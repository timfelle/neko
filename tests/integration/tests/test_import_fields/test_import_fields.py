"""Integration coverage for all case-file users of ``import_fields``.

Flow and scalar initial conditions, sponge baseflows, and non-normal boundary
data all expose field import through the case file.  Exercise each pathway on
the original mesh and through global interpolation.  Interpolation is checked
both with coordinates embedded in the data file and with coordinates taken
from an explicitly named mesh field, the data file having none.

Only the coordinates may be taken from the mesh field: a field requested from
a data file that does not contain it must stop the run, even if the mesh field
contains it.

ALE base-shape import is the only other production caller.  It does not expose
interpolation and its same-mesh path is covered by ``test_ale.py``.
"""

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

import conftest
from testlib import (
    configure_nprocs,
    get_genmeshbox,
    get_neko,
    run_neko,
    which_command,
)

COMPONENTS = ("fluid", "scalar", "non_normal", "sponge")
IMPORT_MODES = ("same_mesh", "embedded_mesh", "separate_mesh")
NPROCS = 2
INTERPOLATION_TOLERANCE = {"dp": 1.0e-8, "sp": 1.0e-4}
VALUE_TOLERANCE = {"dp": 1.0e-8, "sp": 2.0e-4}
# Position of the rdcode, which lists the fields of an fld file, in the
# 132-byte header.
RDCODE = slice(83, 93)
# Pressure of the mesh field, which the data file of the missing-field test
# does not contain.
MESH_PRESSURE = 5.0


def _solver(solver_type):
    """Return a small, robust linear-solver configuration."""
    return {
        "type": solver_type,
        "preconditioner": {"type": "jacobi"},
        "absolute_tolerance": 1.0e-7,
        "max_iterations": 100,
    }


def _field_options(field_file, mode, mesh_field):
    """Build the shared case-file field import options."""
    options = {
        "file_name": str(field_file),
        "interpolate": mode != "same_mesh",
    }
    if mode == "separate_mesh":
        options["mesh_file_name"] = str(mesh_field)
    if mode != "same_mesh":
        options["interpolation"] = {
            "tolerance": INTERPOLATION_TOLERANCE[conftest.RP],
            "padding": 0.05,
        }
    return options


def _write_fields(source, target, codes, pressure=None):
    """Write some of the fields of an fld file to a new fld file.

    `codes` lists the fields to keep in file order, e.g. "XU" for the
    coordinates and the velocity. If `pressure` is given, the pressure field
    takes that value everywhere. The bounding boxes at the end of the source
    file, which Neko does not read, are not copied.
    """
    data = Path(source).read_bytes()
    header = data[:132].decode("ascii")
    word_size, lx, ly, lz, n_elements = map(int, header.split()[1:6])
    rdcode = header[RDCODE].strip()
    assert set(rdcode) <= set("XUPT"), f"Unexpected fields in {source}"

    # The header is followed by a test pattern and the element ids
    start = 132 + 4 + 4 * n_elements
    scalar_size = lx * ly * lz * n_elements * word_size
    vector_size = (3 if lz > 1 else 2) * scalar_size
    fields = {}
    for code in rdcode:
        size = vector_size if code in "XU" else scalar_size
        fields[code] = data[start:start + size]
        start += size
    if pressure is not None:
        dtype = "<f8" if word_size == 8 else "<f4"
        values = np.full(scalar_size // word_size, pressure, dtype=dtype)
        fields["P"] = values.tobytes()

    with open(target, "wb") as f:
        f.write(header[:RDCODE.start].encode("ascii"))
        f.write(codes.ljust(RDCODE.stop - RDCODE.start).encode("ascii"))
        f.write(header[RDCODE.stop:].encode("ascii"))
        f.write(data[132:132 + 4 + 4 * n_elements])
        for code in codes:
            f.write(fields[code])


def _fluid(initial_condition, boundary_conditions, source_terms=None):
    """Build the fluid section shared by source and consumer cases."""
    fluid = {
        "scheme": "pnpn",
        "rho": 1.0,
        "mu": 1.0,
        "freeze": True,
        "initial_condition": initial_condition,
        "velocity_solver": _solver("cg"),
        "pressure_solver": _solver("gmres"),
        "boundary_conditions": boundary_conditions,
        "output_control": "never",
    }
    if source_terms is not None:
        fluid["source_terms"] = source_terms
    return fluid


def _scalar(initial_condition):
    """Build a passive scalar that preserves its imported constant value."""
    return {
        "enabled": True,
        "name": "temperature",
        "advection": False,
        "lambda": 1.0,
        "cp": 1.0,
        "initial_condition": initial_condition,
        "solver": _solver("cg"),
        "boundary_conditions": [
            {
                "type": "neumann",
                "zone_indices": [1, 2, 3, 4, 5, 6],
                "flux": 0.0,
            }
        ],
    }


def _base_case(mesh, output_directory):
    """Return the inexpensive one-step case shared by all runs."""
    return {
        "version": 1,
        "case": {
            "mesh_file": str(mesh),
            "output_directory": str(output_directory),
            "output_at_end": False,
            "output_boundary": False,
            "output_checkpoints": False,
            "time": {"end_time": 1.0e-3, "timestep": 1.0e-3},
            "numerics": {
                "time_order": 1,
                "polynomial_order": 3,
                "dealias": False,
            },
        },
    }


def _probe(fields, output_file, coordinates):
    """Build a one-point probe used to verify imported values."""
    return {
        "type": "probes",
        "compute_control": "nsamples",
        "compute_value": 1,
        "append_output": False,
        "output_file": output_file,
        "fields": fields,
        "points": [{"type": "points", "coordinates": coordinates}],
    }


def _source_case(mesh, workdir):
    """Create velocity and scalar source data with known values."""
    config = _base_case(mesh, workdir)
    case = config["case"]
    case["output_at_end"] = True
    case["fluid"] = _fluid(
        {
            "type": "expression",
            "value": ["1 + x", "2 + y", "3 + z"],
        },
        [
            {
                "type": "velocity_value",
                "zone_indices": [1, 2, 3, 4, 5, 6],
                "value": [1.5, 2.5, 3.5],
            }
        ],
    )
    case["fluid"]["output_filename"] = "source"
    case["fluid"]["output_mesh_in_all_files"] = True
    case["scalar"] = _scalar({"type": "uniform", "value": 4.0})
    return config


def _consumer_case(component, mode, assets, run_dir):
    """Create one consumer case and its expected probe values."""
    mesh = assets["source_mesh"] if mode == "same_mesh" else assets["target_mesh"]
    if mode == "separate_mesh":
        field = assets["field_without_mesh"]
    else:
        field = assets["field"]
    options = _field_options(field, mode, assets["mesh_field"])
    config = _base_case(mesh, run_dir)
    case = config["case"]

    no_slip = {
        "type": "no_slip",
        "zone_indices": [1, 2, 3, 4, 5, 6],
    }
    fluid_ic = {"type": "uniform", "value": [0.0, 0.0, 0.0]}
    boundaries = [no_slip]
    source_terms = None

    if component == "fluid":
        fluid_ic = {"type": "field", **options}
        probe_fields = ["u", "v", "w"]
        probe_point = [0.5, 0.5, 0.5]
        expected = [1.5, 2.5, 3.5]
    elif component == "scalar":
        # The first scalar in fluid output occupies the temperature slot,
        # whose import index is zero.
        scalar_ic = {"type": "field", "target_index": 0, **options}
        case["scalar"] = _scalar(scalar_ic)
        probe_fields = ["temperature"]
        probe_point = [0.5, 0.5, 0.5]
        expected = [4.0]
    elif component == "non_normal":
        boundaries = [
            {"type": "normal_outflow", "zone_indices": [2], **options},
            {"type": "no_slip", "zone_indices": [1, 3, 4, 5, 6]},
        ]
        probe_fields = ["u", "v", "w"]
        probe_point = [1.0, 0.5, 0.5]
        expected = [0.0, 2.5, 3.5]
    elif component == "sponge":
        source_terms = [
            {
                "type": "sponge",
                "amplitudes": [0.0, 0.0, 0.0],
                "fringe_registry_name": "u",
                "baseflow": {"method": "field", **options},
            }
        ]
        probe_fields = ["sponge_bf_u", "sponge_bf_v", "sponge_bf_w"]
        probe_point = [0.5, 0.5, 0.5]
        expected = [1.5, 2.5, 3.5]
    else:  # pragma: no cover - protects additions to COMPONENTS
        raise ValueError(f"Unknown import component: {component}")

    case["fluid"] = _fluid(fluid_ic, boundaries, source_terms)
    if component == "non_normal":
        # A frozen fluid does not apply boundary conditions during the time
        # loop. Take one non-advecting step so the imported tangential values
        # are observable on the boundary probe.
        case["fluid"]["freeze"] = False
        case["fluid"]["advection"] = False
    probe_name = f"probe_{component}_{mode}.csv"
    case["simulation_components"] = [_probe(probe_fields, probe_name, probe_point)]
    return config, run_dir / probe_name, expected


def _resolve_mesh_checker(neko):
    """Locate mesh_checker next to Neko or on PATH."""
    checker = which_command("mesh_checker")
    if checker is not None:
        return Path(checker).resolve()
    candidate = Path(neko).resolve().with_name("mesh_checker")
    assert candidate.is_file(), "The mesh_checker executable could not be found"
    return candidate


def _generate_mesh(genmeshbox, mesh_checker, workdir, name, dimensions):
    """Generate and validate one non-periodic box mesh."""
    result = subprocess.run(
        [
            str(genmeshbox),
            "0",
            "1",
            "0",
            "1",
            "0",
            "1",
            *(str(value) for value in dimensions),
            ".false.",
            ".false.",
            ".false.",
        ],
        cwd=workdir,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout

    generated = workdir / "box.nmsh"
    mesh = workdir / name
    generated.rename(mesh)
    result = subprocess.run(
        [str(mesh_checker), mesh.name],
        cwd=workdir,
        capture_output=True,
        text=True,
        errors="replace",
    )
    assert result.returncode == 0, result.stdout
    return mesh


@pytest.fixture(scope="module")
def import_assets(tmp_path_factory, request):
    """Generate meshes and a single source field shared by the matrix."""
    workdir = tmp_path_factory.mktemp("import_fields")
    neko = Path(get_neko()).resolve()
    genmeshbox = Path(get_genmeshbox()).resolve()
    mesh_checker = _resolve_mesh_checker(neko)
    launcher = Path(request.config.getoption("--launcher-script")).resolve()

    source_mesh = _generate_mesh(
        genmeshbox, mesh_checker, workdir, "source.nmsh", (2, 2, 2)
    )
    target_mesh = _generate_mesh(
        genmeshbox, mesh_checker, workdir, "target.nmsh", (3, 2, 2)
    )

    source_case = workdir / "source.case"
    source_case.write_text(
        json.dumps(_source_case(source_mesh, workdir), indent=2) + "\n",
        encoding="utf-8",
    )
    source_log = workdir / "source.log"
    result = run_neko(
        str(launcher),
        configure_nprocs(NPROCS),
        str(source_case),
        str(neko),
        str(source_log),
    )
    assert result.returncode == 0, source_log.read_text(encoding="utf-8")

    field = workdir / "source0.f00000"
    assert field.is_file(), "The source run did not write source0.f00000"
    mesh_field = workdir / "source_mesh0.f00000"
    shutil.copyfile(field, mesh_field)
    # A data file whose coordinates have to be taken from the mesh field
    field_without_mesh = workdir / "source_data0.f00000"
    _write_fields(field, field_without_mesh, "UPT")

    return {
        "workdir": workdir,
        "neko": neko,
        "launcher": launcher,
        "source_mesh": source_mesh,
        "target_mesh": target_mesh,
        "field": field,
        "field_without_mesh": field_without_mesh,
        "mesh_field": mesh_field,
    }


@pytest.mark.parametrize("component", COMPONENTS)
@pytest.mark.parametrize("mode", IMPORT_MODES)
def test_case_file_import_fields(component, mode, import_assets):
    """Import known data through every case-file pathway and verify values."""
    run_dir = import_assets["workdir"] / f"{component}_{mode}"
    run_dir.mkdir()
    config, probe_file, expected = _consumer_case(
        component, mode, import_assets, run_dir
    )
    case_file = run_dir / "consumer.case"
    case_file.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    log_file = run_dir / "consumer.log"

    result = run_neko(
        str(import_assets["launcher"]),
        configure_nprocs(NPROCS),
        str(case_file),
        str(import_assets["neko"]),
        str(log_file),
    )
    log = log_file.read_text(encoding="utf-8")
    assert result.returncode == 0, log
    assert "Import fields" in log
    assert f"Interpolation : {'T' if mode != 'same_mesh' else 'F'}" in log
    if mode == "separate_mesh":
        assert f"Mesh file     : {import_assets['mesh_field']}" in log

    values = [
        float(value) for value in probe_file.read_text().splitlines()[-1].split(",")
    ]
    assert values[1:] == pytest.approx(expected, abs=VALUE_TOLERANCE[conftest.RP])


def test_field_missing_from_data_file(import_assets):
    """A field missing from the data file must not come from the mesh field."""
    run_dir = import_assets["workdir"] / "missing_field"
    run_dir.mkdir()
    # The mesh field has a pressure, the data file has only the velocity
    mesh_field = run_dir / "pressure_mesh0.f00000"
    field = run_dir / "velocity0.f00000"
    _write_fields(import_assets["field"], mesh_field, "XUP",
                  pressure=MESH_PRESSURE)
    _write_fields(import_assets["field"], field, "U")

    config = _base_case(import_assets["target_mesh"], run_dir)
    options = _field_options(field, "separate_mesh", mesh_field)
    no_slip = {"type": "no_slip", "zone_indices": [1, 2, 3, 4, 5, 6]}
    case = config["case"]
    case["fluid"] = _fluid({"type": "field", **options}, [no_slip])
    probe_file = run_dir / "probe_missing_field.csv"
    case["simulation_components"] = [
        _probe(["p"], probe_file.name, [0.5, 0.5, 0.5])
    ]
    case_file = run_dir / "consumer.case"
    case_file.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    log_file = run_dir / "consumer.log"

    result = run_neko(
        str(import_assets["launcher"]),
        configure_nprocs(NPROCS),
        str(case_file),
        str(import_assets["neko"]),
        str(log_file),
    )
    if result.returncode == 0:
        probe = probe_file.read_text().splitlines()[-1]
        pytest.fail(
            "Importing a pressure the data file does not contain did not fail,"
            f" the pressure was taken from the mesh field (t, p: {probe})"
        )
    log = log_file.read_text(encoding="utf-8")
    assert "Field p requested for import is not in the fld file" in log, log
