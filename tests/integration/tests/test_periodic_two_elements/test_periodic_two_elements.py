"""A periodic direction with 2 elements is rejected when the mesh is set up.

Neko identifies an edge or a face by the ids of its vertices. With 2 elements
in a periodic direction, the two elements of each periodic column have the
same vertices, so their edges and faces cannot be told apart and the
connectivity would be wrong. A single element, or 3 and more, are accepted.
"""

import json
import subprocess

import pytest

import conftest
from testlib import configure_nprocs, get_genmeshbox, get_neko, run_neko


ERROR_MESSAGE = "A periodic direction with 2 elements is not supported"
POLYNOMIAL_ORDER = 3
SOLVER_TOLERANCE = {"dp": 1.0e-10, "sp": 1.0e-5}
TWO_PI = "6.283185307179586"


def generate_box(tmp_path, nel, periodic):
    """Generate a box on [0, 2 pi]^3 with nel elements per direction."""
    flags = [".true." if p else ".false." for p in periodic]
    result = subprocess.run(
        [get_genmeshbox(), "0", TWO_PI, "0", TWO_PI, "0", TWO_PI,
         *(str(n) for n in nel), *flags],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout
    mesh = tmp_path / "box.nmsh"
    assert mesh.is_file()
    return mesh


def write_case(tmp_path, mesh, periodic):
    """Two time steps of a uniform flow, no-slip walls where not periodic."""
    tolerance = SOLVER_TOLERANCE[conftest.RP]
    walls = [zone for d, p in enumerate(periodic) if not p
             for zone in (2 * d + 1, 2 * d + 2)]
    fluid = {
        "scheme": "pnpn",
        "mu": 0.1,
        "rho": 1.0,
        "output_control": "never",
        "initial_condition": {"type": "uniform", "value": [1.0, 0.0, 0.0]},
        "velocity_solver": {
            "type": "cg",
            "preconditioner": {"type": "jacobi"},
            "absolute_tolerance": tolerance,
            "max_iterations": 200,
        },
        "pressure_solver": {
            "type": "gmres",
            "preconditioner": {"type": "jacobi"},
            "absolute_tolerance": tolerance,
            "max_iterations": 200,
        },
    }
    if walls:
        fluid["boundary_conditions"] = [
            {"type": "no_slip", "zone_indices": walls}
        ]

    case = {
        "version": 1,
        "case": {
            "mesh_file": str(mesh),
            "output_directory": str(tmp_path),
            "output_at_end": False,
            "output_boundary": False,
            "output_checkpoints": False,
            "time": {"end_time": 0.02, "timestep": 0.01},
            "numerics": {
                "time_order": 1,
                "polynomial_order": POLYNOMIAL_ORDER,
                "dealias": False,
            },
            "fluid": fluid,
        },
    }
    case_file = tmp_path / "box.case"
    case_file.write_text(json.dumps(case, indent=2))
    return case_file


def unique_points(nel, periodic):
    """The number of distinct GLL points of the box."""
    count = 1
    for n, p in zip(nel, periodic):
        count *= n * POLYNOMIAL_ORDER + (0 if p else 1)
    return count


@pytest.mark.parametrize("nprocs", [1, 2])
@pytest.mark.parametrize(
    "nel, periodic",
    [
        ((2, 2, 2), (True, True, True)),
        ((2, 3, 3), (True, False, False)),
        # On 2 ranks each periodic column of 2 is split between the ranks
        ((3, 3, 2), (False, False, True)),
    ],
    ids=["2x2x2-periodic", "2x3x3-periodic-x", "3x3x2-periodic-z"],
)
def test_two_elements_rejected(
    launcher_script, log_file, tmp_path, nel, periodic, nprocs
):
    """Setting up the mesh stops with an error naming the restriction."""
    mesh = generate_box(tmp_path, nel, periodic)
    case_file = write_case(tmp_path, mesh, periodic)

    result = run_neko(launcher_script, configure_nprocs(nprocs),
                      str(case_file), get_neko(), log_file)

    with open(log_file, errors="replace") as f:
        log = f.read()
    assert result.returncode != 0, "neko ran on a mesh it cannot represent"
    assert ERROR_MESSAGE in log


@pytest.mark.parametrize("nprocs", [1, 2])
@pytest.mark.parametrize(
    "nel, periodic",
    [
        ((3, 3, 3), (True, True, True)),
        ((2, 2, 2), (False, False, False)),
        # One element in a periodic direction, as in quasi-2D cases
        ((3, 3, 1), (False, False, True)),
    ],
    ids=["3x3x3-periodic", "2x2x2-walls", "3x3x1-periodic-z"],
)
def test_other_meshes_run(
    launcher_script, log_file, tmp_path, nel, periodic, nprocs
):
    """Valid meshes run and have the expected number of distinct points."""
    mesh = generate_box(tmp_path, nel, periodic)
    case_file = write_case(tmp_path, mesh, periodic)

    result = run_neko(launcher_script, configure_nprocs(nprocs),
                      str(case_file), get_neko(), log_file)

    with open(log_file, errors="replace") as f:
        log = f.read()
    assert result.returncode == 0, f"neko failed with {result.returncode}"
    assert ERROR_MESSAGE not in log
    assert f"Unique pts.: {unique_points(nel, periodic)}\n" in log
