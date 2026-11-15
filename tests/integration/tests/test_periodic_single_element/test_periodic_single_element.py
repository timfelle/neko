"""Warn about a single element in two or more periodic directions.

Neko identifies a face or an edge by the ids of its vertices. With a single
element in two periodic directions, the faces of an element normal to them
have the same vertices, as have its edges along them, so they are joined.
Only fields that are uniform in these directions, and have no vector
component along them, are then handled correctly, which Neko warns about
when it sets up the mesh. A single element in one periodic direction, as in
pseudo-2D cases, is exact and gives no warning.
"""

import json
import subprocess

import pytest

import conftest
from testlib import configure_nprocs, get_genmeshbox, get_neko, run_neko


WARNING_MESSAGE = "Two or more periodic directions with a single element"
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


def write_case(tmp_path, mesh, nel, periodic):
    """Two time steps of a uniform flow, no-slip walls where not periodic.

    The flow is along the first direction with more than one element. On the
    meshes with the warning, a velocity along a direction with a single element
    is not handled, and the pressure solver would not converge.
    """
    tolerance = SOLVER_TOLERANCE[conftest.RP]
    walls = [zone for d, p in enumerate(periodic) if not p
             for zone in (2 * d + 1, 2 * d + 2)]
    axis = next(d for d, n in enumerate(nel) if n > 1)
    velocity = [1.0 if d == axis else 0.0 for d in range(3)]
    fluid = {
        "scheme": "pnpn",
        "mu": 0.1,
        "rho": 1.0,
        "output_control": "never",
        "initial_condition": {"type": "uniform", "value": velocity},
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


def run_box(launcher_script, log_file, tmp_path, nel, periodic, nprocs):
    """Run the case on the box and return the exit status and the log."""
    mesh = generate_box(tmp_path, nel, periodic)
    case_file = write_case(tmp_path, mesh, nel, periodic)

    result = run_neko(launcher_script, configure_nprocs(nprocs),
                      str(case_file), get_neko(), log_file)

    with open(log_file, errors="replace") as f:
        log = f.read()
    return result.returncode, log


@pytest.mark.parametrize("nprocs", [1, 2])
@pytest.mark.parametrize(
    "nel, periodic",
    [
        ((4, 1, 1), (False, True, True)),
        ((1, 4, 1), (True, False, True)),
        ((1, 1, 4), (True, True, False)),
    ],
    ids=["4x1x1-periodic-yz", "1x4x1-periodic-xz", "1x1x4-periodic-xy"],
)
def test_warning_printed(
    launcher_script, log_file, tmp_path, nel, periodic, nprocs
):
    """The case runs, and the warning is printed once, by the first rank."""
    returncode, log = run_box(launcher_script, log_file, tmp_path, nel,
                              periodic, nprocs)

    assert returncode == 0, f"neko failed with {returncode}"
    assert log.count(WARNING_MESSAGE) == 1


@pytest.mark.parametrize("nprocs", [1, 2])
@pytest.mark.parametrize(
    "nel, periodic",
    [
        # One element in a periodic direction, as in pseudo-2D cases
        ((3, 3, 1), (False, False, True)),
        ((3, 3, 1), (True, True, True)),
        ((4, 1, 1), (False, False, False)),
        ((3, 3, 3), (True, True, True)),
    ],
    ids=["3x3x1-periodic-z", "3x3x1-periodic", "4x1x1-walls",
         "3x3x3-periodic"],
)
def test_no_warning(
    launcher_script, log_file, tmp_path, nel, periodic, nprocs
):
    """Other meshes run without the warning, with every point distinct."""
    returncode, log = run_box(launcher_script, log_file, tmp_path, nel,
                              periodic, nprocs)

    assert returncode == 0, f"neko failed with {returncode}"
    assert WARNING_MESSAGE not in log
    assert f"Unique pts.: {unique_points(nel, periodic)}\n" in log
