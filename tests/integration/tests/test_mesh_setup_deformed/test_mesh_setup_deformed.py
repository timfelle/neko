"""Elements deformed by a user mesh_setup routine are flagged as deformed.

The deformed-element flags decide where the cross terms of the metric are
applied, for example by the libxsmm Helmholtz operator. They are made when the
mesh is read, but the user mesh_setup routine runs later and may move the
vertices. The user file shears the box, y in proportion to x, which deforms
every element, and logs the number of elements flagged once the case is
initialised. Without the shear no element is flagged, which shows that the
count does not come from something that flags every element, as the hsmg
preconditioner does.
"""

import json
import re
import subprocess
from pathlib import Path

import pytest

from testlib import configure_nprocs, get_genmeshbox, get_makeneko, run_neko

USER_FILE = Path("tests", "test_mesh_setup_deformed",
                 "test_mesh_setup_deformed.f90")

# Number of elements of the 3 x 3 x 3 box
NELV = 27


@pytest.fixture(scope="module")
def user_neko(tmp_path_factory):
    """Build the user neko and the mesh once for all the tests."""
    build_dir = tmp_path_factory.mktemp("mesh_setup_deformed")

    result = subprocess.run(
        [get_makeneko(), str(USER_FILE.resolve())],
        cwd=build_dir, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    result = subprocess.run(
        [get_genmeshbox(), "0", "1", "0", "1", "0", "1", "3", "3", "3",
         ".false.", ".false.", ".false."],
        cwd=build_dir, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    return build_dir / "neko", build_dir / "box.nmsh"


def write_case(tmp_path, mesh, shear):
    """One time step of a frozen fluid in a box with no-slip walls. Only
    Jacobi preconditioners, as the hsmg one flags every element."""
    solver = {
        "type": "cg",
        "preconditioner": {"type": "jacobi"},
        "absolute_tolerance": 1.0e-7,
        "max_iterations": 100,
    }
    case = {
        "version": 1,
        "case": {
            "mesh_file": str(mesh),
            "mesh_shear": shear,
            "output_directory": str(tmp_path),
            "output_at_end": False,
            "output_boundary": False,
            "output_checkpoints": False,
            "time": {"end_time": 1.0e-3, "timestep": 1.0e-3},
            "numerics": {
                "time_order": 1,
                "polynomial_order": 3,
                "dealias": False,
            },
            "fluid": {
                "scheme": "pnpn",
                "rho": 1.0,
                "mu": 1.0,
                "freeze": True,
                "initial_condition": {
                    "type": "uniform",
                    "value": [0.0, 0.0, 0.0],
                },
                "velocity_solver": solver,
                "pressure_solver": dict(solver, type="gmres"),
                "boundary_conditions": [
                    {"type": "no_slip", "zone_indices": [1, 2, 3, 4, 5, 6]}
                ],
                "output_control": "never",
            },
        },
    }
    case_file = tmp_path / "mesh_setup_deformed.case"
    case_file.write_text(json.dumps(case, indent=2))
    return case_file


@pytest.mark.parametrize("shear, expected", [(0.0, 0), (0.5, NELV)],
                         ids=["box", "sheared"])
def test_mesh_setup_deformed(launcher_script, log_file, tmp_path, user_neko,
                             shear, expected):
    """Count the elements flagged as deformed after the user mesh_setup."""
    neko, mesh = user_neko
    case_file = write_case(tmp_path, mesh, shear)

    result = run_neko(launcher_script, configure_nprocs(1), str(case_file),
                      str(neko), log_file)
    with open(log_file, errors="replace") as f:
        log = f.read()
    assert result.returncode == 0, \
        f"neko failed with exit code {result.returncode}:\n{log[-3000:]}"

    match = re.search(r"Deformed elements: (\d+) of (\d+)", log)
    assert match, "The user file did not log the deformed elements"
    assert int(match.group(2)) == NELV
    assert int(match.group(1)) == expected, \
        f"{match.group(1)} of {NELV} elements flagged as deformed, " \
        f"expected {expected}"
