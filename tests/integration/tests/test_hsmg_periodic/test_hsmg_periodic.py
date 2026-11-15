"""Regression test for the hsmg preconditioner on a fully periodic box.

Without a pressure Dirichlet condition the coarse problem of hsmg is singular,
with the constants as its null space. Its right-hand side picks up a component
in that null space from round-off in GMRES, and once that component is above
the coarse CG tolerance, the coarse solve diverges and the pressure goes NaN.
"""

import json
import re
import subprocess

import pytest

import conftest
from testlib import configure_nprocs, get_genmeshbox, get_neko, run_neko


# The null space component grows with every GMRES iteration, so with a tight
# pressure tolerance a coarse solve that does not remove it breaks down
# already in step 1.
PRESSURE_TOLERANCE = 1.0e-8
TWO_PI = "6.283185307179586"


def _case(mesh, output_directory, lx):
    """Return three steps of a Taylor-Green vortex on a periodic box."""
    return {
        "version": 1,
        "case": {
            "mesh_file": str(mesh),
            "output_directory": str(output_directory),
            "output_at_end": False,
            "output_boundary": False,
            "output_checkpoints": False,
            "time": {"end_time": 0.003, "timestep": 0.001},
            "numerics": {
                "time_order": 3,
                "polynomial_order": lx - 1,
                "dealias": False,
            },
            "fluid": {
                "scheme": "pnpn",
                "Re": 100.0,
                "output_control": "never",
                "initial_condition": {
                    "type": "expression",
                    "value": [
                        "sin(x)*cos(y)*cos(z)",
                        "-cos(x)*sin(y)*cos(z)",
                        "0",
                    ],
                },
                "velocity_solver": {
                    "type": "cg",
                    "preconditioner": {"type": "jacobi"},
                    "projection_space_size": 0,
                    "absolute_tolerance": 1.0e-8,
                    "max_iterations": 100,
                },
                "pressure_solver": {
                    "type": "gmres",
                    "preconditioner": {"type": "hsmg"},
                    "projection_space_size": 0,
                    "absolute_tolerance": PRESSURE_TOLERANCE,
                    "max_iterations": 100,
                },
            },
        },
    }


def _pressure_solves(log):
    """Iterations and start residual of each pressure solve in the log."""
    return [(int(m.group(1)), float(m.group(2))) for m in
            re.finditer(r"\|\s*Pressure\s+(\d+)\s+(\S+)", log)]


@pytest.mark.parametrize("lx", [8, 16])
def test_hsmg_periodic(launcher_script, log_file, tmp_path, lx):
    """hsmg solves the pure-Neumann pressure problem of a periodic box."""
    if conftest.RP == "sp":
        pytest.skip("In single precision the coarse CG tolerance is at "
                    "round-off, and the singular coarse solve can still "
                    "break down")
    # A periodic direction needs at least 3 elements for a valid mesh.
    result = subprocess.run(
        [
            get_genmeshbox(),
            "0", TWO_PI, "0", TWO_PI, "0", TWO_PI,
            "3", "3", "3",
            ".true.", ".true.", ".true.",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout
    mesh = tmp_path / "box.nmsh"

    case_file = tmp_path / f"hsmg_periodic_lx{lx}.case"
    case_file.write_text(json.dumps(_case(mesh, tmp_path, lx), indent=2))

    result = run_neko(launcher_script, configure_nprocs(2), str(case_file),
                      get_neko(), log_file)
    assert result.returncode == 0, f"Neko failed, see {log_file}"

    with open(log_file) as f:
        log = f.read()
    solves = _pressure_solves(log)
    assert len(solves) == 3, f"Expected 3 pressure solves, see {log_file}"
    for iterations, _ in solves:
        assert iterations <= 50, f"Pressure solve stalled, see {log_file}"
    # The pressure changes little over a step, so with the previous pressure
    # as the initial guess, the later solves start from a much smaller
    # residual. A diverged coarse solve spoils the solution of the first.
    for _, start_residual in solves[1:]:
        assert start_residual < 0.1 * solves[0][1], \
            f"Pressure solve inaccurate, see {log_file}"
    assert "Coarse null space   : constants" in log
