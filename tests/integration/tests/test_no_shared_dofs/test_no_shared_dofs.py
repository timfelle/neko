"""A rank holding elements but sharing no dofs with another rank runs.

A single element on two ranks leaves the rank holding it without any dof
shared with the other rank. On more than one rank, such a rank still runs the
shared part of every gather-scatter operation, which was handed block arrays
that were only allocated when a rank had shared dofs, and the run crashed
while it was set up.

The element is not periodic: a periodic element gives the rank holding it
shared dofs.
"""

import json
import subprocess

import pytest

from testlib import configure_nprocs, get_genmeshbox, get_neko, run_neko

NPROCS = 2


def write_case(tmp_path, mesh):
    """Two time steps of a uniform flow in a box with no-slip walls."""
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
            "output_directory": str(tmp_path),
            "output_at_end": False,
            "output_boundary": False,
            "output_checkpoints": False,
            "time": {"end_time": 2.0e-3, "timestep": 1.0e-3},
            "numerics": {
                "time_order": 1,
                "polynomial_order": 3,
                "dealias": False,
            },
            "fluid": {
                "scheme": "pnpn",
                "rho": 1.0,
                "mu": 1.0,
                "initial_condition": {
                    "type": "uniform",
                    "value": [1.0, 0.0, 0.0],
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
    case_file = tmp_path / "no_shared_dofs.case"
    case_file.write_text(json.dumps(case, indent=2))
    return case_file


def test_no_shared_dofs(launcher_script, log_file, tmp_path):
    """Run one element on two ranks through the setup and the time steps."""
    nprocs = configure_nprocs(NPROCS)
    if nprocs < NPROCS:
        pytest.skip("Leaving a rank without shared dofs needs two ranks")

    result = subprocess.run(
        [get_genmeshbox(), "0", "1", "0", "1", "0", "1", "1", "1", "1",
         ".false.", ".false.", ".false."],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    mesh = tmp_path / "box.nmsh"

    case_file = write_case(tmp_path, mesh)
    result = run_neko(launcher_script, nprocs, str(case_file), get_neko(),
                      log_file)
    with open(log_file, errors="replace") as f:
        log = f.read()

    # The element is on one rank only, so the other one holds none
    assert "has zero elements" in log, \
        "Every rank holds elements, so the test checks nothing"
    # Not the exit status: that also covers the finalisation of the rank
    # without elements, which this test is not about
    assert "Normal end" in log, \
        f"neko failed with exit code {result.returncode}:\n{log[-3000:]}"
