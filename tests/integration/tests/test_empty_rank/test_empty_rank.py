"""A run with an MPI rank that holds no elements must end cleanly.

A single-element mesh on two ranks leaves the second rank without elements, so
its coefficient arrays have zero size.  The finalisation used to free the mass
matrix storage of such a rank three times: once as B and once through each of
the pointers to its lagged copies, which alias B unless a moving mesh gives
them storage of their own.

The element is periodic in all three directions.  That is the configuration in
which the empty rank reaches the finalisation; a non-periodic single element
on two ranks stops earlier, in the gather-scatter setup of the occupied rank,
which this test does not cover.
"""

import json
import subprocess

import pytest

from testlib import configure_nprocs, get_genmeshbox, get_neko, run_neko

NPROCS = 2


def _case(mesh, output_directory):
    """Return a two-step fluid case on the given mesh."""
    solver = {
        "type": "cg",
        "preconditioner": {"type": "jacobi"},
        "absolute_tolerance": 1.0e-7,
        "max_iterations": 100,
    }
    return {
        "version": 1,
        "case": {
            "mesh_file": str(mesh),
            "output_directory": str(output_directory),
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
                "boundary_conditions": [],
                "output_control": "never",
            },
        },
    }


def test_empty_rank(launcher_script, log_file, tmp_path):
    """Run on two ranks with one element and check the run ends cleanly."""
    nprocs = configure_nprocs(NPROCS)
    if nprocs < NPROCS:
        pytest.skip("Leaving a rank without elements needs two ranks")

    result = subprocess.run(
        [get_genmeshbox(), "0", "1", "0", "1", "0", "1", "1", "1", "1",
         ".true.", ".true.", ".true."],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    mesh = tmp_path / "box.nmsh"

    case_file = tmp_path / "empty_rank.case"
    case_file.write_text(json.dumps(_case(mesh, tmp_path), indent=2))

    result = run_neko(launcher_script, nprocs, str(case_file), get_neko(),
                      log_file)
    with open(log_file, errors="replace") as f:
        log = f.read()

    # The check is only meaningful if a rank actually holds no elements
    assert "has zero elements" in log, \
        "No rank was left without elements, so the test checks nothing"
    assert result.returncode == 0, \
        f"neko failed with exit code {result.returncode}:\n{log[-3000:]}"
    assert "Normal end" in log
