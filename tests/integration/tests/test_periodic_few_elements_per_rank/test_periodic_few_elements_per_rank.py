"""A periodic mesh runs on ranks that hold only a few of its elements.

Reading the periodic zones gives the points of a periodic facet the ids of
the points they are matched with. A rank holding one end of a periodic
direction but not the other has not seen these ids before, and used to store
an extra copy of each such point. The point storage of a rank is sized for
the points of its own elements and has no room for the copies, so with one
or two elements per rank they overflowed it, and the run crashed while the
mesh was read.
"""

import json
import subprocess

import pytest

import conftest
from testlib import configure_nprocs, get_genmeshbox, get_neko, run_neko

SOLVER_TOLERANCE = {"dp": 1.0e-10, "sp": 1.0e-5}

# Elements and periodicity per direction, and the number of ranks. The
# elements are shared out in order, so the last rank holds the far end of
# each periodic direction, or part of it, without the near end.
CASES = {
    "3x1x1-periodic-x-2-ranks": ((3, 1, 1), (True, False, False), 2),
    "3x1x1-periodic-x-3-ranks": ((3, 1, 1), (True, False, False), 3),
    "3x3x3-periodic-8-ranks": ((3, 3, 3), (True, True, True), 8),
}


def generate_box(tmp_path, nel, periodic):
    """Generate a box on the unit cube with nel elements per direction."""
    flags = [".true." if p else ".false." for p in periodic]
    result = subprocess.run(
        [get_genmeshbox(), "0", "1", "0", "1", "0", "1",
         *(str(n) for n in nel), *flags],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    mesh = tmp_path / "box.nmsh"
    assert mesh.is_file()
    return mesh


def write_case(tmp_path, mesh, periodic):
    """Two time steps of a uniform flow, no-slip walls where not periodic."""
    tolerance = SOLVER_TOLERANCE[conftest.RP]
    walls = [zone for d, p in enumerate(periodic) if not p
             for zone in (2 * d + 1, 2 * d + 2)]
    solver = {
        "type": "cg",
        "preconditioner": {"type": "jacobi"},
        "absolute_tolerance": tolerance,
        "max_iterations": 200,
    }
    fluid = {
        "scheme": "pnpn",
        "mu": 0.1,
        "rho": 1.0,
        "output_control": "never",
        "initial_condition": {"type": "uniform", "value": [1.0, 0.0, 0.0]},
        "velocity_solver": solver,
        "pressure_solver": dict(solver, type="gmres"),
        "boundary_conditions": [],
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
                "polynomial_order": 3,
                "dealias": False,
            },
            "fluid": fluid,
        },
    }
    case_file = tmp_path / "periodic.case"
    case_file.write_text(json.dumps(case, indent=2))
    return case_file


@pytest.mark.parametrize("nel, periodic, nprocs", CASES.values(),
                         ids=CASES.keys())
def test_periodic_few_elements_per_rank(launcher_script, log_file, tmp_path,
                                        nel, periodic, nprocs):
    """The run ends cleanly with a few elements per rank."""
    if configure_nprocs(nprocs) < nprocs:
        pytest.skip(f"Needs {nprocs} ranks")

    mesh = generate_box(tmp_path, nel, periodic)
    case_file = write_case(tmp_path, mesh, periodic)

    result = run_neko(launcher_script, nprocs, str(case_file), get_neko(),
                      log_file)
    with open(log_file, errors="replace") as f:
        log = f.read()

    assert result.returncode == 0, \
        f"neko failed with exit code {result.returncode}:\n{log[-3000:]}"
    assert "Normal end" in log
