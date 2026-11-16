"""Test that a restart continues the output schedule of the Lagrangian
particle tracking component, without adding or skipping a snapshot.

The particle tracker writes its trajectory snapshots with an output
controller of its own, a copy of the one of the simulation component. The
case runs from 0 to 0.12 with a time step of 0.01, in a frozen flow, and
writes a snapshot every 0.04 and a checkpoint every 0.05. The second run
restarts from the checkpoint at t = 0.05, so that its first step, to
t = 0.06, does not end on a scheduled time.

Only the snapshots written in the time loop are compared. The snapshot
written at the initialisation, with time step 0, is not part of the
schedule.
"""

import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from testlib import configure_nprocs, get_neko, get_neko_dir, run_neko


def _write_case(path, case, output_dir, mesh, restart_file=None):
    """Write the case file, pointing it at the working directory."""
    case = json.loads(json.dumps(case))
    case["case"]["mesh_file"] = str(mesh)
    case["case"]["output_directory"] = str(output_dir)
    if restart_file is not None:
        case["case"]["restart_file"] = str(restart_file)
    with path.open("w", encoding="utf-8") as stream:
        json.dump(case, stream, indent=2)
    return path


def _snapshot_times(directory):
    """The times of the trajectory snapshots written in the time loop."""
    files = sorted(directory.glob("tracers_*.csv"))
    assert files, f"No trajectory output in {directory}"
    times = []
    for path in files:
        with path.open(encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                if int(float(row["tstep"])) < 1:
                    continue
                time = round(float(row["time"]), 8)
                if time not in times:
                    times.append(time)
    return times


def test_lpt_restart(launcher_script, request, tmp_path):
    """Check the snapshots written by a run and by a restart of it."""
    neko = get_neko()
    neko_dir = Path(get_neko_dir()).resolve()
    test_dir = (
        neko_dir / "tests" / "integration" / "tests" / "test_lpt_restart"
    )
    mesh = (
        neko_dir / "tests" / "integration" / "meshes" / "small_test_cyl.nmsh"
    )
    nprocs = configure_nprocs(1)

    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)

    # Keep the restart path short enough for Neko's fixed-size log buffer.
    work_context = TemporaryDirectory(prefix="neko-lpt-", dir="/tmp")
    request.addfinalizer(work_context.cleanup)
    work_dir = Path(work_context.name)

    with (test_dir / "lpt_restart.case").open(encoding="utf-8") as stream:
        case = json.load(stream)

    #
    # Part 1, the whole interval in one run
    #
    part1_dir = work_dir / "part1"
    part1_dir.mkdir()
    part1_case = _write_case(tmp_path / "part1.case", case, part1_dir, mesh)
    part1_log = log_dir / f"{request.node.name}_part1.log"
    result = run_neko(
        launcher_script, nprocs, str(part1_case), neko, str(part1_log)
    )
    assert result.returncode == 0, (
        f"neko process failed with exit code {result.returncode}"
    )

    restart_time = 0.05
    part1_times = _snapshot_times(part1_dir)
    after_restart = [t for t in part1_times if t > restart_time + 1e-8]
    # The snapshots the restarted run has to write.
    assert after_restart == [0.08, 0.12], (
        f"Unexpected snapshots in the whole run: {part1_times}"
    )

    #
    # Part 2, the second half of the same interval, restarted from the
    # checkpoint written at t = 0.05
    #
    part2_dir = work_dir / "part2"
    part2_dir.mkdir()
    part2_case = _write_case(
        tmp_path / "part2.case",
        case,
        part2_dir,
        mesh,
        restart_file=part1_dir / "chkp00000.chkp",
    )
    part2_log = log_dir / f"{request.node.name}_part2.log"
    result = run_neko(
        launcher_script, nprocs, str(part2_case), neko, str(part2_log)
    )
    assert result.returncode == 0, (
        f"neko process failed with exit code {result.returncode}"
    )

    # The schedule of the whole run is continued: no snapshot at the first
    # step after the restart, and none skipped.
    part2_times = _snapshot_times(part2_dir)
    assert part2_times == after_restart, (
        f"Snapshots after the restart: {part2_times}, "
        f"expected {after_restart}"
    )
