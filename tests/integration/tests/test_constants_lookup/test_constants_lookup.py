"""Keys that may name an entry of case.constants instead of a literal.

Each case gives one key the name of a constant whose value differs from the
key's default, runs a few steps, and checks that Neko ends normally and, where
it logs the key, logs the constant's value. Keys that Neko does not log are
checked by the run completing only.
One key per reader family is covered: source-term times and values, gradient
jump penalty, multigrid preconditioners, fluid and scalar boundary
conditions, the wall model, the scalar field initial condition, the restart
interpolation tolerance and the compressible solver.
"""
from pathlib import Path
import json
import re
import shutil
import subprocess

import json5
import pytest
from testlib import get_neko, get_makeneko, run_neko, configure_nprocs

TEST_DIR = Path(__file__).resolve().parent
MESH = str(TEST_DIR.parents[1] / "meshes" / "small_test_cyl.nmsh")


def make_case(directory, name, constants, edit):
    """Write the base case with the given constants and edit applied.

    The constants are a dict, or a list of (name, value) pairs when a name
    is repeated."""
    with open(TEST_DIR / "constants_lookup.case") as f:
        case = json5.load(f)
    case_object = case["case"]
    case_object["mesh_file"] = MESH
    case_object["output_directory"] = str(directory)
    if isinstance(constants, dict):
        constants = constants.items()
    case_object["constants"] = [{"name": key, "value": value}
                                for key, value in constants]
    edit(case_object)
    case_file = str(Path(directory) / (name + ".case"))
    with open(case_file, "w") as f:
        json.dump(case, f, indent=2)
    return case_file


def run_neko_log(launcher_script, case_file, log_file, neko=None):
    """Run neko and return its exit code, the log and its error lines."""
    result = run_neko(launcher_script, configure_nprocs(2), case_file,
                      neko or get_neko(), log_file)
    with open(log_file) as f:
        log = f.read()
    errors = sorted({line.strip() for line in log.splitlines()
                     if "*** ERROR" in line})
    return result, log, errors


def run(launcher_script, case_file, log_file, neko=None):
    """Run neko and return the log, which must show a normal end."""
    result, log, errors = run_neko_log(launcher_script, case_file, log_file,
                                       neko)
    assert result.returncode == 0, \
        f"neko failed with exit code {result.returncode}: {errors}"
    assert "Normal end." in log
    return log


def scalar_source(values, **times):
    """A constant source term on the scalar."""
    def edit(case):
        case["scalar"]["source_terms"] = [
            dict({"type": "constant", "values": values}, **times)]
    return edit


def fluid_source(term):
    def edit(case):
        case["fluid"]["source_terms"] = [term]
    return edit


def pressure_preconditioner(preconditioner):
    def edit(case):
        case["fluid"]["pressure_solver"]["preconditioner"] = preconditioner
    return edit


def fluid_key(key, value):
    def edit(case):
        case["fluid"][key] = value
    return edit


def dong(case):
    case["fluid"]["boundary_conditions"][1].update(
        {"delta": "dong_delta", "velocity_scale": "u_ref"})


def scalar_bc_zones(case):
    case["scalar"]["boundary_conditions"][1]["zone_indices"] = "outlets"


def wall_model(case):
    case["fluid"]["full_stress_formulation"] = True
    case["fluid"]["velocity_solver"]["type"] = "coupled_cg"
    case["fluid"]["boundary_conditions"][4] = {
        "type": "wall_model", "model": "spalding", "zone_indices": [7],
        "kappa": 0.41, "B": 5.2, "h_index": "h_idx"}


# (id, constants, edit, pattern the log must contain or None). A None
# pattern means Neko does not log the key, so only the run is checked:
# source term times, the gradient jump penalty, schwarz_iterations, the Dong
# outflow parameters, the scalar boundary zones, the wall model h_index and
# the phmg eigenvalue refresh frequency.
CASES = [
    ("source_start_time", {"t_on": 0.005, "t_off": 1.0},
     scalar_source([0.5], start_time="t_on", end_time="t_off"), None),
    ("constant_values_number", {}, scalar_source(0.5), None),
    ("constant_values_scalar", {"q": 0.5}, scalar_source("q"), None),
    ("constant_values_vector", {"q": [0.5]}, scalar_source("q"), None),
    ("gjp_scaling_factor", {"gjp_a": 0.5, "gjp_t": 0.005},
     fluid_source({"type": "gradient_jump_penalty",
                   "scaling_factor": "gjp_a", "start_time": "gjp_t"}),
     None),
    ("hsmg_coarse_iterations", {"n_crs": 7},
     pressure_preconditioner({"type": "hsmg",
                              "coarse_grid": {"iterations": "n_crs"}}),
     r"Coarse grid iters\.\s*:\s*7\n"),
    ("phmg_smoother_iterations", {"n_smoo": 2, "n_freq": 10},
     pressure_preconditioner({"type": "phmg",
                              "smoother_iterations": "n_smoo",
                              "update": {"eigs": {"frequency": "n_freq"}}}),
     r"Smoother Iters:\s*2\n"),
    ("schwarz_iterations", {"n_schwarz": 1},
     fluid_key("schwarz_iterations", "n_schwarz"), None),
    ("dong_delta", {"dong_delta": 0.02, "u_ref": 2.0}, dong, None),
    ("scalar_bc_zone_indices", {"outlets": [2.0, 5.0, 6.0]},
     scalar_bc_zones, None),
    ("wall_model_h_index", {"h_idx": 1}, wall_model, None),
]


@pytest.mark.parametrize("name, constants, edit, pattern", CASES,
                         ids=[c[0] for c in CASES])
def test_constant_reference(launcher_script, log_file, tmp_path, name,
                            constants, edit, pattern):
    """The key is resolved from case.constants."""
    log = run(launcher_script, make_case(tmp_path, name, constants, edit),
              log_file)
    if pattern is not None:
        assert re.search(pattern, log), f"{pattern!r} not in the log"


def test_constant_values_name_clash(launcher_script, log_file, tmp_path):
    """A constant source term's values naming both a scalar and a vector
    constant is an error."""
    case_file = make_case(tmp_path, "clash", [("c", 0.5), ("c", [0.5])],
                          scalar_source("c"))
    result, log, errors = run_neko_log(launcher_script, case_file, log_file)
    assert result.returncode != 0
    assert "*** ERROR: Constant source term values: c names both a scalar" \
        " and a vector constant ***" in errors, errors


@pytest.fixture(scope="module")
def previous_run(request, tmp_path_factory):
    """A run writing a field file and a checkpoint to start from."""
    directory = tmp_path_factory.mktemp("previous")

    def edit(case):
        case["output_at_end"] = True
        case["output_checkpoints"] = True
        case["checkpoint_control"] = "nsamples"
        case["checkpoint_value"] = 1

    case_file = make_case(directory, "previous", {}, edit)
    run(request.config.getoption("--launcher-script"), case_file,
        str(directory / "previous.log"))
    return directory


def scalar_field_ic(directory, ic):
    """Add a second scalar, #1, started from scalar #0 of the previous run's
    field file with the given initial condition keys."""
    fields = sorted(directory.glob("field0.f*"))
    assert fields, "the previous run wrote no field file"

    def edit(case):
        second = json.loads(json.dumps(case["scalar"]))
        second["name"] = "tracer2"
        second["initial_condition"] = dict(
            {"type": "field", "file_name": str(fields[-1]),
             "target_index": 0}, **ic)
        case["scalars"] = [case.pop("scalar"), second]
    return edit


# target_index defaults to the scalar's own index, 1 here, and the
# interpolation tolerance to about 2e-13.
@pytest.mark.parametrize("ic, constants, pattern", [
    ({"target_index": "idx"}, {"idx": 0},
     r"Loading scalar #0 into scalar #1\n"),
    ({"interpolate": True, "interpolation": {"tolerance": "itol"}},
     {"itol": 1e-6}, r"Tolerance:\s*0\.1000000E-05\n"),
], ids=["target_index", "interpolation_tolerance"])
def test_scalar_field_ic(launcher_script, log_file, tmp_path, previous_run,
                         ic, constants, pattern):
    """Scalar field initial condition keys are resolved from constants."""
    log = run(launcher_script,
              make_case(tmp_path, "field_ic", constants,
                        scalar_field_ic(previous_run, ic)),
              log_file)
    assert re.search(pattern, log), f"{pattern!r} not in the log"


def test_mesh2mesh_tolerance(launcher_script, log_file, tmp_path,
                             previous_run):
    """case.mesh2mesh_tolerance of a restart is resolved from constants."""
    chkp = sorted(previous_run.glob("*.chkp"))
    assert chkp, "the previous run wrote no checkpoint"

    def edit(case):
        case["restart_file"] = str(chkp[-1])
        case["restart_mesh_file"] = MESH
        case["mesh2mesh_tolerance"] = "m2m"
        case["time"]["end_time"] = 0.04

    log = run(launcher_script,
              make_case(tmp_path, "restart", {"m2m": 1e-5}, edit), log_file)
    assert "Restarting from checkpoint" in log
    # The default is 1e-6; the restart interpolation logs the tolerance.
    assert re.search(r"Tolerance:\s*0\.1000000E-04\n", log)


@pytest.fixture(scope="module")
def compressible_neko(tmp_path_factory):
    """A neko executable built with a uniform compressible start."""
    directory = tmp_path_factory.mktemp("compressible")
    shutil.copy(TEST_DIR / "uniform_state.f90", directory)
    result = subprocess.run([get_makeneko(), "uniform_state.f90"],
                            cwd=directory, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    return str(directory / "neko")


def test_compressible(launcher_script, log_file, tmp_path,
                      compressible_neko):
    """gamma, time_order and boundary zone indices of the compressible
    solver are resolved from constants."""
    with open(TEST_DIR / "compressible.case") as f:
        case = json5.load(f)
    case_object = case["case"]
    case_object["mesh_file"] = MESH
    case_object["output_directory"] = str(tmp_path)
    case_object["constants"] = [
        {"name": "gam", "value": 1.3}, {"name": "rk", "value": 3},
        {"name": "walls", "value": [7.0]}]
    fluid = case_object["fluid"]
    fluid["gamma"] = "gam"
    fluid["boundary_conditions"][1]["zone_indices"] = "walls"
    case_object["numerics"]["time_order"] = "rk"
    case_file = str(tmp_path / "compressible.case")
    with open(case_file, "w") as f:
        json.dump(case, f, indent=2)
    log = run(launcher_script, case_file, log_file, neko=compressible_neko)
    assert re.search(r"gamma\s*:\s*1\.300000E\+00", log)
    assert re.search(r"RK order\s*:\s*3\n", log)
