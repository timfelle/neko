import json
from pathlib import Path
import re
import subprocess
import sys

import json5
import pytest


REPO_ROOT = Path(__file__).resolve().parents[4]
VALIDATOR = REPO_ROOT / "contrib" / "validate_case_schema.py"
EXAMPLES_DIR = REPO_ROOT / "examples"
TESTS_DIR = REPO_ROOT / "tests"
CASE_FILES = (
    sorted(EXAMPLES_DIR.rglob("*.case"))
    + sorted(EXAMPLES_DIR.rglob("*.json"))
    + sorted(TESTS_DIR.rglob("*.case"))
)


@pytest.mark.parametrize(
    "case_file",
    CASE_FILES,
    ids=lambda path: path.relative_to(REPO_ROOT).as_posix(),
)
def test_example_case_schema(case_file):
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, (
        f"schema validation failed for {case_file}\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )


@pytest.mark.parametrize(
    "property_path",
    [
        ("case", "load_balance"),
        ("case", "joblimit"),
        ("case", "output_format"),
        ("case", "fluid", "variable_material_properties"),
    ],
)
def test_rejects_known_ignored_properties(tmp_path, property_path):
    with (EXAMPLES_DIR / "tgv" / "tgv.case").open(encoding="utf-8") as handle:
        data = json5.load(handle)

    target = data
    for key in property_path[:-1]:
        target = target[key]
    target[property_path[-1]] = False

    case_file = tmp_path / "invalid.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0


@pytest.mark.parametrize(
    "property_name", ["c_avisc_low", "c_avisc_entropy"]
)
def test_rejects_legacy_artificial_viscosity_parameters(
    tmp_path, property_name
):
    with (EXAMPLES_DIR / "tgv" / "tgv.case").open(encoding="utf-8") as handle:
        data = json5.load(handle)

    data["case"]["numerics"][property_name] = 0.5
    case_file = tmp_path / "invalid.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert property_name in result.stdout


def test_allows_user_extension_properties(tmp_path):
    with (EXAMPLES_DIR / "tgv" / "tgv.case").open(encoding="utf-8") as handle:
        data = json5.load(handle)
    data["case"]["fluid"]["user_parameter"] = 42

    case_file = tmp_path / "custom.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr


def test_allows_custom_scalar_pnpn_boundary_condition(tmp_path):
    with (EXAMPLES_DIR / "scalar_mms" / "scalar_mms.case").open(
        encoding="utf-8"
    ) as handle:
        data = json5.load(handle)

    boundary = data["case"]["scalar"]["boundary_conditions"][0]
    boundary.clear()
    boundary.update(
        {
            "type": "custom_scalar_bc",
            "zone_indices": [1],
            "custom_parameter": 2.0,
        }
    )

    case_file = tmp_path / "custom_scalar_bc.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr


def test_rejects_custom_properties_for_builtin_scalar_bc(tmp_path):
    with (EXAMPLES_DIR / "scalar_mms" / "scalar_mms.case").open(
        encoding="utf-8"
    ) as handle:
        data = json5.load(handle)

    boundary = data["case"]["scalar"]["boundary_conditions"][0]
    boundary["custom_parameter"] = 2.0

    case_file = tmp_path / "invalid_scalar_bc.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0


@pytest.mark.parametrize(
    "example, property_path, value",
    [
        ("tgv/tgv.case",
         ("case", "fluid", "initial_condition", "amplitude"), 1.0),
        ("scalar_mms/scalar_mms.case",
         ("case", "scalar", "initial_condition", "value"), 1.0),
        ("scalar_mms/scalar_mms.case",
         ("case", "scalar", "boundary_conditions", 0),
         {"type": "user_dirichlet", "zone_indices": [1],
          "custom_parameter": 2.0}),
        ("scalar_mms/scalar_mms.case",
         ("case", "scalar", "boundary_conditions", 1),
         {"type": "user_neumann", "zone_indices": [2],
          "custom_parameter": 2.0}),
        ("scalar_mms/scalar_mms.case",
         ("case", "scalar", "boundary_conditions", 0),
         {"type": "user_dirichlet", "zone_indices": [1],
          "value": "my own meaning", "custom_parameter": 2.0}),
        ("scalar_mms/scalar_mms.case",
         ("case", "scalar", "boundary_conditions", 1),
         {"type": "user_neumann", "zone_indices": [2],
          "flux": [1, 2], "custom_parameter": 2.0}),
        ("cylinder/cylinder.case",
         ("case", "fluid", "boundary_conditions", 1),
         {"type": "outflow+user", "zone_indices": [2],
          "value": {"any": "thing"}, "custom_parameter": 2.0}),
        ("cylinder/cylinder.case",
         ("case", "fluid", "boundary_conditions", 0),
         {"type": "user_velocity", "zone_indices": [1],
          "profile": "parabolic"}),
        ("tgv/tgv.case",
         ("case", "fluid", "source_terms"),
         [{"type": "user", "custom_parameter": 2.0}]),
        ("tgv/tgv.case", ("postprocess",), {"tool": "custom"}),
    ],
    ids=[
        "fluid_user_initial_condition",
        "scalar_user_initial_condition",
        "user_dirichlet",
        "user_neumann",
        "user_dirichlet_value",
        "user_neumann_flux",
        "outflow_user_value",
        "user_velocity",
        "user_source_term",
        "top_level_section",
    ],
)
def test_allows_custom_keys_in_user_components(
    tmp_path, example, property_path, value
):
    """User-defined components take parameters the user file defines, and
    other tools may add sections next to the case."""
    with (EXAMPLES_DIR / example).open(encoding="utf-8") as handle:
        data = json5.load(handle)

    target = data
    for key in property_path[:-1]:
        target = target[key]
    target[property_path[-1]] = value

    case_file = tmp_path / "custom.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr


def test_rejects_custom_keys_in_builtin_initial_condition(tmp_path):
    with (EXAMPLES_DIR / "cylinder" / "cylinder.case").open(
        encoding="utf-8"
    ) as handle:
        data = json5.load(handle)

    data["case"]["fluid"]["initial_condition"]["custom_parameter"] = 2.0

    case_file = tmp_path / "invalid_initial_condition.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0


def test_rejects_real_encoding_for_integer_controller_value(tmp_path):
    with (EXAMPLES_DIR / "rayleigh_benard_cylinder" / "rayleigh.case").open(
        encoding="utf-8"
    ) as handle:
        data = json5.load(handle)
    data["case"]["fluid"]["output_value"] = 10.0

    case_file = tmp_path / "invalid.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0


def test_rejects_integer_encoding_for_real_parameter(tmp_path):
    with (EXAMPLES_DIR / "api" / "c" / "cylinder.case").open(
        encoding="utf-8"
    ) as handle:
        data = json5.load(handle)
    data["case"]["fluid"]["Re"] = 200

    case_file = tmp_path / "invalid.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0


def test_rejects_integer_encoding_in_real_array(tmp_path):
    # The first boundary condition is a velocity_value, whose value is real
    with (EXAMPLES_DIR / "cylinder" / "cylinder.case").open(
        encoding="utf-8"
    ) as handle:
        data = json5.load(handle)
    boundary = data["case"]["fluid"]["boundary_conditions"][0]
    boundary["value"] = [1, 0, 0]

    case_file = tmp_path / "invalid.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0


@pytest.mark.parametrize(
    "boundary_type",
    ["outflow", "normal_outflow"],
)
@pytest.mark.parametrize("property_name", ["delta", "velocity_scale"])
def test_rejects_dong_parameters_for_non_dong_outflow(
    tmp_path, boundary_type, property_name
):
    with (EXAMPLES_DIR / "cylinder" / "cylinder.case").open(
        encoding="utf-8"
    ) as handle:
        data = json5.load(handle)

    boundary = data["case"]["fluid"]["boundary_conditions"][1]
    boundary["type"] = boundary_type
    boundary[property_name] = 0.1

    case_file = tmp_path / "invalid.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0


def test_rejects_value_on_plain_outflow(tmp_path):
    """Plain outflow never reads `value`; only the user variant may carry it."""
    with (EXAMPLES_DIR / "cylinder" / "cylinder.case").open(
        encoding="utf-8"
    ) as handle:
        data = json5.load(handle)

    boundary = data["case"]["fluid"]["boundary_conditions"][1]
    boundary["type"] = "outflow"
    boundary["value"] = [1.0, 0.0, 0.0]

    case_file = tmp_path / "invalid.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0


@pytest.mark.parametrize("boundary_type", ["outflow+dong", "normal_outflow+dong"])
def test_allows_dong_outflow_parameters(tmp_path, boundary_type):
    with (EXAMPLES_DIR / "cylinder" / "cylinder.case").open(
        encoding="utf-8"
    ) as handle:
        data = json5.load(handle)

    boundary = data["case"]["fluid"]["boundary_conditions"][1]
    boundary["type"] = boundary_type
    if boundary_type == "normal_outflow+dong":
        boundary["value"] = [1.0, 0.0, 0.0]
    boundary["delta"] = 0.01
    boundary["velocity_scale"] = 1.0

    case_file = tmp_path / "valid.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr


def collect_const_values(value, property_name):
    values = set()
    if isinstance(value, dict):
        property_schema = value.get(property_name)
        if isinstance(property_schema, dict) and "const" in property_schema:
            values.add(property_schema["const"])
        for child in value.values():
            values.update(collect_const_values(child, property_name))
    elif isinstance(value, list):
        for child in value:
            values.update(collect_const_values(child, property_name))
    return values


def factory_case_values(path):
    source = path.read_text(encoding="utf-8")
    return set(re.findall(r"case \([\"']([^\"']+)[\"']\)", source))


def referenced_definition_consts(schema, union_name, property_name):
    definition = schema["$defs"][union_name]
    property_schema = definition.get("properties", {}).get(property_name, {})
    if "const" in property_schema:
        return {property_schema["const"]}

    values = set()
    for variant in definition.get("oneOf", []):
        definition_name = variant["$ref"].removeprefix("#/$defs/")
        values.update(
            referenced_definition_consts(schema, definition_name, property_name)
        )
    return values


def test_all_builtin_simcomps_have_schemas():
    runtime_types = factory_case_values(
        REPO_ROOT / "src/simulation_components/simulation_component_fctry.f90"
    )
    with (REPO_ROOT / "doc/schemas/simulation-components.schema.json").open(
        encoding="utf-8"
    ) as handle:
        schema = json.load(handle)
    refs = json.dumps(schema["$defs"]["simulationComponent"])
    schema_types = set(re.findall(r"urn:neko:schema:simcomps:([^#]+)#", refs))

    assert schema_types == runtime_types


def test_all_builtin_scalar_pnpn_bc_types_have_schemas():
    runtime_types = factory_case_values(
        REPO_ROOT / "src/scalar/scalar_pnpn_bc_fctry.f90"
    )
    with (REPO_ROOT / "doc/schemas/scalar.schema.json").open(
        encoding="utf-8"
    ) as handle:
        schema = json.load(handle)
    schema_types = collect_const_values(
        schema["$defs"]["scalarBoundaryCondition"], "type"
    )
    custom_schema = schema["$defs"]["scalarBoundaryCondition"]["oneOf"][-1]
    excluded_custom_types = set(
        custom_schema["not"]["properties"]["type"]["enum"]
    )

    assert schema_types == runtime_types == excluded_custom_types


@pytest.mark.parametrize(
    ("factory", "schema", "property_name", "extra_runtime_types"),
    [
        (
            "src/wall_models/wall_model_fctry.f90",
            "doc/schemas/fluid-pnpn-boundary-conditions.schema.json",
            "model",
            set(),
        ),
        (
            "src/source_terms/source_term_fctry.f90",
            "doc/schemas/source-terms.schema.json",
            "type",
            {"user"},
        ),
        (
            "src/mesh/point_zone_fctry.f90",
            "doc/schemas/point-zones.schema.json",
            "geometry",
            {"combine"},
        ),
    ],
)
def test_all_builtin_factory_types_have_schemas(
    factory, schema, property_name, extra_runtime_types
):
    runtime_types = factory_case_values(REPO_ROOT / factory) | extra_runtime_types
    with (REPO_ROOT / schema).open(encoding="utf-8") as handle:
        schema_data = json.load(handle)
    if schema == "doc/schemas/source-terms.schema.json":
        schema_types = referenced_definition_consts(
            schema_data, "fluidSourceTerm", property_name
        )
    else:
        schema_types = collect_const_values(schema_data, property_name)

    assert schema_types == runtime_types


def _validate_modified(tmp_path, example, property_path, value):
    """Validate EXAMPLES_DIR/example with the entry at property_path
    replaced by value (a missing last key is added)."""
    with (EXAMPLES_DIR / example).open(encoding="utf-8") as handle:
        data = json5.load(handle)

    target = data
    for key in property_path[:-1]:
        target = target[key]
    target[property_path[-1]] = value

    case_file = tmp_path / "modified.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    return subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )


# Keys that Neko reads without resolving constants. A schema that accepts a
# constants reference there lets a case through that the solver then rejects.
# Each entry is (id, example, path, literal value, constants reference).
LITERAL_ONLY_KEYS = [
    ("schwarz_iterations", "tgv/tgv.case",
     ("case", "fluid", "schwarz_iterations"), 1, "n_iter"),
    ("mesh2mesh_tolerance", "tgv/tgv.case",
     ("case", "mesh2mesh_tolerance"), 1e-6, "tol"),
    ("compressible_gamma", "euler_tgv/euler_tgv.case",
     ("case", "fluid", "gamma"), 1.4, "gamma"),
    ("compressible_zone_indices", "euler_tgv/euler_tgv.case",
     ("case", "fluid", "boundary_conditions"),
     [{"type": "no_slip", "zone_indices": [1]}],
     [{"type": "no_slip", "zone_indices": "zones"}]),
    ("hsmg_coarse_grid_iterations", "cylinder/cylinder.case",
     ("case", "fluid", "pressure_solver", "preconditioner"),
     {"type": "hsmg", "coarse_grid": {"iterations": 5}},
     {"type": "hsmg", "coarse_grid": {"iterations": "n_iter"}}),
    ("dong_delta", "cylinder/cylinder.case",
     ("case", "fluid", "boundary_conditions", 1),
     {"type": "outflow+dong", "zone_indices": [2], "delta": 0.01},
     {"type": "outflow+dong", "zone_indices": [2], "delta": "delta"}),
    ("source_term_start_time", "tgv/tgv.case",
     ("case", "fluid", "source_terms"),
     [{"type": "constant", "values": [1.0, 0.0, 0.0], "start_time": 1.0}],
     [{"type": "constant", "values": [1.0, 0.0, 0.0],
       "start_time": "t_on"}]),
    ("gradient_jump_penalty_tau", "tgv/tgv.case",
     ("case", "fluid", "source_terms"),
     [{"type": "gradient_jump_penalty", "tau": 0.02}],
     [{"type": "gradient_jump_penalty", "tau": "tau"}]),
    ("scalar_bc_zone_indices", "scalar_mms/scalar_mms.case",
     ("case", "scalar", "boundary_conditions", 0),
     {"type": "dirichlet", "zone_indices": [1], "value": 1.0},
     {"type": "dirichlet", "zone_indices": "zones", "value": 1.0}),
]


@pytest.mark.parametrize(
    "example, property_path, good, bad",
    [entry[1:] for entry in LITERAL_ONLY_KEYS],
    ids=[entry[0] for entry in LITERAL_ONLY_KEYS],
)
def test_literal_only_keys_reject_constant_references(
    tmp_path, example, property_path, good, bad
):
    accepted = _validate_modified(tmp_path, example, property_path, good)
    assert accepted.returncode == 0, accepted.stdout + accepted.stderr

    rejected = _validate_modified(tmp_path, example, property_path, bad)
    assert rejected.returncode != 0


def test_allows_split_velocity_pressure_boundary_condition(tmp_path):
    path = ("case", "fluid", "boundary_conditions", 1)
    split = {
        "zone_indices": [2],
        "velocity": {"type": "non_normal", "value": [1.0, 0.0, 0.0]},
        "pressure": {"type": "expression", "value": "sin(t)"},
    }
    result = _validate_modified(tmp_path, "cylinder/cylinder.case", path, split)
    assert result.returncode == 0, result.stdout + result.stderr

    # A joint type cannot be combined with the split form.
    split["type"] = "outflow"
    result = _validate_modified(tmp_path, "cylinder/cylinder.case", path, split)
    assert result.returncode != 0


def test_allows_compressible_slip_boundary_condition(tmp_path):
    result = _validate_modified(
        tmp_path, "euler_tgv/euler_tgv.case",
        ("case", "fluid", "boundary_conditions"),
        [{"type": "slip", "zone_indices": [1]}])
    assert result.returncode == 0, result.stdout + result.stderr


def test_requires_dealias_for_pnpn(tmp_path):
    with (EXAMPLES_DIR / "tgv" / "tgv.case").open(encoding="utf-8") as handle:
        data = json5.load(handle)
    del data["case"]["numerics"]["dealias"]
    case_file = tmp_path / "no_dealias.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True, text=True)
    assert result.returncode != 0


def test_requires_output_checkpoints(tmp_path):
    with (EXAMPLES_DIR / "tgv" / "tgv.case").open(encoding="utf-8") as handle:
        data = json5.load(handle)
    del data["case"]["output_checkpoints"]
    case_file = tmp_path / "no_checkpoints.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True, text=True)
    assert result.returncode != 0


def test_time_order_four_only_for_compressible(tmp_path):
    path = ("case", "numerics", "time_order")
    result = _validate_modified(tmp_path, "tgv/tgv.case", path, 4)
    assert result.returncode != 0

    result = _validate_modified(tmp_path, "euler_tgv/euler_tgv.case", path, 4)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("points", [None, []], ids=["absent", "empty"])
def test_allows_probes_without_points(tmp_path, points):
    # Neko sets up zero probes when no points are given.
    path = ("case", "simulation_components")
    probes = {"type": "probes", "fields": ["u", "p"],
              "output_file": "probes.csv"}
    if points is not None:
        probes["points"] = points
    result = _validate_modified(tmp_path, "cylinder/cylinder.case", path,
                                [probes])
    assert result.returncode == 0, result.stdout + result.stderr

    # The probes schema still applies: a malformed points list is rejected.
    probes["points"] = "points.csv"
    result = _validate_modified(tmp_path, "cylinder/cylinder.case", path,
                                [probes])
    assert result.returncode != 0


_DELETE = object()


def _validate_edits(tmp_path, example, edits):
    """Validate EXAMPLES_DIR/example after applying edits, a list of
    (property_path, value) pairs. The value _DELETE removes the entry."""
    with (EXAMPLES_DIR / example).open(encoding="utf-8") as handle:
        data = json5.load(handle)

    for property_path, value in edits:
        target = data
        for key in property_path[:-1]:
            target = target[key]
        if value is _DELETE:
            del target[property_path[-1]]
        else:
            target[property_path[-1]] = value

    case_file = tmp_path / "edited.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    return subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )


def _assert_accepted(result):
    assert result.returncode == 0, result.stdout + result.stderr


def _assert_rejected(result):
    assert result.returncode != 0, result.stdout + result.stderr


# Keys and values the code accepts. The schemas must not take away the
# freedom to use them, e.g. for a user-defined component that reads them.
# Each entry is (id, example, edits).
ACCEPTED_BY_THE_CODE = [
    ("fluid_relative_tolerance", "cylinder/cylinder.case",
     [(("case", "fluid", "velocity_solver", "relative_tolerance"), 1e-3)]),
    ("scalar_relative_tolerance", "scalar_mms/scalar_mms.case",
     [(("case", "scalar", "solver", "relative_tolerance"), 1e-3)]),
    # Re is only rejected together with mu or rho when the case does not
    # provide material properties in the user file.
    ("re_with_mu_and_rho", "cylinder/cylinder.case",
     [(("case", "fluid", "mu"), 0.01), (("case", "fluid", "rho"), 1.0)]),
    ("scalar_pe_with_lambda", "scalar_mms/scalar_mms.case",
     [(("case", "scalar", "lambda"), 1.0), (("case", "scalar", "Pe"), 1.0)]),
    ("svv_without_enabled", "cylinder/cylinder.case",
     [(("case", "fluid", "svv"), {"direction": "rst"})]),
    ("no_version", "cylinder/cylinder.case", [(("version",), _DELETE)]),
    ("flow_rate_direction_zero", "cylinder/cylinder.case",
     [(("case", "fluid", "flow_rate_force"),
       {"direction": 0, "value": 1.0, "use_averaged_flow": True})]),
    ("scalar_neumann_flux_one_entry", "scalar_mms/scalar_mms.case",
     [(("case", "scalar", "boundary_conditions", 1),
       {"type": "neumann", "zone_indices": [2], "flux": [1.0]})]),
    ("empty_simcomp_comment", "cylinder/cylinder.case",
     [(("case", "simulation_components"),
       [{"type": "lambda2", "comment": ""}])]),
    ("simcomp_order_unset", "cylinder/cylinder.case",
     [(("case", "simulation_components"),
       [{"type": "lambda2", "order": -1}])]),
    ("hdf5_field_output", "cylinder/cylinder.case",
     [(("case", "simulation_components"),
       [{"type": "curl", "fields": ["u", "v", "w"],
         "output_filename": "curl", "output_format": "hdf5"}])]),
    # Only hsmg has a floor (pc_hsmg.f90); pnpn itself runs at order 2 and 1.
    ("pnpn_order_two_with_hsmg", "cylinder/cylinder.case",
     [(("case", "numerics", "polynomial_order"), 2),
      (("case", "fluid", "pressure_solver", "preconditioner"),
       {"type": "hsmg"})]),
    ("pnpn_order_one_with_jacobi", "cylinder/cylinder.case",
     [(("case", "numerics", "polynomial_order"), 1),
      (("case", "fluid", "pressure_solver", "preconditioner"),
       {"type": "jacobi"}),
      (("case", "fluid", "velocity_solver", "preconditioner"),
       {"type": "jacobi"})]),
]


@pytest.mark.parametrize(
    "example, edits",
    [entry[1:] for entry in ACCEPTED_BY_THE_CODE],
    ids=[entry[0] for entry in ACCEPTED_BY_THE_CODE],
)
def test_allows_what_the_code_accepts(tmp_path, example, edits):
    _assert_accepted(_validate_edits(tmp_path, example, edits))


def test_scalar_initial_condition_not_required_on_restart(tmp_path):
    restart = (("case", "restart_file"), "restart.chkp")
    no_ic = (("case", "scalar", "initial_condition"), _DELETE)
    _assert_rejected(
        _validate_edits(tmp_path, "scalar_mms/scalar_mms.case", [no_ic]))
    _assert_accepted(
        _validate_edits(tmp_path, "scalar_mms/scalar_mms.case",
                        [restart, no_ic]))


def test_field_subsampler_point_zone_needs_no_output_filename(tmp_path):
    # The subsampler removes point_zone before it sets up its writer and
    # names the output after itself, see the example in the simcomps page.
    subsampler = {"type": "field_subsampler", "name": "sub",
                  "source_fields": ["u"], "point_zone": "box",
                  "polynomial_order": 3}
    zone = {"name": "box", "geometry": "box", "x_bounds": [0.0, 1.0],
            "y_bounds": [0.0, 1.0], "z_bounds": [0.0, 1.0]}
    edits = [(("case", "simulation_components"), [subsampler]),
             (("case", "point_zones"), [zone])]
    _assert_accepted(_validate_edits(tmp_path, "cylinder/cylinder.case",
                                     edits))


def test_field_writer_point_zone_needs_fld_and_output_filename(tmp_path):
    zone = (("case", "point_zones"),
            [{"name": "box", "geometry": "box", "x_bounds": [0.0, 1.0],
              "y_bounds": [0.0, 1.0], "z_bounds": [0.0, 1.0]}])
    writer = {"type": "curl", "fields": ["u", "v", "w"],
              "point_zone": "box"}

    def check(**extra):
        return _validate_edits(
            tmp_path, "cylinder/cylinder.case",
            [zone, (("case", "simulation_components"), [{**writer, **extra}])])

    _assert_rejected(check())
    _assert_accepted(check(output_filename="curl"))
    _assert_accepted(check(output_filename="curl", output_format="fld"))
    _assert_rejected(check(output_filename="curl", output_format="vtkhdf"))
    _assert_rejected(check(output_filename="curl", output_format="adios2"))


# Rejected because the code aborts on them. Each entry is
# (id, example, edits).
REJECTED_AS_THE_CODE_ABORTS = [
    ("order_one_with_hsmg_pressure", "cylinder/cylinder.case",
     [(("case", "numerics", "polynomial_order"), 1),
      (("case", "fluid", "pressure_solver", "preconditioner"),
       {"type": "hsmg"})]),
    ("shear_stress_without_full_stress", "cylinder/cylinder.case",
     [(("case", "fluid", "boundary_conditions", 4),
       {"type": "shear_stress", "zone_indices": [7],
        "value": [0.0, 0.0, 0.0]})]),
    ("split_wall_model_without_full_stress", "cylinder/cylinder.case",
     [(("case", "fluid", "boundary_conditions", 4),
       {"zone_indices": [7], "velocity": {"type": "wall_model"},
        "pressure": {"type": "zero_dirichlet"}})]),
    ("org_output_without_nsamples", "cylinder/cylinder.case",
     [(("case", "fluid", "output_control"), "org")]),
    ("is_user_without_type", "cylinder/cylinder.case",
     [(("case", "simulation_components"), [{"is_user": True}])]),
    ("simcomp_order_zero", "cylinder/cylinder.case",
     [(("case", "simulation_components"),
       [{"type": "lambda2", "order": 0}])]),
    ("probes_output_file_suffix", "cylinder/cylinder.case",
     [(("case", "simulation_components"),
       [{"type": "probes", "fields": ["u"], "output_file": "probes.dat"}])]),
    ("probes_line_with_one_point", "cylinder/cylinder.case",
     [(("case", "simulation_components"),
       [{"type": "probes", "fields": ["u"], "output_file": "probes.csv",
         "points": [{"type": "line", "start": [0.0, 0.0, 0.0],
                     "end": [1.0, 0.0, 0.0], "amount": 1}]}])]),
    ("probes_circle_radius_zero", "cylinder/cylinder.case",
     [(("case", "simulation_components"),
       [{"type": "probes", "fields": ["u"], "output_file": "probes.csv",
         "points": [{"type": "circle", "center": [0.0, 0.0, 0.0],
                     "normal": [0.0, 0.0, 1.0], "radius": 0.0,
                     "amount": 4, "axis": "x"}]}])]),
    ("wall_distance_integer", "cylinder/cylinder.case",
     [(("case", "fluid", "full_stress_formulation"), True),
      (("case", "fluid", "boundary_conditions", 4),
       {"type": "wall_model", "model": "spalding", "zone_indices": [7],
        "kappa": 0.41, "B": 5.2,
        "sampling": {"type": "distance", "value": 1}})]),
    ("brinkman_pde_filter_without_radius", "immersed_zones/sphere.case",
     [(("case", "fluid", "source_terms", 0, "filter"), {"type": "PDE"})]),
    ("brinkman_file_index_constant", "immersed_zones/sphere.case",
     [(("case", "fluid", "source_terms", 0, "objects", 0),
       {"type": "file", "file_name": "indicator0.f00000",
        "file_index": "idx"})]),
    ("scalar_source_term_of_fluid_type", "scalar_mms/scalar_mms.case",
     [(("case", "scalar", "source_terms"),
       [{"type": "boussinesq", "g": [0.0, 0.0, -1.0],
         "reference_value": 1.0}])]),
]


@pytest.mark.parametrize(
    "example, edits",
    [entry[1:] for entry in REJECTED_AS_THE_CODE_ABORTS],
    ids=[entry[0] for entry in REJECTED_AS_THE_CODE_ABORTS],
)
def test_rejects_what_the_code_aborts_on(tmp_path, example, edits):
    _assert_rejected(_validate_edits(tmp_path, example, edits))


def test_wall_distance_real_is_accepted(tmp_path):
    # The control of wall_distance_integer: a real passes
    edits = [(("case", "fluid", "full_stress_formulation"), True),
             (("case", "fluid", "boundary_conditions", 4),
              {"type": "wall_model", "model": "spalding",
               "zone_indices": [7], "kappa": 0.41, "B": 5.2,
               "sampling": {"type": "distance", "value": 1.0}})]
    _assert_accepted(_validate_edits(tmp_path, "cylinder/cylinder.case",
                                     edits))


def test_wall_zone_indices_reject_constant_reference(tmp_path):
    path = ("case", "fluid", "wall_zone_indices")
    _assert_accepted(_validate_edits(tmp_path, "cylinder/cylinder.case",
                                     [(path, [7])]))
    _assert_rejected(_validate_edits(tmp_path, "cylinder/cylinder.case",
                                     [(path, "walls")]))


def test_shear_stress_accepted_with_full_stress(tmp_path):
    edits = [(("case", "fluid", "full_stress_formulation"), True),
             (("case", "fluid", "boundary_conditions", 4),
              {"type": "shear_stress", "zone_indices": [7],
               "value": [0.0, 0.0, 0.0]})]
    _assert_accepted(_validate_edits(tmp_path, "cylinder/cylinder.case",
                                     edits))


def test_brinkman_pde_filter_with_radius_is_accepted(tmp_path):
    edits = [(("case", "fluid", "source_terms", 0, "filter"),
              {"type": "PDE", "radius": 0.1})]
    _assert_accepted(_validate_edits(tmp_path, "immersed_zones/sphere.case",
                                     edits))


# Types registered through Neko's extension mechanism define keys of their
# own, so the schemas accept any type name that is not a built-in one.
# Each entry is (id, example, edits) of a valid case with such a type.
REGISTERED_TYPES = [
    ("simcomp", "cylinder/cylinder.case",
     [(("case", "simulation_components"),
       [{"type": "my_simcomp", "my_key": [1, 2], "compute_control": "tsteps",
         "compute_value": 10}])]),
    ("fluid_source_term", "cylinder/cylinder.case",
     [(("case", "fluid", "source_terms"),
       [{"type": "my_source", "my_key": "value"}])]),
    ("scalar_source_term", "scalar_mms/scalar_mms.case",
     [(("case", "scalar", "source_terms"),
       [{"type": "my_source", "my_key": "value"}])]),
    ("point_zone", "cylinder/cylinder.case",
     [(("case", "point_zones"),
       [{"name": "mine", "geometry": "my_geometry", "my_key": 1.0}])]),
    ("point_zone_in_combine", "cylinder/cylinder.case",
     [(("case", "point_zones"),
       [{"name": "mine", "geometry": "my_geometry", "my_key": 1.0},
        {"name": "both", "geometry": "combine", "operator": "OR",
         "subsets": [{"name": "mine"}]}])]),
    ("les_model", "cylinder/cylinder.case",
     [(("case", "simulation_components"),
       [{"type": "les_model", "model": "my_model", "my_key": 1.0}])]),
    ("avm_model", "euler_tgv/euler_tgv.case",
     [(("case", "simulation_components"),
       [{"type": "artificial_viscosity_model", "model": "my_model",
         "my_key": 1.0}])]),
    ("wall_model", "cylinder/cylinder.case",
     [(("case", "fluid", "full_stress_formulation"), True),
      (("case", "fluid", "boundary_conditions", 4),
       {"type": "wall_model", "model": "my_model", "zone_indices": [7],
        "h_index": 1, "my_key": 1.0})]),
    ("preconditioner", "cylinder/cylinder.case",
     [(("case", "fluid", "velocity_solver", "preconditioner"),
       {"type": "my_precon", "my_key": 1.0})]),
]


@pytest.mark.parametrize(
    "example, edits",
    [entry[1:] for entry in REGISTERED_TYPES],
    ids=[entry[0] for entry in REGISTERED_TYPES],
)
def test_allows_registered_types(tmp_path, example, edits):
    _assert_accepted(_validate_edits(tmp_path, example, edits))


# The built-in types are still checked strictly: a key they do not define
# or one with the wrong type is reported, rather than passed on as a custom
# key. (A misspelled type name cannot be told from a registered one.)
BUILTIN_TYPES_STAY_STRICT = [
    ("simcomp_invalid_control_value", "cylinder/cylinder.case",
     [(("case", "simulation_components"),
       [{"type": "lambda2", "compute_control": "tsteps",
         "compute_value": 0}])]),
    ("source_term_unknown_key", "cylinder/cylinder.case",
     [(("case", "fluid", "source_terms"),
       [{"type": "constant", "values": [1.0, 0.0, 0.0], "my_key": 1.0}])]),
    ("point_zone_unknown_key", "cylinder/cylinder.case",
     [(("case", "point_zones"),
       [{"name": "box", "geometry": "box", "x_bounds": [0.0, 1.0],
         "y_bounds": [0.0, 1.0], "z_bounds": [0.0, 1.0], "my_key": 1.0}])]),
    ("preconditioner_unknown_key", "cylinder/cylinder.case",
     [(("case", "fluid", "velocity_solver", "preconditioner"),
       {"type": "jacobi", "my_key": 1.0})]),
    ("les_model_missing_key", "cylinder/cylinder.case",
     [(("case", "simulation_components"),
       [{"type": "les_model", "model": "deardorff"}])]),
    ("wall_model_missing_key", "cylinder/cylinder.case",
     [(("case", "fluid", "full_stress_formulation"), True),
      (("case", "fluid", "boundary_conditions", 4),
       {"type": "wall_model", "model": "spalding", "zone_indices": [7],
        "h_index": 1})]),
    ("avm_wrong_type", "euler_tgv/euler_tgv.case",
     [(("case", "simulation_components"),
       [{"type": "artificial_viscosity_model",
         "model": "entropy_viscosity", "c_avisc_low": "c"}])]),
]


@pytest.mark.parametrize(
    "example, edits",
    [entry[1:] for entry in BUILTIN_TYPES_STAY_STRICT],
    ids=[entry[0] for entry in BUILTIN_TYPES_STAY_STRICT],
)
def test_builtin_types_stay_strict(tmp_path, example, edits):
    _assert_rejected(_validate_edits(tmp_path, example, edits))


def _case_values(path):
    source = (REPO_ROOT / path).read_text(encoding="utf-8")
    return set(re.findall(r"case\s*\(\s*[\"']([^\"']+)[\"']\s*\)", source))


def _excluded_types(node, property_name):
    """The enum of the `not` that keeps the built-in types out of a
    fallback branch for registered types, searched in node."""
    found = set()
    if isinstance(node, dict):
        negated = node.get("not")
        if isinstance(negated, dict):
            enum = negated.get("properties", {}).get(property_name, {}).get(
                "enum")
            if enum:
                found.update(enum)
        for child in node.values():
            found.update(_excluded_types(child, property_name))
    elif isinstance(node, list):
        for child in node:
            found.update(_excluded_types(child, property_name))
    return found


@pytest.mark.parametrize(
    ("factory", "schema", "defs", "property_name", "extra_types"),
    [
        ("src/simulation_components/simulation_component_fctry.f90",
         "doc/schemas/simulation-components.schema.json",
         ["simulationComponent"], "type", set()),
        ("src/source_terms/source_term_fctry.f90",
         "doc/schemas/source-terms.schema.json",
         ["customSourceTerm"], "type", {"user"}),
        ("src/mesh/point_zone_fctry.f90",
         "doc/schemas/point-zones.schema.json",
         ["pointZone"], "geometry", {"combine"}),
        ("src/les/les_model_fctry.f90",
         "doc/schemas/simcomps/les_model.schema.json",
         ["simulationComponent"], "model", set()),
        ("src/viscous_regularization/avm/avm_fctry.f90",
         "doc/schemas/simcomps/artificial_viscosity_model.schema.json",
         ["simulationComponent"], "model", set()),
        ("src/wall_models/wall_model_fctry.f90",
         "doc/schemas/fluid-pnpn-boundary-conditions.schema.json",
         ["fluidPnpnBoundaryCondition"], "model", set()),
        ("src/krylov/precon_fctry.f90",
         "doc/schemas/common.schema.json",
         ["customPreconditioner"], "type", set()),
    ],
    ids=["simcomp", "source_term", "point_zone", "les_model", "avm_model",
         "wall_model", "preconditioner"],
)
def test_registered_type_fallback_excludes_exactly_the_builtins(
    factory, schema, defs, property_name, extra_types
):
    with (REPO_ROOT / schema).open(encoding="utf-8") as handle:
        schema_data = json.load(handle)
    excluded = set()
    for name in defs:
        excluded |= _excluded_types(schema_data["$defs"][name], property_name)

    assert excluded == _case_values(factory) | extra_types
