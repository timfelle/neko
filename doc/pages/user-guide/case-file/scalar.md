# Scalar {#case-file_scalar}
The scalar object allows to add a scalar transport equation to the solution. The
solution variable is called `s` by default, but can be controlled by the
`name` entry in the case file. When several scalars are defined in the `scalars`
array, the default names are `s_1`, `s_2`, and so on. In the fld files, it is
saved as `temperature`. Some properties of the object are inherited from `fluid`: the
value of the density, and the output control.

## Material properties

The scalar equation requires defining additional material properties: the
specific heat capacity and thermal conductivity. These are provided as `cp` and
`lambda`. Similarly to the fluid, one can provide the Peclet number, `Pe`, as an
alternative. In this case, `cp` is set to 1 and `lambda` to the inverse of `Pe`.

Different from the setup in the fluid, turbulence modelling is enabled by setting the `alphat` json entry.

## Turbulence modelling

The user could choose to either relate the eddy diffusivity field to the eddy viscosity
field in the Fluid, or model the eddy diffusivity field by some particular SGS models,
by setting up the `nut_dependency` entry.
If the eddy diffusivity field is associated to the eddy viscosity field by a coefficient
`Pr_t`, the eddy viscosity values will be divided by it to produce eddy diffusivity.
And the corresponding setting could be done by the following:

```json
{
  "alphat": {
    "nut_dependency": true,
    "nut_field": "nut",
    "Pr_t": 0.7
  }
}
```

Otherwise one could have some SGS models providing an eddy diffusivity field, and the
user could set it up by the following manner to include an eddy diffusivity field called
`temperature_alphat`:

```json
{
  "alphat": {
    "nut_dependency": false,
    "alphat_field": "temperature_alphat"
  }
}
```

## Boundary conditions {#case-file_scalar-boundary-conditions}

The boundary conditions for the scalar are specified through the
`boundary_conditions` keyword, which follows the same format as the fluid, for
specifying the type of the condition and where it is applied.
The following types of conditions are available for the scalar:

* `dirichlet`. Sets the value of the scalar, controlled by the `value` keyword.
  ```json
  {
    "type": "dirichlet",
    "value": 1,
    "zone_indices": [1, 2]
  }
  ```
* `expression_dirichlet`. Sets the value of the scalar from a mathematical
  expression, given as a string under the `value` keyword. See
  [expressions](@ref case-file_expressions) for the syntax.
  ```json
  {
    "type": "expression_dirichlet",
    "value": "T_wall + dT*tanh((x - x0)/L)",
    "zone_indices": [1, 2]
  }
  ```
* `neumann`. Sets the flux of the scalar, controlled by the `flux` keyword.
  ```json
  {
    "type": "neumann",
    "flux": 1,
    "zone_indices": [1, 2]
  }
  ```
* `user_dirichlet`. User-updated Dirichlet boundary condition, see the
  [user-file documentation](@ref user-file_field-dirichlet-update).
  Custom keys of your choosing, including `value` or `flux`, are accepted,
  see [user-defined condition keys](@ref user-file_user-bc-keys).
  ```json
  {
    "type": "user_dirichlet",
    "zone_indices": [1, 2]
  }
  ```
* `user_neumann`. User-updated Neumann boundary condition, see the
  [user-file documentation](@ref user-file_field-neumann-update).
  Custom keys of your choosing, including `value` or `flux`, are accepted,
  see [user-defined condition keys](@ref user-file_user-bc-keys).
  ```json
  {
    "type": "user_neumann",
    "zone_indices": [1, 2]
  }
  ```
* `user_neumann`. User-defined scalar flux, configured through the
  [user-file Neumann callback](@ref user-file_field-neumann-update).
  ```json
  {
    "type": "user_neumann",
    "zone_indices": [1, 2]
  }
  ```
* `overset_interface`. Retrieves scalar Dirichlet data from another concurrent
  Neko simulation on an overlapping domain. It accepts `name`, `order`,
  `relaxation`, `interpolation.tolerance`, `interpolation.padding`, and `log`
  with the same meanings and defaults as the
  [fluid overset interface](@ref case-file_overset-interface). The first IEXT
  update of each physical timestep is unrelaxed; relaxation is applied only if
  the boundary is updated again at the same timestep.
  ```json
  {
    "type": "overset_interface",
    "zone_indices": [1, 2],
    "order": 3,
    "relaxation": 0.7,
    "interpolation": {
      "tolerance": 1.0e-8,
      "padding": 0.01
    },
    "log": false
  }
  ```

## Initial conditions

The object `initial_condition` is used to provide initial conditions.
It is mandatory.
The means of prescribing the values are controlled via the `type` keyword:

1. `user`, the values are set inside the compiled user file as explained in the
[user defined initial condition](@ref user-file_user-ic) section of the user
file documentation.
2. `uniform`, the value is a constant scalar, looked up under the `value`
   keyword.
3. `expression`, the `value` keyword holds a single string, a mathematical
   expression prescribing the initial value. See
   [expressions](@ref case-file_expressions) for the syntax. For example
   ~~~~~~~~~~~~~~~{.json}
   {
     "type": "expression",
     "value": "T_hot + dT*tanh((z - z0)/L)"
   }
   ~~~~~~~~~~~~~~~
4. `point_zone`, the values are set to a constant base value, supplied under the
   `base_value` keyword, and then assigned the value given under the
   `zone_value` keyword inside a point zone. The point zone is specified by the
   `zone_name` keyword, and should be defined in the `case.point_zones` object.
   See more about [point zones](@ref point-zones).
5. `field`, where the initial condition is retrieved from a field file. Works
   in the same way as for the fluid. See the
   [fluid section](@ref case-file_fluid-ic) for detailed explanations. In
   addition to the keywords described there, `target_index` selects which
   field in the file is read. It is a non-negative integer, where 0 is the
   temperature. It defaults to the position of the scalar among the enabled
   scalars not named `temperature`, counted from 1; a scalar named
   `temperature` defaults to 0. Scalars with `enabled` set to `false` are not
   counted. When initialising from a file written with a different set of
   enabled scalars, set `target_index` explicitly, otherwise the fields are
   read shifted. The keywords `interpolation.tolerance` and `interpolation.padding` must be given as
   literal reals.

## Source terms

The configuration of source terms is the same as for the fluid. A demonstration
of using source terms for the scalar can be found in the `scalar_mms` example.

## Linear solver configuration

Should be provided as an object under the `solver` keyword. For available
configuration options, see the corresponding documentation for the fliud. A
standard choice would be `"type": "cg"` and `"preconditioner": "jacobi"`.

## Full parameter table

| Name                           | Description                                                           | Admissible values                           | Default value |
| ------------------------------ | --------------------------------------------------------------------- | ------------------------------------------- | ------------- |
| `enabled`                      | Whether to enable the scalar computation. If `false`, `initial_condition` and `solver` are not required. For an entry of `scalars`, the entry is skipped, and the other entries keep their position-based default names `s_<i>`. | `true` or `false`                           | `true`        |
| `scheme`                       | The scalar solver scheme.                                             | `pnpn`                                      | `pnpn`        |
| `name`                         | The name associated to the scalar solver, which is also the name of the solution in the field registry. | A string                                    | `s` for a single scalar, `s_<i>` for entry `i` of `scalars` |
| `freeze`                       | Whether to skip the time step of the scalar equation, leaving the field unchanged. | `true` or `false`                           | `false`       |
| `Pe`                           | The Peclet number.                                                    | Positive real                               | -             |
| `cp`                           | Specific heat capacity.                                               | Positive real                               | -             |
| `lambda`                       | Thermal conductivity.                                                 | Positive real                               | -             |
| `alphat.nut_dependency`        | Whether the eddy diffusivity depends on the eddy kinematic viscosity. | `true` or `false`                           | -             |
| `alphat.alphat_field`          | Name of the turbulent diffusivity field. Required if `nut_dependency` is `false`. | String                          | -             |
| `alphat.nut_field`             | Name of the turbulent kinematic viscosity field. Required if `nut_dependency` is `true`. | String                      | -             |
| `alphat.Pr_t`                  | Turbulent Prandtl number. Required if `nut_dependency` is `true`.     | Positive real                               | -             |
| `boundary_conditions`          | Array of JSON objects, defining the boundary conditions.              | See list of boundary conditions above       | -             |
| `initial_condition.type`       | Initial condition type.                                               | `user`, `uniform`, `expression`, `point_zone`, `field` | -           |
| `initial_condition.value`      | Value of the scalar initial condition.                                | Real, or a string if `"type" = "expression"` | -             |
| `source_terms`                 | Array of JSON objects, defining additional source terms.              | See list of source terms above              | -             |
| `advection`                    | Whether to compute the advetion term.                                 | `true` or `false`                           | `true`        |
| `solver.type`                  | Linear solver for scalar equation.                                    | `cg`, `pipecg`, `bicgstab`, `cacg`, `gmres` | -             |
| `solver.preconditioner.type`   | Linear solver preconditioner for the momentum equation.               | `ident`, `hsmg`, `jacobi`                   | -             |
| `solver.absolute_tolerance`    | Linear solver convergence criterion for the momentum equation.        | Positive real                               | -             |
| `solver.max_iterations`        | Linear solver max iteration count for the scalar equation.            | Positive integer                            | 1000          |
| `solver.space_size`            | Krylov space size of `gmres` for the scalar equation.                 | Positive integer                            | 30            |
| `solver.projection_space_size` | Projection space size for the scalar equation.                        | Positive integer                            | 0             |
| `solver.projection_hold_steps` | Holding steps of the projection for the scalar equation.              | Positive integer                            | 5             |
