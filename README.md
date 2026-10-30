# APC Propeller Tools

Analyze APC propeller public geometry and performance `.dat` files. Fit thrust and power coefficients as a function of advance ratio and Reynolds number, then export a single JSON file per propeller for use in other tools.

## Data files

Point the analyzer at a folder that uses this naming scheme:

| File | Contents |
|------|----------|
| `<name>_performance.dat` | APC CT/CP vs J tables by RPM |
| `<name>_geometry.dat` | APC airfoil / inertia summary |

Example: `5x11E_performance.dat` + `5x11E_geometry.dat`

## Usage (MATLAB command window)

```matlab
% Prompt for folder + selection
analyze_apc_prop

% Folder with one prop -> runs automatically
analyze_apc_prop('/path/to/prop_data')

% Folder with multiple props -> lists them and asks you to choose
analyze_apc_prop('/path/to/prop_data')

% Analyze specific props
analyze_apc_prop('/path/to/prop_data', '5x11E')
analyze_apc_prop('/path/to/prop_data', '5x11E,7x11E')
analyze_apc_prop('/path/to/prop_data', {'5x11E','7x15E'})

% Analyze every prop in the folder
analyze_apc_prop('/path/to/prop_data', 'all')
```

Each run writes:

```text
<dataFolder>/<name>_analysed.json
```

## Fit model

Bivariate polynomial response surface:

```text
Re_hat = Re / 1e5
C(J, Re_hat) = sum_i sum_j a_ij * J^i * Re_hat^j
```

Default degrees: `deg_J = 4`, `deg_Re = 2` (editable near the top of `analyze_apc_prop.m`).

Coefficient ordering is stored explicitly in the JSON (`fit_model.basis` + `coefficient_ordering`).

## JSON contents (high level)

- Prop ID, manufacturer, model, diameter, pitch
- Source filenames/paths and fit timestamp
- Fit type, CT/CP coefficients, basis ordering
- Valid J and Re ranges, units/sign conventions, schema version
- Chord min/max/average and chord at 0.75R, plus mass/MoI when geometry is present
