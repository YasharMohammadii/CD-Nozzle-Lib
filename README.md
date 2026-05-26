# CD-Nozzle-Lib

**1D converging-diverging nozzle flow simulator** – educational & preliminary design tool.  
Built for the *Advanced Gas Dynamics* course, this library computes flow properties (Mach, pressure, temperature, density) along a CD nozzle given geometry, gas properties, and boundary conditions.

## Features

- 5 modular classes:
  - `Constants` – stores input parameters
  - `NozzleGeometry` – computes radius & area along nozzle axis
  - `ThermoRelations` – thermodynamic relations for 1D isentropic / Fanno / Rayleigh flow (as needed)
  - `FlowSolver` – solves flow regimes with three flexible input options
  - `PlotFunc` – visualization of results & animated back‑pressure effect
- Solver outputs:
  - Array with `[x, M, P/P0, T/T0, rho/rho0, P0/Pref]`
  - Flow type classification:
    - Fully subsonic
    - Choked subsonic
    - Normal shock inside nozzle
    - Fully supersonic (under‑/ideal/over‑expanded)

## Inputs
| Category | Parameters |
|----------|------------|
| Nozzle geometry | Converging length `Lc`, diverging length `Ld`, inlet radius `Ri`, throat radius `Rth`, exhaust radius `Re` |
| Gas properties | Specific heat ratio `GAMMA` (γ), gas constant `R` |
| Inlet conditions | Stagnation pressure `P0` |
| Ambient (back) pressure | Static back pressure `Pe` |

## Dependencies
- Python 3.7+
- `numpy`
- `matplotlib`
- `scipy`

Install them with:
```bash
pip install numpy matplotlib scipy
```

## Usage
__See the example script `Demo_CD_Nozzle.py` in the `examples/` directory. It demonstrates how to import classes, define inputs, call the solver, and plot results.__

The `FlowSolver.solver_flow` function accepts three modes:

  `solve_flow(Pe=<user defined value>)` – provide ambient static pressure value

  `solve_flow(shock_position=<user defined value>)` – provide normal shock location along nozzle axis

  `solve_flow(M_throat=<user defined value>)` – provide throat Mach number (typically 1 for choked flow)

Results can be visualized using `PlotFunc.plot_nozzle_results`.
An animation of back‑pressure effects can be generated with `PlotFunc.create_flow_animation` using an array `Pe_range`.

## Outputs
| Column | Quantity | Description |
|--------|----------|-------------|
| `x` | Distance (m) | Along nozzle axis |
| `M` | Mach number | Local Mach number |
| `P/P0` | Pressure ratio | Static pressure / stagnation pressure |
| `T/T0` | Temperature ratio | Static temperature / stagnation temperature |
| `rho/rho0` | Density ratio | Static density / stagnation density |
| `P0/Pref` | Total pressure ratio | Local total pressure / inlet total pressure (indicates shock losses) |

`flow_type` string indicates the regime (e.g., `Fully Supersonic (Under Expanded)`).

## Use Cases
   __Education__ – understand 1D gas dynamics, normal shocks, over/under‑expansion, and back‑pressure effects.

  __Preliminary design__ – evaluate nozzle performance for given geometry and operating conditions.
