# IDL – Ice-Retreat Model with Hydrated Salt

This repository provides an IDL implementation of the model described in Landis et al. (2025), *"Role of Natron in Delaying Retreat of Buried Ice Tables on Ceres."*

## `script_iteration1`

This script estimates the delay in ice table descent as a function of natron abundance in the regolith.

- Regolith parameters (porosity, grain size, temperature) are hard-coded in `salty_ice_loss1.pro`. To modify these, you must edit that file directly.
- The delay is defined as the additional time required for the ice table to retreat to a given depth compared to a salt-free scenario.
- Output: A plot showing delay as a function of depth, with natron weight fractions ranging from 0.1 to 0.9.

There are no input parameters other than a background color toggle for the plot.

## `salty_ice_loss1`

### Description

Models ice retreat under the influence of vapor released by hydrated salt. The presence of water vapor from salt dehydration reduces ice loss by suppressing sublimation.

### Inputs

- `zf` – Ice table depth (m)  
- `wsalt` – Salt weight fraction (g/g); default: `0.05`  
- `pp_profile` – Structure defining partial pressure profile as a function of depth (Pa). To ignore pressure effects:

  ```idl
  pp_profile = {z:[0d,1d], pp:[0d,0d]}
  ```

  This is the default if the keyword is unset.

- `verbose` – If set, enables plots and diagnostic output  

> Note: Regolith parameters and salt composition are hard-coded in the script.

### Outputs

Returned structure includes:

- `.Pp` – Estimated partial pressure from salt dehydration at `zf` (Pa)  
- `.tt` – Vector of recession times (s); `tt[-1]` is the total time to reach depth `zf`

### Usage

See `script_iteration1.pro` for an example call.

---

## `ice_loss_rate`

### Description

Calculates the instantaneous ice retreat rate for a homogeneous, porous planetary surface under vacuum, following the formulation of Schorghofer (2008).

### Inputs  
(*Keywords with default values*)

- `temperature` – Surface/layer temperature (K) (default: `293`)  
- `depth` – Ice table depth (m) (default: `0.001`)  
- `porosity` – Regolith porosity (default: `0.5`)  
- `tortuosity` – Tortuosity of the pore network (default: `2`)  
- `diameter` – Grain or pore diameter (m) (default: `1e-4`)  
- `Pp` – Partial vapor pressure at the ice surface from sources other than ice (e.g., hydrated salt)  

*Note:*

- When `Pp ≠ 0` and `Pp < Ps`, sublimation rate is reduced.  
- If `Pp > Ps`, the result is negative mass flux, indicating ice accumulation.

### Outputs  
(structure)

- `.Ps` – Equilibrium vapor pressure of ice (Pa)  
- `.Pp` – Input partial vapor pressure (Pa)  
- `.vbar` – Mean thermal velocity (m/s)  
- `.Dk` – Effective diffusion coefficient (m²/s)  
- `.MFP` – Molecular mean free path: `2*Dk/vbar` (m)  
- `.DT` – Mean time between molecular collisions: `.MFP/vbar` (s)  
- `.Jbar` – Surface water mass flux (kg/m²/s)  
- `.Ri` – Instantaneous ice recession velocity (m/s)

### Example

Schorghofer (2008) example:  
Assuming particle diameter = 100 µm, tortuosity = 2, porosity = 0.5, depth = 1 m, temperature = 150 K

```idl
r = ice_loss_rate(temperature=150., depth=1)
print, r.Ri * 1e9 * 3.154e+7  ; Convert to meters per billion years
; Expected result: ~8.5 m/Gy
```