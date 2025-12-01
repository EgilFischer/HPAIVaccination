# HPAI_vaccination
Codes for the Public Private Partnership on HPAI vaccination.



---

## Overview of Within-farm scripts 🐔

| Step | Script | Description | Approx. Run Time (*48 scenarios) |
|------|--------|-------------|----------------|
| 1 | `faster_multitypetransitionsSEIR_tleap.R` | SEIR simulation with vaccination & egg dynamics | ~40–50 min |
| 2 | `Scenarios.R` | Run multiple scenarios (farm size, intro time, strain) | ~32–40* hrs |
| 3 | `pre_processing_detModule.R` | Aggregate results to daily level; calculate mortality metrics | ~20* min |
| 4 | `detectionModule.R` | Evaluate passive and active detection protocols | ~5* min |
| 5 | `post_detMod_figures.R` | All plotting scripts | Not long |

---

## Scripts and Their Functions

### 1. `faster_multitypetransitionsSEIR_tleap.R`
Stochastic SEIR model for within-farm HPAI transmission in layers.

**Key Features:**
- Two bird types:
  - Type 1: Unvaccinated
  - Type 2: Vaccinated
- Gamma-distributed latent and infectious periods
- Immunity buildup and waning 
- Daily egg production (healthy vs. infected)
- Eggs picked up at regular intervals

**Key Parameters:**
- `beta`: 2×2 transmission matrix
- `infectious.period`, `latency.period`, `k.*`: Gamma (Erlang) parameters
- `p.protect`, `trans.mean.wane`, `trans.mean.buildup`: Vaccination
- `eh`, `ei`: Egg-laying rates
- `pickup_time`: Shipment interval

**Outputs:**
- S, E, I, R compartments by type
- Deaths by compartment
- Daily and cumulative egg counts (healthy and infected)

---

### 2. `Scenarios.R`
Runs multiple simulations varying:
- Flock size 
- Introduction time
- Waning immunity for homogeneous and heterogeneous strains

**Scenarios:**
- Scenario 0: Baseline (no vaccination)
- Scenario 1: Size and introduction time (100% vx)
- Scenario 2.1: Homogeneous waning and introduction time (100% vx)
- Scenario 2.2: Heterogeneous waning and introduction time (100% vx)

**Outputs:**
- `.RDS` files in `./output/`
- Plots in `./figures/`

---

### 3. `pre_processing_detModule.R`
Processes many timesteps per day simulation output to daily summaries.

**New Columns Added:**
- `day`, `N.dead`, `N.total`
- `daily.death.incidence`, `daily.mortality.rate`
- `N.dead.detectables`, `N.live.detectables`

**Outputs:**
- Daily-aggregated `.RDS` files prefixed with `DailyData.`

---

### 4. `detectionModule.R`
Applies passive and active surveillance protocols to estimate outbreak detection timing and contaminated egg risk.

**Passive Detection:**
- Mortality threshold method (2-day 0.05% mortality)
- Ratio method (2-day vs. previous week)

**Active Detection:**
- Sampling frequency (e.g., every 7, 14, 30 days)
- Dead and live confirmatory tests

**Outputs:**
- Detection timing (`min.det.time`)
- Contaminated eggs on-farm and shipped at detection
- Final outbreak size classification (MAJOR/MINOR)

---

### 5. `post_detMod_figures.R`
Plots for DailyData (faster plotting) and Detection Module outputs

---

## Unit Testing

### `Tests_Model.R`
Validates `faster_multitypetransitionsSEIR_tleap.R` behaviour.  

### `pmaj_calculation.R`
Calculation of `pmajor` using:
- Gamma-distributed infectious period
- 2-type population (vaccinated/unvaccinated)
- Extinction probabilities (`q1`, `q2`)
- Configurable `n`, `k.infectious`, and `p`

