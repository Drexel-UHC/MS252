# Heat and Road Traffic Mortality in 272 LA Cities


## Files

Files names are ordered in their intended execution order. Files
`00_aggregate_data.R`, `01_load_data.R`, and `02_model_setup.R` must be
run before any of the analysis code `03`-`06`.

| File Name                           | File Purpose                                                                                                                                                                                                                                                                                                                                                              |
|-------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `Preliminary_subgroup.ipynb`        | Process (non-imputed) data for use in analysis.                                                                                                                                                                                                                                                                                                                           |
| `00_aggregate_data.R`               | Equivalent to `Preliminary_subgroup.ipynb` (processes non-imputed data), but written in `R`. The final file is called `final_temp_cluster_subgroup.csv`. Each row is a city-year-month-day with the associated mortality, temperature, and other relevant city/time variables.                                                                                            |
| `01_load_data.R`                    | Imports all libraries and loads the data from `final_temp_cluster_subgroup.csv` into memory.                                                                                                                                                                                                                                                                              |
| `02_model_setup.R`                  | Used to specify general choices for the DLNM model, such as using percentile temperature vs absolute temperature, how many lag and variable knots to use for splines, and how to stratify the data. Generates the crossbasis used for modeling and prediction. Stratifies by city-year-month-dow (day of week) and uses natural cubic splines for nonlinear associations. |
| `03_fixed_effect_road_mortality.R`  | Runs a quassipoisson DLNM model on all forms of road mortality.                                                                                                                                                                                                                                                                                                           |
| `04_fixed_effect_mode_stratified.R` | Runs separate quassipoisson DLNM models for car, motorcycle, bicycle, and pedestrian mortality.                                                                                                                                                                                                                                                                           |
| `05_fixed_effect_age_stratified.R`  | Runs separate quassipoisson DLNM models for mortality in specified age ranges (9 or under, 10-19, 20-34, 35-64, and 65 or over).                                                                                                                                                                                                                                          |
| `06_fixed_effect_sex_stratified.R`  | Runs separate quassipoisson DLNM models for male and female mortality.                                                                                                                                                                                                                                                                                                    |

## Data

### Original Data

Description of data generated by `00_aggregate_data.R`.

### Imputed Data

Description of data generated by `10_aggregate_data.R`.

## Manuscript

These files produce the analysis found in MS252, which analyzes
associations of heat exposure with road traffic mortality in 272 Latin
American cities.
