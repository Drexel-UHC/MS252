# Heat and Road Traffic Mortality in 272 LA Cities


## Files

| File Name                          | File Purpose                                                                                                                                                                                                                                                                   |
|------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `Preliminary_subgroup.ipynb`       | Process (non-imputed) data for use in analysis.                                                                                                                                                                                                                                |
| `00_aggregate_data.R`              | Equivalent to `Preliminary_subgroup.ipynb` (processes non-imputed data), but written in `R`. The final file is called `final_temp_cluster_subgroup.csv`. Each row is a city-year-month-day with the associated mortality, temperature, and other relevant city/time variables. |
| `01_load_data.R`                   | Called via `source` in other files. Imports all libraries and loads the data from `final_temp_cluster_subgroup.csv` into memory.                                                                                                                                               |
| `02_fixed_effect_road_mortality.R` | Runs a quassipoisson DLNM model on the loaded data. Stratified by city-year-month-dow (day of week) and uses natural cubic splines for nonlinear associations.                                                                                                                 |

## Manuscript

These files produce the analysis found in MS252, which analyzes
associations of heat exposure with road traffic mortality in 272 Latin
American cities.
