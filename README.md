# Heat and Road Traffic Mortality in 272 LA Cities


## Files

Files names are ordered in their intended execution order. Files
`10` and `11` must be run before any of the analysis code.
`12`-`14` are for main analysis (`12` nonimputed, pooled using median; `13` imputed, using Rubin's rule; `14` simply compares models across different datasets and temp measures)
`15`-`17` for subgroup analyses by sex, age, and mode of transport.
`18` adds interaction for city-level effect modifiers

| File Name                           | File Purpose                                                                                                                                                                                     |
|-------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `10_data wrangling 1.R`             | Turns raw data (272 files, each representing a city with 100 imputations) into 27200 files (each representing a city-imputation)                                                                 |
| `11_data wrangling 2.R`             | Turns 27200 files (each representing a city-imputation) into 100 files (each representing an imputation for 272 cities)                                                                          |
| `12_nonimputed main.R`              | Define strata and model setup, run nonimputed analysis for percentile and °C.                                                                                                                    |
| `13_imputed main rubin.R`           | Run imputed analysis (Rubin's rule) and plot %tile results (non-imputed and imputed) and °C results (imputed and temp distribution)                                                              |
| `14_compare og vs pooled.R`         | Only for internal check. Simply compares across original models (percentiles and °C) and models with full and reduced coefs.                                                                     |
| `15_subgroup sex.R`                 | Runs separate quassipoisson DLNM models for male and female mortality. (1st block: non-imputed; 2nd block: imputed (takes long); 3rd block: reconstruction; fourth block: plotting)              |
| `16_subgroup_age.R`                 | Runs separate quassipoisson DLNM models for mortality in specified age ranges (9 or under, 10-19, 20-34, 35-64, and 65 or over). (1st block: non-imputed; 2nd block: imputed (takes long); 3rd block: reconstruction; fourth block: plotting)   |
| `17_subgroup_mode.R`                | Runs separate quassipoisson DLNM models for car, motorcycle, bicycle, and pedestrian mortality. (1st block: non-imputed; 2nd block: imputed (takes long); 3rd block: reconstruction; fourth block: plotting)                                    |
| `18_interaction.R`                  | Includes interaction with city-level effect modifiers (i.e., temp mean, temp std, street segment length and peak-hour travel time)                                                               |
| `19_climate cluster.R`              | Stratifies imputed analysis by six temperature clusters (cluster_ward_std_6)                                                                                                                     |
| `19B_cluster plotting.R`            | Plot cluster-specific association curve and temperature distribution                                                                                                                             |
| `20_EDF.R`                          | Calculate EDFs                                                                                                                                                                                   |



## Manuscript

These files produce the analysis found in MS252, which analyzes
associations of heat exposure with road traffic mortality in 272 Latin
American cities.
