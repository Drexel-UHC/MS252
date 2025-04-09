# Submission_NatCities
# MS252: Temperature and Road-Traffic Mortality Across Latin American Cities
**"Effects of ambient temperature on road-traffic mortality across Latin America: Individual and city-level variations across 272 cities from the SALURBAL study."**

## 📁 Script Overview

| Script | Description |
|--------|-------------|
| `00_data wrangling_1.R` | Loads and reshapes city-level 100-imputatoins road-traffic mortality data. |
| `01_data_wrangling_2.R` | Aggregates previously processed city-imputation CSV files into imputation-level datasets. |
| `02_median_imputation_main_effects.R` | Runs DLNM models using median-imputed data across all cities; estimates main effects. |
| `03_rubin_imputation_main_effect.R` | Applies Rubin’s Rules to combine results from 100 imputed datasets; produces pooled estimates for main effects. |
| `04_sensitivity_knots.R` | Conducts sensitivity analyses by testing alternative spline knot specifications in DLNM models. |
| `05_subgroup_sex.R` | Stratifies analysis by sex (male/female). |
| `06_subgroup_age.R` | Stratifies analysis by age group (<=9, 10-19, 20–39, 40–64, 65+). |
| `07_subgroup_mode.R` | Stratifies analysis by transport mode (e.g., pedestrian, motorcyclist, bicyclists, motor vehicle). |
| `08_interaction.R` | Explores effect modification by city-level characteristics (e.g., average temperature, peak-hour travel times). |
| `09_temp_cluster.R` | Explores effect modification by temperature clusters. |
| `10_EDF.R` | Generates excess death fractions (EDF) plots for main and subgroup analyses. |
| `11_non-imputation.R` | Runs non-imputed models as a sensitivity check and comparison against imputed results. |
