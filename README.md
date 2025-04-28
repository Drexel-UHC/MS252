# MS252
**"Effects of ambient temperature on road traffic mortality in Latin America: Individual- and city-level variations across 272 SALURBAL cities"**

This repository contains the full modeling pipeline used for our multi-city time-series analysis of temperature and road-traffic deaths across 272 cities in Latin America. The pipeline includes data processing, stratified and interaction analyses, pooled effect estimation, and sensitivity checks.

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
| `11a_sensitivity_outcome.R` | Compares effect estimates between imputed vs. non-imputed outcomes. |

💻 System Requirements
Software:
R version ≥ 4.2.0 (tested with 4.2.2 and 4.3.1)

Operating Systems tested:

macOS Monterey & Ventura

Linux (Ubuntu 22.04)

Windows-compatible with minor path edits

Required R Packages:
r
Copy
Edit
install.packages(c("data.table", "dlnm", "gnm", "ggplot2", "patchwork",
                   "lubridate", "haven", "RColorBrewer", "tidyr"))
Hardware:
Recommended: ≥16GB RAM for running all 100 imputation datasets efficiently

No non-standard hardware required

⚙️ Installation Guide
Clone the repository:

bash
Copy
Edit
git clone https://github.com/Drexel-UHC/MS252.git
Open R or RStudio and set the working directory:

r
Copy
Edit
setwd("path/to/MS252")
Install all required packages (see above).

Update paths in scripts that reference:

r
Copy
Edit
/Volumes/TOSHIBA Kai/MS252/...
Tip: Use here::here() or define a base_path <- "your_project_folder" to simplify.

⏱️ Typical install time:
R package setup: ~5–10 minutes on stable internet

File download time varies by size (raw mortality + imputation files not stored in repo)

🚀 Demo Instructions (Optional)
To test the pipeline structure on a lightweight dataset:

Use demo_run.R (or similar script) to:

Simulate 5 cities × 60 days × 3 imputations

Run DLNM main effect

Perform Rubin’s rule pooling

Generate subgroup and interaction plots

Expected output:

3-panel plot of:

Main effect

Subgroup (e.g., male deaths)

Interaction (e.g., city-level mean temp)

⏱️ Run time: ~1–2 minutes on standard desktop (8GB RAM)

📘 Instructions for Use (Full Pipeline)
Step 1: Data Preparation
Run 00_data wrangling_1.R → processes .sas7bdat files

Run 01_data_wrangling_2.R → aggregates into imputation-level datasets

Step 2: Main Effects & Stratification
Run 02_median_imputation_main_effects.R for median imputation results

Run 03_rubin_imputation_main_effect.R for pooled main effects

Subgroup Analyses:
Sex: 05_subgroup_sex.R

Age: 06_subgroup_age.R

Mode: 07_subgroup_mode.R

Step 3: Sensitivity Analyses
Knot placement: 04_sensitivity_knots.R

Non-imputed data: 11_non-imputation.R + 11a_sensitivity_outcome.R

Step 4: Interaction Effects
Model: 08_interaction.R

Plot: 08a_interaction_plot.R

Formal tests via ANOVA or pooled coefficient inference

Step 5: Temperature Cluster Modeling
Modeling: 09_temp_cluster.R

Plotting: 09a_cluster_plotting.R

Step 6: Attributable Risk (EDF)
EDF estimation: 10_EDF.R (includes subgroups + simulation-based CI)

🔁 Reproducibility
All modeling scripts are modular and reproducible across systems with R ≥ 4.2.0.

All model coefficients, variance-covariance matrices, and prediction objects are saved in .rds format.

Outputs are version-controlled by file name and stored in clearly labeled subfolders:

Derived Data/, Figures/, Pooled results/, Median-imputed/

