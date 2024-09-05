library("tidyverse")
library("lubridate")

# City Road Mortality Data ------------------------------------------------
# Path to data
path_1 <- "//files.drexel.edu/colleges/SOPH/Shared/UHC/Projects/Wellcome_Trust/Manuscripts/MS252_Hsu/Data Build_20240710/Prelim Derived Data/"

# List of data files
data_files_1 <- list.files(path_1) %>%
  str_subset(".sas7bdat")

# Function to add year_month_dow (day of week) column
add_year_month_dow <- function(df) {
  df %>% 
    mutate(age_cat = factor(age_cat),
           year_month_dow = paste(str_extract(allDate, "\\d\\d\\d\\d-\\d\\d"),
                                  wday(allDate, label = TRUE),
                                  sep = "-"))
}

# Function to extract relevant data and aggregate deaths
custom_summary_function <- function(df) {
  df %>% 
    summarize(across(c(country,
                       city_size,
                       salid1, 
                       allDate), first),
              across(c(deaths,
                       median_road_round,
                       median_vehicle_round, 
                       median_motorcycle_round,
                       median_bicycle_round,
                       median_ped_round), sum),
              across(c(L1ADtemp_pw,
                       L1ADtemp_x,
                       tmp_pw_percentile,
                       tmp_x_percentile), mean),
              male_deaths    = sum(median_road*(male == 1)),
              female_deaths  = sum(median_road*(male == 0)),
              deaths_under_9 = sum(median_road*(age_cat == "<=9")),
              deaths_10_19   = sum(median_road*(age_cat == "10-19")),
              deaths_20_34   = sum(median_road*(age_cat == "20-34")),
              deaths_35_64   = sum(median_road*(age_cat == "35-64")),
              deaths_65_plus = sum(median_road*(age_cat == "65+")))
}

# Read in data and create stratified summary
# Each row represents a city-year-month-dow
df <- data_files_1 %>% 
  map_dfr(\(x) {
    haven::read_sas(paste0(path_1, x)) %>% 
      add_year_month_dow() %>% 
      group_by(year_month_dow) %>% 
      custom_summary_function() %>% 
      arrange(allDate)
  })

# Save data
path_write <- "C:/Users/dpw48/OneDrive - Drexel University/git/MS252/Data/"
write_rds(df, paste0(path_write, "combined_subgroup.rds"))

# Read saved data
df <- read_rds(paste0(path_write, "combined_subgroup.rds"))


# City Level Modifier Data ------------------------------------------------
# Path to data
path_2 <- "//files.drexel.edu/colleges/SOPH/Shared/UHC/Projects/Wellcome_Trust/Manuscripts/MS252_Hsu/2024_01_03/"

# Read and combine BEC data
BEC1 <- read_csv(paste0(path_2, "BEC_L1AD_08162023.csv")) %>% 
  select(c('ISO2','SALID1','BECADSTTLGAVGL1AD','BECCZL1AD','BECADSTTDENSL1AD',
           'BECADLRDENSL1AD','BECADINTDENSL1AD','BECPTCHDENSL1AD','BECADINTDENS3L1AD',
           'BECADINTDENS4L1AD','BECADSTTPNODEAVGL1AD','BECADSTTPNODESDL1AD',
           'BECADSTTLGAVGL1AD','BECADCRCTYAVGL1AD','BECSTTPL1AD','BECPCTURBANL1AD',
           'BECGSPCTL1AD', 'BECGSPTCHDENSL1AD','BECMINWAGEL1AD','BECELEVATIONMAXL1AD',
           'BECELEVATIOVEL1AD','BECELEVATIONMEDIANL1AD','BECELEVATIONMINL1AD',
           'BECELEVATIONP25L1AD','BECELEVATIONP75L1AD','BECELEVATIONSTDL1AD',
           'BECSLOPEMAXL1AD','BECSLOPEAVEL1AD','BECSLOPEMEDIANL1AD','BECSLOPEMINL1AD',
           'BECSLOPEP25L1AD','BECSLOPEP75L1AD','BECSLOPESTDL1AD'))

BEC2 <- read_csv(paste0(path_2, "BEC_RESTRICTED_L1AD_08162023.csv")) %>% 
  select(c('SALID1','BECURBTRVDELAYINDEXL1AD','BECURBAVGTRAFTIMEL1AD',
           'BECURBTRVDELAYTIMEL1AD','BECPARKPCTAREAL1AD'))

BEC <- inner_join(BEC1, BEC2, by = "SALID1")

# Read SEC data
SEC <- read_csv(paste0(path_2, 'SEC_INDEXSCORES_L1AD_07102023.csv')) %>% 
  filter(YEAR == max(YEAR), .by = SALID1) %>% 
  select(c('SALID1','CNSSEI_L1AD', 'CNSSE1_L1AD', 'CNSSE2_L1AD','CNSSE3_L1AD'))

# Read REG data
REG <- read_csv(paste0(path_2, 'VehicleRegistration_L1AD_20201027.csv')) %>% 
  filter(YEAR == max(YEAR), .by = SALID1) %>% 
  select(c('SALID1','BECMTRBRATEL1AD','BECPAVRATEL1AD','BECTOTVRRATEL1AD'))

# Combine road mortality and city level modifier data
df_all <- df %>% 
  left_join(BEC, by = join_by(salid1 == SALID1)) %>% 
  left_join(SEC, by = join_by(salid1 == SALID1)) %>% 
  left_join(SEC, by = join_by(salid1 == SALID1))

# Save data
write_rds(df_all, paste0(path_write, "final_subgroup.csv"))

# Read data
df_all <- read_rds(paste0(path_write, "final_subgroup.csv"))

# Temperature Cluster Data ------------------------------------------------
path_3 <- "//files.drexel.edu/colleges/SOPH/Shared/UHC/Projects/Wellcome_Trust/Manuscripts/MS252_Hsu/2024_02_27/"

temp_cluster <- haven::read_sas(paste0(path_3, "city_level_temp_w_clusters.sas7bdat"))

df_temp_cluster <- df_all %>% 
  left_join(temp_cluster, by = join_by(salid1 == nsalid1))

# Save data
write_rds(paste0(path_write, "final_temp_cluster_subgroup.csv"))
