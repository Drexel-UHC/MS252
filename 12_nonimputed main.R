library(haven)     # For reading .sas7bdat files
library(dplyr)     # For data manipulation
library(lubridate) # For date manipulation
library(dlnm)
library(splines)
library(RColorBrewer)
library(tsModel)
library(data.table) # HANDLE LARGE DATASETS
library(dlnm) ; library(gnm) ; library(splines) # MODELLING TOOLS
library(sf) ; library(terra) # HANDLE SPATIAL DATA
library(exactextractr) # FAST EXTRACTION OF AREA-WEIGHTED RASTER CELLS
library(dplyr) ; library(tidyr) # DATA MANAGEMENT TOOLS
library(ggplot2) ; library(patchwork) # PLOTTING TOOLS
library(gnm)
library(stargazer)
library(broom)
#########################################Derek's 00_aggregate_data.R###########################################################
library("tidyverse")
library("lubridate")

# City Road Mortality Data ------------------------------------------------
# Path to data
path_1 <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/Non-imputed/Prelim Derived Data_20240920/"

# List of data files
data_files_1 <- list.files(path_1) %>%
  str_subset(".sas7bdat")

# Function to add year_month_dow (day of week) column
add_year_month_dow <- function(df) {
  df %>% 
    mutate(age_cat = factor(age_cat),
           year = str_extract(allDate, "\\d\\d\\d\\d"),
           year_month = str_extract(allDate, "\\d\\d\\d\\d-\\d\\d"),
           doy = yday(allDate),
           dow = wday(allDate, label = TRUE),
           year_month_dow = paste(year_month, dow, sep = "-"))
}

# Function to extract relevant data and aggregate deaths
custom_summary_function <- function(df) {
  df %>% 
    summarize(across(c(country,
                       city_size,
                       salid1,
                       year_month,
                       dow), first),
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
# Each row represents a city-year-month-day
df <- data_files_1 %>% 
  map_dfr(\(x) {
    haven::read_sas(paste0(path_1, x)) %>% 
      add_year_month_dow() %>% 
      group_by(allDate) %>% 
      custom_summary_function()
  })

# Save data
path_write <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/Non-imputed/"
#write_csv(df, paste0(path_write, "combined_subgroup.csv"))

# Read saved data
df <- read_csv(paste0(path_write, "combined_subgroup.csv"))


# City Level Modifier Data ------------------------------------------------
# Path to data
path_2 <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252 requested data/"

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
write_csv(df_all, paste0(path_write, "final_subgroup.csv"))

# Read data
df_all <- read_csv(paste0(path_write, "final_subgroup.csv"))

# Temperature Cluster Data ------------------------------------------------
path_3 <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/"

temp_cluster <- haven::read_sas(paste0(path_3, "city_level_temp_w_clusters.sas7bdat"))

df_temp_cluster <- df_all %>% 
  left_join(temp_cluster, by = join_by(salid1 == nsalid1))

# Save data
#write_csv(df_temp_cluster, paste0(path_write, "final_temp_cluster_subgroup.csv"))
##################################################
fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}

non_imputed <- read.csv("/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/Non-imputed/final_temp_cluster_subgroup.csv")
non_imputed <- as.data.table(non_imputed)
summary(non_imputed)

# Model set-up (using non-imputed data)
data <- non_imputed %>% filter(!grepl("^2020", `year_month`))

data[, stratum:=factor(paste(salid1, year_month, dow, sep=":"))]
data[,  keep:=sum(median_road_round)>0, by=stratum]

Temp_measure <- "tmp_pw_percentile"
minT <- min(data[[Temp_measure]], na.rm=TRUE)
maxT <- max(data[[Temp_measure]], na.rm=TRUE)
medT <- median(data[[Temp_measure]], na.rm=TRUE)
quan01 <- quantile(data[[Temp_measure]], probs = 0.01, na.rm = TRUE)
quan25 <- quantile(data[[Temp_measure]], probs = 0.25, na.rm = TRUE)
quan75 <- quantile(data[[Temp_measure]], probs = 0.75, na.rm = TRUE)
quan99 <- quantile(data[[Temp_measure]], probs = 0.99, na.rm = TRUE)

Temp_measure_celcius <- "L1ADtemp_pw"
minT_celcius <- min(data[[Temp_measure_celcius]], na.rm=TRUE)
maxT_celcius <- max(data[[Temp_measure_celcius]], na.rm=TRUE)
medT_celcius <- median(data[[Temp_measure_celcius]], na.rm=TRUE)
quan01_celcius <- quantile(data[[Temp_measure_celcius]], probs = 0.01, na.rm = TRUE)
quan25_celcius <- quantile(data[[Temp_measure_celcius]], probs = 0.25, na.rm = TRUE)
quan75_celcius <- quantile(data[[Temp_measure_celcius]], probs = 0.75, na.rm = TRUE)
quan99_celcius <- quantile(data[[Temp_measure_celcius]], probs = 0.99, na.rm = TRUE)

# Modeling (non-imputed; "original model")
lagknots <- logknots(3, df = 3)
knotstmean <- quantile(data[[Temp_measure]], c(10,75,90)/100, na.rm=T)
argvartmean <- list(fun="ns", knots=knotstmean)
cbt <- crossbasis(data[[Temp_measure]], lag = 2, argvar = argvartmean, arglag = list(knots = lagknots), group=data$salid1)
model_road <- gnm(median_road_round ~ cbt, eliminate = stratum, 
                  family = quasipoisson(), data = data, subset=keep)

lagknots_celcius <- logknots(3, df = 3)
knotstmean_celcius <- quantile(data[[Temp_measure_celcius]], c(10,75,90)/100, na.rm=T)
argvartmean_celcius <- list(fun="ns", knots=knotstmean_celcius)
cbt_celcius <- crossbasis(data[[Temp_measure_celcius]], lag = 2, argvar = argvartmean_celcius, arglag = list(knots = lagknots_celcius), group=data$salid1)
model_road_celcius <- gnm(median_road_round ~ cbt_celcius, eliminate = stratum, 
                  family = quasipoisson(), data = data, subset=keep)

pred_road <- crosspred(cbt, model_road, cum=TRUE, cen=medT)
pred_road_celcius <- crosspred(cbt_celcius, model_road_celcius, cum=TRUE, cen=medT_celcius)


# Plotting
col='gray'
par(mfrow = c(1, 2))
par(mar=c(4,5,1,0.5), las=1, mgp=c(2.5,1,0))
plot(pred_road, "overall", lag=0, cumul=TRUE, ylim=c(0.9,1.3), ylab="RR", col=col, lwd=1.5,ci='area',lty=1,
       xlab="Temperature (%tile) - pooled across 100 imputations", ci.arg=list(col=alpha(col, 0.3)))
ind1 <- pred_road$predvar<=medT
ind2 <- pred_road$predvar>=medT
lines(pred_road$predvar[ind1],pred_road$allRRfit[ind1],col=4,lwd=1.5)
lines(pred_road$predvar[ind2],pred_road$allRRfit[ind2],col=2,lwd=1.5)    
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")
  
plot(pred_road_celcius, "overall", lag=0, cumul=TRUE, ylim=c(0.8,1.4), ylab="RR", col=col, lwd=1.5,ci='area',lty=1,
       xlab="Temperature (°C) - pooled across 100 imputations", ci.arg=list(col=alpha(col, 0.3)))
ind1 <- pred_road_celcius$predvar<=medT_celcius
ind2 <- pred_road_celcius$predvar>=medT_celcius
lines(pred_road_celcius$predvar[ind1],pred_road_celcius$allRRfit[ind1],col=4,lwd=1.5)
lines(pred_road_celcius$predvar[ind2],pred_road_celcius$allRRfit[ind2],col=2,lwd=1.5)    
abline(v = c(quan01_celcius,quan25_celcius,medT_celcius,quan75_celcius,quan99_celcius), lty = 4, col = "gray")
  

summary(data)
