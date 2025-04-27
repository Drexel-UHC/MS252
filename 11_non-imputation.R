################################################################################
# Purpose:
# This script processes non-imputed road death data from the SALURBAL project
# and repeated the analysis as previously done, as a sensitivity analysis. 
################################################################################

path_4 <- "/Volumes/TOSHIBA Kai/MS252/Non-imputed_Nature Cities/Derived Data/"

# Get full file paths
data_files_4 <- list.files(path_4, pattern = ".sas7bdat$", full.names = TRUE)

# Print total number of files found
cat("Number of SAS files found:", length(data_files_4), "\n")

# Function to safely process a file
process_file <- function(file_path) {
  cat("Processing:", basename(file_path), "\n")
  
  # Try reading file
  dt <- tryCatch({
    as.data.table(read_sas(file_path))
  }, error = function(e) {
    cat("Error reading:", basename(file_path), ":", conditionMessage(e), "\n")
    return(NULL)
  })
  
  # Skip if failed
  if (is.null(dt) || !"salid1" %in% names(dt)) return(NULL)
  
  # Ensure allDate is Date class
  dt[, allDate := as.Date(allDate)]
  
  # Add time variables
  dt[, `:=`(
    year = year(allDate),
    year_month = format(allDate, "%Y-%m"),
    dow = lubridate::wday(allDate, label = TRUE),
    year_month_dow = paste0(format(allDate, "%Y-%m"), "-", lubridate::wday(allDate, label = TRUE))
  )]
  
  # Collapse by city-day
  dt_agg <- dt[, .(
    country = first(country),
    city_size = first(city_size),
    year_month = first(year_month),
    dow = first(dow),
    deaths = sum(deaths, na.rm = TRUE),
    road = sum(road, na.rm = TRUE),
    tmp_pw_percentile = mean(tmp_pw_percentile, na.rm = TRUE),
    L1ADtemp_pw = mean(L1ADtemp_pw, na.rm = TRUE)
  ), by = .(salid1, allDate, year_month_dow)]
  
  return(dt_agg)
}

# Loop through files and combine
nonimp_road_deaths <- map_dfr(data_files_4, process_file)

nonimp_road_deaths


# Exclude 2020 and prepare strata
data_nonimp_road_deaths <- nonimp_road_deaths %>% filter(!grepl("^2020", year_month)) %>% as.data.table()
data_nonimp_road_deaths[, stratum := factor(paste(salid1, year_month, dow, sep = ":"))]
data_nonimp_road_deaths[, keep := sum(road) > 0, by = stratum]

# Temperature Variables Setup --------------------------------------------
Temp_measure_nonimp_road_deaths <- "tmp_pw_percentile"

# Calculate quantiles (percentiles)
minT_nonimp_road_deaths <- min(data_nonimp_road_deaths[[Temp_measure_nonimp_road_deaths]], na.rm=TRUE)
maxT_nonimp_road_deaths <- max(data_nonimp_road_deaths[[Temp_measure_nonimp_road_deaths]], na.rm=TRUE)
medT_nonimp_road_deaths <- median(data_nonimp_road_deaths[[Temp_measure_nonimp_road_deaths]], na.rm=TRUE)
quan01_nonimp_road_deaths <- quantile(data_nonimp_road_deaths[[Temp_measure_nonimp_road_deaths]], probs = 0.01, na.rm = TRUE)
quan25_nonimp_road_deaths <- quantile(data_nonimp_road_deaths[[Temp_measure_nonimp_road_deaths]], probs = 0.25, na.rm = TRUE)
quan75_nonimp_road_deaths <- quantile(data_nonimp_road_deaths[[Temp_measure_nonimp_road_deaths]], probs = 0.75, na.rm = TRUE)
quan99_nonimp_road_deaths <- quantile(data_nonimp_road_deaths[[Temp_measure_nonimp_road_deaths]], probs = 0.99, na.rm = TRUE)

# Modeling
lagknots_nonimp_road_deaths <- logknots(3, df = 3)
knotstmean_nonimp_road_deaths <- quantile(data_nonimp_road_deaths[[Temp_measure_nonimp_road_deaths]], c(10,75,90)/100, na.rm=T)
argvartmean_nonimp_road_deaths <- list(fun="ns", knots=knotstmean)
cbt_nonimp_road_deaths <- crossbasis(data_nonimp_road_deaths[[Temp_measure_nonimp_road_deaths]], lag = 2, 
                                     argvar = argvartmean_nonimp_road_deaths, arglag = list(knots = lagknots_nonimp_road_deaths), group=data_nonimp_road_deaths$salid1)
model_nonimp_road_deaths <- gnm(road ~ cbt_nonimp_road_deaths, eliminate = stratum, 
                  family = quasipoisson(), data = data_nonimp_road_deaths, subset=keep)

pred_nonimp_road_deaths <- crosspred(cbt_nonimp_road_deaths, model_nonimp_road_deaths, cum=TRUE, cen=quan01, by=0.1)


# Plotting
col <- 'gray'
  par(mfrow = c(1, 1))
  par(mar = c(4, 5, 1, 0.5), las = 1, mgp = c(3, 1, 0))
  
  # Plot a
  plot(pred_nonimp_road_deaths, "overall", lag = 0, cumul = TRUE, ylim = c(0.9, 1.4), 
       ylab = "RR", col = col, lwd = 1.5, ci = 'area', lty = 1,
       xlab = "Temperature (%tile)", ci.arg = list(col = alpha(col, 0.3)))
  ind2 <- pred_nonimp_road_deaths$predvar >= quan01_nonimp_road_deaths & pred_nonimp_road_deaths$predvar <= quan99
  lines(pred_nonimp_road_deaths$predvar[ind2], pred_nonimp_road_deaths$allRRfit[ind2], 
        col = 'firebrick3', lwd = 1.5)
  abline(v = c(quan01_nonimp_road_deaths, quan25_nonimp_road_deaths, medT_nonimp_road_deaths, quan75_nonimp_road_deaths, quan99_nonimp_road_deaths), lty = 4, col = "gray")
  mtext("a", side = 3, adj = -0.3, line = 0, cex = 1)  # Label 'a'
  
  
