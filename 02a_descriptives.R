################################################################################
# Purpose:
# This script summarizes road traffic mortality data across Latin American cities
# based on the *median* of 100 imputations. It generates descriptive statistics
# overall, by city, by country, and by temperature cluster.
#
# Outputs include:
# - Summary tables (mean, SD, median, min, max) across all cities.
# - Mortality, population, and temperature metrics by city and by country.
# - Spatial plots showing city-level mortality and temperature.
# - Temperature percentiles by climate clusters.
#
# Author: Cheng-Kai Hsu
# Date: April 2025
################################################################################

################################################################################
# Purpose:
# This script summarizes road traffic mortality data across Latin American cities
# based on the *median* of 100 imputations. It generates descriptive statistics
# overall, by city, by country, and by temperature cluster.
#
# Outputs include:
# - Summary tables (mean, SD, median, min, max) across all cities.
# - Mortality, population, and temperature metrics by city and by country.
# - Spatial plots showing city-level mortality and temperature.
# - Temperature percentiles by climate clusters.
#
# Author: Cheng-Kai Hsu
# Date: April 2025
################################################################################

# Load libraries
library(dplyr)
library(data.table)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(grid)

# 1. Summary statistics across all cities ----------------------------------
summary_stats <- data %>%
  group_by(salid1) %>%
  summarise(
    across(
      c(median_road_round, median_vehicle_round, median_motorcycle_round, 
        median_bicycle_round, median_ped_round, male_deaths, female_deaths, 
        deaths_under_9, deaths_10_19, deaths_20_34, deaths_35_64, deaths_65_plus),
      ~ 100 * mean(.x, na.rm = TRUE), .names = "mean_{col}"
    ),
    across(c(L1ADtemp_pw, tmp_pw_percentile),
           ~ mean(.x, na.rm = TRUE), .names = "mean_{col}"
    ),
    .groups = "drop"
  )

summary_across_rows <- summary_stats %>%
  summarise(
    across(-salid1,
           list(mean = mean, std = sd, median = median, min = min, max = max),
           .names = "{col}_{fn}"
    )
  )

summary_long <- summary_across_rows %>%
  pivot_longer(cols = everything(),
               names_to = c("variable", "statistic"),
               names_pattern = "(.*)_(mean|std|median|min|max)",
               values_to = "value")

summary_wide <- summary_long %>%
  pivot_wider(names_from = statistic, values_from = value)

print(summary_wide)

# 2. Descriptives by city --------------------------------------------------
descriptive <- data[, .(country, salid1, year_month, median_road_round, L1ADtemp_pw)]
descriptive[, year := as.integer(substr(year_month, 1, 4))]
setDT(descriptive)
setDT(BEC1)

annual_summary <- descriptive[, .(
  annual_deaths = sum(median_road_round, na.rm = TRUE),
  annual_temp = mean(L1ADtemp_pw, na.rm = TRUE)
), by = .(salid1, year)]

pop_long <- melt(BEC1, id.vars = "SALID1",
                 measure.vars = patterns("^BECTPOP\\d{4}L1AD$"),
                 variable.name = "year_col", value.name = "population")
pop_long[, year := as.integer(gsub("BECTPOP(\\d{4})L1AD", "\\1", year_col))]
setnames(pop_long, "SALID1", "salid1")

merged <- merge(annual_summary, pop_long[, .(salid1, year, population)],
                by = c("salid1", "year"), all.x = TRUE)
merged[, mortality_rate := (annual_deaths / population) * 1e5]

city_summary <- merged[, .(
  median_deaths     = as.numeric(median(annual_deaths, na.rm = TRUE)),
  p10_deaths        = as.numeric(quantile(annual_deaths, 0.10, na.rm = TRUE)),
  p90_deaths        = as.numeric(quantile(annual_deaths, 0.90, na.rm = TRUE)),
  
  median_temp       = as.numeric(median(annual_temp, na.rm = TRUE)),
  p10_temp          = as.numeric(quantile(annual_temp, 0.10, na.rm = TRUE)),
  p90_temp          = as.numeric(quantile(annual_temp, 0.90, na.rm = TRUE)),
  
  median_pop_1k     = as.numeric(median(population, na.rm = TRUE)) / 1e3,
  p10_pop_1k        = as.numeric(quantile(population, 0.10, na.rm = TRUE)) / 1e3,
  p90_pop_1k        = as.numeric(quantile(population, 0.90, na.rm = TRUE)) / 1e3,
  
  median_mort_rate  = as.numeric(median(mortality_rate, na.rm = TRUE)),
  p10_mort_rate     = as.numeric(quantile(mortality_rate, 0.10, na.rm = TRUE)),
  p90_mort_rate     = as.numeric(quantile(mortality_rate, 0.90, na.rm = TRUE))
), by = salid1]


# 3. Spatial mapping --------------------------------------------------------
l1ad <- st_read("../../Data/MS252 requested data/ForSALURBALCEC20240927.gdb", layer = "SALURBAL_L1AD")
l1ad$SALID1 <- as.character(l1ad$SALID1)
city_summary[, salid1 := as.character(salid1)]

l1ad_summary <- merge(l1ad, city_summary, by.x = "SALID1", by.y = "salid1")
countries <- ne_countries(scale = "medium", returnclass = "sf")
countries_proj <- st_transform(countries, crs = st_crs(l1ad_summary))
l1ad_points <- st_centroid(l1ad_summary)
l1ad_points_wgs <- st_transform(l1ad_points, 4326)
countries_wgs <- st_transform(countries_proj, 4326)
bbox <- st_bbox(l1ad_points_wgs)

# Mortality map
ggplot() +
  geom_sf(data = countries_wgs, fill = "gray90", color = "gray20", size = 0.3) +
  geom_sf(data = l1ad_points_wgs, aes(color = median_mort_rate), size = 1, stroke = 1) +
  scale_color_viridis_c(option = "B", name = "Mortality Rate\n(per 100k)") +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE, datum = NA) +
  theme_minimal()

# Temperature map
ggplot() +
  geom_sf(data = countries_wgs, fill = "gray90", color = "gray20", size = 0.3) +
  geom_sf(data = l1ad_points_wgs, aes(color = median_temp), size = 1, stroke = 1) +
  scale_color_viridis_c(option = "B", name = "Mean Temp\n(\u00b0C)") +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE, datum = NA) +
  theme_minimal()

# 4. Descriptives by Country -----------------------------------------------
descriptive <- data[, .(country, salid1, year_month, median_road_round, L1ADtemp_pw)]
descriptive[, year := as.integer(substr(year_month, 1, 4))]
setDT(descriptive)

annual_summary <- descriptive[, .(
  annual_deaths = sum(median_road_round, na.rm = TRUE),
  annual_temp = mean(L1ADtemp_pw, na.rm = TRUE)
), by = .(salid1, year)]

pop_long <- melt(BEC1, id.vars = "SALID1",
                 measure.vars = patterns("^BECTPOP\\d{4}L1AD$"),
                 variable.name = "year_col", value.name = "population")
pop_long[, year := as.integer(gsub("BECTPOP(\\d{4})L1AD", "\\1", year_col))]
setnames(pop_long, "SALID1", "salid1")

merged <- merge(annual_summary, pop_long[, .(salid1, year, population)],
                by = c("salid1", "year"), all.x = TRUE)
merged <- merge(merged, unique(descriptive[, .(salid1, country)]), by = "salid1")
merged[, country := fifelse(country %in% c("CR", "SV", "PA"), "Central America", country)]

country_summary <- merged[, .(
  median_deaths = median(annual_deaths, na.rm = TRUE),
  p10_deaths = quantile(annual_deaths, 0.10, na.rm = TRUE),
  p90_deaths = quantile(annual_deaths, 0.90, na.rm = TRUE),
  median_temp = median(annual_temp, na.rm = TRUE),
  p10_temp = quantile(annual_temp, 0.10, na.rm = TRUE),
  p90_temp = quantile(annual_temp, 0.90, na.rm = TRUE),
  median_pop_1k = median(population, na.rm = TRUE) / 1e3,
  p10_pop_1k = quantile(population, 0.10, na.rm = TRUE) / 1e3,
  p90_pop_1k = quantile(population, 0.90, na.rm = TRUE) / 1e3,
  count = uniqueN(salid1)
), by = country]

country_table <- country_summary[, .(
  Country = country,
  Deaths = sprintf("%d (%.1f, %.1f)", round(median_deaths), p10_deaths, p90_deaths),
  Temp_C = sprintf("%.2f (%.2f, %.2f)", median_temp, p10_temp, p90_temp),
  Population_1k = sprintf("%d (%d, %d)", round(median_pop_1k), round(p10_pop_1k), round(p90_pop_1k)),
  Cities = count
)]

print(country_table)

# 5. Descriptives by Temperature Cluster -----------------------------------
temp_percentiles_by_cluster <- data %>%
  group_by(cluster_ward_std_6) %>%
  summarize(
    p01 = quantile(L1ADtemp_pw, 0.01, na.rm = TRUE),
    p05 = quantile(L1ADtemp_pw, 0.05, na.rm = TRUE),
    p25 = quantile(L1ADtemp_pw, 0.25, na.rm = TRUE),
    p50 = quantile(L1ADtemp_pw, 0.50, na.rm = TRUE),
    p75 = quantile(L1ADtemp_pw, 0.75, na.rm = TRUE),
    p95 = quantile(L1ADtemp_pw, 0.95, na.rm = TRUE),
    p99 = quantile(L1ADtemp_pw, 0.99, na.rm = TRUE),
    .groups = "drop"
  )

print(round(temp_percentiles_by_cluster, 1))
