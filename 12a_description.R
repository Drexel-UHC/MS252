library(dplyr)
length(row_number(data))
non_imputed
length(unique(non_imputed[non_imputed$country == "CR"]$year))

sum(data$median_road_round)

summary_stats <- data %>% 
  group_by(salid1) %>% 
  summarise(
    # Columns that should be multiplied by 100
    across(
      c(
        median_road_round, median_vehicle_round, median_motorcycle_round, 
        median_bicycle_round, median_ped_round, male_deaths, female_deaths, 
        deaths_under_9, deaths_10_19, deaths_20_34, deaths_35_64, deaths_65_plus
      ), 
      ~ 100 * mean(.x, na.rm = TRUE), 
      .names = "mean_{col}"
    ),
    
    # Columns that should not be multiplied by 100 (L1ADtemp_pw, tmp_pw_percentile)
    across(
      c(L1ADtemp_pw, tmp_pw_percentile), 
      ~ mean(.x, na.rm = TRUE), 
      .names = "mean_{col}"
    ),
    
    .groups = "drop"
  )


summary_across_rows <- summary_stats %>%
  summarise(
    across(
      -salid1,  # Exclude 'salid1' from the summary
      list(
        mean = ~mean(.x, na.rm = TRUE),
        std = ~sd(.x, na.rm = TRUE),
        median = ~median(.x, na.rm = TRUE),
        min = ~min(.x, na.rm = TRUE),
        max = ~max(.x, na.rm = TRUE)
      ),
      .names = "{col}_{fn}"
    )
  )

# Pivot to long format first
summary_long <- summary_across_rows %>%
  pivot_longer(
    cols = everything(),
    names_to = c("variable", "statistic"),
    names_pattern = "(.*)_(mean|std|median|min|max)$",
    values_to = "value"
  )

# Pivot to wide format: variables as rows, statistics as columns
summary_wide <- summary_long %>%
  pivot_wider(
    names_from = statistic,
    values_from = value
  )

# View the reshaped summary
print(summary_wide)


# write.csv(summary_wide, "summary_across_rows.csv", row.names = FALSE)
