# Plotting (sensitivity by outcome) — styled to match the sensitivity by knot plots
col <- 'gray'
  
# Set up 1 row, 2 columns (you can change to mfrow = c(2, 4) in the combined version)
# par(mfrow = c(1, 2), mar = c(4, 5, 3, 1.5), las = 1, mgp = c(3, 1, 0))

# Plot a — Non-imputed road deaths
plot(pred_nonimp_road_deaths, "overall", lag = 0, cumul = TRUE, ylim = c(0.9, 1.4), 
     ylab = "RR", col = col, lwd = 1.5, ci = 'area', lty = 1,
     xlab = "Temperature (%tile)", ci.arg = list(col = alpha(col, 0.3)))
ind2 <- pred_nonimp_road_deaths$predvar >= quan01_nonimp_road_deaths & 
  pred_nonimp_road_deaths$predvar <= quan99_nonimp_road_deaths
lines(pred_nonimp_road_deaths$predvar[ind2], 
      pred_nonimp_road_deaths$allRRfit[ind2], col = 'firebrick3', lwd = 1.5)
abline(v = c(quan01_nonimp_road_deaths, quan25_nonimp_road_deaths, 
             medT_nonimp_road_deaths, quan75_nonimp_road_deaths, 
             quan99_nonimp_road_deaths), lty = 4, col = "gray")
# mtext("g", side = 3, adj = -0.3, line = 0, cex = 1)
mtext("Non-imputed road deaths", side = 3, line = 1.2, cex = 1, font = 2)

# Plot h — Median-imputation outcome
plot(pred_road, "overall", lag = 0, cumul = TRUE, ylim = c(0.9, 1.4), 
     ylab = "RR", col = col, lwd = 1.5, ci = 'area', lty = 1,
     xlab = "Temperature (%tile)", ci.arg = list(col = alpha(col, 0.3)))
ind2 <- pred_road$predvar >= quan01 & pred_road$predvar <= quan99
lines(pred_road$predvar[ind2], pred_road$allRRfit[ind2], col = 'firebrick3', lwd = 1.5)
abline(v = c(quan01, quan25, medT, quan75, quan99), lty = 4, col = "gray")
# mtext("h", side = 3, adj = -0.3, line = 0, cex = 1)
mtext("Median-imputed road deaths", side = 3, line = 1.2, cex = 1, font = 2)

dev.off()

############################################################
# Calculate standard errors from RR for indices 951 to 991
index_range <- 951:991

coef_std_imp <- (log(pred_road$allRRhigh[index_range]) - log(pred_road$allRRfit[index_range])) / 1.96
coef_std_nonimp <- (log(pred_nonimp_road_deaths$allRRhigh[index_range]) - log(pred_nonimp_road_deaths$allRRfit[index_range])) / 1.96

# Create percentile sequence for x-axis (e.g., from 95th to 99th percentile)
temp_percentile <- seq(95, 99, length.out = length(index_range))

# Create data frame
plot_data <- data.frame(
  temp_percentile = temp_percentile,
  SE_imputed = coef_std_imp,
  SE_nonimputed = coef_std_nonimp
)

# Convert to long format for ggplot2
library(tidyr)
plot_data_long <- pivot_longer(
  plot_data,
  cols = c("SE_imputed", "SE_nonimputed"),
  names_to = "Type",
  values_to = "SE"
)

# Clean up labels
plot_data_long$Type <- factor(plot_data_long$Type,
                              levels = c("SE_imputed", "SE_nonimputed"),
                              labels = c("Imputed", "Non-Imputed"))

# Plot
library(ggplot2)

ggplot(plot_data_long, aes(x = temp_percentile, y = SE, color = Type)) +
  geom_line(size = 1.1) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = seq(1, 99, by = 1)) +
  labs(
    title = "Standard Error of log(RR) Across High Temperature Percentiles",
    x = "Temperature Percentile",
    y = "Standard Error of log(RR)",
    color = "Model Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )
########################################################################################################################
library(ggplot2)
library(tidyr)
library(patchwork)

# Define index range and temperature percentiles
index_range <- 951:991
temp_percentile <- seq(95, 99, length.out = length(index_range))

# Compute SE of log(RR)
coef_std_imp <- (log(pred_road$allRRhigh[index_range]) - log(pred_road$allRRfit[index_range])) / 1.96
coef_std_nonimp <- (log(pred_nonimp_road_deaths$allRRhigh[index_range]) - log(pred_nonimp_road_deaths$allRRfit[index_range])) / 1.96

# SE dataframe
se_df <- data.frame(
  temp_percentile = temp_percentile,
  Imputed = coef_std_imp,
  "Non-Imputed" = coef_std_nonimp,
  check.names = FALSE
)

se_long <- pivot_longer(
  se_df,
  cols = c("Imputed", "Non-Imputed"),
  names_to = "Type",
  values_to = "SE"
)

# SE plot
se_plot <- ggplot(se_long, aes(x = temp_percentile, y = SE, color = Type)) +
  geom_line(size = 1.1) +
  geom_point(size = 1.5) +
  labs(
    title = "Standard Error of log(RR)",
    x = "Temperature Percentile",
    y = "SE of log(RR)",
    color = "Model Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# Compute CI widths
ci_width_imp <- pred_road$allRRhigh[index_range] - pred_road$allRRlow[index_range]
ci_width_nonimp <- pred_nonimp_road_deaths$allRRhigh[index_range] - pred_nonimp_road_deaths$allRRlow[index_range]

# CI dataframe
ci_df <- data.frame(
  temp_percentile = temp_percentile,
  Imputed = ci_width_imp,
  "Non-Imputed" = ci_width_nonimp,
  check.names = FALSE
)

ci_long <- pivot_longer(
  ci_df,
  cols = c("Imputed", "Non-Imputed"),
  names_to = "Type",
  values_to = "CI_width"
)

# CI plot
ci_plot <- ggplot(ci_long, aes(x = temp_percentile, y = CI_width, color = Type)) +
  geom_line(size = 1.1) +
  geom_point(size = 1.5) +
  labs(
    title = "Width of 95% CI for RR",
    x = "Temperature Percentile",
    y = "CI Width (RRhigh - RRlow)",
    color = "Model Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# Combine plots
se_plot 
ci_plot 
