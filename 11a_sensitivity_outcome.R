# Plotting (sensitivity by outcome) — styled to match the sensitivity by knot plots
col <- 'gray'
  
# Set up 1 row, 2 columns (you can change to mfrow = c(2, 4) in the combined version)
par(mfrow = c(1, 2), mar = c(4, 5, 3, 1.5), las = 1, mgp = c(3, 1, 0))

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
