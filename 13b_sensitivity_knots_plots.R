# Main
knotstmean <- quantile(data[["tmp_pw_percentile"]], c(10, 75, 90) / 100, na.rm=TRUE)
argvartmean <- list(fun = "ns", knots = knotstmean)
cbt <- crossbasis(data[[Temp_measure]], lag = 2, argvar = argvartmean, arglag = list(knots = lagknots), group=data$salid1)
path_for_pooled <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/imputed/Pooled results/"
main <- crosspred(cbt, model.link="log",
                                          coef = readRDS(paste0(path_for_pooled, "pooled_main_results_coef.rds")),
                                          vcov = readRDS(paste0(path_for_pooled, "pooled_main_results_vcov.rds")),
                                          cum=TRUE,cen=quan01,by=0.1)

path_for_pooled <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/imputed/Pooled results/sensitivity/"
# 10, 75, 90
knotstmean <- quantile(data[["tmp_pw_percentile"]], c(10, 75, 90) / 100, na.rm=TRUE)
argvartmean <- list(fun = "ns", knots = knotstmean)
cbt <- crossbasis(data[[Temp_measure]], lag = 2, argvar = argvartmean, arglag = list(knots = lagknots), group=data$salid1)
knots_10_75_90 <- crosspred(cbt, model.link="log",
                            coef = readRDS(paste0(path_for_pooled, "knots_10_75_90_pooled_main_results_coef.rds")),
                            vcov = readRDS(paste0(path_for_pooled, "knots_10_75_90_pooled_main_results_vcov.rds")),
                            cum=TRUE,cen=quan01,by=0.1)
# 10, 50, 90
knotstmean <- quantile(data[["tmp_pw_percentile"]], c(10, 50, 90) / 100, na.rm=TRUE)
argvartmean <- list(fun = "ns", knots = knotstmean)
cbt <- crossbasis(data[[Temp_measure]], lag = 2, argvar = argvartmean, arglag = list(knots = lagknots), group=data$salid1)
knots_10_50_90 <- crosspred(cbt, model.link="log",
                            coef = readRDS(paste0(path_for_pooled, "knots_10_50_90_pooled_main_results_coef.rds")),
                            vcov = readRDS(paste0(path_for_pooled, "knots_10_50_90_pooled_main_results_vcov.rds")),
                            cum=TRUE,cen=quan01,by=0.1)
# 25, 50, 75
knotstmean <- quantile(data[["tmp_pw_percentile"]], c(25, 50, 75) / 100, na.rm=TRUE)
argvartmean <- list(fun = "ns", knots = knotstmean)
cbt <- crossbasis(data[[Temp_measure]], lag = 2, argvar = argvartmean, arglag = list(knots = lagknots), group=data$salid1)
knots_25_50_75 <- crosspred(cbt, model.link="log",
                            coef = readRDS(paste0(path_for_pooled, "knots_25_50_75_pooled_main_results_coef.rds")),
                            vcov = readRDS(paste0(path_for_pooled, "knots_25_50_75_pooled_main_results_vcov.rds")),
                            cum=TRUE,cen=quan01,by=0.1)

# 10, 90
knotstmean <- quantile(data[["tmp_pw_percentile"]], c(10, 90) / 100, na.rm=TRUE)
argvartmean <- list(fun = "ns", knots = knotstmean)
cbt <- crossbasis(data[[Temp_measure]], lag = 2, argvar = argvartmean, arglag = list(knots = lagknots), group=data$salid1)
knots_10_90 <- crosspred(cbt, model.link="log",
                         coef = readRDS(paste0(path_for_pooled, "knots_10_90_pooled_main_results_coef.rds")),
                         vcov = readRDS(paste0(path_for_pooled, "knots_10_90_pooled_main_results_vcov.rds")),
                         cum=TRUE,cen=quan01,by=0.1)

# 25, 75
knotstmean <- quantile(data[["tmp_pw_percentile"]], c(25, 75) / 100, na.rm=TRUE)
argvartmean <- list(fun = "ns", knots = knotstmean)
cbt <- crossbasis(data[[Temp_measure]], lag = 2, argvar = argvartmean, arglag = list(knots = lagknots), group=data$salid1)
knots_25_75 <- crosspred(cbt, model.link="log",
                         coef = readRDS(paste0(path_for_pooled, "knots_25_75_pooled_main_results_coef.rds")),
                         vcov = readRDS(paste0(path_for_pooled, "knots_25_75_pooled_main_results_vcov.rds")),
                         cum=TRUE,cen=quan01,by=0.1)

# 33, 67
knotstmean <- quantile(data[["tmp_pw_percentile"]], c(33.3, 66.7) / 100, na.rm=TRUE)
argvartmean <- list(fun = "ns", knots = knotstmean)
cbt <- crossbasis(data[[Temp_measure]], lag = 2, argvar = argvartmean, arglag = list(knots = lagknots), group=data$salid1)

knots_33_67 <- crosspred(cbt, model.link="log",
                         coef = readRDS(paste0(path_for_pooled, "knots_33_67_pooled_main_results_coef.rds")),
                         vcov = readRDS(paste0(path_for_pooled, "knots_33_67_pooled_main_results_vcov.rds")),
                         cum=TRUE,cen=quan01,by=0.1)

# 50
knotstmean <- quantile(data[["tmp_pw_percentile"]], c(50) / 100, na.rm=TRUE)
argvartmean <- list(fun = "ns", knots = knotstmean)
cbt <- crossbasis(data[[Temp_measure]], lag = 2, argvar = argvartmean, arglag = list(knots = lagknots), group=data$salid1)
knot_50 <- crosspred(cbt, model.link="log",
                     coef = readRDS(paste0(path_for_pooled, "knot_50_pooled_main_results_coef.rds")),
                     vcov = readRDS(paste0(path_for_pooled, "knot_50_pooled_main_results_vcov.rds")),
                     cum=TRUE,cen=quan01,by=0.1)



# main, knots_10_75_90, knots_10_50_90,knots_25_50_75, knots_10_90, knots_25_75, knots_33_67,knot_50


# Define panel labels for each knot specification
panel_labels <- c(
  "Knots: 10, 50, 90", "Knots: 25, 50, 75", 
  "Knots: 10, 90", "Knots: 25, 75", 
  "Knots: 33, 67", "Knots: 50"
)


png("/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Submission-LTPH/figa1.png", width = 5, height = 6, units = "in", res = 300)  # 8x8 inches, 300 DPI

# Setup multi-panel layout
par(mfrow = c(3, 2)) 
col = 'gray'
  par(mar = c(4, 5, 4, 1.5), las = 1, mgp = c(3, 1, 0))
  
  # List of models to iterate over
  models <- list(knots_10_50_90, knots_25_50_75, 
                 knots_10_90, knots_25_75, knots_33_67, knot_50)
  
  # Loop through models and plot
  for (i in seq_along(models)) {
    plot(models[[i]], "overall", lag = 0, cumul = TRUE, ylim = c(0.9, 1.4), 
         ylab = "RR", col = col, lwd = 1.5, ci = 'area', lty = 1, 
         xlab = "Temperature (%tile)", ci.arg = list(col = alpha(col, 0.3)))
    
    ind2 <- models[[i]]$predvar >= quan01 & models[[i]]$predvar <= quan99
    lines(models[[i]]$predvar[ind2],models[[i]]$allRRfit[ind2],col='firebrick3',lwd=1.5)    
  
    # Add vertical reference lines
    abline(v = c(quan01, quan25, medT, quan75, quan99), lty = 4, col = "gray")
    
    # Add mtext label indicating knot specification
    mtext(panel_labels[i], side = 3, line = 1, cex = 1.2, font = 2) 
  }
  
  
dev.off()