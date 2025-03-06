path_for_pooled <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/imputed/Pooled results/"

# Reconstruction
pred_male_pooled <- crosspred(cbt, model.link="log",
                              coef = readRDS(paste0(path_for_pooled, "pooled_male_results_coef.rds")),
                              vcov = readRDS(paste0(path_for_pooled, "pooled_male_results_vcov.rds")),
                              cum=TRUE, cen=quan01, by=0.1)

pred_female_pooled <- crosspred(cbt, model.link="log",
                                coef = readRDS(paste0(path_for_pooled, "pooled_female_results_coef.rds")),
                                vcov = readRDS(paste0(path_for_pooled, "pooled_female_results_vcov.rds")),
                                cum=TRUE, cen=quan01, by=0.1)


############################################################################################################################################
png("/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Submission-LTPH/fig3.png", width = 7, height = 2.5, units = "in", res = 300)  # 8x8 inches, 300 DPI

# Plotting
par(mfrow = c(1, 3))
par(mar=c(4,4,2,2), las=1, mgp=c(2.5,1,0))
col <- brewer.pal(n = 11, name = "Blues")[c(9,7)]
plot(pred_male_pooled, "overall", lag=0, cumul=TRUE, ylim=c(0.95,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=1,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.01)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")

lines(pred_female_pooled, "overall", lag=0, cumul=TRUE, ci="bars",col=col[2], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.01)))

legend("topleft", c("male","female"), lty=1, lwd=1, col=col,
       bty="n", inset=0, y.intersp=1.0, cex=1)
mtext("a", side = 3, adj = 0.5, line = 0.5, cex = 0.95)

###########
# Reconstruction
age_groups <- c("age1", "age2", "age3", "age4", "age5")
pred_pooled_list <- setNames(lapply(age_groups, function(age_group) {
  crosspred(cbt, model.link = "log",
            coef = readRDS(paste0(path_for_pooled, "pooled_", age_group, "_results_coef.rds")),
            vcov = readRDS(paste0(path_for_pooled, "pooled_", age_group, "_results_vcov.rds")),
            cum = TRUE, cen = quan01, by = 0.1)
}), age_groups)

# Plotting
col <- c('darkorange2','darkorange4', "steelblue2","steelblue4",'cadetblue4')
col <- heat.colors(10)[3:8]
col <- brewer.pal(n = 21, name = "OrRd")[c(4,5,7,8,9)]

plot(pred_pooled_list[['age1']], "overall", lag=0, cumul=TRUE, ylim=c(0.95,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=1,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.01)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")

lines(pred_pooled_list[['age2']], "overall", lag=0, cumul=TRUE, ci="bars",col=col[2], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.01)))
lines(pred_pooled_list[['age3']], "overall", lag=0, cumul=TRUE, ci="bars",col=col[3], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[3], 0.01)))
lines(pred_pooled_list[['age4']], "overall", lag=0, cumul=TRUE, ci="bars",col=col[4], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[4], 0.01)))
lines(pred_pooled_list[['age5']], "overall", lag=0, cumul=TRUE, ci="bars",col=col[5], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[5], 0.01)))

legend("topleft", c("<=9","10-19","20-34","35-64",">=65"), lty=1, lwd=2, col=col,
       bty="n", inset=0, y.intersp=1.0, cex=1)
mtext("b", side = 3, adj = 0.5, line = 0.5, cex = 0.95)


#####
# Reconstruction
modes <- c("vehicle", "motorcycle", "bicycle", "ped")
pred_pooled_list <- setNames(lapply(modes, function(mode) {
  crosspred(cbt, model.link = "log",
            coef = readRDS(paste0(path_for_pooled, "pooled_", mode, "_results_coef.rds")),
            vcov = readRDS(paste0(path_for_pooled, "pooled_", mode, "_results_vcov.rds")),
            cum = TRUE, cen = quan01, by = 0.1)
}), modes)

# Plotting
col <- c("firebrick2",'darkorange2','cadetblue4', 'steelblue3')
col <- brewer.pal(n = 11, name = "RdYlBu")[c(1,3,9,11)]
plot(pred_pooled_list[['motorcycle']], "overall", lag=0, cumul=TRUE, ylim=c(0.95,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=1,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.01)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")

lines(pred_pooled_list[['bicycle']], "overall", lag=0, cumul=TRUE, ci="bars",col=col[2], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.01)))
lines(pred_pooled_list[['vehicle']], "overall", lag=0, cumul=TRUE, ci="bars",col=col[3], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[3], 0.01)))
lines(pred_pooled_list[['ped']], "overall", lag=0, cumul=TRUE, ci="bars",col=col[4], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[4], 0.01)))
legend("topleft", c("motorcycle","bicycle","vehicle","ped"), lty=1, lwd=1, col=col,
       bty="n", inset=0, y.intersp=1.0, cex=1)
mtext("c", side = 3, adj = 0.5, line = 0.5, cex = 0.95)

dev.off()
############################################################################################################################################