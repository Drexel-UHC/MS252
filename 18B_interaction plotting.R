################################################################################################################################################################
png("/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Submission-LTPH/figa2.png", width = 6, height = 3, units = "in", res = 300)  # 8x8 inches, 300 DPI

# Reconstruction by modifier
modifier = 'mean' #mean #std #BECADSTTLGAVGL1AD #BECURBAVGTRAFTIMEL1AD
coef_int1_pooled <- readRDS(paste0(path_for_pooled, "pooled_int1_", modifier, "_results_coef.rds"))
vcov_int1_pooled <- readRDS(paste0(path_for_pooled, "pooled_int1_", modifier, "_results_vcov.rds"))
coef_int2_pooled <- readRDS(paste0(path_for_pooled, "pooled_int2_", modifier, "_results_coef.rds"))
vcov_int2_pooled <- readRDS(paste0(path_for_pooled, "pooled_int2_", modifier, "_results_vcov.rds"))

pred_int1_pooled <- crosspred(cbt, model.link="log", coef = coef_int1_pooled,vcov = vcov_int1_pooled, cum=TRUE, cen=quan01, by=0.1)
pred_int2_pooled <- crosspred(cbt, model.link="log",coef = coef_int2_pooled, vcov = vcov_int2_pooled, cum=TRUE, cen=quan01, by=0.1)



#Plotting
col <- c("steelblue4", "firebrick3")
par(mfrow = c(1, 2))
par(mar=c(4,3.5,2,1), las=1, mgp=c(2.5,1,0))
plot(pred_int1_pooled, "overall", lag=0, cumul=TRUE, ylim=c(0.95,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=1,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.01)))
lines(pred_int2_pooled, "overall", lag=0, cumul=TRUE, ci="bars",col=col[2], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.01)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")
legend("topleft", c("lower","higher"), lty=1, lwd=1, col=col,
       bty="n", inset=0.01, y.intersp=0.9, cex=1.0)
mtext("a", side = 3, adj = 0.5, line = 0.5, cex = 1.2)



# Reconstruction by modifier
modifier = 'std' #mean #std #BECADSTTLGAVGL1AD #BECURBAVGTRAFTIMEL1AD
coef_int1_pooled <- readRDS(paste0(path_for_pooled, "pooled_int1_", modifier, "_results_coef.rds"))
vcov_int1_pooled <- readRDS(paste0(path_for_pooled, "pooled_int1_", modifier, "_results_vcov.rds"))
coef_int2_pooled <- readRDS(paste0(path_for_pooled, "pooled_int2_", modifier, "_results_coef.rds"))
vcov_int2_pooled <- readRDS(paste0(path_for_pooled, "pooled_int2_", modifier, "_results_vcov.rds"))

pred_int1_pooled <- crosspred(cbt, model.link="log", coef = coef_int1_pooled,vcov = vcov_int1_pooled, cum=TRUE, cen=quan01, by=0.1)
pred_int2_pooled <- crosspred(cbt, model.link="log",coef = coef_int2_pooled, vcov = vcov_int2_pooled, cum=TRUE, cen=quan01, by=0.1)

#Plotting
col <- c("steelblue4", "firebrick3")
plot(pred_int1_pooled, "overall", lag=0, cumul=TRUE, ylim=c(0.95,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=1,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.01)))
lines(pred_int2_pooled, "overall", lag=0, cumul=TRUE, ci="bars",col=col[2], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.01)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")
legend("topleft", c("lower","higher"), lty=1, lwd=1, col=col,
       bty="n", inset=0.01, y.intersp=0.9, cex=1.0)
mtext("b", side = 3, adj = 0.5, line = 0.5, cex = 1.2)
dev.off()
################################################################################################################################################################
################################################################################################################
png("/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Submission-LTPH/fig4.png", width = 6, height = 3, units = "in", res = 300)  # 8x8 inches, 300 DPI
# Reconstruction by modifier
modifier = 'BECADSTTLGAVGL1AD' #mean #std #BECADSTTLGAVGL1AD #BECURBAVGTRAFTIMEL1AD
coef_int1_pooled <- readRDS(paste0(path_for_pooled, "pooled_int1_", modifier, "_results_coef.rds"))
vcov_int1_pooled <- readRDS(paste0(path_for_pooled, "pooled_int1_", modifier, "_results_vcov.rds"))
coef_int2_pooled <- readRDS(paste0(path_for_pooled, "pooled_int2_", modifier, "_results_coef.rds"))
vcov_int2_pooled <- readRDS(paste0(path_for_pooled, "pooled_int2_", modifier, "_results_vcov.rds"))

pred_int1_pooled <- crosspred(cbt, model.link="log", coef = coef_int1_pooled,vcov = vcov_int1_pooled, cum=TRUE, cen=quan01, by=0.1)
pred_int2_pooled <- crosspred(cbt, model.link="log",coef = coef_int2_pooled, vcov = vcov_int2_pooled, cum=TRUE, cen=quan01, by=0.1)



#Plotting
col <- c("steelblue4", "firebrick3")
par(mfrow = c(1, 2))
par(mar=c(4,3.5,2,1), las=1, mgp=c(2.5,1,0))
plot(pred_int1_pooled, "overall", lag=0, cumul=TRUE, ylim=c(0.95,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=1,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.01)))
lines(pred_int2_pooled, "overall", lag=0, cumul=TRUE, ci="bars",col=col[2], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.01)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")
legend("topleft", c("lower","higher"), lty=1, lwd=1, col=col,
       bty="n", inset=0.01, y.intersp=0.9, cex=1.0)
mtext("a", side = 3, adj = 0.5, line = 0.5, cex = 1.2)



# Reconstruction by modifier
modifier = 'BECURBAVGTRAFTIMEL1AD' #mean #std #BECADSTTLGAVGL1AD #BECURBAVGTRAFTIMEL1AD
coef_int1_pooled <- readRDS(paste0(path_for_pooled, "pooled_int1_", modifier, "_results_coef.rds"))
vcov_int1_pooled <- readRDS(paste0(path_for_pooled, "pooled_int1_", modifier, "_results_vcov.rds"))
coef_int2_pooled <- readRDS(paste0(path_for_pooled, "pooled_int2_", modifier, "_results_coef.rds"))
vcov_int2_pooled <- readRDS(paste0(path_for_pooled, "pooled_int2_", modifier, "_results_vcov.rds"))

pred_int1_pooled <- crosspred(cbt, model.link="log", coef = coef_int1_pooled,vcov = vcov_int1_pooled, cum=TRUE, cen=quan01, by=0.1)
pred_int2_pooled <- crosspred(cbt, model.link="log",coef = coef_int2_pooled, vcov = vcov_int2_pooled, cum=TRUE, cen=quan01, by=0.1)

#Plotting
col <- c("steelblue4", "firebrick3")
plot(pred_int1_pooled, "overall", lag=0, cumul=TRUE, ylim=c(0.95,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=1,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.01)))
lines(pred_int2_pooled, "overall", lag=0, cumul=TRUE, ci="bars",col=col[2], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.01)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")
legend("topleft", c("lower","higher"), lty=1, lwd=1, col=col,
       bty="n", inset=0.01, y.intersp=0.9, cex=1.0)
mtext("b", side = 3, adj = 0.5, line = 0.5, cex = 1.2)
dev.off()