# Modeling (non-imputed; "original model")
pred_road <- crosspred(cbt, model_road, cum=TRUE,cen=medT, by=0.1)

coef = coef(pred_road)
vcov = vcov(pred_road)

# Modeling ("reduced original model")
pred_road_reduce <- crossreduce(cbt, model_road,cen=medT, by=0.1)
coef_red = coef(pred_road_reduce)
vcov_red = vcov(pred_road_reduce)

# Reconstruct (reduced coefsl")
## set a cb with lag=0
cbt_reduce <- crossbasis(data[[Temp_measure]], lag = 0, argvar = argvartmean, group=data$salid1)
## reconstruct with reduced coef
pred_road_reconstruct <-  crosspred(cbt_reduce, model.link="log",
                                    coef = coef_red,
                                    vcov = vcov_red,cen=medT, by=0.1)
coef_red_reconstruct = coef(pred_road_reconstruct)
vcov_red_reconstruct = vcov(pred_road_reconstruct)


# Pooled Reconstruct (12 coefs)
pred_road_reconstruct_pooled <- crosspred(cbt, model.link="log",
                                          coef = readRDS(paste0(path_for_pooled, "pooled_main_results_coef.rds")),
                                          vcov = readRDS(paste0(path_for_pooled, "pooled_main_results_vcov.rds")),
                                         cum=TRUE,cen=medT)

# Pooled Reduced-Reconstruct (from 03_imputed: use reduced coef/vcov from pooling 100 imputations)
pred_road_red_pooled <- crosspred(cbt_reduce, model.link="log",
                                  coef = readRDS(paste0(path_for_pooled, "pooled_main_reduced_results_coef.rds")),
                                  vcov = readRDS(paste0(path_for_pooled, "pooled_main_reduced_results_vcov.rds")),
                                          cen=medT, by=0.1)
coef_red_pooled <- coef(pred_road_red_pooled)
vcov_red_pooled <- vcov(pred_road_red_pooled)

# Pooled Reconstruct (12 coefs) celcius
pred_road_reconstruct_celcius_pooled <- crosspred(cbt, model.link="log",
                                          coef = readRDS(paste0(path_for_pooled, "pooled_main_celcius_results_coef.rds")),
                                          vcov = readRDS(paste0(path_for_pooled, "pooled_main_celcius_results_vcov.rds")),
                                          cum=TRUE,cen=medT_celcius)



col='rosybrown'
par(mfrow = c(5, 1))
par(mar=c(4,5,1,0.5), las=1, mgp=c(2.5,1,0))
plot(pred_road, "overall", lag=0, cumul=TRUE, ylim=c(0.95,1.2), ylab="RR", col=col, lwd=1.5,ci='area',lty=1,
     xlab="Original", ci.arg=list(col=alpha(col, 0.3)))

plot(pred_road_reconstruct_celcius_pooled, "overall", lag=0, cumul=TRUE, xlim=c(0, 40),ylim=c(0.95,1.2), ylab="RR", col=col, lwd=1.5,ci='area',lty=1,
     xlab="Original", ci.arg=list(col=alpha(col, 0.3)))

plot(pred_road_reduce, ylab="RR", col=col, ylim=c(0.95,1.2), lwd=1.5,ci='area',lty=1,
     xlab="Reduced", ci.arg=list(col=alpha(col, 0.3)))

plot(pred_road_reconstruct_pooled, "overall", lag=0, cumul=TRUE, ylim=c(0.95,1.2), ylab="RR", col=col, lwd=1.5,ci='area',lty=1,
     xlab="Pooled (12 coef)", ci.arg=list(col=alpha(col, 0.3)))

plot(pred_road_red_pooled, ylab="RR", col=col, ylim=c(0.95,1.2), lwd=1.5,ci='area',lty=1,
     xlab="Pooled (4 coef)", ci.arg=list(col=alpha(col, 0.3)))




