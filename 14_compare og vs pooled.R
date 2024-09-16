# Modeling (non-imputed; "original model")
pred_road <- crosspred(cbt, model_road, cum=TRUE)

coef = coef(pred_road)
vcov = vcov(pred_road)

# Pooled Reconstruct (from 03_imputed: use coef/vcov from pooling 100 imputations)
coef_pooled <- pooled_results$pooled_coef
vcov_pooled <- pooled_results$pooled_vcov

pred_road_reconstruct_pooled <- crosspred(cbt, model.link="log",
                                         coef = coef_pooled,
                                         vcov = vcov_pooled, cum=TRUE)

# Modeling ("reduced original model")
pred_road_reduce <- crossreduce(cbt, model_road)
coef_red = coef(pred_road_reduce)
vcov_red = vcov(pred_road_reduce)

# Reconstruct ("reconstructed original model", same as "reduced original model")
## set a cb with lag=0
cbt_reduce <- crossbasis(data[[Temp_measure]], lag = 0, argvar = argvartmean, group=data$salid1)
## reconstruct with reduced coef
pred_road_reconstruct <-  crosspred(cbt_reduce, model.link="log",
                                    coef = coef_red,
                                    vcov = vcov_red)
coef_red_reconstruct = coef(pred_road_reconstruct)
vcov_red_reconstruct = vcov(pred_road_reconstruct)

# Pooled Reduced-Reconstruct (from 03_imputed: use reduced coef/vcov from pooling 100 imputations)
coef_red_pooled <- pooled_results_reduced$pooled_coef
vcov_red_pooled <- pooled_results_reduced$pooled_vcov

pred_road_red_pooled <- crosspred(cbt_reduce, model.link="log",
                                          coef = coef_red_pooled,
                                          vcov = vcov_red_pooled)


par(mfrow = c(4, 1))
par(mar=c(4,5,1,0.5), las=1, mgp=c(2.5,1,0))
plot(pred_road, "overall", lag=0, cumul=TRUE, ylim=c(0.8,1.3), ylab="RR", col="brown4", lwd=1.5,ci='bars',lty=1,
     xlab="Original", ci.arg=list(col=alpha("black", 0.03)))

plot(pred_road_reconstruct_pooled, "overall", lag=0, cumul=TRUE, ylim=c(0.8,1.3), ylab="RR", col="brown4", lwd=1.5,ci='bars',lty=1,
     xlab="Pooled (12 coef)", ci.arg=list(col=alpha("black", 0.03)))

plot(pred_road_reduce, ylab="RR", col="brown4", ylim=c(0.8,1.3), lwd=1.5,ci='bars',lty=1,
     xlab="Reduced", ci.arg=list(col=alpha("black", 0.03)))

plot(pred_road_red_pooled, ylab="RR", col="brown4", ylim=c(0.8,1.3), lwd=1.5,ci='bars',lty=1,
     xlab="Pooled (4 coef)", ci.arg=list(col=alpha("black", 0.03)))

