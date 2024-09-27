library("tidyverse")
library("data.table")
library("dlnm"); library("gnm"); library("splines") # Modeling tools


# Load data ---------------------------------------------------------------

# Number of imputations to run (max 100)
n_imp <- 100

# Load data
path_data <- "C:/Users/dpw48/git/MS252/Data/"

all <- fread(paste0(path_data, "final_temp_cluster_subgroup_1.csv"))

# Select all temperature clusters for analysis
data <- all

# Select only specific temperature clusters for analysis
# data <- all[all[['cluster_ward_std_6']] %in% c(1)]

# Set up model ------------------------------------------------------------

# Specify temperature measure
Temp_measure <- "tmp_pw_percentile" # percentile
# Temp_measure <- "L1ADtemp_pw"     # C°
minT <- min(data[[Temp_measure]], na.rm=TRUE)
maxT <- max(data[[Temp_measure]], na.rm=TRUE)
medT <- median(data[[Temp_measure]], na.rm=TRUE)
quan01 <- quantile(data[[Temp_measure]], probs = 0.01, na.rm = TRUE)
quan25 <- quantile(data[[Temp_measure]], probs = 0.25, na.rm = TRUE)
quan75 <- quantile(data[[Temp_measure]], probs = 0.75, na.rm = TRUE)
quan99 <- quantile(data[[Temp_measure]], probs = 0.99, na.rm = TRUE)

# Specify knot number for lags
lagknots <- logknots(3, df = 3)

# Specify knot locations for temperature
knotstmean <- quantile(data[[Temp_measure]], c(10,75,90)/100, na.rm=T)
argvartmean <- list(fun="ns", knots=knotstmean)

# Create crossbasis for DLNM
cbt <- crossbasis(data[[Temp_measure]], lag = 2,
                  argvar = argvartmean, arglag = list(knots = lagknots))

# Model data --------------------------------------------------------------

fit_imputed <- 1:n_imp |> 
  map(\(i) {
    all <- fread(paste0(path_data, "final_temp_cluster_subgroup_", i, ".csv"))
    
    # Select all temperature clusters for analysis
    data <- all
    
    # Specify stratification scheme
    # stratify by city-year-month-dow
    data[, stratum:=factor(paste(salid1, year_month, dow, sep=":"))]
    
    # Keep all strata with nonzero number of road deaths
    data[,  keep:=sum(road)>0, by=stratum]
    
    # Fit model
    model_road <- gnm(road ~ cbt, eliminate = stratum, 
                      family = quasipoisson(), data = data, subset=keep)
    
    list(coef = coef(model_road),
         vcov = vcov(model_road))
  })

# Extract coefficients and covariance matrices
coef_imputed <- map(fit_imputed, \(x) x$coef)
vcov_imputed <- map(fit_imputed, \(x) x$vcov)

# Estimated imputed coefs
q_bar <- coef_imputed |> 
  unlist() |> 
  matrix(ncol = n_imp) |> 
  apply(1, mean)

W_bar <- Reduce('+', vcov_imputed)/n_imp

B <- coef_imputed |> 
  unlist() |> 
  matrix(ncol = n_imp) |> 
  t() |> 
  cov()

T_0 <- W_bar + (1+1/n_imp)*B
  

# Using DLNM functions with q_bar and T_0
pred_road_bar <- crosspred(cbt, coef = q_bar, vcov = T_0, model.link = "log",
                           at = seq(minT, maxT, by = .1),
                           cen = medT, bylag = 1, cumul = TRUE, ci.level = 0.95)

# Plot overal RR and temperature histogram
par(mfrow = c(2, 1))
par(mar=c(4,5,1,0.5), las=1, mgp=c(2.5,1,0))

plot(pred_road_bar, "overall", lag=0, cumul=TRUE,
     ylim=c(0.8,1.3), ylab="RR", col="brown4", lwd=1.5,ci='bars',lty=1,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha("black", 0.03)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")

ind1 <- pred_road_bar$predvar<=medT
ind2 <- pred_road_bar$predvar>=medT
lines(pred_road_bar$predvar[ind1],pred_road_bar$allRRfit[ind1],col=4,lwd=1.5)
lines(pred_road_bar$predvar[ind2],pred_road_bar$allRRfit[ind2],col=2,lwd=1.5)


hist(data[[Temp_measure]], breaks = 100, probability = TRUE, col = "grey60",
     main = '', xlab = "Temperature (%tile)", ylab = "Probability Density", yaxt = "n", ylim=c(0,0.25))




# Extra code --------------------------------------------------------------


# Create crossbasis for DLNM
cbt <- crossbasis(data[[Temp_measure]], lag = 2,
                  argvar = argvartmean, arglag = list(knots = lagknots))

# Specify stratification scheme
# stratify by city-year-month-dow
data[, stratum:=factor(paste(salid1, year_month, dow, sep=":"))]

# Keep all strata with nonzero number of road deaths
data[,  keep:=sum(road)>0, by=stratum]

# Fit model
model_road <- gnm(road ~ cbt, eliminate = stratum, 
                  family = quasipoisson(), data = data, subset=keep)

# Estimate association via the DLNM basis
pred_road <- crosspred(cbt, model_road,
                       at = seq(minT, maxT, by = .1),
                       cen = medT, bylag = 1, cumul = TRUE, ci.level = 0.95)

# Plot overal RR and temperature histogram
par(mfrow = c(2, 1))
par(mar=c(4,5,1,0.5), las=1, mgp=c(2.5,1,0))

plot(pred_road, "overall", lag=0, cumul=TRUE,
     ylim=c(0.8,1.3), ylab="RR", col="brown4", lwd=1.5,ci='bars',lty=1,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha("black", 0.03)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")

ind1 <- pred_road$predvar<=medT
ind2 <- pred_road$predvar>=medT
lines(pred_road$predvar[ind1],pred_road$allRRfit[ind1],col=4,lwd=1.5)
lines(pred_road$predvar[ind2],pred_road$allRRfit[ind2],col=2,lwd=1.5)


hist(data[[Temp_measure]], breaks = 100, probability = TRUE, col = "grey60",
     main = '', xlab = "Temperature (%tile)", ylab = "Probability Density", yaxt = "n", ylim=c(0,0.25))


# Pool models -------------------------------------------------------------

# For now work to pool these three models
model_road
model_road2
model_road3


q1 <- coef(model_road)
q2 <- coef(model_road2)
q3 <- coef(model_road3)
q_bar <- apply(matrix(c(q1, q2, q3), ncol = 3), 1, mean)

W1 <- vcov(model_road)
W2 <- vcov(model_road2)
W3 <- vcov(model_road3)
W_bar <- (W1 + W2 + W3)/3

B <- cov(t(matrix(c(q1, q2, q3), ncol = 3)))

T_0 <- W_bar + (1+1/3)*B

# So q_bar is our mean and T_0 is its associated covariance matrix

# Using DLNM functions with q_bar and T_0

pred_road_bar <- crosspred(cbt, coef = q_bar, vcov = T_0, model.link = "log",
                           at = seq(minT, maxT, by = .1),
                           cen = medT, bylag = 1, cumul = TRUE, ci.level = 0.95)

pred_road_bar$model.class <- c("gnm", "glm", "lm")

# Plot overal RR and temperature histogram
par(mfrow = c(2, 1))
par(mar=c(4,5,1,0.5), las=1, mgp=c(2.5,1,0))

plot(pred_road_bar, "overall", lag=0, cumul=TRUE,
     ylim=c(0.8,1.3), ylab="RR", col="brown4", lwd=1.5,ci='bars',lty=1,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha("black", 0.03)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")

ind1 <- pred_road_bar$predvar<=medT
ind2 <- pred_road_bar$predvar>=medT
lines(pred_road_bar$predvar[ind1],pred_road_bar$allRRfit[ind1],col=4,lwd=1.5)
lines(pred_road_bar$predvar[ind2],pred_road_bar$allRRfit[ind2],col=2,lwd=1.5)


hist(data[[Temp_measure]], breaks = 100, probability = TRUE, col = "grey60",
     main = '', xlab = "Temperature (%tile)", ylab = "Probability Density", yaxt = "n", ylim=c(0,0.25))

# Usama input -------------------------------------------------------------

# Estimate = average of estimates
# Standard error = sqrt(var(estimates across iterations) + mean(se^2 across iterations))
# Use estimates + SEs to make tables/figures


# restricted cubic splines (natural splines) with 3 knots (matching our quartiles above)
knots_loc<-c(0.25, 0.50, 0.75)
basis<-onebasis(.x$SEI, fun="ns", knots=quantile(.x$SEI, probs=knots_loc))
f<-as.formula(paste0("deaths~age_cat+basis+offset(log(pop))+(1|ISO2)"))
m<-glmmTMB(f, data=.x, family = nbinom2(), verbose = F)

# since the dlnm package does not have a prediction function for glmmTMB, i'me xtracting the coefficients manually
coef<-fixef(m)$cond
vcov<-vcov(m)$cond

# restricting to the spline basis
coef<-coef[grep("basis", names(coef))]
vcov<-vcov[grep("basis", rownames(vcov)),grep("basis", rownames(vcov))]

# predicting with reference at 0
pred<-crosspred(basis, coef=coef, vcov=vcov,
                cen=0)

result<-data.frame(rr=exp(as.numeric(pred$matfit)),
                   lci=exp(as.numeric(pred$matlow)),
                   uci=exp(as.numeric(pred$mathigh)),
                   SEI=as.numeric(pred$predvar))