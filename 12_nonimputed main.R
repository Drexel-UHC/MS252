library(haven)     # For reading .sas7bdat files
library(dplyr)     # For data manipulation
library(lubridate) # For date manipulation
library(dlnm)
library(splines)
library(RColorBrewer)
library(tsModel)
library(data.table) # HANDLE LARGE DATASETS
library(dlnm) ; library(gnm) ; library(splines) # MODELLING TOOLS
library(sf) ; library(terra) # HANDLE SPATIAL DATA
library(exactextractr) # FAST EXTRACTION OF AREA-WEIGHTED RASTER CELLS
library(dplyr) ; library(tidyr) # DATA MANAGEMENT TOOLS
library(ggplot2) ; library(patchwork) # PLOTTING TOOLS
library(gnm)
library(stargazer)
library(broom)

fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}

non_imputed <- read.csv("/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS 252 data_July 2024/Python Processed Data/final_temp_cluster_subgroup.csv")
non_imputed <- as.data.table(non_imputed)
summary(non_imputed)

# Model set-up (using non-imputed data)
data <- non_imputed %>% filter(!grepl("^2020", `year`))

Temp_measure <- "tmp_pw_percentile"
#Temp_measure <- "L1ADtemp_pw"
minT <- min(data[[Temp_measure]], na.rm=TRUE)
maxT <- max(data[[Temp_measure]], na.rm=TRUE)
medT <- median(data[[Temp_measure]], na.rm=TRUE)
quan01 <- quantile(data[[Temp_measure]], probs = 0.01, na.rm = TRUE)
quan25 <- quantile(data[[Temp_measure]], probs = 0.25, na.rm = TRUE)
quan75 <- quantile(data[[Temp_measure]], probs = 0.75, na.rm = TRUE)
quan99 <- quantile(data[[Temp_measure]], probs = 0.99, na.rm = TRUE)

# Modeling (non-imputed; "original model")
lagknots <- logknots(3, df = 3)
knotstmean <- quantile(data[[Temp_measure]], c(10,75,90)/100, na.rm=T)
argvartmean <- list(fun="ns", knots=knotstmean)
cbt <- crossbasis(data[[Temp_measure]], lag = 2, argvar = argvartmean, arglag = list(knots = lagknots), group=data$salid1)

data[, stratum:=factor(paste(salid1, year.month, dow, sep=":"))]
data[,  keep:=sum(median_road_round)>0, by=stratum]
model_road <- gnm(median_road_round ~ cbt, eliminate = stratum, 
                  family = quasipoisson(), data = data, subset=keep)
pred_road <- crosspred(cbt, model_road, cum=TRUE)

coef = coef(pred_road)
vcov = vcov(pred_road)
