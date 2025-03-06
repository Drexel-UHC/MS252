#EDF, or AF
attrdl <- function(x,basis,cases,model=NULL,coef=NULL,vcov=NULL,model.link=NULL,
                   type="af",dir="back",tot=TRUE,cen,range=NULL,sim=FALSE,nsim=5000) {
  ################################################################################
  #
  # CHECK VERSION OF THE DLNM PACKAGE
  if(packageVersion("dlnm")<"2.2.0") 
    stop("update dlnm package to version >= 2.2.0")
  #
  # EXTRACT NAME AND CHECK type AND dir
  name <- deparse(substitute(basis))
  type <- match.arg(type,c("an","af"))
  dir <- match.arg(dir,c("back","forw"))
  #
  # DEFINE CENTERING
  if(missing(cen) && is.null(cen <- attr(basis,"argvar")$cen))
    stop("'cen' must be provided")
  if(!is.numeric(cen) && length(cen)>1L) stop("'cen' must be a numeric scalar")
  attributes(basis)$argvar$cen <- NULL
  #  
  # SELECT RANGE (FORCE TO CENTERING VALUE OTHERWISE, MEANING NULL RISK)
  if(!is.null(range)) x[x<range[1]|x>range[2]] <- cen
  #
  # COMPUTE THE MATRIX OF
  #   - LAGGED EXPOSURES IF dir="back"
  #   - CONSTANT EXPOSURES ALONG LAGS IF dir="forw"
  lag <- attr(basis,"lag")
  if(NCOL(x)==1L) {
    at <- if(dir=="back") tsModel:::Lag(x,seq(lag[1],lag[2])) else 
      matrix(rep(x,diff(lag)+1),length(x))
  } else {
    if(dir=="forw") stop("'x' must be a vector when dir='forw'")
    if(ncol(at <- x)!=diff(lag)+1) 
      stop("dimension of 'x' not compatible with 'basis'")
  }
  #
  # NUMBER USED FOR THE CONTRIBUTION AT EACH TIME IN FORWARD TYPE
  #   - IF cases PROVIDED AS A MATRIX, TAKE THE ROW AVERAGE
  #   - IF PROVIDED AS A TIME SERIES, COMPUTE THE FORWARD MOVING AVERAGE
  #   - THIS EXCLUDES MISSING ACCORDINGLY
  # ALSO COMPUTE THE DENOMINATOR TO BE USED BELOW
  if(NROW(cases)!=NROW(at)) stop("'x' and 'cases' not consistent")
  if(NCOL(cases)>1L) {
    if(dir=="back") stop("'cases' must be a vector if dir='back'")
    if(ncol(cases)!=diff(lag)+1) stop("dimension of 'cases' not compatible")
    den <- sum(rowMeans(cases,na.rm=TRUE),na.rm=TRUE)
    cases <- rowMeans(cases)
  } else {
    den <- sum(cases,na.rm=TRUE) 
    if(dir=="forw") 
      cases <- rowMeans(as.matrix(tsModel:::Lag(cases,-seq(lag[1],lag[2]))))
  }
  #
  ################################################################################
  #
  # EXTRACT COEF AND VCOV IF MODEL IS PROVIDED
  if(!is.null(model)) {
    cond <- paste0(name,"[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}")
    if(ncol(basis)==1L) cond <- name
    model.class <- class(model)
    coef <- dlnm:::getcoef(model,model.class)
    ind <- grep(cond,names(coef))
    coef <- coef[ind]
    vcov <- dlnm:::getvcov(model,model.class)[ind,ind,drop=FALSE]
    model.link <- dlnm:::getlink(model,model.class)
    if(!model.link %in% c("log","logit"))
      stop("'model' must have a log or logit link function")
  }
  #
  # IF REDUCED ESTIMATES ARE PROVIDED
  typebasis <- ifelse(length(coef)!=ncol(basis),"one","cb")
  #
  ################################################################################
  #
  # PREPARE THE ARGUMENTS FOR TH BASIS TRANSFORMATION
  predvar <- if(typebasis=="one") x else seq(NROW(at))
  predlag <- if(typebasis=="one") 0 else dlnm:::seqlag(lag)
  #  
  # CREATE THE MATRIX OF TRANSFORMED CENTRED VARIABLES (DEPENDENT ON typebasis)
  if(typebasis=="cb") {
    Xpred <- dlnm:::mkXpred(typebasis,basis,at,predvar,predlag,cen)
    Xpredall <- 0
    for (i in seq(length(predlag))) {
      ind <- seq(length(predvar))+length(predvar)*(i-1)
      Xpredall <- Xpredall + Xpred[ind,,drop=FALSE]
    }
  } else {
    basis <- do.call(onebasis,c(list(x=x),attr(basis,"argvar")))
    Xpredall <- dlnm:::mkXpred(typebasis,basis,x,predvar,predlag,cen)
  }
  #  
  # CHECK DIMENSIONS  
  if(length(coef)!=ncol(Xpredall))
    stop("arguments 'basis' do not match 'model' or 'coef'-'vcov'")
  if(any(dim(vcov)!=c(length(coef),length(coef)))) 
    stop("arguments 'coef' and 'vcov' do no match")
  if(typebasis=="one" && dir=="back")
    stop("only dir='forw' allowed for reduced estimates")
  #
  ################################################################################
  #
  # COMPUTE AF AND AN 
  af <- 1-exp(-drop(as.matrix(Xpredall%*%coef)))
  an <- af*cases
  #
  # TOTAL
  #   - SELECT NON-MISSING OBS CONTRIBUTING TO COMPUTATION
  #   - DERIVE TOTAL AF
  #   - COMPUTE TOTAL AN WITH ADJUSTED DENOMINATOR (OBSERVED TOTAL NUMBER)
  if(tot) {
    isna <- is.na(an)
    af <- sum(an[!isna])/sum(cases[!isna])
    an <- af*den
  }
  #
  ################################################################################
  #
  # EMPIRICAL CONFIDENCE INTERVALS
  if(!tot && sim) {
    sim <- FALSE
    warning("simulation samples only returned for tot=T")
  }
  if(sim) {
    # SAMPLE COEF
    k <- length(coef)
    eigen <- eigen(vcov)
    X <- matrix(rnorm(length(coef)*nsim),nsim)
    coefsim <- coef + eigen$vectors %*% diag(sqrt(eigen$values),k) %*% t(X)
    # RUN THE LOOP
    # pre_afsim <- (1 - exp(- Xpredall %*% coefsim)) * cases # a matrix
    # afsim <- colSums(pre_afsim,na.rm=TRUE) / sum(cases[!isna],na.rm=TRUE)
    afsim <- apply(coefsim,2, function(coefi) {
      ani <- (1-exp(-drop(Xpredall%*%coefi)))*cases
      sum(ani[!is.na(ani)])/sum(cases[!is.na(ani)])
    })
    ansim <- afsim*den
  }
  #
  ################################################################################
  #
  res <- if(sim) {
    if(type=="an") ansim else afsim
  } else {
    if(type=="an") an else af    
  }
  #
  return(res)
}

edf <- attrdl(data$tmp_pw_percentile,cbt,data$median_road_round, dir="forw",cen=minT,
       coef = readRDS(paste0(path_for_pooled, "pooled_main_results_coef.rds")),
       vcov = readRDS(paste0(path_for_pooled, "pooled_main_results_vcov.rds")),
       range=c(95,100))

edf_ci <- quantile(attrdl(data$tmp_pw_percentile,cbt,data$median_road_round, dir="forw",cen=minT,
                coef = readRDS(paste0(path_for_pooled, "pooled_main_results_coef.rds")),
                vcov = readRDS(paste0(path_for_pooled, "pooled_main_results_vcov.rds")),
                range=c(95, 100), sim=T,nsim=1000), c(0.025,0.975))

# temperature cluster
sapply(1:6, function(n) length(unique(data[data$cluster_ward_std_6 == n, ]$salid1))) 
# 1) mild N=21; 2) warm/con, N=82; 3) hot/con N=103; 4) warm/dis N=30; 5) hot/dis N=19; 6) cool N=17

# Initialize lists to store results
edf_cluster <- list()
edf_ci_cluster <- list()

# Loop through cluster numbers 1 to 6
for (n in 1:6) {
  # Subset data dynamically
  subset_data <- data[data$cluster_ward_std_6 == n, ]
  
  # Compute EDF
  edf_cluster[[n]] <- attrdl(
    subset_data$tmp_pw_percentile, cbt, subset_data$median_road_round, dir = "forw", cen = minT,
    coef = readRDS(paste0(path_for_pooled, "pooled_results_cluster", n, "_coef.rds")),
    vcov = readRDS(paste0(path_for_pooled, "pooled_results_cluster", n, "_vcov.rds")),
    range = c(95, 100)
  )
  
  # Compute EDF confidence interval
  edf_ci_cluster[[n]] <- quantile(
    attrdl(
      subset_data$tmp_pw_percentile, cbt, subset_data$median_road_round, dir = "forw", cen = minT,
      coef = readRDS(paste0(path_for_pooled, "pooled_results_cluster", n, "_coef.rds")),
      vcov = readRDS(paste0(path_for_pooled, "pooled_results_cluster", n, "_vcov.rds")),
      range = c(95, 100), sim = TRUE, nsim = 1000
    ), 
    c(0.025, 0.975)
  )
}

edf_cluster
edf_ci_cluster

# sex
edf_male <- attrdl(data$tmp_pw_percentile,cbt,data$male_deaths, dir="forw",cen=minT,
                   coef = readRDS(paste0(path_for_pooled, "pooled_male_results_coef.rds")),
                   vcov = readRDS(paste0(path_for_pooled, "pooled_male_results_vcov.rds")),
                   range=c(95,100))

edf_ci_male <- quantile(attrdl(data$tmp_pw_percentile,cbt,data$male_deaths, dir="forw",cen=minT,
                               coef = readRDS(paste0(path_for_pooled, "pooled_male_results_coef.rds")),
                               vcov = readRDS(paste0(path_for_pooled, "pooled_male_results_vcov.rds")),
                               range=c(95, 100), sim=T,nsim=1000), c(0.025,0.975))

edf_female <- attrdl(data$tmp_pw_percentile,cbt,data$female_deaths, dir="forw",cen=minT,
                   coef = readRDS(paste0(path_for_pooled, "pooled_female_results_coef.rds")),
                   vcov = readRDS(paste0(path_for_pooled, "pooled_female_results_vcov.rds")),
                   range=c(95,100))

edf_ci_female <- quantile(attrdl(data$tmp_pw_percentile,cbt,data$female_deaths, dir="forw",cen=minT,
                               coef = readRDS(paste0(path_for_pooled, "pooled_female_results_coef.rds")),
                               vcov = readRDS(paste0(path_for_pooled, "pooled_female_results_vcov.rds")),
                               range=c(95, 100), sim=T,nsim=1000), c(0.025,0.975))
edf_female
edf_ci_female
# age
edf_age1 <- attrdl(data$tmp_pw_percentile,cbt,data$deaths_under_9, dir="forw",cen=minT,
                         coef = readRDS(paste0(path_for_pooled, "pooled_age1_results_coef.rds")),
                         vcov = readRDS(paste0(path_for_pooled, "pooled_age1_results_vcov.rds")),
                         range=c(95,100))

edf_ci_age1 <- quantile(attrdl(data$tmp_pw_percentile,cbt,data$deaths_under_9, dir="forw",cen=minT,
                                     coef = readRDS(paste0(path_for_pooled, "pooled_age1_results_coef.rds")),
                                     vcov = readRDS(paste0(path_for_pooled, "pooled_age1_results_vcov.rds")),
                                     range=c(95, 100), sim=T,nsim=1000), c(0.025,0.975))

edf_age2 <- attrdl(data$tmp_pw_percentile,cbt,data$deaths_10_19, dir="forw",cen=minT,
                   coef = readRDS(paste0(path_for_pooled, "pooled_age2_results_coef.rds")),
                   vcov = readRDS(paste0(path_for_pooled, "pooled_age2_results_vcov.rds")),
                   range=c(95,100))

edf_ci_age2 <- quantile(attrdl(data$tmp_pw_percentile,cbt,data$deaths_10_19, dir="forw",cen=minT,
                               coef = readRDS(paste0(path_for_pooled, "pooled_age2_results_coef.rds")),
                               vcov = readRDS(paste0(path_for_pooled, "pooled_age2_results_vcov.rds")),
                               range=c(95, 100), sim=T,nsim=1000), c(0.025,0.975))

edf_age3 <- attrdl(data$tmp_pw_percentile,cbt,data$deaths_20_34, dir="forw",cen=minT,
                   coef = readRDS(paste0(path_for_pooled, "pooled_age3_results_coef.rds")),
                   vcov = readRDS(paste0(path_for_pooled, "pooled_age3_results_vcov.rds")),
                   range=c(95,100))

edf_ci_age3 <- quantile(attrdl(data$tmp_pw_percentile,cbt,data$deaths_20_34, dir="forw",cen=minT,
                               coef = readRDS(paste0(path_for_pooled, "pooled_age3_results_coef.rds")),
                               vcov = readRDS(paste0(path_for_pooled, "pooled_age3_results_vcov.rds")),
                               range=c(95, 100), sim=T,nsim=1000), c(0.025,0.975))


edf_age4 <- attrdl(data$tmp_pw_percentile,cbt,data$deaths_35_64, dir="forw",cen=minT,
                   coef = readRDS(paste0(path_for_pooled, "pooled_age4_results_coef.rds")),
                   vcov = readRDS(paste0(path_for_pooled, "pooled_age4_results_vcov.rds")),
                   range=c(95,100))

edf_ci_age4 <- quantile(attrdl(data$tmp_pw_percentile,cbt,data$deaths_35_64, dir="forw",cen=minT,
                               coef = readRDS(paste0(path_for_pooled, "pooled_age4_results_coef.rds")),
                               vcov = readRDS(paste0(path_for_pooled, "pooled_age4_results_vcov.rds")),
                               range=c(95, 100), sim=T,nsim=1000), c(0.025,0.975))

edf_age5 <- attrdl(data$tmp_pw_percentile,cbt,data$deaths_65_plus, dir="forw",cen=minT,
                   coef = readRDS(paste0(path_for_pooled, "pooled_age5_results_coef.rds")),
                   vcov = readRDS(paste0(path_for_pooled, "pooled_age5_results_vcov.rds")),
                   range=c(95,100))

edf_ci_age5 <- quantile(attrdl(data$tmp_pw_percentile,cbt,data$deaths_65_plus, dir="forw",cen=minT,
                               coef = readRDS(paste0(path_for_pooled, "pooled_age5_results_coef.rds")),
                               vcov = readRDS(paste0(path_for_pooled, "pooled_age5_results_vcov.rds")),
                               range=c(95, 100), sim=T,nsim=1000), c(0.025,0.975))

edf_age5
edf_ci_age5
# mode
edf_motorcycle <- attrdl(data$tmp_pw_percentile,cbt,data$median_motorcycle_round, dir="forw",cen=minT,
                         coef = readRDS(paste0(path_for_pooled, "pooled_motorcycle_results_coef.rds")),
                         vcov = readRDS(paste0(path_for_pooled, "pooled_motorcycle_results_vcov.rds")),
                         range=c(95,100))

edf_ci_motorcycle <- quantile(attrdl(data$tmp_pw_percentile,cbt,data$median_motorcycle_round, dir="forw",cen=minT,
                                     coef = readRDS(paste0(path_for_pooled, "pooled_motorcycle_results_coef.rds")),
                                     vcov = readRDS(paste0(path_for_pooled, "pooled_motorcycle_results_vcov.rds")),
                                     range=c(95, 100), sim=T,nsim=1000), c(0.025,0.975))

edf_bicycle <- attrdl(data$tmp_pw_percentile,cbt,data$median_bicycle_round, dir="forw",cen=minT,
              coef = readRDS(paste0(path_for_pooled, "pooled_bicycle_results_coef.rds")),
              vcov = readRDS(paste0(path_for_pooled, "pooled_bicycle_results_vcov.rds")),
              range=c(95,100))

edf_ci_bicycle <- quantile(attrdl(data$tmp_pw_percentile,cbt,data$median_bicycle_round, dir="forw",cen=minT,
                          coef = readRDS(paste0(path_for_pooled, "pooled_bicycle_results_coef.rds")),
                          vcov = readRDS(paste0(path_for_pooled, "pooled_bicycle_results_vcov.rds")),
                          range=c(95, 100), sim=T,nsim=1000), c(0.025,0.975))

edf_ped <- attrdl(data$tmp_pw_percentile,cbt,data$median_ped_round, dir="forw",cen=minT,
                         coef = readRDS(paste0(path_for_pooled, "pooled_ped_results_coef.rds")),
                         vcov = readRDS(paste0(path_for_pooled, "pooled_ped_results_vcov.rds")),
                         range=c(95,100))

edf_ci_ped <- quantile(attrdl(data$tmp_pw_percentile,cbt,data$median_ped_round, dir="forw",cen=minT,
                                     coef = readRDS(paste0(path_for_pooled, "pooled_ped_results_coef.rds")),
                                     vcov = readRDS(paste0(path_for_pooled, "pooled_ped_results_vcov.rds")),
                                     range=c(95, 100), sim=T,nsim=1000), c(0.025,0.975))

edf_vehicle <- attrdl(data$tmp_pw_percentile,cbt,data$median_vehicle_round, dir="forw",cen=minT,
                  coef = readRDS(paste0(path_for_pooled, "pooled_vehicle_results_coef.rds")),
                  vcov = readRDS(paste0(path_for_pooled, "pooled_vehicle_results_vcov.rds")),
                  range=c(95,100))

edf_ci_vehicle <- quantile(attrdl(data$tmp_pw_percentile,cbt,data$median_vehicle_round, dir="forw",cen=minT,
                              coef = readRDS(paste0(path_for_pooled, "pooled_vehicle_results_coef.rds")),
                              vcov = readRDS(paste0(path_for_pooled, "pooled_vehicle_results_vcov.rds")),
                              range=c(95, 100), sim=T,nsim=1000), c(0.025,0.975))
