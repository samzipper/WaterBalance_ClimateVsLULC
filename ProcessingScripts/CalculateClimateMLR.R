CalculateClimateMLR <- function(df, mo, flux.name, var.options, 
                                yr.baseline.start, yr.baseline.end, n.val.yr,
                                p.thres=0.10, n.perm=250, seed=1,
                                neg.allowed=F,  n.vars.keep=1, write.vars.keep=T, write.perm=T,
                                plot.dir=NULL){
  #' This script is intended to use multiple linear regression (MLR) to 
  #' build statistical relationships between meteorological variables and a
  #' hydrological flux of interest.
  #' 
  #' required inputs:
  #'  -df                = data frame with input data (described below)
  #'  -mo                = month you want to build MLR for [numeric]
  #'  -flux.name         = name of column containing flux of interest [character]
  #'  -var.options       = name of candidate predictor variables [character vector]
  #'  
  #'  -yr.baseline.start = start year of baseline period (inclusive) [numeric]
  #'  -yr.baseline.end   = end year of baseline period (inclusive) [numeric]
  #'  -n.val.yr          = number of years to use for model validation (rest of baseline period will be calibration) [numeric]
  #'  
  #'  -p.thres           = significance threshold for pruning [numeric]
  #'  -cum.var           = cumulative variance explained used for selecting retained PCs [numeric, 0-1]
  #'  -min.var           = minimum variance explained used for selecting retained PCs [numeric, 0-1]
  #'  -n.perm            = number of permutations to iterate for uncertainty calculations [numeric]
  #'  -seed              = seed value for random number generator [numeric]
  #'  
  #'  -neg.allowed       = are negative predictions allowed? [logical] - should be F for streamflow and ET, T for drainage
  #'  -n.vars.keep       = number of predictors to retain [numeric]
  #'  -write.vars.keep   = write out vars retained? [logical]
  #'  -write.perm        = write out permutation results? [logical]
  #'  -plot.dir          = path to save plots
  #' 
  #' df is a data frame containing the following columns:
  #'  -year
  #'  -month
  #'  -[flux.name] = hydrological flux of interest
  #'  -[var.options] = one column for each candidate predictor variable
  #'  All other columns will be ignored.
  #'  
  #' output returned is:
  #'  -df.out.all = data frame with the following columns
  #'    -year  = year [numeric]
  #'    -month = month [numeric]
  #'    -perm  = permutation number [numeric]
  #'    -flux  = flux of interest [numeric]
  #'    -MLR   = predicted flux of interest [numeric]
  #'    -group = group (cal, val, or prediction) [character]
  
  ## define other useful functions
  scale.baseline <- function(x,year,yr.start,yr.end){
    # this function is intended to take in a vector x of data and associated 
    # year timestamps, and and scale x based on the mean and standard deviation only 
    # within the baseline period, defined by yr.start and yr.end (inclusive)
    
    # find measurements within baseline period
    i.baseline <- (year >= yr.start & year <= yr.end)
    
    # find mean and standard deviation of baseline period
    mean.baseline <- mean(x[i.baseline], na.rm=T)
    sd.baseline <- sd(x[i.baseline], na.rm=T)
    
    # scale whole record based on baseline period
    x.scale <- (x-mean.baseline)/sd.baseline
    
    return(x.scale)
  }
  lmp <- function (modelobject) {
    # calculate the p-value for a linear regression
    # i got this off the internet but can't remember where... (StackOverflow I think)
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }
  
  # rename flux of interest
  colnames(df)[colnames(df)==flux.name] <- "flux"
  
  # get data frame for this month
  df.mo <- subset(df, month==mo & is.finite(flux))
  
  # significance pruning- go through and select variables with significant linear relationships in historial period
  vars.keep <- character(0)
  R2.keep <- numeric(0)
  for (v in var.options){
    # scale to mean and sd
    col.v <- which(colnames(df.mo)==v)
    df.mo[,col.v] <- scale.baseline(df.mo[,col.v], df.mo$year, yr.baseline.start, yr.baseline.end)
    if (length(unique(df.mo[df.mo$year<=yr.baseline.end,v]))>1){
      # only check this potential variable if it has variability
      fmla.v <- paste0("flux ~ ", v)
      fit.v <- lm(fmla.v, data=subset(df.mo, year<=yr.baseline.end))
      if (lmp(fit.v) <= p.thres){
        vars.keep <- c(vars.keep, v)
        R2.keep <- c(R2.keep, summary(fit.v)$r.squared)
      }
    }
  }
  
  # pick only the best vars
  vars.keep <- vars.keep[order(-R2.keep)][1:n.vars.keep]
  
  # report error is <= 1 var retained
  if (length(vars.keep)<=1){
    stop(paste0(flux.name, " mo ", mo, ": <= 1 var retained (", vars.keep, ")"))
  }
  
  # write out variables to retain
  if (write.vars.keep){
    write.csv(data.frame(var=vars.keep), paste0("GHCN_MonthlyRegressions_MLR_vars.keep_", sprintf("%02d", mo), ".csv"),
              row.names=F, quote=F)
  }
  
  # trim df.mo to vars only
  df.mo.vars <- df.mo[, c("year", "flux", vars.keep)]
  
  ## 5. Build MLR
  fmla.vars.keep <- as.formula(paste0("flux~ ", paste0(vars.keep, collapse="+")))
  
  ## permutation-based approach to uncertainty analysis
  set.seed(seed)

  # empty data frame to hold output
  df.perm <- data.frame(perm = seq(1,n.perm),
                        proj.RMSE.cal = NaN,
                        proj.RMSE.val = NaN,
                        proj.NRMSE.cal = NaN,
                        proj.NRMSE.val = NaN,
                        proj.R2.cal = NaN,
                        proj.R2.val = NaN,
                        fit.intercept = NaN)  # intercept
  # make a column to store coefficient for each PC
  for (var in vars.keep){
    df.perm[,dim(df.perm)[2]+1] <- NaN
    colnames(df.perm)[dim(df.perm)[2]] <- var
  }
  
  # store indices of validation cells
  perm.val <- matrix(ncol=n.val.yr, nrow=n.perm)
  
  for (p in 1:n.perm){
    # randomly assign to calibration, validation, or prediction
    df.mo.vars$group <- "prediction"  # set all years to prediction category
    df.mo.vars$group[df.mo.vars$year <= yr.baseline.end] <- "cal"  
    df.mo.vars$group[sample(seq(1,sum(df.mo.vars$group=="cal")), n.val.yr)] <- "val"   # random sample of years for validation
    
    # save indices of validation
    perm.val[p, ] <- which(df.mo.vars$group=="val")
    
    # fit model
    fit.var <- lm(fmla.vars.keep, data=subset(df.mo.vars, group=="cal"))
    
    # estimate output using model
    df.mo.vars$estimate <- predict(fit.var, newdata=df.mo.vars)
    
    if (!neg.allowed){
      # don't allow negative estimates
      df.mo.vars$estimate[df.mo.vars$estimate<0] <- 0
      
    }
    
    # get rid of missing values
    df.mo.vars <- subset(df.mo.vars, is.finite(estimate) & is.finite(flux))
    
    #fit metrics
    df.perm$proj.RMSE.val[p] <- RMSE(df.mo.vars$estimate[df.mo.vars$group=="val"], df.mo.vars$flux[df.mo.vars$group=="val"])
    df.perm$proj.NRMSE.val[p] <- NRMSE(df.mo.vars$estimate[df.mo.vars$group=="val"], df.mo.vars$flux[df.mo.vars$group=="val"])
    df.perm$proj.R2.val[p] <- R2(df.mo.vars$estimate[df.mo.vars$group=="val"], df.mo.vars$flux[df.mo.vars$group=="val"])
    df.perm$proj.RMSE.cal[p] <- RMSE(df.mo.vars$estimate[df.mo.vars$group=="cal"], df.mo.vars$flux[df.mo.vars$group=="cal"])
    df.perm$proj.NRMSE.cal[p] <- NRMSE(df.mo.vars$estimate[df.mo.vars$group=="cal"], df.mo.vars$flux[df.mo.vars$group=="cal"])
    df.perm$proj.R2.cal[p] <- R2(df.mo.vars$estimate[df.mo.vars$group=="cal"], df.mo.vars$flux[df.mo.vars$group=="cal"])
    
    # collect all output data
    df.out <- data.frame(year=df.mo.vars$year,
                         month=mo,
                         perm=p,
                         flux = df.mo.vars$flux,
                         MLR = df.mo.vars$estimate,
                         group = df.mo.vars$group)
    if (exists("df.out.all")){
      df.out.all <- rbind(df.out.all, df.out)
    } else {
      df.out.all <- df.out
    }
    
    # fill in coefficients
    n.coef <- length(vars.keep)+1
    df.perm[p, (dim(df.perm)[2]-n.coef+1):dim(df.perm)[2]] <- as.numeric(coef(fit.var))
  }
  
  if (write.perm){
    # save df.perm data frame
    write.csv(df.perm, paste0("GHCN_MonthlyRegressions_MLR_permutations_", sprintf("%02d", mo), ".csv"),
              row.names=F, quote=F)
    
    # save df.perm data frame
    write.csv(perm.val, paste0("GHCN_MonthlyRegressions_MLR_perm.val_", sprintf("%02d", mo), ".csv"),
              row.names=F, quote=F)
  }
  
  return(df.out.all)
}