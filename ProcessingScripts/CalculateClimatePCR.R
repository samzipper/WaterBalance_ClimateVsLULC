CalculateClimatePCR <- function(df, mo, flux.name, var.options, 
                                yr.baseline.start, yr.baseline.end, n.val.yr,
                                p.thres=0.10, cum.var=0.80, min.var=0.01, n.perm=250, seed=1,
                                neg.allowed=F,  write.vars.keep=T, write.PC.keep=T, write.perm=T,
                                PC.plot=F, plot.dir=NULL){
  #' This script is intended to use principal components regression (PCR) to 
  #' build statistical relationships between meteorological variables and a
  #' hydrological flux of interest.
  #' 
  #' required inputs:
  #'  -df                = data frame with input data (described below)
  #'  -mo                = month you want to build PCR for [numeric]
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
  #'  -write.vars.keep   = write out vars retained? [logical]
  #'  -write.PC.keep     = write out retained principal components? [logical]
  #'  -write.perm        = write out permutation results? [logical]
  #'  -PC.plot           = save plots showing principal component loading? [logical]
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
  #'    -PCR   = predicted flux of interest [numeric]
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
      }
    }
  }
  
  # report error is <= 1 var retained
  if (length(vars.keep)<=1){
    stop(paste0(flux.name, " mo ", mo, ": <= 1 var retained (", vars.keep, ")"))
  }
  
  # write out variables to retain
  if (write.vars.keep){
    write.csv(data.frame(var=vars.keep), paste0("GHCN_MonthlyRegressions_vars.keep_", sprintf("%02d", mo), ".csv"),
              row.names=F, quote=F)
  }
  
  # trim df.mo to vars only
  df.mo.vars <- df.mo[, vars.keep]
  
  ## 5. Calculate principal components
  # calculate PC rotations based on 1986-2013 data only
  fmla.vars.keep <- as.formula(paste0("~ ", paste0(vars.keep, collapse="+")))
  PCA.fit <- prcomp(fmla.vars.keep, data=df.mo.vars[df.mo$year <= yr.baseline.end, ], na.action=na.omit)
  
  # apply PCA.fit to entire period
  df.mo.PCs <- as.data.frame(predict(PCA.fit, newdata=df.mo.vars))
  
  # build a model
  df.mo.PCs$flux <- df.mo$flux     # add in your desired predictor variable
  df.mo.PCs$year <- df.mo$year
  
  ## select number of principal components to use for final relationship
  # first, select those that explain a cumulative 95% of the total variance
  PC.95 <- min(which(summary(PCA.fit)$importance["Cumulative Proportion", ] >= cum.var))  # find which PC explains >95% of cumulative proportion
  PC.keep <- colnames(summary(PCA.fit)$importance)[1:PC.95]
  
  # now, look through remaining PCs and see if any are significantly correlated with output variable, and retain those as well
  # only select from those that explain >0.5% of total variance observed in predictors
  if (sum(summary(PCA.fit)$importance["Proportion of Variance", ]<min.var) > 0){
    PC.max <- min(which(summary(PCA.fit)$importance["Proportion of Variance", ] < min.var))  # max allowed PC
  } else {
    # if all PCs are above min.var, add 1 to PC.max
    PC.max <- length(summary(PCA.fit)$importance["Proportion of Variance", ])+1
  }
  
  if (PC.max > PC.95){
    PC.options <-paste0("PC", (PC.95+1):(PC.max-1))  # number of PCs will correspond to length of vars.keep
    for (PC in PC.options){
      fmla.PC <- paste0("flux ~ ", PC)
      fit.PC <- lm(fmla.PC, data=df.mo.PCs)
      if (lmp(fit.PC) <= p.thres){
        PC.keep <- c(PC.keep, PC)
      }
    }
  }
  
  if (PC.plot){
    # extract components: first 10 only
    PCA.comps <- PCA.fit$rotation[,PC.keep]
    
    # prepare data to plot
    PCA.comps.plot <- as.data.frame(PCA.comps)
    PCA.comps.plot$varName <- rownames(PCA.comps.plot)
    PCA.comps.plot.long <- gather(PCA.comps.plot, "PC", "loading", starts_with("PC"))
    
    # plot singular values
    p.magnitudes <-
      dotplot_identity(frame = data.frame(pc=1:length(PCA.fit$sdev), 
                                          magnitude=PCA.fit$sdev), 
                       xvar="pc",yvar="magnitude") +
      geom_hline(yintercept=0, color="gray65") +
      scale_x_continuous(name="PC", breaks=seq(1,length(vars.keep), 1)) +
      ylab("Magnitude") +
      ggtitle("Magnitudes of singular values") +
      theme_bw() +
      theme(panel.grid=element_blank())
    
    # plot variable loadings
    p.loadings <-
      dotplot_identity(PCA.comps.plot.long, "varName", "loading") + 
      facet_wrap(~PC,nrow=1) + coord_flip() + 
      geom_hline(yintercept=0, color="gray65") +
      xlab("Variable") +
      ylab("Loading") +
      ggtitle("x scaled variable loadings") +
      theme_bw() +
      theme(panel.grid=element_blank())
    
    # save plots
    ggsave(paste0(plot.dir, "GHCN_MonthlyRegressions_PCRloadings_", sprintf("%02d", mo), ".png"),
           arrangeGrob(p.magnitudes, p.loadings, ncol=1, heights=c(0.5,1)), width=12, height=12, units="in")
  }
  
  # write out PCs to retain
  if (write.PC.keep){
    write.csv(data.frame(PC=PC.keep), paste0("GHCN_MonthlyRegressions_PC.keep_", sprintf("%02d", mo), ".csv"),
              row.names=F, quote=F)
  }
  
  # build formula based on retained PCs
  ncomp <- length(PC.keep)  # number of principal components you want to use
  fmla <- paste("flux ~", paste(PC.keep, collapse='+'))  # build your formula
  
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
  for (PC in PC.keep){
    df.perm[,dim(df.perm)[2]+1] <- NaN
    colnames(df.perm)[dim(df.perm)[2]] <- PC
  }
  
  # store indices of validation cells
  perm.val <- matrix(ncol=n.val.yr, nrow=n.perm)
  
  for (p in 1:n.perm){
    # randomly assign to calibration, validation, or prediction
    df.mo.PCs$group <- "prediction"  # set all years to prediction category
    df.mo.PCs$group[df.mo.PCs$year <= yr.baseline.end] <- "cal"  
    df.mo.PCs$group[sample(seq(1,sum(df.mo.PCs$group=="cal")), n.val.yr)] <- "val"   # random sample of years for validation
    
    # save indices of validation
    perm.val[p, ] <- which(df.mo.PCs$group=="val")
    
    # fit model
    fit.var <- lm(fmla, data=subset(df.mo.PCs, group=="cal"))
    
    # estimate output using model
    df.mo.PCs$estimate <- predict(fit.var, newdata=df.mo.PCs)
    
    # don't allow negative estimates
    df.mo.PCs$estimate[df.mo.PCs$estimate<0] <- 0
    
    # get rid of missing values
    df.mo.PCs <- subset(df.mo.PCs, is.finite(estimate) & is.finite(flux))
    
    #fit metrics
    df.perm$proj.RMSE.val[p] <- RMSE(df.mo.PCs$estimate[df.mo.PCs$group=="val"], df.mo.PCs$flux[df.mo.PCs$group=="val"])
    df.perm$proj.NRMSE.val[p] <- NRMSE(df.mo.PCs$estimate[df.mo.PCs$group=="val"], df.mo.PCs$flux[df.mo.PCs$group=="val"])
    df.perm$proj.R2.val[p] <- R2(df.mo.PCs$estimate[df.mo.PCs$group=="val"], df.mo.PCs$flux[df.mo.PCs$group=="val"])
    df.perm$proj.RMSE.cal[p] <- RMSE(df.mo.PCs$estimate[df.mo.PCs$group=="cal"], df.mo.PCs$flux[df.mo.PCs$group=="cal"])
    df.perm$proj.NRMSE.cal[p] <- NRMSE(df.mo.PCs$estimate[df.mo.PCs$group=="cal"], df.mo.PCs$flux[df.mo.PCs$group=="cal"])
    df.perm$proj.R2.cal[p] <- R2(df.mo.PCs$estimate[df.mo.PCs$group=="cal"], df.mo.PCs$flux[df.mo.PCs$group=="cal"])
    
    # collect all output data
    df.out <- data.frame(year=df.mo.PCs$year,
                         month=mo,
                         perm=p,
                         flux = df.mo.PCs$flux,
                         PCR = df.mo.PCs$estimate,
                         group = df.mo.PCs$group)
    if (exists("df.out.all")){
      df.out.all <- rbind(df.out.all, df.out)
    } else {
      df.out.all <- df.out
    }
    
    # fill in coefficients
    n.coef <- length(PC.keep)+1
    df.perm[p, (dim(df.perm)[2]-n.coef+1):dim(df.perm)[2]] <- as.numeric(coef(fit.var))
  }
  
  if (write.perm){
    # save df.perm data frame
    write.csv(df.perm, paste0("GHCN_MonthlyRegressions_permutations_", sprintf("%02d", mo), ".csv"),
              row.names=F, quote=F)
    
    # save df.perm data frame
    write.csv(perm.val, paste0("GHCN_MonthlyRegressions_perm.val_", sprintf("%02d", mo), ".csv"),
              row.names=F, quote=F)
  }
  
  return(df.out.all)
}