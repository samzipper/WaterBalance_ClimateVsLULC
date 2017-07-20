## PenmanMonteith_Daily.R
# Goal: calculate Penman-Monteith reference ET (ETo) as described in FAO Paper #56 at an hourly timestep.
# This is based on the DailyPET.m script by Eric Booth, adapted for R and vectorized.
#
# Reference: Allen, R., Pereira, L.S., Raes, D. and Smith, M., 1998. Crop 
# evapotranspiration: Guidelines for computing crop requirements, FAO, 
# Rome, Italy. http://www.fao.org/docrep/X0490E/x0490e00.htm

PenmanMonteith <- function(metData){
  # metData should be a vector of: 
  #    DOY = day, 
  #    srad = [W m-2], 
  #    TMAX [C], 
  #    TMIN [C], 
  #    es_mean [kPa],
  #    ea_mean [kPa],
  #    prec [mm], 
  #    wind [m/s],
  #    LAT [degN],
  #    ELEV [m]
  # field names are intended to line up with U of Idaho Metdata daily gridded data
  
  require(lubridate)
  
  # initialize some matrices
  ETo <- data.frame(DOY=metData$DOY, ETo=NaN)
  
  # process your input data as necessary
  TMEAN <- (metData$TMAX+metData$TMIN)/2
  srad <- metData$srad*0.0864  # convert from [W m-2] to [MJ m-2 d-1]
  
  # convert LAT to radians
  LAT.rad <- (pi/180)*metData$LAT        # convert latitude to radians
  
  # calculate pressure based on elevation
  P <- 101.3*((293-0.0065*metData$ELEV)/293)^5.26    # Eq. 7
  
  # Determine psychrometric constant
  psy <- 0.665e-3*P
  
  ## Determine air humidity values
  # actual vapor pressure      # eqn 17
  # slope of saturation vapor pressure curve at air temperature
  delta <- 4098*0.6108*exp(17.27*TMEAN/(TMEAN+237.3))/(TMEAN+237.3)^2; # eqn 13
  
  ## Determination of extraterrestrial radiation, R_a
  # inverse relative distance Earth-Sun
  d_r <- 1 + 0.033*cos(2*pi*metData$DOY/365); # eqn 23
  # solar decimation
  del <- 0.409*sin(2*pi*metData$DOY/365-1.39); # eqn 24
  # sunset hour angle
  w_s <- acos(-tan(LAT.rad)*tan(del)); # eqn 25
  # extraterrestrial radiation
  R_a <- (24*60/pi)*0.0820*d_r*(w_s*sin(LAT.rad)*sin(del) + 
                                  cos(LAT.rad)*cos(del)*sin(w_s)); # eqn 21
  
  ## Determination of clear-sky solar radiation, R_so
  R_so <- (0.75 + 2e-5*metData$ELEV)*R_a; # eqn 37
  # make sure that srad is less than or equal to R_so
  srad[srad > R_so] <- R_so[srad > R_so]

  
  ## Determination of net incoming solar radiation, R_ns
  albedo <- 0.23; # value for hypothetical grass reference crop
  R_ns <- (1 - albedo)*srad; # eqn 38
  
  ## Determination of net outgoing longwave radiation, R_nl
  s <- 4.903e-9; # Stefan-Boltzmann constant
  TmaxK4 <- (metData$TMAX + 273.16)^4;
  TminK4 <- (metData$TMIN + 273.16)^4;
  R_nl <- s*(TmaxK4+TminK4)/2*(0.34-0.14*metData$ea_mean^0.5)*(1.35*srad/R_so-0.35); # eqn 39
  
  ## Determination of net radiation, R_n
  R_n <- R_ns - R_nl; # eqn 40
  
  ## Calculate grass reference potential evapotranspiration, ETo
  ETo$ETo <- (0.408*delta*R_n+psy*(900/(TMEAN+273))*metData$wind*(metData$es_mean-metData$ea_mean))/
    (delta + psy*(1 + 0.34*metData$wind)) # [mm/day] eqn 6
  
  ## if negative ETo values calculated, set equal to 0
  # this can happen at high latitudes in winter when incoming solar
  # radiation is very low, leading to a negative net radiation term.
  ETo$ETo[ETo$ETo<0] <- 0
  
  return(ETo)
}

