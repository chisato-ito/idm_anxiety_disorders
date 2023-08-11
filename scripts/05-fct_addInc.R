####################### Run after 04-covid_2020-2022 ###########################
# Peaks of COVID-19 waves (in weeks relative to 1st week in 2019):--------------
# Vary w.in, w.out, delay, lambda, incidence spike h0

fct_addInc <- function(t, a, w.in, w.out, delay, lambda, h0){
 
  # description of COVID-waves (measured in weeks relative to 1st week in 2019)
  waveMidTime  <- 2019 + (52 + peaks.position)/52  # where are midpoints of waves located
  waveDuration <-              fwhm/52  # how long do waves last (FWHM)
  # end description of waves
  
  t.pts <- c()
  vals  <- c()
  nWaves      <- length(waveMidTime)
  
  for(wNr in 1:nWaves){
    thisP   <- waveMidTime[wNr]
    thisDur <- waveDuration[wNr]
    thisDel <- delay * thisDur
    washIn  <- w.in  * thisDur
    washOut <- w.out * thisDur
    mhIncH  <- h0 * exp(-lambda*(wNr-1)) # incidence increase with damping effect

    t.pts   <- c(t.pts, thisP + c(-0.5*thisDur, -0.5*thisDur + washIn, 0.5*thisDur + thisDel, 0.5*thisDur + thisDel + washOut))
    vals    <- c( vals,         c(     0      ,        mhIncH        ,          mhIncH      ,                 0              ))
  }
  fac_      <- 1 + approx(t.pts, vals, xout = t, method = "linear", rule =2, ties = "ordered")$y
  return(fac_)
}
