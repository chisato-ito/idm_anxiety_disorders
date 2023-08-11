####################### Run after 01-prep_data #################################
# Incidence --------------------------------------------------------------------
# Incidence rate data frame by sex
idat <- gbd[gbd$sex_name == sexString &
              gbd$metric_name == "Rate" &
              gbd$measure_name == "Incidence",
            c("year", "age", "val", "upper", "lower")]
summary(idat$val)

# Fit model                                                                               
i.mod <- lm(log(1e-5*val) ~ ns(age, df = 10)*ns(year, df = 2), data = idat)
summary(i.mod)

fct_inc <- function(t, a){
  exp(predict(i.mod, newdata = data.frame(year = rep(t, length(a)), age = a)))
}

# Check fit
ages.e <- 2.5 + 5*0:17
this.year <- 1995
matplot(ages.e, fct_inc(this.year, ages.e), type = "l", pch = 16, las = 1, ylim = c(0, .1))
with(idat, matplot(age[year == this.year], 1e-5*val[year == this.year], type = "l", col = "blue", add = TRUE))
for (ii in 1:length(unique(idat$age))) {
  aa <- (idat$age[idat$year == this.year]) [ii]
  up <- (idat$upper[idat$year == this.year]) [ii]
  lw <- (idat$lower[idat$year == this.year]) [ii]
  lines(rep(aa,2), 1e-5*c(lw,up), col = "blue")
}

# Prevalence -------------------------------------------------------------------
pdat <- gbd[gbd$sex_name == sexString &
              gbd$metric_name == "Percent" &
              gbd$measure_name == "Prevalence",
            c("year", "age", "val", "upper", "lower")]
summary(pdat$val)

# Fit model                                                                     
p.mod <- lm(logit(val) ~ ns(age, df = 10)*ns(year, df = 2), data = pdat, subset = val > 0)
summary(p.mod)  

fct_p <- function(t, a){
  expit(predict(p.mod, newdata = data.frame(year = rep(t, length(a)), age = a)))
}

# Check fit
ages.e <- 2.5 + 5*0:17
this.year <- 1995
matplot(ages.e, fct_p(this.year, ages.e), type = "l", pch = 16, las = 1, ylim = c(0, .2))
with(pdat, matplot(age[year == this.year], val[year == this.year], type = "l", col = "blue", add = TRUE))
for (ii in 1:length(unique(pdat$age))) {
  aa <- (pdat$age[pdat$year == this.year]) [ii]
  up <- (pdat$upper[pdat$year == this.year]) [ii]
  lw <- (pdat$lower[pdat$year == this.year]) [ii]
  lines(rep(aa,2), c(lw, up), col = "blue")
}

# Mortality --------------------------------------------------------------------
# We assume a mortality rate ratio of 1.4 based on Walker et al. 2015
# Mortality among the non-diseased, with a log-term for annual decrease         
fct_m0 <- function(t, a){
  exp(-10.7 + 0.1*a + (t-1990)*log(1 - 0.005))
}
# Mortality among the diseased
fct_m1 <- function(t, a){
  exp(-10.4 + 0.095*a + (t-1990)*log(1 - 0.007))
}

# Check
ages.e <- 2.5 + 5*0:17
this.year <- 1995
matplot(ages.e, fct_m0(this.year, ages.e), type = "l", log = "y", las = 1, ylim = c(1e-5, 0.15))
matplot(ages.e, fct_m1(this.year, ages.e), type = "l", add = TRUE, col = "red")

# Additional time-scale: period ------------------------------------------------
years   <- 1990:2019
ages    <- 1:95
nYears  <- length(years)
nAges   <- length(ages)

p.mat   <- matrix(data = 0, nrow = nAges, ncol = nYears)
po.mat  <- matrix(data = 0, nrow = nAges, ncol = nYears)
for (yIdx in 1:nYears) {
  pp_ <- fct_p(years[yIdx], ages)
  p.mat [,yIdx] <- pp_
  po.mat[,yIdx] <- pp_/(1-pp_)
}

inc.mat <- matrix(data = 0, nrow = nAges, ncol = nYears)
for (yIdx in 1:nYears) {
  inc.mat[,yIdx] <- fct_inc(years[yIdx], ages)
}

m0.mat  <- matrix(data = 0, nrow = nAges, ncol = nYears)
for (yIdx in 1:nYears) {
  m0.mat[,yIdx] <- fct_m0(years[yIdx], ages)
}

m1.mat  <- matrix(data = 0, nrow = nAges, ncol = nYears)
for (yIdx in 1:nYears) {
  m1.mat[,yIdx] <- fct_m1(years[yIdx], ages)
}
