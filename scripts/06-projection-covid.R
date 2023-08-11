####################### Run after 05-fct_addInc ###############################
# Specify projection time frame-------------------------------------------------
# Model 2019 to 2030
yearsToModel  <- 2019:2030
nProjYears    <- length(yearsToModel)

# Select relevant years from Destatis population projection
###
yearsPast <- 2019:2021 # In the 15th projection 2019-2021 are historical data (Variant0)
yearsFuture <- 2022:2030
N0             <- t(as.matrix(pop[pop$Variante == 0 & # Variant0:
                                   pop$mw == sexString2 &
                                   pop$Simulationsjahr %in% yearsPast, 5+ages]))
N1             <- t(as.matrix(pop[pop$Variante == 2 & # Variant2:
                                    pop$mw == sexString2 &
                                    pop$Simulationsjahr %in% yearsFuture, 5+ages]))
N <- cbind(N0,N1)
###
dt            <- 1/52 #temporal resolution in units per year: weekly
timeSteps     <- seq(from = yearsToModel[1], to = yearsToModel[nProjYears], by = dt)
nTimeSteps    <- length(timeSteps)
ageSteps      <- seq(from = ages[1], to = ages[nAges], by = dt)
nAgeSteps     <- length(ageSteps)

S.hd <- matrix(data = 0, nrow = nAgeSteps, ncol = nTimeSteps) # susceptible
C.hd <- matrix(data = 0, nrow = nAgeSteps, ncol = nTimeSteps) # cases

# Upsample population with an implicit assumption that births are uniformly distributed
# over year
# People entering population
N0.hd <- dt*approx(yearsToModel, N[1,], xout = timeSteps)$y 
# Age-distribution at the start year
N.hd <- dt*approx(ages, N[,1], xout = ageSteps)$y

#-------------------------------------------------------------------------------

# Initial point in time
p_        <- fct_p(timeSteps[1], ageSteps)
S.hd[,1]  <- N.hd*(1 - p_)
C.hd[,1]  <- N.hd*p_

#-------------------------------------------------------------------------------
# Without COVID-19
#-------------------------------------------------------------------------------

for (tIdx in 2:nTimeSteps) {
  ii_ <- fct_inc(timeSteps[tIdx], ageSteps[-1])
  rr_ <- fct_rem(timeSteps[tIdx], ageSteps[-1])
  m0  <- fct_m0 (timeSteps[tIdx], ageSteps[-1])
  m1  <- fct_m1 (timeSteps[tIdx], ageSteps[-1])
  
  # S, C from previous step
  Sp_ <- S.hd[-nAgeSteps, tIdx-1]
  Cp_ <- C.hd[-nAgeSteps, tIdx-1]
  
  # S, C in actual step
  S_  <- Sp_ * (1 - (ii_ + m0)*dt) + Cp_ * rr_*dt
  C_  <- Cp_ * (1 - (rr_ + m1)*dt) + Sp_ * ii_*dt
  
  S.hd[,tIdx] <- c(N0.hd[tIdx], S_)
  C.hd[,tIdx] <- c(          0, C_)
}

# case number
plot(timeSteps, colSums(C.hd), pch = '.', ylim = c(0, 5000))
# prevalence
plot(timeSteps, colSums(C.hd)/(colSums(C.hd)+colSums(S.hd)), pch = '.', ylim = c(0,0.2))

# Total cases & susceptible for each time step
df.C <- data.frame("year" = timeSteps, 
                   "susceptible" = colSums(S.hd),
                   "case" = colSums(C.hd), 
                   "prev" = colSums(C.hd) / (colSums(C.hd) + colSums(S.hd)))
assign(paste0("df_",sexString2), df.C)
# Export to CSV
write_csv(df.C, file = paste0("out_data/","df_",sexString2,"_cov0.csv"))

#-------------------------------------------------------------------------------
# With COVID-19
#-------------------------------------------------------------------------------

# parameters
w.in   <- c(0.1, 0.2 , 0.3 ) # w.in*100  = % of wave duration
w.out  <- c(0.1, 0.2 , 0.3 ) # w.out*100 = % of wave duration 
delay  <- c(0  , 0.25, 0.5 ) # delay*100 = % of wave duration
lambda <- c(0.1, 0.2 , 0.3 )
h0     <- c(1.0, 5.0 , 10.0) # h0*100   = % increase in incidence               

var <- data.frame(w.in, w.out, delay, lambda, h0)
param <- expand.grid(var) %>% filter(w.in == w.out)

# create empty lists
Smatlist <- list()
Cmatlist <- list()

for (i in 1:nrow(param)) {
  for (tIdx in 2:nTimeSteps) {
    ii_ <- fct_inc(timeSteps[tIdx], ageSteps[-1]) * 
           fct_addInc(timeSteps[tIdx], ageSteps[-1], param[i,1], param[i,2], param[i,3], param[i,4], param[i,5])
    rr_ <- fct_rem(timeSteps[tIdx], ageSteps[-1])
    m0_ <- fct_m0 (timeSteps[tIdx], ageSteps[-1])
    m1_ <- fct_m1 (timeSteps[tIdx], ageSteps[-1])
    
    # S, C from previous step
    Sp_ <- S.hd[-nAgeSteps, tIdx-1]
    Cp_ <- C.hd[-nAgeSteps, tIdx-1]
    
    # S, C in actual step
    S_  <- Sp_ * (1 - (ii_ + m0_)*dt) + Cp_ * rr_*dt
    C_  <- Cp_ * (1 - (rr_ + m1_)*dt) + Sp_ * ii_*dt
    
    S.hd[,tIdx] <- c(N0.hd[tIdx], S_)
    C.hd[,tIdx] <- c(          0, C_)

  }
  Smatlist[[i]] <- assign(paste0("S.hd_out_", i), S.hd)
  Cmatlist[[i]] <- assign(paste0("C.hd_out_", i), C.hd)
}

# data frames for all scenarios
#for (i in 1:27) {
for (i in 1:81) {
  df <- data.frame("year" = timeSteps,
             "susceptible" = colSums(Smatlist[[i]]),
             "case" = colSums(Cmatlist[[i]]),
             "prev" = colSums(Cmatlist[[i]])/(colSums(Smatlist[[1]])+colSums(Cmatlist[[i]])))
  assign(paste0("df_",sexString2,"_cov",i), df)
  write_csv(df, file = paste0("out_data/","df_",sexString2,"_cov",i,".csv")) 
}

#-------------------------------------------------------------------------------

