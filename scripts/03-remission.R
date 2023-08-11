####################### Run after 02-models ####################################
# Now calculate remission via recursion formula --------------------------------
r.mat   <- matrix(data = 0, nrow = nAges, ncol = nYears)                  
for (yIdx in 1:(nYears-1)) {
  for (aIdx in 1:(nAges-1)) {
    r.mat[aIdx,yIdx] <- ( po.mat[aIdx  ,yIdx  ] * (1 - m1.mat[aIdx  ,yIdx  ]) -
                          po.mat[aIdx+1,yIdx+1] * (1 - m0.mat[aIdx  ,yIdx  ]) +
                          inc.mat[aIdx ,yIdx  ] * (1 + po.mat[aIdx+1,yIdx+1])
                         )/(
                          po.mat[aIdx  ,yIdx  ] * (po.mat[aIdx+1,yIdx+1] + 1)
                          )
  }
}
summary(as.numeric(r.mat))

# Check values
hist(as.numeric(r.mat), breaks = 50)
validRvals <- 0 < as.numeric(r.mat) & as.numeric(r.mat) < 1
table(validRvals)

# Regression model for remission -----------------------------------------------
t.mat <- matrix(data = rep(years, each = nAges),                nrow = nAges, ncol = nYears)
a.mat <- matrix(data = rep(ages,  each = nYears), byrow = TRUE, nrow = nAges, ncol = nYears)

df.r <- data.frame(rem = as.numeric(r.mat[validRvals]),
                     t = as.numeric(t.mat[validRvals]),
                     a = as.numeric(a.mat[validRvals])
                   )

r.mod <- lm(logit(rem) ~ ns(a, df = 6)*ns(t, df = 4), data = df.r) #            
summary(r.mod)

fct_rem <- function(t, a){
  expit(predict(r.mod, newdata = data.frame(t = rep(t, length(a)), a = a)))
}

# Visual check
ages.e <- 2.5 + 5*0:17
this.year <- 1995
matplot(ages.e, fct_rem(this.year, ages.e), type = "l", las = 1)
points(ages, r.mat[,6], pch = 20, col = "blue") # select 1995 in r.mat
