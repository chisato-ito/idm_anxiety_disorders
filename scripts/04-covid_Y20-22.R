####################### Run after 03-remission #################################
# COVID-19 data: Germany 2020-2022 (7 waves)
#library(ggpmisc)

# Read in SurvStat@RKI COVID data
cov <- read.csv(here("raw_data", "survstat_wk_20230127", "Data1.csv"), 
                 col.names = c("wk", "2020", "2021","2022","2023"),            
                 na.strings=c(""), skip = 1)
cov$X2022 <- as.numeric(gsub(",","",cov$X2022)) # remove "," for 1000s
cov <- cov[1:4] # Select until the end of 2022
str(cov)  # 53 weeks in 2020, 52 weeks in 2021 & 2022

cov <- cov %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%  # replace NA with 0
  pivot_longer(cols = 2:4, names_to = "year", names_prefix = "Y", values_to = "val") %>% 
  filter(!((wk == 53 & year == "X2021") | (wk == 53 & year == "X2022"))) %>%    
  # Remove 53rd week in 2021 & 2022 (there is no wk 53)
  mutate(week = ifelse(year == "X2020", wk, ifelse(year == "X2021", wk+53, wk+105))) %>% 
  arrange(week)
str(cov)

# Visualize local maxima with ggpmisc package
p1 <- ggplot(cov, aes(week, val)) +
  geom_line() +
  stat_peaks(col="red",
             ignore_threshold = .01,      # .01 to identify the first wave
             geom = "text_s",
             position = position_nudge_to(y = 300),
             arrow = arrow(length = grid::unit(1.5, "mm")),
             hjust = 0.2,
             x.label.fmt = "%.0f")
p1

# Get local maxima to give starting values for modeling
(wk.maxima <- cov$week[ggpmisc:::find_peaks(cov$val, ignore_threshold = .01)]) 
(wk.peaks <- wk.maxima[c(1,3,6,8,10,11,12)]) # visually select peaks - 7 waves in total
(val.peaks <- cov$val[wk.peaks])

# Fit a model - multiple Gaussians ---------------------------------------------
wk <- cov$week
val <- cov$val
df <- data.frame(wk, val)

cov.mod <- nls(val ~ (a1*exp(-((wk-b1)^2)/(2*c1^2))) + # wave 1
                     (a2*exp(-((wk-b2)^2)/(2*c2^2))) + # wave 2
                     (a3*exp(-((wk-b3)^2)/(2*c3^2))) + # wave 3
                     (a4*exp(-((wk-b4)^2)/(2*c4^2))) + # wave 4 
                     (a5*exp(-((wk-b5)^2)/(2*c5^2))) + # wave 5
                     (a6*exp(-((wk-b6)^2)/(2*c6^2))) + # wave 6
                     (a7*exp(-((wk-b7)^2)/(2*c7^2))), # wave 7
              data = df,
              start = list(a1=val.peaks[1], b1=wk.peaks[1], c1=1.4,
                           a2=val.peaks[2], b2=wk.peaks[2], c2=5.5,
                           a3=val.peaks[3], b3=wk.peaks[3], c3=3,
                           a4=val.peaks[4], b4=wk.peaks[4], c4=4,
                           a5=val.peaks[5], b5=wk.peaks[5], c5=4,
                           a6=val.peaks[6], b6=wk.peaks[6], c6=3,
                           a7=val.peaks[7], b7=wk.peaks[7], c7=3),
              # trace = TRUE
)
(model.fit <- summary(cov.mod))                                                 

# Predict the fitted model 
dffit <- data.frame(wk=seq(0, 160, .1))
dffit$val <- predict(cov.mod, newdata=dffit)

# Plot
print(ggplot(df, aes(wk, val)) + geom_point() +
        geom_smooth(data=dffit, stat="identity", color="red")) 

# Model coefficients for fct_addInc
coef <- as.data.frame(model.fit$coefficients[,1])
peaks.height <- coef[c(1,4,7,10,13,16,19),]
peaks.position <- coef[c(2,5,8,11,14,17,20),]
sd <- coef[c(3,6,9,12,15,18,21),]

# Full with at tenth maximum (FWTM)
fct_fwtm <- function(x){
  2*sqrt(2*log(10))*x
}

fwtm <- fct_fwtm(sd)

# Full with at half maximum (FWHM)
fct_fwhm <- function(x){
  2*sqrt(2*log(2))*x
}

fwhm <- fct_fwhm(sd)
