rm(list = ls(all=TRUE))

# Load packages
library(here)
library(splines)
library(tidyverse)
library(ggpmisc)
library(gt)

# Use this to toggle between sexes
sexString <- "Female" #"Female" or "Male"
sexString2 <- "w" #"w" or "m"

# Aux functions
logit <- function(x){
  log(x/(1-x))
}

expit <- function(x){
  exp(x)/(1+exp(x))
}

# Prepare data set -------------------------------------------------------------
# Read in GBD 2019 data set (comma-separated)
gbd <- read.csv(here("raw_data", "IHME-GBD_2019_DATA-69583287-1", "IHME-GBD_2019_DATA-69583287-1.csv"),
                stringsAsFactors = TRUE)
str(gbd)

# Check age groups in the data set
table(gbd$age_name)

# Add numeric age variable according to age groups
ageVals <- c(3, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 7.5, #          
             52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 82.5, 87.5, 92.5)
age <- ageVals[as.numeric(gbd$age_name)]
gbd$age <- age
# Check
with(gbd, table(age, age_name, useNA = "alw"))

# Read in population projection from Destatis
pop <- read.csv2(here("raw_data", "15_bevoelkerungsvorausberechnung_daten.csv"), header = TRUE) # NEW 15th coordinated projection
str(pop)
