install.packages("readcsvcolumns", repos="http://193.190.10.42/jori/")
install.packages("tidyverse")
install.packages("data.table")
install.packages(c("exactci",
                   "fBasics",
                   "VarianceGamma",
                   "fitdistrplus",
                   "EnvStats",
                   "gtools",
                   "lhs",
                   "GGally",
                   "lmtest",
                   "mice",
                   "miceadds",
                   "randtoolbox",
                   "EasyABC"))

library(RSimpactCyan)
library(RSimpactHelper)
library(ggplot2)
library(GGally)
library(lmtest)
library(mice)
library(miceadds)
library(parallel)
library(randtoolbox)
library(RSimpactHelper)
library(EasyABC)
library(dplyr)

# First we set up a configuration of Simpact parameters to run and get target statistics out

