setwd("/Users/wimdelva/Documents/MiceABC/R")
source("simpact.wrapper.R")
source("simpact.parallel.R")
source("MaC.incremental.R")

targets.empirical <- c(4, 4, 0.7, 3, 1.5, -1, 1.29, 0.66, 0.134, 0.2, 0.1, 0.42, 0.17, 0.37, 0.24, 1.01)



library(dplyr)
library(MASS)
library(splines)
library(boot)
library(haven)
library(ggplot2)
library(GGally)
library(fitdistrplus)
library(lmtest)
library(mclust)
library(depth)
library(devtools)
#install_github("wdelva/RSimpactHelp")
library(RSimpactCyan)
library(RSimpactHelper)
library(lmtest)
library(mice)
#library(miceadds)
library(parallel)
library(randtoolbox)
library(EasyABC)
library(dplyr)
library(tidyr)
library(nlme)
library(lme4)
library(boot)
library(data.table)

test.MaC.incremental <- MaC.incremental(targets.empirical = targets.empirical,
                                                RMSD.tol.max = 1,
                                                min.givetomice = 64,
                                                n.experiments = 256,
                                                lls = c(1,  0.12, -0.3, 2.5, 0.1, 0.1, 20, 20, -0.8, 2, -0.35, -0.35, -3.6, -0.8, -0.16),
                                                uls = c(1.2, 0.37, 0.3, 3.5, 0.4, 0.4, 66, 66, -0.25, 3.9, -0.1, -0.1, -1.4, -0.3,  -0.001),
                                                model = simpact.wrapper,
                                                maxit = 50,
                                                maxwaves = 8)
