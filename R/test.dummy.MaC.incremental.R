setwd("/Users/delvaw/Documents/MiceABC/R")
source("simpact.wrapper.R")
source("dummy.wrapper.R")
source("simpact.parallel.R")
source("MaC.incremental.R")
source("dummy.MaC.incremental.R")


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
library(pcaPP)
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


dummy.input.vector <- c(1000, 2, 3, 4)
n.experiments <- 10000

dummy.master <- simpact.parallel(model = dummy.wrapper,
                             actual.input.matrix = matrix(rep(dummy.input.vector, each = n.experiments), nrow = n.experiments),
                             seed_count = 0,
                             n_cluster = 8)
dummy.master <- dummy.master %>%
  as.data.frame() %>%
  dplyr::filter(complete.cases(.))


dummy.targets.empirical <- l1median(dummy.master)


test.MaC.incremental <- dummy.MaC.incremental(targets.empirical = dummy.targets.empirical,
                                        RMSD.tol.max = 0.5,
                                        min.givetomice = 200, # 400
                                        n.experiments = 2000, # 1000
                                        lls = c(590,  0, 0, 2),
                                        uls = c(2000, 4, 6, 6),
                                        model = dummy.wrapper,
                                        maxit = 20,
                                        maxwaves = 10) # 6
(round(l1median(head(test.MaC.incremental$selected.experiments[[length(test.MaC.incremental$selected.experiments)]]), 1), 99)[5:8] - dummy.targets.empirical[1:4]) / dummy.targets.empirical[1:4]
round(l1median(head(test.MaC.incremental$selected.experiments[[length(test.MaC.incremental$selected.experiments)]]), 1), 2)
round(l1median(head(test.MaC.incremental$selected.experiments[[length(test.MaC.incremental$selected.experiments)]]), 20), 2)[1:4]
