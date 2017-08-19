setwd("/Users/wimdelva/Documents/MiceABC/R")
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
#library(depth)
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
x.offset <- length(dummy.input.vector)
n.experiments <- 10000

dummy.master <- simpact.parallel(model = dummy.wrapper,
                             actual.input.matrix = matrix(rep(dummy.input.vector, each = n.experiments), nrow = n.experiments),
                             seed_count = 0,
                             n_cluster = 8)
dummy.master <- dummy.master %>%
  as.data.frame() %>%
  dplyr::filter(complete.cases(.))


dummy.targets.empirical <- l1median(dummy.master)


predictorMatrix <- (1 - diag(1, length(c(dummy.input.vector, dummy.targets.empirical)))) # This is the default matrix.
# # Let's now modify the first 15 rows of this matrix, corresponding to the indicators of predictor variables for the input variables. In brackets the values for the master model.

predictorMatrix[1:x.offset, ] <- 0 # First we "empty" the relevant rows, then we refill them.
# We are currently not allowing input variables to be predicted by other predictor variables. Only via output variables. We could change this at a later stage.
predictorMatrix[1, x.offset + c(3, 4)] <- 1
predictorMatrix[2, x.offset + 2:4] <- 1
predictorMatrix[3, x.offset + 2:4] <- 1
predictorMatrix[4, x.offset + 2:4] <- 1
predictorMatrix[5, c(2, 3)] <- 1
predictorMatrix[6, 2:4] <- 1
predictorMatrix[7, 2:4] <- 1
predictorMatrix[8, 2:4] <- 1


test.MaC.incremental <- dummy.MaC.incremental(targets.empirical = dummy.targets.empirical,
                                        RMSD.tol.max = 0.5,
                                        min.givetomice = 200, # 400
                                        n.experiments = 800, # 1000
                                        lls = c(590,  0, 0, 2),
                                        uls = c(2000, 4, 6, 6),
                                        model = dummy.wrapper,
                                        strict.positive.params = c(1:4),
                                        predictorMatrix = predictorMatrix,
                                        maxit = 20,
                                        maxwaves = 20) # 6
(round(l1median(head(test.MaC.incremental$selected.experiments[[length(test.MaC.incremental$selected.experiments)]]), 1), 99)[5:8] - dummy.targets.empirical[1:4]) / dummy.targets.empirical[1:4]
round(l1median(head(test.MaC.incremental$selected.experiments[[length(test.MaC.incremental$selected.experiments)]]), 1), 2)
test.MaC.incremental$secondspassed
