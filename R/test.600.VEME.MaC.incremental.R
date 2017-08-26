#setwd("/Users/wimdelva/Documents/MiceABC/R")
#source("simpact.wrapper.R")
#source("VEME.wrapper.R")
#source("dummy.wrapper.R")
#source("simpact.parallel.R")
#source("MaC.incremental.R")
#source("dummy.MaC.incremental.R")

#source("/user/data/gent/vsc400/vsc40070/phylo/scripts/VEME.wrapper.R")
#source("/user/data/gent/vsc400/vsc40070/phylo/scripts/simpact.parallel.R")
#source("/user/data/gent/vsc400/vsc40070/phylo/scripts/dummy.MaC.incremental.R")

source("/Users/delvaw/Documents/MiceABC/R/VEME.wrapper.R")
source("/Users/delvaw/Documents/MiceABC/R/VEME.wrapper2.R")
source("/Users/delvaw/Documents/MiceABC/R/00-Functions.R")

source("/Users/delvaw/Documents/MiceABC/R/simpact.parallel.R")
source("/Users/delvaw/Documents/MiceABC/R/dummy.MaC.incremental.R")

library(dplyr)
library(MASS)
library(splines)
library(boot)
#library(haven)
library(ggplot2)
library(GGally)
library(fitdistrplus)
library(lmtest)
library(mclust)
#library(depth)
library(pcaPP)
#library(devtools)
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


dummy.input.vector <- c(1.1, 0.25, 0, 3, 0.23, 0.23,
                        45, 45, #45, 45,
                        -0.5, 2.8, -0.2, -0.2, -2.5, -0.52, -0.05)# c(1000, 2, 3, 4)
x.offset <- length(dummy.input.vector)
n.experiments <- 120
dummy.master2 <- simpact.parallel(model = VEME.wrapper2,
                             actual.input.matrix = matrix(rep(dummy.input.vector, each = n.experiments), nrow = n.experiments),
                             seed_count = 0,
                             n_cluster = 8)

#####
# The output of the master model
#####
inc.master.vector <- dummy.master2[, 17]
#save(inc.master.vector, file = "/Users/delvaw/Documents/MiceABC/inc.master.vector.RData")
hist(inc.master.vector, 8) # The distribution of HIV incidence after 10 years
mean(inc.master.vector) # The mean HIV incidence after 10 years
median(inc.master.vector) # The median HIV incidence after 10 years
quantile(inc.master.vector, c(0.025, 0.975))

new.infect.vector <- dummy.master2[, 18]
hist(new.infect.vector[!is.na(new.infect.vector) & new.infect.vector < Inf], 20)

plot(new.infect.vector[!is.na(new.infect.vector) & new.infect.vector < Inf],
     inc.master.vector[!is.na(new.infect.vector) & new.infect.vector < Inf])

cor.test(inc.master.vector[!is.na(new.infect.vector) & new.infect.vector < Inf],
         new.infect.vector[!is.na(new.infect.vector) & new.infect.vector < Inf])

    sim <- sequence.simulation(transtree = tree0, seedSeq = hiv.seq.env.short, alpha = 0.90,
                               rate.list = rate, base.freq = freq)
    saveAlignment.PhyloSim(sim,
                           file = paste("/Users/delvaw/Documents/MiceABC/HIVSeq_fullNetwork.fasta",
                                        sep = ""),
                           skip.internal = TRUE,
                           paranoid = TRUE)
    # to handle sequences, let use phangorn and ape packages and read the sequence data and build the phylogenetic tree
    seq.sim.size_full <- read.FASTA("/Users/delvaw/Documents/MiceABC/HIVSeq_fullNetwork.fasta")
    tree.dat.full <- phyDat(seq.sim.size_full, type = "DNA")
    tree.ml.full <- dist.ml(tree.dat.full)
    tree.sim.full <- upgma(tree.ml.full)



dummy.master2 <- dummy.master2 %>%
  as.data.frame() %>%
  dplyr::filter(complete.cases(.))


dummy.targets.empirical <- l1median(dummy.master2)


predictorMatrix <- (1 - diag(1, length(c(dummy.input.vector, dummy.targets.empirical)))) # This is the default matrix.
# # Let's now modify the first 15 rows of this matrix, corresponding to the indicators of predictor variables for the input variables. In brackets the values for the master model.

predictorMatrix[1:x.offset, ] <- 0 # First we "empty" the relevant rows, then we refill them.
# We are currently not allowing input variables to be predicted by other predictor variables. Only via output variables. We could change this at a later stage.

predictorMatrix[1, x.offset + 10:11] <- 1 # relative susceptibility in young women is predicted by HIV prevalence in young men and women
      predictorMatrix[2, x.offset + 3] <- 1 # agescale predicted by slope
      predictorMatrix[3, x.offset + c(1, 3, 6)] <- 1 # mean of the person-specific age gap preferences is predicted by slope, intercept and AAD
      predictorMatrix[4, x.offset + c(2, 4, 5)] <- 1 # sd of the person-specific age gap preferences is predicted by SD, WSD, BSD
      predictorMatrix[5, x.offset + c(7, 8, 9, 13, 16)] <- 1 # man gamma a predicted by gamma shape.male, scale.male, pp.cp, hiv.prev.25.34.men, exp(growthrate)
      predictorMatrix[6, x.offset + c(7, 8, 9, 12, 16)] <- 1 # woman gamma a predicted by gamma shape.male, scale.male, pp.cp, hiv.prev.25.34.women, exp(growthrate)
      predictorMatrix[7, x.offset + c(7, 8, 9, 13, 16)] <- 1 # man gamma b predicted by gamma shape.male, scale.male, pp.cp, hiv.prev.25.34.men, exp(growthrate)
      predictorMatrix[8, x.offset + c(7, 8, 9, 12, 16)] <- 1 # woman gamma b predicted by gamma shape.male, scale.male, pp.cp, hiv.prev.25.34.men, exp(growthrate)
      predictorMatrix[9, x.offset + c(2, 4, 5, 7, 8, 14, 15, 16)] <- 1 # formation.hazard.agegapry.gap_factor_x_exp is predicted by population growth, age gap variance, hiv prevalence,
      predictorMatrix[10, x.offset + c(7, 8, 9, 12, 13, 16)] <- 1 # baseline formation hazard predicted by HIV prevalence, cp, degree distrib. HIV prevalence.
      predictorMatrix[11, x.offset + c(7, 8, 9, 12, 13, 16)] <- 1 # numrel man penalty is predicted by degree distrib, cp, prev, popgrowth
      predictorMatrix[12, x.offset + c(7, 8, 9, 12, 13, 16)] <- 1 # # numrel woman penalty is predicted by degree distrib, cp, prev, popgrowth
      predictorMatrix[13, x.offset + 16] <- 1 # conception.alpha_base is predicted by popgrowth
      predictorMatrix[14, x.offset + c(7, 8, 9, 16)] <- 1 # baseline dissolution hazard predicted by degree distrib, cp, popgrowth
      predictorMatrix[15, x.offset + c(7, 8, 9, 16)] <- 1 # age effect on dissolution hazard predicted by degree distrib, cp, popgrowth, HIV prev in older people (maybe?)

      # NOTE: As it stands, each output statistic is predicted by ALL input and ALL other output statistics. That may not be a great idea, or even possible, if there is collinearity.

test.600.VEME.MaC.incremental <- dummy.MaC.incremental(targets.empirical = dummy.targets.empirical,
                                        RMSD.tol.max = 0.95,
                                        min.givetomice = 32, # 400
                                        n.experiments = 96, # 1000
                                        lls = c(1,  0.12, -0.3, 2.5, 0.1, 0.1, 20, 20, -0.8, 2, -0.35, -0.35, -3.6, -0.8, -0.16),
                                        uls = c(1.2, 0.37, 0.3, 3.5, 0.4, 0.4, 66, 66, -0.25, 3.9, -0.1, -0.1, -1.4, -0.3,  -0.001),
                                        model = VEME.wrapper,
                                        strict.positive.params = c(4:8),
                                        predictorMatrix = predictorMatrix,
                                        maxit = 10,
                                        maxwaves = 2,
                                        n_cluster = 8) # 6
#(round(l1median(head(test.MaC.incremental$selected.experiments[[length(test.MaC.incremental$selected.experiments)]]), 1), 99)[5:8] - dummy.targets.empirical[1:4]) / dummy.targets.empirical[1:4]
#round(l1median(head(test.MaC.incremental$selected.experiments[[length(test.MaC.incremental$selected.experiments)]]), 1), 2)
#test.MaC.incremental$secondspassed

test.600.VEME.MaC.incremental$secondspassed
test.600.VEME.MaC.incremental$max.RMSD
test.600.VEME.MaC.incremental$n.close.to.targets
head(test.600.VEME.MaC.incremental$selected.experiments[[length(test.600.VEME.MaC.incremental$selected.experiments)]])

save(dummy.master, dummy.targets.empirical, test.600.VEME.MaC.incremental, file = "/user/data/gent/vsc400/vsc40070/phylo/output/test.600.VEME.MaC.incremental.RData")
