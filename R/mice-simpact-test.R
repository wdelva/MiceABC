# Install packages if not already done so
# install.packages("readcsvcolumns", repos="http://193.190.10.42/jori/")
# install.packages("tidyverse")
# install.packages("data.table")
# install.packages(c("exactci",
#                    "fBasics",
#                    "VarianceGamma",
#                    "fitdistrplus",
#                    "EnvStats",
#                    "gtools",
#                    "lhs",
#                    "GGally",
#                    "lmtest",
#                    "mice",
#                    "miceadds",
#                    "randtoolbox",
#                    "EasyABC"))

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

# First we define a simpact wrapper function that takes as argument a vector of model parameter values and returns a vector of summary statistics






simpact.wrapper <- function(inputvector = input.vector){
  cfg <- cfg.list
  cfg["person.agegap.man.dist.fixed.value"] <- inputvector[1]
  cfg["person.agegap.woman.dist.fixed.value"] <- inputvector[1]
  cfg["formation.hazard.agegapry.gap_factor_man_exp"] <- inputvector[2]
  cfg["formation.hazard.agegapry.gap_factor_woman_exp"] <- inputvector[2]
  results <- simpact.run(configParams = cfg,
                         destDir = destDir,
                         agedist = age.distr)
  datalist.agemix <- readthedata(results)
  agegapsd <- sd(datalist.agemix$rtable$AgeGap)
  agegapmean <-mean(datalist.agemix$rtable$AgeGap)
  outputvector <- c(agegapmean, agegapsd)
  return(outputvector)
}


sumsmstatsvect <- simpact.wrapper(inputvector = c(-6, -6))


# Next we set up a configuration of Simpact parameters to run and get target statistics out
destDir <- "/Users/wimdelva/Downloads"
age.distr <- agedistr.creator(shape = 5, scale = 65)
cfg.list <- input.params.creator(population.eyecap.fraction = 0.21,#1,
                                 population.simtime = 40,
                                 population.nummen = 5000,
                                 population.numwomen = 5000,
                                 hivseed.time = 10,
                                 hivseed.type = "amount",
                                 hivseed.amount = 30,
                                 hivseed.age.min = 20,
                                 hivseed.age.max = 25,
                                 hivtransmission.param.a = -1,
                                 hivtransmission.param.b = -90,
                                 hivtransmission.param.c = 0.5,
                                 hivtransmission.param.f1 = log(2),
                                 hivtransmission.param.f2 = log(log(1.4) / log(2)) / 5,
                                 formation.hazard.agegapry.gap_factor_man_age = -0.072653928518528523251061,
                                 formation.hazard.agegapry.gap_factor_woman_age = -0.0726539285185285232510561,
                                 formation.hazard.agegapry.meanage = -0.05,
                                 formation.hazard.agegapry.gap_factor_man_const = 0,
                                 formation.hazard.agegapry.gap_factor_woman_const = 0,
                                 formation.hazard.agegapry.gap_factor_man_exp = -6,#-1.5,
                                 formation.hazard.agegapry.gap_factor_woman_exp = -6,#-1.5,
                                 formation.hazard.agegapry.gap_agescale_man = 0.3,
                                 formation.hazard.agegapry.gap_agescale_woman = 0.3,#-0.30000007,#-0.03,
                                 debut.debutage = 15,
                                 person.art.accept.threshold.dist.fixed.value = 0
)
cfg.list["formation.hazard.agegapry.baseline"] <- 3.5
cfg.list["mortality.aids.survtime.C"] <- 65
cfg.list["mortality.aids.survtime.k"] <- -0.2

cfg.list["person.agegap.man.dist.type"] <- "fixed"
cfg.list["person.agegap.man.dist.fixed.value"] <- -6
cfg.list["person.agegap.woman.dist.type"] <- "fixed"
cfg.list["person.agegap.woman.dist.fixed.value"] <- -6

cfg.list["mortality.aids.survtime.C"] <- 65
cfg.list["mortality.aids.survtime.k"] <- -0.2
cfg.list["monitoring.cd4.threshold"] <- NULL
cfg.list["person.agegap.man.dist.normal.mu"] <- NULL
cfg.list["person.agegap.woman.dist.normal.mu"] <- NULL
cfg.list["person.agegap.man.dist.normal.sigma"] <- NULL
cfg.list["person.agegap.woman.dist.normal.sigma"] <- NULL

cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1

