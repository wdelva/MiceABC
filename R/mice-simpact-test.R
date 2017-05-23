## OLD verison before SIMID visit

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




# A wrapper function to run Simpact once
simpact.wrapper <- function(inputvector = input.vector){
  library(RSimpactHelper)
  library(RSimpactCyan)

  destDir <- "/Users/wimdelva/Downloads"
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.21,#1,
                                   population.simtime = 40,
                                   population.nummen = 2500,
                                   population.numwomen = 2500,
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

  cfg <- cfg.list
  seedid <- inputvector[1]
  cfg["person.agegap.man.dist.fixed.value"] <- inputvector[2]
  cfg["person.agegap.woman.dist.fixed.value"] <- inputvector[2]
  cfg["formation.hazard.agegapry.gap_factor_man_exp"] <- inputvector[3]
  cfg["formation.hazard.agegapry.gap_factor_woman_exp"] <- inputvector[3]
  results <- simpact.run(configParams = cfg,
                         destDir = destDir,
                         agedist = age.distr,
                         seed = seedid)
  datalist.agemix <- readthedata(results)
  agegap.mean <-mean(datalist.agemix$rtable$AgeGap)
  agegap.sd <- sd(datalist.agemix$rtable$AgeGap)
  outputvector <- c(agegap.mean, agegap.sd)
  return(outputvector)
}

# A function to run Simpact in parallel
simpact.parallel <- function(model = simpact.wrapper,
                             actual.input.matrix = matrix(rep(c(-6, -6), 16), nrow = 16),
                             #nb_simul = 16,
                             seed_count = 0,
                             n_cluster = 8){
  cl <- makeCluster(getOption("cl.cores", n_cluster))
  tab_simul_summarystat = NULL
  list_param <- list(NULL)
  tab_param <- NULL
  paramtemp <- NULL
  simultemp <- NULL

  nb_simul <- nrow(actual.input.matrix)

  for (i in 1:nb_simul) {
    l <- ncol(actual.input.matrix)
    param <- c((seed_count + i), actual.input.matrix[i, ])
    list_param[[i]] <- param
    tab_param <- rbind(tab_param, param[2:(l + 1)])
    paramtemp <- rbind(paramtemp, param[2:(l + 1)])
  }
  list_simul_summarystat = parLapplyLB(cl, list_param,
                                       model)
  tab_simul_summarystat <- do.call(rbind, list_simul_summarystat)
  stopCluster(cl)
  return(tab_simul_summarystat)
}



sumsmstatsvect <- simpact.wrapper(inputvector = c(1, -6, -6))
sumsmstatsvect # Mean age gap of 4.3 years and sd of 4.6 years

sumsmstats.df <- simpact.parallel(model = simpact.wrapper,
                                  actual.input.matrix = matrix(rep(c(-6, -6), 32), nrow = 32),
                                  seed_count = 0,
                                  n_cluster = 8)
summary(sumsmstats.df)
# Min.   :3.419   Min.   :4.265
# 1st Qu.:3.685   1st Qu.:4.387
# Median :3.743   Median :4.482
# Mean   :3.750   Mean   :4.479
# 3rd Qu.:3.875   3rd Qu.:4.528
# Max.   :3.971   Max.   :4.684

miceabc.output <- MiceABC(maxwaves = 1)

miceabc.output2 <- MiceABC(targets = c(3.75, 4.48),
                          n.experiments = 64,
                          lls = c(-7, -7),
                          uls = c(-2, -2),
                          model = simpact.wrapper,
                          maxit = 20,
                          maxwaves = 5,#0,
                          #reps = 5,
                          alpha = 0.25,
                          saturation.crit = 0)


save(simpact.wrapper, simpact.parallel, MiceABC, sumsmstats.df, miceabc.output, miceabc.output2, rejectabc.output, lenormandabc.output, file = "/Users/wimdelva/Documents/miceabc.comparison.RData")

# Plotting all the input parameter sets and then overlaying this with the final calibration results (16 best)
all.results.with.design <- do.call(rbind, miceabc.output$sim.results.with.design)

plot(miceabc.output$inputoutput[[5]][,1],
     miceabc.output$inputoutput[[5]][,2],
     xlab = "x1",
     ylab = "x2",
     xlim = c(-7, -2),
     ylim = c(-7, -2),
     pch = 20)
points(-6, -6, pch = 3, col = "red3", lwd = 3)


plot(miceabc.output$inputoutput[[5]][,3],
     miceabc.output$inputoutput[[5]][,4],
     xlab = "y1",
     ylab = "y2",
     xlim = c(3.5, 4),
     ylim = c(4.25, 4.75),
     pch = 20)
points(3.75, 4.48, pch = 3, col = "red3", lwd = 3)




# As comparison, let's run a basic ABC accept rejection scheme
help("ABC_rejection")
rejectabc.output <- ABC_rejection(model = simpact.wrapper,
                                  prior = list(c("unif", -7, -2),
                                               c("unif", -7, -2)),
                                  nb_simul = 250,
                                  #summary_stat_target = c(3.75, 4.48),   # this and tol is commented out, so that we can see ALL the output, not just the retained fraction.
                                  #tol = 0.064, # That is 16/250 to be in line with the output of miceabc
                                  use_seed = TRUE,
                                  seed_count = 0,
                                  n_cluster = 8)
# Calculating squared relative distances and sum thereof
targets <- c(3.75, 4.48)
rejectabc.output$y.1.sq.rel.dist <- ((rejectabc.output$stats[,1] - targets[1]) / targets[1])^2
rejectabc.output$y.2.sq.rel.dist <- ((rejectabc.output$stats[,2] - targets[2]) / targets[2])^2
rejectabc.output$sum.sq.rel.dist <- rejectabc.output$y.1.sq.rel.dist + rejectabc.output$y.2.sq.rel.dist

abcreject.sim.results.with.design.df <- as.data.frame(cbind(rejectabc.output$stats, rejectabc.output$y.1.sq.rel.dist, rejectabc.output$y.2.sq.rel.dist, rejectabc.output$sum.sq.rel.dist))
names(abcreject.sim.results.with.design.df) <- c("y1", "y2", "y.1.sq.rel.dist", "y.2.sq.rel.dist", "sum.sq.rel.dist")
# Now subsetting the 16 shortest distances
dist.order <- order(abcreject.sim.results.with.design.df$sum.sq.rel.dist) # Ordering the squared distances from small to big. The last observation (targets) should be ranked first
last.one.selected <- 16
selected.distances <- dist.order[1:last.one.selected]
abcreject.sim.results.with.design.df$shortest16 <- 0
abcreject.sim.results.with.design.df$shortest16[selected.distances] <- 1


# Plotting all the input parameter sets and then overlaying this with the final calibration results (16 best)
plot(abcreject.sim.results.with.design.df$y1,
     abcreject.sim.results.with.design.df$y2,
     xlab = "y1",
     ylab = "y2",
     xlim = c(min(abcreject.sim.results.with.design.df$y1), max(abcreject.sim.results.with.design.df$y1)),
     ylim = c(4,5.2), #c(min(abcreject.sim.results.with.design.df$y2), max(abcreject.sim.results.with.design.df$y2)),
     pch = 20 - abcreject.sim.results.with.design.df$shortest16,
     col = abcreject.sim.results.with.design.df$shortest16 * 5 + 1,
     cex = 0.75)

# And overlaying the target statistics
points(3.75, 4.48, pch = 3, col = "blue3", lwd = 3, cex = 2)

# And now we overlay the best 16 from rejectabc with the 16 best from miceabc
plot(abcreject.sim.results.with.design.df$y1[abcreject.sim.results.with.design.df$shortest16 == 1],
     abcreject.sim.results.with.design.df$y2[abcreject.sim.results.with.design.df$shortest16 == 1],
     xlab = "y1",
     ylab = "y2",
     #xlim = c(min(abcreject.sim.results.with.design.df$y1), max(abcreject.sim.results.with.design.df$y1)),
     #ylim = c(4,5.2), #c(min(abcreject.sim.results.with.design.df$y2), max(abcreject.sim.results.with.design.df$y2)),
     pch = 20, # - abcreject.sim.results.with.design.df$shortest16,
     col = 6, #abcreject.sim.results.with.design.df$shortest16 * 5 + 1,
     cex = 1)

points(miceabc.output$inputoutput[[5]][,3],
     miceabc.output$inputoutput[[5]][,4],
     xlab = "y1",
     ylab = "y2",
     pch = 20,
     col = "blue2",
     cex = 1)

# And overlaying the target statistics
points(3.75, 4.48, pch = 3, col = "blue3", lwd = 3, cex = 2)


# And the input parameters
plot(rejectabc.output$param[abcreject.sim.results.with.design.df$shortest16 == 1, 1],
     rejectabc.output$param[abcreject.sim.results.with.design.df$shortest16 == 1, 2],
     xlab = "x1",
     ylab = "x2",
     xlim = c(-7, -2), # c(min(abcreject.sim.results.with.design.df$y1), max(abcreject.sim.results.with.design.df$y1)),
     ylim = c(-7, -2), #c(min(abcreject.sim.results.with.design.df$y2), max(abcreject.sim.results.with.design.df$y2)),
     pch = 20, # - abcreject.sim.results.with.design.df$shortest16,
     col = 6, #abcreject.sim.results.with.design.df$shortest16 * 5 + 1,
     cex = 1)
points(miceabc.output$inputoutput[[5]][,1],
            miceabc.output$inputoutput[[5]][,2],
            xlab = "x1",
            ylab = "x2",
       col = "blue2",
       pch = 20)
points(-6, -6, pch = 3, col = "darkgrey", lwd = 3, cex = 2)


# Second comparison: let's run Lenormand's adaptive population ABC scheme
help("ABC_sequential")
lenormandabc.output <- ABC_sequential(method = "Lenormand",
                                      model = simpact.wrapper,
                                      prior = list(c("unif", -7, -2),
                                                   c("unif", -7, -2)),
                                      nb_simul = 64,
                                      summary_stat_target = c(3.75, 4.48),   # this and tol is commented out, so that we can see ALL the output, not just the retained fraction.
                                      alpha = 0.25,
                                      p_acc_min = 0.05,
                                      use_seed = TRUE,
                                      seed_count = 0,
                                      n_cluster = 8)
