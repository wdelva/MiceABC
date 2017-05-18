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

## Code as copied from the iMac before leaving for Belgium, on 14 May 2017

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
library(nlme)


# Let's create a colour palette that is colourblind-friendly
# The palette with black:
cb.Palette <- c("#F0E442", "#E69F00", "#0072B2", "#D55E00", "#999999")
#                yellow     orange      blue       red        grey

# First we define a simpact wrapper function that takes as argument a vector of model parameter values and returns a vector of summary statistics

# A wrapper function to run Simpact once
input.vector <- c(0, 3, -0.5)
simpact.wrapper <- function(inputvector = input.vector){
  library(RSimpactHelper)
  library(RSimpactCyan)

  destDir <- "/Users/delvaw/Downloads"
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.21,#1,
                                   population.simtime = 20,
                                   population.nummen = 2500,
                                   population.numwomen = 2500,
                                   hivseed.time = 10,
                                   hivseed.type = "amount",
                                   hivseed.amount = 30,
                                   hivseed.age.min = 20,
                                   hivseed.age.max = 50,
                                   hivtransmission.param.a = -1,
                                   hivtransmission.param.b = -90,
                                   hivtransmission.param.c = 0.5,
                                   hivtransmission.param.f1 = log(2),
                                   hivtransmission.param.f2 = log(log(1.4) / log(2)) / 5,
                                   formation.hazard.agegapry.gap_factor_man_age = -0.01, #-0.01472653928518528523251061,
                                   formation.hazard.agegapry.gap_factor_woman_age = -0.01, #-0.0726539285185285232510561,
                                   formation.hazard.agegapry.meanage = -0.025,
                                   formation.hazard.agegapry.gap_factor_man_const = 0,
                                   formation.hazard.agegapry.gap_factor_woman_const = 0,
                                   formation.hazard.agegapry.gap_factor_man_exp = -1, #-6,#-1.5,
                                   formation.hazard.agegapry.gap_factor_woman_exp = -1, #-6,#-1.5,
                                   formation.hazard.agegapry.gap_agescale_man = 0.25,
                                   formation.hazard.agegapry.gap_agescale_woman = 0.25,#-0.30000007,#-0.03,
                                   debut.debutage = 15,
                                   conception.alpha_base = -2.5#,
                                   #person.art.accept.threshold.dist.fixed.value = 0
  )
  # cfg.list["diagnosis.baseline"] <- 100
  # cfg.list["diagnosis.HSV2factor"] <- 0
  # cfg.list["diagnosis.agefactor"] <- 0
  # cfg.list["diagnosis.beta"] <- 0
  # cfg.list["diagnosis.diagpartnersfactor"] <- 0
  # cfg.list["diagnosis.genderfactor"] <- 0
  # cfg.list["diagnosis.isdiagnosedfactor"] <- 0
  # cfg.list["diagnosis.t_max"] <- 200
  cfg.list["formation.hazard.agegapry.baseline"] <- 2
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.fraction.log_viralload"] <- 0.3
  cfg.list["dropout.interval.dist.uniform.min"] <- 100
  cfg.list["dropout.interval.dist.uniform.max"] <- 200

  cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
  #cfg.list["person.agegap.man.dist.fixed.value"] <- -6
  cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
  #cfg.list["person.agegap.woman.dist.fixed.value"] <- -6

  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  #cfg.list["monitoring.cd4.threshold"] <- 0
  cfg.list["person.agegap.man.dist.normal.mu"] <- 0
  cfg.list["person.agegap.woman.dist.normal.mu"] <- 0
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[2] ######### 3
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[2] ######### 3

  cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
  cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
  cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1

  cfg <- cfg.list
  seedid <- inputvector[1]
  #cfg["person.agegap.man.dist.fixed.value"] <- -2 # inputvector[2]
  #cfg["person.agegap.woman.dist.fixed.value"] <- -2 # inputvector[2]
  cfg["formation.hazard.agegapry.gap_factor_man_exp"] <- inputvector[3] ######### -0.5
  cfg["formation.hazard.agegapry.gap_factor_woman_exp"] <- inputvector[3] ######### -0.5

  # Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly
  art.intro <- list()
  art.intro["time"] <- 25
  art.intro["person.art.accept.threshold.dist.fixed.value"] <- 0.75
  art.intro["diagnosis.baseline"] <- 100
  art.intro["monitoring.cd4.threshold"] <- 100 # 1200
  #art.intro["monitoring.interval.piecewise.cd4s"] <- "0,1300"


  # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500

  art.intro2 <- list()
  art.intro2["time"] <- 30
  art.intro2["monitoring.cd4.threshold"] <- 200

  art.intro3 <- list()
  art.intro3["time"] <- 33
  art.intro3["monitoring.cd4.threshold"] <- 350

  art.intro4 <- list()
  art.intro4["time"] <- 36
  art.intro4["monitoring.cd4.threshold"] <- 500



  intervention <- list(art.intro, art.intro2, art.intro3, art.intro4)


  results <- simpact.run(configParams = cfg,
                         destDir = destDir,
                         agedist = age.distr,
                         seed = seedid, #, Introducing ART has helped to keep the prevalence high
                         intervention = intervention)

  datalist.agemix <- readthedata(results)
  agemix.df <- agemix.df.maker(datalist.agemix)

  agemix.model <- pattern.modeller(dataframe = agemix.df,
                                   agegroup = c(15, 50),
                                   timepoint = datalist.agemix$itable$population.simtime[1],
                                   timewindow = 3)

  men.lme <- lme(pagerelform ~ agerelform0,
                 data = dplyr::filter(agemix.model[[1]], Gender =="male"),
                 control=lmeControl(maxIter=200, returnObject=TRUE, opt = "optim"),
                 random = ~1 | ID,
                 method = "REML",
                 weight = varPower(value = 0.5, form = ~agerelform0 + 1))

  AAD <- mean(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap)
  SDAD <- sd(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap)
  powerm <- as.numeric(attributes(men.lme$apVar)$Pars["varStruct.power"])
  slopem <- summary(men.lme)$tTable[2, 1]
  WVAD.base <- men.lme$sigma^2
  BVAD = getVarCov(men.lme)[1,1]

  #agegap.mean <- mean(datalist.agemix$rtable$AgeGap)
  #agegap.sd <- sd(datalist.agemix$rtable$AgeGap)
  hivprev.15.50 <- prevalence.calculator(datalist = datalist.agemix,
                        agegroup = c(15, 50),
                        timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[3]
  growthrate <- pop.growth.calculator(datalist = datalist.agemix,
                        timewindow = c(0, datalist.agemix$itable$population.simtime[1]))

  outputvector <- c(AAD, SDAD, powerm, slopem, WVAD.base, BVAD, hivprev.15.50, growthrate)
  outputvector
  return(outputvector)
}

# Visually checking the age-mixing pattern
agemix.df <- agemix.df.maker(datalist.agemix)
# Scatter plot of age-mixing pattern
ggplot(data = dplyr::filter(agemix.df, Gender == "male", agerelform >=18, agerelform < 50),
       aes(x = agerelform, y = pagerelform)) +
  geom_point(alpha = 0.35, size = 0.5) +
  geom_abline(size = 1,
              aes(intercept = 0, slope = 1, linetype = "Same age"),
              show.legend = FALSE) +
  #facet_grid(. ~ Gender) +
  scale_y_continuous(name = "Partner's ages") +
  # scale_linetype_manual('Lines',
  #                        values = c("Population mean" = 1, "Same age" = 2)) +
  xlab("Individual in population") +
  guides(linetype = guide_legend(keywidth = 2, keyheight = 1)) +
  coord_fixed(ratio = 1) #+ theme

#  Histogram of age differences
ggplot(data = dplyr::filter(agemix.df, Gender == "male", agerelform >=18, agerelform < 50),
       aes(AgeGap, ..density..)) +
  geom_histogram()


# HIV prevalence over time
prevalence.plotter(datalist = datalist.agemix, agegroup = c(15, 50))

datalist.agemix$ttable


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
                                  actual.input.matrix = matrix(rep(c(3, -0.5), 32), nrow = 32),
                                  seed_count = 0,
                                  n_cluster = 8)
summary(sumsmstats.df)
# Min.   :3.419   Min.   :4.265
# 1st Qu.:3.685   1st Qu.:4.387
# Median :3.743   Median :4.482
# Mean   :3.750   Mean   :4.479
# 3rd Qu.:3.875   3rd Qu.:4.528
# Max.   :3.971   Max.   :4.684

miceabc.output <- MiceABC(maxwaves = 1) # This was overwritten by a function call, identical to miceabc.output2

# The miceabc output that I will work with for the SIMID presentation
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

# The rejection abc output that I will work with for the SIMID presentation
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

# The Lenormand abc output that I will work with for the SIMID presentation
# Second comparison: let's run Lenormand's adaptive population ABC scheme
help("ABC_sequential")
lenormandabc.output2 <- ABC_sequential(method = "Lenormand",
                                      model = simpact.wrapper,
                                      prior = list(c("unif", -7, -2),
                                                   c("unif", -7, -2)),
                                      nb_simul = 64,
                                      summary_stat_target = c(3.75, 4.48),   # this and tol is commented out, so that we can see ALL the output, not just the retained fraction.
                                      alpha = 0.25,
                                      p_acc_min = 0.08, # 0.05,
                                      use_seed = TRUE,
                                      seed_count = 0,
                                      n_cluster = 8,
                                      verbose = TRUE)


# On the iMac I saved these outputs, and I will now load them and work from there.

# save(simpact.wrapper, simpact.parallel, MiceABC, sumsmstats.df, miceabc.output, miceabc.output2, rejectabc.output, lenormandabc.output, file = "/Users/wimdelva/Documents/miceabc.comparison.RData")

abc.comparison <- load(file = "/Users/delvaw/Downloads/miceabc.comparison.RData")


# First, let's get all the model output in the same format:
# x1, x2, y1, y2, best16
# 1. For rejectabc
# Calculating squared relative distances and sum thereof for the ABC rejection method
targets <- c(3.75, 4.48)
rejectabc.output$y.1.sq.rel.dist <- ((rejectabc.output$stats[,1] - targets[1]) / targets[1])^2
rejectabc.output$y.2.sq.rel.dist <- ((rejectabc.output$stats[,2] - targets[2]) / targets[2])^2
rejectabc.output$sum.sq.rel.dist <- rejectabc.output$y.1.sq.rel.dist + rejectabc.output$y.2.sq.rel.dist
dist.order <- order(rejectabc.output$sum.sq.rel.dist) # Ordering the squared distances from small to big.
last.one.selected <- 16
selected.distances <- dist.order[1:last.one.selected]
best16.rejection <- rep(FALSE, length(dist.order))
best16.rejection[selected.distances] <- TRUE

abcreject.df <- as.data.frame(cbind(rejectabc.output$param,
                                    rejectabc.output$stats))
abcreject.df$best16 <- best16.rejection
abcreject.df$method <- "rejection"
names(abcreject.df) <- c("x1", "x2", "y1", "y2", "best16", "method")

# 2. For Lenormandabc
# Reading in the model_step1 to model_step5 files with input and output values
Lenormand_inout_all <- vector("list", 5)
for (step in 1:5){
  filename <- paste0("model_step", step)
  Lenormand_inout <- read.delim(file = filename,
                                header = FALSE,
                                sep = " ")
  Lenormand_inout_all <- rbind(Lenormand_inout_all, Lenormand_inout)
}
Lenormand.best16.y1 <- read.delim(file = "output_step5", header = FALSE, sep = " ")#lenormandabc.output2$intermediary[[5]]$posterior[,4]
best16.Lenormand <- round(as.numeric(unlist(Lenormand_inout_all["V4"])), 6) %in% round(as.numeric(unlist(Lenormand.best16.y1["V4"])), 6)
abcLenorm.df <- as.data.frame(cbind(Lenormand_inout_all[, 2:5],
                                    best16.Lenormand))
abcLenorm.df$method <- "Lenormand"
names(abcLenorm.df) <- c("x1", "x2", "y1", "y2", "best16", "method")

# 3. For miceabc
miceabc_inout_all <- do.call(rbind, miceabc.output2$sim.results.with.design)[,1:4]
best16.mice <- miceabc_inout_all$y.1 %in% miceabc.output2$inputoutput[[5]]$y.1
abcmice.df <- as.data.frame(cbind(miceabc_inout_all,
                                  best16.mice))
abcmice.df$method <- "mice"
names(abcmice.df) <- c("x1", "x2", "y1", "y2", "best16", "method")

# Adding the target statistics and their inputs
target.df <- data.frame(x1 = -6, x2 = -6, y1 = 3.75, y2 = 4.48, best16 = TRUE, method = "target")

# Now we can rbind these 4 datasets into one dataset that we can plot with ggplot2
abc.comparison.df <- rbind(abcreject.df,
                           abcLenorm.df,
                           abcmice.df,
                           target.df)
abc.comparison.df$method <- factor(abc.comparison.df$method, labels = c("Lenormand",
                                                                        "MICE",
                                                                        "Rejection",
                                                                        "Target"))

abc.comparison.df$method <- factor(abc.comparison.df$method,levels(abc.comparison.df$method)[c(3, 1, 2, 4)]) # We want the levels in a particular order, for plotting and colour purposes
levels(abc.comparison.df$method)

input.plot <- ggplot(data = abc.comparison.df,
                     aes(x = x1,
                         y = x2,
                         col = method,
                         alpha = best16,
                         shape = method == "Target",
                         size = best16)) +
  scale_colour_manual(values = cb.Palette) +
  scale_alpha_discrete(breaks = NULL,
                       range = c(0.4, 1)) +
  scale_shape_manual(breaks = NULL,
                     values = c(20, 10)) +
  scale_size_manual(breaks = NULL,
                    values = c(2, 4)) +
  xlab("Preferred Age Gap at age 0") +
  ylab("Penalty of deviation from Preferred Age Gap") +
  geom_point() +
  theme_bw() +
  guides(colour = guide_legend(title = "Method",
                               override.aes = list(shape = c(20, 20, 20, 10),
                                                   size = c(4, 4, 4, 4))))
plot(input.plot)


output.plot <- ggplot(data = abc.comparison.df,
                     aes(x = y1,
                         y = y2,
                         col = method,
                         alpha = best16,
                         shape = method == "Target",
                         size = best16)) +
  scale_colour_manual(values = cb.Palette) +
  scale_alpha_discrete(breaks = NULL,
                       range = c(0.4, 1)) +
  scale_shape_manual(breaks = NULL,
                     values = c(20, 10)) +
  scale_size_manual(breaks = NULL,
                    values = c(2, 4)) +
  xlab("Average Age Gap") +
  ylab("Standard Deviation of Age Gaps") +
  geom_point() +
  theme_bw() +
  guides(colour = guide_legend(title = "Method",
                               override.aes = list(shape = c(20, 20, 20, 10),
                                                   size = c(4, 4, 4, 4))))
plot(output.plot)


bestoutput.plot <- ggplot(data = filter(abc.comparison.df, best16 == TRUE),
                      aes(x = y1,
                          y = y2,
                          col = method,
                          shape = method == "Target",
                          size = best16)) +
  scale_colour_manual(values = cb.Palette) +
  scale_shape_manual(breaks = NULL,
                     values = c(20, 10)) +
  scale_size_manual(breaks = NULL,
                    values = 4) +
  xlab("Average Age Gap") +
  ylab("Standard Deviation of Age Gaps") +
  geom_point() +
  theme_bw() +
  guides(colour = guide_legend(title = "Method",
                               override.aes = list(shape = c(20, 20, 20, 10),
                                                   size = c(4, 4, 4, 4))))
plot(bestoutput.plot)



# Plotting all the input parameter sets and then overlaying this with the final calibration results (16 best)

all.results.with.design <- do.call(rbind, miceabc.output2$sim.results.with.design)
dim(all.results.with.design)
plot(all.results.with.design[,1],
     all.results.with.design[,2],
     xlab = "x1",
     ylab = "x2",
     xlim = c(-7, -2),
     ylim = c(-7, -2),
     pch = 20)
points(-6, -6, pch = 3, col = "red3", lwd = 3)

# miceabc, rejection and Lenormand wrt proximity to target statistics.
# 1. rejection abc
plot(abcreject.sim.results.with.design.df$y1[abcreject.sim.results.with.design.df$shortest16 == 1],
     abcreject.sim.results.with.design.df$y2[abcreject.sim.results.with.design.df$shortest16 == 1],
     xlab = "y1",
     ylab = "y2",
     xlim = c(3.5,4.1), #c(min(abcreject.sim.results.with.design.df$y1), max(abcreject.sim.results.with.design.df$y1)),
     ylim = c(4.3,4.65), #c(min(abcreject.sim.results.with.design.df$y2), max(abcreject.sim.results.with.design.df$y2)),
     pch = 20, # - abcreject.sim.results.with.design.df$shortest16,
     col = 6, #abcreject.sim.results.with.design.df$shortest16 * 5 + 1,
     cex = 1)

# 2. Lenormand
points(lenormandabc.output2$intermediary[[5]]$posterior[,4],
       lenormandabc.output2$intermediary[[5]]$posterior[,5],
       xlab = "y1",
       ylab = "y2",
       xlim = c(3.5, 4),
       ylim = c(4.25, 4.75),
       col = "red",
       pch = 20)

# 3. miceabc
points(miceabc.output2$inputoutput[[5]][,3],
       miceabc.output2$inputoutput[[5]][,4],
       xlab = "y1",
       ylab = "y2",
       xlim = c(3.5, 4),
       ylim = c(4.25, 4.75),
       col = "yellow2",
       pch = 20)

points(3.75, 4.48, pch = 3, col = "black", lwd = 5)


plot(miceabc.output$inputoutput[[5]][,3],
     miceabc.output$inputoutput[[5]][,4],
     xlab = "y1",
     ylab = "y2",
     xlim = c(3.5, 4.1),
     ylim = c(4.25, 4.75),
     pch = 20)
points(3.75, 4.48, pch = 3, col = "red3", lwd = 3)

points(miceabc.output2$inputoutput[[5]][,3],
     miceabc.output2$inputoutput[[5]][,4],
     xlab = "y1",
     ylab = "y2",
     xlim = c(3.5, 4),
     ylim = c(4.25, 4.75),
     col = "yellow2",
     pch = 20)



# Calculating squared relative distances and sum thereof for the ABC rejection method
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






















### From here, it's old code

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

