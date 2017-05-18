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
  library(magrittr)
  library(dplyr)
  library(nlme)

  agemixing.lme.errFunction <- function(e)
  {
    return(list())
  }

  agemixing.lme.fitter <- function(data = dplyr::filter(agemix.model[[1]], Gender =="male")){
    men.lme <- lme(pagerelform ~ agerelform0,
                   data = data,
                   control=lmeControl(maxIter=10000, returnObject=TRUE, opt = "optim"),
                   random = ~1 | ID,
                   method = "REML",
                   weight = varPower(value = 0.5, form = ~agerelform0 + 1))
  }




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

  men.lme <- tryCatch(agemixing.lme.fitter(data = dplyr::filter(agemix.model[[1]], Gender =="male")),
                      error = agemixing.lme.errFunction) # Returns an empty list if the lme model can't be fitted

  bignumber <- 999

  AAD <- ifelse(length(men.lme) > 0, mean(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
  SDAD <- ifelse(length(men.lme) > 0, sd(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
  powerm <- ifelse(length(men.lme) > 0, as.numeric(attributes(men.lme$apVar)$Pars["varStruct.power"]), bignumber)
  slopem <- ifelse(length(men.lme) > 0, summary(men.lme)$tTable[2, 1], bignumber)
  WVAD.base <- ifelse(length(men.lme) > 0, men.lme$sigma^2, bignumber)
  BVAD <- ifelse(length(men.lme) > 0, getVarCov(men.lme)[1,1], bignumber)

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
  scale_y_continuous(name = "Female partner age") +
  # scale_linetype_manual('Lines',
  #                        values = c("Population mean" = 1, "Same age" = 2)) +
  xlab("Male partner age") +
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



sumsmstatsvect <- simpact.wrapper(inputvector = c(1, 3, -0.5))
sumsmstatsvect # Mean age gap of 4.3 years and sd of 4.6 years

sumsmstats.df <- simpact.parallel(model = simpact.wrapper,
                                  actual.input.matrix = matrix(rep(c(3, -0.5), 32), nrow = 32, byrow = TRUE),
                                  seed_count = 0,
                                  n_cluster = 8)
summary(sumsmstats.df)
# V1              V2              V3                 V4               V5              V6               V7
# Min.   :5.480   Min.   :3.653   Min.   :-0.08131   Min.   :0.6135   Min.   :3.051   Min.   :0.8859   Min.   :0.01293
# 1st Qu.:5.741   1st Qu.:3.847   1st Qu.: 0.02042   1st Qu.:0.6549   1st Qu.:4.082   1st Qu.:1.2559   1st Qu.:0.11263
# Median :5.844   Median :3.925   Median : 0.06849   Median :0.6623   Median :4.968   Median :1.4600   Median :0.13981
# Mean   :5.816   Mean   :3.924   Mean   : 0.05916   Mean   :0.6628   Mean   :5.229   Mean   :1.4985   Mean   :0.12573
# 3rd Qu.:5.884   3rd Qu.:4.008   3rd Qu.: 0.10639   3rd Qu.:0.6758   3rd Qu.:5.973   3rd Qu.:1.6189   3rd Qu.:0.15104
# Max.   :6.267   Max.   :4.250   Max.   : 0.17067   Max.   :0.6904   Max.   :9.306   Max.   :2.4328   Max.   :0.17041
# V8
# Min.   :0.001797
# 1st Qu.:0.002722
# Median :0.003125
# Mean   :0.003311
# 3rd Qu.:0.003871
# Max.   :0.005951

miceabc.output.8V <- MiceABC.V8(maxwaves = 1) # This was overwritten by a function call, identical to miceabc.output2

# The miceabc output that I will work with for the SIMID presentation
ptm <- proc.time()
miceabc.output2.8V <- MiceABC.V8(targets = c(5.816, 3.924, 0.05916, 0.6628, 5.229, 1.4985, 0.12573, 0.003311),
                           n.experiments = 64,
                           lls = c(0.1, -3),
                           uls = c(5, -0.1),
                           model = simpact.wrapper,
                           maxit = 20,
                           maxwaves = 5, # and we can also try 10 to do the same as the 10 waves of Lenormand
                           #reps = 5,
                           alpha = 0.25,
                           saturation.crit = 0)
mice.timeused <- proc.time() - ptm
mice.timeused # 1745 seconds

ptm <- proc.time()
miceabc.output2.8V.10waves <- MiceABC.V8(targets = c(5.816, 3.924, 0.05916, 0.6628, 5.229, 1.4985, 0.12573, 0.003311),
                                 n.experiments = 64,
                                 lls = c(0.1, -3),
                                 uls = c(5, -0.1),
                                 model = simpact.wrapper,
                                 maxit = 20,
                                 maxwaves = 10, # and we can also try 10 to do the same as the 10 waves of Lenormand
                                 #reps = 5,
                                 alpha = 0.25,
                                 saturation.crit = 0)
mice.timeused <- proc.time() - ptm
mice.timeused # 2974.431 seconds

# The rejection abc output that I will work with for the SIMID presentation
# As comparison, let's run a basic ABC accept rejection scheme
help("ABC_rejection")

ptm <- proc.time()
rejectabc.output <- ABC_rejection(model = simpact.wrapper,
                                  prior = list(c("unif", 0.1, 5),
                                               c("unif", -3, -0.1)),
                                  nb_simul = 254, # and we can also try 496 to do the same as the 10 waves of Lenormand
                                  #summary_stat_target = c(3.75, 4.48),   # this and tol is commented out, so that we can see ALL the output, not just the retained fraction.
                                  #tol = 0.064, # That is 16/250 to be in line with the output of miceabc
                                  use_seed = TRUE,
                                  seed_count = 0,
                                  n_cluster = 8)
rejection.timeused <- proc.time() - ptm
rejection.timeused # 1248.939 seconds


ptm <- proc.time()
rejectabc.output.10waves <- ABC_rejection(model = simpact.wrapper,
                                  prior = list(c("unif", 0.1, 5),
                                               c("unif", -3, -0.1)),
                                  nb_simul = 496, # and we can also try 496 to do the same as the 10 waves of Lenormand
                                  #summary_stat_target = c(3.75, 4.48),   # this and tol is commented out, so that we can see ALL the output, not just the retained fraction.
                                  #tol = 0.064, # That is 16/250 to be in line with the output of miceabc
                                  use_seed = TRUE,
                                  seed_count = 0,
                                  n_cluster = 8)
rejection.timeused.10waves <- proc.time() - ptm
rejection.timeused.10waves # 2244.733 seconds


# The Lenormand abc output that I will work with for the SIMID presentation
# Second comparison: let's run Lenormand's adaptive population ABC scheme
help("ABC_sequential")

ptm <- proc.time()
lenormandabc.output.8V <- ABC_sequential(method = "Lenormand",
                                      model = simpact.wrapper,
                                      prior = list(c("unif", 0.1, 5),
                                                   c("unif", -3, -0.1)),
                                      nb_simul = 64,
                                      summary_stat_target = c(5.816, 3.924, 0.05916, 0.6628, 5.229, 1.4985, 0.12573, 0.003311),
                                      alpha = 0.25,
                                      p_acc_min = 0.08, # 0.05,
                                      use_seed = TRUE,
                                      seed_count = 0,
                                      n_cluster = 8,
                                      verbose = TRUE)
lenormand.timeused <- proc.time() - ptm
lenormand.timeused # 5608.776 seconds


# On the iMac I saved these outputs, and I will now load them and work from there.

save(simpact.wrapper, simpact.parallel, MiceABC, MiceABC.V8, sumsmstats.df, miceabc.output.8V, miceabc.output2.8V, miceabc.output2.8V.10waves, rejectabc.output, rejectabc.output.10waves, lenormandabc.output.8V, file = "/Users/delvaw/Documents/miceabc.V8.comparison.RData")

# save(simpact.wrapper, simpact.parallel, MiceABC, sumsmstats.df, miceabc.output, miceabc.output2, rejectabc.output, lenormandabc.output, file = "/Users/wimdelva/Documents/miceabc.comparison.RData")

abc.comparison.V8 <- load(file = "/Users/delvaw/Documents/miceabc.V8.comparison.RData")


# First, let's get all the model output in the same format:
# x1, x2, y1, y2, ..., y8, best16
# 1. For rejectabc
# Calculating squared relative distances and sum thereof for the ABC rejection method
targets <- c(5.816, 3.924, 0.05916, 0.6628, 5.229, 1.4985, 0.12573, 0.003311)
sum.sq.rel.dist <- rep(0, nrow(rejectabc.output$stats))
for(i in 1:length(targets)) {
  nam <- paste0("y.", i, ".sq.rel.dist")
  val <- ((rejectabc.output$stats[,i] - targets[i]) / targets[i])^2
  assign(nam, val)
  sum.sq.rel.dist <- sum.sq.rel.dist + get(nam)
}

#rejectabc.output$y.1.sq.rel.dist <- ((rejectabc.output$stats[,1] - targets[1]) / targets[1])^2
#rejectabc.output$y.2.sq.rel.dist <- ((rejectabc.output$stats[,2] - targets[2]) / targets[2])^2
#rejectabc.output$sum.sq.rel.dist <- rejectabc.output$y.1.sq.rel.dist + rejectabc.output$y.2.sq.rel.dist
dist.order <- order(sum.sq.rel.dist) # Ordering the squared distances from small to big.
last.one.selected <- 16
selected.distances <- dist.order[1:last.one.selected]
best16.rejection <- rep(FALSE, length(dist.order))
best16.rejection[selected.distances] <- TRUE

abcreject.df <- as.data.frame(cbind(rejectabc.output$param,
                                    rejectabc.output$stats))
abcreject.df$best16 <- best16.rejection
abcreject.df$method <- "rejection"

xvars.names <- paste0("x", 1:2) # This will also be used by the other methods
yvars.names <- paste0("y", 1:8) # This will also be used by the other methods
names(abcreject.df) <- c(xvars.names, yvars.names, "best16", "method")

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
abcLenorm.df <- as.data.frame(cbind(Lenormand_inout_all[, 2:ncol(Lenormand_inout_all)],
                                    best16.Lenormand))
abcLenorm.df$method <- "Lenormand"
names(abcLenorm.df) <- c(xvars.names, yvars.names, "best16", "method")

# And now for the full 10 waves
Lenormand_inout_all.10waves <- vector("list", 10)
for (step in 1:10){
  filename <- paste0("model_step", step)
  Lenormand_inout <- read.delim(file = filename,
                                header = FALSE,
                                sep = " ")
  Lenormand_inout_all.10waves <- rbind(Lenormand_inout_all.10waves, Lenormand_inout)
}
Lenormand.best16.y1 <- read.delim(file = "output_step10", header = FALSE, sep = " ")#lenormandabc.output2$intermediary[[5]]$posterior[,4]
best16.Lenormand <- round(as.numeric(unlist(Lenormand_inout_all.10waves["V4"])), 6) %in% round(as.numeric(unlist(Lenormand.best16.y1["V4"])), 6)
abcLenorm.df <- as.data.frame(cbind(Lenormand_inout_all.10waves[, 2:ncol(Lenormand_inout_all.10waves)],
                                    best16.Lenormand))
abcLenorm.df$method <- "Lenormand"
names(abcLenorm.df) <- c(xvars.names, yvars.names, "best16", "method")


# 3. For miceabc
# First for the 5 waves
xplusyvars.length <- length(xvars.names) + length(yvars.names)
miceabc_inout_all <- do.call(rbind, miceabc.output.8V$sim.results.with.design)[,1:xplusyvars.length]
best16.mice <- miceabc_inout_all$y.1 %in% miceabc.output.8V$inputoutput[[5]]$y.1
abcmice.df <- as.data.frame(cbind(miceabc_inout_all,
                                  best16.mice))
abcmice.df$method <- "mice"
names(abcmice.df) <- c(xvars.names, yvars.names,"best16", "method")

# Then for the 10 waves
xplusyvars.length <- length(xvars.names) + length(yvars.names)
miceabc_inout_all <- do.call(rbind, miceabc.output2.8V.10waves$sim.results.with.design)[,1:xplusyvars.length]
best16.mice <- miceabc_inout_all$y.1 %in% miceabc.output2.8V.10waves$inputoutput[[10]]$y.1
abcmice.df <- as.data.frame(cbind(miceabc_inout_all,
                                  best16.mice))
abcmice.df$method <- "mice"
names(abcmice.df) <- c(xvars.names, yvars.names,"best16", "method")


# Adding the target statistics and their inputs
target.df <- data.frame(matrix(c(3, -0.5, 5.816, 3.924, 0.05916, 0.6628, 5.229, 1.4985, 0.12573, 0.003311), nrow = 1), # The inputs and targets
                        best16 = TRUE, method = "target")
names(target.df) <- c(xvars.names, yvars.names, "best16", "method")


# Now we can rbind these 4 datasets into one dataset that we can plot with ggplot2
abc.comparison.V8.df <- rbind(abcreject.df,
                           abcLenorm.df,
                           abcmice.df,
                           target.df)
abc.comparison.V8.df$method <- factor(abc.comparison.V8.df$method, labels = c("Lenormand",
                                                                        "MICE",
                                                                        "Rejection",
                                                                        "Target"))

abc.comparison.V8.df$method <- factor(abc.comparison.V8.df$method,levels(abc.comparison.V8.df$method)[c(3, 1, 2, 4)]) # We want the levels in a particular order, for plotting and colour purposes
levels(abc.comparison.V8.df$method)

input.plot <- ggplot(data = abc.comparison.V8.df,
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
  xlab("Standard deviation of preferred Age Gap distribution") +
  ylab("Penalty of deviation from Preferred Age Gap") +
  geom_point() +
  theme_bw() +
  guides(colour = guide_legend(title = "Method",
                               override.aes = list(shape = c(20, 20, 20, 10),
                                                   size = c(4, 4, 4, 4))))
plot(input.plot)

input.best.plot <- ggplot(data = filter(abc.comparison.V8.df,
                                        best16 == TRUE),
                     aes(x = x1,
                         y = x2,
                         col = method,
                        # alpha = best16,
                         shape = method == "Target",
                         size = best16)) +
  scale_colour_manual(values = cb.Palette) +
  #scale_alpha_discrete(breaks = NULL,
  #                     range = c(0.4, 1)) +
  scale_shape_manual(breaks = NULL,
                     values = c(20, 10)) +
  scale_size_manual(breaks = NULL,
                    values = 4) +
  xlab("Standard deviation of preferred Age Gap distribution") +
  ylab("Penalty of deviation from Preferred Age Gap") +
  geom_point() +
  theme_bw() +
  guides(colour = guide_legend(title = "Method",
                               override.aes = list(shape = c(20, 20, 20, 10),
                                                   size = c(4, 4, 4, 4))))
plot(input.best.plot)


####################
### OUTPUT PLOTS
####################
# As a reminder:   outputvector <- c(AAD, SDAD, powerm, slopem, WVAD.base, BVAD, hivprev.15.50, growthrate)


# First we remove all rows that contain a value 999 (Lenormand method has this)
abc.comparison.V8.df[abc.comparison.V8.df==999] <- NA
abc.comparison.V8.clean.df <- na.omit(abc.comparison.V8.df)

output.plot <- ggplot(data = abc.comparison.V8.clean.df,
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


bestoutput.plot <- ggplot(data = filter(abc.comparison.V8.clean.df, best16 == TRUE),
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


bestoutput.plot <- ggplot(data = filter(abc.comparison.V8.clean.df, best16 == TRUE),
                          aes(x = y5,
                              y = y6,
                              col = method,
                              shape = method == "Target",
                              size = best16)) +
  scale_colour_manual(values = cb.Palette) +
  scale_shape_manual(breaks = NULL,
                     values = c(20, 10)) +
  scale_size_manual(breaks = NULL,
                    values = 4) +
  xlab("Within-person variance of Age Gaps") +
  ylab("Between-person variance of Age Gaps") +
  geom_point() +
  theme_bw() +
  guides(colour = guide_legend(title = "Method",
                               override.aes = list(shape = c(20, 20, 20, 10),
                                                   size = c(4, 4, 4, 4))))
plot(bestoutput.plot)

# Plots are saved like output.V8.10waves.y1y2.pdf
bestoutput.plot <- ggplot(data = filter(abc.comparison.V8.clean.df, best16 == TRUE),
                          aes(x = y7,
                              y = y8,
                              col = method,
                              shape = method == "Target",
                              size = best16)) +
  scale_colour_manual(values = cb.Palette) +
  scale_shape_manual(breaks = NULL,
                     values = c(20, 10)) +
  scale_size_manual(breaks = NULL,
                    values = 4) +
  xlab("HIV prevalence among 15-50 year olds") +
  ylab("Population growth rate") +
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


time.df <- data.frame(Seconds = c(rep("Rejection", 2245), rep("Lenormand", 5608), rep("MICE", 2974)))

g <- ggplot(time.df, aes(Seconds))
g + geom_bar() + xlab("ABC method") + ylab("Computation time (seconds)") + theme_bw()


df <- data.frame(x = rep(c(2.9, 3.1, 4.5), c(5, 10, 4)))
ggplot(df, aes(x)) + geom_bar()

output.plot <- ggplot(data = abc.comparison.V8.clean.df,
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


