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

#source("/Users/delvaw/Documents/MiceABC/R/VEME.wrapper.R")
source("/Users/delvaw/Documents/MiceABC/R/VEME.wrapper2.R")
source("/Users/delvaw/Documents/MiceABC/R/mice.wrapper.R")
source("/Users/delvaw/Documents/MiceABC/R/00-Functions.R")

source("/Users/delvaw/Documents/MiceABC/R/simpact.parallel.R")
source("/Users/delvaw/Documents/MiceABC/R/mice.parallel.R")

source("/Users/delvaw/Documents/MiceABC/R/dummy.MaC.incremental.R")
source("/Users/delvaw/Documents/MiceABC/R/dummy.MaC.incremental.parallel.mice.R")

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
library(parallel)


dummy.input.vector <- c(1.1, 0.25, 0, 3, 0.23, 0.23,  # what if 1.1 becomes 1.4
                        45, 45, #45, 45,              # what if 45 becomes 60
                        -0.5, 2.8, -0.2, -0.2, -2.5, -0.52, -0.05)# c(1000, 2, 3, 4)
x.offset <- length(dummy.input.vector)
n.experiments <- 80
dummy.master2 <- simpact.parallel(model = VEME.wrapper2,
                             actual.input.matrix = matrix(rep(dummy.input.vector, each = n.experiments), nrow = n.experiments),
                             seed_count = 0,
                             n_cluster = 8)

save(valid.dummy.master2, file = "/Users/delvaw/Documents/MiceABC/valid.dummy.master2.RData")






#####
# The output of the master model
#####
#inc.master.vector only exists for the validation version of the master model

# inc.master.vector <- dummy.master2[, 17]
# #save(inc.master.vector, file = "/Users/delvaw/Documents/MiceABC/inc.master.vector.RData")
# #save(inc.master.vector.inflated, file = "/Users/delvaw/Documents/MiceABC/inc.master.vector.inflated.RData")
#
# hist(inc.master.vector[!is.na(inc.master.vector)], 14) # The distribution of HIV incidence after 10 years
# mean(inc.master.vector) # The mean HIV incidence after 10 years: 0.01527118
# median(inc.master.vector) # The median HIV incidence after 10 years: 0.01532902
# quantile(inc.master.vector, c(0.025, 0.975)) # 0.002462209 0.026417416
head(dummy.master2)
inc.master.vector <- dummy.master2[, 17]
new.infect.vector <- dummy.master2[, 18]
recent.ratio.vector <- dummy.master2[, 19]

mean(new.infect.vector, na.rm = TRUE) # The mean HIV incidence after 10 years
median(new.infect.vector, na.rm = TRUE) # The median HIV incidence after 10 years
quantile(new.infect.vector, c(0.025, 0.975), na.rm = TRUE)


# mean.br.len.vector <- dummy.master2[, 19]


hist(new.infect.vector[!is.na(new.infect.vector) & new.infect.vector < Inf], 20)
hist(recent.ratio.vector[!is.na(recent.ratio.vector) & recent.ratio.vector < Inf], 20)

plot(new.infect.vector[!is.na(new.infect.vector) & new.infect.vector < Inf],
     inc.master.vector[!is.na(new.infect.vector) & new.infect.vector < Inf])

plot(recent.ratio.vector[inc.master.vector > 0.0142 & inc.master.vector < 0.0162],
     inc.master.vector[inc.master.vector > 0.0142 & inc.master.vector < 0.0162])

cor.test(inc.master.vector[!is.na(new.infect.vector) & new.infect.vector < Inf],
         new.infect.vector[!is.na(new.infect.vector) & new.infect.vector < Inf])

cor.test(inc.master.vector[!is.na(new.infect.vector) & new.infect.vector < Inf],
         new.infect.vector[!is.na(new.infect.vector) & new.infect.vector < Inf])

#plot(mean.br.len.vector[!is.na(mean.br.len.vector) & mean.br.len.vector < Inf],
#     inc.master.vector[!is.na(mean.br.len.vector) & mean.br.len.vector < Inf])

#cor.test(mean.br.len.vector[!is.na(mean.br.len.vector) & mean.br.len.vector < Inf],
#     inc.master.vector[!is.na(mean.br.len.vector) & mean.br.len.vector < Inf])

#plot(mean.br.len.vector,
#     new.infect.vector)


dummy.master2 <- dummy.master2 %>%
  as.data.frame() %>%
  dplyr::filter(complete.cases(.))


dummy.targets.empirical <- l1median(dummy.master2)
# dummy.targets.empirical (based on 400 repeats):
# 5.79782475  4.16132170  0.63770526  2.69246630  1.21359926 -1.56183331  0.25472579  2.09719359  0.02891330  0.13671790
# 0.06816378  0.34554390  0.31086783  0.37050293  0.41706055  1.01682475  7.07415631

dummy.targets.empirical <- c(5.79782475,
                             4.16132170,
                             0.63770526,
                             2.69246630,
                             1.21359926,
                             -1.56183331,
                             0.25472579,
                             2.09719359,
                             0.02891330,
                             0.13671790,
                             0.06816378,
                             0.34554390,
                             0.31086783,
                             0.37050293,
                             0.41706055,
                             1.01682475,
                             7.07415631)
# round(colMeans(dummy.master2), 3) # For interest sake, what are the marginal means?


predictorMatrix <- (1 - diag(1, length(c(dummy.input.vector, dummy.targets.empirical)))) # This is the default matrix.
# # Let's now modify the first 15 rows of this matrix, corresponding to the indicators of predictor variables for the input variables. In brackets the values for the master model.

predictorMatrix[1:x.offset, ] <- 0 # First we "empty" the relevant rows, then we refill them.
# We are currently not allowing input variables to be predicted by other predictor variables. Only via output variables. We could change this at a later stage.

predictorMatrix[1, x.offset + c(10, 11, 17)] <- 1 # relative susceptibility in young women is predicted by HIV prevalence in young men and women, and recent infections (~ incidence)
      predictorMatrix[2, x.offset + 3] <- 1 # agescale predicted by slope
      predictorMatrix[3, x.offset + c(1, 3, 6)] <- 1 # mean of the person-specific age gap preferences is predicted by slope, intercept and AAD
      predictorMatrix[4, x.offset + c(2, 4, 5)] <- 1 # sd of the person-specific age gap preferences is predicted by SD, WSD, BSD
      predictorMatrix[5, x.offset + c(7, 8, 9, 13, 16)] <- 1 # man gamma a predicted by gamma shape.male, scale.male, pp.cp, hiv.prev.25.34.men, exp(growthrate)
      predictorMatrix[6, x.offset + c(7, 8, 9, 12, 16)] <- 1 # woman gamma a predicted by gamma shape.male, scale.male, pp.cp, hiv.prev.25.34.women, exp(growthrate)
      predictorMatrix[7, x.offset + c(7, 8, 9, 13, 16, 17)] <- 1 # man gamma b predicted by gamma shape.male, scale.male, pp.cp, hiv.prev.25.34.men, exp(growthrate), and recent infections (~ incidence)
      predictorMatrix[8, x.offset + c(7, 8, 9, 12, 16, 17)] <- 1 # woman gamma b predicted by gamma shape.male, scale.male, pp.cp, hiv.prev.25.34.men, exp(growthrate), and recent infections (~ incidence)
      predictorMatrix[9, x.offset + c(2, 4, 5, 7, 8, 14, 15, 16, 17)] <- 1 # formation.hazard.agegapry.gap_factor_x_exp is predicted by population growth, age gap variance, hiv prevalence, and recent infections (~ incidence)
      predictorMatrix[10, x.offset + c(7, 8, 9, 12, 13, 16, 17)] <- 1 # baseline formation hazard predicted by HIV prevalence, cp, degree distrib. HIV prevalence, and recent infections (~ incidence)
      predictorMatrix[11, x.offset + c(7, 8, 9, 12, 13, 16, 17)] <- 1 # numrel man penalty is predicted by degree distrib, cp, prev, popgrowth, and recent infections (~ incidence)
      predictorMatrix[12, x.offset + c(7, 8, 9, 12, 13, 16, 17)] <- 1 # # numrel woman penalty is predicted by degree distrib, cp, prev, popgrowth, and recent infections (~ incidence)
      predictorMatrix[13, x.offset + 16] <- 1 # conception.alpha_base is predicted by popgrowth
      predictorMatrix[14, x.offset + c(7, 8, 9, 16)] <- 1 # baseline dissolution hazard predicted by degree distrib, cp, popgrowth
      predictorMatrix[15, x.offset + c(7, 8, 9, 16)] <- 1 # age effect on dissolution hazard predicted by degree distrib, cp, popgrowth, HIV prev in older people (maybe?)

      # NOTE: As it stands, each output statistic is predicted by ALL input and ALL other output statistics. That may not be a great idea, or even possible, if there is collinearity.

# Test dummy.MaC.incremental, and also dummy.MaC.incremental.parallel.mice

test.VEME2.MaC.incremental <- dummy.MaC.incremental.parallel.mice(targets.empirical = dummy.targets.empirical,
                                        RMSD.tol.max = 0.95,
                                        min.givetomice = 20, # 400
                                        n.experiments = 80, # 1000
                                        lls = c(1,  0.12, -0.3, 2.5, 0.1, 0.1, 20, 20, -0.8, 2, -0.35, -0.35, -3.6, -0.8, -0.16),
                                        uls = c(1.2, 0.37, 0.3, 3.5, 0.4, 0.4, 66, 66, -0.25, 3.9, -0.1, -0.1, -1.4, -0.3,  -0.001),
                                        model = VEME.wrapper2,
                                        strict.positive.params = c(4:8),
                                        predictorMatrix = predictorMatrix,
                                        maxit = 5,
                                        maxwaves = 6,
                                        n_cluster = 8) # 6
#(round(l1median(head(test.MaC.incremental$selected.experiments[[length(test.MaC.incremental$selected.experiments)]]), 1), 99)[5:8] - dummy.targets.empirical[1:4]) / dummy.targets.empirical[1:4]
#round(l1median(head(test.MaC.incremental$selected.experiments[[length(test.MaC.incremental$selected.experiments)]]), 1), 2)
#test.MaC.incremental$secondspassed

test.VEME2.MaC.incremental$secondspassed
test.VEME2.MaC.incremental$max.RMSD
test.VEME2.MaC.incremental$n.close.to.targets
head(test.VEME2.MaC.incremental$selected.experiments[[length(test.VEME2.MaC.incremental$selected.experiments)]])

save(dummy.targets.empirical, test.VEME2.MaC.incremental, file = "/Users/delvaw/Documents/MiceABC/test.VEME2.MaC.incremental.RData")


### Now we simulate for the 1 (or 5?) best fitting model(s)
inputs.calib <- as.numeric(test.VEME2.MaC.incremental$selected.experiments[[length(test.VEME2.MaC.incremental$selected.experiments)]][1, 1:15])


calib.dummy.master2 <- simpact.parallel(model = VEME.wrapper2,
                                        actual.input.matrix = matrix(rep(inputs.calib, each = n.experiments), nrow = n.experiments),
                                        seed_count = 0,
                                        n_cluster = 8)

save(calib.dummy.master2, file = "/Users/delvaw/Documents/MiceABC/calib.dummy.master2.RData")



## AND THE OUTPUT LOOKS LIKE:
valid.inc.wide <- as.data.frame(valid.dummy.master2[, 17:31])
names(valid.inc.wide) <- 1:15 #paste0("incid.", 1:15)
valid.inc.wide$rep <- 1:nrow(valid.inc.wide)
valid.inc.wide$type <- "Master model"

calib.inc.wide <- as.data.frame(calib.dummy.master2[, 17:31])
names(calib.inc.wide) <- 1:15 #paste0("incid.", 1:15)
calib.inc.wide$rep <- (1 + nrow(valid.inc.wide)):(nrow(valid.inc.wide) + nrow(calib.inc.wide))
calib.inc.wide$type <- "Calibrated model"

inc.both <- rbind(valid.inc.wide, calib.inc.wide)
inc.both.long <- gather(inc.both, year, incidence, 1:15)
inc.both.long$year <- as.numeric(inc.both.long$year)
inc.both.long$smoother <- factor(paste0("smoother.", inc.both.long$type))


###
# Plotting the result
###
library(metafolio)
n.colours <- 2
cols <- gg_color_hue(n.colours, l=65)
darkcols <- gg_color_hue(n.colours, l=40)

incplot <- ggplot(filter(inc.both.long,
                         year >=10),
                  aes(year, incidence, group = rep, colour = type)) +
  geom_line() +
  facet_wrap(~ type) +
  stat_smooth(se=FALSE,
              method="loess",
              span=1,
              aes(year, incidence, group = type, colour=smoother)) +
  xlab("Time (years)") +
  ylab("HIV incidence (cases/person-year)") +
  #  scale_colour_hue(l = c(rep(65, 5), rep(10, 5))) +
  scale_color_manual(values = c("Master model"=cols[1],
                                "Calibrated model"=cols[2],
                                "smoother.Master model"=darkcols[1],
                                "smoother.Calibrated model"=darkcols[2]),
                     guide = FALSE)
plot(incplot)


