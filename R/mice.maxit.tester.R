# mice maxit tester

load(file = "/Users/delvaw/Documents/MiceABC/test.VEME2.MaC.incremental.RData")
vsc.data.loaded <- load(file = "/Users/wimdelva/Documents/MiceABC/test.600.4waves.VEME2.MaC.increm.vsc.RData")
dummy.targets.empirical
test.VEME2.MaC.incremental <- test.VEME2.MaC.increm.vsc

test.VEME2.MaC.increm.vsc$max.RMSD
test.VEME2.MaC.increm.vsc$n.close.to.targets

head(test.VEME2.MaC.increm.vsc$selected.experiments[[4]])
# holds dummy.targets.empirical and test.VEME2.MaC.incremental
library(mice)
library(dplyr)
library(tidyr)
# help(mice)

length(test.VEME2.MaC.increm.vsc$selected.experiments)

summary(test.VEME2.MaC.incremental$selected.experiments[[1]][, 1:15])
hist(test.VEME2.MaC.incremental$selected.experiments[[1]][, 1])

pairs(test.VEME2.MaC.incremental$selected.experiments[[20]][, 1:15])


hist(mice.guesses[, 1])

off.target.factor <- 1

data4mice <- test.VEME2.MaC.incremental$selected.experiments[[1]] %>% dplyr::select(., -RMSD) %>% rbind(., c(rep(NA, 15), off.target.factor * dummy.targets.empirical))

str(data4mice)
dummy.input.vector <- c(1.1, 0.25, 0, 3, 0.23, 0.23,  # what if 1.1 becomes 1.4
                        45, 45, #45, 45,              # what if 45 becomes 60
                        -0.5, 2.8, -0.2, -0.2, -2.5, -0.52, -0.05)# c(1000, 2, 3, 4)
x.offset <- length(dummy.input.vector)

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

# The distribution of values for Parameters 1 to 15 in the training dataset:
percentiles.data4mice <- data.frame(t(apply(head(dplyr::select(data4mice,
                                                                       x.1:x.15), -1), 2, quantile, c(0.05, 0.5, 0.95))),
                                    row.names = NULL)

percentiles.data4mice$range5.95 <- percentiles.data4mice$X95. - percentiles.data4mice$X5.
percentiles.data4mice$Var <- 1:15
percentiles.data4mice$Type <- "Training"

# Now the range and medians for the mice imputations
mice.guesses <- mice.parallel(mice.model = mice.wrapper,
                                df.give.to.mice = data4mice,
                                m = 1,
                                predictorMatrix = predictorMatrix,
                                method = "norm",
                                defaultMethod = "norm",
                                maxit = 5,
                                printFlag = FALSE,
                                seed_count = 0,
                                n_cluster = 8,
                                nb_simul = 1000)
mice.guesses.df <- data.frame(mice.guesses)


percentiles.mice.guesses <- data.frame(t(apply(mice.guesses, 2, quantile, c(0.05, 0.5, 0.95))),
                                    row.names = NULL)

percentiles.mice.guesses$range5.95 <- percentiles.mice.guesses$X95. - percentiles.mice.guesses$X5.
percentiles.mice.guesses$Var <- 1:15
percentiles.mice.guesses$Type <- "Imputations"


TrainingAndImputations <- rbind(percentiles.data4mice, percentiles.mice.guesses)

###
# Plotting the result
###


compar.1 <- ggplot(data = TrainingAndImputations,
                     aes(Var, range5.95, colour = Type)) +
  geom_point() +
  xlab("Parameter") +
  ylab("Distance between p5 and p95")
plot(compar.1)

compar.2 <- ggplot(data = TrainingAndImputations,
                   aes(Var, X50., colour = Type)) +
  geom_point() +
  xlab("Parameter") +
  ylab("Median")
plot(compar.2)












output <- data.frame(p5 = double(),
                     p25 = double(),
                     p50 = double(),
                     p75 = double(),
                     p95 = double(),
                     m.i = double())
m.i.vector <- c(1:5, 10, 15, 20, 25)
for (m.i in 1:length(m.i.vector)){
  print(m.i.vector[m.i])
  mice.guesses <- mice.parallel(mice.model = mice.wrapper,
                        df.give.to.mice = data4mice,
                        m = 1,
                        predictorMatrix = predictorMatrix,
                        method = "norm",
                        defaultMethod = "norm",
                        maxit = m.i.vector[m.i],
                        printFlag = FALSE,
                        seed_count = 0,
                        n_cluster = 8,
                        nb_simul = 1000)

  # Summary of output:
  maxit <- m.i.vector[m.i]
  percentiles <- cbind(t(apply(mice.guesses, 2, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95))), maxit)
  output <- rbind(output, percentiles)
}

output$Var <- as.factor(paste0("Parameter ", rep(1:15, times = length(m.i.vector))))

output.long <- gather(data = output,
                      key = Percentile,
                      value = value,
                      1:5)


###
# Plotting the result
###
library(metafolio)
n.colours <- 5
cols <- gg_color_hue(n.colours, l=65)
darkcols <- gg_color_hue(n.colours, l=40)

# Making sure that the order of the variables is ok
output.long <- within(output.long, Var <- factor(Var, levels = paste0("Parameter ", 1:15)))
output.long <- within(output.long, Percentile <- factor(Percentile, levels = c("95%",
                                                                               "75%",
                                                                               "50%",
                                                                               "25%",
                                                                               "5%")))

with(output.long, levels(Var))


maxit.plot <- ggplot(data = output.long,
                     aes(maxit, value, group = Percentile, colour = Percentile)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ Var, nrow = 3, scales="free") +

  xlab("Iterations") +
  ylab("Percentile value") +
  #  scale_colour_hue(l = c(rep(65, 5), rep(10, 5))) +
  scale_color_manual(values = c("5%"=cols[1],
                                "25%"=cols[2],
                                "50%"=cols[3],
                                "75%"=cols[4],
                                "95%"=cols[5]))
plot(maxit.plot)

# l1median example
library(pcaPP)
x.test <- runif(1000, min = -10, max = 10)
y.test <- x.test^2 + rnorm(1000, mean = 0, sd = 10)
test.mat <- matrix(cbind(x.test, y.test), ncol = 2)
l1med <- l1median(test.mat)
l1med
univ.med <- apply(test.mat, 2, median)

plot(x.test, y.test, pch = 16, cex = 0.5)
points(l1med[1], l1med[2], col = "green3", pch = 16)
points(univ.med[1], univ.med[2], col = "orange3", pch = 16)


library(mvtnorm)
x <- rbind(rmvnorm(200, rep(0, 4), diag(c(1, 1, 2, 2))),
           rmvnorm( 50, rep(3, 4), diag(rep(2, 4))))
l1median(x, trace = -1)
# compare with coordinate-wise median:
apply(x,2,median)
