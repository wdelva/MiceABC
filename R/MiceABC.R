#' Adaptive Population ABC with MICE-based supervised learning
#'
#' This function implements an adaptive population ABC scheme with supervised learning through Multivariate Imputation by Chained Equations (MICE) to produce new particles (i.e. samples from the parameter space and their associated output statistics)
#'
#' @param targets The vector of target statistics to calibrate the simulation model against
#' @param n.experiments The number of initial multivariate samples from the parameter space
#' @param lls Lower limits of the parameter ranges to be explored
#' @param uls Upper limits of the parameter ranges to be explored
#' @param probs Indicate in the lls and uls vectors, which elements represent probabilities strictly between 0 and 1
#' @param model Handle to the wrapper function that runs the simulation model
#' @param maxit The number of iterations through the chained equations to settle on imputed values
#' @param maxwaves Maximum number of iterations (waves of simulations) before the algorithm is forced to stop
#' @param alpha Minimum fraction of n.experiments that gets used as input data for MICE
#' @param saturation.crit Saturation criterion. Not used at the moment.
#' @return A list with the following components: blah
#'
#' @examples
#' calibration.list <- MiceABC(targets = c(logit(0.15), log(0.015)),
#'                             n.experiments = 100,
#'                             lls = c(0, 0),
#'                             uls = c(10, 1),
#'                             probs = c(3, 7),
#'                             model = model.wrapper,
#'                             maxit = 20,
#'                             maxwaves = 10,
#'                             reps = 5,
#'                             alpha = 0.5,
#'                             saturation.crit = 0)
#' @import boot
#' @import randtoolbox
#' @import dplyr

# lls = c(1.01, 0.15, -1, 2, 0.05, 0.05, 10,  10,  -1,   0, -1,   -1,   -5, -2,   -0.2), #c(2, -1,   0.1, 2, 1, 1, 0.3), #c(0.1, -3, 0.1,  0,  0,  0, 0.2),
# uls = c(3,    0.35,  1, 5, 1,       1, 100, 100, -0.2, 5, -0.1, -0.1, -1, -0.2,    0), #c(4, -0.2, 0.6, 8, 6, 6, 0.95), #c(5, -0.1, 1.0, 10, 10, 10, 1.0),

MICE_ABC <- function(targets = targets.empirical,
                     n.experiments = 16,
                     lls = c(1.01, 0.15, -1, 2, 0.05, 0.05, 10,  10,  -1,   0, -1,   -1,   -5, -2,   -0.2), #c(2, -1,   0.1, 2, 1, 1, 0.3), #c(0.1, -3, 0.1,  0,  0,  0, 0.2),
                     uls = c(3,    0.35,  1, 5, 1,       1, 100, 100, -0.2, 5, -0.1, -0.1, -1, -0.2,    0), #c(4, -0.2, 0.6, 8, 6, 6, 0.95), #c(5, -0.1, 1.0, 10, 10, 10, 1.0),
                     # probs = c(3, 7), # indicate in the lls and uls vectors, which elements represent probabilities strictly between 0 and 1
                     model = simpact.wrapper,
                     maxit = 10,
                     maxwaves = 2, # and we can also try 10 to do the same as the 10 waves of Lenormand
                     #reps = 5,
                     alpha = 0.25,
                     saturation.crit = 0){
  ptm <- proc.time() # Start the clock

  range.width <- uls - lls
  ll.mat <- matrix(rep(lls, n.experiments), nrow = n.experiments, byrow = TRUE)
  range.width.mat <- matrix(rep(range.width, n.experiments), nrow = n.experiments, byrow = TRUE)
  sobol.seq.0.1 <- sobol(n = n.experiments, dim = length(lls), init = TRUE, scrambling = 1, seed = 1, normal = FALSE)
  experiments <- ll.mat + sobol.seq.0.1 * range.width.mat

  calibration.list <- list() # initiating the list where all the output of MiceABC will be stored

  wave <- 1 # initiating the loop of waves of simulations (one iteration is one wave)
  rel.dist.cutoff <- Inf # initially it is infinitely large, but in later iterations it shrinks
  saturation <- 0
  sim.results.with.design.df.selected <- NULL


  while (wave <= maxwaves && saturation == 0){
    print(c("wave", wave), quote = FALSE) # only for local testing. To be deleted after testing is completed
    # These large and small enough tests could be commented out because the "true" parameter value may lie outside the initial parameter range.
    # They may be useful to prevent infeasible MICE-derived proposals from running (e.g. positive parameters values that don't make sense)
    large.enough <- t(t(experiments) >= lls)
    small.enough <- t(t(experiments) <= uls)
    fine <- large.enough & small.enough
    fine.flag <- as.logical(rowSums(fine) == ncol(fine))
    experiments <- experiments[fine.flag, ]
    #try(if(nrow(experiments) == 0) stop("no experiments in the range"))
    if (nrow(experiments) < 2){
      wave <- maxwaves
    } else {
    print(c(nrow(experiments), "experiments to be run"), quote = FALSE) # only for local testing. To be deleted after testing is completed

    calibration.list$experiments.executed[[wave]] <- nrow(experiments)

    # Running the experiments
    # nl.path <- "/Applications/NetLogo 6.0/app/"
    # NLStart(nl.path, gui=FALSE, nl.obj="agemix.nl1", nl.jarname = 'netlogo-6.0.0.jar')
    # model.path <- "/Users/delvaw/Google Drive/IBM course/tutorials/AgeMixing/Calibration/dummy.0.1.nlogo"
    # NLLoadModel(model.path, nl.obj="agemix.nl1")
    # source(file = "/Users/delvaw/Google Drive/IBM course/tutorials/AgeMixing/Calibration/agemix.simulate.R")

    # sim.results.simple <- t(apply(experiments, 1, model, no.repeated.sim=reps)) # Maybe this apply() can be replaced by a faster foreach loop?
    # Replacing NetLogo model with parallel Simpact model

    # 1. Initially: n.experiments from sobol sequences. 7. Later, (1-alpha) fraction of n.experiments from MICE proposals
    sim.results.simple <- simpact.parallel(model = model,
                                           actual.input.matrix = experiments,
                                           seed_count = 0,
                                           n_cluster = 8)

    #experiment.number <- 1:nrow(experiments)
    sim.results.with.design.df <- as.data.frame(cbind(#experiment.number,
      experiments,
      sim.results.simple))
    x.names <- paste0("x.", seq(1:ncol(experiments)))
    y.names <- paste0("y.", seq(1:ncol(sim.results.simple)))
    names(sim.results.with.design.df) <- c(x.names, y.names)

    # 2. All n.experiments from sobol sequences get summarised measure of distance from target. 8. Distance measures for new proposals
    # NEEDS TO BE AUTOMATED

    # for(i in 1:6) { #-- Create objects  'r.1', 'r.2', ... 'r.6' --
    #   nam <- paste("r", i, sep = ".")
    #   assign(nam, 1:i)
    # }

    #browser()

    sum.sq.rel.dist <- rep(0, nrow(sim.results.with.design.df))
    for(i in 1:length(targets)) {
      nam <- paste0("y.", i, ".sq.rel.dist")
      val <- ((sim.results.simple[,i] - targets[i]) / targets[i])^2
      assign(nam, val)
      sum.sq.rel.dist <- sum.sq.rel.dist + get(nam)
    }

    sim.results.with.design.df$sum.sq.rel.dist <- sum.sq.rel.dist

    # sim.results.with.design.df$y.1.sq.rel.dist <- ((sim.results.with.design.df$y.1 - targets[1]) / targets[1])^2
    # sim.results.with.design.df$y.2.sq.rel.dist <- ((sim.results.with.design.df$y.2 - targets[2]) / targets[2])^2
    # sim.results.with.design.df$y.3.sq.rel.dist <- ((sim.results.with.design.df$y.3 - targets[3]) / targets[3])^2
    # sim.results.with.design.df$y.4.sq.rel.dist <- ((sim.results.with.design.df$y.4 - targets[4]) / targets[4])^2
    # sim.results.with.design.df$y.5.sq.rel.dist <- ((sim.results.with.design.df$y.5 - targets[5]) / targets[5])^2
    # sim.results.with.design.df$y.6.sq.rel.dist <- ((sim.results.with.design.df$y.6 - targets[6]) / targets[6])^2
    # sim.results.with.design.df$y.7.sq.rel.dist <- ((sim.results.with.design.df$y.7 - targets[7]) / targets[7])^2
    # sim.results.with.design.df$y.8.sq.rel.dist <- ((sim.results.with.design.df$y.8 - targets[8]) / targets[8])^2
    # sim.results.with.design.df$y.9.sq.rel.dist <- ((sim.results.with.design.df$y.9 - targets[9]) / targets[9])^2
    #
    #
    # sim.results.with.design.df$sum.sq.rel.dist <- sim.results.with.design.df$y.1.sq.rel.dist + sim.results.with.design.df$y.2.sq.rel.dist + sim.results.with.design.df$y.3.sq.rel.dist + sim.results.with.design.df$y.4.sq.rel.dist + sim.results.with.design.df$y.5.sq.rel.dist + sim.results.with.design.df$y.6.sq.rel.dist + sim.results.with.design.df$y.7.sq.rel.dist + sim.results.with.design.df$y.8.sq.rel.dist + sim.results.with.design.df$y.9.sq.rel.dist

    # 2b. Writing to list all the input and output of the executed experiments, so that we can plot it later
    calibration.list$sim.results.with.design[[wave]] <- sim.results.with.design.df

    # 10. Calculate fraction of new (1-alpha frac *n.experiments) distances that are below "old" distance threshold. NOT CURRENTLY USED BECAUSE saturation.crit = 0.
    below.old.treshold <- sim.results.with.design.df$sum.sq.rel.dist < rel.dist.cutoff
    frac.below.old.threshold <- sum(below.old.treshold %in% TRUE) / round(n.experiments * (1-alpha))
    if(frac.below.old.threshold < saturation.crit) saturation <- 1 # If less than the fraction saturation.crit of the new experiments a closer fit than the previous batch of retained experiments, the loop must be terminated.


    # 9. Merge with previously kept alpha fraction shortest distances
    sim.results.with.design.df <- rbind(sim.results.with.design.df.selected, sim.results.with.design.df) # Initially sim.results.with.design.df.selected = NULL

    # 3. Keeping alpha fraction shortest distances

    dist.order <- order(sim.results.with.design.df$sum.sq.rel.dist) # Ordering the squared distances from small to big. The last observation (targets) should be ranked first
    last.one.selected <- round(alpha * n.experiments) # OR IT COULD ALSO BE WRITTEN AS round(alpha * length(dist.order)) BUT IF RANGE CHECKS ARE ON, THE nrow IS NOT NECESSARILY = n.experiments
    selected.distances <- dist.order[1:last.one.selected]
    sim.results.with.design.df.selected <- sim.results.with.design.df[selected.distances, ] %>% filter(complete.cases(.)) # Here we filter out any experiments that gave NA output

    # 4. Set distance threshold at alpha-percentile best distance. 11. Update distance threshold
    rel.dist.cutoff <- max(sim.results.with.design.df.selected$sum.sq.rel.dist)  # We update the rel.dist.cutoff, based on the alpha fraction best parametersets

    # 5. Write alpha fraction of input-output-distances and distance threshold to list. 12. Write new best alpha fraction and updated threshold to list
    calibration.list$inputoutput[[wave]] <- sim.results.with.design.df.selected
    calibration.list$rel.dist.cutoff[[wave]] <- rel.dist.cutoff

    # 6. Feed MICE the best alpha fraction and ask for (1-alpha) fraction proposals based on targets
    targets.df <- as.data.frame(matrix(targets, ncol = length(targets))) # Putting target in dataframe format

    names(targets.df) <- y.names
    df.give.to.mice <- full_join(dplyr::select(sim.results.with.design.df.selected,
                                               -contains("sq.rel.dist")), # adding target to training dataset
                                 targets.df,
                                 by = names(targets.df)) # "by" statement added to avoid printing message of the variables were used for joining

    ### Before actually giving it to MICE, we need to logit transform the input params that represent probabilities, so that we can treat these variables as continuous variables
    # "probs" shows that it is the 3rd and 7th input variable (not counting the random seed variable)
    # df.give.to.mice[, probs] <- car::logit(df.give.to.mice[, probs]) # Commented out because none of the input parameters are probabilities for this age-mixing calibration
    # INSTEAD, we are transforming parameters that are necessarily strictly positive: sigma, gamma.a, gamma.b.
    # We could also consider a similar transformation for input parameters that we think should be negative (e.g. formation.hazard.agegapry.gap_factor_man_exp) but for now not yet
    strict.positive.params <- c(4:8)
    df.give.to.mice[, strict.positive.params] <- log(df.give.to.mice[, strict.positive.params])

    # Let's think a bit more carefully about which variables should be allowed as input for which input parameters.
    # IN THE FUTURE THIS COULD BE AUTOMATED WITH VARIABLE SELECTION ALGORITHMS.
    predictorMatrix <- (1 - diag(1, ncol(df.give.to.mice))) # This is the default matrix.
    # Let's now modify the first 15 rows of this matrix, corresponding to the indicators of predictor variables for the input variables. In brackets the values for the master model.
    # 1. hivtransmission.param.f1 (2)
    # 2. formation.hazard.agegapry.gap_agescale_man = formation.hazard.agegapry.gap_agescale_woman (0.25)
    # 3. person.agegap.man.dist.normal.mu = person.agegap.woman.dist.normal.mu (0)
    # 4. person.agegap.man.dist.normal.sigma = person.agegap.woman.dist.normal.sigma (3)
    # 5. person.person.eagerness.man.dist.gamma.a (0.23)
    # 6. person.person.eagerness.woman.dist.gamma.a (0.23)
    # 7. person.person.eagerness.man.dist.gamma.b (45)
    # 8. person.person.eagerness.woman.dist.gamma.b (45)
    # 9. formation.hazard.agegapry.gap_factor_man_exp = formation.hazard.agegapry.gap_factor_woman_exp (-0.5)
    # 10. formation.hazard.agegapry.baseline (0)
    # 11. formation.hazard.agegapry.numrel_man (-0.5)
    # 12. formation.hazard.agegapry.numrel_woman (-0.5)
    # 13. conception.alpha_base (-2.5)
    # 14. dissolution.alpha_0 (-0.52)
    # 15. dissolution.alpha_4 (-0.05))

    # And a reminder of the output vector
    # 1. AAD.male (4)
    # 2. SDAD.male (4)
    # 3. slope.male (0.7)
    # 4. WSD.male (3)
    # 5. BSD.male (1.5)
    # 6. intercept.male (-1)
    # 7. shape.nb.male (1.29)
    # 8. scale.nb.male (0.66)
    # 9. meandegree.male (1)
    # 10. pp.cp.6months.male (0.134)
    # 11. hiv.prev.lt25.women (0.2)
    # 12. hiv.prev.lt25.men (0.1)
    # 13. hiv.prev.25.34.women (0.42)
    # 14. hiv.prev.25.34.men (0.17)
    # 15. hiv.prev.35.44.women (0.37)
    # 16. hiv.prev.35.44.men (0.24)
    # 17. exp(growthrate)) (1.01)

    x.offset <- length(x.names)


    predictorMatrix[1:length(x.names), ] <- 0 # First we "empty" the relevant rows, then we refill them.
    # We are currently not allowing input variables to be predicted by other predictor variables. Only via output variables. We could change this at a later stage.
    predictorMatrix[1, x.offset + 11:12] <- 1 # relative susceptibility in young women is predicted by HIV prevalence in young men and women
    predictorMatrix[2, x.offset + 3] <- 1 # agescale predicted by slope
    predictorMatrix[3, x.offset + c(1, 3, 6)] <- 1 # mean of the person-specific age gap preferences is predicted by slope, intercept and AAD
    predictorMatrix[4, x.offset + c(2, 4, 5)] <- 1 # sd of the person-specific age gap preferences is predicted by SD, WSD, BSD
    predictorMatrix[5, x.offset + c(7, 8, 10, 14, 17)] <- 1 # man gamma a predicted by gamma shape.male, scale.male, pp.cp, hiv.prev.25.34.men, exp(growthrate)
    predictorMatrix[6, x.offset + c(7, 8, 10, 13, 17)] <- 1 # woman gamma a predicted by gamma shape.male, scale.male, pp.cp, hiv.prev.25.34.women, exp(growthrate)
    predictorMatrix[7, x.offset + c(7, 8, 10, 14, 17)] <- 1 # man gamma b predicted by gamma shape.male, scale.male, pp.cp, hiv.prev.25.34.men, exp(growthrate)
    predictorMatrix[8, x.offset + c(7, 8, 10, 13, 17)] <- 1 # woman gamma b predicted by gamma shape.male, scale.male, pp.cp, hiv.prev.25.34.men, exp(growthrate)
    predictorMatrix[9, x.offset + c(2, 4, 5, 7, 8, 15, 16, 17)] <- 1 # formation.hazard.agegapry.gap_factor_x_exp is predicted by population growth, age gap variance, hiv prevalence,
    predictorMatrix[10, x.offset + c(7, 8, 10, 13, 14, 17)] <- 1 # baseline formation hazard predicted by HIV prevalence, cp, degree distrib. HIV prevalence.
    predictorMatrix[11, x.offset + c(7, 8, 10, 13, 14, 17)] <- 1 # numrel man penalty is predicted by degree distrib, cp, prev, popgrowth
    predictorMatrix[12, x.offset + c(7, 8, 10, 13, 14, 17)] <- 1 # # numrel woman penalty is predicted by degree distrib, cp, prev, popgrowth
    predictorMatrix[13, x.offset + 17] <- 1 # conception.alpha_base is predicted by popgrowth
    predictorMatrix[14, x.offset + c(7, 8, 10, 17)] <- 1 # baseline dissolution hazard predicted by degree distrib, cp, popgrowth
    predictorMatrix[15, x.offset + c(7, 8, 10, 17)] <- 1 # age effect on dissolution hazard predicted by degree distrib, cp, popgrowth, HIV prev in older people (maybe?)

    # NOTE: As it stands, each output statistic is predicted by ALL input and ALL other output statistics. That may not be a great idea, or even possible, if there is collinearity.

    # Because there are some many interaction terms, let's first try without any
    #df.give.to.mice <- cbind(df.give.to.mice, data.frame(y.1.y.2 = df.give.to.mice$y.1 * df.give.to.mice$y.2)) # adding interaction term
    print(c(nrow(df.give.to.mice), "nrows to give to mice"), quote = FALSE)

    mice.test <- mice(data = df.give.to.mice, # the dataframe with missing data
                      m = round(n.experiments * (1-alpha)), # number of imputations
                      predictorMatrix = predictorMatrix,
                      method = "norm",
                      defaultMethod = "norm",
                      maxit = maxit,
                      printFlag = FALSE,
                      data.init = NULL)



    # experiments <- cbind(unlist(mice.test$imp$x.1), # these are the (1-alpha) fraction of n.experiments proposals
    #                      unlist(mice.test$imp$x.2))

    experiments <- unlist(mice.test$imp) %>% matrix(., byrow = FALSE, ncol = length(x.names))
    # Before we check the suitability of the new experimental input parameter values, we must backtransform the logits to probabilities
    # experiments[, probs] <- boot::inv.logit(experiments[, probs])
    # Before we check the suitability of the new experimental input parameter values, we must backtransform the log values to natural values
    experiments[, strict.positive.params] <- exp(experiments[, strict.positive.params])
    }
    wave <- wave + 1
  }

  calibration.list$secondspassed <- proc.time() - ptm # Stop the clock
  return(calibration.list)
}






























MiceABC <- function(targets = c(3.75, 4.48),
                    n.experiments = 16,
                    lls = c(-6.5, -6.5),
                    uls = c(-5.5, -5.5),
                    model = simpact.wrapper,
                    maxit = 20,
                    maxwaves = 2,#0,
                    #reps = 5,
                    alpha = 0.5,
                    saturation.crit = 0.1){
  ptm <- proc.time() # Start the clock

  range.width <- uls - lls
  ll.mat <- matrix(rep(lls, n.experiments), nrow = n.experiments, byrow = TRUE)
  range.width.mat <- matrix(rep(range.width, n.experiments), nrow = n.experiments, byrow = TRUE)
  sobol.seq.0.1 <- sobol(n = n.experiments, dim = length(lls), init = TRUE, scrambling = 1, seed = 1, normal = FALSE)
  experiments <- ll.mat + sobol.seq.0.1 * range.width.mat

  calibration.list <- list() # initiating the list where all the output of MiceABC will be stored

  wave <- 1 # initiating the loop of waves of simulations (one iteration is one wave)
  rel.dist.cutoff <- Inf # initially it is infinitely large, but in later iterations it shrinks
  saturation <- 0
  sim.results.with.design.df.selected <- NULL


  while (wave <= maxwaves && saturation == 0){
    print(c("wave", wave), quote = FALSE) # only for local testing. To be deleted after testing is completed
    # These large and small enough tests could be commented out because the "true" parameter value may lie outside the initial parameter range.
    # They may be useful to prevent infeasible MICE-derived proposals from running (e.g. positive parameters values that don't make sense)
    large.enough <- t(t(experiments) >= lls)
    small.enough <- t(t(experiments) <= uls)
    fine <- large.enough & small.enough
    fine.flag <- as.logical(rowSums(fine) == ncol(fine))
    experiments <- experiments[fine.flag, ]
    try(if(nrow(experiments) == 0) stop("no experiments in the range"))
    print(c(nrow(experiments), "experiments to be run"), quote = FALSE) # only for local testing. To be deleted after testing is completed

    calibration.list$experiments.executed[[wave]] <- nrow(experiments)

    # Running the experiments
    # nl.path <- "/Applications/NetLogo 6.0/app/"
    # NLStart(nl.path, gui=FALSE, nl.obj="agemix.nl1", nl.jarname = 'netlogo-6.0.0.jar')
    # model.path <- "/Users/delvaw/Google Drive/IBM course/tutorials/AgeMixing/Calibration/dummy.0.1.nlogo"
    # NLLoadModel(model.path, nl.obj="agemix.nl1")
    # source(file = "/Users/delvaw/Google Drive/IBM course/tutorials/AgeMixing/Calibration/agemix.simulate.R")

    # sim.results.simple <- t(apply(experiments, 1, model, no.repeated.sim=reps)) # Maybe this apply() can be replaced by a faster foreach loop?
    # Replacing NetLogo model with parallel Simpact model

    # 1. Initially: n.experiments from sobol sequences. 7. Later, (1-alpha) fraction of n.experiments from MICE proposals
    sim.results.simple <- simpact.parallel(model = model,
                                           actual.input.matrix = experiments,
                                           seed_count = 0,
                                           n_cluster = 8)

    #experiment.number <- 1:nrow(experiments)
    sim.results.with.design.df <- as.data.frame(cbind(#experiment.number,
                                                      experiments,
                                                      sim.results.simple))
    x.names <- paste0("x.", seq(1:ncol(experiments)))
    y.names <- paste0("y.", seq(1:ncol(sim.results.simple)))
    names(sim.results.with.design.df) <- c(x.names, y.names)

    # 2. All n.experiments from sobol sequences get summarised measure of distance from target. 8. Distance measures for new proposals
    sim.results.with.design.df$y.1.sq.rel.dist <- ((sim.results.with.design.df$y.1 - targets[1]) / targets[1])^2
    sim.results.with.design.df$y.2.sq.rel.dist <- ((sim.results.with.design.df$y.2 - targets[2]) / targets[2])^2
    sim.results.with.design.df$sum.sq.rel.dist <- sim.results.with.design.df$y.1.sq.rel.dis + sim.results.with.design.df$y.2.sq.rel.dis

    # 2b. Writing to list all the input and output of the executed experiments, so that we can plot it later
    calibration.list$sim.results.with.design[[wave]] <- sim.results.with.design.df

    # 10. Calculate fraction of new (1-alpha frac *n.experiments) distances that are below "old" distance threshold
    below.old.treshold <- sim.results.with.design.df$sum.sq.rel.dist < rel.dist.cutoff
    frac.below.old.threshold <- sum(below.old.treshold) / round(n.experiments * (1-alpha))
    if(frac.below.old.threshold < saturation.crit) saturation <- 1 # If less than the fraction saturation.crit of the new experiments a closer fit than the previous batch of retained experiments, the loop must be terminated.


    # 9. Merge with previously kept alpha fraction shortest distances
    sim.results.with.design.df <- rbind(sim.results.with.design.df.selected, sim.results.with.design.df) # Initially sim.results.with.design.df.selected = NULL

    # 3. Keeping alpha fraction shortest distances
    dist.order <- order(sim.results.with.design.df$sum.sq.rel.dist) # Ordering the squared distances from small to big. The last observation (targets) should be ranked first
    last.one.selected <- round(alpha * n.experiments) # OR IT COULD ALSO BE WRITTEN AS round(alpha * length(dist.order)) BUT IF RANGE CHECKS ARE ON, THE nrow IS NOT NECESSARILY = n.experiments
    selected.distances <- dist.order[1:last.one.selected]
    sim.results.with.design.df.selected <- sim.results.with.design.df[selected.distances, ] %>% filter(complete.cases(.)) # Here we filter out any experiments that gave NA output

    # 4. Set distance threshold at alpha-percentile best distance. 11. Update distance threshold
    rel.dist.cutoff <- max(sim.results.with.design.df.selected$sum.sq.rel.dist)  # We update the rel.dist.cutoff, based on the alpha fraction best parametersets

    # 5. Write alpha fraction of input-output-distances and distance threshold to list. 12. Write new best alpha fraction and updated threshold to list
    calibration.list$inputoutput[[wave]] <- sim.results.with.design.df.selected
    calibration.list$rel.dist.cutoff[[wave]] <- rel.dist.cutoff

    # 6. Feed MICE the best alpha fraction and ask for (1-alpha) fraction proposals based on targets
    targets.df <- as.data.frame(matrix(targets, ncol = length(targets))) # Putting target in dataframe format
    names(targets.df) <- y.names
    df.give.to.mice <- full_join(dplyr::select(sim.results.with.design.df.selected,
                                               -contains("sq.rel.dist")), # adding target to training dataset
                                 targets.df)
    df.give.to.mice <- cbind(df.give.to.mice, data.frame(y.1.y.2 = df.give.to.mice$y.1 * df.give.to.mice$y.2)) # adding interaction term
    print(c(nrow(df.give.to.mice), "nrows to give to mice"), quote = FALSE)
    mice.test <- mice(data = df.give.to.mice, # the dataframe with missing data
                      m = round(n.experiments * (1-alpha)), # number of imputations
                      method = "norm",
                      defaultMethod = "norm",
                      maxit = maxit,
                      printFlag = FALSE,
                      data.init = NULL)

    experiments <- cbind(unlist(mice.test$imp$x.1), # these are the (1-alpha) fraction of n.experiments proposals
                         unlist(mice.test$imp$x.2))

    wave <- wave + 1
  }

  calibration.list$secondspassed <- proc.time() - ptm # Stop the clock
  return(calibration.list)
}






MiceABC.V8 <- function(targets = c(5.816, 3.924, 0.05916, 0.6628, 5.229, 1.4985, 0.12573, 0.003311),
                       n.experiments = 16,
                       lls = c(0.1, -3),
                       uls = c(5, -0.1),
                       model = simpact.wrapper,
                       maxit = 20,
                       maxwaves = 2,#0,
                       #reps = 5,
                       alpha = 0.25,
                       saturation.crit = 0){
  ptm <- proc.time() # Start the clock

  range.width <- uls - lls
  ll.mat <- matrix(rep(lls, n.experiments), nrow = n.experiments, byrow = TRUE)
  range.width.mat <- matrix(rep(range.width, n.experiments), nrow = n.experiments, byrow = TRUE)
  sobol.seq.0.1 <- sobol(n = n.experiments, dim = length(lls), init = TRUE, scrambling = 1, seed = 1, normal = FALSE)
  experiments <- ll.mat + sobol.seq.0.1 * range.width.mat

  calibration.list <- list() # initiating the list where all the output of MiceABC will be stored

  wave <- 1 # initiating the loop of waves of simulations (one iteration is one wave)
  rel.dist.cutoff <- Inf # initially it is infinitely large, but in later iterations it shrinks
  saturation <- 0
  sim.results.with.design.df.selected <- NULL


  while (wave <= maxwaves && saturation == 0){
    print(c("wave", wave), quote = FALSE) # only for local testing. To be deleted after testing is completed
    # These large and small enough tests could be commented out because the "true" parameter value may lie outside the initial parameter range.
    # They may be useful to prevent infeasible MICE-derived proposals from running (e.g. positive parameters values that don't make sense)
    large.enough <- t(t(experiments) >= lls)
    small.enough <- t(t(experiments) <= uls)
    fine <- large.enough & small.enough
    fine.flag <- as.logical(rowSums(fine) == ncol(fine))
    experiments <- experiments[fine.flag, ]
    try(if(nrow(experiments) == 0) stop("no experiments in the range"))
    print(c(nrow(experiments), "experiments to be run"), quote = FALSE) # only for local testing. To be deleted after testing is completed

    calibration.list$experiments.executed[[wave]] <- nrow(experiments)

    # Running the experiments
    # nl.path <- "/Applications/NetLogo 6.0/app/"
    # NLStart(nl.path, gui=FALSE, nl.obj="agemix.nl1", nl.jarname = 'netlogo-6.0.0.jar')
    # model.path <- "/Users/delvaw/Google Drive/IBM course/tutorials/AgeMixing/Calibration/dummy.0.1.nlogo"
    # NLLoadModel(model.path, nl.obj="agemix.nl1")
    # source(file = "/Users/delvaw/Google Drive/IBM course/tutorials/AgeMixing/Calibration/agemix.simulate.R")

    # sim.results.simple <- t(apply(experiments, 1, model, no.repeated.sim=reps)) # Maybe this apply() can be replaced by a faster foreach loop?
    # Replacing NetLogo model with parallel Simpact model

    # 1. Initially: n.experiments from sobol sequences. 7. Later, (1-alpha) fraction of n.experiments from MICE proposals
    sim.results.simple <- simpact.parallel(model = model,
                                           actual.input.matrix = experiments,
                                           seed_count = 0,
                                           n_cluster = 8)

    #experiment.number <- 1:nrow(experiments)
    sim.results.with.design.df <- as.data.frame(cbind(#experiment.number,
      experiments,
      sim.results.simple))
    x.names <- paste0("x.", seq(1:ncol(experiments)))
    y.names <- paste0("y.", seq(1:ncol(sim.results.simple)))
    names(sim.results.with.design.df) <- c(x.names, y.names)

    # 2. All n.experiments from sobol sequences get summarised measure of distance from target. 8. Distance measures for new proposals
    # NEEDS TO BE AUTOMATED
    sim.results.with.design.df$y.1.sq.rel.dist <- ((sim.results.with.design.df$y.1 - targets[1]) / targets[1])^2
    sim.results.with.design.df$y.2.sq.rel.dist <- ((sim.results.with.design.df$y.2 - targets[2]) / targets[2])^2
    sim.results.with.design.df$y.3.sq.rel.dist <- ((sim.results.with.design.df$y.3 - targets[3]) / targets[3])^2
    sim.results.with.design.df$y.4.sq.rel.dist <- ((sim.results.with.design.df$y.4 - targets[4]) / targets[4])^2
    sim.results.with.design.df$y.5.sq.rel.dist <- ((sim.results.with.design.df$y.5 - targets[5]) / targets[5])^2
    sim.results.with.design.df$y.6.sq.rel.dist <- ((sim.results.with.design.df$y.6 - targets[6]) / targets[6])^2
    sim.results.with.design.df$y.7.sq.rel.dist <- ((sim.results.with.design.df$y.7 - targets[7]) / targets[7])^2
    sim.results.with.design.df$y.8.sq.rel.dist <- ((sim.results.with.design.df$y.8 - targets[8]) / targets[8])^2

    sim.results.with.design.df$sum.sq.rel.dist <- sim.results.with.design.df$y.1.sq.rel.dist + sim.results.with.design.df$y.2.sq.rel.dist + sim.results.with.design.df$y.3.sq.rel.dist + sim.results.with.design.df$y.4.sq.rel.dist + sim.results.with.design.df$y.5.sq.rel.dist + sim.results.with.design.df$y.6.sq.rel.dist + sim.results.with.design.df$y.7.sq.rel.dist + sim.results.with.design.df$y.8.sq.rel.dist

    # 2b. Writing to list all the input and output of the executed experiments, so that we can plot it later
    calibration.list$sim.results.with.design[[wave]] <- sim.results.with.design.df

    # 10. Calculate fraction of new (1-alpha frac *n.experiments) distances that are below "old" distance threshold
    below.old.treshold <- sim.results.with.design.df$sum.sq.rel.dist < rel.dist.cutoff
    frac.below.old.threshold <- sum(below.old.treshold %in% TRUE) / round(n.experiments * (1-alpha))
    if(frac.below.old.threshold < saturation.crit) saturation <- 1 # If less than the fraction saturation.crit of the new experiments a closer fit than the previous batch of retained experiments, the loop must be terminated.


    # 9. Merge with previously kept alpha fraction shortest distances
    sim.results.with.design.df <- rbind(sim.results.with.design.df.selected, sim.results.with.design.df) # Initially sim.results.with.design.df.selected = NULL

    # 3. Keeping alpha fraction shortest distances
    dist.order <- order(sim.results.with.design.df$sum.sq.rel.dist) # Ordering the squared distances from small to big. The last observation (targets) should be ranked first
    last.one.selected <- round(alpha * n.experiments) # OR IT COULD ALSO BE WRITTEN AS round(alpha * length(dist.order)) BUT IF RANGE CHECKS ARE ON, THE nrow IS NOT NECESSARILY = n.experiments
    selected.distances <- dist.order[1:last.one.selected]
    sim.results.with.design.df.selected <- sim.results.with.design.df[selected.distances, ] %>% filter(complete.cases(.)) # Here we filter out any experiments that gave NA output

    # 4. Set distance threshold at alpha-percentile best distance. 11. Update distance threshold
    rel.dist.cutoff <- max(sim.results.with.design.df.selected$sum.sq.rel.dist)  # We update the rel.dist.cutoff, based on the alpha fraction best parametersets

    # 5. Write alpha fraction of input-output-distances and distance threshold to list. 12. Write new best alpha fraction and updated threshold to list
    calibration.list$inputoutput[[wave]] <- sim.results.with.design.df.selected
    calibration.list$rel.dist.cutoff[[wave]] <- rel.dist.cutoff

    # 6. Feed MICE the best alpha fraction and ask for (1-alpha) fraction proposals based on targets
    targets.df <- as.data.frame(matrix(targets, ncol = length(targets))) # Putting target in dataframe format

    names(targets.df) <- y.names
    df.give.to.mice <- full_join(dplyr::select(sim.results.with.design.df.selected,
                                               -contains("sq.rel.dist")), # adding target to training dataset
                                 targets.df)

    # Because there are some many interaction terms, let's first try without any
    #df.give.to.mice <- cbind(df.give.to.mice, data.frame(y.1.y.2 = df.give.to.mice$y.1 * df.give.to.mice$y.2)) # adding interaction term
    print(c(nrow(df.give.to.mice), "nrows to give to mice"), quote = FALSE)
    mice.test <- mice(data = df.give.to.mice, # the dataframe with missing data
                      m = round(n.experiments * (1-alpha)), # number of imputations
                      method = "norm",
                      defaultMethod = "norm",
                      maxit = maxit,
                      printFlag = FALSE,
                      data.init = NULL)

    experiments <- cbind(unlist(mice.test$imp$x.1), # these are the (1-alpha) fraction of n.experiments proposals
                         unlist(mice.test$imp$x.2))

    wave <- wave + 1
  }

  calibration.list$secondspassed <- proc.time() - ptm # Stop the clock
  return(calibration.list)
}














    df <- dplyr::filter(sim.results.with.design.df,
                        sum.sq.rel.dist <= rel.dist.cutoff) # We force the new batch to be at least as good as the previous one. Initially, rel.dist.cutoff is Inf

    if(nrow(df) < saturation.crit * nrow(sim.results.with.design.df)) saturation <- 1 # If less than the fraction saturation.crit of the new experiments a closer fit than the previous batch of retained experiments, the loop must be terminated.


    calibration.list$inputoutput[[wave]] <- df # sim.results.with.design.df
    calibration.list$inputoutput[[wave]] <- df.selected # sim.results.with.design.df

    all.sim.results.with.design.df <- do.call(rbind, calibration.list$inputoutput) %>% filter(complete.cases(.)) # Some of the runs may have produced NA values.



    dist.order <- order(df$sum.sq.rel.dist) # Ordering the squared distances from small to big. The last observation (targets) should be ranked first
    shortest.dist.percentile <- alpha # 0.25
    last.one.selected <- round(shortest.dist.percentile * length(dist.order))
    selected.distances <- dist.order[1:(1+last.one.selected)] # we are adding one because the targets row is also included.
    df.selected <- df[selected.distances, ]


    rel.dist.cutoff <- max(df.selected$sum.sq.rel.dist)  # We update the rel.dist.cutoff, based on the alpha fraction best parametersets
    print(c("relative distance cut-off (sum of squared relative distances)", rel.dist.cutoff), quote = FALSE)

    calibration.list$rel.dist.cutoff[[wave]] <- rel.dist.cutoff







    all.sim.results.with.design.df <- do.call(rbind, calibration.list$inputoutput) %>% filter(complete.cases(.)) # Some of the runs may have produced NA values.

    df.select.for.mice <- dplyr::filter(all.sim.results.with.design.df,
                                 sum.sq.rel.dist <= rel.dist.cutoff)


    ##
    # Now we combine this new set of simulation inputs and outputs with the target statistics

    targets.df <- as.data.frame(matrix(targets, ncol = length(targets)))
    names(targets.df) <- y.names

    df.give.to.mice <- full_join(df.select.for.mice, #all.sim.results.with.design.df,
                                 targets.df)
    df.give.to.mice <- cbind(df.give.to.mice, data.frame(y.1.y.2 = df.give.to.mice$y.1 * df.give.to.mice$y.2))

    print(c(nrow(df.give.to.mice), "nrows to give to mice"), quote = FALSE)

    mice.test <- mice(data = df.give.to.mice, # the dataframe with missing data
                      m = round(n.experiments * (1-alpha)), # number of imputations
                      method = "norm",
                      defaultMethod = "norm",
                      maxit = maxit,
                      printFlag = FALSE,
                      data.init = NULL)

    experiments <- cbind(unlist(mice.test$imp$x.1),
                         unlist(mice.test$imp$x.2))
    ##

    wave <- wave + 1

  }
  calibration.list$secondspassed <- proc.time() - ptm # Stop the clock
  return(calibration.list)
}

    # # Suppose we first filter the data: We only retain observations within x% from the targets
    # rel.tol <- rel.dist.tol
    # # We also calculate sum of squared relative distances
    # df$y.1.sq.rel.dist <- ((df$y.1 - targets[1]) / targets[1])^2
    # df$y.2.sq.rel.dist <- ((df$y.2 - targets[2]) / targets[2])^2
    # df$sum.sq.rel.dist <- df$y.1.sq.rel.dis + df$y.2.sq.rel.dis

    df <- dplyr::filter(df,
                        sum.sq.rel.dist <= rel.dist.cutoff) # We force the new batch to be at least as good as the previous one. Initially, rel.dist.cutoff is Inf

    dist.order <- order(df$sum.sq.rel.dist) # Ordering the squared distances from small to big. The last observation (targets) should be ranked first
    shortest.dist.percentile <- alpha # 0.25
    last.one.selected <- round(shortest.dist.percentile * length(dist.order))
    selected.distances <- dist.order[1:(1+last.one.selected)] # we are adding one because the targets row is also included.
    df.selected <- df[selected.distances, ]


    rel.dist.cutoff <- max(df.selected$sum.sq.rel.dist)
    print(c("relative distance cut-off (sum of squared relative distances)", rel.dist.cutoff), quote = FALSE)

    calibration.list$rel.dist.cutoff[[wave]] <- rel.dist.cutoff

    # df.selected <- dplyr::filter(df,
    #                              sum.sq.rel.dist <= rel.dist.cutoff)


    #df.filtered <- filter(df,
    #                      abs(df$y.1 - targets[1]) / targets[1] <= rel.tol,
    #                      abs(df$y.2 - targets[2]) / targets[2] <= rel.tol)

    #df.choices <-list(df.selected, df.filtered)
    #df.choice.index <- ifelse(nrow(df.selected) > nrow(df.filtered),
    #                          1,
    #                          2)
    df.give.to.mice <- df.selected #df.choices[[df.choice.index]]

    interactions.df <- data.frame(y.1.y.2 = df.give.to.mice$y.1 * df.give.to.mice$y.2)

    df.give.to.mice <- cbind(df.give.to.mice[, 1:ncol(all.sim.results.with.design.df)],
                             interactions.df)
    print(c(nrow(df.give.to.mice), "nrows to give to mice"), quote = FALSE)

    mice.test <- mice(data = df.give.to.mice, # the dataframe with missing data
                      m = round(n.experiments * (1-alpha)), # number of imputations
                      method = "norm",
                      defaultMethod = "norm",
                      maxit = maxit,
                      printFlag = FALSE,
                      data.init = NULL)

    experiments <- cbind(unlist(mice.test$imp$x.1),
                         unlist(mice.test$imp$x.2))
    ##

    wave <- wave + 1

  }
  calibration.list$secondspassed <- proc.time() - ptm # Stop the clock
  return(calibration.list)
}




    age.gap.tol.intercept <- experiments[, 1] # guess.1.df$age.gap.tol.intercept[x];
    age.gap.tol.coef <- experiments[, 2] #guess.1.df$age.gap.tol.coef[x];

    experiment <- 1:nrow(experiments)

    sim.results.with.design.df <- as.data.frame(cbind(sim.results.simple, experiment, age.gap.tol.intercept, age.gap.tol.coef))

    calibration.list$inputoutput[[iteration]] <- sim.results.with.design.df
    calibration.list$intercept[[iteration]] <- summary(experiments[ , 1])
    calibration.list$coef[[iteration]] <- summary(experiments[, 2])

    all.sim.results.with.design.df <- do.call(rbind, calibration.list$inputoutput) # sim.results.with.design.df # do.call(rbind, calibration.list$inputoutput)


    # Now we combine this new set of simulation inputs and outputs with the target statistics
    df <- dplyr::select(full_join(all.sim.results.with.design.df, data.frame(cumsum.transmissions = targets[1],
                                                                         mean.age.hivpos = targets[2],
                                                                         var.age.hivpos = targets[3])),
                        age.gap.tol.intercept,
                        age.gap.tol.coef,
                        cumsum.transmissions,
                        mean.age.hivpos,
                        var.age.hivpos)

    # Suppose we first filter the data: We only retain observations within x% from the targets
    rel.tol <- rel.dist.tol
    # We also calculate sum of squared relative distances
    df$cumsum.transmissions.sq.rel.dist <- ((df$cumsum.transmissions - targets[1]) / targets[1])^2
    df$mean.age.hivpos.sq.rel.dist <- ((df$mean.age.hivpos - targets[2]) / targets[2])^2
    df$var.age.hivpos.sq.rel.dist <- ((df$var.age.hivpos - targets[3]) / targets[3])^2
    df$sum.sq.rel.dist <- df$cumsum.transmissions.sq.rel.dist + df$mean.age.hivpos.sq.rel.dist +df$var.age.hivpos.sq.rel.dist

    df <- dplyr::filter(df,
                        sum.sq.rel.dist <= rel.dist.cutoff) # We force the new batch to be at least as good as the previous one

    dist.order <- order(df$sum.sq.rel.dist) # Ordering the squared distances from small to big. The last observation (targets) should be ranked first
    shortest.dist.percentile <- alpha # 0.25
    last.one.selected <- round(shortest.dist.percentile * length(dist.order))
    selected.distances <- dist.order[1:last.one.selected]
    df.selected <- df[selected.distances, ]

    rel.dist.cutoff <- df$sum.sq.rel.dist[dist.order[last.one.selected]]
    print(c("relative distance cut-off (sum of squared relative distances)", rel.dist.cutoff), quote = FALSE)

    calibration.list$rel.dist.cutoff[[iteration]] <- rel.dist.cutoff

    df.selected <- dplyr::filter(df,
                        sum.sq.rel.dist <= rel.dist.cutoff)


    df.filtered <- filter(df,
                          abs(df$cumsum.transmissions - targets[1]) / targets[1] <= rel.tol,
                          abs(df$mean.age.hivpos - targets[2]) / targets[2] <= rel.tol,
                          abs(df$var.age.hivpos - targets[3]) / targets[3] <= rel.tol)

    df.choices <-list(df.selected, df.filtered)
    df.choice.index <- ifelse(nrow(df.selected) > nrow(df.filtered),
                              1,
                              2)
    df.give.to.mice <- df.choices[[df.choice.index]]

    interactions.df <- data.frame(x1x2 = df.give.to.mice$cumsum.transmissions * df.give.to.mice$mean.age.hivpos,
                                  x1x3 = df.give.to.mice$cumsum.transmissions * df.give.to.mice$var.age.hivpos,
                                  x2x3 = df.give.to.mice$mean.age.hivpos * df.give.to.mice$var.age.hivpos)

    df.give.to.mice <- cbind(dplyr::select(df.give.to.mice,
                                           age.gap.tol.intercept,
                                           age.gap.tol.coef,
                                           cumsum.transmissions,
                                           mean.age.hivpos,
                                           var.age.hivpos),
                             interactions.df)
    print(c(nrow(df.give.to.mice), "nrows to give to mice"), quote = FALSE)

    mice.test <- mice(data = df.give.to.mice, # the dataframe with missing data
                      m = n.experiments, # number of imputations
                      method = "norm",
                      defaultMethod = "norm",
                      maxit = maxit,
                      printFlag = FALSE,
                      data.init = NULL)

    experiments <- cbind(unlist(mice.test$imp$age.gap.tol.intercept),
                         unlist(mice.test$imp$age.gap.tol.coef))
    wave <- wave + 1
  }
  calibration.list$secondspassed <- proc.time() - ptm # Stop the clock
  return(calibration.list)
}
