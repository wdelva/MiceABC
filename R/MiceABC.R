#' Adaptive Population ABC with MICE-based supervised learning
#'
#' This function implements an adaptive population ABC scheme with supervised learning through Multivariate Imputation by Chained Equations (MICE) to produce new particles (i.e. samples from the parameter space and their associated output statistics)
#'
#' @param targets The vector of target statistics to calibrate the simulation model against
#' @param n.experiments The number of initial multivariate samples from the parameter space
#' @param lls Lower limits of the parameter ranges to be explored
#' @param uls Upper limits of the parameter ranges to be explored
#' @param model Handle to the wrapper function that runs the simulation model
#' @param maxit The number of iterations through the chained equations to settle on imputed values
#' @param maxwaves Maximum number of iterations (waves of simulations) before the algorithm is forced to stop
#' @param reps Number of times each parameter set is repeated
#' @param alpha Minimum fraction of n.experiments that gets used as input data for MICE
#' @param rel.dist.tol Maximum tolerated relative distance from target statistics (may be overwritten by alpha)
#' @return A list with the following components: blah
#'
#' @examples
#' calibration.list <- MiceABC(targets = c(logit(0.15), log(0.015)),
#'                             n.experiments = 100,
#'                             lls = c(0, 0),
#'                             uls = c(10, 1),
#'                             model = model.wrapper,
#'                             maxit = 20,
#'                             maxwaves = 10,
#'                             reps = 5,
#'                             alpha = 0.5,
#'                             rel.dist.tol = 0.1)
#' @import boot
#' @import randtoolbox
#' @import dplyr

MiceABC <- function(targets = c(4.3, 4.7),
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
