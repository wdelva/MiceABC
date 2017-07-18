MaC.incremental <- function(targets.empirical = targets.empirical,
                            RMSD.tol.max = 2,
                            min.givetomice = 64,
                            n.experiments = 256,
                            lls = c(1,  0.12, -0.3, 2.5, 0.1, 0.1, 20, 20, -0.8, 2, -0.35, -0.35, -3.6, -0.8, -0.16),
                            uls = c(1.2, 0.37, 0.3, 3.5, 0.4, 0.4, 66, 66, -0.25, 3.9, -0.1, -0.1, -1.4, -0.3,  -0.001),
                            model = simpact.wrapper,
                            maxit = 50,
                            maxwaves = 4){
  # 0. Start the clock
  ptm <- proc.time()
  calibration.list <- list() # initiating the list where all the output of MiceABC will be stored
  wave <- 1 # initiating the loop of waves of simulations (one iteration is one wave)
  rel.dist.cutoff <- Inf # initially it is infinitely large, but in later iterations it shrinks
  sim.results.with.design.df <- NULL # Will be growing with each wave (appending)
  sim.results.with.design.df.selected <- NULL
  final.intermediate.features <- NULL

  # 1. Start loop of waves, based on comparing intermediate features with targets.empirical
  while (wave <= maxwaves & !identical(final.intermediate.features, targets.empirical)){
    if (wave == 1){
      # 2. Initial, naive results, based on Sobol sequences
      range.width <- uls - lls
      ll.mat <- matrix(rep(lls, n.experiments), nrow = n.experiments, byrow = TRUE)
      range.width.mat <- matrix(rep(range.width, n.experiments), nrow = n.experiments, byrow = TRUE)
      sobol.seq.0.1 <- sobol(n = n.experiments, dim = length(lls), init = TRUE, scrambling = 1, seed = 1, normal = FALSE)
      experiments <- ll.mat + sobol.seq.0.1 * range.width.mat
    }

    sim.results.simple <- simpact.parallel(model = model,
                                           actual.input.matrix = experiments,
                                           seed_count = 0,
                                           n_cluster = 8)
    # save(sim.results.simple, file = "/Users/delvaw/Documents/MiceABC/sim.results.simple.RData")
    # load(file = "sim.results.simple.RData")

    new.sim.results.with.design.df <- as.data.frame(cbind(experiments,
                                                          sim.results.simple))
    x.names <- paste0("x.", seq(1:ncol(experiments)))
    y.names <- paste0("y.", seq(1:ncol(sim.results.simple)))
    x.offset <- length(x.names)
    names(new.sim.results.with.design.df) <- c(x.names, y.names)

    new.sim.results.with.design.df <- new.sim.results.with.design.df %>% dplyr::filter(complete.cases(.))

    if (is.null(sim.results.with.design.df)){ # TRUE for the first wave only
      sim.results.with.design.df <- rbind(sim.results.with.design.df,
                                          new.sim.results.with.design.df)
    } else {
      sim.results.with.design.df <- rbind(dplyr::select(sim.results.with.design.df,
                                                        -contains("RMSD")),
                                          new.sim.results.with.design.df)
    }


    # save(sim.results.with.design.df, file = "/Users/delvaw/Documents/MiceABC/sim.results.with.design.df.RData")
    # load(file = "/Users/delvaw/Documents/MiceABC/sim.results.with.design.df.RData")

    experim.median.features <- med(dplyr::select(sim.results.with.design.df, contains("y.")))$median

    # 3. Find intermediate features and RMSD.tol for which n.close.to.targets >= min.givetomice
    targets.diff <- targets.empirical - experim.median.features # First we determine how far the empirical targets are away from the median features of the executed experiments
    candidate.RMSD.tol <- Inf # Initially, we assume that the RMSD cut-off needs to be infinitely large to have sufficient observations to give to mice.

    n.steps.targets <- 100
    n.steps.RMSD.tol <- 100

    n.close.to.targets.mat <- matrix(NA, nrow = n.steps.targets+1, ncol = n.steps.RMSD.tol+1)
    final.RMSD.tol <- NA
    shift.fraction <- NA



    for (steps.intermediate.targets in 0:n.steps.targets){ # We make 10% steps from empirical targets to experimental median features
      candidate.intermediate.features <- targets.empirical - (targets.diff * steps.intermediate.targets) / n.steps.targets
      for (steps.RMSD.tol in 0:n.steps.RMSD.tol){ # We make 10% steps from RMSD.tol.max to RMSD.tol = 0
        RMSD.tol <- RMSD.tol.max * (n.steps.RMSD.tol - steps.RMSD.tol) / n.steps.RMSD.tol

        sum.sq.rel.dist <- rep(0, nrow(sim.results.with.design.df))
        for (i in 1:length(candidate.intermediate.features)) {
          name.dist <- paste0("y.", i, ".sq.rel.dist")
          value.dist <- ((sim.results.with.design.df[,i + x.offset] - candidate.intermediate.features[i]) / candidate.intermediate.features[i])^2
          assign(name.dist, value.dist)
          sum.sq.rel.dist <- sum.sq.rel.dist + get(name.dist)
        }
        RMSD <- sqrt(sum.sq.rel.dist / length(candidate.intermediate.features))
        n.close.to.targets <- sum(RMSD < RMSD.tol, na.rm = TRUE)
        n.close.to.targets.mat[(1+steps.intermediate.targets), (1+steps.RMSD.tol)] <- n.close.to.targets
        large.enough.training.df <- n.close.to.targets >= min.givetomice
        if (large.enough.training.df == TRUE){
          if (RMSD.tol < candidate.RMSD.tol){
            candidate.RMSD.tol <- RMSD.tol # Updating the candidate.RMSD.tol
            final.RMSD.tol <- RMSD.tol
            final.intermediate.features <- candidate.intermediate.features
            shift.fraction <- steps.intermediate.targets/n.steps.targets # 0 means targets.empirical; 1 means experim.median.features
            #calibration.list$intermediate.features[[wave]] <- final.intermediate.features
            #calibration.list$steps.intermediate.targets[[wave]] <- steps.intermediate.targets
            #calibration.list$RMSD.tol[[wave]] <- RMSD.tol
            #calibration.list$n.close.to.targets[[wave]] <- n.close.to.targets
          }
        }
      }
    }

    # n.close.to.targets.mat
    # final.RMSD.tol
    # final.intermediate.features
    # shift.fraction


    # 4. Calculate RMSD values for these intermediate.features and the RMSD.tol found
    sum.sq.rel.dist <- rep(0, nrow(sim.results.with.design.df))
    for (i in 1:length(final.intermediate.features)) {
      name.dist <- paste0("y.", i, ".sq.rel.dist")
      value.dist <- ((sim.results.with.design.df[,i + x.offset] - final.intermediate.features[i]) / final.intermediate.features[i])^2
      assign(name.dist, value.dist)
      sum.sq.rel.dist <- sum.sq.rel.dist + get(name.dist)
    }
    RMSD <- sqrt(sum.sq.rel.dist / length(final.intermediate.features))
    sim.results.with.design.df$RMSD <- RMSD
    n.close.to.targets <- sum(RMSD < final.RMSD.tol, na.rm = TRUE)

    # 5. Select n.close.to.targets shortest distances
    dist.order <- order(RMSD) # Ordering the squared distances from small to big.
    selected.distances <- dist.order[1:n.close.to.targets]
    sim.results.with.design.df.selected <- sim.results.with.design.df[selected.distances, ]
    # 5.b. Record highest RMSD value for that the selected experiments
    calibration.list$max.RMSD[[wave]] <- max(sim.results.with.design.df.selected$RMSD)

    # 6. Record selected experiments to give to mice for this wave
    calibration.list$selected.experiments[[wave]] <- sim.results.with.design.df.selected

    # 7. Put intermediate features in dataframe format
    final.intermediate.features.df <- as.data.frame(matrix(final.intermediate.features, ncol = length(final.intermediate.features)))
    names(final.intermediate.features.df) <- y.names

    # 8. Prepare dataframe to give to mice: selected experiments plus intermediate features
    df.give.to.mice <- dplyr::full_join(dplyr::select(sim.results.with.design.df.selected,
                                                      -contains("RMSD")), # adding target to training dataset
                                        final.intermediate.features.df,
                                        by = names(final.intermediate.features.df)) # "by" statement added to avoid printing message of the variables were used for joining

    # We are transforming parameters that are necessarily strictly positive: sigma, gamma.a, gamma.b.
    # We could also consider a similar transformation for input parameters that we think should be negative (e.g. formation.hazard.agegapry.gap_factor_man_exp) but for now not yet
    strict.positive.params <- c(4:8)
    df.give.to.mice[, strict.positive.params] <- log(df.give.to.mice[, strict.positive.params])

    # 9. Override default predictorMatrix with a sparser matrix
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
    # 9. pp.cp.6months.male (0.134)
    # 10. hiv.prev.lt25.women (0.2)
    # 11. hiv.prev.lt25.men (0.1)
    # 12. hiv.prev.25.34.women (0.42)
    # 13. hiv.prev.25.34.men (0.17)
    # 14. hiv.prev.35.44.women (0.37)
    # 15. hiv.prev.35.44.men (0.24)
    # 16. exp(growthrate)) (1.01)


    predictorMatrix[1:length(x.names), ] <- 0 # First we "empty" the relevant rows, then we refill them.
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

    print(c(nrow(df.give.to.mice), "nrows to give to mice"), quote = FALSE)

    # 10. Let mice propose parameter values that are predicted by the intermediate features
    mice.test <- mice(data = df.give.to.mice, # the dataframe with missing data
                      m = n.experiments, # number of imputations
                      predictorMatrix = predictorMatrix,
                      method = "norm",
                      defaultMethod = "norm",
                      maxit = maxit,
                      printFlag = FALSE,
                      data.init = NULL)

    # 11. Turn mice proposals into a new matrix of experiments
    experiments <- unlist(mice.test$imp) %>% matrix(., byrow = FALSE, ncol = length(x.names))
    # Before we check the suitability of the new experimental input parameter values, we must backtransform the log values to natural values
    experiments[, strict.positive.params] <- exp(experiments[, strict.positive.params])

    # 12. Update wave count
    wave <- wave + 1
  }

  # 13. Check if intermediate features have converged to empirical targests (if not, then maxwaves is reached)
  if (final.intermediate.features == targets.empirical){
    # 14. Do more "traditional mice optimisation" for the remaining waves
    while (wave <= maxwaves){
      sim.results.simple <- simpact.parallel(model = model,
                                             actual.input.matrix = experiments,
                                             seed_count = 0,
                                             n_cluster = 8)
      # save(sim.results.simple, file = "/Users/delvaw/Documents/MiceABC/sim.results.simple.RData")
      # load(file = "sim.results.simple.RData")

      new.sim.results.with.design.df <- as.data.frame(cbind(experiments,
                                                            sim.results.simple))
      x.names <- paste0("x.", seq(1:ncol(experiments)))
      y.names <- paste0("y.", seq(1:ncol(sim.results.simple)))
      x.offset <- length(x.names)
      names(new.sim.results.with.design.df) <- c(x.names, y.names)

      new.sim.results.with.design.df <- new.sim.results.with.design.df %>% dplyr::filter(complete.cases(.))

      if (is.null(sim.results.with.design.df)){ # TRUE for the first wave only
        sim.results.with.design.df <- rbind(sim.results.with.design.df,
                                            new.sim.results.with.design.df)
      } else {
        sim.results.with.design.df <- rbind(dplyr::select(sim.results.with.design.df,
                                                          -contains("RMSD")),
                                            new.sim.results.with.design.df)
      }


      # save(sim.results.with.design.df, file = "/Users/delvaw/Documents/MiceABC/sim.results.with.design.df.RData")
      # load(file = "/Users/delvaw/Documents/MiceABC/sim.results.with.design.df.RData")

      #experim.median.features <- med(dplyr::select(sim.results.with.design.df, contains("y.")))$median

      # 3. Find intermediate features and RMSD.tol for which n.close.to.targets >= min.givetomice
      # targets.diff <- targets.empirical - experim.median.features # First we determine how far the empirical targets are away from the median features of the executed experiments
      # candidate.RMSD.tol <- Inf # Initially, we assume that the RMSD cut-off needs to be infinitely large to have sufficient observations to give to mice.

      # n.steps.targets <- 100
      # n.steps.RMSD.tol <- 100

      # n.close.to.targets.mat <- matrix(NA, nrow = n.steps.targets+1, ncol = n.steps.RMSD.tol+1)
      # final.RMSD.tol <- NA
      # shift.fraction <- NA



      # for (steps.intermediate.targets in 0:n.steps.targets){ # We make 10% steps from empirical targets to experimental median features
      #   candidate.intermediate.features <- targets.empirical - (targets.diff * steps.intermediate.targets) / n.steps.targets
      #   for (steps.RMSD.tol in 0:n.steps.RMSD.tol){ # We make 10% steps from RMSD.tol.max to RMSD.tol = 0
      #     RMSD.tol <- RMSD.tol.max * (n.steps.RMSD.tol - steps.RMSD.tol) / n.steps.RMSD.tol
      #
      #     sum.sq.rel.dist <- rep(0, nrow(sim.results.with.design.df))
      #     for (i in 1:length(candidate.intermediate.features)) {
      #       name.dist <- paste0("y.", i, ".sq.rel.dist")
      #       value.dist <- ((sim.results.with.design.df[,i + x.offset] - candidate.intermediate.features[i]) / candidate.intermediate.features[i])^2
      #       assign(name.dist, value.dist)
      #       sum.sq.rel.dist <- sum.sq.rel.dist + get(name.dist)
      #     }
      #     RMSD <- sqrt(sum.sq.rel.dist / length(candidate.intermediate.features))
      #     n.close.to.targets <- sum(RMSD < RMSD.tol, na.rm = TRUE)
      #     n.close.to.targets.mat[(1+steps.intermediate.targets), (1+steps.RMSD.tol)] <- n.close.to.targets
      #     large.enough.training.df <- n.close.to.targets >= min.givetomice
      #     if (large.enough.training.df == TRUE){
      #       if (RMSD.tol < candidate.RMSD.tol){
      #         candidate.RMSD.tol <- RMSD.tol # Updating the candidate.RMSD.tol
      #         final.RMSD.tol <- RMSD.tol
      #         final.intermediate.features <- candidate.intermediate.features
      #         shift.fraction <- steps.intermediate.targets/n.steps.targets # 0 means targets.empirical; 1 means experim.median.features
      #         #calibration.list$intermediate.features[[wave]] <- final.intermediate.features
      #         #calibration.list$steps.intermediate.targets[[wave]] <- steps.intermediate.targets
      #         #calibration.list$RMSD.tol[[wave]] <- RMSD.tol
      #         #calibration.list$n.close.to.targets[[wave]] <- n.close.to.targets
      #       }
      #     }
      #   }
      # }

      # n.close.to.targets.mat
      # final.RMSD.tol
      # final.intermediate.features
      # shift.fraction

      # We can calculate the alpha fraction of close-enough data, based on the final.RMSD.tol from the previous wave, and use that as a starting point?


      # 4. Calculate RMSD values for these intermediate.features and the RMSD.tol found
      sum.sq.rel.dist <- rep(0, nrow(sim.results.with.design.df))
      for (i in 1:length(final.intermediate.features)) { # These are now equal to targets.empirical
        name.dist <- paste0("y.", i, ".sq.rel.dist")
        value.dist <- ((sim.results.with.design.df[,i + x.offset] - final.intermediate.features[i]) / final.intermediate.features[i])^2
        assign(name.dist, value.dist)
        sum.sq.rel.dist <- sum.sq.rel.dist + get(name.dist)
      }
      RMSD <- sqrt(sum.sq.rel.dist / length(final.intermediate.features))
      sim.results.with.design.df$RMSD <- RMSD
      n.close.to.targets <- sum(RMSD < final.RMSD.tol, na.rm = TRUE) # From hereon, we will keep using this same n.close.to.targets
      # alpha.fraction <- n.close.to.targets / nrow(sim.results.with.design.df) # Here we automatically set the alpha level.


      # 2b. Writing to list all the input and output of the executed experiments, so that we can plot it later
      # calibration.list$sim.results.with.design[[wave]] <- sim.results.with.design.df

      # 10. Calculate fraction of new (1-alpha frac *n.experiments) distances that are below "old" distance threshold. NOT CURRENTLY USED BECAUSE saturation.crit = 0.
      # below.old.treshold <- sim.results.with.design.df$sum.sq.rel.dist < rel.dist.cutoff
      # frac.below.old.threshold <- sum(below.old.treshold %in% TRUE) / round(n.experiments * (1-alpha))
      # if(frac.below.old.threshold < saturation.crit) saturation <- 1 # If less than the fraction saturation.crit of the new experiments a closer fit than the previous batch of retained experiments, the loop must be terminated.


      # 9. Merge with previously kept alpha fraction shortest distances
      # sim.results.with.design.df <- rbind(sim.results.with.design.df.selected, sim.results.with.design.df) # Initially sim.results.with.design.df.selected = NULL

      # 3. Keeping alpha fraction shortest distances


      # 5. Select n.close.to.targets shortest distances
      dist.order <- order(RMSD) # Ordering the squared distances from small to big.
      selected.distances <- dist.order[1:n.close.to.targets]
      sim.results.with.design.df.selected <- sim.results.with.design.df[selected.distances, ]
      # 5.b. Record highest RMSD value for that the selected experiments
      calibration.list$max.RMSD[[wave]] <- max(sim.results.with.design.df.selected$RMSD)

      # 6. Record selected experiments to give to mice for this wave
      calibration.list$selected.experiments[[wave]] <- sim.results.with.design.df.selected

      # 7. Put intermediate features in dataframe format
      final.intermediate.features.df <- as.data.frame(matrix(final.intermediate.features, ncol = length(final.intermediate.features)))
      names(final.intermediate.features.df) <- y.names

      # 8. Prepare dataframe to give to mice: selected experiments plus intermediate features
      df.give.to.mice <- dplyr::full_join(dplyr::select(sim.results.with.design.df.selected,
                                                        -contains("RMSD")), # adding target to training dataset
                                          final.intermediate.features.df,
                                          by = names(final.intermediate.features.df)) # "by" statement added to avoid printing message of the variables were used for joining

      # We are transforming parameters that are necessarily strictly positive: sigma, gamma.a, gamma.b.
      # We could also consider a similar transformation for input parameters that we think should be negative (e.g. formation.hazard.agegapry.gap_factor_man_exp) but for now not yet
      strict.positive.params <- c(4:8)
      df.give.to.mice[, strict.positive.params] <- log(df.give.to.mice[, strict.positive.params])

      # 9. Override default predictorMatrix with a sparser matrix
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
      # 9. pp.cp.6months.male (0.134)
      # 10. hiv.prev.lt25.women (0.2)
      # 11. hiv.prev.lt25.men (0.1)
      # 12. hiv.prev.25.34.women (0.42)
      # 13. hiv.prev.25.34.men (0.17)
      # 14. hiv.prev.35.44.women (0.37)
      # 15. hiv.prev.35.44.men (0.24)
      # 16. exp(growthrate)) (1.01)


      predictorMatrix[1:length(x.names), ] <- 0 # First we "empty" the relevant rows, then we refill them.
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

      print(c(nrow(df.give.to.mice), "nrows to give to mice"), quote = FALSE)

      # 10. Let mice propose parameter values that are predicted by the intermediate features
      mice.test <- mice(data = df.give.to.mice, # the dataframe with missing data
                        m = n.experiments, # number of imputations
                        predictorMatrix = predictorMatrix,
                        method = "norm",
                        defaultMethod = "norm",
                        maxit = maxit,
                        printFlag = FALSE,
                        data.init = NULL)

      # 11. Turn mice proposals into a new matrix of experiments
      experiments <- unlist(mice.test$imp) %>% matrix(., byrow = FALSE, ncol = length(x.names))
      # Before we check the suitability of the new experimental input parameter values, we must backtransform the log values to natural values
      experiments[, strict.positive.params] <- exp(experiments[, strict.positive.params])

      wave <- wave + 1
    }
  }

  # 15. Stop clock and return calibration list
  calibration.list$secondspassed <- proc.time() - ptm # Stop the clock
  return(calibration.list)
}
