#' Adaptive Population ABC with MICE-based supervised learning
#'
#' This function implements an adaptive population ABC scheme with supervised learning through Multivariate Imputation by Chained Equations (MICE) to produce new particles (i.e. samples from the parameter space and their associated output statistics)
#'
#' @param targets The vector of target statistics to calibrate the simulation model against
#' @param n.experiments The number of initial multivariate samples from the parameter space
#' @param lls Lower limits of the parameter ranges to be explored
#' @param uls Upper limits of the parameter ranges to be explored
#' @param model Handle to the wrapper function that runs the simulation model
#' @param maxiter Maximum number of iterations (waves of simulations) before the algorithm is forced to stop
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
#'                             maxiter = 10,
#'                             alpha = 0.5,
#'                             rel.dist.tol = 0.1)
#' @import boot
#' @import randtoolbox
#' @import dplyr

MiceABC <- function(targets = c(logit(0.15), log(0.015)),
                    n.experiments = 100,
                    lls = c(0, 0),
                    uls = c(10, 1),
                    model = model.wrapper,
                    maxiter = 10,
                    alpha = 0.5,
                    rel.dist.tol = 0.1){
  range.width <- uls - lls
  ll.mat <- matrix(rep(lls, n.experiments), nrow = n.experiments, byrow = TRUE)
  range.width.mat <- matrix(rep(range.width, n.experiments), nrow = n.experiments, byrow = TRUE)
  sobol.seq.0.1 <- sobol(n = n.experiments, dim = length(lls), init = TRUE, scrambling = 1, seed = 1, normal = FALSE)
  experiments <- ll.mat + sobol.seq.0.1 * range.width.mat

  calibration.list <- list() # initiating the list where all the output of MiceABC will be stored

  iteration <- 1 # initiating the loop of waves of simulations (one iteration is one wave)

  while (iteration <= maxiter){
    print(c("iteration", iteration), quote = FALSE) # only for local testing. To be deleted after testing is completed
    large.enough <- t(t(experiments) >= lls)
    small.enough <- t(t(experiments) <= uls)
    fine <- large.enough & small.enough
    fine.flag <- as.logical(rowSums(fine) == ncol(fine))
    experiments <- experiments[fine.flag, ]
    try(if(nrow(experiments) == 0) stop("no experiments in the range"))
    print(c(nrow(experiments), "experiments to be run"), quote = FALSE) # only for local testing. To be deleted after testing is completed

    # Running the experiments
    nl.path <- "/Applications/NetLogo 6.0/app/"
    NLStart(nl.path, gui=FALSE, nl.obj="agemix.nl1", nl.jarname = 'netlogo-6.0.0.jar')
    model.path <- "/Users/delvaw/Google Drive/IBM course/tutorials/AgeMixing/Calibration/dummy.0.1.nlogo"
    NLLoadModel(model.path, nl.obj="agemix.nl1")
    source(file = "/Users/delvaw/Google Drive/IBM course/tutorials/AgeMixing/Calibration/agemix.simulate.R")
    sim.results.simple <- apply(experiments, 1, agemix.simulate) # Maybe this apply() can be replaced by a faster foreach loop?
    if (ncol(sim.results.simple) >= 2) {sim.results.simple <- t(sim.results.simple)}
    # numCores <- detectCores()
    # sobol.experiments.list <- split(c[1:10, ], seq(nrow(sobol.experiments[1:10, ]))) # Turning the dataframe into a list
    # sim.results.paral <- mclapply(sobol.experiments.list, agemix.simulate, mc.cores = numCores)
    lapply.limit <- ifelse(is.vector(sim.results.simple), length(sim.results.simple), nrow(sim.results.simple))

    #sim.results.with.design <- lapply(1:lapply.limit, function(x) {experiment <- x;
    age.gap.tol.intercept <- experiments[, 1] # guess.1.df$age.gap.tol.intercept[x];
    age.gap.tol.coef <- experiments[, 2] #guess.1.df$age.gap.tol.coef[x];
    #cbind(experiment, sim.results.simple[x, ], age.gap.tol.intercept, age.gap.tol.coef)})
    # First we combine all the elements of the list into one long dataset
    #sim.results.with.design.df <- as.data.frame(do.call(rbind, sim.results.with.design))

    experiment <- 1:nrow(experiments)

    sim.results.with.design.df <- as.data.frame(cbind(sim.results.simple, experiment, age.gap.tol.intercept, age.gap.tol.coef))

    calibration.list$inputoutput[[iteration]] <- sim.results.with.design.df
    calibration.list$intercept[[iteration]] <- summary(experiments[ , 1])
    calibration.list$coef[[iteration]] <- summary(experiments[, 2])

    # Now we combine this new set of simulation inputs and outputs with the target statistics
    df <- dplyr::select(full_join(sim.results.with.design.df, data.frame(cumsum.transmissions = targets[1],
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
    dist.order <- order(df$sum.sq.rel.dist) # Ordering the squared distances from small to big. The last observation (targets) should be ranked first
    shortest.dist.percentile <- alpha # 0.25
    last.one.selected <- round(shortest.dist.percentile * length(dist.order))
    selected.distances <- dist.order[1:last.one.selected]
    df.selected <- df[selected.distances, ]
    df.filtered <- filter(df,
                          abs(df$cumsum.transmissions - targets[1]) / targets[1] <= rel.tol,
                          abs(df$mean.age.hivpos - targets[2]) / targets[2] <= rel.tol,
                          abs(df$var.age.hivpos - targets[3]) / targets[3] <= rel.tol)

    df.give.to.mice <- df.filtered
    if (as.numeric(last.one.selected) > nrow(df.filtered)){df.give.to.mice <- df.selected}

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

    mice.test <- mice(data = df.give.to.mice, # the dataframe with missing data
                      m = n.experiments, # number of imputations
                      method = "norm",
                      defaultMethod = "norm",
                      maxit = 20,
                      printFlag = FALSE,
                      data.init = NULL)

    experiments <- cbind(unlist(mice.test$imp$age.gap.tol.intercept),
                         unlist(mice.test$imp$age.gap.tol.coef))
    iteration <- iteration + 1

  }

  proc.time() - ptm

  df <- do.call(rbind, guesses$inputoutput[1:4])

  rel.tol <- 0.002 # 0.25

  # We also calculate sum of squared relative distances
  df$cumsum.transmissions.sq.rel.dist <- ((df$cumsum.transmissions - targets[1]) / targets[1])^2
  df$mean.age.hivpos.sq.rel.dist <- ((df$mean.age.hivpos - targets[2]) / targets[2])^2
  df$var.age.hivpos.sq.rel.dist <- ((df$var.age.hivpos - targets[3]) / targets[3])^2
  df$sum.sq.rel.dist <- df$cumsum.transmissions.sq.rel.dist + df$mean.age.hivpos.sq.rel.dist +df$var.age.hivpos.sq.rel.dist
  dist.order <- order(df$sum.sq.rel.dist) # Ordering the squared distances from small to big. The last observation (targets) should be ranked first
  shortest.dist.percentile <- 0.1344086 #50/372 # 0.25
  last.one.selected <- round(shortest.dist.percentile * length(dist.order))
  selected.distances <- dist.order[1:last.one.selected]
  df.selected <- df[selected.distances, ]
  df.filtered <- filter(df,
                        abs(df$cumsum.transmissions - targets[1]) / targets[1] <= rel.tol,
                        abs(df$mean.age.hivpos - targets[2]) / targets[2] <= rel.tol,
                        abs(df$var.age.hivpos - targets[3]) / targets[3] <= rel.tol)

  df.give.to.mice <- df.filtered
  if (as.numeric(last.one.selected) > nrow(df.filtered)){df.give.to.mice <- df.selected}

  summary(df.give.to.mice)
  pairs(df.give.to.mice[,6:7])


}


