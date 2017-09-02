# mice.model (mice.wrapper) is the function that takes as arguments a dataframe to be given to mice,
# as well as the arguments for the mice function:m = n.experiments (1), # number of imputations
# predictorMatrix = predictorMatrix
# method = "norm"
# defaultMethod = "norm"
# maxit = maxit
# printFlag = FALSE
# The output of mice.model (mice.guess) is vector of proposed values for the missing input parameters of the IBM:
# experiments <- unlist(mice.imputation$imp)



mice.parallel <- function(mice.model = mice.wrapper,
                          df.give.to.mice = df.give.to.mice,
                          m = 1,
                          predictorMatrix = predictorMatrix,
                          method = "norm",
                          defaultMethod = "norm",
                          maxit = 5,
                          printFlag = FALSE,
                          seed_count = 0,
                          n_cluster = 8,
                          nb_simul = n.experiments){
  cl <- makeCluster(getOption("cl.cores", n_cluster))


  mice.input.lists <- vector("list", length = nb_simul) #list(NULL)

  mice.input.list <- list() # An element of the mice.input.lists object. Used as input for mice
  mice.input.list$data <- df.give.to.mice # Not dependent on i
  mice.input.list$m <- m # number of imputations
  mice.input.list$predictorMatrix <- predictorMatrix
  mice.input.list$method <- method
  mice.input.list$defaultMethod <- defaultMethod
  mice.input.list$maxit <- maxit
  mice.input.list$printFlag <- printFlag

  tab_simul_guess = NULL

  # nb_simul <- nrow(actual.input.matrix)

  for (i in 1:nb_simul) {
    mice.input.list$seed <- seed_count + i
    mice.input.lists[[i]] <- mice.input.list
  }
  list_simul_guesses = parLapplyLB(cl, mice.input.lists,
                                       mice.model)
  tab_simul_guess <- do.call(rbind, list_simul_guesses)
  stopCluster(cl)
  return(tab_simul_guess)
}

