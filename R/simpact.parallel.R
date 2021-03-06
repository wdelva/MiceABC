simpact.parallel <- function(model = simpact.wrapper,
                             actual.input.matrix = matrix(rep(c(1:15), 16), nrow = 16),
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
