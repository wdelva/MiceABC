library(parallel)
# Using all cores can slow down the computer
# significantly, I therefore try to leave one
# core alone in order to be able to do something
# else during the time the code runs


# creating a training dataset:
n.training <- 100
mean.x1 <- 100
mean.x2 <- 2
mean.x3 <- 3
mean.x4 <- 4

x1 <- rnorm(n = n.training, mean = mean.x1, sd = mean.x1/10)
x2 <- rnorm(n = n.training, mean = mean.x2, sd = mean.x2/10)
x3 <- rnorm(n = n.training, mean = mean.x3, sd = mean.x3/10)
x4 <- rnorm(n = n.training, mean = mean.x4, sd = mean.x4/10)

train.df <- simpact.parallel(model = dummy.wrapper,
                             actual.input.matrix = matrix(cbind(x1, x2, x3, x4), nrow = n.training),
                             #nb_simul = 16,
                             seed_count = 0,
                             n_cluster = 8) %>% cbind(., x1, x2, x3, x4) %>% as.data.frame() %>% dplyr::filter(complete.cases(.))

targets.row <- c(200, 4000, 24, 17200, NA, NA, NA, NA)
formice.df <- rbind(train.df, targets.row)

# The sequential experiment
ptm <- proc.time() # Start the clock
mice.sequential <- mice(formice.df, m = 200, printFlag = FALSE, maxit = 20)
mice.results.sequential <- unlist(mice.sequential$imp) %>% matrix(., byrow = FALSE, ncol = 4)
secondspassed <- proc.time() - ptm # Stop the clock
secondspassed


# The partial parallel experiment
ptm <- proc.time() # Start the clock
# cores_2_use <- detectCores()# - 1
cl <- makeCluster(8)# cores_2_use)
clusterSetRNGStream(cl, 9956)
clusterExport(cl, "formice.df")
no.output<-clusterEvalQ(cl, library(mice))
imp_pars <-
  parLapply(cl = cl, X = 1:cores_2_use, fun = function(no){
    mice(formice.df, m = 50, printFlag = FALSE, maxit = 20)
  })
stopCluster(cl)
imp_merged <- imp_pars[[1]]
for (n in 2:length(imp_pars)){
  imp_merged <-
    ibind(imp_merged,
          imp_pars[[n]])
}
mice.results.mostly.parallel <- unlist(imp_merged$imp) %>% matrix(., byrow = FALSE, ncol = 4)
secondspassed <- proc.time() - ptm # Stop the clock
secondspassed


# The full parallel experiment
ptm <- proc.time() # Start the clock
library(foreach)
library(doParallel)
cl <- makeCluster(cores_2_use)
clusterSetRNGStream(cl, 9956)
registerDoParallel(cl)

library(mice)
imp_merged <-
  foreach(no = 1:cores_2_use,
          .combine = ibind,
          .export = "formice.df",
          .packages = "mice") %dopar%
          {
            mice(formice.df, m = 50, printFlag = FALSE, maxit = 20)
          }
stopCluster(cl)
mice.results.full.parallel <- unlist(imp_merged$imp) %>% matrix(., byrow = FALSE, ncol = 4)
secondspassed <- proc.time() - ptm # Stop the clock
secondspassed



