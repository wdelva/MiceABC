mice.wrapper <- function(mice.input.list){
  library(mice)
  library(magrittr)
  mice.imputation <- mice(data = mice.input.list$data, # the dataframe with missing data
                          m = mice.input.list$m, # number of imputations
                          predictorMatrix = mice.input.list$predictorMatrix,
                          method = mice.input.list$method,
                          defaultMethod = mice.input.list$defaultMethod,
                          maxit = mice.input.list$maxit,
                          printFlag = mice.input.list$printFlag,
                          seed = mice.input.list$seed)
  mice.guess <- unlist(mice.imputation$imp) %>% matrix(., byrow = FALSE, ncol = length(unlist(mice.imputation$imp))) %>% as.numeric()
  return(mice.guess)
}
