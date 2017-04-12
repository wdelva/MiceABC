# a function to handle a simulation
# gets a set of parameters
# returns results of evaluation criteria
# NOTE: It runs repeated simulations for stochastic models.
#       To control stochasticity it runs replicated simulations for current parameter combination
#       and calculates the mean simulation output.
#       If your model is deterministic, just set no.repeated.sim to 1.
agemix.simulate <- function(param.set, parameter.names = c("age.gap.tol.intercept", "age.gap.tol.coef"), no.repeated.sim=10, nl.obj="agemix.nl1") {#, trace.progress, iter.length, function.name) {
  # some security checks
  if (length(param.set) != length(parameter.names))
  { stop("Wrong length of param.set!") }
  if (no.repeated.sim <= 0)
  { stop("Number of repetitions must be > 0!") }
  if (length(parameter.names) <= 0)
  { stop("Length of parameter.names must be > 0!") }
  
  # an empty list to save the simulation results
  eval.values <- NULL
  
  # repeated simulations (to control stochasticity)
  for (i in 1:no.repeated.sim)
  {
    # create a random-seed for NetLogo from R, based on min/max of NetLogo's random seed
    # for NetLogo 4:
    #NLCommand("random-seed",runif(1,-9007199254740992,9007199254740992), nl.obj=nl.obj)
    # since NetLogo 5:
    NLCommand("random-seed",runif(1,-2147483648,2147483647), nl.obj=nl.obj)
    
    # TODO: adapt the following to your simulation model
    # This is the stuff for one simulation
    NLCommand("setup", nl.obj=nl.obj)
    
    # set NetLogo parameters to current parameter values
    lapply(seq(1:length(parameter.names)), function(x) {NLCommand("set ",parameter.names[x], param.set[x], nl.obj=nl.obj)})
    
    # run simulation
    # TODO: adapt to your simulation process
    # warm-up (1299)
    NLDoCommand(99, "go", nl.obj=nl.obj)
    
    cal.crit <- NLDoReport(1,
                           "go",
                           c("ticks",
                             "cumsum.transmissions",
                             "mean.age.infecteds",
                             "var.age.infecteds"),
                           as.data.frame = TRUE,
                           df.col.names = c("timestep",
                                            "cumsum.transmissions",
                                            "mean.age.hivpos",
                                            "var.age.hivpos"),
                           nl.obj=nl.obj)
    
    # append to former results
    eval.values <- rbind(eval.values, cal.crit) #calibration.criteria)
  }
  
  
  # return the mean of the repeated simulation results
  if (no.repeated.sim > 1)
  {
    return(colMeans(eval.values, na.rm = TRUE))
  }
  else {
  return(colMeans(eval.values, na.rm = TRUE))
  }
}
