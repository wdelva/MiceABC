VEME.wrapper2 <- function(inputvector = input.vector){ # This function does not produce the incidence and does not have an input for the model scenario
  #source("/user/data/gent/vsc400/vsc40070/phylo/scripts/00-Functions.R")
  source("/Users/delvaw/Documents/MiceABC/R/00-Functions.R")





  cfg.list["hivtransmission.param.f1"] = log(inputvector[2])
  cfg.list["hivtransmission.param.f2"] = log(log(sqrt(inputvector[2])) / log(inputvector[2])) / 5
  cfg.list["formation.hazard.agegapry.gap_agescale_man"] = inputvector[3]
  cfg.list["formation.hazard.agegapry.gap_agescale_woman"] = inputvector[3]
  cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[4]
  cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[4]
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[5]
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[5]
  cfg.list["person.eagerness.man.dist.gamma.a"] <- inputvector[6]
  cfg.list["person.eagerness.woman.dist.gamma.a"] <- inputvector[7]
  cfg.list["person.eagerness.man.dist.gamma.b"] <- inputvector[8]
  cfg.list["person.eagerness.woman.dist.gamma.b"] <- inputvector[9]




  cfg <- cfg.list

  cfg["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3
  cfg["monitoring.fraction.log_viralload"] <- 0.3
  cfg["person.vsp.toacute.x"] <- 5 # See Bellan PLoS Medicine

  seedid <- inputvector[1]
  #cfg["person.agegap.man.dist.fixed.value"] <- -2 # inputvector[2]
  #cfg["person.agegap.woman.dist.fixed.value"] <- -2 # inputvector[2]
  cfg["formation.hazard.agegapry.gap_factor_man_exp"] <- inputvector[10] ######### -0.5
  cfg["formation.hazard.agegapry.gap_factor_woman_exp"] <- inputvector[10] ######### -0.5
  cfg["formation.hazard.agegapry.baseline"] <- inputvector[11]

  cfg["formation.hazard.agegapry.numrel_man"] <- inputvector[12]
  cfg["formation.hazard.agegapry.numrel_woman"] <- inputvector[13]
  # inputvector[14] is conception.alpha.base (higher up)
  cfg["dissolution.alpha_0"] <- inputvector[15]
  cfg["dissolution.alpha_4"] <- inputvector[16]



  results <- tryCatch(simpact.run(configParams = cfg,
                                  destDir = destDir,
                                  agedist = age.distr,
                                  seed = seedid, #, Introducing ART has helped to keep the prevalence high
                                  intervention = intervention),
                      error = simpact.errFunction)

  if (length(results) == 0){
    outputvector <- rep(NA, 18)
  } else {
    if (as.numeric(results["eventsexecuted"]) >= (as.numeric(cfg["population.maxevents"]) - 1)) {
      outputvector <- rep(NA, 18)
    } else {
      datalist.agemix <- readthedata(results)
      agemix.df <- agemix.df.maker(datalist.agemix)

      agemix.model <- pattern.modeller(dataframe = agemix.df,
                                       agegroup = c(15, 50),
                                       timepoint = datalist.agemix$itable$population.simtime[1],
                                       timewindow = 1)#3)

      # men.lme <- tryCatch(agemixing.lme.fitter(data = dplyr::filter(agemix.model[[1]], Gender =="male")),
      #                     error = agemixing.lme.errFunction) # Returns an empty list if the lme model can't be fitted

      men.lmer <- tryCatch(ampmodel(data = dplyr::filter(agemix.model[[1]], Gender =="male")),
                           error = agemixing.lme.errFunction) # Returns an empty list if the lme model can't be fitted

      bignumber <- NA # let's try if NA works (instead of 9999 for example)

      AAD.male <- ifelse(length(men.lmer) > 0, mean(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
      SDAD.male <- ifelse(length(men.lmer) > 0, sd(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
      #powerm <- ifelse(length(men.lme) > 0, as.numeric(attributes(men.lme$apVar)$Pars["varStruct.power"]), bignumber)
      slope.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$coefficients[2, 1], bignumber) #summary(men.lmer)$tTable[2, 1], bignumber)
      WSD.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$sigma, bignumber) #WVAD.base <- ifelse(length(men.lme) > 0, men.lme$sigma^2, bignumber)

      BSD.male <- ifelse(length(men.lmer) > 0, bvar(men.lmer), bignumber) # Bad name for the function because it actually extracts between subject standard deviation # BVAD <- ifelse(length(men.lmer) > 0, getVarCov(men.lme)[1,1], bignumber)

      intercept.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$coefficients[1,1] - 15, bignumber)

      #agegap.mean <- mean(datalist.agemix$rtable$AgeGap)
      #agegap.sd <- sd(datalist.agemix$rtable$AgeGap)


      ###
      # Now we fit the negative binomial distribution
      ###
      degrees.df <- degree.df.maker(dataframe.df = agemix.df,
                                    agegroup = c(18, 50),
                                    hivstatus = 2,
                                    survey.time = datalist.agemix$itable$population.simtime[1],
                                    window.width = 1,
                                    gender.degree = "male",
                                    only.new = TRUE)


      # If we want to know % of men who had 0 partners in the past 12 months,
      # We need to compare nrow(degree.df) with the number of 18-50 year old men #HIV negative women
      # that were alive at the time of the survey
      allsurveymen <- dplyr::filter(datalist.agemix$ptable,
                                    Gender == 1, # Male
                                    TOD > datalist.agemix$itable$population.simtime[1], # Still alive at the time of the survey
                                    TOB <= datalist.agemix$itable$population.simtime[1] - 18, # Not younger than 18 at the time of the survey
                                    TOB > datalist.agemix$itable$population.simtime[1] - 50)#, # Not older than 50 at the time of the survey
      #InfectTime > datalist.agemix$itable$population.simtime[1]) # HIV negative at the time of the survey

      # Creating the vector of degrees
      degree.vector <- c(rep(0, times = (nrow(allsurveymen) - nrow(degrees.df))), degrees.df$Degree)
      meandegree.male <- mean(degree.vector)
      # hist(degree.vector, 10)
      # degree.vector <- rnbinom(n = 100, size =1.28, mu = 0.66) # PLACEHOLDER FOR NOW
      # Fitting the negative binomial distribution to this vector

      fit.negbin <- tryCatch(fitdist(degree.vector, "nbinom"), error = agemixing.lme.errFunction)
      shape.nb.male <- ifelse(length(fit.negbin) > 0, as.numeric(fit.negbin$estimate[2]), bignumber)
      scale.nb.male <- ifelse(length(fit.negbin) > 0, as.numeric(fit.negbin$estimate[1]), bignumber) #(theta = p/(1-p))


      # Concurrency point prevalence 6 months before a survey, among men
      pp.cp.6months.male <- concurr.pointprev.calculator(datalist = datalist.agemix,
                                                         timepoint = datalist$itable$population.simtime[1] - 0.5)


      hiv.prev.lt25.women <- prevalence.calculator(datalist = datalist.agemix,
                                                   agegroup = c(15, 25),
                                                   timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[2]
      hiv.prev.lt25.men <- prevalence.calculator(datalist = datalist.agemix,
                                                 agegroup = c(15, 25),
                                                 timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[1]
      hiv.prev.25.34.women <- prevalence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(25, 35),
                                                    timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[2]
      hiv.prev.25.34.men <- prevalence.calculator(datalist = datalist.agemix,
                                                  agegroup = c(25, 35),
                                                  timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[1]
      hiv.prev.35.44.women <- prevalence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(35, 45),
                                                    timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[2]
      hiv.prev.35.44.men <- prevalence.calculator(datalist = datalist.agemix,
                                                  agegroup = c(35, 45),
                                                  timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[1]




      growthrate <- pop.growth.calculator(datalist = datalist.agemix,
                                          timewindow = c(0, datalist.agemix$itable$population.simtime[1]))


      #######
      # And now we add HIV incidence at two selected time points
      #######
      # inc.times <- seq(from = datalist.agemix$itable$population.simtime[1] - 5,
      #                  by = 5,
      #                  to = datalist.agemix$itable$population.simtime[1])
      # inc.vector <- rep(NA, length(inc.times))
      # for (inc.time in (1:length(inc.times))){
      #
      #   inc.tibble <- incidence.calculator(datalist = datalist.agemix,
      #                                                agegroup = c(0, 500), # Essentially incidence in the entire population
      #                                                timewindow = c(inc.times[inc.time] - 1, inc.times[inc.time]),
      #                                                only.active = "No")
      #   inc.vector[inc.time] <- inc.tibble$incidence[3]
      # }
      # inc.ratio <- inc.vector[2] / inc.vector[1]
      #
      # inc.end <- inc.vector[2]
      inc.time <- datalist.agemix$itable$population.simtime[1]
      inc.tibble <- incidence.calculator(datalist = datalist.agemix,
                                         agegroup = c(0, 500), # Essentially incidence in the entire population
                                         timewindow = c(inc.time - 1, inc.time),
                                         only.active = "No")
      inc.end <- inc.tibble$incidence[3]


      #######
      # And transmission tree
      #######

      recent.vs.old.nodetimes.ratio <- NA # Will be overwritten if it can be computed

      transm.ls <- transmNetworkBuilder.baseline(datalist = datalist.agemix,
                                                 endpoint = datalist.agemix$itable$population.simtime[1])

      transm.ls.attrib <- attributes.trans.network(datalist = datalist.agemix,
                                                   endpoint = datalist.agemix$itable$population.simtime[1])

      ID.vect.alive <- alive.infected(datalist = datalist.agemix,
                                      timepoint = datalist.agemix$itable$population.simtime[1]) %>% filter(., Infected == "TRUE") %>% dplyr::select(., ID) %>% unlist() %>% as.numeric() #datalist.agemix$itable$population.simtime[1])  #### REPLACE 15 with datalist.agemix$itable$population.simtime[1])   !!!!




      internal.node.times <- numeric() # Initiation, to append more branching times, for all seeds

      for (transm.ls.element in 1:length(transm.ls)){
        transnetwork <- transm.ls[[transm.ls.element]] # 1 seed considered at a time
        transnetw.attr <- transm.ls.attrib[[transm.ls.element]]

        if (length(ID.vect.alive) > 1){

          alive.selected.vector <- transnetw.attr$id.orig.recip %in% ID.vect.alive
          sequences.keep <- transnetw.attr$id[alive.selected.vector] # Only the sequences of these tip labels (i.e. new IDs) can be used to build the phylo tree
          numbered.tips <- 1:length(alive.selected.vector)


          if (length(transnetwork$parent) > 1){ # Let's start with only small networks
            epi.tree <- epi2tree2(transnetwork)
            tree.dat.full <- simSeq(epi.tree,
                                    l = 100,
                                    bf = freq,
                                    #rootseq = hiv.seq.env.short,
                                    type = "DNA",
                                    rate = overall.rate)
            #
            #     #sim <- sequence.simulation(transtree = epi.tree,
            #     #                           seedSeq = hiv.seq.env.short,
            #     #                           base.freq = freq)
            #
            #     #saveAlignment.PhyloSim(sim,
            #     #                       file = paste("/Users/delvaw/Documents/MiceABC/HIVSeq_fullNetwork.fasta",
            #     #                                    sep = ""),
            #     #                       skip.internal = TRUE,
            #     #                       paranoid = TRUE)
            #     # to handle sequences, let use phangorn and ape packages and read the sequence data and build the phylogenetic tree
            #     #seq.sim.size_full <- read.FASTA("/Users/delvaw/Documents/MiceABC/HIVSeq_fullNetwork.fasta")
            #     #tree.dat.full <- phyDat(seq.sim.size_full, type = "DNA")
            #
            numbered.keep <- numbered.tips[as.numeric(names(tree.dat.full)) %in% sequences.keep]
            #
            tree.dat.pruned <- subset(tree.dat.full,
                                      subset = numbered.keep)# Only keep sequences of people that were alive at datalist.agemix$itable$population.simtime[1]
            if (length(tree.dat.pruned) > 1){
              tree.ml.pruned <- dist.ml(tree.dat.pruned, model = "F81", bf = freq) # F81
              tree.sim.pruned <- upgma(tree.ml.pruned)

              internal.node.times.element <- branching.times(tree.sim.pruned) # datalist.agemix$itable$hivseed.time[1] + branching.times(tree.sim.full) * (datalist.agemix$itable$population.simtime[1] - datalist.agemix$itable$hivseed.time[1])
              internal.node.times.unsorted <- c(internal.node.times, internal.node.times.element)
              internal.node.times.sorted <- sort(as.numeric(internal.node.times.unsorted))
              recent.transmissions <- sum(internal.node.times.sorted < 1 * overall.rate)
              #less.recent.transmissions <- sum(internal.node.times.sorted >= 2 * overall.rate & internal.node.times.sorted < 4 * overall.rate)
              #recent.vs.old.nodetimes.ratio <- recent.transmissions / less.recent.transmissions
              # This can be a real number >=0, NaN (0/0) or Inf (num/0)
            }
          }
        }
      }






      outputvector <- c(AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male,
                        shape.nb.male, scale.nb.male,
                        #meandegree.male,
                        pp.cp.6months.male,
                        hiv.prev.lt25.women,
                        hiv.prev.lt25.men,
                        hiv.prev.25.34.women,
                        hiv.prev.25.34.men,
                        hiv.prev.35.44.women,
                        hiv.prev.35.44.men,
                        exp(growthrate),
                        inc.end,
                        recent.transmissions#recent.vs.old.nodetimes.ratio
                        )
    }
  }
  return(outputvector)#,
           #   transm.ls = transm.ls,
          #    tree.sim.full = tree.sim.full,
          #    internal.node.times.sorted = internal.node.times.sorted)) #return(list(outputvector, transm.ls))
}
