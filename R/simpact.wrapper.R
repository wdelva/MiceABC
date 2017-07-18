simpact.wrapper <- function(inputvector = input.vector){ # This function does not produce the incidence and does not have an input for the model scenario
  library(RSimpactHelper)
  library(RSimpactCyan)
  library(magrittr)
  library(dplyr)
  library(nlme)
  library(lme4)
  library(fitdistrplus)
  library(lmtest)
  library(tidyr)
  library(data.table)

  agemixing.lme.errFunction <- function(e)
  {
    return(list())
  }

  simpact.errFunction <-function(e){
    if (length(grep("NaN",e$message)) != 0){
      return(list())
    }
  }


  ampmodel <- function(data = dplyr::filter(agemix.model[[1]], Gender =="male")) {
    lmer(pagerelform ~ agerelform0 + (1 | ID),
         data = data,
         REML = TRUE)
  }



  bvar <- function(model) {

    # Outputs a df with between-subject variance, upr & lwr limits

    # Must take an merMod  object

    bsd <- as.numeric(as.data.frame(VarCorr(model))[1,5])


  }



  degree.df.maker <- function (dataframe.df, agegroup = c(15, 30), hivstatus = 0,
                               survey.time = 40, window.width = 1, gender.degree = "female",
                               only.new = TRUE)
  {
    dataframe.rels.df <- dataframe.df
    dataframe.rels.df <- dplyr::select(dataframe.rels.df, -episodeorder, -FormTime, -DisTime)
    dataframe.rels.df <- unique.data.frame(dataframe.rels.df)
    rels.form.dis.df <- dplyr::summarise(dplyr::group_by(dataframe.df,
                                                         relid), FormTime = min(FormTime), DisTime = max(DisTime))
    dfnew <- dplyr::left_join(x = dataframe.rels.df, y = rels.form.dis.df,
                              by = "relid")
    dfnew <- dplyr::filter(dfnew, TOD > survey.time)
    {
      if (hivstatus == 0) {
        dfnew <- dplyr::filter(dfnew, InfectTime > survey.time)
      }
      else if (hivstatus == 1) {
        dfnew <- dplyr::filter(dfnew, InfectTime <= survey.time)
      }
    }
    {
      if (only.new) {
        dfnew <- dplyr::filter(dfnew, FormTime >= survey.time -
                                 window.width, FormTime < survey.time, DisTime >
                                 survey.time - window.width, Gender == gender.degree,
                               survey.time - TOB >= agegroup[1], survey.time -
                                 TOB < agegroup[2])
      }
      else {
        dfnew <- dplyr::filter(dfnew, FormTime < survey.time,
                               DisTime > survey.time - window.width, Gender ==
                                 gender.degree, survey.time - TOB >= agegroup[1],
                               survey.time - TOB < agegroup[2])
      }
    }
    uniqueOut <- dfnew %>% dplyr::select(ID, relid) %>% distinct %>%
      rename(Degree = relid)
    if (dim(uniqueOut)[1] != 0) {
      degreedata.df <- aggregate(Degree ~ ID, data = uniqueOut,
                                 length)
    } else {
      degreedata.df <- data.frame(ID = NA, Degree = NA)
    }
    return(degreedata.df)
  }





  concurr.pointprev.calculator <- function(datalist,
                                           timepoint = datalist$itable$population.simtime[1] - 0.5){

    output <- data.table()
    DTalive.infected <- alive.infected(datalist = datalist,
                                       timepoint = timepoint, site = "All") # First we only take the data of people who were alive at time_i
    agemix.df <- agemix.df.maker(datalist)
    degrees.df <- degree.df.maker(dataframe.df = agemix.df,
                                  agegroup = c(15, 50),
                                  hivstatus = 2,
                                  survey.time = timepoint,
                                  window.width = 0,
                                  gender.degree = "male",
                                  only.new = FALSE)

    number.people.with.cps <- sum(degrees.df$Degree > 1)
    popsize <- nrow(DTalive.infected)
    concurr.pointprevalence <- number.people.with.cps / popsize

    return(concurr.pointprevalence)
  }




  cd4.atARTinit <- function(datalist = datalist,
                            agegroup = c(15, 30),
                            timewindow = c(15, 30),
                            cd4count=350, site="All"){

    cd4.atARTinit <- age.group.time.window(datalist = datalist,
                                           agegroup = agegroup,
                                           timewindow = timewindow, site="All")

    cd4.atARTinit <- subset(cd4.atARTinit, TreatTime !=Inf) #HIV positive individuals

    raw.df <- data.frame(cd4.atARTinit)
    art.df <- subset(datalist$ttable, ID %in% cd4.atARTinit$ID &
                       TStart > timewindow[1] & TStart < timewindow[2])

    ##What if the person dropped out and come back again?
    art.df <- data.frame(dplyr::summarise(dplyr::group_by(art.df, ID, Gender),
                                          CD4atARTstart = min(CD4atARTstart)))

    #indicate those who started their treatment when their CD4 count was below a given threshold
    art.df <- art.df %>% dplyr::mutate(ART.start.CD4 = CD4atARTstart < cd4count)

    # Now we apply the left_join dplyr function to add the ART status to raw.df.
    raw.df <- dplyr::left_join(x = raw.df, y = art.df, by = c("ID", "Gender"))

    #provide a summary of those that are on treatment and those that started below a threshold
    cd4count.atARTInit <- data.frame(dplyr::summarise(dplyr::group_by(raw.df, Gender),
                                                      TotalCases = n(),
                                                      LessCD4initThreshold =sum(ART.start.CD4)))

    return(cd4count.atARTInit)
  }
  root.path <- dirname(getwd())
  destDir <- paste0(root.path, "/temp")
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2, #0.21,#1,
                                   population.simtime = 35, #40,
                                   population.nummen = 2500,
                                   population.numwomen = 2500,
                                   hivseed.time = 10,
                                   hivseed.type = "amount",
                                   hivseed.amount = 25, #30,
                                   hivseed.age.min = 20,
                                   hivseed.age.max = 50,
                                   hivtransmission.param.a = -1,
                                   hivtransmission.param.b = -90,
                                   hivtransmission.param.c = 0.5,
                                   hivtransmission.param.f1 = log(inputvector[2]) , #log(2),
                                   hivtransmission.param.f2 = log(log(sqrt(inputvector[2])) / log(inputvector[2])) / 5, #log(log(1.4) / log(2)) / 5,
                                   formation.hazard.agegapry.gap_factor_man_age = -0.01, #-0.01472653928518528523251061,
                                   formation.hazard.agegapry.gap_factor_woman_age = -0.01, #-0.0726539285185285232510561,
                                   formation.hazard.agegapry.meanage = -0.025,
                                   formation.hazard.agegapry.gap_factor_man_const = 0,
                                   formation.hazard.agegapry.gap_factor_woman_const = 0,
                                   formation.hazard.agegapry.gap_factor_man_exp = -1, #-6,#-1.5,
                                   formation.hazard.agegapry.gap_factor_woman_exp = -1, #-6,#-1.5,
                                   formation.hazard.agegapry.gap_agescale_man = inputvector[3], # 0.25,
                                   formation.hazard.agegapry.gap_agescale_woman = inputvector[3], # 0.25,#-0.30000007,#-0.03,
                                   debut.debutage = 15,
                                   conception.alpha_base = inputvector[14]#-2.5#,
                                   #person.art.accept.threshold.dist.fixed.value = 0
  )


  cfg.list["formation.hazard.agegapry.baseline"] <- 2
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.fraction.log_viralload"] <- 0.3
  cfg.list["dropout.interval.dist.uniform.min"] <- 100
  cfg.list["dropout.interval.dist.uniform.max"] <- 200


  cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
  #cfg.list["person.agegap.man.dist.fixed.value"] <- -6
  cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
  #cfg.list["person.agegap.woman.dist.fixed.value"] <- -6

  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.cd4.threshold"] <- 0
  cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[4] # 0
  cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[4] # 0
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[5] ######### 3
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[5] ######### 3
  cfg.list["person.eagerness.man.dist.gamma.a"] <- inputvector[6] ######### 0.23
  cfg.list["person.eagerness.woman.dist.gamma.a"] <- inputvector[7] ######### 0.23
  cfg.list["person.eagerness.man.dist.gamma.b"] <- inputvector[8] ######### 45
  cfg.list["person.eagerness.woman.dist.gamma.b"] <- inputvector[9] ######### 45


  cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
  cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
  cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1

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

  # Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly
  art.intro <- list()
  art.intro["time"] <- 25
  art.intro["person.art.accept.threshold.dist.fixed.value"] <- 0.5 # inputvector[4] ######### 0.5
  art.intro["diagnosis.baseline"] <- 0#100
  art.intro["monitoring.cd4.threshold"] <- 100 # 1200
  #art.intro["monitoring.interval.piecewise.cd4s"] <- "0,1300"


  # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500

  art.intro2 <- list()
  art.intro2["time"] <- 25 + 5 # inputvector[5] ######### 30
  art.intro2["monitoring.cd4.threshold"] <- 200

  art.intro3 <- list()
  art.intro3["time"] <- 33 # inputvector[4] + inputvector[5] + inputvector[6] ########### 33
  art.intro3["monitoring.cd4.threshold"] <- 350

  art.intro4 <- list()
  art.intro4["time"] <- 3 # inputvector[4] + inputvector[5] + inputvector[6] + inputvector[7] ########### 36
  art.intro4["monitoring.cd4.threshold"] <- 500

  art.intro5 <- list()
  art.intro5["time"] <- 38
  art.intro5["monitoring.cd4.threshold"] <- 5000 # This is equivalent to immediate access
  art.intro5["person.art.accept.threshold.dist.fixed.value"] <- 0.5 # inputvector[8] ########### 0.75

  # tasp.indicator <- inputvector[9] # 1 if the scenario is TasP, 0 if the scenario is current status

  interventionlist <- list(art.intro, art.intro2, art.intro3, art.intro4, art.intro5)

  intervention <- interventionlist # scenario(interventionlist, tasp.indicator)

  results <- tryCatch(simpact.run(configParams = cfg,
                                  destDir = destDir,
                                  agedist = age.distr,
                                  seed = seedid, #, Introducing ART has helped to keep the prevalence high
                                  intervention = intervention),
                      error = simpact.errFunction)

  if (length(results) == 0){
    outputvector <- rep(NA, 16)
  } else {
    if (as.numeric(results["eventsexecuted"]) >= (as.numeric(cfg["population.maxevents"]) - 1)) {
      outputvector <- rep(NA, 16)
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
                        exp(growthrate))
    }
  }
  return(outputvector)
}
