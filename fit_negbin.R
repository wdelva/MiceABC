library(dplyr)
library(MASS)
library(splines)
library(boot)
library(ggplot2)

data<-read_dta("C:/Users/Christianah/Google Drive/PTOR/Sex Network Master (participant level) 09.08.2012.dta") #
names(data)
data.needed <- dplyr::select(data, c(id,
                                     age,
                                     race,
                                     gender,
                                     cpcount,
                                     ompcount,
                                     hivtestresult,
                                     totalpartners))
data.needed$cpcount[is.na(data.needed$cpcount)] <- 0
data.needed$ompcount[is.na(data.needed$ompcount)] <- 0

#########################################################################################
# data.filtered is the dataset for the partner turnover analysis

data.filtered <- dplyr::filter(data.needed,
                               race == 2,
                               gender != "2",
                               age > 15,
                               age < 41)

data.filtered$totalpartnerslstyr <- data.filtered$cpcount +
  data.filtered$ompcount # Total NEW relationships formed in the last year

# We create 4 subsets, one for each race-gender stratum
md.wideBM <- dplyr::filter(data.filtered,
                           race == 2,
                           gender == "0")
md.wideBW <- dplyr::filter(data.filtered,
                           race == 2,
                           gender == "1")
##########################################################################################

diffest<- function(md.wide,indices){
  d <- md.wide[indices,] # taking a random sample of the full dataset
  final<-glm.nb(totalpartnerslstyr~ns(age,4),
                data=d,
                control=glm.control(maxit=100000, trace = 3),
                init.theta=0.5) # fitting the PLY model to this sample
  final2<-glm.nb(totalpartners~ns(age,3),
                 data=d,
                 control=glm.control(maxit=100000)) # fitting the LTP model to this sample
  
  made <- data.frame(age=16:40) # creating artificial dataset for model predictions
  
  try <- predict(final, newdata=made, type="response") # model predictions for PLY
  s <- cumsum(try) # cumulating PLY to obtain a model-based estimate of LTP
  
  tryL <- predict(final2, newdata=made, type="response") # model predictions for LTP
  
  diffvect <- s - tryL # Difference between the cumulated estimates and the model-prediction estimates
  return(diffvect)
}
###################################
#Negative binomial model
###################################
# Model for black men in the last year

finalBm<-glm.nb(totalpartnerslstyr~ns(age,4),
                data=md.wideBM,
                control=glm.control(maxit=100000, trace = 3),
                init.theta=0.5)

summary(finalBm)

# Model for black women in the last year

finalBw<-glm.nb(totalpartnerslstyr~ns(age,4),
                data=md.wideBW,
                control=glm.control(maxit=100000, trace = 3),
                init.theta=0.5) # fitting the PLY model to this sample
summary(finalBw)

# Model for the lifetime number of partners for black men 
final2Bm<-glm.nb(totalpartners~ns(age,4),
                 data=md.wideBM,
                 control=glm.control(maxit=100000, trace = 3),
                 init.theta=0.5)

# Model for the lifetime number of partners for black women
final2Bw<-glm.nb(totalpartners~ns(age,4),
                 data=md.wideBW,
                 control=glm.control(maxit=100000, trace = 3),
                 init.theta=0.5) # fitting the LTP model to this sample

###############################################################################
# Made up data set to predict the expected number of sexual partners per age
##############################################################################
madeforBM<-data.frame(age = 16:40, gender = "man", race = "African")
madeforBW<-data.frame(age=16:40, gender="woman", race="African")

madeforBM$totalpartnerslstyr <- NA
madeforBW$totalpartnerslstyr <- NA
########################################################################################
tryBM<-predict(finalBm,newdata=madeforBM,type="response")
tryBW<-predict(finalBw,newdata=madeforBW,type="response")
madeforBM$try<-tryBM;madeforBW$try<-tryBW

madeforBM$s<-cumsum(madeforBM$try)
madeforBW$u<-cumsum(madeforBW$try)

s<-cumsum(tryBM)
u<-cumsum(tryBW)
#######################################################################################
# BML- Black men lifetime
######################################
madeforBML<-data.frame(age=16:40, gender = "man", race = "African")
madeforBWL<-data.frame(age=16:40, gender="woman", race="African")
madeforBML$totalpartners <- NA
madeforBWL$totalpartners <- NA

tryBML<-predict(final2Bm,newdata=madeforBML,type="response") # model predictions for LTP
tryBWL<-predict(final2Bw,newdata=madeforBWL,type="response")
madeforBML$tryL<-tryBML;madeforBWL$tryL<-tryBWL

madeforBML$s<-s
madeforBWL$u<-u

madeforBML$YL <- (madeforBML$s - madeforBML$tryL)
madeforBWL$YL <- (madeforBWL$u - madeforBWL$tryL)

diffvectBM<-s-tryBML # Difference between the cumulated estimates and the model-prediction estimates
diffvectBW<-u-tryBWL
diffvectALL <- c(diffvectBM, diffvectBW)

############################
# We construct confidence band around our values.
# We do this 4 times: once for each stratum

# 1. For Black Men
t0 <- proc.time()
results<-boot(data=md.wideBM, statistic = diffest, R=50)
proc.time() - t0

lowerbandBM <-rep(NA,length(results$t0))
upperbandBM <- lowerbandBM
for (i in 1:length(results$t0) ) {
  BCI <- boot.ci(results, type="perc", index=i)
  lowerbandBM[i] <- BCI$percent[4]
  upperbandBM[i] <- BCI$percent[5]
}


# 2. For Black Women
t0 <- proc.time()
results<-boot(data=md.wideBW, statistic = diffest, R=50)
proc.time() - t0

lowerbandBW <-rep(NA,length(results$t0))
upperbandBW <- lowerbandBW
for (i in 1:length(results$t0) ) {
  BCI <- boot.ci(results, type="perc", index=i)
  lowerbandBW[i] <- BCI$percent[4]
  upperbandBW[i] <- BCI$percent[5]
}

### Final step: plotting the results
plot(seq(16,40,length=25),diffvectBM, type="l", col="navyblue",main="Black men",xlab="Age",ylab="Difference",ylim=c(-20,75))
lines(seq(16,40,length=25),lowerbandBM[1:25], lty=2)
lines(seq(16,40,length=25),upperbandBM[1:25], lty=2)
lines(c(16,40), c(0,0), lty=3)
#dev.copy(jpeg,"diffBM.jpeg")
#dev.off()
plot(seq(16,40,length=25),diffvectBW, type="l", col="red3",main="Black women",xlab="Age",ylab="Difference")#,ylim=c(-20,75), col="red3")
lines(seq(16,40,length=25),lowerbandBW, lty=2)
lines(seq(16,40,length=25),upperbandBW, lty=2)
lines(c(16,40), c(0,0), lty=3)