# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_121.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
require(rJava)
library(ggplot2)
library(GGally)
library(lmtest)
library(mice)
library(miceadds)
library(RNetLogo)
library(parallel)
library(randtoolbox)
library(RSimpactCyan)
library(RSimpactHelper)
library(EasyABC)
library(dplyr)

nl.path <- "/Applications/NetLogo 6.0/app/"
NLStart(nl.path, gui=FALSE, nl.obj="agemix.nl1", nl.jarname = 'netlogo-6.0.0.jar')
model.path <- "/Users/delvaw/Google Drive/IBM course/tutorials/AgeMixing/Calibration/dummy.0.1.nlogo"
NLLoadModel(model.path, nl.obj="agemix.nl1")
setwd("/Users/delvaw/Documents/MiceABC/R")
master.data <- agemix.simulate(param.set = c(8, 0.2), no.repeated.sim=25)

targets <- c(48.41, 27.91908, 54.22174)
n.experiments <- 30#60
maxit <- 10
maxiter <- 2#6
reps <- 2#5
alpha <- 0.7
rel.dist.tol <- 0.000001
rel.dist.tol.final <- 0.005

test.list <- MiceABC(targets = targets,
                    n.experiments = n.experiments,
                    lls = c(2, 0.05),
                    uls = c(10, 0.65),
                    model = model.wrapper,
                    maxit = maxit,
                    maxiter = maxiter,
                    reps = reps,
                    alpha = alpha,
                    rel.dist.tol = rel.dist.tol)


all.data <- do.call(rbind, test.list$inputoutput)

# We also calculate sum of squared relative distances
all.data$cumsum.transmissions.sq.rel.dist <- ((all.data$cumsum.transmissions - targets[1]) / targets[1])^2
all.data$mean.age.hivpos.sq.rel.dist <- ((all.data$mean.age.hivpos - targets[2]) / targets[2])^2
all.data$var.age.hivpos.sq.rel.dist <- ((all.data$var.age.hivpos - targets[3]) / targets[3])^2
all.data$sum.sq.rel.dist <- all.data$cumsum.transmissions.sq.rel.dist + all.data$mean.age.hivpos.sq.rel.dist +all.data$var.age.hivpos.sq.rel.dist
dist.order <- order(all.data$sum.sq.rel.dist) # Ordering the squared distances from small to big. The last observation (targets) should be ranked first
shortest.dist.percentile <- (n.experiments * alpha) / nrow(all.data)  # 0.1344086 #50/372 # 0.25
last.one.selected <- round(shortest.dist.percentile * length(dist.order))
selected.distances <- dist.order[1:last.one.selected]
all.data.selected <- all.data[selected.distances, ]
all.data.filtered <- dplyr::filter(all.data,
                                   sum.sq.rel.dist <= 0.00346296164344922)#test.list$rel.dist.cutoff[length(test.list$rel.dist.cutoff)])
                            # abs(all.data$cumsum.transmissions - targets[1]) / targets[1] <= rel.dist.tol.final,
                            # abs(all.data$mean.age.hivpos - targets[2]) / targets[2] <= rel.dist.tol.final,
                            # abs(all.data$var.age.hivpos - targets[3]) / targets[3] <= rel.dist.tol.final)

all.data.choices <-list(all.data.selected, all.data.filtered)
all.data.choice.index <- ifelse(nrow(all.data.selected) > nrow(all.data.filtered),
                                1,
                                2)
all.data.final <- all.data.choices[[2]] # all.data.choices[[all.data.choice.index]]
test.list$estimates <- all.data.final

summary(all.data.final)
pairs(all.data.final[,c(2:4, 6:7)])
pairs(all.data.final[,6:7])
var(all.data.final[,6:7])


dim(all.data.final[,6:7])

save.image(file = "/Users/delvaw/Google Drive/IBM course/tutorials/AgeMixing/Calibration/10May.RData")

