root <- getwd()
source(file.path(root, "Rpgm/Simulation/simulateFuns.R"))
library("JMbayes")

##########################################################################################
##########################################################################################

cutoff <- 0.8
gap <- 2
Dt <- 5

# simulate data
Data <- simulateJoint()

# Fit joint model
lmeFit <- lme(y ~ ns(time, k = c(2.1, 5.5), B = c(0, 19.5)), data = Data$DF,
              random = ~ ns(time, k = c(2.1, 5.5), B = c(0, 19.5)) | id,
              control = lmeControl(opt = "optim", niterEM = 45))
coxFit <- coxph(Surv(Time, event) ~ group, data = Data$DF.id, x = TRUE)
dForm <- list(fixed = ~ 0 + dns(time, k = c(2.1, 5.5), B = c(0, 19.5)), 
              random = ~ 0 + dns(time, k = c(2.1, 5.5), B = c(0, 19.5)),
              indFixed = 2:4, indRandom = 2:4)
jointFit <- jointModelBayes(lmeFit, coxFit, timeVar = "time", extraForm = dForm,
                            param = "td-both")
summary(jointFit)

true_opt_times <- true_optimal_times(Data$b, Data$predDat, Data$trueValues, Dt, cutoff)
system.time(opt_times_fixed <- optimalTime("fixed", Data, jointFit, cutoff, gap, Dt))
system.time(opt_times_pers <- optimalTime("personalized", Data, jointFit, cutoff, gap, Dt))

Effect1 <- abs(opt_times_pers$optimal_times - true_opt_times)
Effect0 <- abs(opt_times_fixed$optimal_times - true_opt_times)
Diff_trueTime1 <- (opt_times_pers$optimal_times - Data$predDat$Time[!duplicated(Data$predDat$id)])
Diff_trueTime0 <- (opt_times_fixed$optimal_times - Data$predDat$Time[!duplicated(Data$predDat$id)])
Cost1 <- sapply(opt_times_pers$visit_times, length)
Cost0 <- sapply(opt_times_fixed$visit_times, length)
ICER <- (Cost1 - Cost0) / (Effect1 - Effect0)
boxplot(list(Effect1, Effect0))
boxplot(list(Diff_trueTime1, Diff_trueTime0))
boxplot(list(Cost1, Cost0))
summary(Effect1)
summary(Effect0)
summary(Diff_trueTime1)
summary(Diff_trueTime0)
summary(Cost1)
summary(Cost0)



cbind("true" = true_opt_times, "Event" = Data$predDat$Time[!duplicated(Data$predDat$id)], 
      "Pers" = opt_times_pers$optimal_times, "Fixed" = opt_times_fixed$optimal_times)


fff <- function (i) {
    list(fixed = opt_times_fixed$visit[[i]], personalized = opt_times_pers$visit[[i]])
}

##########################################################################################

rootOut <- file.path(root, "LaTex Files/vers2/Figures")
pdf(file.path(rootOut, "Simulation.pdf"), height = 5)
op <- par(mfrow = c(1, 2))
boxplot(list("Fixed" = Effect0, "Personalized" = Effect1), 
        ylab = "Absolute Error Optimal Intervention Point")
boxplot(list("Fixed" = Cost0, "Personalized" = Cost1), ylab = "Number of Screenings")
par(op)
dev.off()

##########################################################################################


library("parallel")

runSim <- function (inds) {
    K <- length(inds)
    res <- vector("list", K)
    for (k in 1:K) {
        set.seed(max(inds) + k)
        test <- try({
            root <- "C:/Users/Dimitris/Documents/Papers/Paper23"
            source(file.path(root, "Rpgm/Simulation/simulateFuns.R"))
            library("JMbayes")
            
            ##########################################################################################
            ##########################################################################################
            
            cutoff <- 0.8
            gap <- 2
            Dt <- 5
            
            # simulate data
            Data <- simulateJoint()
            
            # Fit joint model
            lmeFit <- lme(y ~ ns(time, k = c(2.1, 5.5), B = c(0, 19.5)), data = Data$DF,
                          random = ~ ns(time, k = c(2.1, 5.5), B = c(0, 19.5)) | id,
                          control = lmeControl(opt = "optim", niterEM = 45))
            coxFit <- coxph(Surv(Time, event) ~ group, data = Data$DF.id, x = TRUE)
            dForm <- list(fixed = ~ 0 + dns(time, k = c(2.1, 5.5), B = c(0, 19.5)), 
                          random = ~ 0 + dns(time, k = c(2.1, 5.5), B = c(0, 19.5)),
                          indFixed = 2:4, indRandom = 2:4)
            jointFit <- jointModelBayes(lmeFit, coxFit, timeVar = "time", extraForm = dForm,
                                        param = "td-both")
            
            true_opt_times <- true_optimal_times(Data$b, Data$predDat, Data$trueValues, Dt, cutoff)
            system.time(opt_times_fixed <- optimalTime("fixed", Data, jointFit, cutoff, gap, Dt))
            system.time(opt_times_pers <- optimalTime("personalized", Data, jointFit, cutoff, gap, Dt))
            
            Effect1 <- abs(opt_times_pers$optimal_times - true_opt_times)
            Effect0 <- abs(opt_times_fixed$optimal_times - true_opt_times)
            Diff_trueTime1 <- (opt_times_pers$optimal_times - Data$predDat$Time[!duplicated(Data$predDat$id)])
            Diff_trueTime0 <- (opt_times_fixed$optimal_times - Data$predDat$Time[!duplicated(Data$predDat$id)])
            Cost1 <- sapply(opt_times_pers$visit_times, length)
            Cost0 <- sapply(opt_times_fixed$visit_times, length)
            ICER <- (Cost1 - Cost0) / (Effect1 - Effect0)
            list(Effect1 = Effect1, Effect0 = Effect0, Cost1 = Cost1, Cost0 = Cost0, 
                 ICER = ICER)
        }, TRUE)
        res[[k]] <- if (!inherits(test, "try-error")) test
    }
    res
}

splits <- matrix(1:16, 2, 8)
splits <- split(splits, col(splits))
cl <- makeCluster(8)
res <- parLapply(cl, splits, runSim)
stopCluster(cl)

res. <- unlist(res, recursive = FALSE, use.names = FALSE)
Effect1 <- c(sapply(res., function (x) x$Effect1))
Effect0 <- c(sapply(res., function (x) x$Effect0))
Cost1 <- c(sapply(res., function (x) x$Cost1))
Cost0 <- c(sapply(res., function (x) x$Cost0))
Cost1 <- Cost1[Cost1 >= 3]
Cost0 <- Cost0[Cost0 >= 3]



#####

rootOut <- file.path(root, "LaTex Files/vers2/Figures")
pdf(file.path(rootOut, "Simulation.pdf"), height = 5)
op <- par(mfrow = c(1, 2))
boxplot(list("Fixed" = Effect0, "Personalized" = Effect1), 
        ylab = "Absolute Error Optimal Intervention Point")
boxplot(list("Fixed" = Cost0, "Personalized" = Cost1), ylab = "Number of Screenings")
par(op)
dev.off()

