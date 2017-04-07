args <- commandArgs(TRUE)
nsub <- as.numeric(args[1L])

rootIn <- file.path('..', 'Rpgms')
rootOut <- file.path('..', 'Results')

#########################################################################

tt <- try({
    set.seed(2015 + nsub)
    source(file.path(rootIn, 'simulateFuns.R'))
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
                                param = "td-both", seed = 2015 + nsub)
    
    true_opt_times <- true_optimal_times(Data$b, Data$predDat, Data$trueValues, Dt, cutoff)
    opt_times_fixed_1y <- optimalTime("fixed", Data, jointFit, cutoff, gap = 1, Dt)
    opt_times_fixed_2y <- optimalTime("fixed", Data, jointFit, cutoff, gap = 2, Dt)
    opt_times_fixed_3y <- optimalTime("fixed", Data, jointFit, cutoff, gap = 4, Dt)
    opt_times_pers <- optimalTime("personalized", Data, jointFit, cutoff, gap = 3, Dt)
    
    Effect_pers <- abs(opt_times_pers$optimal_times - true_opt_times)
    Effect_fixed_1y <- abs(opt_times_fixed_1y$optimal_times - true_opt_times)
    Effect_fixed_2y <- abs(opt_times_fixed_2y$optimal_times - true_opt_times)
    Effect_fixed_3y <- abs(opt_times_fixed_3y$optimal_times - true_opt_times)
    eventTime <- Data$predDat$Time[!duplicated(Data$predDat$id)]
    zero_length <- function (x) if (is.null(x)) NA else length(x)
    Cost_pers <- sapply(opt_times_pers$visit_times, length)
    Cost_fixed_1y <- sapply(opt_times_fixed_1y$visit_times, zero_length) - 3
    Cost_fixed_2y <- sapply(opt_times_fixed_2y$visit_times, zero_length) - 3
    Cost_fixed_3y <- sapply(opt_times_fixed_3y$visit_times, zero_length) - 3
    

    list(Effect_pers = Effect_pers, Effect_fixed_1y = Effect_fixed_1y, 
         Effect_fixed_2y = Effect_fixed_2y,  Effect_fixed_3y = Effect_fixed_3y, 
         Cost_pers = Cost_pers, Cost_fixed_1y = Cost_fixed_1y, Cost_fixed_2y = Cost_fixed_2y,
         Cost_fixed_3y = Cost_fixed_3y,
         trueTimes = true_opt_times, persTimes = opt_times_pers$optimal_times,
         fixedTimes_1y = opt_times_fixed_1y$optimal_times, 
         fixedTimes_2y = opt_times_fixed_2y$optimal_times, 
         fixedTimes_3y = opt_times_fixed_3y$optimal_times, eventTime = eventTime)
}, TRUE)


out <- if (!inherits(tt, "try-error")) tt else {
    list(Effect1 = NULL, Effect0 = NULL, Cost1 = NULL, Cost0 = NULL, 
         ICER = NULL, trueTimes = NULL, persTimes = NULL, fixedTimes = NULL, 
         eventTime = NULL)
    
}

nam <- paste0('results_', nsub, '.RData')
save(out, file = file.path(rootIn, nam))


boxplot(list(Fixed1 = Effect_fixed_1y, Fixed2 = Effect_fixed_2y, 
             Fixed4 = Effect_fixed_3y, Personalized = Effect_pers))
boxplot(list(Fixed1 = Cost_fixed_1y, Fixed2 = Cost_fixed_2y, Fixed3 = Cost_fixed_3y,
             Personalized = Cost_pers))








