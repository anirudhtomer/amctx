library(ggplot2)
library(ggmcmc)
library(coda)
library(parallel)
library(doParallel)
library(survival)
library(splines)
library(nlme)
library(JMbayes)

ticksX = function(from=0, max, by, labels=waiver()){
  scale_x_continuous(breaks = seq(from, max, by = by), labels = labels)
}

ticksY = function(from=0, max, by, labels = waiver()){
  scale_y_continuous(breaks = round(seq(from, max, by = by),3), labels=waiver())
}

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

source("../JMBayes/Anirudh/dev/replaceMCMCContentsLight.R")
source("../JMBayes/Anirudh/dev/multiplot.R")

plotRandomProfile = function(count=1, fitted=F, creatinine=T){
  
  dataset = if(creatinine==T){amctx_creatinine}else{amctx_pcr}
  
  pid_sample = sample(x = unique(dataset$amctx), size = count)
  plot<-ggplot(data=dataset[dataset$amctx %in% pid_sample,], aes(x=tx_s_years, y=log(value))) + 
    geom_line(aes(group=amctx))
  if(fitted==T){
    plot = plot + geom_line(aes(y=fitted, x=tx_s_years, color=amctx, group=amctx)) 
  }else{
    plot
  }
  plot = plot + xlab("Time (years)") +  ggtitle(ifelse(creatinine==T, "Creatinine", "PCR"))
  print(plot)
}

getFittedSlopes = function(object){
  postRandEff = object$statistics$postMeans$b[, 6:10]
  Zmatrix = object$model_info$mvglmer_components$Z2
  Ids = mvJoint_pcr_creatinine_tdboth_complex$model_info$mvglmer_components$id2
  rowSums(Zmatrix * postRandEff[as.numeric(Ids),,drop=F])
}

#######################################################################
# the following function creates the predicted values
# and the 95% CIs
##################################################
effectPlotData <- function (object, newdata, orig_data) {
  form <- formula(object)
  namesVars <- all.vars(form)
  betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
  V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
  orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
  Terms <- delete.response(terms(form))
  mfX <- model.frame(Terms, data = orig_data)
  Terms_new <- attr(mfX, "terms")
  mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
  X <- model.matrix(Terms_new, mfX_new)
  pred <- c(X %*% betas)
  ses <- sqrt(diag(X %*% V %*% t(X)))
  newdata$pred <- pred
  newdata$low <- pred - 1.96 * ses
  newdata$upp <- pred + 1.96 * ses
  newdata
}

plotPCRFittedCurve = function(models, transform=F, individually=T){
  
  ds = amctx_merged[!is.na(amctx_merged$pcr),]
  newDF <- with(ds, expand.grid(tx_s_years = seq(0, 10, length.out = 30),
                                rec_gender = c("F","M"),
                                d_age = median(amctx.id$d_age),
                                d_cadaveric=c("no", "yes")))
  
  
  plotData = lapply(models, function(model){
    temp = effectPlotData(model, newDF, ds)
    
    if(transform==T){
      temp[,c("pred","low","upp")] = exp(temp[,c("pred","low","upp")])
    }
    temp
  })
  
  if(individually==T){
    lapply(plotData, function(data){
      plot = ggplot(data=data)
      if(transform==F){
        plot = plot + geom_ribbon(aes(x=tx_s_years, ymin=low, ymax=upp), fill = "grey", alpha=0.4)
      } 
      print(plot + geom_line(aes(y=pred, x=tx_s_years), color="black") +
              facet_grid(rec_gender~d_cadaveric, labeller = label_both) + 
              ticksX(from=0, max = 10, by = 1, labels = c(0:10)) + 
              xlab("Follow-up time (Years)") + 
              ylab("PCR"))
    })
  }else{
    newPlotData = do.call(rbind, plotData)
    newPlotData$model = rep(paste("Model",1:length(plotData)), each=nrow(newDF))
    
    plot = ggplot(data=newPlotData)
    if(transform==F){
      plot = plot + geom_ribbon(aes(x=tx_s_years, ymin=low, ymax=upp, group=model), 
                                fill = "grey", alpha=0.4)
    } 
    print(plot + geom_line(aes(y=pred, x=tx_s_years, group=model, color=model)) +
            facet_grid(rec_gender~d_cadaveric, labeller = label_both) + 
            ticksX(from=0, max = 10, by = 1, labels = c(0:10)) + 
            xlab("Follow-up time (Years)") + 
            ylab("PCR"))
  }
}

fitUnivariatePCRModel = function(fixedSplineKnots=c(50, 200, 365)/365, 
                                 randomSplineKnots=c(50, 200)/365, 
                                 boundaryKnots=range(amctx_pcr$tx_s_years), method="ML"){
  
  spline = paste("ns(tx_s_years, knots=c(", paste(fixedSplineKnots, collapse=", "), 
                 "), Boundary.knots=c(", paste(boundaryKnots, collapse=", "),"))", sep="")
  
  fixedFormula = as.formula(paste("log(pcr) ~ d_age + rec_gender + d_cadaveric + ",
                                  spline, "* d_cadaveric + ", 
                                  spline, "* rec_gender", sep=""))
  
  randomFormula = as.formula(paste("~ns(tx_s_years, knots=c(", paste(randomSplineKnots, collapse=", "), 
                                   "), Boundary.knots=c(", paste(boundaryKnots, collapse=", "),"))|amctx", sep=""))
  
  model = lme(fixed=fixedFormula, random = randomFormula, 
              data=amctx_merged[!is.na(amctx_merged$pcr),],
              control = lmeControl(opt = "optim"), 
              method = method)
  
  model$call$fixed = fixedFormula
  model$call$random = randomFormula
  
  print(anova(model, type="marginal"))
  plot = qplot(x=model$fitted[,2], y=model$residuals[,2], xlab="Fitted", 
               ylab="Residuals", geom=c("point", "smooth"))
  print(plot)
  
  return(model)
}


plotCreatinineFittedCurve = function(models, transform=F, individually=T){
  
  ds = amctx_creatinine[!is.na(amctx_creatinine$creatinine),]
  newDF <- data.frame(tx_s_years = seq(0, 10, length.out = 30),
                      amctx.id[1,])
  effectPlotData(models[[1]], newDF, ds)  
  plotData = lapply(models, function(model){
    temp = effectPlotData(model, newDF, ds)
    
    if(transform==T){
      temp[,c("pred","low","upp")] = exp(temp[,c("pred","low","upp")])
    }
    temp
  })
  
  if(individually==T){
    lapply(plotData, function(data){
      plot = ggplot(data=data)
      if(transform==F){
        plot = plot + geom_ribbon(aes(x=tx_s_years, ymin=low, ymax=upp), fill = "grey", alpha=0.4)
      } 
      print(plot + geom_line(aes(y=pred, x=tx_s_years), color="black") +
              facet_grid(tx_dgf~d_cadaveric, labeller = label_both) + 
              ticksX(from=0, max = 10, by = 1, labels = c(0:10)) + 
              xlab("Follow-up time (Years)") + 
              ylab("Creatinine"))
    })
  }else{
    newPlotData = do.call(rbind, plotData)
    newPlotData$model = rep(paste("Model",1:length(plotData)), each=nrow(newDF))
    
    plot = ggplot(data=newPlotData)
    if(transform==F){
      plot = plot + geom_ribbon(aes(x=tx_s_years, ymin=low, ymax=upp, group=model), 
                                fill = "grey", alpha=0.4)
    } 
    print(plot + geom_line(aes(y=pred, x=tx_s_years, group=model, color=model)) +
            facet_grid(tx_dgf~d_cadaveric, labeller = label_both) + 
            ticksX(from=0, max = 10, by = 1, labels = c(0:10)) + 
            xlab("Follow-up time (Years)") + 
            ylab("Creatinine"))
  }
}

fitUnivariateCreatinineModel = function(fixedSplineKnots=c(50, 100, 900)/365, 
                                        randomSplineKnots=c(50, 100)/365, 
                                        boundaryKnots=range(amctx_creatinine$tx_s_years), method="ML"){
  
  spline = paste("ns(tx_s_years, knots=c(", paste(fixedSplineKnots, collapse=", "), 
                 "), Boundary.knots=c(", paste(boundaryKnots, collapse=", "),"))", sep="")
  
  fixedFormula = as.formula(paste("log(creatinine) ~ rec_age_fwp1 + rec_gender + d_age + tx_pra + ah_nr + tx_dm +",
                                  spline, "* d_cadaveric + ", 
                                  spline, "* tx_dgf", sep=""))
  
  randomFormula = as.formula(paste("~ns(tx_s_years, knots=c(", paste(randomSplineKnots, collapse=", "), 
                                   "), Boundary.knots=c(", paste(boundaryKnots, collapse=", "),"))|amctx", sep=""))
  
  model = lme(fixed=fixedFormula, random = randomFormula, 
              data=amctx_creatinine,
              control = lmeControl(opt = "optim"), 
              method = method)
  
  model$call$fixed = fixedFormula
  model$call$random = randomFormula
  
  print(anova(model, type="marginal"))
  plot = qplot(x=model$fitted[,2], y=model$residuals[,2], xlab="Fitted", 
               ylab="Residuals", geom=c("point", "smooth"))
  print(plot)
  
  return(model)
}

plotAUCOverTime = function(auc_roc_list){
  auc = unlist(lapply(auc_roc_list, function(fold){
    lapply(fold, function(auc_roc_types){
      lapply(auc_roc_types[[1]], function(auc_roc_types_instance){
        auc_roc_types_instance$auc
      })
    })
  }))
  tstart = unlist(lapply(auc_roc_list, function(fold){
    lapply(fold, function(auc_roc_types){
      lapply(auc_roc_types[[1]], function(auc_roc_types_instance){
        auc_roc_types_instance$Tstart
      })
    })
  }))
  
  aucMean = c(by(auc, tstart, function(x){mean(x, na.rm = T)}))
  aucLow = c(by(auc, tstart, function(x){quantile(x, probs = 0.025, na.rm = T)}))
  aucUp = c(by(auc, tstart, function(x){quantile(x, probs = 0.975, na.rm = T)}))
  aucTimes = unique(tstart)
  
  plotData = data.frame(aucTimes, aucMean, aucUp, aucLow)
  plot = ggplot(data=plotData) + geom_line(aes(y=aucMean, x=aucTimes)) + 
    geom_ribbon(aes(x=aucTimes, ymin = aucLow, ymax = aucUp), fill="grey", alpha=0.5)
  
  print(plot)
  return(plotData)
}

plotDynamicRiskProb = function(pid, fittedJointModel, origdata, maxLongTimeHorizon=Inf, 
                               maxFutureTime=NA, futureTimeDt = NA, riskLabelProb = NA,
                               RISK_COLOR="firebrick2", ymax=NA,
                               POINT_COLOR="dodgerblue4", FONT_SIZE=14, POINT_SIZE=3){
  
  if(is.na(futureTimeDt) & is.na(maxFutureTime)){
    stop("No future max time or time window supplied")
  }
  
  patientDs = origdata[origdata$amctx %in% pid & origdata$tx_s_years<=maxLongTimeHorizon &
                         !is.na(origdata$creatinine),]
  lastCreatinineTime = max(patientDs$tx_s_years)
  
  if(is.na(maxFutureTime) & !is.na(futureTimeDt)){
    maxFutureTime = lastCreatinineTime + futureTimeDt
  }
  
  futureTimes = seq(lastCreatinineTime, maxFutureTime, length.out = 10)
  
  sfit = survfitJM(fittedJointModel, patientDs, idVar="amctx", survTimes = futureTimes)
  
  longprof = predict(fittedJointModel, patientDs, type = "Subject",
                     interval = "confidence", return = TRUE, idVar="amctx", FtTimes = futureTimes)
  
  longprof[longprof$tx_s_years>lastCreatinineTime,]$creatinine=NA
  longprof$riskMean = rep(NA, nrow(longprof))
  longprof$riskLow = rep(NA, nrow(longprof))
  longprof$riskUp = rep(NA, nrow(longprof))
  
  ymin = min(c(longprof[longprof$tx_s_years<=lastCreatinineTime,]$pred, log(longprof[longprof$tx_s_years<=lastCreatinineTime,]$creatinine)))
  
  if(is.na(ymax)){
    ymax = max(c(longprof[longprof$tx_s_years<=lastCreatinineTime,]$pred, log(longprof[longprof$tx_s_years<=lastCreatinineTime,]$creatinine)))
  }
  riskProbRes = 1 - rbind(c(1,1,1), sfit$summaries[[1]][, c("Mean", "Lower", "Upper")])
  
  #subsetting twice because there are two rows for the last time, and -1 to remove one of those two rows
  longprof[longprof$tx_s_years>=lastCreatinineTime, c("riskMean", "riskLow", "riskUp")][-1,] =
    (riskProbRes * (ymax-ymin) + ymin)
  
  p=ggplot() +
    geom_line(data = longprof[longprof$tx_s_years<=lastCreatinineTime,], aes(x = tx_s_years, y=pred), color=POINT_COLOR) +
    geom_point(data = longprof[longprof$tx_s_years<=lastCreatinineTime,], aes(y=log(creatinine), x=tx_s_years),
               color=POINT_COLOR, size=POINT_SIZE) +
    geom_vline(xintercept = lastCreatinineTime, linetype="dashed") +
    geom_line(data = longprof, aes(x=tx_s_years, y=riskMean), color=RISK_COLOR) +
    geom_ribbon(data = longprof, aes(ymin=riskLow, ymax=riskUp, x= tx_s_years), fill=RISK_COLOR, alpha=0.2) +
    xlab("Time since transplantation (years)") + ylab(expression('Serum creatinine (log umol/l)')) +
    theme_bw() +
    theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE), 
          axis.text.y.right = element_text(color=RISK_COLOR), axis.title.y.right = element_text(color=RISK_COLOR),
          plot.title = element_text(hjust = 0.5, size=FONT_SIZE)) +
    scale_x_continuous(breaks = seq(0, maxFutureTime, 1)) + 
    scale_y_continuous(limits = c(ymin, ymax),breaks = round(seq(ymin, ymax, length.out = 5),2), 
                       sec.axis = sec_axis(~(.-ymin)/(ymax-ymin), name = "Dynamic risk of graft failure (%)", 
                                           breaks = c(0, 0.05, 0.5, 1),
                                           labels = paste0(c(0, 0.05, 0.5, 1)*100, "%")))
  
  if(!is.na(riskLabelProb) & riskLabelProb < max(riskProbRes)){
    riskProbResTime = futureTimes[which.min(abs(riskProbRes[,1]-riskLabelProb))]
    print(paste0("Time at which risk prob = ", riskLabelProb, " is ", riskProbResTime))
    print(paste0("The scaled risk prob = ", riskLabelProb * (ymax-ymin) + ymin))
  }
  
  print(p)
  return(p)
}

plotInformationGain = function(pid, fittedJointModel, origdata, maxLongTimeHorizon, 
                               futureTimeDt, INFO_COLOR="forestgreen",
                               POINT_COLOR="dodgerblue4", FONT_SIZE=14, POINT_SIZE=3){
  
  if(is.na(futureTimeDt)){
    stop("No future max time or time window supplied")
  }
  
  patientDs = origdata[origdata$amctx %in% pid & origdata$tx_s_years<=maxLongTimeHorizon &
                         !is.na(origdata$creatinine),]
  lastCreatinineTime = max(patientDs$tx_s_years)
  
  futureTimes = seq(lastCreatinineTime, lastCreatinineTime + futureTimeDt, length.out = 16)
  
  longprof = predict(fittedJointModel, patientDs, type = "Subject",
                     interval = "confidence", return = TRUE, idVar="amctx", FtTimes = futureTimes)
  
  dynInfoRes = dynInfoPar(fittedJointModel, newdata = patientDs, Dt = futureTimeDt,
                          K=16, seed = 2018, idVar="amctx")
  
  informationGain = dynInfoRes$summary$Info - min(dynInfoRes$summary$Info)
  
  longprof[longprof$tx_s_years>lastCreatinineTime,]$creatinine=NA
  longprof$informationGain = rep(NA, nrow(longprof))
  
  ymin = min(c(longprof[longprof$tx_s_years<=lastCreatinineTime,]$pred, log(longprof[longprof$tx_s_years<=lastCreatinineTime,]$creatinine)))
  ymax = max(c(longprof[longprof$tx_s_years<=lastCreatinineTime,]$pred, log(longprof[longprof$tx_s_years<=lastCreatinineTime,]$creatinine)))
  
  scaledInformationGain = (informationGain - min(informationGain)) * (ymax-ymin)/(max(informationGain)-min(informationGain)) + ymin
  #subsetting twice because there are two rows for the last time, and -1 to remove one of those two rows
  longprof[longprof$tx_s_years>=lastCreatinineTime, c("informationGain")][-1] = scaledInformationGain
  
  p=ggplot() +
    geom_line(data = longprof[longprof$tx_s_years<=lastCreatinineTime,], aes(x = tx_s_years, y=pred), color=POINT_COLOR) +
    geom_point(data = longprof[longprof$tx_s_years<=lastCreatinineTime,], aes(y=log(creatinine), x=tx_s_years),
               color=POINT_COLOR, size=POINT_SIZE) +
    geom_vline(xintercept = lastCreatinineTime, linetype="dashed") +
    geom_line(data = longprof, aes(x=tx_s_years, y=informationGain), color=INFO_COLOR) + 
    scale_y_continuous(limits = c(ymin, ymax + 0.2 * (ymax - ymin)),labels = round(seq(ymin, ymax + 0.2 * (ymax - ymin), length.out = 5),2),
                       breaks = seq(ymin, ymax + 0.2 * (ymax - ymin), length.out = 5),
                       sec.axis = sec_axis(~(.-ymin)*(max(informationGain)-min(informationGain))/(ymax-ymin) + min(informationGain),
                                           breaks  = seq(min(informationGain), max(informationGain), length.out = 5),
                                           labels = round(seq(min(informationGain), max(informationGain), length.out = 5),2),
                                           name = "Information from new measurement")) + 
    theme_bw() +
    theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE), 
          axis.text.y.right = element_text(color=INFO_COLOR), 
          axis.title.y.right = element_text(color=INFO_COLOR),
          plot.title = element_text(hjust = 0.5, size=FONT_SIZE)) +
    xlab("Time since transplantation (years)") + ylab(expression('Serum creatinine (log umol/l)'))
  
  print(paste0("Time at which information gain is maximum = ", futureTimes[which.max(informationGain)]))
  print(paste0("The scaled max info gain is = ", max(scaledInformationGain)))
  
  print(p)
  return(p)
}

#p1 = plotDynamicRiskProb(190, joint_creatinine_tdboth_complex_replaced, maxLongTimeHorizon=3.5, maxFutureTime = 6) + ggtitle("Using creatinine measurements up to 3.5 years")
#p2 = plotDynamicRiskProb(190, joint_creatinine_tdboth_complex_replaced, maxLongTimeHorizon=4, maxFutureTime = 6) + ggtitle("Using creatinine measurements up to 4 years")
#ggsave(multiplot(p1, p2, cols = 1), filename = "report/hessel/images/dynriskProb.eps", width=8.27, height=8.27, device=cairo_ps, dpi=500)
