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
  scale_y_continuous(breaks = seq(from, max, by = by), labels=waiver())
}

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

source("../JMBayes/Anirudh/dev/replaceMCMCContents.R")

plotRandomProfile = function(count=1, fitted=F, creatinine=T){
  
  dataset = if(creatinine==T){amctx_creatinine}else{amctx_pcr}
  
  pid_sample = sample(x = unique(dataset$amctx), size = count)
  plot<-ggplot(data=dataset[dataset$amctx %in% pid_sample,], aes(x=tx_s_years, y=log(value))) + 
    geom_line(aes(group=amctx))
  if(fitted==T){
    plot + geom_line(aes(y=fitted, x=tx_s_years, color=amctx, group=amctx)) 
  }else{
    plot
  }
  plot + xlab("Time (years)") +  ggtitle(ifelse(creatinine==T, "Creatinine", "PCR"))

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

fitUnivaritePCRModel = function(fixedSplineKnots=c(50, 200, 365)/365, 
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

  ds = amctx_merged[!is.na(amctx_merged$creatinine),]
  newDF <- with(ds, expand.grid(tx_s_years = seq(0, 10, length.out = 30),
                                   rec_gender = "M",
                                   rec_age_fwp1 = median(amctx.id$rec_age),
                                   d_age = median(amctx.id$d_age),
                                   tx_pra = median(amctx.id$tx_pra),
                                   tx_dm = "no", ah_nr=median(amctx.id$ah_nr),
                                   tx_dgf = c("no", "yes"), 
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

fitUnivariteCreatinineModel = function(fixedSplineKnots=c(50, 100, 900)/365, 
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
              data=amctx_merged[!is.na(amctx_merged$creatinine),],
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

# 
# idList = droplevels(unique(amctx.id[amctx.id$gl_failure==0,]$amctx))
# 
# idList = c(73, 94, 195, 209)
# for(id in idList){
#   ND = amctx_creatinine[amctx_creatinine$amctx %in% id,]
#   futureTimes = seq(max(ND$tx_s_years), (max(ND$tx_s_years) + 3), 0.1)
#   
#   ####################
#   longprof = predict(jmbayes_creatinine_orig, ND, type = "Subject",
#           interval = "confidence", return = TRUE, idVar="amctx", FtTimes = futureTimes)
#   last.time <- with(longprof, tx_s_years[!is.na(low)][1])
#   longprof[longprof$tx_s_years>last.time,]$value=NA
#   ggplot(data = longprof, aes(x = tx_s_years, y=pred)) + geom_line() + 
#     geom_ribbon(aes(ymin=low, ymax=upp), fill="grey", alpha=0.5) + 
#     geom_point(aes(y=log(value)), colour="red", alpha=0.4) + 
#     geom_vline(xintercept = last.time, linetype="dotted") + 
#     xlab("Time (years)") + ylab("Predicted log(serum creatinine)") + 
#     ggtitle(paste("amctx =",id))
#   
#   #####################
#   sfit.patient2 = survfitJM(jmbayes_creatinine_orig, ND, idVar="amctx", survTimes = futureTimes)
#   plot(sfit.patient2, estimator="mean", include.y=T, conf.int=T, fill.area=T, col.area="lightgrey", main=paste("amctx =",id))
#   
#   longprof$survMean = rep(NA, nrow(longprof))
#   longprof$survLow = rep(NA, nrow(longprof))
#   longprof$survUp = rep(NA, nrow(longprof))
#   
#   ymin = min(c(longprof[longprof$tx_s_years<=last.time,]$pred, log(longprof[longprof$tx_s_years<=last.time,]$value)))
#   ymax = max(c(longprof[longprof$tx_s_years<=last.time,]$pred, log(longprof[longprof$tx_s_years<=last.time,]$value)))
#   
#   longprof[longprof$tx_s_years>last.time, c("survMean", "survLow", "survUp")] = 
#     (sfit.patient2$summaries[[1]][, c("Mean", "Lower", "Upper")] * (ymax-ymin) + ymin)
#   
#   ggplot() + 
#     geom_line(data = longprof[longprof$tx_s_years<=last.time,], aes(x = tx_s_years, y=pred)) + 
#     geom_point(data = longprof, aes(y=log(value), x=tx_s_years), colour="red", alpha=0.4) + 
#     geom_vline(xintercept = last.time, linetype="dotted") + 
#     geom_line(data = longprof, aes(x=tx_s_years, y=survMean)) + 
#     geom_ribbon(data = longprof, aes(ymin=survLow, ymax=survUp, x= tx_s_years), fill="grey", alpha=0.5) + 
#     xlab("Time (years)") + ylab("log(serum creatinine)") + 
#     scale_y_continuous(limits = c(ymin, ymax), sec.axis = sec_axis(~(.-ymin)/(ymax-ymin)))
#   
#   plotSurv = ggplot(data=data.frame(sfit.patient2$summaries[[1]])) + 
#     geom_line(aes(x=times, y=Mean)) + 
#     geom_ribbon(aes(ymin=Lower, ymax=Upper, x= times),  fill="grey", alpha=0.5) + 
#     scale_y_continuous(position = "right")
#   
#   multiplot(plotLong, plotSurv, cols = 2)
#   multiplot(templot, plotSurv, cols = 2)
# }
# 
# rocJM(jmbayes_creatinine_orig, dt = c(1, 2, 4), data = ND,
#       M = 1000, burn.in = 500)
