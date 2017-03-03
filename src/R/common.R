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


idList = droplevels(unique(amctx.id[amctx.id$gl_failure==0,]$amctx))

idList = c(73, 94, 195, 209)
for(id in idList){
  ND = amctx_creatinine[amctx_creatinine$amctx %in% id,]
  futureTimes = seq(max(ND$tx_s_years), (max(ND$tx_s_years) + 3), 0.1)
  
  ####################
  longprof = predict(jmbayes_creatinine_orig, ND, type = "Subject",
          interval = "confidence", return = TRUE, idVar="amctx", FtTimes = futureTimes)
  last.time <- with(longprof, tx_s_years[!is.na(low)][1])
  longprof[longprof$tx_s_years>last.time,]$value=NA
  ggplot(data = longprof, aes(x = tx_s_years, y=pred)) + geom_line() + 
    geom_ribbon(aes(ymin=low, ymax=upp), fill="grey", alpha=0.5) + 
    geom_point(aes(y=log(value)), colour="red", alpha=0.4) + 
    geom_vline(xintercept = last.time, linetype="dotted") + 
    xlab("Time (years)") + ylab("Predicted log(serum creatinine)") + 
    ggtitle(paste("amctx =",id))
  
  #####################
  sfit.patient2 = survfitJM(jmbayes_creatinine_orig, ND, idVar="amctx", survTimes = futureTimes)
  plot(sfit.patient2, estimator="mean", include.y=T, conf.int=T, fill.area=T, col.area="lightgrey", main=paste("amctx =",id))
  
  longprof$survMean = rep(NA, nrow(longprof))
  longprof$survLow = rep(NA, nrow(longprof))
  longprof$survUp = rep(NA, nrow(longprof))
  
  ymin = min(c(longprof[longprof$tx_s_years<=last.time,]$pred, log(longprof[longprof$tx_s_years<=last.time,]$value)))
  ymax = max(c(longprof[longprof$tx_s_years<=last.time,]$pred, log(longprof[longprof$tx_s_years<=last.time,]$value)))
  
  longprof[longprof$tx_s_years>last.time, c("survMean", "survLow", "survUp")] = 
    (sfit.patient2$summaries[[1]][, c("Mean", "Lower", "Upper")] * (ymax-ymin) + ymin)
  
  ggplot() + 
    geom_line(data = longprof[longprof$tx_s_years<=last.time,], aes(x = tx_s_years, y=pred)) + 
    geom_point(data = longprof, aes(y=log(value), x=tx_s_years), colour="red", alpha=0.4) + 
    geom_vline(xintercept = last.time, linetype="dotted") + 
    geom_line(data = longprof, aes(x=tx_s_years, y=survMean)) + 
    geom_ribbon(data = longprof, aes(ymin=survLow, ymax=survUp, x= tx_s_years), fill="grey", alpha=0.5) + 
    xlab("Time (years)") + ylab("log(serum creatinine)") + 
    scale_y_continuous(limits = c(ymin, ymax), sec.axis = sec_axis(~(.-ymin)/(ymax-ymin)))
  
  plotSurv = ggplot(data=data.frame(sfit.patient2$summaries[[1]])) + 
    geom_line(aes(x=times, y=Mean)) + 
    geom_ribbon(aes(ymin=Lower, ymax=Upper, x= times),  fill="grey", alpha=0.5) + 
    scale_y_continuous(position = "right")
  
  multiplot(plotLong, plotSurv, cols = 2)
  multiplot(templot, plotSurv, cols = 2)
}

rocJM(jmbayes_creatinine_orig, dt = c(1, 2, 4), data = ND,
      M = 1000, burn.in = 500)


# # the data frame that contains the combination of values to
# # create the plot
# newDF <- with(amctx, expand.grid(tx_s_years = seq(0, 15, length.out = 30),
#                                 rec_gender = "M", 
#                                 rec_age_fwp1 = median(amctx.id$rec_age), 
#                                 d_age = median(amctx.id$d_age),
#                                 d_bmi = median(amctx.id$d_bmi),
#                                 rec_bmi = median(amctx.id$rec_bmi),
#                                 tx_dial_days = median(amctx.id$tx_dial_days),
#                                 tx_cit = median(amctx.id$tx_cit),
#                                 tx_pra = median(amctx.id$tx_pra),
#                                 tx_dm = "no", tx_previoustx = "no",
#                                 tx_hla = 3, tx_dgf = "no", d_type="HBD", amctx = "100"))
# 
# # the effects plot
# xyplot(pred + low + upp ~ tx_s_years | rec_gender, 
#        data = effectPlotData(model_pcr, newDF, amctx), 
#        lty = c(1, 2, 2), col = c(2, 1, 1), lwd = 2, type = "l",
#        xlab = "Follow-up time (years)",
#        ylab = "log (pcr)")
