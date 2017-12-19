source("src/R/common.R")

library(MASS)

#Further feedback came back from Clinicians.
# Some considerations to leave out covariates on the basis of knowledge:
#   - DGF can be left out, is not a risk factor for graft failure 
# (at least it could only be a risk factor for very old recipients)
# - HLA can be treated continous
# - diuretics it could be left out as well to my opinion (Further preferably no medication use included in the model)
# - I think donorage is a measure to predict the number of (working) nephrones in the kidney 
# (higher age, less working nephrones), and this is resembled in 
# the intercept of log(creat) at t=0 of the receiver. 
# Interesting to see that donorage is not significant anymore 
# (while this is the most prominent predictor for graft failure in the literature).
# Maybe they are to much correlated? As aformentioned, donorage could be left out.

kmfit = survfit(Surv(years_tx_gl, gl_failure)~1, conf.type="log-log", data=amctx.id)

plotdf = data.frame(time=kmfit$time, n.risk=kmfit$n.risk, surv=kmfit$surv, 
                    upper=kmfit$upper, lower=kmfit$lower)

timeIndices = sapply(0:12, function(x){which.min(abs(plotdf$time-x))})

nrisktimes = c(min(plotdf$time), 1:12)
nrisk = plotdf$n.risk[timeIndices]

ggplot(data=plotdf) + geom_line(aes(x=time, y=surv)) +
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper),fill = "grey", alpha=0.6) +
  scale_y_continuous(breaks=seq(0,1,0.1),limits = c(0,1)) + 
  scale_x_continuous(breaks = seq(0,12,1), limits=c(0,12)) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0, y = 0.1, xend=12, yend = 0.1)) + 
  annotate("text",x=6, y=0.13, size=4.5, label="Number of at-risk patients") +
  annotate("text", x = nrisktimes, y = 0.05, label = paste(nrisk)) +
  theme(text = element_text(size=13), axis.text=element_text(size=13)) + 
  ylab("Survival Probability") + xlab("Time (years)")

ggsave(filename = "report/hessel/images/km.eps", width=8.27, height=9.69/2, device=cairo_ps)

modelNull = coxph(Surv(days_tx_gl, gl_failure) ~ 1,
                  data = amctx.id)

cox_All = coxph(Surv(days_tx_gl, gl_failure) ~ rec_age_fwp1 + 
                   rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
                   tx_hla + tx_pra + tx_dgf + ah_nr + tx_dm + tx_hvdis + 
                   rr_sys + rr_dias + tx_cit + d_cadaveric + 
                   rr_map + tx_dial_days ,
               data = amctx.id, x=T, model=T)

#All 3 lead to the same conclusion
backward_elim = stepAIC(cox_All, direction = "backward")
forward_selec = stepAIC(modelNull,direction="forward",scope=list(upper=cox_All,lower=modelNull))
stepwise_selec = stepAIC(modelNull,direction="both",scope=list(upper=cox_All,lower=modelNull))

###########################################################
# Using glmnet
###########################################################
cv.fit = cv.glmnet(model.matrix(cox_All), Surv(amctx.id$years_tx_gl, amctx.id$gl_failure), 
                   family = "cox", alpha=1, nfolds = 10, lambda = exp(seq(from = -10, to = 8,by = 0.01)),
                   standardize = T)
plot(cv.fit)
log(cv.fit$lambda.min)

fit = glmnet(model.matrix(cox_All), Surv(amctx.id$years_tx_gl, amctx.id$gl_failure), 
             family = "cox", alpha=1, standardize = T)

lasso_coef = coef(fit, s = cv.fit$lambda.min)
attributes(lasso_coef)$Dimnames[[1]][abs(as.matrix(lasso_coef)) > 0]

#Final cox model.
#Chosen on the basis of expert advice. besides hardly any good difference between various models
coxModel = coxph(Surv(years_tx_gl, gl_failure) ~ d_age + tx_previoustx + d_gender + rec_bmi + tx_pra + 
                   rec_age_fwp1 + I(tx_dial_days/365),
                 data = amctx.id, x=T, model=T)

amctx.id_scaled = amctx_merged_scaled[!duplicated(amctx_merged_scaled$amctx),]

coxModel_clinical = coxph(Surv(years_tx_gl, gl_failure) ~ tx_hla + tx_previoustx + 
                            tx_cit + tx_dial_days,
                 data = amctx.id_scaled, x=T, model=T)


save.image("Rdata/feedbackmeeting.Rdata")

