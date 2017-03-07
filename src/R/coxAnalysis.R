source("src/R/common.R")

library(MASS)
library(glmnet)

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
survminer::ggsurvplot(kmfit, risk.table = T,break.time.by = 1, 
                      xlab = "Time(years)", ylim = c(0.5,1), 
                      conf.int = T)


kmfit = survfit(Surv(years_tx_gl, gl_failure)~I(rec_age < tx_dial_days), conf.type="log-log", data=amctx.id)
survminer::ggsurvplot(kmfit, risk.table = T,break.time.by = 1, 
                      xlab = "Time(years)", ylim = c(0.5,1), 
                      conf.int = T)

modelNull = coxph(Surv(days_tx_gl, gl_failure) ~ 1,
                  data = amctx.id)

cox_All = coxph(Surv(days_tx_gl, gl_failure) ~ rec_age + 
                   rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
                   tx_hla + tx_pra + tx_dgf + ah_nr + tx_dm + tx_hvdis + 
                   rr_sys + rr_dias +
                   rr_map + tx_dial_days + d_cadaveric,
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
cv.fit$lambda.min

fit = glmnet(model.matrix(cox_All), Surv(amctx.id$years_tx_gl, amctx.id$gl_failure), 
             family = "cox", alpha=1, standardize = T)

lasso_coef = coef(fit, s = cv.fit$lambda.min)
attributes(lasso_coef)$Dimnames[[1]][abs(as.matrix(lasso_coef)) > 0]

#Final cox model.
#rec_age and tx_dial_days are chosen on the basis of KM curves
coxModel = coxph(Surv(years_tx_gl, gl_failure) ~ rec_age + 
                   d_age + tx_previoustx + d_gender + 
                   rec_bmi + tx_pra + I(tx_dial_days/365),
                 data = amctx.id, x=T, model=T)
save.image("Rdata/feedbackmeeting.Rdata")
