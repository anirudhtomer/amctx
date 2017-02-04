library(survival)
library(MASS)
library(glmnet)

colsOfInterest = c("rec_age", "rec_gender",
                   "tx_previoustx", "d_age", "d_gender", "d_bmi", "rec_bmi",
                   "tx_hla", "tx_pra", "tx_dgf", "tx_cit", "is_nr", "is_pred", 
                   "is_aza", "is_cni", "is_mmf", "is_mtor", "ah_nr", "ah_diur",      
                   "ah_ace", "ah_arb", "ah_raasi", "ah_bb", "ah_ccb", "statin", 
                   "dm_oad", "dm_insulin", "tx_dm","tx_hvdis", "rr_sys",
                   "rr_dias", "rr_map", "tx_dial_days", "d_type", "d_cadaveric")

# new suggestion: tx_hla tx_pra tx_dgf rec_age 
# old suggestion: d_age d_bmi d_type tx_cit tx_hla rec_age tx_previoustx tx_dial_days tx_dm rec_bmi tx_pra
# avoid as they have almost empty categories: is_pred is_aza is_mtor is_cni statin
# avoid as they have almost empty categories in 2x2 table against failure: dm_oad is_mmf ah_arb ah_ace 

zeroVar = c("statin", "is_pred", "is_mtor")
#What to do with ah_nr(1 col with 0), dm_insulin(infinite beta for no reason), d_cadaveric(singular)

modelNull = coxph(Surv(days_tx_gl, gl_failure) ~ 1,
                  data = amctx.id)

modelAll = coxph(Surv(days_tx_gl, gl_failure) ~ rec_age  + rec_gender +
                 tx_previoustx + d_age + d_gender + d_bmi +
                 rec_bmi + tx_hla + tx_pra + tx_dgf + 
                 tx_cit+ is_nr  +  
                 ah_diur +   
                 ah_raasi + ah_bb + ah_ccb  + 
                dm_insulin + tx_dm + tx_hvdis + 
                 rr_sys + rr_dias + rr_map + 
                 tx_dial_days + d_type,
               data = amctx.id, x=T, model=T)

colnames(model.matrix(model1))

#All 3 lead to the same conclusion
backward_elim = stepAIC(modelAll, direction = "backward")
forward_selec = stepAIC(modelNull,direction="forward",scope=list(upper=modelAll,lower=modelNull))
stepwise_selec = stepAIC(modelNull,direction="both",scope=list(upper=modelAll,lower=modelNull))

###########################################################
# Using glmnet
###########################################################
cv.fit = cv.glmnet(model.matrix(modelAll), Surv(amctx.id$years_tx_gl, amctx.id$gl_failure), 
                   family = "cox", alpha=1, nfolds = 10, lambda = exp(seq(from = -10, to = 4,by = 0.01)),
                   standardize = T)
plot(cv.fit)

fit = glmnet(model.matrix(modelAll), Surv(amctx.id$years_tx_gl, amctx.id$gl_failure), 
             family = "cox", alpha=1, standardize = T)

lasso_coef = coef(fit, s = cv.fit$lambda.min)
attributes(lasso_coef)$Dimnames[[1]][abs(as.matrix(lasso_coef)) > 0]

#model without the subject 346
amctx.id = amctx.id[amctx.id$amctx != 346, ]
amctx.id$amctx = droplevels(amctx.id$amctx)

coxModel = coxph(Surv(years_tx_gl, gl_failure) ~ d_age + tx_previoustx + rec_bmi +
                   d_type + ah_diur,
                 data = amctx.id, x=T, model=T)
coxModelConfInt = summary(coxModel)$conf.int[,c(3,4)]

time.dep.zph <- cox.zph(coxModel, transform = function(x){x})
for(i in 1:ncol(model.matrix(coxModel))){
  plot(time.dep.zph[i])
  abline(h = 0, lty=3)
  abline(h=coef(coxModel)[i], col="red")
}

#Age, rec_bmi, d_type are pretty much remaining at the same level...
#for tx_previoustx I am confused what to do
#for ah_diur we have a strong effect against time...like spline

#Can't export with time-transform tt..,how will I use in joint model??

#At the end I have decided....keep what clinician says as well

coxModel = coxph(Surv(years_tx_gl, gl_failure) ~ d_age + tx_previoustx + rec_bmi +
                   d_type + ah_diur + 
                   tx_dgf + tx_hla +rec_age + tx_pra,
               data = amctx.id, x=T, model=T)

plot(survfit(Surv(years_tx_gl, gl_failure)~1, data=amctx.id))


coxSnellRes=-log(amctx.id$gl_failure-resid(coxModel,type="martingale"))
fitres=survfit(coxph(Surv(coxSnellRes,amctx.id$gl_failure)~1,method='breslow'),type='aalen')
plot(fitres$time,-log(fitres$surv),type='s',xlab='Cox-Snell Residuals', 
     ylab='Estimated Cumulative Hazard Function')
abline(0,1,col='red',lty=2)
