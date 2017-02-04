library(ggplot2)
library(JMbayes)
library(nlme)
library(splines)
library(glmmLasso)

ticks = function(by){
  scale_x_continuous(breaks = round(seq(0, max(amctx_creatinine$tx_s_years), by = by),0))
}

#####################################################
# Longitudinal analysis for PCR
#####################################################
ggplot(data=amctx_pcr, aes(x=tx_s_days,y=value)) + geom_point() + stat_smooth() + ticks(200) + 
  ylim(0,250)

ggplot(data=amctx_pcr[1:5000,], aes(x=tx_s_years,y=log(value))) + 
  geom_line(aes(group=amctx)) + stat_smooth() + ticks(1) + ylab("log(pcr)") +xlab("tx_s_years")

#Firsly log transform is better. 
#Secondly I could see some striations on that plot at the bottom, which are marked now
ggplot(data=amctx_pcr, aes(x=tx_s_days,y=log(value), group=amctx, color=factor(log(value)))) + geom_point() + ticks(200) +
  guides(color=F) + 
  geom_abline(intercept = log(2), slope=0) +
  geom_abline(intercept = log(3), slope=0) +
  geom_abline(intercept = log(4), slope=0) +
  geom_abline(intercept = log(5), slope=0) 

#What if I add a bit of noise to get rid of striations
noise = rnorm(sum(log(amctx_pcr$value) <= 2.5), mean = 0.5, sd = 0.2)
range(exp(noise))
amctx_pcr$noise_val = amctx_pcr$value
k = 1
for(i in 1:length(amctx_pcr$noise_val)){
  if(log(amctx_pcr$noise_val[i]) <= 2.5){
    amctx_pcr$noise_val[i] = amctx_pcr$noise_val[i] + exp(noise[k])
    k = k + 1
  }    
}

ggplot(data=amctx_pcr, aes(x=tx_s_days,y=log(noise_val), group=amctx, color=factor(log(noise_val)))) + 
  geom_point() + ticks(200) + guides(color=F)
  
# take residuals and repeat what you did above
linearmodel=lm(log(value)~rec_age_fwp1 + rec_gender +
                 tx_previoustx + d_age + d_gender + d_bmi +
                 rec_bmi + tx_hla + tx_pra + tx_dgf + 
                 tx_cit+ is_nr + is_aza + 
                 is_cni + is_mmf  + ah_nr + 
                 ah_diur + ah_ace + ah_arb + 
                 ah_raasi + ah_bb + ah_ccb  + 
                 dm_oad + dm_insulin + tx_dm + tx_hvdis + 
                 rr_sys + rr_dias + rr_map + 
                 tx_dial_days + d_type, data=amctx_pcr)

amctx_pcr$residuals = linearmodel$residuals
amctx_pcr$fitted = linearmodel$fitted.values

ggplot(data=amctx_pcr, aes(x=tx_s_days,y=residuals)) + geom_point() + stat_smooth() + 
  ticks(150) + facet_grid(.~rec_gender)
#Striations gone in above plot but can be seen in below plot
ggplot(data=amctx_pcr, aes(x=fitted,y=residuals)) + geom_point() + stat_smooth()
ggplot(data=amctx_pcr, aes(x=fitted,y=residuals^2)) + geom_point() + stat_smooth()


# Put knots at 90, 450 and make an additive model
model_rand_nsslope = lme(data=amctx_pcr, fixed=log(value)~rec_age + rec_gender +
                           tx_previoustx + d_age + d_gender + d_bmi +
                           rec_bmi + tx_hla + tx_pra + tx_dgf + 
                           tx_cit+ is_nr + is_aza + 
                           is_cni + is_mmf  + ah_nr + 
                           ah_diur + ah_ace + ah_arb + 
                           ah_raasi + ah_bb + ah_ccb  + 
                           dm_oad + dm_insulin + tx_dm + tx_hvdis + 
                           rr_sys + rr_dias + rr_map + 
                           tx_dial_days + d_type + 
                           ns(tx_s_years,knots=c(90, 450)/365),
                         random = ~ns(tx_s_years,knots=c(90)/365)|amctx, method="ML")

anova.lme(model_rand_nsslope, type = "marginal", adjustSigma = F)

model_rand_nsslope_2 = lme(data=amctx_pcr, fixed=log(value)~d_age + is_aza + ah_ace + ah_arb + 
                             ah_raasi  + rec_bmi + 
                             ns(tx_s_years,knots=c(90, 450)/365),
                           random = ~ns(tx_s_years,knots=c(90)/365)|amctx, method="ML")

anova.lme(model_rand_nsslope_2, type = "marginal", adjustSigma = F)
anova.lme(model_rand_nsslope_2, model_rand_nsslope)

amctx_pcr$residuals = model_rand_nsslope_2$residuals[,2]
amctx_pcr$fitted = model_rand_nsslope_2$fitted[,2]

ggplot(data=amctx_pcr, aes(x=tx_s_days,y=residuals)) + geom_point() + stat_smooth() + 
  ticks(200)

ggplot(data=amctx_pcr, aes(x=fitted,y=residuals)) + geom_point() + stat_smooth()
ggplot(data=amctx_pcr, aes(x=fitted,y=residuals^2)) + geom_point() + stat_smooth()

#Dunno what to do with striations in residuals. Now i'll check the knots
model_rand_nsslope_3 = lme(data=amctx_pcr, fixed=log(value)~d_age + is_aza + ah_ace + ah_arb + 
                             ah_raasi  + rec_bmi + 
                             ns(tx_s_years,knots=c(100, 200, 350)/365),
                           random = ~ns(tx_s_years,knots=c(100)/365)|amctx, method="ML"
)

anova.lme(model_rand_nsslope_3, model_rand_nsslope_2)

amctx_pcr$residuals = model_rand_nsslope_3$residuals[,2]
amctx_pcr$fitted = model_rand_nsslope_3$fitted[,2]

ggplot(data=amctx_pcr, aes(x=tx_s_days,y=residuals)) + geom_point() + stat_smooth() + 
  ticks(200)

ggplot(data=amctx_pcr, aes(x=fitted,y=residuals)) + geom_point() + stat_smooth()
ggplot(data=amctx_pcr, aes(x=fitted,y=residuals^2)) + geom_point() + stat_smooth()

#######with more complicated random effects
model_rand_nsslope_4 = lme(data=amctx_pcr, fixed=log(value)~ rec_age_fwp1 + d_age + d_gender +
                             d_bmi + rec_bmi + tx_hla + tx_pra + tx_dgf + tx_cit + 
                             is_nr + is_mmf + ah_nr + ah_diur + ah_ace + ah_arb + ah_raasi + 
                             ah_bb + ah_ccb + dm_oad + dm_insulin + tx_dm + tx_hvdis + 
                             rr_sys + rr_dias + rr_map + tx_dial_days + d_type +
                             ns(tx_s_years,knots=c(100, 200, 350)/365),
                           random = ~ns(tx_s_years,knots=c(100, 200)/365)|amctx, method="ML",
                           control = lmeControl(opt = "optim"))
anova.lme(model_rand_nsslope_4, type = "marginal", adjustSigma = F)

# glmmLasso(points ~ transfer.spendings + ave.unfair.score
#           + ball.possession + tackles
#           + ave.attend + sold.out, rnd = list(team=~1),
#           lambda=10, data = soccer)

lambdaVec = seq(from=50, to = 4000, by = 200)
bicVec = vector("numeric", length(lambdaVec))
aicVec = vector("numeric", length(lambdaVec))
for(i in 1:length(lambdaVec)){
lasso_pcr = glmmLasso(fix = log(value)~ rec_age_fwp1 + d_age + d_gender +
            d_bmi + rec_bmi + tx_hla + tx_pra + tx_dgf + tx_cit + 
            is_nr + is_mmf + ah_nr + ah_diur + ah_ace + ah_arb + ah_raasi + 
            ah_bb + ah_ccb + dm_oad + dm_insulin + tx_dm + tx_hvdis + 
            rr_sys + rr_dias + rr_map + tx_dial_days + d_type +
            tx_s_years,
          rnd = list(amctx=~1 + tx_s_years), lambda=lambdaVec[i], switch.NR=T,final.re=T,
          data=amctx_pcr)
bicVec[i] = lasso_pcr$bic
aicVec[i] = lasso_pcr$aic
}
#lambda around 800 to 1000 seems ok


amctx_pcr$residuals = model_rand_nsslope_4$residuals[,2]
amctx_pcr$fitted = model_rand_nsslope_4$fitted[,2]

ggplot(data=amctx_pcr, aes(x=tx_s_days,y=residuals)) + geom_point() + stat_smooth() + 
  ticks(200)

ggplot(data=amctx_pcr[amctx_pcr$value==2,], aes(y=tx_s_days,x=fitted)) + geom_point() + stat_smooth() + 
  ticks(2)

ggplot(data=amctx_pcr[amctx_pcr$value>0,], aes(x=fitted,y=residuals)) + geom_point() + stat_smooth()
ggplot(data=amctx_pcr, aes(x=fitted,y=residuals^2)) + geom_point() + stat_smooth()

#Final choice
model_pcr = lme(data=amctx_pcr, fixed=log(value)~d_age + ah_ace + ah_arb + 
                  ah_raasi  + rec_bmi + d_type + d_bmi + tx_cit+ 
                  tx_hla+ rec_age_fwp1 + tx_previoustx + tx_dial_days + 
                  tx_dm + tx_pra + ah_nr + 
                  ns(tx_s_years,knots=c(100, 200, 350)/365),
                random = ~ns(tx_s_years,knots=c(100, 200)/365)|amctx,
                control = lmeControl(opt = "optim"))
anova.lme(model_pcr, type = "marginal", adjustSigma = F)
