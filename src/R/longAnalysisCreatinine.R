source("src/R/common.R")

ticks = function(from=0, to, by){
  scale_x_continuous(breaks = seq(from, to, by = by))
}

#################################
# longitudinal analysis for creatinine
#################################

# First check the trend overall
ggplot(data=amctx_creatinine, aes(x=tx_s_days,y=value, group=amctx)) + geom_line() + 
  facet_grid(.~rec_gender)

ggplot(data=amctx_creatinine, aes(x=tx_s_days,y=value)) + geom_point() + stat_smooth() + ticks(200)

#For d_cadaveric or d_type difference evolution
ggplot(data=amctx_creatinine, aes(x=tx_s_years,y=log(value))) + 
  geom_line(aes(group=amctx)) + stat_smooth() + ticks(0, 10, 0.5) + ylab("log(serum creatinine)") + 
  xlab("tx_s_years") + facet_grid(.~d_cadaveric)

#Tiny bit of difference for males and females in terms of where to put the knots
ggplot(data=amctx_creatinine[amctx_creatinine$tx_s_years<=1,], aes(x=tx_s_years,y=log(value))) + 
  geom_line(aes(group=amctx)) + stat_smooth() + ticks(0, 1, 0.05) + ylab("log(serum creatinine)") + 
  xlab("tx_s_years") + facet_grid(.~rec_gender)

#Tiny bit of difference for males and females in terms of where to put the knots
ggplot(data=amctx_creatinine, aes(x=tx_s_years,y=log(value))) + 
  geom_line(aes(group=amctx)) + stat_smooth() + ticks(0, 2, 0.5) + ylab("log(serum creatinine)") + 
  xlab("tx_s_years") + facet_grid(.~tx_dgf) + xlim(0,2)

idList = unique(amctx_creatinine$amctx)
ggplot(data=amctx_creatinine[amctx$amctx==idList[1],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[2],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[3],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[4],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[5],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[6],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[7],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[8],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[9],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[10],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[11],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[12],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[13],], aes(x=tx_s_days,y=value)) + geom_line()

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
                 tx_dial_days + d_type, data=amctx_creatinine)

amctx_creatinine$residuals = linearmodel$residuals

ggplot(data=amctx_creatinine, aes(x=tx_s_days,y=residuals)) + geom_point() + stat_smooth() + 
   ticks(200)

# This plot gives much better idea of where to put knots. females have some outliers which may influence knot selection
ggplot(data=amctx_creatinine, aes(x=tx_s_days,y=residuals)) + geom_point() + stat_smooth() + 
  facet_grid(rec_gender~.) + ticks(200)

# Put knots at 150, 400, 1100 and make an additive model
model_rand_nsslope = lme(data=amctx_creatinine, fixed=log(value)~rec_age + rec_gender +
              tx_previoustx + d_age + d_gender + d_bmi +
              rec_bmi + tx_hla + tx_pra + tx_dgf + 
              tx_cit+ is_nr + is_aza + 
              is_cni + is_mmf  + ah_nr + 
              ah_diur + ah_ace + ah_arb + 
              ah_raasi + ah_bb + ah_ccb  + 
              dm_oad + dm_insulin + tx_dm + tx_hvdis + 
              rr_sys + rr_dias + rr_map + 
              tx_dial_days + d_type + 
              ns(tx_s_years,knots=c(150, 400, 1000)/365),
            random = ~ns(tx_s_years,knots=c(150)/365)|amctx, method="ML")

anova.lme(model_rand_nsslope, type = "marginal", adjustSigma = F)

model_rand_nsslope_2 = lme(data=amctx_creatinine, fixed=log(value)~rec_age + rec_gender +
                             d_age +  tx_dgf + 
                             is_cni +
                             ns(tx_s_years,knots=c(150, 400, 1000)/365),
                           random = ~ns(tx_s_years,knots=c(150)/365)|amctx, method = "ML")
anova.lme(model_rand_nsslope_2, type = "marginal", adjustSigma = F)
anova(model_rand_nsslope, model_rand_nsslope_2)

amctx_creatinine$residuals = model_rand_nsslope_2$residuals[,2]
amctx_creatinine$fitted = model_rand_nsslope_2$fitted[,2]

ggplot(data=amctx_creatinine, aes(x=tx_s_days,y=residuals)) + geom_point() + stat_smooth() + 
  ticks(200)

ggplot(data=amctx_creatinine, aes(x=fitted,y=residuals)) + geom_point() + stat_smooth()
ggplot(data=amctx_creatinine, aes(x=fitted,y=residuals^2)) + geom_point() + stat_smooth()

#change knot positions, 100, 400, 1000
model_rand_nsslope_3 = lme(data=amctx_creatinine, fixed=log(value)~rec_age + rec_gender +
                             d_age +  tx_dgf + 
                             is_cni +
                             ns(tx_s_years,knots=c(100, 400, 1000)/365),
                           random = ~ns(tx_s_years,knots=c(100)/365)|amctx, method = "ML")
anova(model_rand_nsslope_2, model_rand_nsslope_3)

amctx_creatinine$residuals = model_rand_nsslope_3$residuals[,2]
amctx_creatinine$fitted = model_rand_nsslope_3$fitted[,2]

ggplot(data=amctx_creatinine, aes(x=tx_s_days,y=residuals)) + geom_point() + stat_smooth() + 
  ticks(200)

ggplot(data=amctx_creatinine, aes(x=fitted,y=residuals)) + geom_point() + stat_smooth()
ggplot(data=amctx_creatinine, aes(x=fitted,y=residuals^2)) + geom_point() + stat_smooth()

#change knot positions, 100, 300, 1000: fairly well supported by the residual plot from linear model
model_rand_nsslope_4 = lme(data=amctx_creatinine, fixed=log(value)~rec_age + rec_gender +
                             d_age +  tx_dgf + 
                             is_cni +
                             ns(tx_s_years,knots=c(100, 300, 1000)/365),
                           random = ~ns(tx_s_years,knots=c(100)/365)|amctx, method = "ML")
anova(model_rand_nsslope_4, model_rand_nsslope_3)

amctx_creatinine$residuals = model_rand_nsslope_4$residuals[,2]
amctx_creatinine$fitted = model_rand_nsslope_4$fitted[,2]

ggplot(data=amctx_creatinine, aes(x=tx_s_days,y=residuals)) + geom_point() + stat_smooth() + 
  ticks(150)

ggplot(data=amctx_creatinine, aes(x=fitted,y=residuals)) + geom_point() + stat_smooth()
ggplot(data=amctx_creatinine, aes(x=fitted,y=residuals^2)) + geom_point() + stat_smooth()

# can see tiny bit heterogeneity in residuals...change random effect structure
model_rand_nsslope_5 = lme(data=amctx_creatinine, fixed=log(value)~rec_age + rec_gender +
                             d_age +  tx_dgf + is_cni +
                             ns(tx_s_years,knots=c(100, 300, 1000)/365),
                           random = ~ns(tx_s_years,knots=c(100, 300)/365)|amctx, method="ML",
                           control = lmeControl(opt = "optim"))

anova(model_rand_nsslope_5, model_rand_nsslope_4)

amctx_creatinine$residuals = model_rand_nsslope_5$residuals[,2]
amctx_creatinine$fitted = model_rand_nsslope_5$fitted[,2]

ggplot(data=amctx_creatinine, aes(x=tx_s_days,y=residuals)) + geom_point() + stat_smooth() + 
  ticks(200)

ggplot(data=amctx_creatinine, aes(x=fitted,y=residuals)) + geom_point() + stat_smooth()
ggplot(data=amctx_creatinine, aes(x=fitted,y=residuals^2)) + geom_point() + stat_smooth()

#Final choice
model_creatinine = lme(data=amctx_creatinine, fixed=log(value)~rec_age_fwp1 + rec_gender +
                             d_age +  tx_dgf + is_cni + d_bmi + tx_hla + tx_previoustx + 
                         tx_pra + tx_cit + tx_dial_days + tx_dm + rec_bmi + 
                             ns(tx_s_years,knots=c(100, 300, 1000)/365),
                           random = ~ns(tx_s_years,knots=c(100, 300)/365)|amctx,
                           control = lmeControl(opt = "optim"))

# ah_ace + ah_arb = ah_raasi, perhaps it’s better to not include all three in the model. 
# Preferably no medication use is included, since they have not been collected in a dynamic way 
# (just baseline) and without doses. 
# tx_hla can be modelled as a continuous parameter
# (hessel: maybe not statistically correct, but this is often done also to reduce parameters and error)  
# 30 to 60 days for people to change trend
model_creatinine_feedback1 = lme(data=amctx_creatinine, fixed=log(value) ~ rec_age_fwp1 + 
                                   rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
                                   tx_hla + tx_pra + tx_dgf + ah_nr + tx_dm + tx_hvdis + 
                                   rr_sys + rr_dias +
                                   rr_map + tx_dial_days + d_cadaveric +
                                   ns(tx_s_years,knots=c(50, 100, 900)/365) * d_cadaveric + 
                                   ns(tx_s_years,knots=c(50, 100, 900)/365) * tx_dgf,
                                 random = ~ns(tx_s_years,knots=c(50, 100)/365)|amctx,
                                 control = lmeControl(opt = "optim"), method="ML")

#Despite the fact that BIC advises against interaction of time with tx_dgf, the graphs make it clear that such an interaction is possible
lme_creatinine_feedback2 = lme(data=amctx_merged[!is.na(amctx_merged$creatinine),], 
                               fixed=log(creatinine) ~ rec_age_fwp1 + 
                                   rec_gender + d_age + 
                                   tx_pra + ah_nr + tx_dm + 
                                   ns(tx_s_years,knots=c(50, 100, 900)/365) * d_cadaveric + 
                                   ns(tx_s_years,knots=c(50, 100, 900)/365) * tx_dgf,
                                 random = ~ns(tx_s_years,knots=c(50, 100)/365)|amctx,
                                 control = lmeControl(opt = "optim"), method="ML")

jmbayes_creatinine_orig = jointModelBayes(lmeObject = lme_creatinine_feedback2, survObject = coxModel, timeVar = "tx_s_years", control = list(n.iter=1000))
jmbayes_creatinine_replaced = replaceMCMCContents(mvJoint_creatinine_tdval, jmbayes_creatinine_orig)
