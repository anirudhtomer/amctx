library(ggplot2)
library(JMbayes)
library(nlme)
library(splines)

ticks = function(by){
  scale_x_continuous(breaks = round(seq(0, max(amctx_creatinine$tx_s_days), by = by),0))
}

#################################
# longitudinal analysis for creatinine
#################################

# First check the trend overall
ggplot(data=amctx_creatinine, aes(x=tx_s_days,y=value, group=amctx)) + geom_line() + 
  facet_grid(.~rec_gender)

ggplot(data=amctx_creatinine, aes(x=tx_s_days,y=value)) + geom_point() + stat_smooth() + ticks(200)

ggplot(data=amctx_creatinine[1:5000,], aes(x=tx_s_years,y=log(value))) + 
  geom_line(aes(group=amctx)) + stat_smooth() + ticks(1) + ylab("log(serum creatinine)") + 
  xlab("tx_s_years")

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

