source("src/R/common.R")

#################################
# Just to check the schedule
##################################
ggplot(data=amctx_creatinine[amctx_creatinine$tx_s_years<1,], aes(factor(visit_num), tx_s_days)) + 
  geom_boxplot() + stat_summary(fun.data = function(x){
    return(c(y = 1.2, label = length(x))) 
  }, geom = "text", fun.y = median)

#################################
# longitudinal analysis for creatinine
#################################

# First check the trend overall
ds = amctx_creatinine[amctx_creatinine$tx_s_days>0 & amctx_creatinine$tx_s_days<=200,]
ds$amctx = droplevels(ds$amctx)
ggplot(data=ds, aes(x=tx_s_years,y=log(value))) + geom_line(aes(group=amctx)) + stat_smooth() + 
  facet_grid(tx_dgf~d_cadaveric, labeller = label_both)

ggplot(data=amctx_creatinine, aes(x=tx_s_years,y=value)) + geom_point() + stat_smooth()

#For d_cadaveric or d_type difference evolution
ggplot(data=amctx_creatinine, aes(x=tx_s_years,y=log(value))) + 
  geom_line(aes(group=amctx)) + stat_smooth() + ticksX(0, 10, 1) + ylab("log(serum creatinine)") + 
  xlab("Follow up time (years)") + facet_grid(.~d_cadaveric, labeller = label_both)

#Tiny bit of difference for males and females in terms of where to put the knots
ggplot(data=amctx_creatinine[amctx_creatinine$tx_s_years<=1,], aes(x=tx_s_years,y=log(value))) + 
  geom_line(aes(group=amctx)) + stat_smooth() + ticksX(0, 1, 0.05) + ylab("log(serum creatinine)") + 
  xlab("tx_s_years") + facet_grid(.~rec_gender)

ggplot(data=amctx_creatinine, aes(x=tx_s_years,y=log(value))) + 
  geom_line(aes(group=amctx)) + stat_smooth() + ticksX(0, 10, 1) + ylab("log(serum creatinine)") + 
  xlab("Follow up time (years)") + facet_grid(.~tx_dgf, labeller = label_both)

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

# ah_ace + ah_arb = ah_raasi, perhaps its better to not include all three in the model. 
# Preferably no medication use is included, since they have not been collected in a dynamic way 
# (just baseline) and without doses. 
# tx_hla can be modelled as a continuous parameter
# (hessel: maybe not statistically correct, but this is often done also to reduce parameters and error)  
# 30 to 60 days for people to change trend
basic_model = lme(data=amctx_creatinine,
                  fixed=log(creatinine) ~ ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6)),
                  random = ~ns(tx_s_years,knots=c(50, 100)/365, Boundary.knots = c(0.03917808, 6))|amctx,
                  control = lmeControl(opt = "optim"), method="ML")

main_effect_complex_model = lme(data=amctx_creatinine,
                        fixed=log(creatinine) ~ rec_age_fwp1 + 
                          rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
                          tx_hla + tx_pra + tx_dgf + ah_nr + tx_dm + tx_hvdis + 
                          rr_sys + rr_dias + tx_cit + 
                          rr_map + tx_dial_days + d_cadaveric +
                          ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6)),
                        random = ~ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))|amctx,
                        control = lmeControl(opt = "optim"), method="ML")

main_effect_complex_model2 = lme(data=amctx_creatinine,
                                fixed=log(creatinine) ~ rec_age_fwp1 + 
                                  rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
                                  tx_hla + tx_pra + tx_dgf + ah_nr + tx_dm + tx_hvdis + 
                                  rr_sys + rr_dias + tx_cit + 
                                  rr_map + tx_dial_days + d_cadaveric +
                                  ns(tx_s_years,knots=c(30, 70, 900)/365, Boundary.knots = c(0.03917808, 6)),
                                random = ~ns(tx_s_years,knots=c(30, 70, 900)/365, Boundary.knots = c(0.03917808, 6))|amctx,
                                control = lmeControl(opt = "optim"), method="ML")

main_effect_complex_model3 = lme(data=amctx_creatinine,
                                 fixed=log(creatinine) ~ rec_age_fwp1 + 
                                   rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
                                   tx_hla + tx_pra + tx_dgf + ah_nr + tx_dm + tx_hvdis + 
                                   rr_sys + rr_dias + tx_cit + 
                                   rr_map + tx_dial_days + d_cadaveric +
                                   ns(tx_s_years,knots=c(30, 80, 900)/365, Boundary.knots = c(0.03917808, 6)),
                                 random = ~ns(tx_s_years,knots=c(30, 80, 900)/365, Boundary.knots = c(0.03917808, 6))|amctx,
                                 control = lmeControl(opt = "optim"), method="ML")

main_effect_complex_model4 = lme(data=amctx_creatinine,
                                 fixed=log(creatinine) ~ rec_age_fwp1 + 
                                   rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
                                   tx_hla + tx_pra + tx_dgf + ah_nr + tx_dm + tx_hvdis + 
                                   tx_cit + 
                                   tx_dial_days + d_cadaveric +
                                   ns(tx_s_years,knots=c(30, 80, 800)/365, Boundary.knots = c(0.03917808, 6)),
                                 random = ~ns(tx_s_years,knots=c(30, 80, 800)/365, Boundary.knots = c(0.03917808, 6))|amctx,
                                 control = lmeControl(opt = "optim"), method="ML")

main_effect_complex_model5 = lme(data=amctx_creatinine,
                                 fixed=log(creatinine) ~ rec_age_fwp1 + 
                                   rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
                                   tx_hla + tx_pra + tx_dgf + ah_nr + tx_dm + tx_hvdis + 
                                   tx_cit + 
                                   tx_dial_days + d_cadaveric +
                                   ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                                 random = ~ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6))|amctx,
                                 control = lmeControl(opt = "optim"), method="ML")


backward_selec = stepAIC(full_model, direction = "backward")

forward_selec = stepAIC(basic_model, direction = "forward", scope=log(creatinine) ~ rec_age_fwp1 + 
                          rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
                          tx_hla + tx_pra + tx_dgf + ah_nr + tx_dm + tx_hvdis + 
                          rr_sys + rr_dias + tx_cit + 
                          rr_map + tx_dial_days + d_cadaveric +
                          rec_age_fwp1*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          rec_gender*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          tx_previoustx*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          d_age*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          d_gender*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          d_bmi*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          rec_bmi*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          tx_hla*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          tx_pra*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          tx_dgf*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          ah_nr*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          tx_dm*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          tx_hvdis*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          rr_sys*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          rr_dias*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          tx_cit*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          rr_map*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          tx_dial_days*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          d_cadaveric*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6)))

both_selec = stepAIC(basic_model, direction = "both", scope=log(creatinine) ~ rec_age_fwp1 + 
                          rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
                          tx_hla + tx_pra + tx_dgf + ah_nr + tx_dm + tx_hvdis + 
                          rr_sys + rr_dias + tx_cit + 
                          rr_map + tx_dial_days + d_cadaveric +
                          rec_age_fwp1*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          rec_gender*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          tx_previoustx*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          d_age*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          d_gender*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          d_bmi*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          rec_bmi*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          tx_hla*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          tx_pra*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          tx_dgf*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          ah_nr*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          tx_dm*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          tx_hvdis*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          rr_sys*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          rr_dias*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          tx_cit*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          rr_map*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          tx_dial_days*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6))+
                          d_cadaveric*ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0.03917808, 6)))


#Despite the fact that BIC advises against interaction of time with tx_dgf, the graphs make it clear that such an interaction is possible
lme_creatinine_feedback2 = lme(data=amctx_merged[!is.na(amctx_merged$creatinine),],
                               fixed=log(creatinine) ~ rec_age_fwp1 + 
                                   rec_gender + d_age + tx_dgf + 
                                   tx_pra + ah_nr + tx_dm + d_cadaveric + 
                                   ns(tx_s_years,knots=c(50, 100, 900)/365) * d_cadaveric + 
                                   ns(tx_s_years,knots=c(50, 100, 900)/365) * tx_dgf,
                                 random = ~ns(tx_s_years,knots=c(50, 100)/365)|amctx,
                                 control = lmeControl(opt = "optim"), method="ML")

fixedSplines = lapply(1:50, function(i){
  second = runif(n=1, 50, 365)
  third = runif(n=1, second, 900)
  c(50, second, third)/365
})

model_splines = foreach(i=1:50, .packages = c("splines", "nlme", "ggplot2"), 
                        .export = c("amctx_merged")) %dopar%{
  fitUnivariteCreatinineModel(fixedSplineKnots = fixedSplines[[i]],
                               randomSplineKnots = fixedSplines[[i]][1:2],
                               boundaryKnots = c(0, 6))
                        }

lapply(model_splines, BIC)

model0 = fitUnivariteCreatinineModel(boundaryKnots = c(0,4))
model1 = fitUnivariteCreatinineModel(boundaryKnots = c(0,6), 
                                     fixedSplineKnots = c(50, 100, 1000)/365, 
                                     randomSplineKnots = c(50, 100)/365)

model2 = fitUnivariteCreatinineModel(boundaryKnots = c(0,6), 
                                    fixedSplineKnots = c(50, 80, 1000)/365, 
                                    randomSplineKnots = c(50, 80)/365)
model3 = fitUnivariteCreatinineModel(boundaryKnots = c(0,6), 
                                     fixedSplineKnots = c(30, 70, 1000)/365, 
                                     randomSplineKnots = c(30, 70)/365)

#After trying various models from the above list, the best fit is by model3
lme_creatinine_final = model3
save.image("Rdata/feedbackmeeting.Rdata")

plotCreatinineFittedCurve(list(model1, model2, model3), 
                          individually = F, transform = T)
