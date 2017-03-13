#library(parallel)

#cluster = makeCluster(5)
#clusterExport(cl=cluster, varlist=c("replaceMCMCContents", "amctx_merged", "amctx.id"))

auc_roc = function(testIds){
  
  training = amctx_merged[!(amctx_merged$amctx %in% testIds),]
  training$amctx = droplevels(training$amctx)
  training.id = amctx.id[!(amctx.id$amctx %in% testIds),]
  training.id$amctx = droplevels(training.id$amctx)
   
  test = amctx_merged[amctx_merged$amctx %in% testIds,]
  test$amctx = as.numeric(test$amctx)
  test = test[!is.na(test$creatinine),]
   
  coxModel_training = coxph(Surv(years_tx_gl, gl_failure) ~ rec_age + 
                      d_age + tx_previoustx + d_gender + 
                      rec_bmi + tx_pra + I(tx_dial_days/365),
                    data = training.id, x=T, model=T)
   
  mvglmer_training = mvglmer(list(log(creatinine) ~ rec_age_fwp1 +
                                     rec_gender + d_age +
                                     tx_pra + ah_nr + tx_dm +
                                     ns(tx_s_years,knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6)) * d_cadaveric +
                                     ns(tx_s_years,knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6)) * tx_dgf +
                                     (ns(tx_s_years, knots=c(30, 70)/365, Boundary.knots = c(0, 6))|amctx)),
                              data = training, families = list(gaussian))

  mvJoint_training = mvJointModelBayes(mvglmer_training, coxModel_training, timeVar = "tx_s_years",
                                       Formulas = list("log(creatinine)" = "value",
                                                       "log(creatinine)" = list(fixed = ~ 0 + dns(tx_s_years, knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6)) +
                                                                                  I(dns(tx_s_years, knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6))*(as.numeric(d_cadaveric)-1)) +
                                                                                  I(dns(tx_s_years, knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6))*(as.numeric(tx_dgf)-1)),
                                                                                random = ~ 0 + dns(tx_s_years, knots = c(30, 70)/365, Boundary.knots = c(0, 6)),
                                                                                indFixed = c(8:11, 14:17, 18:21), indRandom = 2:4,
                                                                                name = "slope")),
                                       priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

  jmbayes_training_replaced = replaceMCMCContents(mvJoint_training, jmbayes_creatinine_tdboth)

  timeRange = seq(1, 9, by = 0.5)
  
  auc_halfyear = vector("list", length(timeRange))
  auc_1year = vector("list", length(timeRange))
  roc_halfyear = vector("list", length(timeRange))
  roc_1year = vector("list", length(timeRange))
  k=1
  for(tstart in timeRange){
    auc_halfyear[[k]] = aucJM(jmbayes_training_replaced, test, Tstart=tstart, Dt=0.5, idVar = "amctx")
    auc_1year[[k]] = aucJM(jmbayes_training_replaced, test, Tstart=tstart, Dt=1, idVar = "amctx")
    roc_halfyear[[k]] = rocJM(jmbayes_training_replaced, test, Tstart=tstart, Dt=0.5, idVar = "amctx")
    roc_1year[[k]] = rocJM(jmbayes_training_replaced, test, Tstart=tstart, Dt=1, idVar = "amctx")
    
    k = k + 1
  }

  list(auc_halfyear, auc_1year, roc_halfyear, roc_1year)
}

################################
#Run simulations
##############################
set.seed(1000)

nrep = 30
folds = 3
n = nrow(amctx.id)

rep_splits = lapply(1:nrep, function(i){
  splits <- split(seq_len(n), sample(rep(seq_len(folds), length.out = n)))
  splits = lapply(splits, function(index){
    amctx.id$amctx[index]
  })
})

tStart = Sys.time()
# auc_roc_all = lapply(rep_splits, function(splits){
#   lapply(splits, auc_roc)
# })

auc_roc_all = vector("list", nrep)
for(i in 1:nrep){
  auc_roc_all[[i]] = vector("list", folds)
  for(j in 1:folds){
    auc_roc_all[[i]][[j]] = auc_roc(rep_splits[[i]][[j]])
  }
  save.image("Rdata/auc.Rdata")
}

tEnd = Sys.time()
save.image("Rdata/auc.Rdata")
