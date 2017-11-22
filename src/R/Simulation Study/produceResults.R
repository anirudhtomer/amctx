#############################################
#Load the results from Rdata files, 6 months 5% risk
#############################################
simTestDs.id$diff_fixed_new = simTestDs.id$stopTime_fixed_new - simTestDs.id$stoptime_True

i = 1
for(rdataname in c("pt05risk/u1_dynInfo_mod_6mo_nFix_6_k50.Rdata", 
                   "pt05risk/u1_dynInfo_mod_6mo_nFix_10_k50.Rdata",
                   "pt05risk/u1_dynInfo_mod_6mo_nFix_15_k50.Rdata",
                   "pt05risk/u1_dynInfoPar_6mo_nFix_6_k100.Rdata",
                   "pt05risk/u1_dynInfoPar_6mo_nFix_10_k100.Rdata",
                   "pt05risk/u1_dynInfoPar_6mo_nFix_15_k50.Rdata")){
  load(paste("Rdata/old/", rdataname, sep=""))
  
  simTestDs.id[,paste("nObs_",i, sep="")] = sapply(patientDsList, nrow)
  simTestDs.id[,paste("stopTime_",i, sep="")] = sapply(patientDsList, function(x){max(x$tx_s_years)})
  simTestDs.id[,paste("diff_",i, sep="")] = simTestDs.id[,paste("stopTime_",i, sep="")] - simTestDs.id$stoptime_True 
  
  i = i + 1
}


startCol = which(colnames(simTestDs.id)=="nObs_fixed_new")
simTestDs.id_long=reshape(simTestDs.id, direction='long', idvar='amctx', timevar = "methodNumber",
                          varying=list(seq(startCol, ncol(simTestDs.id), 3), 
                                       seq(startCol+1, ncol(simTestDs.id), 3),
                                       seq(startCol+2, ncol(simTestDs.id), 3)),
                          v.names=c('nObs', 'stopTime', 'offset'))
simTestDs.id_long = simTestDs.id_long[order(simTestDs.id_long$amctx, simTestDs.id_long$methodNumber, na.last = T), ]
simTestDs.id_long$offset_fail = simTestDs.id_long$stopTime - simTestDs.id_long$years_tx_gl

simTestDs.id_long$methodName = factor(simTestDs.id_long$methodNumber, 
                                      labels =  c("Fixed", "Max-KL_6obs", "Max-KL_10obs",
                                                  "Max-KL_15obs", "Max_Entropy_6obs", 
                                                  "Max_Entropy_10obs", "Max_Entropy_15obs"))
################# In depth comparison Fixed vs Pers
betterPers_2Obs = simTestDs.id$amctx[simTestDs.id$nObs_fixed_new > simTestDs.id$nObs_1]
betterPers_8Obs = simTestDs.id$amctx[simTestDs.id$nObs_fixed_new > simTestDs.id$nObs_2]
better_Fixed_2Obs = simTestDs.id$amctx[simTestDs.id$nObs_fixed_new <= simTestDs.id$nObs_1]
better_Fixed_8Obs = simTestDs.id$amctx[simTestDs.id$nObs_fixed_new <= simTestDs.id$nObs_2]

betterPers_2offset = simTestDs.id$amctx[simTestDs.id$diff_fixed_new > simTestDs.id$diff_1]
betterPers_8offset= simTestDs.id$amctx[simTestDs.id$diff_fixed_new > simTestDs.id$diff_2]
better_Fixed_2offset = simTestDs.id$amctx[simTestDs.id$diff_fixed_new <= simTestDs.id$diff_1]
better_Fixed_8offset = simTestDs.id$amctx[simTestDs.id$diff_fixed_new <= simTestDs.id$diff_2]

temp = simTestDs.id_long
temp = simTestDs.id_long[simTestDs.id_long$amctx %in% betterPers_8Obs,]

summaryDs = data.frame(matrix(nrow=max(temp$methodNumber), ncol=4))
rownames(summaryDs) = levels(temp$methodName)
colnames(summaryDs) = c("methodName", "nObs", "abs_offset", "abs_offset_fail")
summaryDs$methodName = levels(temp$methodName)
summaryDs$nObs = c(by(temp$nObs, temp$methodName, median, na.rm=T))
summaryDs$abs_offset = c(by(abs(temp$offset), temp$methodName, median, na.rm=T))
summaryDs$abs_offset_fail = c(by(abs(temp$offset_fail), temp$methodName, median, na.rm=T))

#ggplot(data=simTestDs.id[simTestDs.id$amctx %in% better_Fixed_2offset,]) + 
#  geom_boxplot(aes(x="", y=stoptime_True*12))
  
ggplot(data=summaryDs, aes(nObs, abs_offset * 12, label = methodName)) + 
  geom_label(size=4.5, nudge_y = -0.25) + geom_point() + 
  xlab("Median Number of observations") + ylab("Median |Stop Time - True threshold time| (months)") +
  scale_x_continuous(breaks = seq(0, ceiling(max(summaryDs$nObs)) + 5, 5), limits = c(0,ceiling(max(summaryDs$nObs)) + 5)) +
  scale_y_continuous(breaks = seq(0, max(summaryDs$abs_offset)*12 + 1, 1), limits = c(0,max(summaryDs$abs_offset)*12 + 1)) + 
  theme(text = element_text(size=12), axis.text=element_text(size=12))
  
ggplot(data=summaryDs, aes(nObs, abs_offset_fail * 12, label = methodName)) + 
  geom_label(size=4.5, nudge_y = -1) + geom_point() +
  xlab("Median Number of observations") + ylab("Median |Stop Time - True failure time| (months)") +
  scale_x_continuous(breaks = seq(0, ceiling(max(summaryDs$nObs)) + 5, 5), limits = c(0,ceiling(max(summaryDs$nObs)) + 5)) +
  scale_y_continuous(breaks = seq(0, max(summaryDs$abs_offset_fail)*12 + 1, 5), limits = c(0, max(summaryDs$abs_offset_fail)*12 + 1)) + 
  theme(text = element_text(size=12), axis.text=element_text(size=12))

a = ggplot(data=temp) + 
  geom_boxplot(aes(reorder(methodName, offset_fail, FUN=median), offset_fail*12)) + 
  xlab("Method") + ylab("Stop time - True fail time; [months]") + 
  theme(text = element_text(size=12), axis.text=element_text(size=12)) 
plot(a)

b = ggplot(data=temp) + 
  geom_boxplot(aes(reorder(methodName, offset, FUN=median), offset*12)) + 
  scale_y_continuous(breaks = seq(-100, 40, 10), limits = c(-100, 40)) + 
  xlab("Method") + ylab("Stop time - True threshold time; [months]") + 
  theme(text = element_text(size=12), axis.text=element_text(size=12)) 
plot(b)
  
c = ggplot(data=temp) + 
  geom_boxplot(aes(reorder(methodName, nObs, FUN=median), nObs)) + 
  scale_y_continuous(breaks=seq(0, max(temp$nObs), by = 10), limits = c(0,max(temp$nObs))) +
  xlab("Method") + ylab("Number of observations") + 
  theme(text = element_text(size=12), axis.text=element_text(size=12))
plot(c)

multiplot(a,b,c)

