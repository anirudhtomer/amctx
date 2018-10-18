#############################################
#Load the results from Rdata files, 6 months 5% risk
#############################################
simTestDs.id$diff_fixed_new = simTestDs.id$stopTime_fixed_new - simTestDs.id$stoptime_True

dir = "Rdata/new/pt05risk/"
i = 1

#Index is 1 for both 2.5% and 5% risk
fileNames = list.files(path = dir)[1]
for(rdataname in fileNames){
  load(paste(dir, rdataname, sep=""))
  
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
                                      labels =  c("Fixed", paste(1:length(fileNames))))
simTestDs.id_long$methodName = factor(simTestDs.id_long$methodNumber, 
                                      labels =  c("Fixed", "Personalized"))

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
  geom_boxplot(aes(methodName, offset_fail*12), outlier.shape = NA) + 
  xlab("Method") + ylab("Failure offset (months)") + 
  theme(text = element_text(size=13), axis.text=element_text(size=13)) +
  geom_hline(yintercept=0, linetype="dashed") + scale_y_continuous(breaks = seq(-200, 200, 15)) +
  coord_cartesian(ylim= c(-125,15))
plot(a)
ggsave(filename = "report/hessel/images/truestoptimept025.eps", width=8.27, height=9.69/2, device=cairo_ps)

b = ggplot(data=temp) + 
  geom_boxplot(aes(methodName, offset*12), outlier.shape = NA) + 
  scale_y_continuous(breaks = seq(-16, 12.5, 2)) + 
  coord_cartesian(ylim= c(-16,12.5)) +
  xlab("Method") + ylab("Intervention offset (months)") + 
  theme(text = element_text(size=13), axis.text=element_text(size=13)) +
  geom_hline(yintercept=0, linetype="dashed")
plot(b)
ggsave(filename = "report/hessel/images/truethrestimept025.eps", width=8.27, height=9.69/2, device=cairo_ps)

c = ggplot(data=temp) + 
  geom_boxplot(aes(methodName, nObs), outlier.shape = NA) + 
  scale_y_continuous(breaks=seq(3, max(temp$nObs), by = 3)) +
  coord_cartesian(ylim= c(3,max(temp$nObs))) +
  xlab("Schedule") + ylab("Number of measurements") + 
  theme(text = element_text(size=13), axis.text=element_text(size=13)) 
plot(c)

ggsave(filename = "report/hessel/images/nObspt025.eps", width=8.27, height=9.69/2, device=cairo_ps)

