ND = amctx_merged[!is.na(amctx_merged$creatinine),]
ND$amctx = as.numeric(ND$amctx)

roc_creatinine = rocJM(jmbayes_creatinine_replaced, ND, Tstart=3.5, Dt=c(1), idVar = "amctx")

qplot(y=roc_creatinine$TP, x=roc_creatinine$FP,geom="line") + 
  geom_abline(aes(slope=1, intercept=0),linetype = "dotted") + ylab("Sensitivity") + 
  xlab("1-Specificity") + xlim(0,1) + ylim(0,1)
