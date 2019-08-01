FONT_SIZE = 14
POINT_SIZE = 2
THEME_COLOR = "dodgerblue4"
RISK_COLOR = "firebrick2"
INFO_COLOR = "forestgreen"

tempdata_scaled = testDs[testDs$amctx==198 & testDs$tx_s_years<=3,]
curVisit = max(tempdata_scaled$tx_s_years)

plot_problem_explanation=ggplot() + 
  geom_ribbon(aes(x=seq(curVisit,5,length.out = 4), ymin=-Inf, ymax=Inf, 
                  fill="Region with risk of graft failure"), alpha=0.2) + 
  geom_vline(xintercept = curVisit, linetype="dashed") +
  geom_vline(aes(xintercept = curVisit + (5-curVisit)/2), 
               linetype="dotted", color=THEME_COLOR, size=0.6) +
  geom_point(data=tempdata_scaled,
             aes(x=tx_s_years, y=creatinine), size=POINT_SIZE, color = THEME_COLOR) +
  geom_line(data=tempdata_scaled,
             aes(x=tx_s_years, y=creatinine), color = THEME_COLOR, alpha=0.2) + 
  theme_bw() + theme(text = element_text(size=FONT_SIZE), legend.position = "bottom",
                     legend.title = element_blank(),
                     axis.text.y = element_text(color=THEME_COLOR), 
                     axis.title.y = element_text(color=THEME_COLOR)) + 
  geom_label(aes(x=curVisit + (5-curVisit)/2, y=180, label="Optimal time of next\n creatinine measurement?"), size=4, color=THEME_COLOR)+
  scale_fill_manual(name="", values=RISK_COLOR) +
  scale_x_continuous(breaks = c(0,1.5, curVisit,curVisit + (5-curVisit)/2, 5),
                     labels = c("0","1.5","Current visit \n(3 years)",
                                "Future visit \n(u years)",
                                "\u221E"), limits = c(0,5)) +
  ylim(110, 260) +
  xlab("Time since transplantation (years)") + ylab("Serum creatinine (µmol/L)") 

ggsave(plot_problem_explanation, file="report/hessel/images/new_2018/plot_problem_explanation.eps", device=cairo_ps)


#JM explanation plot
plotJM_fitValue = ggplot(testDs[1:timesPerSubject,]) + 
  geom_point(aes(x=tx_s_years, y=rLogCreatinine(6, tx_s_years, F)), size=POINT_SIZE, color=THEME_COLOR) +
  geom_line(aes(x=tx_s_years, y=rLogCreatinine(6, tx_s_years, T)), color=THEME_COLOR) +
  theme_bw() + theme(text = element_text(size=FONT_SIZE), legend.position = "bottom",
                     legend.title = element_blank(),
                     axis.text.y = element_text(color=THEME_COLOR), 
                     axis.title.y = element_text(color=THEME_COLOR),
                     axis.text.x = element_blank(), axis.title.x = element_blank()) +
  xlim(0,10) + 
  xlab("Time since transplantation (years)") + ylab("Fitted\nserum creatinine\n levels (log µmol/L)") 

plotJM_fitVelocity = ggplot(testDs[1:timesPerSubject,]) + 
  geom_line(aes(x=tx_s_years, y=trueLogCreatinineVelocity(6, tx_s_years)), color=THEME_COLOR) +
  theme_bw() + theme(text = element_text(size=FONT_SIZE), legend.position = "bottom",
                     legend.title = element_blank(),
                     axis.text.y = element_text(color=THEME_COLOR), 
                     axis.title.y = element_text(color=THEME_COLOR),
                     axis.text.x = element_blank(), axis.title.x = element_blank()) +
  xlim(0,10) + 
  xlab("Time since transplantation (years)") + ylab("Fitted\nserum creatinine\n velocity") 

plot_hazard = ggplot() +
  geom_line(aes(x=seq(0.3, 10, 0.1), y=hazardFunc(seq(0.3, 10, 0.1), 6, testDs.id[6,],wGamma[6], b_creatinine[6,], noBaseline = T)),
            color=RISK_COLOR) +
  theme_bw() + theme(text = element_text(size=FONT_SIZE), legend.position = "bottom",
                     legend.title = element_blank(),
                     axis.text.y = element_text(color=THEME_COLOR), 
                     axis.title.y = element_text(color=THEME_COLOR)) +
  xlim(0,10) +
  geom_point(aes(x=1, y=1, shape="Observed serum creatinine\n (log µmol/L)", 
                 color="Observed serum creatinine\n (log µmol/L)"), size=POINT_SIZE) +
  geom_line(aes(x=1:2, y=1:2, linetype="Fitted serum creatinine\n (log µmol/L)"), color=THEME_COLOR) +
  scale_color_manual(name="", values=THEME_COLOR) +
  scale_shape_manual(name="", values=16) + 
  scale_linetype_manual(name="", values="solid") + 
  xlab("Time since transplantation (years)") + ylab("Hazard of \ngraft failure")

plot_survival = ggplot() +
  geom_line(aes(x=seq(0.3, 10, 0.1), y=1-sapply(seq(0.3, 10, 0.1), survivalFunc, 6, testDs.id[6,],wGamma[6], b_creatinine[6,])),
            color=RISK_COLOR) +
  theme_bw() + theme(text = element_text(size=FONT_SIZE), legend.position = "bottom",
                     legend.title = element_blank(),
                     axis.text.y = element_text(color=RISK_COLOR), 
                     axis.title.y = element_text(color=RISK_COLOR)) +
  xlim(0,10) + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                                  labels = c("0%", "25%", "50%", "75%", "100%"), limits = c(0,1)) +  
  geom_point(aes(x=1, y=2, shape="Observed serum creatinine\n (log µmol/L)", 
                 color="Observed serum creatinine\n (log µmol/L)"), size=POINT_SIZE) +
  geom_line(aes(x=1:2, y=1:2, linetype="Fitted serum creatinine\n (log µmol/L)"), color=THEME_COLOR) +
  scale_color_manual(name="", values=THEME_COLOR) +
  scale_shape_manual(name="", values=16) + 
  scale_linetype_manual(name="", values="solid") + 
  xlab("Time since transplantation (years)") + ylab("Risk of \ngraft failure (%)")


jmPlot = ggpubr::ggarrange(plotJM_fitValue, plotJM_fitVelocity, plot_survival,
                                   ncol = 1, nrow=3,align = "v",
                  common.legend = T,legend = "bottom", labels = "AUTO")

ggsave(jmPlot, file="report/hessel/images/new_2018/jmPlot.eps", device=cairo_ps, 
       width = 6, height = 7)

#STEP 1: personalized approach
plotJMFitted = plotDynamicRiskProb(198, simJointModel_replaced, origdata = testDs, ymax=5.87,
                                   maxLongTimeHorizon=3, maxFutureTime = 5, riskLabelProb = 0.05) + 
  geom_segment(aes(x=3.9, xend=3.9, y=-Inf, yend=Inf), 
               linetype="dotted", color=RISK_COLOR, size=0.65) +
  geom_segment(aes(x=3, xend=Inf, y=4.83368854033109, yend=4.83368854033109), 
               linetype="dotted", color=RISK_COLOR, size=0.65) + 
  geom_label(aes(x=3.9, y=5.1), label="Time of 5%\n cumulative risk\n (3.9 years)", color=RISK_COLOR) + 
  geom_point(aes(x=3.9, y=4.83368854033109), color=RISK_COLOR) + 
  theme(axis.text.y.left = element_text(color=THEME_COLOR), 
        axis.title.y.left = element_text(color=THEME_COLOR), 
        title = element_text(color=INFO_COLOR),
        axis.title.x = element_text(color="black")) + 
  scale_x_continuous(breaks=c(0,1.5, curVisit,3.9, 5), labels = c("0","1.5","Current visit\n(3 years)",
                                                                  "u=3.9", "5"), limits = c(0,5))

ggsave(plotJMFitted, file="report/hessel/images/new_2018/plotJMFitted.eps", device=cairo_ps)

#STEP 2: Calculate information gain at future
plotInfoGain = plotInformationGain(198, simJointModel_replaced, testDs, 3, 0.9)
optimalTime = 3.42
plotInfoGain = plotInfoGain +
  geom_segment(aes(x=optimalTime, xend=optimalTime, y=-Inf, yend=Inf), 
               linetype="dotted", color=INFO_COLOR) + 
  geom_label(aes(x=optimalTime, y=5.65), 
             label=paste0("Optimal time of\n next measurement\n(",optimalTime,"years)"), 
             color=INFO_COLOR, nudge_x = 0.125) + 
  geom_point(aes(x=optimalTime, y=5.55275544483135), color=INFO_COLOR, size=POINT_SIZE) + 
  theme(axis.text.y.left = element_text(color=THEME_COLOR), 
        axis.title.y.left = element_text(color=THEME_COLOR)) + 
  scale_x_continuous(breaks=c(0,1.5, curVisit, 3.9), 
                     labels = c("0","1.5","Current visit\n(3 years)", "u=3.9"))

ggsave(plotInfoGain, file="report/hessel/images/new_2018/plotInfoGain.eps", device=cairo_ps)


#STEP 3: take a new measurement and calculate risk profile again
newRow = tempdata_scaled[1, ]
newRow$tx_s_years = optimalTime
newRow$logCreatinine = 5.87
newRow$creatinine = exp(newRow$logCreatinine)
tempdata_scaled = rbind(tempdata_scaled, newRow)

1 - survfitJM(simJointModel_replaced, tempdata_scaled, idVar = "amctx", survTimes = optimalTime + 0.5)$summaries[[1]][,"Mean"]

plotJMFitted_newmeasurement = plotDynamicRiskProb(198, simJointModel_replaced, origdata = tempdata_scaled, ymax=5.87,
                                   maxLongTimeHorizon=optimalTime, maxFutureTime = 5, riskLabelProb = 0.05) + 
  geom_segment(aes(x=3.82, xend=3.82, y=-Inf, yend=Inf), 
               linetype="dotted", color=RISK_COLOR, size=0.65) +
  geom_segment(aes(x=3.42, xend=Inf, y=4.84955076808953, yend=4.84955076808953), 
               linetype="dotted", color=RISK_COLOR, size=0.65) + 
  geom_label(aes(x=3.82, y=5.1), label="Time of 5%\n cumulative risk\n (3.8 years)", color=RISK_COLOR) + 
  geom_point(aes(x=3.82, y=4.84955076808953), color=RISK_COLOR) + 
  
  geom_point(aes(x=optimalTime, y=newRow$logCreatinine), color=THEME_COLOR, shape=8, size=POINT_SIZE+1)+
  theme(axis.text.y.left = element_text(color=THEME_COLOR), 
        axis.title.y.left = element_text(color=THEME_COLOR),
        title = element_text(color=RISK_COLOR), 
        axis.title.x = element_text(color="black")) + 
  scale_x_continuous(breaks=c(0,1.5, optimalTime,4, 5),
                     labels = c("0","1.5","New visit\n(3.42 years)","4","5"),
                     limits = c(0,5))

ggsave(plotJMFitted_newmeasurement, file="report/hessel/images/new_2018/plotJMFitted_newmeasurement.eps", 
       device=cairo_ps)


dynSurvDemo = ggpubr::ggarrange(plotJMFitted + xlab("") + ggtitle("Time of 5% risk (3.9 years) is after 6 months of latest visit"), 
                  plotJMFitted_newmeasurement + ggtitle("Updated time of 5% risk (3.8 years) is within 6 months of latest visit"), ncol = 1, nrow=2,
                  labels = "AUTO")
ggsave(dynSurvDemo, file="report/hessel/images/new_2018/dynSurvDemo.eps", 
       device=cairo_ps, height = 9)
