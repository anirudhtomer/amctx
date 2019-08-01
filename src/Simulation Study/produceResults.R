startCol = which(colnames(testDs.id)=="nObs_fixed_pt05")
testDs.id_long=reshape(testDs.id, direction='long', idvar='amctx', timevar = "methodNumber",
                          varying=list(seq(startCol, ncol(testDs.id), 2), 
                                       seq(startCol+1, ncol(testDs.id), 2)),
                          v.names=c('nObs', 'stopTime'))
testDs.id_long = testDs.id_long[order(testDs.id_long$amctx, testDs.id_long$methodNumber, na.last = T), ]
testDs.id_long$offset_fail = testDs.id_long$years_tx_gl - testDs.id_long$stopTime

testDs.id_long$methodName = factor(testDs.id_long$methodNumber, 
                                      labels =  c("Fixed_pt05", "Personalized_pt05", "Fixed_pt025", "Personalized_pt025"))

FONT_SIZE = 14

pt05_nObs_glT = ggplot(data=testDs.id_long[testDs.id_long$methodName %in% c("Fixed_pt05", "Personalized_pt05") & testDs.id_long$gl_failure==T,]) + 
  geom_boxplot(aes(methodName, nObs), outlier.shape = NA) + scale_y_continuous(breaks = seq(0, 70, 15), limits = c(0, 70)) + 
  scale_x_discrete(labels=c("Fixed", "Personalized")) +
  xlab("Schedule") + ylab("Nr. of creatinine measurements") + theme_bw() + 
  theme(text = element_text(size = FONT_SIZE), plot.title = element_text(size=FONT_SIZE)) + coord_flip() +
  ggtitle("For (54%) patients who observe \ngraft failure")

pt05_offset_glT = ggplot(data=testDs.id_long[testDs.id_long$methodName %in% c("Fixed_pt05", "Personalized_pt05") & testDs.id_long$gl_failure==T,]) + 
  geom_boxplot(aes(methodName, offset_fail), outlier.shape = NA) + 
  scale_y_continuous(breaks = c(0,seq(-10, 8, 4)), limits = c(-10, 8)) + 
  scale_x_discrete(labels=c("Fixed", "Personalized")) + coord_flip() +
  xlab("Schedule") + ylab("Time available for \nproactive treatment (years)") + theme_bw() + 
  theme(text = element_text(size = FONT_SIZE), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank()) 

pt05_nObs_glF = ggplot(data=testDs.id_long[testDs.id_long$methodName %in% c("Fixed_pt05", "Personalized_pt05") & testDs.id_long$gl_failure==F,]) + 
  geom_boxplot(aes(methodName, nObs), outlier.shape = NA) + scale_y_continuous(breaks = seq(0, 70, 15), limits = c(0, 70)) + 
  scale_x_discrete(labels=c("Fixed", "Personalized")) +
  xlab("Schedule") + ylab("Nr. of creatinine measurements") + theme_bw() + 
  theme(text = element_text(size = FONT_SIZE), plot.title = element_text(size=FONT_SIZE)) + coord_flip() +
  ggtitle("For (46%) patients who are \nright censored")

pt05_offset_glF = ggplot() + 
  geom_text(aes(x=1, y=-1), size=4,
            label="These patients are right censored\n at the 10 year follow-up period mark. \n \nTime available for proactive treatment \ncannot be calculated due to right censoring.") + 
  scale_y_continuous(breaks = c(0,seq(-10, 8, 4)), limits = c(-10, 8)) + 
  scale_x_discrete(labels=c("Fixed", "Personalized")) +
  xlab("Schedule") + ylab("Time available for \nproactive treatment (years)") + theme_bw() + 
  theme(text = element_text(size = FONT_SIZE), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank())+ coord_flip()

combined_plot_pt05 = ggpubr::ggarrange(ggpubr::ggarrange(pt05_nObs_glT, pt05_offset_glT, widths = c(1.275,1), align = "h"),
                                       ggpubr::ggarrange(pt05_nObs_glF, pt05_offset_glF, widths = c(1.275,1), align = "h"),
                                       ncol = 1, nrow=2, labels = "AUTO")
ggsave(combined_plot_pt05, filename = "report/hessel/images/new_2018/combined_plot_pt05.eps",
       width=8, height=7)

#######################

pt025_nObs_glT = ggplot(data=testDs.id_long[testDs.id_long$methodName %in% c("Fixed_pt025", "Personalized_pt025") & testDs.id_long$gl_failure==T,]) + 
  geom_boxplot(aes(methodName, nObs), outlier.shape = NA) + scale_y_continuous(breaks = seq(0, 70, 15), limits = c(0, 70)) + 
  scale_x_discrete(labels=c("Fixed", "Personalized")) +
  xlab("Schedule") + ylab("Nr. of creatinine measurements") + theme_bw() + 
  theme(text = element_text(size = FONT_SIZE), plot.title = element_text(size=FONT_SIZE)) + coord_flip() +
  ggtitle("For (54%) patients who observe \ngraft failure")

pt025_offset_glT = ggplot(data=testDs.id_long[testDs.id_long$methodName %in% c("Fixed_pt025", "Personalized_pt025") & testDs.id_long$gl_failure==T,]) + 
  geom_boxplot(aes(methodName, offset_fail), outlier.shape = NA) + 
  scale_y_continuous(breaks = c(0,seq(-6, 10, 4)), limits = c(-6, 10)) + 
  scale_x_discrete(labels=c("Fixed", "Personalized")) + coord_flip() +
  xlab("Schedule") + ylab("Time available for \nproactive treatment (years)") + theme_bw() + 
  theme(text = element_text(size = FONT_SIZE), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank()) 

pt025_nObs_glF = ggplot(data=testDs.id_long[testDs.id_long$methodName %in% c("Fixed_pt025", "Personalized_pt025") & testDs.id_long$gl_failure==F,]) + 
  geom_boxplot(aes(methodName, nObs), outlier.shape = NA) + scale_y_continuous(breaks = seq(0, 70, 15), limits = c(0, 70)) + 
  scale_x_discrete(labels=c("Fixed", "Personalized")) +
  xlab("Schedule") + ylab("Nr. of creatinine measurements") + theme_bw() + 
  theme(text = element_text(size = FONT_SIZE), plot.title = element_text(size=FONT_SIZE)) + coord_flip() +
  ggtitle("For (46%) patients who are \nright censored")

pt025_offset_glF = ggplot() + 
  geom_text(aes(x=1, y=2), size=4,
            label="These patients are right censored\n at the 10 year follow-up period mark. \n \nTime available for proactive treatment \ncannot be calculated due to right censoring.") + 
  scale_y_continuous(breaks = c(0,seq(-6, 10, 4)), limits = c(-6, 10)) + 
  scale_x_discrete(labels=c("Fixed", "Personalized")) +
  xlab("Schedule") + ylab("Time available for \nproactive treatment (years)") + theme_bw() + 
  theme(text = element_text(size = FONT_SIZE), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank())+ coord_flip()

combined_plot_pt025 = ggpubr::ggarrange(ggpubr::ggarrange(pt025_nObs_glT, pt025_offset_glT, widths = c(1.275,1), align = "h"),
                                       ggpubr::ggarrange(pt025_nObs_glF, pt025_offset_glF, widths = c(1.275,1), align = "h"),
                                       ncol = 1, nrow=2, labels = "AUTO")
ggsave(combined_plot_pt025, filename = "report/hessel/images/new_2018/combined_plot_pt025.eps",
       width=8, height=7)
