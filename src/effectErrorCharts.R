cbbPalette <- c("#000000","#D55E00", "#0072B2", "#E69F00", "#56B4E9", "#009E73", "#F0E442",  "#CC79A7")

rows = 200:1000

names = colnames(mvJoint_pcr_creatinine_tdboth_complex$mcmc$alphas)
names = c(c("log(PCR) Value", "log(PCR) Velocity", 
            "log(creatinine) Value", "log(creatinine) Velocity"),
          colnames(mvJoint_pcr_creatinine_tdboth_complex$mcmc$gammas))
means = c(apply(mvJoint_pcr_creatinine_tdboth_complex$mcmc$alphas[rows,], 2, mean),
          apply(mvJoint_pcr_creatinine_tdboth_complex$mcmc$gammas[rows,], 2, mean))
upper = c(apply(mvJoint_pcr_creatinine_tdboth_complex$mcmc$alphas[rows,], 2, 
              function(x){HPDinterval(as.mcmc(x))})[2,],
          apply(mvJoint_pcr_creatinine_tdboth_complex$mcmc$gammas[rows,], 2, 
                function(x){HPDinterval(as.mcmc(x))})[2,])
lower = c(apply(mvJoint_pcr_creatinine_tdboth_complex$mcmc$alphas[rows,], 2, 
              function(x){HPDinterval(as.mcmc(x))})[1,],
          apply(mvJoint_pcr_creatinine_tdboth_complex$mcmc$gammas[rows,], 2, 
                function(x){HPDinterval(as.mcmc(x))})[1,])
type = rep(c("Biomarker", "Other"), each=4)

plotdf = data.frame(names, means, upper, lower, type)

ggplot(data=plotdf) + geom_errorbar(aes(x=names, ymin=lower, ymax=upper, color=type), width=0.25) + 
  geom_point(aes(x=names, y=means, color=type)) + geom_hline(yintercept = 0, linetype="dashed") + 
  xlab("Effect name") + ylab("Mean and 95% HPDI log(HR)") +
  scale_fill_manual(values=cbbPalette) + scale_y_continuous(breaks = seq(0,10,0.5)) + 
scale_colour_manual(values=cbbPalette) +
  coord_flip() + theme(text = element_text(size=14),
                       axis.text.x = element_text(size=14)) 

##################################################
# Betas
#################################################
rows = 200:1000

names = colnames(mvJoint_pcr_creatinine_tdboth_complex$mcmc$betas1)
means = apply(mvJoint_pcr_creatinine_tdboth_complex$mcmc$betas1[rows,], 2, mean)
upper = apply(mvJoint_pcr_creatinine_tdboth_complex$mcmc$betas1[rows,], 2, 
                    function(x){HPDinterval(as.mcmc(x))})[2,]
lower = apply(mvJoint_pcr_creatinine_tdboth_complex$mcmc$betas1[rows,], 2, 
                    function(x){HPDinterval(as.mcmc(x))})[1,]

plotdf = data.frame(names, means, upper, lower)[-c(1,18:21),]

ggplot(data=plotdf) + geom_errorbar(aes(x=names, ymin=lower, ymax=upper), width=0.25) + 
  geom_point(aes(x=names, y=means)) + geom_hline(yintercept = 0, linetype="dashed") + 
  xlab("Effect name") + ylab("Mean and 95% HPDI") +
  scale_fill_manual(values=cbbPalette) + scale_y_continuous(breaks = seq(-5,10,0.1)) + 
  scale_colour_manual(values=cbbPalette) +
  coord_flip() + theme(text = element_text(size=14),
                       axis.text.x = element_text(size=14)) 

##########################################
# effect plot creatinine, and pcr
##########################################
xTimes = seq(0,8, 0.05)
nsTime = ns(xTimes,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6))
evolution = t(apply(mvJoint_pcr_creatinine_tdboth_complex$mcmc$betas1[rows,c(1,18:21)], MARGIN = 1, 
                    function(x){cbind(1, nsTime) %*% x}))

means = apply(evolution, 2, mean)
upper = apply(evolution, 2, 
              function(x){
                quantile(as.mcmc(x), c(0.025, 0.975))
                #HPDinterval(as.mcmc(x))
                })[2,]
lower = apply(evolution, 2, 
              function(x){
                quantile(as.mcmc(x), c(0.025, 0.975))
                #HPDinterval(as.mcmc(x))
                })[1,]

plotdf = data.frame(means, upper, lower)
ggplot(data=plotdf, aes(x=xTimes)) + geom_line(aes(y=means)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = "grey", alpha=0.6) +
  xlab("Time (years)") + ylab("log(PCR)") +
  #geom_vline(xintercept = c(30, 80)/365, linetype="dashed") + 
  theme(text = element_text(size=13),
        axis.text.x = element_text(size=13)) 

ggsave(filename = "report/hessel/images/pcr.eps", width=8.27, height=9.69/2, device=cairo_ps)

############################################
# Dynamic survival prob
############################################
newds = amctx_creatinine[amctx_creatinine$amctx==3,][c(1:60),]
lasttime = tail(newds$tx_s_years,1)
plot(survfitJM(mvJoint_creatinine_tdboth_complex, newds, idVar="amctx", survTimes = seq(lasttime, 5, 0.1)))
