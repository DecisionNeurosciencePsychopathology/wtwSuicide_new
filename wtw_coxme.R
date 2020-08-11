
setwd("~/Box Sync/Project_wtw/R")
library(readr)
library(ggplot2)
library(tidyr)
library(tibble)
library(Hmisc)
library(reshape2)
library(corrplot)
library(emmeans)
# library(factoextra)
# library(ggfortify)
library(compareGroups)
# library(RColorBrewer)
# library(MASS)
library(effects)
# library(readr)
# library(VIM)
# library(mice)
library(multcompView)
library(stargazer)
library(dplyr)
library(lme4)
library(survival)
library(coxme)
library(survminer)
# library(OIsurv)
library(ggpubr)
library(sjPlot)
library(sjlabelled)
library(sjmisc)

source('~/code/R/vif.lme.R')

load(file = "wtw.RData")

# #  for collinearity diagnostics
# vif.lme <- function (fit) {
#   ## adapted from rms::vif
#   v <- vcov(fit)
#   nam <- names(fixef(fit))
#   ## exclude intercepts
#   ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
#   if (ns > 0) {
#     v <- v[-(1:ns), -(1:ns), drop = FALSE]
#     nam <- nam[-(1:ns)] }
#   d <- diag(v)^0.5
#   v <- diag(solve(v/(d %o% d)))
#   names(v) <- nam
#   v }


# make a na.omit df for plotting predicted values
new_tddf <- na.omit(tddf[,c("ID", "latency", "quit","lat_lag", "initialTime", "win_lag", "waitValue", "Group")])

cf0 <- coxph(Surv(latency,quit)~1, new_tddf)
# plot(cf0)
c0 <- coxme(Surv(latency, quit) ~ (1|ID), new_tddf)
new_tddf$predicted_lp_c0 <- predict(c0, type = "risk")
pdf(file = 'predicted risk c0.pdf')
ggplot(new_tddf[new_tddf$latency<4,],aes(x = latency, y = predicted_lp_c0, color = Group)) + geom_smooth() 
# + facet_wrap(~ID)
dev.off()

c0v <- coxme(Surv(latency, quit) ~ scale(waitValue) + (1|ID), tddf)
summary(c0v)


# build model systematically
c1 <- coxme(Surv(latency, quit) ~ scale(lat_lag) + (1|ID), new_tddf)
new_tddf$predicted_h_c1 <- predict(c1, type = "risk")
p1 <- ggplot(new_tddf[new_tddf$latency<4,],aes(x = latency, y = predicted_h_c1, color = Group)) + geom_smooth() 
new_tddf$predicted_logh_c1 <- predict(c1, type = "lp")
p2 <- ggplot(new_tddf[new_tddf$latency<4,],aes(x = latency, y = predicted_logh_c1, color = Group)) + geom_smooth() 

pdf(file = 'predicted c1.pdf')
ggarrange(p1,p2,nrow = 2)# + facet_wrap(~ID)
dev.off()


summary(c1)
c2 <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + (1|ID), tddf)
summary(c2)
c3 <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + (1|ID), tddf)
summary(c3)
vif.lme(c3)
c4 <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(waitValue) + (1|ID), tddf)
summary(c4)
vif.lme(c4)
################
# more complex models fit better, but may not add much to interpretation beyond this
################
c4g <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(waitValue) * Group + (1|ID), tddf)
summary(c4g)
vif.lme(c4g)

c4g_new <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(waitValue) * Group + (1|ID), new_tddf)
summary(c4g_new)

new_tddf$predicted_lp_c4g <- predict(c4g_new, type = "risk")
pdf(file = 'predicted risk by group.pdf')
ggplot(new_tddf[new_tddf$latency<4,],aes(x = latency, y = predicted_lp_c4g, color = Group)) + geom_smooth() + facet_wrap(~win_lag)
dev.off()

# remove lagged latency
c4g_new_nolag <- coxme(Surv(latency, quit) ~ scale(initialTime) + win_lag + scale(waitValue) * Group + (1|ID), new_tddf)
summary(c4g_new_nolag)
new_tddf$predicted_h_c4g_new_nolag <- predict(c4g_new_nolag, type = "risk")
p1 <- ggplot(new_tddf[new_tddf$latency<4,],aes(x = latency, y = predicted_h_c4g_new_nolag, color = Group)) + geom_smooth() 
new_tddf$predicted_logh_c4g_new_nolag <- predict(c4g_new_nolag, type = "lp")
p2 <- ggplot(new_tddf[new_tddf$latency<4,],aes(x = latency, y = predicted_logh_c4g_new_nolag, color = Group)) + geom_smooth() 

pdf(file = 'predicted_logh_c4g_new_nolag.pdf')
ggarrange(p1,p2,nrow = 2)# + facet_wrap(~ID)
dev.off()

# interaction with win -- NS
c4g1 <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag * scale(waitValue) * Group + (1|ID), tddf)
summary(c4g1)
vif.lme(c4g1)
anova(c4g,c4g1)

c4gdemo <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(waitValue) * Group + 
                   scale(waitValue) * scale(age) + scale(waitValue) * scale(education) + (1|ID), tddf)
summary(c4gdemo)

# add immediate quits -- maybe not, coefficients get really crazy...
# ic4g <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(waitValue) * Group + (1|ID), tdf)
# summary(ic4g)
# ic4gdemo <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(waitValue) * Group + 
#                    scale(waitValue) * scale(age) + scale(waitValue) * scale(education) + (1|ID), tdf)
# summary(ic4gdemo)

# forward-looking value
c5 <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(sv) + (1|ID), tddf)
summary(c5)
vif.lme(c5)
c5g <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(sv)*Group + (1|ID), tddf)
summary(c5g)
vif.lme(c5g)

# including both trial-level and forward-looking value improves fits
c6 <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(waitValue) + scale(sv) + (1|ID), tddf)
summary(c6)
vif.lme(c6)

c7 <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(initialTime)*scale(waitValue) + scale(initialTime)*scale(sv) + (1|ID), tddf)
summary(c7)
vif.lme(c7)

######################
# "best-fitting" model #
######################
c7g <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(initialTime)*scale(waitValue)*Group + scale(initialTime)*scale(sv)*Group + (1|ID), tddf)
summary(c7g)
vif.lme(c7g)

anova(c4,c4g,c5,c5g,c6,c7,c7g)
# trial-level value as time-varying covariate

##########
# individual differences on waitValue (trial-level)
##########

## cognition
# no effect of WTAR
c4gdemo_wtar <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(waitValue) * Group + 
                        scale(waitValue) * scale(WTARSS) +
                        scale(waitValue) * scale(age) + scale(waitValue) * scale(education) + (1|ID), tddf)
summary(c4gdemo_wtar)

# exit disrupts, group, education
c4gdemo_exit <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(waitValue) * Group + 
                        scale(waitValue) * scale(EXITtot) +
                        scale(waitValue) * scale(age) + scale(waitValue) * scale(education) + (1|ID), tddf)
summary(c4gdemo_exit)

## impulsivity
# SPSI ICS increases value effect...
c4gdemo_ICS <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(waitValue) * Group + 
                        scale(waitValue) * scale(ICSSUB) +
                        scale(waitValue) * scale(age) + scale(waitValue) * scale(education) + (1|ID), tddf)
summary(c4gdemo_ICS)

# BIS nonplanning increases
c4gdemo_nonplan <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(waitValue) * Group + 
                        scale(waitValue) * scale(NONPLAN) +
                        scale(waitValue) * scale(age) + scale(waitValue) * scale(education) + (1|ID), tddf)
summary(c4gdemo_nonplan)

# UPPS neg. urgency disrupts
c4gdemo_neg_urg <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(waitValue) * Group + 
                        scale(waitValue) * scale(UPPSPNEGURGENCY) +
                        scale(waitValue) * scale(age) +
                          # scale(waitValue) * scale(education) + 
                          (1|ID), tddf)
summary(c4gdemo_neg_urg)

# pos urgency increases
c4gdemo_pos_urg <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(waitValue) * Group + 
                           scale(waitValue) * scale(UPPSPPOSURGENCY) +
                           scale(waitValue) * scale(age) + scale(waitValue) * scale(education) + (1|ID), tddf)
summary(c4gdemo_pos_urg)

# premeditation did not converge correctly with age and education, but with each individually the effect is similar to persev
c4gdemo_premed <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + 
                          # scale(waitValue) * Group + 
                           scale(waitValue) * scale(UPPSPLACKOFPREMED) +
                           scale(waitValue) * scale(age) +
                          scale(waitValue) * scale(education) +
                          (1|ID), tddf)
summary(c4gdemo_premed)
vif.lme(c4gdemo_premed)

# perseveration increases, PARTLY EXPLAINS ATT-DEP DIFFERENCE
c4gdemo_persev <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(waitValue) * Group + 
                           scale(waitValue) * scale(UPPSPLACKOFPERSEV) +
                           scale(waitValue) * scale(age) + scale(waitValue) * scale(education) + (1|ID), tddf)
summary(c4gdemo_persev)



##########
# individual differences on sv (forward-looking): weak, inconsistent effects
##########

## cognition
c5gdemo_wtar <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(sv) * Group + 
                        scale(sv) * scale(WTARSS) +
                        scale(sv) * scale(age) + scale(sv) * scale(education) + (1|ID), tddf)
summary(c5gdemo_wtar)

c5gdemo_exit <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(sv) * Group + 
                        scale(sv) * scale(EXITtot) +
                        scale(sv) * scale(age) + scale(sv) * scale(education) + (1|ID), tddf)
summary(c5gdemo_exit)

## impulsivity
c5gdemo_ICS <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(sv) * Group + 
                       scale(sv) * scale(ICSSUB) +
                       scale(sv) * scale(age) + scale(sv) * scale(education) + (1|ID), tddf)
summary(c5gdemo_ICS)

c5gdemo_nonplan <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(sv) * Group + 
                           scale(sv) * scale(NONPLAN) +
                           scale(sv) * scale(age) + scale(sv) * scale(education) + (1|ID), tddf)
summary(c5gdemo_nonplan)

c5gdemo_neg_urg <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(sv) * Group + 
                           scale(sv) * scale(UPPSPNEGURGENCY) +
                           scale(sv) * scale(age) +
                           # scale(sv) * scale(education) + 
                           (1|ID), tddf)
summary(c5gdemo_neg_urg)

c5gdemo_pos_urg <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(sv) * Group + 
                           scale(sv) * scale(UPPSPPOSURGENCY) +
                           scale(sv) * scale(age) + scale(sv) * scale(education) + (1|ID), tddf)
summary(c5gdemo_pos_urg)

c5gdemo_premed <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + 
                          # scale(sv) * Group + 
                          scale(sv) * scale(UPPSPLACKOFPREMED) +
                          scale(sv) * scale(age) +
                          scale(sv) * scale(education) +
                          (1|ID), tddf)
summary(c5gdemo_premed)
vif.lme(c5gdemo_premed)

c5gdemo_persev <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(sv) * Group + 
                          scale(sv) * scale(UPPSPLACKOFPERSEV) +
                          scale(sv) * scale(age) + scale(sv) * scale(education) + (1|ID), tddf)
summary(c5gdemo_persev)


# interactions with time -- these don't seem plausible, let's not over-complicate things
# 
# mt00 <- coxph(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag +  Group  
#              , tddf)
# summary(mt00)
# # zp <- cox.zph(mt0, transform= function(latency) log(latency +1))
# zp <- cox.zph(mt00)
# 
# plot(zp[1]) 
# abline(0,0, col=2)
# 
# 
# mt0 <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(waitValue) + scale(waitValue)*Group + scale(waitValue):scale(tPeriod) + scale(waitValue):scale(tPeriod):Group
#             + (1|ID) , tddf)
# summary(mt0)
# vif.lme(mt0)
# zp <- cox.zph(mt0, transform= function(latency) log(latency +1))
# zp <- cox.zph(mt0)
# 
# plot(zp[4]) 
# abline(0,0, col=2)
# abline(h= mt0$coef[4], col=3, lwd=2, lty=2)
save(list = ls(all.names = TRUE), file = "wtw_results.RData")
 load(file = "wtw_results.RData")

