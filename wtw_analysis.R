
setwd("~/Box Sync/Project_wtw/R")
library(readr)
# library(lmerTest)
library(ggplot2)
library(tidyr)
library(tibble)
# library(xtable)
library(Hmisc)
# library(nnet)
library(reshape2)
# library(ggbiplot)
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


load(file = "wtw.RData")

#  for collinearity diagnostics
vif.lme <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)] }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v }



# what predicts immediate quits?

# im00trial <- glmer(
#   immQuit ~  scale(trial) + immQuit_lag  +  win_lag  +
#     (1 | ID),
#   family = binomial(),
#   data = lagddf,
#   glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
# # vif.lme(im00)

## these analyses are now valid, no multicollinearity once we censored post-immediate-quit trials

im00 <- glmer(
  immQuit ~  scale(initialTime) + win_lag  +
    (1 | ID),
  family = binomial(),
  data = lagddf,
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

vif.lme(im00)
anova(im00,im00trial)
# time predicts better than trial

im0 <- glmer(
  immQuit ~  scale(initialTime) * Group + win_lag * Group +
    (1 | ID),
  family = binomial(),
  data = lagddf,
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(im0)
car::Anova(im0, type = 'III')
# vif.lme(im0)

im1 <- glmer(
  immQuit ~  scale(initialTime)* Group +   win_lag * Group +
    scale(initialTime)* scale(education) +  win_lag * scale(education) +
    scale(initialTime)* scale(age) +   win_lag * scale(age) +
    (1 | ID),
  family = binomial(),
  data = lagddf,
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
# lagged immediate quit/lagged win collinearity problematic -- will need to only look at post-non-immediate-quit trials
vif.lme(im1)
summary(im1)
car::Anova(im1, '3')
stargazer(im00, im0, im1, type="html", out="immQuits.htm", digits = 3,single.row=TRUE,omit.stat = "bic",
          column.labels = c("Trial-level variables", "+ group", "+ demographics"),
          star.char = c("+", "*", "**", "***"),
          star.cutoffs = c(0.1, 0.05, 0.01, 0.001),
          notes = c("+ p<0.1; * p<0.05; ** p<0.01; *** p<0.001"), 
          notes.append = F)


plot_models(im00,im0,im1)
im0good <- glmer(
  immQuit ~  scale(initialTime) * Group + win_lag * Group +
    (1 | ID),
  family = binomial(),
  data = lagddf[!lagddf$bad,],
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(im0good)
car::Anova(im0good, type = 'III')


pdf("immediate_quits_by_group.pdf", height = 4, width = 6 )
plot(emmeans(im0, ~ win_lag | Group , p.kr = TRUE), horiz = F)
dev.off()

# sensitivity analysis: exclude bad people

# sensitivity analyses: exclude trials post-immediate quit to get rid of collinearity with win/quit, findings hold
im0_sens <- glmer(
  immQuit ~  scale(initialTime) * Group + win_lag * Group +
    (1 | ID),
  family = binomial(),
  data = lagddf,
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(im0_sens)
car::Anova(im0_sens, type = 'III')
vif.lme(im0_sens)
plot(emmeans(im0_sens, ~ win_lag | Group , p.kr = FALSE), horiz = F)

# control for previous scheduled wait time

im00w <- glmer(
  immQuit ~  scale(initialTime) + win_lag * scale(designatedWait_lag) +
    (1 | ID),
  family = binomial(),
  data = lagddf,
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(im00w)
vif.lme(im00w)
# time predicts better than trial

# no 2-way (wait * quit) or 3-way (...* group) interaction
im0w <- glmer(
  immQuit ~  scale(initialTime)  + win_lag * scale(designatedWait_lag)  * Group + 
    (1 | ID),
  family = binomial(),
  data = lagddf,
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(im0w)
car::Anova(im0w, type = 'III')
vif.lme(im0w)
plot(emmeans(im0w, ~ win_lag | Group , p.kr = TRUE), horiz = F)
# figure for the grant
s <- summary(im0w)
coef <- s$coefficients
terms1 <- labels(coef)[[1]]
terms1 <- terms1[c(3,9:11)]
labels1 <- rev(c("Previous failure vs. successful wait", "Fail.*Controls (n=68)", "Fail.*Depressed (n=63)", "Fail.*Ideators(n=59)"));
p <- plot_model(im0w, terms = terms1,show.values = T, show.p = T,axis.labels = labels1,  
           axis.title = "odds ratio  Wait <-> Quit",
           value.offset = 0.35,vline.color = "slategray3", axis.lim = c(.01,10),
           title = "Quitting after a failure of persistence.  Reference group: attempters (n=87)")
p <- p + font_size(labels.y = 12) + theme_minimal()
ggsave("Immediate_quits_attempters.pdf", width = 5,height = 3)
im1w <- glmer(
  immQuit ~  scale(initialTime)  + win_lag * scale(designatedWait_lag) + win_lag * Group + scale(designatedWait_lag) * Group +
    win_lag * scale(EDUCATION) + scale(designatedWait_lag) * scale(EDUCATION) +
    win_lag * scale(BASELINEAGE) + scale(designatedWait_lag) * scale(BASELINEAGE) +
    (1 | ID),
  family = binomial(),
  data = lagddf,
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
# lagged immediate quit/lagged win collinearity problematic -- will need to only look at post-non-immediate-quit trials
vif.lme(im1w)

pdf("immediate_quits_by_group.pdf", height = 4, width = 6 )
plot(emmeans(im0, ~ win_lag | Group , p.kr = TRUE), horiz = F)
dev.off()


stargazer(im00w, im0w, im1w, type="html", out="design_immQuits.htm", digits = 3,single.row=TRUE,omit.stat = "bic",
          column.labels = c("Trial-level variables", "+ group", "+ demographics"),
          star.char = c("+", "*", "**", "***"),
          star.cutoffs = c(0.1, 0.05, 0.01, 0.001),
          notes = c("+ p<0.1; * p<0.05; ** p<0.01; *** p<0.001"), 
          notes.append = F)



# individual differences in immediate quits: puzzling/intriguing that all psychopathology measures go in the opposite direction to suicide attempt

im2 <- glmer(
  immQuit ~  scale(initialTime)* Group +   win_lag * Group +
    scale(initialTime)* scale(education) +  win_lag * scale(education) +
    scale(initialTime)* scale(age) +   win_lag * scale(age) +
    scale(initialTime)* scale(UPPSPNEGURGENCY) +   win_lag * scale(UPPSPNEGURGENCY) +
    (1 | ID),
  family = binomial(),
  data = lagddf,
  # data = lagddf[!lagddf$group1245==1,],
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
# lagged immediate quit/lagged win collinearity problematic -- will need to only look at post-non-immediate-quit trials
vif.lme(im2)

stargazer(im0, im1, im2, type="html", out="immQuits_impuls.htm", digits = 2,single.row=TRUE,omit.stat = "bic",
          column.labels = c("Trial-level variables + group", "+ demographics", "+ impulsivity"),
          star.char = c("+", "*", "**", "***"),
          star.cutoffs = c(0.1, 0.05, 0.01, 0.001),
          notes = c("+ p<0.1; * p<0.05; ** p<0.01; *** p<0.001"), 
          notes.append = F)

em2 <- emmeans(im2, ~ win_lag |  Group * UPPSPNEGURGENCY, at = list(UPPSPNEGURGENCY = c(14,25,37)),  p.kr = TRUE)
pdf("immediate_quits_grp_neg_urg.pdf", height = 6, width = 18 )
plot(em2, horiz = F)
dev.off()

# what predicts latency
# hist(qdf$latency[qdf$GROUP1245==1],100)
# hist(qdf$latency[qdf$GROUP1245==2],100)
# hist(qdf$latency[qdf$GROUP1245==4],100)
# hist(qdf$latency[qdf$GROUP1245==5],100)

# quits only
ggplot(qdf, aes(x = outcomeTime, y = log(latency), color = Group)) + stat_smooth(method = "gam")

# delayed quits and 20s wins
ggplot(wdf, aes(x = outcomeTime, y = log(latency), color = Group)) + stat_smooth(method = "gam")

#  all latencies
ggplot(df, aes(x = outcomeTime, y = log(latency), color = Group)) + stat_smooth(method = "gam")

ggplot(qdf, aes(x = latency, color = Group, lty = period)) + geom_density(adjust = 2) #+ facet_wrap( ~ period)
ggplot(qdf, aes(x = latency, color = Group, lty = period)) + geom_density(adjust = 1)

# ggplot(qdf, aes(x = latency, color = GroupLeth)) + geom_density(adjust = 1.5) + facet_wrap( ~ period)
ggplot(qdf, aes(x = latency, color = GroupLeth, lty = period)) + geom_density(adjust = 1.5)


ggplot(df, aes(x = outcomeTime, y = as.numeric(immQuit), color = Group)) + stat_smooth()


# reward rate
# ggplot(ddf[ddf$trial>1,], aes(x = initialTime, y = rewardRate, color = Group)) + facet_wrap(~win_lag) + stat_smooth(method = "gam")
# ggplot(df[df$trial>1,], aes(x = initialTime, y = rewardRate, color = GroupLeth)) + stat_smooth(method = "gam")
# 

pdf("latencies_over_time_by_group.pdf", width = 6, height = 4)
ggplot(qdf, aes(x = initialTime, y = latency, color = Group)) + stat_smooth(method = "loess")
dev.off()

pdf("wtw_over_time_by_group.pdf", width = 6, height = 4)
ggplot(wdf, aes(x = initialTime, y = latency, color = Group)) + stat_smooth(method = "loess")
dev.off()


# # reward rate by group on WTW trials
# rrm0 <- lme4::lmer(rewardRate ~ scale(outcomeTime) + win_lag + (1 | ID), data = ddf)
# rrm1 <- lme4::lmer(rewardRate ~ scale(outcomeTime) * Group +  win_lag * Group + (1 | ID), data = ddf)
# rrm2 <- lme4::lmer(rewardRate ~ scale(outcomeTime)  * Group + win_lag * Group +
#               scale(outcomeTime)*scale(EDUCATION) + scale(outcomeTime)*scale(BASELINEAGE) + 
#               scale(outcomeTime)* GENDERTEXT +
#                 win_lag*scale(EDUCATION) + win_lag*scale(BASELINEAGE) + 
#                 win_lag* GENDERTEXT +
#                 (1 | ID), data = ddf)
# stargazer(rrm0, rrm1, rrm2, type="html", out="reward_rate_by_group.htm", digits = 3,single.row=TRUE,omit.stat = "bic",
#           column.labels = c("Trial-level variables", "+ group", "+ demographics"),
#           star.char = c("+", "*", "**", "***"),
#           star.cutoffs = c(0.1, 0.05, 0.01, 0.001),
#           notes = c("+ p<0.1; * p<0.05; ** p<0.01; *** p<0.001"), 
#           notes.append = F)
# 
# # by lethality: LL have a lower RR
# rrm0 <- lme4::lmer(rewardRate ~ scale(outcomeTime) +  (1 | ID), data = wdf)
# rrm1 <- lme4::lmer(rewardRate ~ scale(outcomeTime) * GroupLeth +  (1 | ID), data = wdf)
# rrm2 <- lme4::lmer(rewardRate ~ scale(outcomeTime)  * GroupLeth + 
#                      scale(outcomeTime)*scale(EDUCATION) + scale(outcomeTime)*scale(BASELINEAGE) + 
#                      scale(outcomeTime)* GENDERTEXT +
#                      (1 | ID), data = wdf)
# stargazer(rrm0, rrm1, rrm2, type="html", out="reward_rate_by_group_leth.htm", digits = 2,single.row=TRUE,omit.stat = "bic",
#           column.labels = c("Trial-level variables", "+ group", "+ demographics"),
#           star.char = c("+", "*", "**", "***"),
#           star.cutoffs = c(0.1, 0.05, 0.01, 0.001),
#           notes = c("+ p<0.1; * p<0.05; ** p<0.01; *** p<0.001"), 
#           notes.append = F)
# 

# no group differences in final earnings
fm1 <- lm(finalEarnings ~ Group, data = fdf)
summary(fm1)

# build a better model of reward rate predictors
summary(rr0 <- lme4::lmer(rewardRate ~ scale(outcomeTime) + win_lag + lat_lag + rewardRate_lag + (1 | ID), data = ddf))


################################
df$`Previous win` <- df$win_lag

km0 <- survfit(Surv(latency, !cens) ~ Group, data=ddf)
p0 <- ggsurvplot(km0, data = df, risk.table = FALSE,  conf.int = TRUE, xlim = c(0,19), break.time.by = 2, ggtheme = theme_minimal())
pdf("km_survival_by_group_simple.pdf", width = 6, height = 4)
print(p0)
dev.off()

# KM survival: all data
km1 <- survfit(Surv(latency, !cens) ~ Group + df$win_lag, data=df)
p1 <- ggsurvplot(km1, data = df, risk.table = FALSE,  facet.by = c("win_lag"), conf.int = TRUE, xlim = c(0,19), break.time.by = 2, ggtheme = theme_minimal())

# KM survival excluding preceding immediate quits
km2 <- survfit(Surv(latency, !cens) ~ win_lag + Group, data=lagddf)
p2 <- ggsurvplot(km2, data = lagddf, risk.table = FALSE, 
                 facet.by = c("win_lag"), #group.by = "Group",
                 conf.int = TRUE, xlim = c(0,19), break.time.by = 2, ggtheme = theme_minimal())
# p2 <- p2 + theme(legend.position="none")

p2 <- p2 + scale_color_manual(values=c( "#C77CFF","#F8766D", "#7CAE00", "#00BFC4")) + 
  scale_fill_manual(values=c("#C77CFF","#F8766D", "#7CAE00", "#00BFC4"))

pdf("km_survival_by_group_and_previous_win.pdf", width = 6, height = 4)
p2
dev.off()


km3 <- survfit(Surv(latency, !cens) ~ Group + win_lag, data=ddf)
# ggsurvplot(km1a, data = ddf, risk.table = FALSE, conf.int = FALSE, xlim = c(0,20))
p3 <- ggsurvplot(km3, data = ddf, risk.table = FALSE, facet.by = c("win_lag"), conf.int = TRUE, xlim = c(0,19), break.time.by = 2, ggtheme = theme_minimal())


pdf("km_survival_by_group.pdf", width = 12, height = 8)
# ggarrange(p1,p2,p3,  labels = c("All trials","Excluding trials after immediate quits", "Excluding immediate quits"), hjust = c(-.2),legend = "bottom")
ggarrange(p2,p3,  labels = c("Excluding trials after immediate quits", "Excluding immediate quits"), 
          hjust = c(-.2),legend = "bottom")
dev.off()

pdf("km_survival_by_group_main_graph.pdf", width = 12, height = 8)
ggarrange(p2,  labels = c("Excluding trials after immediate quits"), hjust = c(-1.2),legend = "bottom")
dev.off()

# without bad subjects
# KM survival: all data
km1 <- survfit(Surv(latency, !cens) ~ Group + win_lag, data=df[!df$bad,])
p1 <- ggsurvplot(km1, data = df, risk.table = FALSE, facet.by = c("win_lag"), conf.int = TRUE, xlim = c(0,19), break.time.by = 2, ggtheme = theme_minimal())

# KM survival excluding immediate quits
km2 <- survfit(Surv(latency, !cens) ~ Group + win_lag, data=lagddf[!lagddf$bad,])
p2 <- ggsurvplot(km2, data = lagddf, risk.table = FALSE, facet.by = c("win_lag"), conf.int = TRUE, xlim = c(0,19), break.time.by = 2, ggtheme = theme_minimal())
# p2 <- p2 + theme(legend.position="none")

km3 <- survfit(Surv(latency, !cens) ~ Group + win_lag, data=ddf[!ddf$bad,])
# ggsurvplot(km1a, data = ddf, risk.table = FALSE, conf.int = FALSE, xlim = c(0,20))
p3 <- ggsurvplot(km3, data = ddf, risk.table = FALSE, facet.by = c("win_lag"), conf.int = TRUE, xlim = c(0,19), break.time.by = 2, ggtheme = theme_minimal())


pdf("km_survival_by_group_good.pdf", width = 12, height = 8)
ggarrange(p1,p2,p3,  labels = c("All trials","Excluding trials after immediate quits", "Excluding immediate quits"), hjust = c(-.2),legend = "bottom")
dev.off()

pdf("km_survival_by_group_main_graph_good.pdf", width = 12, height = 8)
ggarrange(p2,  labels = c("Excluding trials after immediate quits"), hjust = c(-1.2),legend = "bottom")
dev.off()

# plot hazard functions
h1 <- ggsurvplot(km1, data = df, risk.table = FALSE, fun = "cumhaz", facet.by = c("win_lag"), conf.int = TRUE, xlim = c(0,19), break.time.by = 2, ggtheme = theme_minimal())
h2 <- ggsurvplot(km2, data = lagddf, risk.table = FALSE, fun = "cumhaz", facet.by = c("win_lag"), conf.int = TRUE, xlim = c(0,19), ylim = c(0,3),break.time.by = 2, ggtheme = theme_minimal())
h3 <- ggsurvplot(km3, data = ddf, risk.table = FALSE, fun = "cumhaz", facet.by = c("win_lag"), conf.int = TRUE, xlim = c(0,19), ylim = c(0,3), break.time.by = 2, ggtheme = theme_minimal())


# difference in survival
survdiff(Surv(latency, !cens) ~ Group, data=lagddf)

# survival regression with Weibull dist.
w00 <- survreg(Surv(latency, !cens) ~ lat_lag + trial + win_lag * designatedWait_lag,dist = "weibull", data=lagddf)
w0 <- survreg(Surv(latency, !cens) ~ lat_lag + trial + win_lag * designatedWait_lag,dist = "weibull", data=lagddf)
# best model:
w1 <- survreg(Surv(latency, !cens) ~ lat_lag +  trial + win_lag * designatedWait_lag + Group * win_lag  + Group * designatedWait_lag,dist = "weibull", data=lagddf)
w1.5 <- survreg(Surv(latency, !cens) ~ lat_lag + Group * trial + win_lag * designatedWait_lag + Group * win_lag  + Group * designatedWait_lag,dist = "weibull", data=lagddf)
w2 <- survreg(Surv(latency, !cens) ~ lat_lag + Group * trial + Group * win_lag * designatedWait_lag,dist = "weibull", data=lagddf)
anova(w00,w0,w1,w1.5,w2)
summary(w1)

# without immediate quits
w1d <- survreg(Surv(latency, !cens) ~ lat_lag +  trial * Group + win_lag * designatedWait_lag + Group * win_lag  + Group * designatedWait_lag,dist = "weibull", data=ddf)
summary(w1d)
# pct <- 1:98/100 
pw1 <- predict(w1, newdata=list(win_lag = TRUE),type="quantile",p=seq(.01,.99,by=.01))
lines(pw1,seq(.99,.01,by=-.01),col="red")

library(flexsurv)
sWei  <- flexsurvreg(Surv(latency, !cens) ~ lat_lag + trial + win_lag * designatedWait_lag, dist = "gengamma", data=lagddf[,])
plot(sWei, type = 'hazard', ci = TRUE)

sWeiGroup  <- flexsurvreg(Surv(latency, !cens) ~ lat_lag + trial + win_lag * designatedWait_lag + win_lag * Group ,dist = "weibull", data=lagddf)
plot(sWeiGroup, newdata = lagddf, type = 'hazard', ci = TRUE)


## proportional hazards models

## Cox regression model censoring immediate quits and post-immediate quit trials
c1 <- coxph(Surv(latency, !cens) ~ Group, ddf1)

## mixed-effects Cox regression models on the same
c2 <- coxme(Surv(latency, !cens) ~ scale(lat_lag) + scale(initialTime) + win_lag + (1|ID), ddf1)
c3 <- coxme(Surv(latency, !cens) ~ scale(lat_lag) + scale(initialTime) * Group + win_lag * Group + (1|ID), ddf1)
print(c3)
car::Anova(c3)
anova(c2,c3)

# same, adding current immediate quits: the effect seems to be driven by post-loss immediate quits in attempters
c4 <- coxme(Surv(latency, !cens) ~ scale(lat_lag) + scale(initialTime) * Group + win_lag * Group + (1|ID), lagddf)
print(c4)
car::Anova(c4)
vif.lme(c4)
anova(c2,c3)

# add individual differences; nothing for EXIT, education, BIS non-planning
# negative effect for ICS
c5 <- coxme(Surv(latency, !cens) ~ scale(lat_lag) + scale(initialTime) * scale(ICSSUB) + (1|ID), ddf)
print(c5)
car::Anova(c5, type = '3')
c5g <- coxme(Surv(latency, !cens) ~ scale(lat_lag) + scale(initialTime) * Group + scale(initialTime) * scale(ICSSUB) + (1|ID), ddf)
print(c5g)
car::Anova(c5g, type = '3')


# NEG, POS urgency
c6 <- coxme(Surv(latency, !cens) ~ scale(lat_lag) + scale(initialTime) * scale(UPPSPNEGURGENCY) + (1|ID), ddf)
print(c6)
car::Anova(c6, type = '3')
# positive interaction (acceleration) with Pos urgency
c7 <- coxme(Surv(latency, !cens) ~ scale(lat_lag) + scale(initialTime) * scale(UPPSPPOSURGENCY) + (1|ID), ddf)
print(c7)
car::Anova(c7, type = '3')
# stands controlling for group
c7g <- coxme(Surv(latency, !cens) ~ scale(lat_lag) + scale(initialTime) * Group + scale(initialTime) * scale(UPPSPPOSURGENCY) + (1|ID), ddf)
print(c7g)
car::Anova(c7g, type = '3')


# weak for PERSEV
c8 <- coxme(Surv(latency, !cens) ~ scale(lat_lag) + scale(initialTime) * scale(UPPSPLACKOFPERSEV) + (1|ID), ddf)
print(c8)
car::Anova(c8, type = '3')
# negative for premed.
c9 <- coxme(Surv(latency, !cens) ~ scale(lat_lag) + scale(initialTime) * scale(UPPSPLACKOFPREMED) + (1|ID), ddf)
print(c9)
car::Anova(c9, type = '3')

# substance: nothing
c10 <- coxme(Surv(latency, !cens) ~ scale(lat_lag) + scale(initialTime) * SubstanceLifetime + scale(initialTime) * Group +  (1|ID), ddf)
print(c10)
car::Anova(c10, type = '3')

# anxiety: nothing
c11 <- coxme(Surv(latency, !cens) ~ scale(lat_lag) + scale(initialTime) * AnxietyLifetime + scale(initialTime) * Group +  (1|ID), ddf)
print(c11)
car::Anova(c11, type = '3')


# remove bad subjects
c4good <- coxme(Surv(latency, !cens) ~ scale(lat_lag) + scale(initialTime) * Group + win_lag * Group + (1|ID), lagddf[!lagddf$bad,])
print(c4good)


# mixPHM package seems useless because it does not deal with censoring

library(muhaz)
# piecewise <- pehaz(ddf$latency, delta=ddf$quit, width=NA, min.time=0, max.time=20.1)
# plot(piecewise, xlab="Time, seconds", ylab="Quit Rate")

pwC <- pehaz(na.omit(ddf$latency[ddf$group1245==1]), delta=ddf$quit[ddf$group1245==1 & !is.na(ddf$latency)], width=.5, min.time=0, max.time=20)
pwD <- pehaz(na.omit(ddf$latency[ddf$group1245==2]), delta=na.omit(ddf$quit[ddf$group1245==2]), width=.5, min.time=0, max.time=20)
pwI <- pehaz(na.omit(ddf$latency[ddf$group1245==4]), delta=na.omit(ddf$quit[ddf$group1245==4]), width=.5, min.time=0, max.time=20)
pwA <- pehaz(na.omit(ddf$latency[ddf$group1245==5]), delta=na.omit(ddf$quit[ddf$group1245==5]), width=.5, min.time=0, max.time=20)

h <- c(pwC$Hazard,pwD$Hazard,pwI$Hazard,pwA$Hazard)
times <- c(pwC$Cuts[2:length(pwC$Cuts)],pwD$Cuts[2:length(pwD$Cuts)],pwI$Cuts[2:length(pwI$Cuts)],pwA$Cuts[2:length(pwA$Cuts)])
g <- c(rep('Controls',length(pwC$Hazard)), rep('Depressed',length(pwD$Hazard)),rep('Ideators',length(pwI$Hazard)),rep('Attempters',length(pwA$Hazard)))
H <- as.tibble(cbind(h,times,g))
H$h <- as.numeric(H$h)
H$times <- as.numeric(H$times)
pdf("piecewise_hazard_fx_all_groups.pdf", width = 12, height = 4)
ggplot(H,aes(times, h, color = g)) + geom_line()
dev.off()
# 
# pdf("piecewise_hazard_fx_controls.pdf", width = 6, height = 4)
# plot(pwC,ylim = c(0,0.15), xlab="Time, seconds", ylab="Quit Rate", col = "dark green")
# dev.off()
# 
# pdf("piecewise_hazard_fx_dep.pdf", width = 6, height = 4)
# plot(pwD,ylim = c(0,0.15), xlab="Time, seconds", ylab="Quit Rate", col = "blue")
# dev.off()
# 
# pdf("piecewise_hazard_fx_ideators.pdf", width = 6, height = 4)
# plot(pwI, ylim = c(0,0.15),xlab="Time, seconds", ylab="Quit Rate", col = "purple")
# dev.off()
# pdf("piecewise_hazard_fx_att.pdf", width = 6, height = 4)
# plot(pwA, ylim = c(0,0.15),xlab="Time, seconds", ylab="Quit Rate", col = "red")
# dev.off()
# 
# ggarrange(p1,p2,p3, p4,  labels = c("Controls","Depressed", "Ideators", "Attempters"), hjust = c(-.2))
# dev.off()


# trial-level value as time-varying covariate
cV1 <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag  + scale(waitValue)*Group + (1|ID), tddf)
summary(cV1)

# # time * value interactions do not improve fits
# cV1Tasktime <- coxme(Surv(latency, quit) ~ scale(lat_lag) + win_lag  + scale(initialTime) * scale(waitValue) * Group + (1|ID), tddf)
# summary(cV1Tasktime)

# forward-looking value estimate: better on AIC
cVf <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag  + scale(sv)*Group + (1|ID), tddf)
summary(cVf)

# cVfTasktime <- coxme(Surv(latency, quit) ~ scale(lat_lag) + win_lag  + scale(initialTime) * scale(sv) + 
#                        scale(sv) * Group + scale(initialTime) * Group + (1|ID), tddf)
# summary(cVfTasktime)



anova(cV1,cVf)

# same with flexsurvreg
cVwei  <- flexsurvreg(Surv(latency, quit) ~ scale(lat_lag) + scale(trial) + win_lag * designatedWait_lag + scale(waitValue) * Group ,dist = "weibull", data=tddf)
plot(cVwei, type = 'hazard', ci = TRUE)


save(list = ls(all.names = TRUE), file = "wtw_results.RData")
load(file = "wtw_results.RData")

