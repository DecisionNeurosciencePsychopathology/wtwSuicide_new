
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

# plot contingency
ggplot(na.omit(df),aes(x=designatedWait)) +  geom_histogram(aes(y = (..count..)/sum(..count..)),bins = 40) + 
  labs(x = "Designated wait, seconds", y = "Probability of 20c reward") + theme_minimal()
ggsave("wtw_discrete_contingency.pdf", width = 4, height = 2.2)
ggplot(v[v$t2>.4,],aes(t2,sv)) + geom_line() + theme_minimal() + labs(x = "Designated wait, seconds", y = "Value of waiting")
ggsave("wtw_sv_discrete_contingency.pdf", width = 4, height = 2.2)

# inspect hazard functions using raw data
library(muhaz)
# piecewise <- pehaz(ddf$latency, delta=ddf$quit, width=NA, min.time=0, max.time=20.1)
# plot(piecewise, xlab="Time, seconds", ylab="Quit Rate")

pwC <- pehaz(na.omit(df$latency), delta=df$quit[!is.na(df$latency)], width=.5, min.time=0, max.time=40)
h <- pwC$Hazard
times <- pwC$Cuts[2:length(pwC$Cuts)]
H <- as.tibble(cbind(h,times))
H$h <- as.numeric(H$h)
H$logh <- log(H$h)
H$times <- as.numeric(H$times)
H$n <- pwC$At.Risk
# pdf("piecewise_hazard_fx_all_groups.pdf", width = 12, height = 4)
p1 <- ggplot(H,aes(times, h, alpha = n)) + geom_line(size = 1.5)+ 
  theme_minimal() + labs(x = "Time, seconds", y = "Hazard of quitting")
p1s<- ggplot(H,aes(times, h, color = group)) + geom_smooth(method = "loess", span = .33, se = F)
# dev.off()

# pdf("piecewise_log_hazard_all_groups.pdf", width = 12, height = 4)
p2 <- ggplot(H,aes(times, logh, color = g)) + geom_line()+ geom_smooth(method = "loess", span = .33, se = F, size = .5, linetype = 3)
p2s<- ggplot(H,aes(times, logh, color = group)) + geom_smooth(method = "loess", span = .33, se = F)

# dev.off()
p3 <- ggplot(v[v$t2>.4,],aes(t2,waitValue)) + geom_line()
p4 <- ggplot(v[v$t2>.4,],aes(t2,-sv)) + geom_line() + theme_minimal() + labs(x = "Time, seconds", y = "Value, reversed")

ggarrange(p1,p2,p3,p4, common.legend = T,nrow = 4)
ggsave("piecewise_hazard_all_groups_value.pdf", width = 10, height = 10)

ggarrange(p1,p4, common.legend = T,nrow = 2)
ggsave("piecewise_hazard_all_groups_negSV.pdf", width = 8, height = 5)


ggarrange(p1s,p2s,p3,p4, common.legend = T,nrow = 4)
ggsave("smoothed_piecewise_hazard_all_groups_value.pdf", width = 10, height = 10)

# alternate visualization:

p5 <- ggplot(tddf2,aes(x = latency,y = as.numeric(quit), color = Group)) + geom_smooth(method = "loess")

# to understand the difference in the first second, add the immediate quits:
p6 <- ggplot(tdf2,aes(x = latency,y = as.numeric(quit), color = Group)) + geom_smooth(method = "loess")


#  add value
vCoarse <- subset()
vCoarse <- v[is.element(v$t2,unique(H$times)),]


summary(m1 <- lme4::glmer(quit ~ scale(sv) + latency + (1|ID), tddf, family = binomial,
                  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))

# one reason for the changing sign for the waitValue effect is its confounding with the 'baseline' hazard function
# the model below estimates a completely general piecewise hazard function in 2hz time bins, and the sign for the 
# waitValue coefficient does indeed change (still misses the tolerance by a small margin even with these
# optimizer settings)
# m1t_factor <- lme4::glmer(quit ~ as.factor(t2hz) + scale(waitValue) + (1|ID), tddf, family = binomial,
#                   glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
# summary(m1t_factor)

# summary(m2 <- lme4::glmer(quit ~  scale(sv) + (as.factor(t2hz)|ID), tddf, family = binomial,
#                   glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))

# summary(m3 <- lme4::glmer(quit ~ scale(latency)*Group + scale(waitValue)*Group + 
#                             (1|ID), tddf, family = binomial,
#                           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))
# summary(m4 <- lme4::glmer(quit ~ scale(latency)*Group + scale(sv)*Group + (1|ID), tddf, family = binomial,
#                           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))

# without latency
summary(m1_2 <- lme4::glmer(quit ~ scale(sv)*scale(initialTime)*Group + (1|ID), tddf2, family = binomial,
                            glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))
plot_model(m1_2,show.p = TRUE, show.values = T,  vline.color = "slategray3", value.offset = 0.4,axis.lim = c(0.5,5),
           axis.title = "Wait  < - >  Quit")


summary(m4_2 <- lme4::glmer(quit ~ scale(latency) + scale(sv)*Group + (1|ID), tddf2, family = binomial,
                          glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))

summary(m5_2 <- lme4::glmer(quit ~ scale(latency) + scale(sv)*scale(initialTime)*Group + (1|ID), tddf2, family = binomial,
                            glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))

summary(m5_2l <- lme4::glmer(quit ~ scale(latency) + scale(sv)*scale(initialTime)*GroupLeth + (1|ID), tddf2, family = binomial,
                            glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))
plot_model(m5_2l,show.p = TRUE, show.values = T,  vline.color = "slategray3", value.offset = 0.4,axis.lim = c(0.5,5),
           axis.title = "Wait  < - >  Quit")

# education: sanity check
summary(m6_2 <- lme4::glmer(quit ~ scale(latency)*scale(education) + scale(sv)*scale(initialTime)*scale(education) +
                              scale(latency)*scale(age) + scale(sv)*scale(initialTime)*scale(age) + 
                              (1|ID), tddf2, family = binomial,
                            glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))
plot_model(m6_2,show.p = TRUE, show.values = T,  vline.color = "slategray3",
           axis.title = "Wait  < - >  Quit")

# add age
summary(m7_2 <- lme4::glmer(quit ~ scale(latency)*Group + scale(sv)*scale(initialTime)*Group + 
                              scale(latency)*scale(age) + scale(sv)*scale(initialTime)*scale(age) + 
                              (1|ID), tddf2, family = binomial,
                            glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))
plot_model(m7_2,show.p = TRUE, show.values = T,  vline.color = "slategray3", value.offset = 0.4,axis.lim = c(0.5,5),
           axis.title = "Wait  < - >  Quit")


# impulsivity: SPSI IC
summary(m8_2 <- lme4::glmer(quit ~ scale(latency) + scale(sv)*scale(initialTime)*Group + 
                              scale(sv)*scale(initialTime)*scale(NONPLAN) + 
                              scale(sv)*scale(initialTime)*scale(EXITtot) + 
                              scale(sv)*scale(initialTime)*scale(age) + 
                              (1|ID), tddf2, family = binomial,
                            glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))
plot_model(m8_2,show.p = TRUE, show.values = T,  vline.color = "slategray3", value.offset = 0.4,axis.lim = c(0.5,5),
           axis.title = "Wait  < - >  Quit")

# stronger effect of SV in high-EXIT subjects is very counter-intuitive and to me suggests that sv is not the right variable

# add 2s lagged value to account for exploration
summary(m20 <- lme4::glmer(quit ~ scale(latency) + scale(sv)*scale(initialTime) +
                              (1|ID), tddf2[!is.na(tddf2$sv_lag3s),], family = binomial,
                            glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))

summary(m21 <- lme4::glmer(quit ~ scale(latency) + scale(sv)*scale(initialTime) +
                              scale(sv_lag2s)*scale(initialTime) +
                              (1|ID), tddf2[!is.na(tddf2$sv_lag3s),], family = binomial,
                            glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))
summary(m22 <- lme4::glmer(quit ~ scale(latency) + scale(sv)*scale(initialTime) +
                             scale(sv_lag2s)*scale(initialTime) +
                             scale(sv_lag3s)*scale(initialTime) +
                             (1|ID), tddf2[!is.na(tddf2$sv_lag3s),], family = binomial,
                           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))
summary(m23 <- lme4::glmer(quit ~ scale(latency) + scale(sv_lag2s)*scale(initialTime) +
                             (1|ID), tddf2[!is.na(tddf2$sv_lag3s),], family = binomial,
                           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))

# the addition of a 3s lag does not account for more variance
# 2s shifted SV alone explains behavior best
anova(m20,m21,m22,m23)

# simple analysis of waiting over time
# remove bad subjects to make sure they were not driving the effect
summary(m24 <- lme4::glmer(quit ~ scale(latency)*Group*scale(initialTime) +
                             (1|ID), tddf2[!tddf2$bad,], family = binomial,
                           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))
s <- summary(m24)
coef <- s$coefficients
terms1 <- labels(coef)[[1]]
terms1 <- terms1[c(14:16)]
labels1 <- rev(c("Time * Learning * Controls vs. attempters","Time * Learning * Depressed vs. attempters","Time * Learning * Ideators vs. attempters" ));

p <- plot_model(m24,show.p = TRUE, terms = terms1,axis.labels = labels1, show.values = T,  vline.color = "slategray3", value.offset = 0.3,axis.lim = c(0.75,1.5),
           axis.title = "Wait  < - >  Quit")
p <- p + theme_minimal()
ggsave("time_learning_attempters.pdf", width = 3,height = 2)

# post-hoc -- only first 1s
summary(m24first <- lme4::glmer(quit ~ Group*scale(initialTime) +
                             (1|ID), tddf2[!tddf2$bad & tddf2$latency<3,], family = binomial,
                           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))

# change 
summary(m24firstBIS <- lme4::glmer(quit ~ NONPLAN*scale(initialTime) +
                                  (1|ID), tddf2[!tddf2$bad & tddf2$latency<2,], family = binomial,
                                glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))

summary(m24firstNEG <- lme4::glmer(quit ~ UPPSPNEGURGENCY*scale(initialTime) +
                                     (1|ID), tddf2[!tddf2$bad & tddf2$latency<2,], family = binomial,
                                   glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))


summary(m25 <- lme4::glmer(quit ~ scale(latency) +
                             scale(sv_lag2s)*Group*scale(initialTime) +
                             (1|ID), tddf2[!is.na(tddf2$sv_lag2s),], family = binomial,
                           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))

# plot hazard functions by group over time
p <- ggplot(tddf2[!tddf2$bad,],aes(latency,quit,color = Group)) + geom_smooth(method = "loess") + facet_wrap(Group~initialTime>60,ncol = 2) 
ggsave("smoothed_piecewise_hazard_over_time_group.pdf", width = 6, height = 10)
library(splines)
p <- ggplot(tddf[!tddf$bad,],aes(latency,quit,color = Group)) + 
  geom_smooth(   method = "glm",
                 method.args = list(family = binomial),
                 formula = y ~ splines::ns(x, df = 3) ) + 
  facet_wrap(~initialTime>60,ncol = 2) + theme_minimal() + labs(x = "Time, seconds", y = "Hazard of quitting")
ggsave("10hz_smoothed_piecewise_hazard_over_time_group.pdf", width = 6, height = 4)


'summary(m9_2 <- lme4::glmer(quit ~ scale(latency) + scale(sv)*scale(initialTime)*Group +
                              scale(sv_lag2s)*scale(initialTime)*Group +
                              (1|ID), tddf2, family = binomial,
                            glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))

summary(m9_2x <- lme4::glmer(quit ~ scale(latency) + scale(sv)*scale(initialTime)*Group +
                              scale(sv_lag2s)*scale(initialTime)*Group +
                               scale(sv)*scale(initialTime)*scale(EXITtot) + 
                               scale(sv_lag2s)*scale(initialTime)*scale(EXITtot) + 
                              (1|ID), tddf2, family = binomial,
                            glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))


# plot coefficients

s <- summary(m5_2)
coef <- s$coefficients
terms1 <- labels(coef)[[1]]
# p1 <- 
p7 <- plot_model(m5_2,  p.kr = T, terms = terms1, order.terms = c(15:17,8,9:11,3,12:14,4,5:7,2),
                 # show.p = TRUE, show.values = TRUE,  group.terms = c(1,1,1,1,2,1),vline.color = "slategray3",
                 show.p = T, show.values = T,  vline.color = "slategray3",axis.lim = c(.7,2),
                 axis.title = "Wait  <    -    >  Quit", value.offset = 0.4,colors = c( "gray47", "red3", "green4", "navy"),
                 title = "Suicide attempters fail to learn optimal timing")
p8 <- plot_model(m5_2,  p.kr = T, terms = terms1, order.terms = c(15:17,8,9:11,3,12:14,4,5:7,2),
                 group.terms = c(1,2,4,1,1,1,3,2,2,2,4,4,4,3,3,3),
                 show.p = T, show.values = T,  vline.color = "slategray3",axis.lim = c(.7,2),
                 axis.title = "Wait  <    -    >  Quit", value.offset = 0.4,colors = c( "gray47", "red3", "green4", "navy"),
                 title = "Effects of subjective value (SV) and learning (SV*time), reference group: attempters")
ggsave("discrete_time_wtw_m5_2.pdf", width = 10, height = 6)


'plot_model(m5_2)

# summary(m5a_2 <- lme4::glmer(quit ~ scale(latency)*Group + scale(sv)*scale(log(initialTime))*Group + (1|ID), tddf2, family = binomial,
                            # glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))
save(list = ls(all.names = TRUE), file = "wtw_discrete_results.RData")
 load(file = "wtw_discrete_results.RData")

