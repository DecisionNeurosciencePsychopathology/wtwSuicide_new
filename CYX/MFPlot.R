library("tidyverse")
library("ggpubr")
library("latex2exp")
source("subFxs/loadFxs.R")
source("subFxs/plotThemes.R")
source('expSchematics.R')
source('MFAnalysis.R')

# output dir
dir.create('figures/MFplot')

# load experiment parameters
load("expParas.RData")

# normative analysis
iti = 2
normResults = expSchematics(iti, F)
optimWaitThresholds = normResults$optimWaitThresholds

# load data
allData = loadAllData() # partipants who never moved the cursor has been excluded
trialData = allData$trialData  
hdrData = allData$hdrData

######################### plot WTW timecourses in two environments ######################
MFResults = MFAnalysis()
sumStats = MFResults[['sumStats']]
sumStats$group = factor(hdrData$group, levels = c("Controls", "Depressed", "Attempters", "Ideators"))
timeWTW_ = MFResults[['timeWTW_']]
nSub = nrow(sumStats)
## use background color to distinguish used and excluded data 
data.frame(wtw = unlist(timeWTW_),
           time = rep(seq(0, blockSec * nBlock -1, by = 1), nSub),
           group = rep(hdrData$group, each = length(tGrid))) %>%
  group_by(time, group) %>%
  dplyr::summarise(mu = mean(wtw, na.rm = F), se = sd(wtw, na.rm = F) / sqrt(sum(!is.na(wtw))),
                   min = mu- se, max = mu + se) %>%
  ggplot(aes(time, mu, color = group))  + geom_line(size = 0.5) +
  geom_ribbon(aes(ymin=min, ymax=max, fill = group, color = NULL), alpha = 0.2) +
  xlab("Task time (s)") + ylab("WTW (s)") + myTheme +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) +
  theme(legend.position = "none")
ggsave("figures/MFPlot/wtw_timecourse.eps", width = 5, height = 3) 
ggsave("figures/MFPlot/wtw_timecourse.png", width = 5, height = 3) 

# plot AUC by group
sumStats %>% mutate(group = factor(group, level = c("Controls", "Depressed", "Ideators", "Attempters"))) %>%
  ggplot(aes(group, muWTW)) + geom_boxplot() + myTheme + ylab("AUC (s)") + xlab("")
ggsave("figures/MFPlot/auc_group.eps", width = 5, height = 3) 
ggsave("figures/MFPlot/auc_group.png", width = 5, height = 3) 

# check the time later
################## plot AUCs in two environments ###################
## controls vs non-controls, sig
## controls vs attempters, sig
wTest = wilcox.test( sumStats[sumStats$group == "Controls", "muWTW"],
                     sumStats[sumStats$group == "Attempters", "muWTW"],paired = F)
sumStats %>% group_by(group) %>% summarise(median = median(muWTW))
##  not significant
aovRes <- aov(muWTW ~ (group), contrasts = contr.sum, data = sumStats)
summary(aovRes)
sumStats %>% ggplot(aes(x = group, y = muWTW)) + 
  geom_dotplot(binaxis='y', stackdir='center', aes(fill = group)) +
  xlab("") + ylab("AUC (s)") +
  myTheme + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) +
  theme(legend.position = "none")
ggsave("figures/MFPlot/muWTW_comparison.eps", width = 6, height = 4)
ggsave("figures/MFPlot/muWTW_comparison.png", width = 6, height = 4)


################### plot CIPs in two environments ###################
# test
sumStats %>% group_by(group) %>% summarise(median = median(stdWTW))
wTest = wilcox.test( sumStats[sumStats$group == "Controls", "stdWTW"],
                     sumStats[sumStats$group == "Attempters", "stdWTW"],paired = F)
aovRes <- aov(stdWTW ~ group, contrasts = contr.sum, data = sumStats)
summary(aovRes)
# plot
sumStats %>% ggplot(aes(x = group, y = stdWTW)) + 
  geom_dotplot(binaxis='y', stackdir='center', aes(fill = group)) +
  xlab("") + ylab(TeX("CIP (s^2)")) +
  myTheme + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) +
  theme(legend.position = "none")
ggsave("figures/MFPlot/stdWTW_comparison.eps", width = 6, height = 4)
ggsave("figures/MFPlot/stdWTW_comparison.png", width = 6, height = 4)


# plot CIP by group
sumStats %>% mutate(group = factor(group, level = c("Controls", "Depressed", "Ideators", "Attempters"))) %>%
  ggplot(aes(group, stdWTW)) + geom_boxplot() + myTheme + ylab(TeX("CIP (s^2)")) + xlab("")
ggsave("figures/MFPlot/CIP_group.eps", width = 5, height = 3) 
ggsave("figures/MFPlot/CIP_group.png", width = 5, height = 3) 


################################ plot survival curves #####################
# stats test
survCurve_ = MFResults$survCurve_
plotData = data.frame(survCurve = unlist(survCurve_),
                      time = rep(kmGrid, nSub),
                      group = rep(hdrData$group, each = length(kmGrid))) %>% 
  group_by(group, time) %>%
  dplyr::summarise(mu = mean(survCurve, na.rm = F), se = sd(survCurve, na.rm = F) / sqrt(sum(!is.na(survCurve))),
          min = mu- se, max = mu + se)

plotData %>%
  ggplot(aes(time, mu, fill = group)) + geom_line(aes(color = group)) +
  geom_ribbon(aes(time, ymin = min, ymax = max, color = NAN), alpha = 0.2, color = NA) +
  xlab("Elapsed time (s)") + ylab("Survival rate") + myTheme +
  theme(legend.position = "none") +
  geom_vline(aes(xintercept = 3))
ggsave("figures/MFPlot/survival_curve.eps", width = 4, height = 4)
ggsave("figures/MFPlot/survival_curve.png", width = 4, height = 4) 


########################### plot wtw timecourse ##################
timeWTW_ = MFResults$timeWTW_
plotData = data.frame(wtw = unlist(timeWTW_),
                      time = rep(tGrid, nSub),
                      group = rep(hdrData$group, each = length(tGrid))) %>% 
  group_by(group, time) %>%
  dplyr::summarise(mu = mean(wtw, na.rm = F), se = sd(wtw, na.rm = F) / sqrt(sum(!is.na(wtw))), min = mu- se, max = mu + se) %>%
  ungroup()

plotData %>%
  mutate(group = factor(group, level = c("Controls", "Depressed", "Ideators", "Attempters"))) %>%
  ggplot(aes(time, mu, fill = group)) + geom_line(aes(color = group)) +
  geom_ribbon(aes(time, ymin = min, ymax = max, color = NAN), alpha = 0.2, color = NA) +
  xlab("Task time (s)") + ylab("WTW (s)") + myTheme 
ggsave("figures/MFPlot/WTW.eps", width = 6, height = 4)
ggsave("figures/MFPlot/WTW.png", width = 6, height = 4)

