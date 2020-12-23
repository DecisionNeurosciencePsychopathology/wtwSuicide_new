load("expParas.RData")
library(latex2exp)
library("ggplot2"); library("Hmisc"); library("coin")
library("dplyr"); library("tidyr")
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R") # load blockData and expPara
source("subFxs/helpFxs.R") # getparaNames
source("subFxs/analysisFxs.R") # plotCorrelation and getCorrelation
source('MFAnalysis.R')
library("ggpubr")

# model Name
modelName = "QL1"
paraNames = getParaNames(modelName)
nPara = length(paraNames)

# output directories
dir.create("figures/expParaAnalysis")
saveDir = sprintf("figures/expParaAnalysis/%s", modelName)
dir.create(saveDir)

# 
MFResults = MFAnalysis()
sumStats = MFResults[['sumStats']]
# load expPara
paraNames = getParaNames(modelName)
nPara = length(paraNames)
parentDir = "genData/expModelFit"
dirName = sprintf("%s/%s",parentDir, modelName)
expPara = loadExpPara(paraNames, dirName)
passCheck = checkFit(paraNames, expPara)

# load exp data
allData = loadAllData()
hdrData = allData$hdrData 
trialData = allData$trialData       
ids = hdrData$ID
nSub = length(ids)


# plot hist 
# paraNames = c("LR", "LP", expression(tau), expression(gamma), "P")
# paraNames = c("LR", "LP", expression(tau), "P")
paraNames = paraNames
expPara$condition = sumStats$condition[1 : nrow(expPara)]
plotData = expPara %>% filter(passCheck ) %>% select(c(paraNames, "condition")) %>%
  gather(-c("condition"), key = "para", value = "value") 
plotData$para = factor(recode(plotData$para, alpha = "Learning Rate",
              gamma = "Discounting",
              tau = "Inverse Temp",
              prior = "Prior WTW"),
       level = c("Learning Rate", "Inverse Temp", "Discounting", "Prior WTW"))

plotData %>%
  ggplot(aes(value)) + geom_histogram(bins = 8) +
  facet_grid(condition ~ para, scales = "free") + 
  myTheme + xlab(" ") + ylab(" ")
fileName = sprintf("%s/%s/hist.pdf", "figures/expParaAnalysis", modelName)
ggsave(fileName, width = 8, height = 4)


# 
data = cbind(hdrData, expPara)
# data$alpha = log(data$alpha)
# data$tau = log(data$tau)
plotData = data %>% select(c("ID", "group", "alpha", "tau", "gamma", "prior")) %>%
  filter(passCheck) %>%
  mutate(group = factor(group, level = c("Controls", "Depressed", "Ideators", "Attempters"))) %>%
  gather("para", "value", -c("ID", "group")) 

plotData$para = factor(recode(plotData$para, alpha = "Learning Rate",
                       gamma = "Discounting",
                       tau = "Inverse Temp",
                       prior = "Prior WTW"),
                       level = c("Learning Rate", "Inverse Temp", "Discounting", "Prior WTW"))
p = ggboxplot(plotData, x = "group", y = "value",
                       color = "group", 
                       add = "jitter", facet.by = "para", scales = "free")
my_comparisons <- list( c("Controls", "Depressed"), c("Controls", "Ideators"), c("Controls", "Attempters") )
p + stat_compare_means(comparisons = my_comparisons, paired = F, method = "wilcox")
fileName = sprintf("%s/%s/groupDiff.pdf", "figures/expParaAnalysis", modelName)
ggsave(fileName, width = 12, height = 8)


##### plot correlations 
data$auc = sumStats$muWTW
data$cip = sumStats$stdWTW
plotData = data %>% select(c("ID", "group", "alpha", "tau", "gamma", "prior", "auc", "cip")) %>%
  filter(passCheck) %>%
  mutate(group = factor(group, level = c("Controls", "Depressed", "Ideators", "Attempters"))) 
names(plotData) = c("ID", "group", "Learning Rate", "Inverse Temp", "Discount", "Prior", "AUC", "CIP")
corMx = cor(plotData[, 3 : ncol(plotData)])
library("corrplot")
p = corrplot(corMx, method="number", type="upper")
fileName = sprintf("%s/%s/corr.pdf", "figures/expParaAnalysis", modelName)
ggsave(fileName, p, width = 6, height = 6)


