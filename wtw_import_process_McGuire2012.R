# import, check and prepare for analyses data from two bandit samples

setwd("~/Box Sync/Project_wtw/R")
# setwd("~/code/wtw/data/")
# setwd("../data/")

library(readr)
library(lme4)
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
# library(lsmeans)
library(emmeans)
# library(factoextra)
# library(ggfortify)
library(compareGroups)
# library(RColorBrewer)
library(MASS)
library(effects)
library(readr)
library(VIM)
library(mice)
library(multcompView)
library(stargazer)
library(dplyr)
library(readxl)
# detach(package:plyr)
###################
# merge and check
df <- read_csv("wtw_df1_fMRI.csv")

summary(df)
# number trials
df = df %>% group_by(ID) %>% 
  dplyr::mutate(
    trial = 1:n()
  ) %>% ungroup()

# inspect:  looks like we have one subjects who had a different contingency
hist(df$latency[df$trialResult=='high'],100)

# basic qc
pdf(file = "wtw_qc_fMRI.pdf", width = 20, height = 20)
ggplot(df, aes(x = outcomeTime, y = latency)) + geom_line() + facet_wrap(~ID)
dev.off()


###########################
# preprocess, compute vars

df$win <- as.factor(df$trialResult=="high")
# df$immQuit <- as.factor(df$initialPos=="optSmall" & df$latency<0.1)


# df$outcome_type <- NA
# df$outcome_type[df$immQuit] <- 'immQuit'
# df$outcome_type[!df$immQuit & df$win] <- 'immQuit'

# reward rate


df = df %>% as_tibble %>% arrange(ID, trial)

# identify non-quitters
df = df %>% arrange(ID, trial) %>% group_by(ID) %>% 
  mutate(
    Nquits = sum(trialResult=='quit' & latency>.1)
  ) %>% ungroup()



df$cens <- df$trialResult=='high' 
df$quit <- !df$trialResult=='high' 

# censor at 40s
df$latency[df$latency>41] <- NA
# and <1s
df$latency[df$latency<1] <- NA


# define end of interval
df$time2 <- 41; 

#not worrying about immediate quits
df$trialID <- as.numeric(factor(unique(df$ID)))*1000 + df$trial


# tGrid <- sort(unique(ddf$latency))  # this is insane and breaks Joe's code
# convert to counting process (trial-period) form
tdf <- survSplit(Surv(latency, quit) ~., df, cut = (1:80)/2)
tdf <- as.tibble(tdf)

# write values as time-varying covariate
tdf = tdf %>% as_tibble %>% arrange(trialID, latency)

# interactions with time
tdf$t2hz <- round(tdf$latency*2)/2
# only non-immediate-quits

save(list = ls(all.names = TRUE), file = "wtw_fMRI.RData")
