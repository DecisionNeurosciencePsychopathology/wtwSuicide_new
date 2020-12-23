# import, check and prepare for analyses data from two bandit samples

setwd("~/Box/Project_wtw/R")

library(tidyverse)

# load Joe's value functions
load(file = "vFunc.RData")
# downsample to 1s bins
vs <- as_tibble(v) %>% mutate(t1s = floor(t1)) %>% 
  group_by(t1s) %>% summarise(svs = mean(sv)) %>% ungroup() %>%
  mutate(t2s = t1s+1)



v# downsample for Aliona

# make a new dataframe following Therneau, Crowson & Atkinson, 2018
ids <- unique(df$ID)
# unique trial ID
df$trialID <- df$ID*1000 + df$trial

# tGrid <- sort(unique(ddf$latency))  # this is insane and breaks Joe's code
# convert to counting process (trial-period) form
tdf <- survSplit(Surv(latency, quit) ~., df, cut=v$t2)
tdf <- as.tibble(tdf)

# write values as time-varying covariate
v$tstart <- v$t1


tdf <- merge(tdf,v)
tdf = tdf %>% as_tibble %>% arrange(trialID, latency)

# interactions with time
tdf$t2hz <- round(tdf$latency*2)/2
# only non-immediate-quits
tddf <- tdf[!as.logical(tdf$immQuit),]

# exclude trials after an immediate quit to get rid of immediate quit/win corrrelation
tlagddf <- tdf[!as.logical(tdf$immQuit_lag),]

# exclude immediate quits AND trials after an immediate quit
tddf1 <- tdf[!as.logical(tdf$immQuit_lag) & !as.logical(tdf$immQuit),]

# only delayed quits reflecting what the subject has learned for sure
tqdf <- tdf[!as.logical(tdf$immQuit) & !as.logical(tdf$win),]

# delayed quits and 20s wins (actual WTW)
twdf <- tdf[!as.logical(tdf$immQuit) & (!as.logical(tdf$win) | df$latency >19.9),]

# make a downsampled, 2hz version for faster analyses

v2 <- v[is.element(v$t1,c((0:40)/2)),]
tdf2 <- survSplit(Surv(latency, quit) ~., df, cut=v2$t1)
tdf2 <- as.tibble(tdf2)
tdf2 <- merge(tdf2,v2)
tdf2 = tdf2 %>% as_tibble %>% arrange(trialID, latency)
tddf2 <- tdf2[!as.logical(tdf2$immQuit),]

# add lags of SV
tddf2 = tddf2 %>% arrange(ID, trial,latency) %>% group_by(ID,trial) %>% 
  mutate(
    sv_lag1s = lag(sv,2),
    sv_lag2s = lag(sv,4),
    sv_lag3s = lag(sv,6),
    sv_lag4s = lag(sv,8)
  ) %>% ungroup()




# sanity check
ggplot(tdf,aes(x = latency, y = sv)) + geom_line(lty = 'dotted')
ggplot(tdf,aes(x = latency, y = waitValue)) + geom_line()
pdf(file = 'hazard by group and period.pdf')
ggplot(tddf[!is.na(tddf$initialTime),],aes(x = latency, y = waitValue/100)) + geom_line() + 
  geom_line(aes(y = sv), lty = 'dotted') + 
  geom_smooth(aes(y = quit*500, color = Group), se = F) + facet_wrap(~period)
# turned off SE because they cannot be trusted with GAM
dev.off()
# binomial_smooth <- function(...) {
#   geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
# }
# spline fit does not work well with these data
# ggplot(tddf[!is.na(tddf$initialTime),],aes(x = latency, y = quit, color = Group)) + 
#   geom_smooth(method = "glm", method.args = list(family = "binomial"),formula = y ~ splines::ns(x, 4)) + facet_wrap(~period)


save(list = ls(all.names = TRUE), file = "wtw.RData")
