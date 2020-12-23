# import, check and prepare for analyses data from two bandit samples

setwd("~/Box/Project_wtw/R")
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
trial_df <-
  read_csv("wtw_df1.csv")

#View(trial_df)
# sub_df <- read_csv("wtw_df2.csv")

sub_df <- read_excel("WTW DATA 05-15-18.xlsx")


# View(sub_df)
sub_df = sub_df %>% as_tibble %>% arrange(ID)

# sub_df <- sub_df[,c("ID", "BASELINEAGE", "GROUP1245", "GROUP12467", "GENDERTEXT", "EDUCATION", "RACETEXT", "MARITALTEXT")]

table(sub_df$group1245)
table(sub_df$group12467)

# just missing one education
missing_ind_chars = aggr(
  sub_df,
  col = mdc(1:2),
  numbers = TRUE,
  sortVars = TRUE,
  labels = names(sub_df),
  cex.axis = .7,
  gap = 3,
  ylab = c("Proportion of missingness", "Missingness Pattern")
)


sub_df[,c("ID",  "group1245", "group12467", "sex", "race", "marital_status", "SubstanceLifetime", "SubstanceCurrent", "AnxietyLifetime", "AnxietyCurrent", "ANTIPSYCHOTICS", "SEDATIVES", "OPIATES")] <- 
  lapply(sub_df[,c("ID",  "group1245", "group12467", "sex", "race", "marital_status", "SubstanceLifetime", "SubstanceCurrent", "AnxietyLifetime", "AnxietyCurrent", "ANTIPSYCHOTICS", "SEDATIVES", "OPIATES")], factor)


sub_df$Group <-
  dplyr::recode(
    sub_df$group1245,
    `1` = "Controls",
    `2` = "Depressed",
    `4` = "Ideators",
    `5` = "Attempters"
  )
contrasts(sub_df$Group) <- contr.treatment(levels(sub_df$Group),
                                           base = which(levels(sub_df$Group) == 'Attempters'))

# recode WTAR misscode
sub_df$WTARSS[sub_df$WTARSS>200] <- NA

# mark subjects who never moved the cursor
bad <- c(45709, 207789, 211245, 212926, 213704, 214542, 215842, 217283, 217293, 219675, 220244, 220947, 221273, 881137)

sub_df$bad <- is.element(sub_df$ID, bad)
library(formattable)
formattable(table(sub_df$Group,sub_df$bad))
# sample characteristics: attempters less educated, HC a bit older
chars <- as.data.frame(sub_df[!sub_df$bad,c("sex", "race", "marital_status", "age", "education", "WTARSS","EXITtot", 
"UPPSPNEGURGENCY","UPPSPPOSURGENCY","UPPSPLACKOFPREMED","UPPSPLACKOFPERSEV" , "NONPLAN",        
"SubstanceLifetime", "AnxietyLifetime")])
c1 <-
  compareGroups(
    chars,
    y = sub_df$Group[!sub_df$bad],
    # bivar = TRUE,
    include.miss = FALSE
  )
t1 <-
  createTable(
    c1,
    # hide = c(sex = "FEMALE", list(race = c(
    #   "WHITE", "ASIAN PACIFIC"
    # ))),
    hide.no = 0,
    digits = 1,
    show.n = TRUE
    ,show.p.mul = TRUE
  )
export2html(t1, "wtw_by_group_mult.html")

# merge trial-by-trial and subject-level data
df <- merge(trial_df, sub_df)

summary(df)
# number trials
df = df %>% group_by(ID) %>% 
  dplyr::mutate(
    trial = 1:n()
  ) %>% ungroup()

# inspect:  looks like we have one subjects who had a different contingency
hist(df$latency[df$trialResult=='win'],100)
unique(df$ID[df$trialResult=='win' & df$latency>4 & df$latency<20])

max(df$outcomeTime[df$ID==207989])


# basic qc
# pdf(file = "wtw_qc.pdf", width = 20, height = 20)
# ggplot(df, aes(x = outcomeTime, y = latency)) + geom_line() + facet_wrap(~ID)
# dev.off()


###########################
# preprocess, compute vars

# df$Group <-
#   dplyr::recode(
#     df$GROUP1245,
#     `1` = "Controls",
#     `2` = "Depressed",
#     `4` = "Ideators",
#     `5` = "Attempters"
#   )
# contrasts(df$Group) <- contr.treatment(levels(df$Group),
#                                         base = which(levels(df$Group) == 'Attempters'))

df$GroupLeth <-
  dplyr::recode(
    df$group12467,
    `1` = "Controls",
    `2` = "Depressed",
    `4` = "Ideators",
    `6` = "LL Attempters",
    `7` = "HL Attempters"
  )
contrasts(df$GroupLeth) <-
  contr.treatment(levels(df$GroupLeth),
                  base = which(levels(df$GroupLeth) == 'HL Attempters'))

df$win <- as.factor(df$trialResult=="win")
df$immQuit <- as.factor(df$initialPos=="optSmall" & df$latency<0.1)


# df$outcome_type <- NA
# df$outcome_type[df$immQuit] <- 'immQuit'
# df$outcome_type[!df$immQuit & df$win] <- 'immQuit'

# reward rate


df = df %>% as_tibble %>% arrange(ID, trial)
df = df %>% arrange(ID, trial) %>% group_by(ID) %>% 
  mutate(
    win_lag = lag(win),
    immQuit_lag = lag(immQuit),
    lat_lag = lag(latency),
    initialTime_lead = lead(initialTime),
    totalEarned_lag = lag(totalEarned),
    designatedWait_lag = lag(designatedWait)
        ) %>% ungroup()


# identify non-quitters
df = df %>% arrange(ID, trial) %>% group_by(ID) %>% 
  mutate(
    Nquits = sum(trialResult=='quit' & latency>.1)
  ) %>% ungroup()


df$rewardRate <- df$payoff/(df$initialTime_lead - df$initialTime)
# hist(df$rewardRate, 100)
df = df %>% arrange(ID, trial) %>% group_by(ID) %>% 
  mutate(
    rewardRate_lag = lag(rewardRate)
  ) %>% ungroup()


lag_test <- df[,c(1:9, 20:28)]
View(lag_test)

df$period <- NA
df$period[df$outcomeTime<100] <- 'early'
df$period[df$outcomeTime>100 & df$outcomeTime<200] <- 'middle'
df$period[df$outcomeTime>200] <- 'late'
df$period <- as.factor(df$period)
levels(df$period) <- c("early", "middle", "late")


# make censoring and event variable
df$cens <- df$trialResult=='win' 
df$quit <- !df$trialResult=='win' 

# remove one abnormally long latency
df$latency[df$latency>21.1] <- NA

# define end of interval
df$time2 <- 20.1; 

# split time into early/late learning
df$period <- NA
df$period[df$initialTime<100] <- "first 100s"
df$period[df$initialTime>100] <- "last 200s"

# only non-immediate-quits
ddf <- df[!as.logical(df$immQuit),]

# exclude trials after an immediate quit to get rid of immediate quit/win corrrelation
lagddf <- df[!as.logical(df$immQuit_lag),]

# exclude immediate quits AND trials after an immediate quit
ddf1 <- df[!as.logical(df$immQuit_lag) & !as.logical(df$immQuit),]


# only delayed quits reflecting what the subject has learned for sure
qdf <- df[!as.logical(df$immQuit) & !as.logical(df$win),]

# delayed quits and 20s wins (actual WTW)
wdf <- df[!as.logical(df$immQuit) & (!as.logical(df$win) | df$latency >19.9),]

# final earnings
fdf <- df %>%
  group_by(ID) %>%
  arrange(initialTime) %>%
  slice(n()) %>%
  ungroup
fdf$finalEarnings <- fdf$totalEarned/fdf$outcomeTime

# sanity check on latency distributions
# ggplot(df,aes(latency)) + facet_wrap(~trialResult) + geom_histogram()
# ggplot(ddf,aes(latency)) + facet_wrap(~trialResult) + geom_histogram()


# incorporate value of waiting as a time-dependent covariate

# load Joe's value functions
load(file = "vFunc.RData")
v <- as_tibble(v)
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
