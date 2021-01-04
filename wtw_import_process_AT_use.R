# import, check and prepare for analyses data from two bandit samples

#setwd("~/Box/Project_wtw/R_2020")--save files here
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
library(tidyverse)
library(stringr)
library(readr)
library(coxme)
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
    ##bivar = TRUE,
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

# merge trial-by-trial (wtw_df1.csv) and subject-level (WTW DATA 05-15-18.xlsx) data
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

#####################################################################################
# incorporate value of waiting as a time-dependent covariate

# load Joe's value functions
load(file = "vFunc.RData")
v <- as_tibble(v)

###AT: SAVE THE FILES
setwd("~/Box/Project_wtw/wtwSuicide_new-master/WTW_files.Rdata")
save(file = "WTW_files.Rdata",list = ls(all = T))

##START HERE
###AT: LOAD ALL THE FILES GENERATED FROM THE STEPS ABOVE:
load("~/Box/Project_wtw/wtwSuicide_new-master/WTW_files.Rdata")

##AT: LOAD YIXIN'S DATA, create multiple rows (trials) per person dataset 
setwd("~/Box/Project_wtw/wtwSuicide_new-master/CYX/genData/medianParaValueFun")

library(tidyverse)
library(stringr)
library(readr)

list_of_files <- list.files(path = "~/Box/Project_wtw/wtwSuicide_new-master/CYX/genData/medianParaValueFun", pattern = "*.csv",
                full.names=T)

df_Yvalues <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read_csv(.x, col_types = cols(), col_names = FALSE), .id = "file_name")  

####Delete every first row (1st of 21 rows) per participant---dropping the first row as 
##the values are aligned with timepoints rather than bins, no such thing as waiting until second 0

 ## Delete every 21st row starting from 21
df_Yvalues<-df_Yvalues[-seq(21,NROW(df_Yvalues), by = 21),]

# Add t2 variable 
df_Yvalues = df_Yvalues %>% group_by(file_name) %>% 
  dplyr::mutate(
    t2s = 1:n()
  ) %>% ungroup()

# Add t1 variable 
df_Yvalues = df_Yvalues %>% group_by(file_name) %>% 
  dplyr::mutate(
    t1s = 1:n()-1
  ) %>% ungroup()

###AT SAVE in r
setwd("~/Box/Project_wtw/wtwSuicide_new-master/")
save(file = 'df_Yvalues.Rda', df_Yvalues)

###AT load Yixin's multiple rows (time bins) and columns (trials) per person dataset with t1 and t2 variables:
load("~/Box/Project_wtw/wtwSuicide_new-master/df_Yvalues.Rda")

dfy <- df_Yvalues
dfy$ID <- as.numeric(unlist(regmatches(dfy$file_name, gregexpr("[[:digit:]]+", dfy$file_name))))
# extract ID

###Creating file with file_name, trials as rows, and values

 dfyl<-dfy %>% 
pivot_longer (!c(file_name, ID, t1s, t2s), names_to = "trial", names_prefix = "X", values_to = "value") 

 #get rid of the file_name column (redundant with ID in there already)
 dfyl=subset(dfyl, select=-c(file_name))
 
 dfyl<-dfyl %>% arrange(ID, trial, t1s)
 
 #sanity check
ggplot(dfyl, aes(t1s, value)) + geom_smooth()

##save this multiple rows per person values data file
setwd("~/Box/Project_wtw/wtwSuicide_new-master/")
save(file = 'dfyl.Rda', dfyl)
########




 # make a new dataframe following Therneau, Crowson & Atkinson, 2018
 ids <- unique(df$ID)
 # unique trial ID
 df$trialID <- df$ID*1000 + df$trial
 
cuts <- unique(dfyl$t2s)
# tGrid <- sort(unique(ddf$latency))  # this is insane and breaks Joe's code
# convert to counting process (trial-period) form
tdf <- survSplit(Surv(latency, quit) ~., df, cut=cuts)
tdf <- as_tibble(tdf)

# write values as time-varying covariate
dfyl$tstart <- dfyl$t1s


tdf <- as_tibble(merge(tdf,dfyl))

tdf = tdf %>% arrange(trialID, tstart)
ggplot(tdf, aes(tstart, value)) + geom_smooth()

###############

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

# sanity check
ggplot(tdf,aes(x = latency, y = value)) + geom_smooth()

###AT: SAVE ALL THEse FILES
setwd("~/Box/Project_wtw/wtwSuicide_new-master/WTW_files_all.Rdata")
save(file = "WTW_files_all.Rdata",list = ls(all = T))

#####AT: Separately also save tddf (only non-immediate quits) and tddf1 (exclude immediate quits AND trials after an immediate quit):
setwd("~/Box/Project_wtw/wtwSuicide_new-master/")
save(file = 'tddf.Rda', tddf)
save(file = 'tddf1.Rda', tddf1)


#####
# downsample Joe's values to 1s bins
vs <- as_tibble(v) %>% mutate(t1s = floor(t1)) %>%
  group_by(t1s) %>% summarise(svs = mean(sv)) %>% ungroup() %>%
  mutate(t2s = t1s+1)

###Merge tddf and vs files on t1s, to add svs (ideal values) to the tddf
# merge two data frames in r
# r merge by rownames---CHECK
tddf_vs <- merge(tddf, vs, by = 't1s')

##check that these are the same variables (they are)
plot(tddf_vs$t2s.x, tddf_vs$t2s.y)

save(file = 'tddf_vs.Rda', tddf_vs)
##plot Joe's idealized values

ggplot(tddf_vs, aes(t1s, svs)) + geom_line()


##see how Joe's idealized values align with Yixin's values
ggplot(tddf_vs, aes(t1s, svs)) + geom_line() + geom_smooth(aes(t1s, value))

ggplot(tddf_vs, aes(t1s, svs)) + geom_line() + geom_smooth(aes(t1s, value*2 - 4)) + geom_smooth(aes(t1s, as.numeric(quit)*40, color= "r"))


# pdf(file = 'hazard by group and period.pdf')
# ggplot(tddf[!is.na(tddf$initialTime),],aes(x = latency, y = waitValue/100)) + geom_line() + 
#   geom_line(aes(y = sv), lty = 'dotted') + 
#   geom_smooth(aes(y = quit*500, color = Group), se = F) + facet_wrap(~period)
# # turned off SE because they cannot be trusted with GAM
# dev.off()
# binomial_smooth <- function(...) {
#   geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
# }
# spline fit does not work well with these data
# ggplot(tddf[!is.na(tddf$initialTime),],aes(x = latency, y = quit, color = Group)) + 
#   geom_smooth(method = "glm", method.args = list(family = "binomial"),formula = y ~ splines::ns(x, 4)) + facet_wrap(~period)

# try a coxme for validation
library(coxme)

## Yixin's values--effects of value negative (z= - 5.71)
## Data: tddf_vs, events, n = 1423, 29152 (16928 observations deleted due to missingness)

c2 <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(value) + (1|ID), tddf_vs)
summary(c2)

##Joe's idealized values, effects of value slightly negative but not significantly (z= -0.08)
#Data: tddf_vs events, n = 2302, 44657 (1423 observations deleted due to missingness)
c3 <- coxme(Surv(latency, quit) ~ scale(lat_lag) + scale(initialTime) + win_lag + scale(svs) + (1|ID), tddf_vs)
summary(c3)

##Older files (Joe's value): Data: tddf events, n = 1847, 304328 (9757 observations deleted due to missingness)

##########################################
##########################################
##########################################
##########################################
##From the code wtw_coxme.R:

# make a na.omit df for plotting predicted values
new_tddf <- na.omit(tddf[,c("ID", "latency", "quit","lat_lag", "initialTime", "win_lag", "value", "Group")])

cf0 <- coxph(Surv(latency,quit)~1, new_tddf)
# plot(cf0)
c0 <- coxme(Surv(latency, quit) ~ (1|ID), new_tddf)
new_tddf$predicted_lp_c0 <- predict(c0, type = "risk")
pdf(file = 'predicted risk c0.pdf')
ggplot(new_tddf[new_tddf$latency<4,],aes(x = latency, y = predicted_lp_c0, color = Group)) + geom_smooth() 
# + facet_wrap(~ID)
dev.off()

####c0v <- coxme(Surv(latency, quit) ~ scale(waitValue) + (1|ID), tddf)
##summary(c0v)


# build model systematically
c1a <- coxme(Surv(latency, quit) ~ scale(lat_lag) + (1|ID), new_tddf)
new_tddf$predicted_h_c1 <- predict(c1a, type = "risk")
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

###From wtw_import_process.R line 318
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

###From wtw_import_process.R line 318 MODIFIED FOR YIXIN'S DATA---
ggplot(tdf,aes(x = latency, y = value)) + geom_line(lty = 'dotted')



pdf(file = 'hazard by group and period.pdf')

#instead of waitValue use sv; panel plots by both subject and early/late in learning; print to 30x30 or larger pdf
ggplot(tddf_vs[!is.na(tddf_vs$initialTime),],aes(x = latency, y = svs/100)) + geom_line() + 
  geom_line(aes(y = value), lty = 'dotted') + 
  geom_smooth(aes(y = quit*500, color = Group), se = F) + facet_wrap(~period)

# turned off SE because they cannot be trusted with GAM
dev.off()



