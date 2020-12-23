# this script plots delay distributions and reward rates in two environments
expSchematics = function(iti, isPlot){
# load experiment parameters 
load("expParas.RData")
seqHP = rewardDelays$HP
seqLP = rewardDelays$LP

## for display purposes, all variables on continous time
## are discretized into 0.1 secondition time bins
bin = 0.1 # width of a time bin
time = list(
  LP = seq(bin, delayMax, by = bin)
) 

## delays in both HP and LP are drawn repeatedly from 8 possible values,
## In HP, they are uniformly distributed from 0 to 16, 
## in LP, they are log-spaced from 0 to 32, therefore forming a heavy tail distribution

# delay PDFs
LP = rep(0, length(time$LP))
LP[sapply(1 :4, function(i) which.min(abs(time$LP - seqLP[i])))] = 1 / 4
rewardDelayPDFs =list('LP' = LP)

# delay CDFs
rewardDelayCDFs = list('LP' = cumsum(rewardDelayPDFs$LP))

# average waiting durations given different policies
## Here we assume rewards occur at the middle of each time bin
LP = cumsum((time[['LP']] ) * rewardDelayPDFs$LP) / cumsum(rewardDelayPDFs$LP)
meanRewardDelays = list('LP' = LP)

# rewardRates given different policies
LP = (largeReward * rewardDelayCDFs$LP  + smallReward * (1 - rewardDelayCDFs$LP))/
  ((meanRewardDelays$LP * rewardDelayCDFs$LP + time[['LP']] * (1 - rewardDelayCDFs$LP)) + iti)
rewardRates = list('LP' = LP)


# optimal raward rates and optimal policies
optimWaitThresholds = list()
optimWaitThresholds$LP = time$LP[which.max(LP)]
optimRewardRates = list()
optimRewardRates$LP = max(LP, na.rm = T)

# calculate theoretical value of waiting as a function of elapsed time 
subjectValues = list()
condition = "LP"
pdf = rewardDelayPDFs[[condition]]
cdf = rewardDelayCDFs[[condition]]
thisTime = time[[condition]]
Rstar = optimRewardRates[[condition]] # opportunity cost
ts = seq(0, delayMax, by = 0.1) # elapsed time
# initialize 
thisSubjectValues = vector(length = length(ts)) # waiting value given the elapsed time value
Tstars = vector(length = length(ts))  # optimal waiting policy given the elapsed time value
# loop over different elapsed time values
for(i in 1 : length(ts)){
  t = ts[i] 
  trctTime = thisTime[thisTime > t]
  # given the elapsed time value, t, loop over different waiting policies to find Tstar 
  if(t == delayMax){
    Tstar = t
    gt_max = largeReward
  }else{
    Tstar = t
    gt_max = -100
    for(T in seq(t, delayMax, by = 0.1)){
      trctPDF = pdf[thisTime > t] / sum(pdf[thisTime > t])
      at = largeReward * sum(trctPDF[trctTime <= T]) + smallReward * (1 - sum(trctPDF[trctTime <= T]))
      trctPDF[trctTime == T] = trctPDF[trctTime == T]
      bt = sum((trctTime[trctTime <= T] - 0.5 * bin - t) * trctPDF[trctTime <= T]) + 
        + (T - t) * sum(trctPDF[trctTime > T]) 
      gt = at - bt * Rstar
      if(gt > gt_max){
        gt_max = gt
        Tstar = T
      }
    }
  }
  thisSubjectValues[i] = gt_max 
  Tstars[i] = Tstar
}
subjectValues[[condition]] = thisSubjectValues


  

# plot 
if(isPlot){
  # plot CDFs 
  library('ggplot2')
  source('subFxs/plotThemes.R')
  library("latex2exp")
  library("tidyr"); library('dplyr')
  dir.create("figures")
  dir.create('figures/expSchematics')
  
  data.frame(
    CDF = rewardDelayCDFs$LP,
    time = time$LP
  ) %>% ggplot(aes(time, CDF)) + geom_line(size = 3, color = conditionColors[2]) +
    ylim(c(0,1)) + scale_y_continuous(breaks = c(0,0.5,1)) + 
    scale_x_continuous(breaks = c(0, max(delayMax)/ 2, max(delayMax)),limits = c(0, max(delayMax) * 1.1)) +
    myTheme + xlab('Delay duration (s)') + ylab('CDF') + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  ggsave('figures/expSchematics/CDF.eps', width =4, height = 3)
  ggsave('figures/expSchematics/CDF.png', width =4, height = 3)
  
  # plot reward rates
  data.frame(rewardRate = rewardRates$LP,
             time = time$LP) %>%
    ggplot(aes(time, rewardRate)) + 
    geom_line(size = 2, color = conditionColors[2])  + myTheme +
    ylab(TeX("Reward rate $\\rho_T$ (¢ $s^{-1}$)")) + xlab("Give-up time (s)") +
    theme(plot.title = element_text(hjust = 0.5, color = themeColor)) +
    theme(legend.position = "none")
  ggsave("figures/expSchematics/reward_rate.eps", width = 4, height = 3)
  ggsave("figures/expSchematics/reward_rate.png", width = 4, height = 3)
  
  # plot subjective value of waiting 
  data.frame(
    value = subjectValues$LP,
    t = seq(0, delayMax, by = 0.1)
  ) %>%
    ggplot(aes(t, value)) + geom_line(color = conditionColors[2], size = 2) +
    myTheme  +
    xlab("Elapsed time (s)") + ylab("Subjective value (¢)")  + 
    theme(legend.position = "none")
  ggsave("figures/expSchematics/subjective.eps", width = 4, height = 3)
  ggsave("figures/expSchematics/subjective.png", width = 4, height = 3)       
}

# return outputs 
outputs = list(
  "optimWaitThresholds" = optimWaitThresholds,
  "optimRewardRates" = optimRewardRates,
  "subjectValues" = subjectValues
)
return(outputs)
}
