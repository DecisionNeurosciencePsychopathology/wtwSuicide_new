# our reinfocement learning generative models simulate adapative persistence behavior as wait-or-quit choices 
# QL1: Q-learning model with a single learning rate
# QL2: Q-learning model with separate learning rates for rewards and non-rewards
# RL1: R-learning model with a single learning rate 
# RL2: R-learning model with separate learning rates for rewards and non-rewards

# inputs:
# paras: learning parameters
# Rs: trialwise earnings 
# Ts: trialwise delays (measured from the ITI onset)
# normaResults: normative analysis outputs 

# outputs
# trialNum : [nTrialx1 int] 1 : nTrial
# Qwaits_ : [20/40 x nTrial num] value of waiting at each second in each trial
# V0_ : [nTrialx1 num] value of entering a pre-trial iti, namely t = 0

QL1 = function(paras, Rs, Ts, nMadeActions, normResults){
  # default settings 
  iti = 2
  
  # normative analysis 
  optimRewardRates = normResults$optimRewardRates
  optimWaitThresholds = normResults$optimWaitThresholds
  
  # learning parameters
  alpha = paras[1]; tau = paras[2]; gamma = paras[3]; prior = paras[4]
  
  # num of trials
  nTrial = length(Rs) 
  
  # duration of a sampling interval 
  stepSec = 1
    
  # initialize action values 
  V0 = mean(unlist(optimRewardRates)) / (1 - 0.85) # state value for t = 0
  tWaits = seq(iti, delayMax + iti, by = stepSec) # decision points 
  tMax = max(tWaits) #  time point for the last decision point
  Qwaits = -0.1 * (tWaits) + prior + V0 + smallReward # value of waiting at each decision points
  
  # initialize output variables
  Qwaits_ = matrix(NA, length(tWaits), nTrial); Qwaits_[,1] = Qwaits 
  V0_ = vector(length = nTrial); V0_[1] = V0
  
  # log likelihood
  logLik_ = vector(length = nTrial, mode = "list")
    
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # current R and T
    T = Ts[tIdx]
    R = Rs[tIdx]
    nMadeAction = nMadeActions[tIdx]
    
    # update log liklihood
    thisLogLik = rep(NA, length = length(tWaits))
    for(i in 1 : nMadeAction){
      if(R == smallReward & i == nMadeAction){
        thisLogLik[i] = exp((V0 + smallReward) * tau) / (exp((V0 + smallReward) * tau) + exp(Qwaits[i] * tau))
      }else{
        thisLogLik[i] = log(exp(Qwaits[i] * tau) / (exp((V0 + smallReward) * tau) + exp(Qwaits[i] * tau)))
      }
    }
    logLik_[[tIdx]] = thisLogLik
  
    # when the trial endes, update value functions for all time points before T in the trial
    if(tIdx < nTrial){
      
      # calculate expected returns for t >= 2
      Gts =  gamma ^ (T - tWaits) * (R + V0)
      # only update value functions before time t = T
      updateFilter = tWaits <= T 
      # update Qwaits
      Qwaits[updateFilter] = Qwaits[updateFilter] + alpha * (Gts[updateFilter] - Qwaits[updateFilter])
      
      # calculate expected returns for t == 0 and update V0
      Gt =  gamma ^ T * (R + V0)
      V0 = V0 + alpha * (Gt - V0)
      
      # record updated values
      Qwaits_[,tIdx + 1] = Qwaits
      V0_[tIdx + 1] = V0
    }# end of the loop within a trial 
  }  
  
  # calculate the total log likelihood 
  tempt = unlist(logLik_)
  totalLogLik = sum(tempt[!is.na(tempt)])
  # return outputs
  outputs = list( 
    "trialNum" = 1 : nTrial, 
    "Qwaits_" = Qwaits_, 
    "V0_" = V0_,
    "totalLogLik" = totalLogLik
  )
  return(outputs)
}