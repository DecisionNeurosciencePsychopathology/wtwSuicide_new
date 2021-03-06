# fit a reinforcement learning model for a single participant in Rstan 

# inputs:
# thisTrialData: behavioral data for this participant
# fileName: the name of the output file
# modelName: the name of   the model 
# paraNames: parameters for the model
# model: the Bayesian model 
# config: a list containing the Rstab configuration 
modelFitSingle = function(id, thisTrialData, modelName, paraNames, model, config, outputFile){
    # load experiment paras
    load('expParas.RData')
  
    # parse the stan configuration
    nChain = config[['nChain']] # number of MCMC chains
    nIter = config[['nIter']] # number of iterations on each chain
    controlList = list(adapt_delta = config[['adapt_delta']],
                       max_treedepth = config[['max_treedepth']] )
    warningFile = config[['warningFile']] # output file for stan warnings and errors
    
    # duration of one time step (namely one temporal state) 
    stepSec = 1
    
    # duration of iti
    iti = 2
    
    # normative analysis
    normResults = expSchematics(iti, F)
    optimRewardRates = normResults$optimRewardRates
    optimWaitThresholds = normResults$optimWaitThresholds
    subjectValue = normResults$subjectValue
    
    # prepare inputs for fitting the model
    condition = "LP"
    
    # excluded
    immQuit = thisTrialData$initialPos == "optSmall"
    thisTrialData = thisTrialData[!immQuit, ]
    
    ## ensure timeWaited = scheduledWait on rewarded trials
    thisTrialData = within(thisTrialData, {timeWaited[trialEarnings!= smallReward] = scheduledWait[trialEarnings!= smallReward]})
    
    ## wait-or-quit decision time points in a trial
    tWaits = seq(iti, delayMax + iti, by = stepSec)
    ## number of wait-or-quit decision time points in a trial
    nWaitOrQuit = length(tWaits) 
    ## number of made actions in each trial 
    nMadeActions = floor(thisTrialData$timeWaited / stepSec) + 1
    ## when a trial ends 
    Ts = thisTrialData$timeWaited + iti
    Ts[thisTrialData$trialEarnings == smallReward] = floor(Ts[thisTrialData$trialEarnings == smallReward] / stepSec) * stepSec
    ## the theoretic present value of the awaited reward sampled at 1 hz
    LPtheoreticValues = subjectValue$LP[seq(0, delayMax, by = 0.1) %in% seq(0, delayMax, by = stepSec)]
    
    ## orgianze inputs into a list
    inputs <- list(
      iti = iti,
      stepSec = stepSec,
      nWaitOrQuit = nWaitOrQuit,
      tWaits = tWaits,
      N = length(thisTrialData$trialEarnings), # number of trials
      Rs = thisTrialData$trialEarnings, # rewards on each trial
      Ts = Ts,
      nMadeActions = nMadeActions,
      LPvalues = LPtheoreticValues,
      smallReward = smallReward,
      largeReward = largeReward
    )
    if(modelName %in% c("QL1", "QL2")){
      V0_ini = mean(unlist(optimRewardRates)) * stepSec / (1 - 0.85)
      inputs$V0_ini = V0_ini
    }else{
      rewardRate_ini = mean(unlist(optimRewardRates)) * stepSec 
      inputs$rewardRate_ini = rewardRate_ini
    }
   
   # strip the path in outputFile
   outputFile_clean = sub(pattern = sprintf("genData/(exp|sim)*ModelFit(CV)*/[A-Z0-9]*/*%s/", modelName),
                      replacement = "", outputFile)
    
   # fit the model
    withCallingHandlers({
      fit = sampling(object = model, data = inputs, cores = 1, chains = nChain,
                     iter = nIter, control = controlList) 
      print(sprintf("Finish %s !", outputFile_clean))
      write(sprintf("Finish %s !", outputFile_clean), warningFile, append = T, sep = "\n")
    }, warning = function(w){
      warnText = paste(modelName, outputFile_clean, w)
      write(warnText, warningFile, append = T, sep = "\n")
    })
  
  # extract posterior samples
  samples = fit %>%
    rstan::extract(permuted = F, pars = c(paraNames, "totalLL"))
  
  # calculate WAIC and Efficient approximate leave-one-out cross-validation (LOO)
  log_lik = extract_log_lik(fit) 
  WAIC = waic(log_lik)
  LOO = loo(log_lik)
  save("WAIC", "LOO", file = sprintf("%s_waic.RData", outputFile))
  
  # summarise posterior parameters and totalLL
  fitSummary <- summary(fit, pars = c(paraNames, "totalLL"), use_cache = F)$summary
  
  # check ESS and Rhat
  # detect participants with low ESSs and high Rhats 
  ESSCols = which(str_detect(colnames(fitSummary), "Effe")) # columns recording ESSs
  if(any(fitSummary[,ESSCols] < nChain * 100)){
    warnText = paste(modelName, id, "Low ESS")
    write(warnText, warningFile, append = T, sep = "\n")
  }
  RhatCols = which(str_detect(colnames(fitSummary), "Rhat")) # columns recording ESSs
  if(any(fitSummary[,RhatCols] > 1.01)){
    warnText = paste(modelName, id, "High Rhat")
    write(warnText, warningFile, append = T, sep = "\n")
  } 
  
  # check divergent transitions
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  nDt = sum(divergent)
  fitSummary = cbind(fitSummary, nDt = rep(nDt, nrow(fitSummary)))
  
  # write outputs  
  write.table(fitSummary, file = sprintf("%s_summary.txt", outputFile), 
              sep = ",", col.names = F, row.names=FALSE)
}

