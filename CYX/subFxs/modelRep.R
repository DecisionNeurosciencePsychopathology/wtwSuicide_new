modelRep = function(trialData, ids, nRep, modelName){
  # 
  nSub = length(ids)
  
  # get the generative model 
  source(sprintf("subFxs/gnrModels/%s.R", modelName))
  gnrModel = get(modelName)
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  
  # normative analysis 
  iti = 2
  normResults = expSchematics(iti, F)
  
  # initialize outputs
  repTrialData = vector(length = nSub * nRep, mode ='list')
  repNo = matrix(1 : (nSub * nRep), nrow = nRep, ncol = nSub)
  
  # loop over participants
  for(sIdx in 1 : nSub){
    # prepare empirical data 
    ID = ids[sIdx]
    thisTrialData = trialData[[ID]]
    # load individually fitted paramters 
    fitSummary = read.table(sprintf("genData/expModelFit/%s/s%s_summary.txt",  modelName, ID),sep = ",", row.names = NULL)
    paras =  fitSummary[1 : nPara,1]
    # simulate nRep times
    for(rIdx in 1 : nRep){
      immQuit = thisTrialData$initialPos == "optSmall"
      thisTrialData = thisTrialData[!immQuit, ] 
      tempt = gnrModel(paras, thisTrialData$condition, thisTrialData$scheduledWait, normResults)
      repTrialData[[repNo[rIdx, sIdx]]] = tempt
    }
  }
  
  # initialize 
  muWTWRep_ = matrix(NA, nrow = nRep , ncol = nSub)
  stdWTWRep_ = matrix(NA, nrow = nRep, ncol = nSub)
  timeWTW_ =  matrix(NA, nrow = length(tGrid), ncol = nSub)
  for(sIdx in 1 : nSub){
    ID = ids[sIdx]
    condition = "LP"
    timeWTW = matrix(NA, nrow = length(tGrid), ncol = nRep)
    for(rIdx in 1 : nRep){
      thisRepTrialData = repTrialData[[repNo[rIdx, sIdx]]]
      thisRepTrialData = do.call(cbind.data.frame, thisRepTrialData[c(1, 3: 6)])
      kmscResults = kmsc(thisRepTrialData, delayMax, F, kmGrid)
      muWTWRep_[rIdx, sIdx] = kmscResults$auc
      stdWTWRep_[rIdx, sIdx] = kmscResults$stdWTW
      wtwResults = wtwTS(thisRepTrialData, tGrid, delayMax, F)
      timeWTW[,rIdx] = wtwResults$timeWTW
    }
    timeWTW_[,sIdx] = apply(timeWTW, 1, mean, na.rm = T)
  }
  
  
  ## summarise WTW across simulations for replicated data 
  muWTWRep_mu = apply(muWTWRep_, MARGIN = 2, mean) # mean of average willingness to wait
  muWTWRep_std = apply(muWTWRep_, MARGIN = 2, sd) # std of average willingess to wait
  stdWTWRep_mu = apply(stdWTWRep_, MARGIN = 2, mean) # mean of std willingness to wait
  stdWTWRep_std = apply(stdWTWRep_, MARGIN = 2, sd) # std of std willingess to wait
  
  outputs = list(
    muWTWRep_mu = muWTWRep_mu,
    muWTWRep_std = muWTWRep_std,
    stdWTWRep_mu = stdWTWRep_mu,
    stdWTWRep_std = stdWTWRep_std,
    "timeWTW_" = timeWTW_,
    repTrialData = repTrialData,
    'repNo' = repNo
  )
  return(outputs)
}
