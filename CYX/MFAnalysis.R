# model-free analysis

# outputs (summarised stats for each participant, nSub):
# sumStats = {
  # ID : [nSubx1 ID]
  # condition : [nSubx1 fac]
  # nExcl : [nSubx1 int] # total number of excluded trials 
  # muWTWs : [nSubx1 num] # average willingness to wait (WTW), measured by area under the Kaplan-Meier survival curve
  # stdWTWs : [nSubx1 num] # standard deviation of WTW, measured in Kaplan-Meier survival analysis
  # totalEarnings_s :  [nSubx1 num] 
# }
# timeWTW_ : list(84x1) # wtw timecourse, each element is a vector
# trialWTW_ : list(84x1) # trial-wise WTW, each element is a vector
# survCurve_ : list(84x1) # Kaplan-Meier survival curve, each element is a vector

MFAnalysis = function(){
  # load libraries
  source('subFxs/loadFxs.R') 
  source("subFxs/plotThemes.R")
  source('subFxs/analysisFxs.R') 
  library('tidyverse')
  
  # create the output directory
  dir.create("genData")
  dir.create("genData/MFAnalysis")
  
  # load experiment parameters
  load("expParas.RData")
  
  # load exp data
  allData = loadAllData() # partipants who never moved the cursor has been excluded
  trialData = allData$trialData   
  ids = names(trialData)
  nSub = length(ids)                    # n
  cat('Analyzing data for',nSub,'subjects.\n')
  
  # initialize output variables 
  nExcls = numeric(length = nSub)
  muWTWs = numeric(length = nSub) 
  stdWTWs = numeric(length = nSub) 
  totalEarnings_s =  numeric(length = nSub) 
  timeWTW_ = vector(mode = "list", length = nSub) 
  trialWTW_ = vector(mode = "list", length = nSub) 
  survCurve_ = vector(mode = "list", length = nSub) 
  
  # loop over individuals
  for (sIdx in 1 : nSub) {
    # load trialData 
    ID = ids[sIdx]
    thisTrialData = trialData[[ID]]

    # calcualte totalEarnings
    totalEarnings_s[sIdx] =  sum(thisTrialData$trialEarnings)
    
    # plot trial-wise behavior
    # trialPlots(thisTrialData)
    
    # exclude trials where the agent quit immediately to get rid of autocorrrelation
    # immQuit = thisTrialData$timeWaited < 0.1 & thisTrialData$initialPos == "optSmall"
    immQuit = thisTrialData$initialPos == "optSmall"
    
    # excluded = immQuit & tmp
    excluded = immQuit
    nExcls[sIdx] = sum(excluded, rm.na = T)
    
    # survival analysis
    kmscResults = kmsc(thisTrialData[!excluded,], delayMax, F, kmGrid)

    # 
    muWTWs[sIdx] = kmscResults[['auc']]
    survCurve_[[sIdx]] = kmscResults$kmOnGrid
    stdWTWs[[sIdx]] = kmscResults$stdWTW
    
    # WTW timecourse
    wtwtsResults = wtwTS(thisTrialData[!excluded,], tGrid, min(delayMax), F)
    timeWTW_[[sIdx]] = wtwtsResults$timeWTW
    trialWTW_[[sIdx]] = wtwtsResults$trialWTW
  }
  
  # return outputs
  sumStats = data.frame(
    ID = ids,
    condition = rep("LP", nSub),
    nExcl = nExcls,
    totalEarnings = totalEarnings_s,
    muWTW = muWTWs,
    stdWTW = stdWTWs
  )
  outputs = list(
    sumStats = sumStats,
    survCurve_ = survCurve_,
    trialWTW_ = trialWTW_,
    timeWTW_ = timeWTW_ 
  )
  return(outputs)
}
