# generate value functions with individual fitted parameters
medianParaValueFun = function(modelName){
  set.seed(123)
  # create output directories
  dir.create("genData/medianParaValueFun")
  
  # load experiment parameters
  load('expParas.RData')
  iti = 2
  stepSec  = 1
  tWaits = seq(iti, iti + delayMax, by = stepSec) - iti
  
  # load packages and sub functions 
  library("tidyverse")
  source("expSchematics.R")
  source("subFxs/plotThemes.R")
  source("subFxs/helpFxs.R") 
  source("subFxs/loadFxs.R") # 
  source("subFxs/analysisFxs.R") 
  source("subFxs/modelRep.R")
  source("subFxs/plotValueAndAction.R")
  
  
  # check fit
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  expPara = loadExpPara(paraNames, sprintf("genData/expModelFit/%s", modelName))
  passCheck = checkFit(paraNames, expPara)
  expPara$group = hdrData$group
  
  # calculate the group-median parameters
  tempt = expPara %>% filter(passCheck) %>% select(c(paraNames, "group")) %>%
    gather("para", "value", -c("group")) %>%
    mutate(para = factor(para, levels = paraNames)) %>%
    group_by(para) %>%
    summarise(median = median(value)) 
  medianParas = tempt$median

  # load empirical data 
  allData = loadAllData()
  hdrData = allData$hdrData 
  trialData = allData$trialData       
  ids = hdrData$ID
  nSub = length(ids)
  
  # model 
  source(sprintf("subFxs/vrfModels/%s.R", modelName))
  vrfModel = get(modelName)
  
  # normative analysis
  normResults = expSchematics(iti, F)
  
  # 
  for(sIdx in 1 : nSub){
    ID = ids[sIdx]
    thisTrialData = trialData[[ID]]
    immQuit = thisTrialData$initialPos == "optSmall"
    thisTrialData = thisTrialData[!immQuit, ] 

    ## when a trial ends 
    Ts = thisTrialData$timeWaited + iti
    Ts[thisTrialData$trialEarnings == smallReward] = floor(Ts[thisTrialData$trialEarnings == smallReward] / stepSec) * stepSec 
    ## trial-wise earnings
    Rs = thisTrialData$trialEarnings
    ## number of made actions in each trial 
    nMadeActions = floor(thisTrialData$timeWaited / stepSec) + 1
    
    
    ## compute value functions
    tempt = vrfModel(medianParas, Rs, Ts, nMadeActions, normResults)
    Qwaits_ = tempt$Qwaits_
    V0_ = tempt$V0_
    decValues_ = Qwaits_ - matrix(rep(V0_ + smallReward, each = nrow(Qwaits_)), ncol = length(Ts))
    
    # check the value functions
    # plotValueAndAction(tWaits, decValues_, thisTrialData, 2) 
    
    # sabe the value functions, each col is a value function of a trial
    fileName = sprintf("genData/medianParaValueFun/%s.csv", ID)
    write.table(decValues_, file = fileName, row.names = F, col.names = F, sep = ",")
  }

  
}

