# generate value functions with individual fitted parameters
genGroupValueFun = function(modelName){
  set.seed(123)
  # create output directories
  dir.create("genData/valueFun")
  
  # load experiment parameters
  load('expParas.RData')
  iti = 2
  stepSec  = 1
  
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
  # medianParas = expPara %>% filter(passCheck) %>% select(c(paraNames, "group")) %>%
  #   gather("para", "value", -c("group")) %>%
  #   group_by(group, para) %>%
  #   summarise(median = median(value))
  

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
  
  # groups
  groups = c("Controls", "Depressed", "Ideators", "Attempters")
  nGroup = length(groups)
  
  # initialize the output 
  decValues_group.mean = list()
  # loop over groups 
  for(i in 1 : nGroup){
    ## emp data 
    group = groups[i]
    trialDataSelect = trialData[hdrData$group == group & passCheck]
    idsSelect = ids[hdrData$group == group & passCheck]
    nSub = length(trialDataSelect)
    
    ## initialize
    decValues_resampled_ = array(NA, dim = c(delayMax / stepSec + 1, 301, nSub))
  
    for(sIdx in 1 : nSub){
      ID = idsSelect[sIdx]
      thisTrialData = trialDataSelect[[ID]]
      immQuit = thisTrialData$initialPos == "optSmall"
      thisTrialData = thisTrialData[!immQuit, ] 
      
      
      ## when a trial ends 
      Ts = thisTrialData$timeWaited + iti
      Ts[thisTrialData$trialEarnings == smallReward] = floor(Ts[thisTrialData$trialEarnings == smallReward] / stepSec) * stepSec 
      ## trial-wise earnings
      Rs = thisTrialData$trialEarnings
      ## number of made actions in each trial 
      nMadeActions = floor(thisTrialData$timeWaited / stepSec) + 1
      
      # load individually fitted paramters 
      fitSummary = read.table(sprintf("genData/expModelFit/%s/s%s_summary.txt",  modelName, ID),sep = ",", row.names = NULL)
      paras =  fitSummary[1 : nPara,1]

      ## compute value functions
      tempt = vrfModel(paras, Rs, Ts, nMadeActions, normResults)
      Qwaits_ = tempt$Qwaits_
      V0_ = tempt$V0_
      decValues_ = Qwaits_ - matrix(rep(V0_ + smallReward, each = nrow(Qwaits_)), ncol = length(Ts))
      
      ## plot decision values in favor of waiting and waiting durations 
      # plotValueAndAction(tWaits, decValues_, thisTrialData, 3) 
      
      ## interpolate 
      x = c(0, head(thisTrialData$sellTime, -1)) # when the value function is computed
      xout = 0 : 300
      tempt = lapply(1 : nrow(Qwaits_), function(i) {
        yout = approx(x, decValues_[i,], xout = xout)$y
        yout[is.na(yout) & xout > 150] = tail(yout[!is.na(yout)], 1)
        return(yout)
      })
      decValues_resampled_[ , , sIdx] = matrix(unlist(tempt), byrow =T, ncol = 301)
    }
    
    # calculate the mean
    decValues_group.mean[[group]] = apply(decValues_resampled_, MARGIN = c(1, 2), function(x) mean(x, na.rm = T))
  }
  
  # save data
  tWaits = seq(iti, iti + delayMax, by = stepSec) - iti
  dir.create("genData/valueFun")
  for(i in 1 : nGroup){
    group = groups[i]
    output = decValues_group.mean[[group]]
    colnames(output) = as.character(seq(0, 300, by = stepSec))
    fileName = sprintf("genData/valueFun/%s.csv", group)
    write.csv(output, file = fileName, row.names = tWaits)
  }
  

  plotTimes = seq(0, 300, by = 60)
  data.frame(
    decValue = unlist(lapply(1 : nGroup, function(i) as.vector(decValues_group.mean[[i]][,which(0 : 300 %in% plotTimes)]))),
    time = rep(tWaits - iti, length(plotTimes) * nGroup),
    taskTime = rep(plotTimes, each = length(tWaits), nGroup),
    group = rep(groups, each = length(tWaits) * length(plotTimes))
  ) %>% 
    mutate(group = factor(group, levels = c("Controls", "Depressed", "Ideators", "Attempters")),
           taskTime = factor(paste0(taskTime, "s"), levels = paste0(unique(taskTime), "s")))%>%
    ggplot(aes(time, decValue, color = group)) + geom_line() + facet_grid(~taskTime) +
    myTheme + ylab("Relative value of waiting") + xlab("Elapsed time (s)")
  
}

