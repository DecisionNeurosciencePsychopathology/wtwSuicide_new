plotValueAndAction = function(tWaits, decValues_, thisTrialData, tIdx){
  
  nTrial = nrow(thisTrialData)
  if(tIdx == nTrial){
    
  }
  
  actionData = data.frame(
    trialEarnings = thisTrialData$trialEarnings[c(tIdx, tIdx + 1)],
    timeWaited = thisTrialData$timeWaited[c(tIdx, tIdx + 1)],
    isQuit = (thisTrialData$trialEarnings[c(tIdx, tIdx + 1)] == smallReward),
    trial = paste0("Trial ", c(tIdx, tIdx + 1)),
    text = ifelse(thisTrialData$trialEarnings[c(tIdx, tIdx + 1)] == smallReward, "Q", "R")
  )
  valueData = data.frame(t = rep(tWaits, 2) - iti,
             trial = rep(paste0("Trial ", c(tIdx, tIdx + 1)), each = length(tWaits)),
             decValue = as.vector(decValues_[,c(tIdx, tIdx + 1)])) 
  
  upperLimit = max(valueData$decValue) * 1.1
  
  p = valueData %>% 
    ggplot(aes(t, decValue)) + geom_point() + facet_grid(~trial) +
    myTheme +
    xlab("Elapsed time (s)") + ylab("Decision value in favor of waiting") 
  
  p + geom_point(aes(x = timeWaited, y = upperLimit, color = isQuit), size = 3, data = actionData, inherit.aes = F) +
    scale_color_manual(values=c("red", "red"))
    
    
    
}
