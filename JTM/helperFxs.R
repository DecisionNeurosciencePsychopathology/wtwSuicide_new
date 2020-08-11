
# analysis helper functions



# check the distribution of scheduled delays
# ...as measured in number of key presses (for the instrumental version of the task)
scheduledDelays <- function(blockData,blockLabel) {
  cat(sprintf('Scheduled delays for %s\n',blockLabel))
  bkDelays = blockData$designatedWait
  print(summary(bkDelays))
  # empirical cumulative distribution of scheduled delays
  fn <- ecdf(blockData$designatedWait)
  plot(fn, main = sprintf('Scheduled delays: %s',blockLabel), xlab='Scheduled delay (s)',
       ylab='Cumulative proportion')
  # autocorrelation function
  acfOutput <- acf(bkDelays, lag.max=20, main = sprintf('Scheduled delays: %s',blockLabel))
}


# plot trialwise responses in detail
trialPlots <- function(blockData,blockLabel) {
  # vectors to be plotted
  rwdIdx = blockData$trialResult == "win"
  quitIdx = blockData$trialResult == "quit"
  rwdTrialNo = blockData$trialNum[rwdIdx]
  quitTrialNo = blockData$trialNum[quitIdx]
  rwdSchedDelay = blockData$designatedWait[rwdIdx]
  quitSchedDelay = blockData$designatedWait[quitIdx]
  waitDuration = blockData$latency
  quitTime = waitDuration[quitIdx]
  # other parameters
  nTrials = nrow(blockData)
  maxDelay = max(blockData$designatedWait)
  # make the plot and add series
  plot(1, type='n', xlim=c(1,nTrials), ylim=c(0,maxDelay), bty='n',
       xlab='Trial', ylab='Wait duration (s)', main=sprintf('Trial data: %s',blockLabel))
  lines(rwdTrialNo, rwdSchedDelay, col='blue', type='o', lwd=2, pch=16)
  lines(quitTrialNo, quitTime, col='red', type='o', lwd=2, pch=16)
  lines(quitTrialNo, quitSchedDelay, col='black', type='o', lwd=2, lty=0, pch=16)
}


# calculate the number of cursor moves between the 'quit' and 'wait' boxes
nCursorMoves <- function(blockData) {
  #   +1 if a trial is not a fast quit but was quit eventually. (wait -> quit during the trial)
  #   +1 if a trial is not a fast quit and the previous trial was quit. (quit -> wait b/w trials)
  #   +1 if a trial *is* a fast quit and the previous trial was rewarded. (wait -> quit b/w trials)
  immedQuitIdx = (blockData$initialPos == 'optSmall')
  rwdIdx = blockData$trialResult == "win"
  quitIdx = blockData$trialResult == "quit"
  priorRwdIdx = c(FALSE, head(rwdIdx, n=-1)) # TRUE if the previous trial was rewarded
  priorQuitIdx = c(FALSE, head(quitIdx, n=-1)) # FALSE if the previous trial was quit
  nMoves = sum(!immedQuitIdx & quitIdx) + sum(!immedQuitIdx & priorQuitIdx) + sum(immedQuitIdx & priorRwdIdx)
  return(nMoves)
}



# calculate kaplan-meier and area under the curve
kmsc <- function(blockData, tMax, blockLabel, makePlot=FALSE, grid=0) {
  
  if (nrow(blockData) == 0) {
    return(list(kmT=NA, kmF=NA, auc=NA))
  }
  
  waitDuration = blockData$latency
  quitIdx = blockData$trialResult == "quit"
  # for rewarded trials, base the duration on the reward delivery time
  waitDuration[!quitIdx] <- blockData$designatedWait[!quitIdx]
  # fit the survival function
  kmfit <- survfit(Surv(waitDuration, quitIdx, type='right') ~ 1, 
                 type='kaplan-meier', conf.type='none', start.time=0, se.fit=FALSE)
  # extract elements of the survival curve object (?survfit.object)
  kmT = kmfit$time
  kmF = kmfit$surv
  # add a point at zero
  kmT = c(0, kmT)
  kmF = c(1, kmF)
  # keep only points up through tMax
  keepIdx = kmT<=tMax
  kmT <- kmT[keepIdx]
  kmF <- kmF[keepIdx]
  # extend the last value to exactly tMax
  kmT <- c(kmT, tMax)
  kmF <- c(kmF, tail(kmF,1))
  # calculate auc
  auc <- sum(diff(kmT) * head(kmF,-1))
  # plot if requested
  if (makePlot) {
    plot(kmT, kmF, type='s', frame.plot=FALSE, xlab='Delay (s)', ylab='Survival rate',
         main=sprintf('KMSC: subject %s (AUC = %1.1f)',blockLabel,auc), ylim=c(0,1), xlim=c(0,tMax),
         xaxp=c(0,180,9))
  }
  # put the survival curve on a standard grid
  kmOnGrid = vector()
  for (gIdx in 1:length(grid)) {
    g = grid[gIdx]
    # use the last point where t is less than or equal to the current grid value
    kmOnGrid[gIdx] = kmF[max(which(kmT<=g))]
  }
  return(list(kmT=kmT, kmF=kmF, auc=auc, kmOnGrid=kmOnGrid))
}



# willingness to wait time-series
wtwTS <- function(blockData, tGrid) {
  trialWTW = numeric(length = nrow(blockData)) # initialize the per-trial estimate of WTW
  quitIdx = blockData$trialResult == "quit"
  # use either the the designated wait (for reward trials) or actual time waited (for quit trials)
  timeWaited = blockData$designatedWait
  timeWaited[quitIdx] = blockData$latency[quitIdx]
  # find the longest time waited up through the first quit trial
  #   (or, if there were no quit trials, the longest time waited at all)
  #   that will be the WTW estimate for all trials prior to the first quit
  firstQuit = which(quitIdx)[1]
  if (is.na(firstQuit)) {firstQuit = nrow(blockData)} # if no quit, set to the last trial
  currentWTW = max(timeWaited[1:firstQuit])
  thisTrialIdx = firstQuit - 1
  trialWTW[1:thisTrialIdx] = currentWTW
  # iterate through the remaining trials, updating currentWTW
  while (thisTrialIdx < nrow(blockData)) {
    thisTrialIdx = thisTrialIdx + 1
    if (quitIdx[thisTrialIdx]) {currentWTW = timeWaited[thisTrialIdx]}
    else {currentWTW = max(currentWTW, timeWaited[thisTrialIdx])}
    trialWTW[thisTrialIdx] = currentWTW
  }
  # convert from per-trial to per-second over the course of the block
  timeWTW = numeric(length = length(tGrid)) # initialize output
  binStartIdx = 1
  thisTrialIdx = 0
  while (thisTrialIdx < nrow(blockData)) {
    thisTrialIdx = thisTrialIdx + 1
    binEndTime = blockData$outcomeTime[thisTrialIdx]
    binEndIdx = max(which(tGrid < binEndTime)) # last grid point that falls within this trial
    timeWTW[binStartIdx:binEndIdx] = trialWTW[thisTrialIdx]
    binStartIdx = binEndIdx + 1
  }
  # extend the final value to the end of the vector
  timeWTW[binStartIdx:length(timeWTW)] = trialWTW[thisTrialIdx]
  return(timeWTW)
}






