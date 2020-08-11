
library("glmnet")
library("beeswarm")

# decision every 1 s, starting at either 0 or 1. 
# exclude immediate quits (defined based on initial position)

# load theoretical forward-looking SV
sv_df = read.csv('output/sv.csv') # has fields t (time) and sv 

# consider decision points every 1 s
# get the value of the SV function at each point
svGrid = data.frame(t=0:19, sv=NA)
svGrid$sv = as.numeric(lapply(svGrid$t, function(tNow) sv_df$sv[sv_df$t==tNow]))

# load all data
source('loadData.R')
allData = loadData()
trialData = allData$trialData         # unpack trial-level data
subjectData = allData$subjectData     # unpack subject-level data
allIDs = names(trialData)             # subject IDs
n = length(allIDs)                    # n
cat('Analyzing data for n','=',n,'subjects.\n')

# initialize
nWith2To7Quits = 0

# model fit for individual subjects
for (sIdx in 1:n) {
  
  # pull this subject's data
  thisID = allIDs[sIdx]
  thisTrialData = trialData[[thisID]]
  subjectRowIdx = (subjectData$ID == thisID)
  
  # exclude immediate quits
  immedQuitIdx = (thisTrialData$initialPos == 'optSmall')
  td_nonImmed = thisTrialData[!immedQuitIdx,] # td for trial data
  
  # for rewarded trials, base the duration on the reward delivery time
  # (this prevents a reward at 3 s from being coded as a choice to wait
  # after passing the 3 s mark)
  td_nonImmed$waitDuration = td_nonImmed$latency
  quitIdx = td_nonImmed$trialResult == "quit"
  td_nonImmed$waitDuration[!quitIdx] <- td_nonImmed$designatedWait[!quitIdx]
  
  # format as a sequence of pseudo choice points every 1 s 
  # loop over all the initial-wait trials
  nTrials = nrow(td_nonImmed)
  pseudoChoices = c() # initialize
  for (trialIdx in 1:nTrials) {
    thisTrial = td_nonImmed[trialIdx,]
    nPoints = ceiling(thisTrial$waitDuration) # number of pseudo-choices
    thisBlock = svGrid[1:nPoints,] # collect the relevant latencies and SVs
    # add pseudo-choice outcomes
    thisBlock$choseWait = rep(1, nPoints) # fill with ones
    if (thisTrial$trialResult == "quit") {
      # if quit, set the last pseudo-choice outcome to zero
      thisBlock$choseWait[nPoints] = 0
    }
    # add a column filled with trial onset time
    thisBlock$trialOnsetTime = rep(1, nPoints) * thisTrial$initialTime
    pseudoChoices = rbind(pseudoChoices, thisBlock)
    
  } # loop over trials
  
  # skip subjects with one or fewer quits
  # this results in skipping 53 of the 273 subjects
  # ** unsure if it's better to treat these as missing or as zero.
  if (sum(pseudoChoices$choseWait==0) <= 1) {
    subjectData[subjectRowIdx,'t_sv_corr'] = NA
    subjectData[subjectRowIdx,'model1_coef_SV'] = NA
    subjectData[subjectRowIdx,'model1Ridge_coef_SV'] = NA
    next
  }
  
  # glmnet produces a warning if there are fewer than 8 quits. 
  # track the number of these cases (not counting those already skipped)
  if (sum(pseudoChoices$choseWait==0) < 8) {
    nWith2To7Quits = nWith2To7Quits + 1
  }
  
  # check the correlation between sv and t (limited to subjects included for model fitting)
  subjectData[subjectRowIdx,'t_sv_corr'] = cor(pseudoChoices$t, pseudoChoices$sv)
  
  # unregularized logistic regression fit
  result = glm(choseWait ~ 1 + sv, family=binomial(link="logit"), data=pseudoChoices)
  subjectData[subjectRowIdx,'model1_coef_SV'] = result$coefficients['sv'] # store the coef for this subject
  
  # the coef is positive for about 67% of subjects
  # however, there are a number of warnings and extreme values.
  # this result as-is is amenable to an overall wilcox.test, but not suitable for testing group differences. 
  # alternative: fit as a mixed-effects model and/or add regularization. 
  # also still to do: add interaction with total time in session; check effect of time within trial. 
  
  # regularized logistic regression fit using glmnet. model 1, SV only
  pseudoChoiceFactor = as.factor(pseudoChoices$choseWait)
  pseudoChoices$zeros = rep.int(0, nrow(pseudoChoices))
  xMat1 = data.matrix(pseudoChoices[c("zeros", "sv")])
  fit1 = glmnet(x=xMat1, y=pseudoChoiceFactor, family="binomial", alpha=0, nlambda=100, intercept=TRUE)
    # ridge regression (alpha=0) for 100 distinct lambda values
  # plot(fit1, xvar="lambda", label=TRUE)
  # coef(fit1, s=0.01)
  subjectData[subjectRowIdx,'model1Ridge_coef_SV'] = coef(fit1, s=0.01)[3] # using regularization level lambda=0.01
  
  # model 2, explanatory variables are SV and t (time within trial)
  xMat2 = data.matrix(pseudoChoices[c("t", "sv")])
  fit2 = glmnet(x=xMat2, y=pseudoChoiceFactor, family="binomial", alpha=0, nlambda=100, intercept=TRUE)
  subjectData[subjectRowIdx,'model2Ridge_coef_time'] = coef(fit2, s=0.01)[2]
  subjectData[subjectRowIdx,'model2Ridge_coef_SV'] = coef(fit2, s=0.01)[3]
  
  # model 3, like model 2 but using log-transformed time. 
  pseudoChoices$logt = log(pseudoChoices$t + 1)
  xMat3 = data.matrix(pseudoChoices[c("logt", "sv")])
  fit3 = glmnet(x=xMat3, y=pseudoChoiceFactor, family="binomial", alpha=0, nlambda=100, intercept=TRUE)
  subjectData[subjectRowIdx,'model3Ridge_coef_logTime'] = coef(fit3, s=0.01)[2]
  subjectData[subjectRowIdx,'model3Ridge_coef_SV'] = coef(fit3, s=0.01)[3]
  
  # model 4: interaction of trialOnsetTime and SV
  pseudoChoices$trialOnsetNorm = pseudoChoices$trialOnsetTime / 300 # normalize to 0-1 (fraction of the session)
  # next set up the interaction term. 
  pseudoChoices$trialOnset_x_SV = pseudoChoices$trialOnsetNorm * pseudoChoices$sv
  xMat4 = data.matrix(pseudoChoices[c("trialOnsetNorm", "sv", "trialOnset_x_SV")])
  fit4 = glmnet(x=xMat4, y=pseudoChoiceFactor, family="binomial", alpha=0, nlambda=100, intercept=TRUE)
  subjectData[subjectRowIdx,'model4Ridge_coef_sessTime'] = coef(fit4, s=0.01)[2]
  subjectData[subjectRowIdx,'model4Ridge_coef_SV'] = coef(fit4, s=0.01)[3]
  subjectData[subjectRowIdx,'model4Ridge_coef_sessTime.x.SV'] = coef(fit4, s=0.01)[4]
  
  # model 5: interaction of trialOnsetTime and log-time-within-trial
  pseudoChoices$trialOnset_x_logt = pseudoChoices$trialOnsetNorm * pseudoChoices$logt
  xMat5 = data.matrix(pseudoChoices[c("trialOnsetNorm", "logt", "trialOnset_x_logt")])
  fit5 = glmnet(x=xMat5, y=pseudoChoiceFactor, family="binomial", alpha=0, nlambda=100, intercept=TRUE)
  subjectData[subjectRowIdx,'model5Ridge_coef_timeInSess'] = coef(fit5, s=0.01)[2]
  subjectData[subjectRowIdx,'model5Ridge_coef_logTimeInTrial'] = coef(fit5, s=0.01)[3]
  subjectData[subjectRowIdx,'model5Ridge_coef_timeInSess.x.logTimeInTrial'] = coef(fit5, s=0.01)[4]
  
  # model 6: both interactions from models 4 and 5
  xMat6 = data.matrix(pseudoChoices[c("trialOnsetNorm", "logt", "sv", "trialOnset_x_logt", "trialOnset_x_SV")])
  fit6 = glmnet(x=xMat6, y=pseudoChoiceFactor, family="binomial", alpha=0, nlambda=100, intercept=TRUE)
  subjectData[subjectRowIdx,'model6Ridge_coef_timeInSess'] = coef(fit6, s=0.01)[2]
  subjectData[subjectRowIdx,'model6Ridge_coef_logTimeInTrial'] = coef(fit6, s=0.01)[3]
  subjectData[subjectRowIdx,'model6Ridge_coef_SV'] = coef(fit6, s=0.01)[4]
  subjectData[subjectRowIdx,'model6Ridge_coef_timeInSess.x.logTimeInTrial'] = coef(fit6, s=0.01)[5]
  subjectData[subjectRowIdx,'model6Ridge_coef_sessTime.x.SV'] = coef(fit6, s=0.01)[6]
  
} # loop over subjects

nNA = sum(is.na(subjectData$model1Ridge_coef_SV))
cat('Number of subjects excluded for <=1 quit:', nNA, '\n')
cat('Number remaining:', n - nNA, '\n')
cat('Number of retained subjects with fewer than 8 quits:', nWith2To7Quits, '\n')

### plots of results

# model 1: SV only. 
plot(ecdf(subjectData$model1_coef_SV), bty="n", xlab="SV coef", ylab="Cumulative proportion", 
     main="Model 1: Unregularized SV coefs"); lines(c(0,0), c(0,1), col='red', lty="dashed")
plot(ecdf(subjectData$model1Ridge_coef_SV), bty="n", xlab="SV coef", ylab="Cumulative proportion", 
     main="Model 1: Ridge-regularized SV coefs"); lines(c(0,0), c(0,1), col='red', lty="dashed")
plot(subjectData$model1_coef_SV, subjectData$model1Ridge_coef_SV, xlab="Unregularized", ylab="Ridge",
     main="Model 1: Unregularized versus ridge-regularized SV coefs")
  lines(c(0,0), c(-10,10), col='red', lty="dashed")
  lines(c(-10,10), c(0,0), col='red', lty="dashed")
  lines(c(-10,10), c(-10,10), col='red', lty="dashed")
beeswarm(model1Ridge_coef_SV ~ group1245, subjectData, cex=1.2, pch=16, 
         xlab="Group", ylab="SV coef", main="Model 1: Ridge-regularized SV coefs",
         bty="n", cex.lab=1.5, cex.axis=1.2)

# model 2: x1 = time within trial, x2 = SV
plot(ecdf(subjectData$t_sv_corr), bty="n", xlab="Pearson r value", ylab="Cumulative proportion", 
     main="Correlation between time-within-trial and SV"); lines(c(0,0), c(0,1), col='red', lty="dashed")
plot(ecdf(subjectData$model2Ridge_coef_time), bty="n", xlab="Time-within-trial coef", ylab="Cumulative proportion", 
     main="Model 2: Ridge-regularized time-within-trial coef"); lines(c(0,0), c(0,1), col='red', lty="dashed")
plot(ecdf(subjectData$model2Ridge_coef_SV), bty="n", xlab="SV coef", ylab="Cumulative proportion", 
     main="Model 2: Ridge-regularized SV coef"); lines(c(0,0), c(0,1), col='red', lty="dashed")
# plot the time coef against the SV coef
r = cor(subjectData$model2Ridge_coef_time, subjectData$model2Ridge_coef_SV, use="complete")
plot(subjectData$model2Ridge_coef_time, subjectData$model2Ridge_coef_SV,
     xlab="Time-within-trial coef", ylab="SV coef", main=sprintf('Model 2 coefs (r = %1.2f)',r))
  lines(c(0,0), c(-10,10), col='red', lty="dashed")
  lines(c(-10,10), c(0,0), col='red', lty="dashed")

# model 3: x1 = log-transformed time within trial, x2 = SV
plot(ecdf(subjectData$model3Ridge_coef_logTime), bty="n", xlab="Log-time-within-trial coef", ylab="Cumulative proportion", 
     main="Model 3: Ridge-regularized log-time-within-trial coef"); lines(c(0,0), c(0,1), col='red', lty="dashed")
plot(ecdf(subjectData$model3Ridge_coef_SV), bty="n", xlab="SV coef", ylab="Cumulative proportion", 
     main="Model 3: Ridge-regularized SV coef"); lines(c(0,0), c(0,1), col='red', lty="dashed")
# plot the log-time coef against the SV coef
r = cor(subjectData$model3Ridge_coef_logTime, subjectData$model3Ridge_coef_SV, use="complete")
plot(subjectData$model3Ridge_coef_logTime, subjectData$model3Ridge_coef_SV,
     xlab="Log-time-within-trial coef", ylab="SV coef", main=sprintf('Model 3 coefs (r = %1.2f)',r))
lines(c(0,0), c(-10,10), col='red', lty="dashed")
lines(c(-10,10), c(0,0), col='red', lty="dashed")

# model 4: x1 = time within session, x2 = SV, x3 = interaction
plot(ecdf(subjectData$model4Ridge_coef_sessTime), bty="n", xlab="Time-within-session coef", ylab="Cumulative proportion", 
     main="Model 4: Ridge-regularized time-within-session coef"); lines(c(0,0), c(0,1), col='red', lty="dashed")
plot(ecdf(subjectData$model4Ridge_coef_SV), bty="n", xlab="SV coef", ylab="Cumulative proportion", 
     main="Model 4: Ridge-regularized SV coef"); lines(c(0,0), c(0,1), col='red', lty="dashed")
plot(ecdf(subjectData$model4Ridge_coef_sessTime.x.SV), bty="n", xlab="Interaction coef", ylab="Cumulative proportion", 
     main="Model 4: Interaction of time-within-session X SV"); lines(c(0,0), c(0,1), col='red', lty="dashed")

# model 5: x1 = time within session, x2 = log-transformed time within trial, x3 = interaction
plot(ecdf(subjectData$model5Ridge_coef_timeInSess), bty="n", xlab="Time-within-session coef", ylab="Cumulative proportion", 
     main="Model 5: Ridge-regularized time-within-session coef"); lines(c(0,0), c(0,1), col='red', lty="dashed")
plot(ecdf(subjectData$model5Ridge_coef_logTimeInTrial), bty="n", xlab="Log-time-within-trial coef", ylab="Cumulative proportion", 
     main="Model 5: Ridge-regularized log-time-within-trial coef"); lines(c(0,0), c(0,1), col='red', lty="dashed")
plot(ecdf(subjectData$model5Ridge_coef_timeInSess.x.logTimeInTrial), bty="n", xlab="Interaction coef", ylab="Cumulative proportion", 
     main="Model 5: Interaction of time-within-session X log-time-within-trial"); lines(c(0,0), c(0,1), col='red', lty="dashed")

# model 6: x1 = time within session, x2 = log-transformed time within trial, x3 = SV, 
#   x4 = x1-by-x2 interaction, x5 = x1-by-x3 interaction
plot(ecdf(subjectData$model6Ridge_coef_timeInSess), bty="n", xlab="Time-within-session coef", ylab="Cumulative proportion", 
     main="Model 6: Ridge-regularized time-within-session coef"); lines(c(0,0), c(0,1), col='red', lty="dashed")
plot(ecdf(subjectData$model6Ridge_coef_logTimeInTrial), bty="n", xlab="Log-time-within-trial coef", ylab="Cumulative proportion", 
     main="Model 6: Ridge-regularized log-time-within-trial coef"); lines(c(0,0), c(0,1), col='red', lty="dashed")
plot(ecdf(subjectData$model6Ridge_coef_SV), bty="n", xlab="SV coef", ylab="Cumulative proportion", 
     main="Model 6: Ridge-regularized SV coef"); lines(c(0,0), c(0,1), col='red', lty="dashed")
plot(ecdf(subjectData$model6Ridge_coef_timeInSess.x.logTimeInTrial), bty="n", xlab="Interaction coef", ylab="Cumulative proportion", 
     main="Model 6: Interaction of time-within-session X log-time-within-trial"); lines(c(0,0), c(0,1), col='red', lty="dashed")
plot(ecdf(subjectData$model6Ridge_coef_sessTime.x.SV), bty="n", xlab="Interaction coef", ylab="Cumulative proportion", 
     main="Model 6: Interaction of time-within-session X SV"); lines(c(0,0), c(0,1), col='red', lty="dashed")





