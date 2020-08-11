
# anGrp.R
# Dombrovski WTW data

library("survival")
library("beeswarm")
library("matrixStats")

# analysis sub-functions

# still need to adapt paths
if (Sys.info()[6]!='Alex')
{source('helperFxs.R')
source('loadData.R')
  }
if (Sys.info()[6]=='Alex')
{setwd('~/code/wtwSuicide/JTM')
  source('~/code/wtwSuicide/JTM/helperFxs.R')
  source('~/code/wtwSuicide/JTM/loadData.R')}


# load all data
allData = loadData()
trialData = allData$trialData         # unpack trial-level data
subjectData = allData$subjectData     # unpack subject-level data
allIDs = names(trialData)             # subject IDs
n = length(allIDs)                    # n
cat('Analyzing data for n','=',n,'subjects.\n')

# control which individual-level plots to generate
plotScheduledDelays = FALSE
plotTrialwiseData = FALSE
plotKMSC = FALSE
plotKMSC_nonImmed = TRUE

# initialize matrices to hold subject-level time series data
kmGrid = seq(0, 20, by=0.1) # grid on which to average survival curves.
subjectKMSC = matrix(NA, nrow=n, ncol=length(kmGrid)) # structure to hold all subjects' survival curves
subjectKMSC_nonImmed = matrix(NA, nrow=n, ncol=length(kmGrid)) # version for the non-immediate-quit analysis
tsGrid = seq(0, 5*60, by=1) # grid on which to average whole-block WTW time series
subjectWTWTS = matrix(NA, nrow=n, ncol=length(tsGrid)) # structure to hold all subjects' WTW time series

# initialize new subject-level columns to hold group data
# subjectData$grpAUC = NA
# subjectData$grpAUC_nonImmed = NA
# subjectData$grpEarnings = NA
# subjectData$nImmedQuits = NA

# descriptive statistics for individual subjects
for (sIdx in 1:n) {

  # pull this subject's data
  thisID = allIDs[sIdx]
  thisTrialData = trialData[[thisID]]
  subjectRowIdx = (subjectData$ID == thisID)
  subjectData[subjectRowIdx,'grpEarnings'] = max(thisTrialData$totalEarned)
  
  # plot and summarize the distribution of scheduled delays
  if (plotScheduledDelays) {
    scheduledDelays(thisTrialData,thisID)
  }

  # plot trial-by-trial data
  if (plotTrialwiseData) {
    trialPlots(thisTrialData,thisID)
  }

  # survival analysis
  tMax = 20 # time window for the survival analysis (the longest wait time was 20 s)
  kmscResults = kmsc(thisTrialData, tMax, thisID, plotKMSC, kmGrid)
  subjectData[subjectRowIdx,'grpAUC'] = kmscResults[['auc']]
  subjectKMSC[subjectRowIdx,] = kmscResults[['kmOnGrid']]
  
  # survival analysis excluding immediate quits. 
  immedQuitIdx = (thisTrialData$initialPos == 'optSmall')
  subjectData[subjectRowIdx,'nImmedQuits'] = sum(immedQuitIdx) # record the number of immediate quits
  kmscResults_nonImmed = kmsc(thisTrialData[!immedQuitIdx,], tMax, thisID, plotKMSC_nonImmed, kmGrid)
  subjectData[subjectRowIdx,'grpAUC_nonImmed'] = kmscResults_nonImmed[['auc']]
  subjectKMSC_nonImmed[subjectRowIdx,] = kmscResults_nonImmed[['kmOnGrid']]

  # note: it's not feasible to compute the KMSC separately for post-long-delay and post-short-delay trials
  # because for many subjects, the set of post-long-delay trials includes no 20 s scheduled delays
  # (long delays have only 25% frequency and scheduled delays are anticorrelated)
  
  # rate of immediate quits, post-long and post-short delay
  isLongDelay = (thisTrialData$designatedWait == 20)
  postLongDelay = c(FALSE, head(isLongDelay, n=-1))
  postShortDelay = c(FALSE, head(!isLongDelay, n=-1))
  subjectData[subjectRowIdx,'propImmedQuits_postLongDelay'] = mean(immedQuitIdx[postLongDelay])
  subjectData[subjectRowIdx,'propImmedQuits_postShortDelay'] = mean(immedQuitIdx[postShortDelay])
  
  # calculate an index of flexibility
  # i.e. how many times the cursor moved between the two boxes, either within or between trials
  #   +1 if this trial is not a fast quit but was quit eventually. (wait -> quit during the trial)
  #   +1 if this trial is not a fast quit and the previous trial was quit. (quit -> wait b/w trials)
  #   +1 if this trial *is* a fast quit and the previous trial was rewarded. (wait -> quit b/w trials)
  subjectData[subjectRowIdx,'nCursorMoves'] = nCursorMoves(thisTrialData)
  # for testing: can uncomment the line below and turn on plotTrialwiseData
  # cat('Number of cursor moves: ',subjectData[subjectRowIdx,'nCursorMoves'],"\n")
  
  # WTW time series
  subjectWTWTS[subjectRowIdx,] = wtwTS(thisTrialData, tsGrid)
  # ***incude a measure of consistency for the WTW time series
  
  # wait for input before continuing. 
  if (any(plotScheduledDelays, plotTrialwiseData, plotKMSC)) {
    readline(prompt = paste('subject',thisID,'(hit ENTER to continue)'))
  }

  # temporary: inspect only a few subjects
  # if (sIdx>2) {break}
  
}

##################################################
##### summarize and save group-level results #####

cat('Group Ns by the group1245 field:')
print(table(subjectData$group1245))

# save group summary statistics to a csv file
outfname = sprintf('output/grpSummary_n=%d.csv', n)
write.csv(x=subjectData, file=outfname, row.names=FALSE)
cat('Saved group summary output to:',outfname,'\n')

# plot and summarize AUC results
# cat('Distribution of AUC values:\n')
# print(summary(subjectData$grpAUC)) # print a summary
fn <- ecdf(subjectData$grpAUC) 
plot(fn, main = sprintf('AUC values, n = %d',n), xlab='AUC (s)',
     ylab='Cumulative proportion') # plot the empirical CDF
hist(subjectData$grpAUC, breaks=16, freq=TRUE, main = sprintf('AUC values (n = %d)',n), 
     xlab='AUC (s)', xaxp=c(0,180,9)) # plot in histogram form

### beeswarm plots and Kruskal-Wallis tests by clinical group

# set up the grouping variable
# subjectData$PATTYPE = as.factor(subjectData$PATTYPE) 

# total earnings by group
beeswarm(grpEarnings ~ group1245, subjectData, cex=1.2, pch=16, 
         xlab="Group", ylab="Earnings", main="Task earnings",
         bty="n", cex.lab=1.5, cex.axis=1.2)
print(kruskal.test(grpEarnings ~ group1245, subjectData))

# AUC by group
beeswarm(grpAUC ~ group1245, subjectData, cex=1.2, pch=16, 
         xlab="Group", ylab="AUC (s)", main="Willingness to wait (AUC)",
         ylim=c(0,20), bty="n", cex.lab=1.5, cex.axis=1.2)
print(kruskal.test(grpAUC ~ group1245, subjectData))

# nImmedQuits by group
beeswarm(nImmedQuits ~ group1245, subjectData, cex=1.2, pch=16, 
         xlab="Group", ylab="Number of immediate quits", main="Immediate quits",
         bty="n", cex.lab=1.5, cex.axis=1.2)
print(kruskal.test(nImmedQuits ~ group1245, subjectData))

# AUC_nonImmed by group
beeswarm(grpAUC_nonImmed ~ group1245, subjectData, cex=1.2, pch=16, 
         xlab="Group", ylab="AUC (s)", main="Willingness to wait (AUC) w/o immediate quits",
         ylim=c(0,20), bty="n", cex.lab=1.5, cex.axis=1.2)
print(kruskal.test(grpAUC_nonImmed ~ group1245, subjectData))

# difference in propImmedQuits after long versus short delays on the previous trial
subjectData$prevDelayDiff = subjectData$propImmedQuits_postLongDelay - subjectData$propImmedQuits_postShortDelay
beeswarm(prevDelayDiff ~ group1245, subjectData, cex=1.2, pch=16, 
         xlab="Group", ylab="Difference of rates (post-long minus post-short)", 
         main="Rate of fast quits after long vs. short prior delays",
         bty="n", cex.lab=1.5, cex.axis=1.2)
print(kruskal.test(prevDelayDiff ~ group1245, subjectData))
print(wilcox.test(subjectData$propImmedQuits_postLongDelay - subjectData$propImmedQuits_postShortDelay))

# number of cursor movements
beeswarm(nCursorMoves ~ group1245, subjectData, cex=1.2, pch=16, 
         xlab="Group", ylab="Cursor movements", main="Number of cursor movements",
         bty="n", cex.lab=1.5, cex.axis=1.2)
print(kruskal.test(nCursorMoves ~ group1245, subjectData))

# plot group-average KMSC
plot(1, type="l", main="Subgroup mean KMSC", 
     xlab="Time (s)", ylab="Survival rate", bty="n", xlim=c(0,20), ylim=c(0,1))
lines(kmGrid, colMeans(subjectKMSC_nonImmed[subjectData$group1245==1,]), col='black', type='l', lwd=3, pch=16)
lines(kmGrid, colMeans(subjectKMSC_nonImmed[subjectData$group1245==2,]), col='red', type='l', lwd=3, pch=16)
lines(kmGrid, colMeans(subjectKMSC_nonImmed[subjectData$group1245==4,]), col='green', type='l', lwd=3, pch=16)
lines(kmGrid, colMeans(subjectKMSC_nonImmed[subjectData$group1245==5,]), col='blue', type='l', lwd=3, pch=16)
legend("bottomleft", c("Grp 1", "Grp 2", "Grp 4", "Grp 5"), bty="n", lty=1, lwd=3, 
       col=c("black", "red", "green", "blue"))

# plot group-average WTW time series
# can test the slope of this (and whether it differs by group)
# can tabulate the absolute difference from optimal over time.
# use matplot instead? put each variable in 1 col. 
plot(1, type="l", main="Subgroup WTW time series", 
     xlab="Time in block (s)", ylab="WTW (s)", bty="n", xlim=c(0,300), ylim=c(0,20))
lines(tsGrid, colMeans(subjectWTWTS[subjectData$group1245==1,]), col='black', type='l', lwd=3, pch=16)
lines(tsGrid, colMeans(subjectWTWTS[subjectData$group1245==2,]), col='red', type='l', lwd=3, pch=16)
lines(tsGrid, colMeans(subjectWTWTS[subjectData$group1245==4,]), col='green', type='l', lwd=3, pch=16)
lines(tsGrid, colMeans(subjectWTWTS[subjectData$group1245==5,]), col='blue', type='l', lwd=3, pch=16)
legend("bottomleft", c("Grp 1", "Grp 2", "Grp 4", "Grp 5"), bty="n", lty=1, lwd=3, 
       col=c("black", "red", "green", "blue"))





