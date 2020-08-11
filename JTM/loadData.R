
# function to load the Dombrovski data set.
#   trial data for all subjects are in wtw_df1.csv
#   subject-level data are in wtw_df2.csv

loadData = function() {

  ### load trial-level data
  
  # load the concatenated trial data
  cat('Loading trial-level data.\n')
  d = read.csv('~/Box Sync/Project_wtw/R/wtw_df1.csv', stringsAsFactors=FALSE)
  # d = read.csv('../data/wtw_df1.csv', stringsAsFactors=FALSE)
  d$ID = as.character(d$ID) # store IDs as character
  
  # divide trial data by subject into a named list
  allIDs = unique(d$ID)
  n = length(allIDs)
  trialData = list() # initialize
  excludedIDs = c() # initialize
  for (sIdx in 1:n) {
    
    # pull this subject's trials
    thisID = allIDs[sIdx]
    rowIdx = (d$ID == thisID)
    thisTrials = d[rowIdx,]
    
    # add a column for the trial number
    thisTrials$trialNum = 1:nrow(thisTrials)
    
    # check for data completeness. intended session duration is 300 s
    # flag if the latest outcome is not later than 270 s. 
    lastOutcomeTime = max(thisTrials$outcomeTime)
    if (lastOutcomeTime > 270) {
      trialData[[thisID]] = thisTrials
    }
    else {
      cat('  WARNING: Subject',thisID,'has last outcome at',lastOutcomeTime,'s (session should be 300 s).\n')
      excludedIDs = append(excludedIDs, thisID)
    }
  }
    
  # determine n after exclusions
  n = length(trialData)
  cat(n,'complete data records loaded.\n')
  
  ### load subject-level data
  library(readr)
  library(readxl)
  
  # load the subject-level data
  # d = read_excel("../data/WTW DATA 05-15-18.xlsx")
  d = read_excel("~/Box Sync/Project_wtw/R/WTW DATA 05-15-18.xlsx")
  # d = read.csv('~/Box Sync/Project_wtw/R/wtw_df2.csv', stringsAsFactors=FALSE)
  d$ID = as.character(d$ID) # store IDs as character
  subjectData = d[,c('ID', 'group1245')]
  
  # remove rows for subjects who lacked complete trial data above
  for (i in 1:length(excludedIDs)) {
    thisExclID = excludedIDs[i]
    subjectData = subjectData[subjectData$ID != thisExclID,]
  }
  
  # return the 2 outputs in a named list
  allData = list(trialData=trialData, subjectData=subjectData)
  return(allData)
  
} # end of function




