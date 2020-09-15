loadAllData = function() {
  library("dplyr")
  library("gtools")
  # load hdrData
  hdrData = read.csv("../wtw_participantGrouping.csv")
  hdrData$group <-
    dplyr::recode(
     hdrData$GROUP1245,
      `1` = "Controls",
      `2` = "Depressed",
      `4` = "Ideators",
      `5` = "Attempters"
    )
  hdrData = hdrData[c('ID', 'group')]
  hdrData= hdrData %>% arrange(ID)
  hdrData$ID = as.character(hdrData$ID)
  
  
  # load trialData
  data = read.csv(file = "../wtw_df1.csv", header = T)
  names(data) = c("blockNum", "trialStartTime", "initialPos", "scheduledWait",
                  "trialResult", "timeWaited", "trialEarnings", "totalEarnings",
                  "timeLeft", "sellTime", "ID")
  ids = mixedsort(unique(data$ID))
  nSub = length(ids)
  

  excludedIDs = c()
  # loop over
  trialData = list()
  for(i in 1 : nSub){
    id = as.character(ids[i])
    tmp = data[data$ID == ids[i], ]
    # flag if the latest outcome is not later than 270 s. 
    lastOutcomeTime = max(tmp$sellTime)
    if (lastOutcomeTime > 270) {
      tmp$trialNum = 1 : nrow(tmp)
      trialData[[id]] = tmp
    }
    else {
      tmp$trialNum = 1 : nrow(tmp)
      trialData[[id]] = tmp
      cat('  WARNING: Subject',id,'has last outcome at',lastOutcomeTime,'s (session should be 300 s).\n')
      excludedIDs = append(excludedIDs, id)
    }
  }
  # exclude one participant with different experiment contigencies
  excludedIDs = append(excludedIDs, "207989")
  # exclude participants who didn't move the cursor at all
  # check 
  bad <- as.character(c(45709, 207789, 211245, 212926, 213704, 214542, 215842, 217283, 217293, 219675, 220244, 220947, 221273, 881137))
  excludedIDs = append(excludedIDs, bad)
  

  # apply the filter to hdrData
  trialData = trialData[!names(trialData) %in% excludedIDs]
  hdrData = hdrData[!hdrData$ID %in% excludedIDs,]
  # output
  outputData = list(hdrData = hdrData, trialData=trialData)
  return(outputData)
} 

loadExpPara = function(paraNames, dirName){
  # number of paraNames 
  nE = length(paraNames) + 1
  # number of files
  fileNames = list.files(path= dirName, pattern=("*_summary.txt"))
  library("gtools")
  fileNames = mixedsort(sort(fileNames))
  n = length(fileNames) 
  sprintf("load %d files", n)
  
  # initialize the outout variable 
  expPara = matrix(NA, n, nE * 4 + 1)
  idList = vector(length = n)
  # loop over files
  for(i in 1 : n){
    fileName = fileNames[[i]]
    address = sprintf("%s/%s", dirName, fileName)
    junk = read.csv(address, header = F)
    idIndexs = str_locate(fileName, "s[0-9]+")
    idList[i] = substr(fileName, idIndexs[1]+1, idIndexs[2])
    # delete the lp__ in the old version
    if(nrow(junk) > nE){
      junk = junk[1:nE,]
    }
    expPara[i, 1:nE] = junk[,1]
    expPara[i, (nE + 1) : (2 * nE)] = junk[,3]
    expPara[i, (2*nE + 1) : (3 * nE)] = junk[,9]
    expPara[i, (3 * nE + 1) : (4 * nE)] = junk[,10]
    expPara[i, nE * 4 + 1] = junk[1,11]
  }
  # transfer expPara to data.frame
  expPara = data.frame(expPara)
  junk = c(paraNames, "LL_all")
  colnames(expPara) = c(c(junk, paste0(junk, "SD"), paste0(junk, "Effe"), paste0(junk, "Rhat")), "nDt")
  expPara$id = idList # ensure the levels are consistent, usually not that problematic though
  return(expPara)
}

# I also need to load 2.5% and 97.5%
loadCVPara = function(paraNames, dirName, pattern){
  # number of paraNames 
  nE = length(paraNames) + 1
  # number of files
  fileNames = list.files(path= dirName, pattern= pattern)
  library("gtools")
  fileNames = mixedsort(sort(fileNames))
  n = length(fileNames) 
  sprintf("load %d files", n)
  
  # initialize the outout variable 
  expPara = matrix(NA, n, nE * 6)
  idList = vector(length = n)
  # loop over files
  for(i in 1 : n){
    fileName = fileNames[[i]]
    address = sprintf("%s/%s", dirName, fileName)
    junk = read.csv(address, header = F)
    sIndexs = str_locate(fileName, "s[0-9]+")
    s = substr(fileName, sIndexs[1]+1, sIndexs[2])
    fIndexs = str_locate(fileName, "f[0-9]+")
    f = substr(fileName, fIndexs[1]+1, fIndexs[2])
    idList[i] = sprintf("s%s_f%s", s, f)
    # delete the lp__ in the old version
    if(nrow(junk) > nE){
      junk = junk[1:nE,]
    }
    expPara[i, 1:nE] = junk[,1]
    expPara[i, (nE + 1) : (2 * nE)] = junk[,3]
    expPara[i, (2*nE + 1) : (3 * nE)] = junk[,9]
    expPara[i, (3 * nE + 1) : (4 * nE)] = junk[,10]
    expPara[i, (4*nE + 1) : (5 * nE)] = junk[,4]
    expPara[i, (5 * nE + 1) : (6 * nE)] = junk[,8]
  }
  # transfer expPara to data.frame
  expPara = data.frame(expPara)
  junk = c(paraNames, "LL_all")
  colnames(expPara) = c(junk, paste0(junk, "SD"), paste0(junk, "Effe"), paste0(junk, "Rhat"),
                        paste0(junk, "2.5"),paste0(junk, "97.5"))
  expPara$id = idList
  return(expPara)
}

