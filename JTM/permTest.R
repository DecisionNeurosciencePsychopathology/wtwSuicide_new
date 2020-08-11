
library("permute")

runPermTest <- function(d, grp, tag='blank', nIter=5000) {
  # d is a vector of data values; grp is vector of group labels; tag is a descriptive string
  cat('Permutation test by group:',tag,'\n')
  n = length(d)
  if (n != length(grp)) {cat('*** data and labels are different lengths ***\n'); stop()}
  cat(n, 'data values,', length(unique(grp)), 'unique group labels,', nIter, 'permutations.\n')
  grp = as.factor(grp)
  dRank = rank(d) # rank-transformation
  
  set.seed(NULL)
  r2 = numeric(length=nIter)
  h = lm(dRank ~ grp) # unpermuted results
  r2[1] = summary(h)$r.squared # unpermuted results
  for (i in 2:nIter) {
    h = lm(dRank ~ grp[shuffle(n)])
    r2[i] = summary(h)$r.squared
  }
  
  pval = mean(r2 >= r2[1])
  plot(ecdf(r2), main=sprintf("Perm test: %s (p = %1.5f)",tag,pval), 
       ylab="Cumulative null distribution", xlab="r^2")
  lines(r2[1]*c(1,1), c(0,1), col='red', type='l', lwd=2)
  
}



runPermTest(subjectData$nCursorMoves, subjectData$GROUP1245, 'Cursor movements', 5000)

runPermTest(subjectData$grpAUC_nonImmed, subjectData$GROUP1245, 'AUC nonImmed', 5000)

runPermTest(subjectData$nImmedQuits, subjectData$GROUP1245, 'nImmedQuits', 5000)

runPermTest(subjectData$grpEarnings, subjectData$GROUP1245, 'Total earnings', 5000)

