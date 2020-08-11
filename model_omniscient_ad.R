omniscient <- function(tGrid,plots)

reward_cdf = data.frame(time=c(0, 1, 2, 3, 20), cdf=c(0, 0.25, 0.50, 0.75, 1)) # ground truth
# tGrid = seq(0, 20, by=0.1) # grid of quitting times to evaluate
iti = 2
rwd = list("large"=10, "small"=1)

### expected total earnings ###

# initialize outputs
expectedRwdPerTrial = numeric(length = length(tGrid))
expectedTimePerTrial = numeric(length = length(tGrid))



# loop over time points
for (idx in 1:length(tGrid)) {
  tNow = tGrid[idx]
  
  # expected reward
  earlierRwdTimes = reward_cdf$time <= tNow
  pRwd = max(c(0, reward_cdf$cdf[earlierRwdTimes]))
  expectedRwdPerTrial[idx] = pRwd*rwd$large + (1 - pRwd)*rwd$small
  
  # expected time
  delay_cdf = data.frame(time = c(reward_cdf$time[earlierRwdTimes], tNow),
                         cdf = c(reward_cdf$cdf[earlierRwdTimes], 1))
  expectedTimePerTrial[idx] = sum(diff(delay_cdf$cdf) * tail(delay_cdf$time,-1)) + iti
}

# expected total reward in 5 min
expectedTotalReward = 5 * 60 * expectedRwdPerTrial / expectedTimePerTrial
if (plots)
{plot(tGrid, expectedTotalReward, type="l", main="Expected total earnings", 
     xlab="Quitting time (s)", ylab="Earnings (cents)", bty="n", xlim=c(0,20), ylim=c(0,600),
     col='red', lwd=3)
}
# long-run average reward rate for the reward-maximizing policy
# (this will be used below for the optimal opportunity cost)
LRARR = max(expectedRwdPerTrial / expectedTimePerTrial)
cat("Optimal long-run average reward rate (LRARR) =",LRARR,"cents per s\n")

### function to calculate the value of waiting at one point in time ###

# inputs: 
#   the forward-looking reward-time distribution (dependent on time elapsed); 
#     this is a data from with fields $time and $cdf, as above;
#   the estimated LRARR (this factors in the ITI);
#   the reward magnitudes.
waitValue <- function(updated_cdf, LRARR, rwd) {
  # define the grid of possible waiting times to evaluate
  tGrid = seq(0, round(max(updated_cdf$time),1), by=0.1)
  # upsample the cdf to the time grid (could probably be done more efficiently)
  upsampled_cdf = numeric(length = length(tGrid))
  for (idx in 1:length(upsampled_cdf)) {
    earlierRwdTimes = updated_cdf$time <= tGrid[idx]
    upsampled_cdf[idx] = max(c(0, updated_cdf$cdf[earlierRwdTimes]))
  }
  # convert the cdf to expected reward for each possible quitting time
  expectedRwd = upsampled_cdf*rwd$large + (1 - upsampled_cdf)*rwd$small
  # expected delay under each quitting policy
  upsampled_survFx = c(1, 1 - head(upsampled_cdf, n=-1)) # survival function
  marginal_time_cost = c(tGrid[1], diff(tGrid)) * upsampled_survFx
  expectedTime = cumsum(marginal_time_cost)
  # net value of waiting at each point
  # n.b. where there is a non-zero immediate reward, this is part of the OC of waiting
  netValue = expectedRwd - expectedTime * LRARR - rwd$small
  # potential value of waiting
  potential = max(netValue)
  return(potential)
}

### iterate through time steps, calculating the forward-looking value of waiting at each step ###

tGrid = head(tGrid, n=-1) # examine all but the last time point (no decision to be made there)
sv = numeric(length = length(tGrid)) # initialize results
sv_loLRARR = numeric(length = length(tGrid)) # if LRARR were underestimated
for (idx in 1:length(tGrid)) {
  tNow = tGrid[idx]
  # update the reward-time cdf based on time elapsed
  # n.b. includes rounding b/c imprecision messes up logical tests.
  futurePoints = reward_cdf$time > tNow
  pElapsed = max(c(0, reward_cdf$cdf[!futurePoints])) # for renormalizing the future part of the function
  updated_cdf = data.frame(time = round(reward_cdf$time[futurePoints] - tNow, 1),
                           cdf = (reward_cdf$cdf[futurePoints] - pElapsed) / (1 - pElapsed))
  sv[idx] = waitValue(updated_cdf, LRARR, rwd)
  sv_loLRARR[idx] = waitValue(updated_cdf, LRARR/2, rwd)
}
# plot results
if (plots)
{plot(tGrid, sv, type="l", main="Theoretical value of waiting", 
     xlab="Time elapsed (s)", ylab="Net value of waiting (cents)", bty="n", xlim=c(0,20), 
     col='blue', lwd=3)
lines(tGrid, sv_loLRARR, col='gray', type='l', lwd=3, lty=2)
legend("topleft", c("Optimal LRARR", "Underestimated LRARR"), bty="n", lty=c(1, 2), lwd=3, 
       col=c("blue", "gray"))
}
# save value functions in survival format
# v <- as.data.frame(sv)
# v$waitValue <- head(expectedTotalReward, n=-1)
# v$t1 <- tGrid
# v$t2 <- tGrid + .1
# save(v,file = 'vFunc.RData')



