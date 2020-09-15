condition = "LP"
delayMax = 20 # max trial durations in secs
nBlock = 1
blockMin = 5 # block duration in mins
blockSec = blockMin * 60 # block duration in secs
largeReward = 10 # value of the matured token
smallReward = 1 # value of the immature token
# possible values of reward dlays
rewardDelays = list(LP = c(1, 2, 3, 20)) 
# analyses parameters
tGrid = seq(0, blockSec-1, by = 1) # time grid for wtw time courses
kmGrid = seq(0, delayMax, by = 0.1) # time grid for Kaplan-Meier survival curves
save("condition" = condition,
     "delayMax" = delayMax,
     "blockMin" = blockMin,
     "blockSec" = blockSec,
     "nBlock" = nBlock,
     "smallReward" = smallReward,
     "largeReward" = largeReward,
     "rewardDelays" = rewardDelays,
     'tGrid' = tGrid,
     'kmGrid' = kmGrid,
     file = "expParas.RData")
