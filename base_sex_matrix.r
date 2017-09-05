###
### Marco Smolla, Manchester University
### 15-05-2014


#########################################
##     #####   ##  ##  ##  ##   ###    ##
##     ##  ##  ##  ##  ### ##   ###    ##
##     #####   ##  ##  ######   ###    ##
##     ## ##   ##  ##  ## ###          ##
##     ##  ##   ####   ##  ##   ###    ##
#########################################

### Set up
# determine the queue number of this process (necessary when saving the results)
### Load world parameters

# BEGIN CONDOR Computer Cluster
### Set-up packages
args <- commandArgs(TRUE)
if(length(args)!=0){ # CSF hands over a process number, if there is a queue variable we are in CONDOR, otherwise we work locally; Note: in CSF the varaiables are already in the right shape
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }


  ### Load methods
  if(file.exists('methods.tar')) untar(tarfile="methods.tar", exdir=".")
  setwd('./methods')
  files<-list.files()
  lapply(files, source)
  setwd('../')
  if(file.exists('methods.tar'))	unlink('methods',recursive=TRUE)
  ### load additional packages
  library(dplyr)

} else { #END CONDOR
  ### For local simulations:
  ### Load methods
  setwd('~/Dropbox/Programming/R/learning_sexes/methods/')
  l <- lapply(list.files(), source)
  ### set queue to 1
  queue <- 1
  ### load additional packages
  library(parallel)
  library(dplyr)
}


#### Setting parameters for simulations on cluster >>>>>>>>>>>>>>>>>
variances <- c(1,100) # payoff variance
mean <- 4 # the average payoff in the world
k_seq <- (mean^2)/variances # k parameter for the gamma distribution

# sequence for the environmental payoffChangeRate values
e_seq <- 10^seq(-3,0,by=.5) # environmental turnover

# create all possible combinations of parameters
grid <- expand.grid(k=k_seq,e=e_seq, x=c(0.01,0.02,0.04,0.08,seq(from=0.1,to=1,by=.1)), p=c(0.5), h=c(T,F), r=c(.20), c=c(1), rep=1:10)

# choose individual parameters for the current queue
K <- grid[queue, 'k']
PCR <- grid[queue, 'e']
X <- grid[queue, 'x']
HFIE <- grid[queue, 'h']
PROP <- grid[queue, 'p']
R <- grid[queue, 'r']
COMP <- grid[queue, 'c']


### Defining run routine
run <- function(setup, ceed){
  # Create and store seed information
  ceed <- round(runif(1,0,99999999))
  set.seed(ceed)
  setup$seed <- ceed

  ### Defining the data objects
  patDF <- as.matrix(data.frame(id=0:(setup$nPat), payoff=0, agentsPresent=NA, competitionPayOff=NA))#;patDF #formerly known as world@patch ### Patches in the world
  indDF <- as.matrix(data.frame(id=1:setup$nInd, nRound=1, nExploit=0, innovateProp=NA, atPatch=0, learn=NA, fitness=0, yield=0, yield2=0, S=NA, fitness_variance=0))#;indDF ### Single data to be stored for every living individual (row-wise)

  ### Multiple data stored for every individual (column-wise)
  indDF_repertoir <- matrix(NA,nrow=setup$nInd, ncol=(setup$nPat+1)) #nrow is max 100 because all individuals die
  indDF_collectedPayoff <- matrix(NA, ncol=setup$maxAge, nrow=setup$nInd)#columns are rounds, rows are indviduals

  patDF <- changeWorld(patDF=patDF, pat=nrow(patDF), create = T, payoffDis = as.character(setup$payoffDis), mod = '', setup=setup)# ;patDF
  indDF[,'S'] <- 0
  indDF <- setInnovationProp(stratProp=setup$strategyProportion, nInd=setup$nInd, indDF=indDF)# ;indDF

  time <- 0
  allSame <- length(unique(indDF[,"innovateProp"]))==1

  while(!allSame){
    time <- time+1
    patDF           <- changeWorld(patDF=patDF, pat=nrow(patDF), create = F, payoffDis = as.character(setup$payoffDis), mod='', setup=setup) #;patDF
    observable 			<- observation(indDF) #;observable
    indDF 					<- whichAction(indDF = indDF, nInd = setup$nInd, indDF_repertoir = indDF_repertoir, patDF=patDF, setup=setup) #;indDF
    indDF 					<- whichLearning(nInd = setup$nInd, indDF = indDF) #;indDF
    indDF_repertoir <- formingMemory(setup = setup, patDF = patDF, indDF = indDF, indDF_repertoir = indDF_repertoir, observable = observable, mod='') #;indDF_repertoir
    patDF           <- individualsAtPatch(patDF = patDF, indDF = indDF, setup=setup) #;patDF
    indDF           <- whichPayoff(setup = setup, patDF = patDF, indDF = indDF, mod = '') #;indDF
    indDF_repertoir <- updateMemory(indDF_repertoir = indDF_repertoir, indDF = indDF, mod=setup$mod) #;indDF_repertoir
      alive <- which(indDF[,"nRound"]>0)
      learn <- which( !is.na(indDF[,'learn']) & (indDF[,"nRound"] > 0) )
      indDF_collectedPayoff[ cbind(alive,indDF[alive,"nRound"]) ] <- indDF[alive,"yield"]
      indDF_collectedPayoff[ cbind(learn,indDF[learn,"nRound"]) ] <- NA #;indDF_collectedPayoff[,1:10] # use NAs to indicate learning turns; otherwise they overly increase the variance
	indDF           <- updatePayoffs(indDF = indDF, indDF_collectedPayoff=indDF_collectedPayoff, setup = setup) #;indDF
    indDF           <- ageAndDie(setup = setup, indDF = indDF, mod='') #;indDF
    if(sum(indDF[,'nRound']<0) > 0){
      for(h in which(indDF[,'nRound']<0)){
        indDF_repertoir[h,1:(setup$nPat+1)] <- NA
        indDF_collectedPayoff[h,] <- NA
      }
        indDF <- reproduction(setup = setup, indDF = indDF, mod='') #;indDF
    }

    allSame <- length(unique(indDF[,"innovateProp"]))==1
    if(time == setup$nRound) allSame <- T

  }

  ### Report additional values from the simulation
  patchesIL <- which(colSums(!is.na(indDF_repertoir[indDF[,"innovateProp"]==1, ]))>0) # average known patch values
  patchesSL <- which(colSums(!is.na(indDF_repertoir[indDF[,"innovateProp"]==0, ]))>0) # average known patch values
  payIL <- ifelse(length(patchesIL)>0, mean(patDF[patchesIL, "payoff"]), 0)
  paySL <- ifelse(length(patchesSL)>0, mean(patDF[patchesSL, "payoff"]), 0)
  diffpatIL <- ifelse(length(patchesIL)>0, length(patchesIL), 0)
  diffpatSL <- ifelse(length(patchesSL)>0, length(patchesSL), 0)

  return(cbind(indDF, setup, payIL, paySL, diffpatIL, diffpatSL))
} # END RUN

### Define simulation parameters
setup <- data.frame(
  nRound = 5,#5000,
  nInd = 10,
  nPat = 10,
  payoffChangeRate = PCR,
  payoffDis = "gamma", #expo or normal = uneven or even.
  strategyProportion = PROP, #how many Individual Learner/Learning
  competition = TRUE,
  deathRate = 1/100,#0.01,
  maxAge = 100, #((1/50)^-1)*2
  propStable = FALSE,#FALSE, #Freely evolving = false.
  sexRatioStable = TRUE, #Freely evolving = false
  mutationRate = 0.001,#2,
  repetition = 1,
  k = K, # shape
  theta = 4/K, # scale
  lethalILearning = 0,
  competitionStrength = COMP,
  mod = '',
  philopatry = FALSE,
  centralIntel = FALSE,
  #share = SHARE, #exclude when not in use !
  SexStratProp = 0, #1 == 100% female strategy, 0 == 100% male strategy
  x = X, #proportion of individuals with highest fitness that can actually reproduce (for D strat)
  reprorule = "fit",
  highestFitnessIsEverything = HFIE,
  # sexualreproduction = FALSE,
  temporalDiscRate = R,
  stringsAsFactors = FALSE # thats an argument to data.frame()
)

# Call simulation to run (repeatedly)
res <- bind_rows(lapply(1:10, run, setup=setup)); table(res[,"innovateProp"]); sum(res[,"innovateProp"]==1)/nrow(res)
# Move working directory to output
if(!file.exists('output')) dir.create(path = 'output'); setwd('output')
# Store results as individual files
save(res, setup, file=paste(paste(format(Sys.time(), '%y%m%d'),'learning_sexes',queue, sep="_"), sep='/'))
