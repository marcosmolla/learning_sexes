### Define reproducing function
reproducing <- function(WHO, indDF, id, setup, rule){
  if(length(WHO)==0){ #This can happen if individuals die but reproductive constraints prevent reproduction.
    indDF[id,'yield'] <- indDF[id,'yield2'] <- indDF[id,'fitness'] <- NA #Keeps individuals clear from means and sums if na.rm=TRUE.
    indDF[id,'atPatch'] <- 0
    indDF[id,'nRound'] <- -100 #Dead individuals remain in the data.frame so that their space can be later taken by a reproducing individual's offspring.
  }else{

    if(setup$propStable){
      deadStrats <- indDF[,'innovateProp'][id]
      deadSex <- indDF[,'S'][id]
      newStrats <- deadStrats
      newSex <- deadSex
    } else {

      if(length(WHO)==1){ #sampling from a bag with one number, x number of times, does not give that number x times (constraint of function). Therefore do a different strategy if just one individual reproducing.
        newStrats <- rep(indDF[WHO, 'innovateProp'],times=length(id)) #repeat innovateProp of sole reproducer across number that died in last round.
        newSex <- rep(indDF[WHO, 'S'],times=length(id)) #repeat Sex of sole reproducer across number that died in last round
      } else {
        # ### EITHER: sample according to fitness (relative income) >>>>>
        if(setup$highestFitnessIsEverything){
          realWHO <- mySample(WHO, size=length(id), replace=TRUE, prob=indDF[WHO, 'fitness']+0.00000001) #control for where everyone can reproduce if old enough and relative to their fitness
          } else {
        # ### <<<<<
        # ### OR: sample randomly from those how are eligible for reproduction (minimum income) >>>>>>>>
          realWHO <- mySample(WHO, size=length(id), replace=TRUE)
          }
        # ### <<<<<
        newStrats <- indDF[realWHO,'innovateProp']
        newSex <- indDF[realWHO,'S']
      }
      }
  ### MUTATION
  if( setup$mutationRate > 0 ){
    mut <- runif(n = length(id), min = 0, max = 1) <= setup$mutationRate
    if( any(mut) ){
      if(setup$strategyProportion==2){
        newStrats[mut] <- unlist(lapply(newStrats[mut], rnorm, n=1, sd=.1))
        newStrats[newStrats < 0] <- 0
        newStrats[newStrats > 1] <- 1
      } else {
        newStrats[mut] <- 1-newStrats[mut]
      }
    }
  }

  ### Feeding back into indDF
  indDF[id,] <- do.call(rbind, lapply(1:length(id), function(h) {
    tmp <- rep(0, ncol(indDF))
    tmp[which(colnames(indDF)%in%c('id','innovateProp','S'))] <- c(id[h],newStrats[h],newSex[h]) #feeds these new values into list of zeros (tmp) applied to reset dead individuals in turn (using lapply)
    tmp[which(colnames(indDF)%in%c('nRound'))] <- 1 # starting age with 1 instead of 0
    return(tmp)
  }))

  }
  return(indDF)
}



# Define reproduction function
setGeneric('reproduction', function(setup, indDF, mod) standardGeneric('reproduction'))
setMethod(f = 'reproduction', signature = c(setup='data.frame', indDF='matrix', mod='ANY'), definition = function(setup, indDF, mod){
	ID <- which(indDF[,'nRound']<0) # select individuals that were marked as dead by ageAndDie()
 	alive <- indDF[,'nRound']>0 # individauls which are still alive
  m <- indDF[,"S"]==0 # indicating males (note, this is legacy)

  # Reproductive Constraints
  reprul <- setup$reprorule

  if(reprul=='fit'){ # selecting by highest variance
    ordered <- order(indDF[m & alive,'fitness'], decreasing=TRUE, na.last = NA) #orders by fitness values, in decreasing order because NA's are put at the bottom, largest . NA's are now removed.#This includes dead individuals!
    minn <- round(setup$nInd*setup$x)
    if(minn>length(ordered)) minn <- length(ordered) # in case there are individuals included with fitnes NA, reduce it to the number of above individuals with fitness != NA
    lowestreproducingmale <- indDF[m&alive,'id'][ordered][ifelse(minn<1,1,minn)] #returns position of x-th lowest male.
    male <- (indDF[,"fitness"] >= indDF[lowestreproducingmale, 'fitness']) & m & alive #returns true if fitness exceeds lowest reproducing male and that male is alive!
  }

  if(reprul=='var'){ # selecting by highest variance
    ordered <- order(indDF[m & alive,'fitness_variance'], decreasing=FALSE, na.last = NA) #orders by fitness values, in decreasing order because NA's are put at the bottom, largest . NA's are now removed.#This includes dead individuals!
    minn <- round(sum(m&alive)*setup$x)
    if(minn>length(ordered)) minn <- length(ordered) # in case there are individuals included with fitnes NA, reduce it to the number of above individuals with fitness != NA
    lowestreproducingmale <- indDF[m&alive,'id'][ordered][ifelse(minn<1,1,minn)] #returns position of x-th lowest male.
    male <- (indDF[,"fitness_variance"] <= indDF[lowestreproducingmale, 'fitness_variance']) & m & alive #returns true if fitness exceeds lowest reproducing male and that male is alive!
    }
  if(reprul!="fit"&reprul!="var") {stop("At the moment there is only fit and var implemented!")}

  who <- which(alive&male) # this is/are the one/ones who can potentially reproduce

  indDF <- reproducing(WHO=who, indDF=indDF, rule=reprul, id=ID, setup=setup)


return(indDF)
})
