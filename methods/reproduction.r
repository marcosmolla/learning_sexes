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
  m <- indDF[,"S"]==0 # 0 indicating males

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
    stop("Rule 'var' not defined. Check reproduction.r code!")
    # ordered <- order(indDF[m & alive,'fitness_variance'], decreasing=FALSE, na.last = NA) #orders by fitness values, in decreasing order because NA's are put at the bottom, largest . NA's are now removed.#This includes dead individuals!
    # minn <- round(sum(m&alive)*setup$x)
    # if(minn>length(ordered)) minn <- length(ordered) # in case there are individuals included with fitnes NA, reduce it to the number of above individuals with fitness != NA
    # lowestreproducingmale <- indDF[m&alive,'id'][ordered][ifelse(minn<1,1,minn)] #returns position of x-th lowest male.
    # male <- (indDF[,"fitness_variance"] <= indDF[lowestreproducingmale, 'fitness_variance']) & m & alive #returns true if fitness exceeds lowest reproducing male and that male is alive!
    }
  if(reprul!="fit"&reprul!="var") {stop("At the moment there is only fit and var implemented!")}

  who <- which(alive&male) # this is/are the one/ones who can potentially reproduce

  indDF <- reproducing(WHO=who, indDF=indDF, rule=reprul, id=ID, setup=setup)


return(indDF)
})


eligible <- function(x, threshold){
  if(length(x)==0){
    return(integer(0))
  } else {
    if(is.null(nrow(x))){ # that is the case when only one value was handed over
      return( as.numeric(x["id"]) )
    } else {
      # x[,"fitness"] %>% percent_rank() >= 1-threshold -> who
      x[,"fitness"] %>% order(decreasing=T) -> o
      id <- x[o,"id"][ 1:( length(x[,"id"])*threshold ) ]
      # stop("check this ranking thingy again")
      return( id )
    }
  }
}

# Define sexual reproduction function
# This function is specifically written for the sexual reproduction version of the model
setGeneric('sexreproduction', function(setup, indDF) standardGeneric('sexreproduction'))
setMethod(f = 'sexreproduction', signature = c(setup='data.frame', indDF='matrix'), definition = function(setup, indDF){
  # Individauls which are still alive
  alive <- indDF[,'nRound']>0
  
  # Determine eligible individuals for reproduction
  females <- eligible(x=indDF[alive & !is.na(indDF[,"fitness"]) & rowSums(indDF[,c("yield","yield2")])>=setup$minIncome & indDF[,"S"]==1, c("id", "fitness")], threshold=setup$xf)
    males <- eligible(x=indDF[alive & !is.na(indDF[,"fitness"]) &                                                         indDF[,"S"]==0, c("id", "fitness")], threshold=setup$xm) 
    # males <- eligible(x=indDF[alive & !is.na(indDF[,"fitness"]) & indDF[,"yield"]>=setup$minIncome & indDF[,"S"]==0, c("id", "fitness")], threshold=setup$xm) 

  # For reproduction, at least one male and one female is required
  if(length(males)>0 & length(females)>0){
    # In this case there can be reproduction
    # Choose Mothers (either with realtive to fitness, or with equal probabilities)
    if(setup$highestFitnessIsEverything) probF <- indDF[females, "fitness"] else probF <- NULL
    mothers <- mySample(n=females, size=sum(!alive), replace=T, prob=probF)
    # Choose Fathers (either with realtive to fitness, or with equal probabilities)
    if(setup$highestFitnessIsEverything) probM <- indDF[  males, "fitness"] else probM <- NULL
    fathers <- mySample(n=  males, size=sum(!alive), replace=T, prob=probM)
    
    # Now we will choose the genes inherited by the two parents, for each new individual 
    for(i in 1:sum(!alive)){
      ID <- which(!alive)[i]
      # Determine female genes (one from the mother, one from the father)
      indDF[ID, c("f1","f2")] <- c( mySample(indDF[mothers[i], c("f1","f2")], size=1), mySample(indDF[fathers[i], c("f1","f2")], size=1) )
      # Determine   male genes (one from the mother, one from the father)
      indDF[ID, c("m1","m2")] <- c( mySample(indDF[mothers[i], c("m1","m2")], size=1), mySample(indDF[fathers[i], c("m1","m2")], size=1) )
      
      # # Add some noise (ONLY FOR FREE INNOVATEPROP)
      # indDF[ID, c("f1","f2","m1","m2")] <- round(rnorm(n=4, mean=indDF[ID, c("f1","f2","m1","m2")], sd=0.1),2)
      # indDF[ID, c("f1","f2","m1","m2")][indDF[ID, c("f1","f2","m1","m2")]>1] <- 1
      # indDF[ID, c("f1","f2","m1","m2")][indDF[ID, c("f1","f2","m1","m2")]<0] <- 0
      # Mutation (ONLY FOR FIXED INNOVATEPROP)
      indDF[ID, c("f1","f2","m1","m2")] <- abs((runif(4,0,1)<setup$mutationRate) - indDF[ID, c("f1","f2","m1","m2")])
      
      # Determine sex of the newborn
      indDF[ID, "S"] <- sample(c(0,1), size=1)
      
      # Determine innovation proportion based on sex
      # if(indDF[ID, "S"]==1) iP <- round(mean(indDF[ ID , c("f1","f2")]), 2)  else iP <- round(mean(indDF[ ID , c("m1","m2")]), 2) # for free innnovateProp 
      if(indDF[ID, "S"]==1) iP <- round(mean(indDF[ ID , c("f1","f2")]), 1)  else iP <- round(mean(indDF[ ID , c("m1","m2")]), 1) # for fixed innovateProp
      indDF[ ID ,"innovateProp"] <- iP
      
      # Reset data frame to accomodate newborn
      indDF[ ID , which(!colnames(indDF)%in%c('id','innovateProp','S',"m1","m2","f1","f2")) ] <- 0
      indDF[ ID ,"nRound" ] <- 1
    }
    
  } else {
    indDF[ !alive ,"nRound" ] <- -100
  }
  
  return(indDF)
})
