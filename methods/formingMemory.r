setPayoffs <- function(setup, patDF){
	if(setup$competition && sum(patDF[,'competitionPayOff'],na.rm = T)!=0) {
		patPayoff <- patDF[,'competitionPayOff']
	} else {
		patPayoff <- patDF[,'payoff']
	}
	return(patPayoff)
}

choosePatches <- function(setup, patPayoff, observable, ID2, patDF){
	# Standard or Preferential Learning?
	# If setup$preferentialLearning does not exist Standard is used by default
	ifelse(test=!is.null(setup$preferentialLearning), yes=prefL<-setup$preferentialLearning, no=FALSE)
	observed <- unlist(lapply(observable[ID2], function(what){
		if(!prefL){
			## STANDARD LEARNING ===
			## the follwoing line is the standard sampling of any patch informatin that is not 0
			mySample(what[what!=0], size = 1, replace = T)
		} else {
			## PREFERENTIAL LEARNING ===
			## with the following line I alter the standard selection process so that individual always choose the individual that provides patch information with the higher payoff
			# First, extract payoffs (+1, because first patch is patch 0, which doesn't count)
			# pay<-patPayoff[what] # this doesn't work when ID's are mixed in the patDF, I need to find the right row in the data.frame that corresponds to the patch ID the learner wants to learn about
			pay <- patPayoff[match(what,patDF[,'id'])]
			# Second, if there are any patches known by neighbours and they weren't previously learning, choose the the patch with the highest payoff, if several are equal in payoff pick one at random
			ifelse(test=length( what[what!=0] )!=0, yes=mySample(what[which(pay==max(pay,na.rm=T))],1), no=NA)
		}
	}))
	# in the next step this becomes 0, this is the case when there was nothing to observe and mySample returns an NA, than this individual learned about patch 0
	observed[is.na(observed)] <- 0
	return(observed)
}






setGeneric('formingMemory', function(setup, patDF, indDF, indDF_repertoir, indDF_repertoir_source, observable, mod, net) standardGeneric('formingMemory'))
### Version for base and base_network
setMethod(f = 'formingMemory', signature = c(setup='data.frame', patDF='matrix', indDF='matrix', indDF_repertoir='ANY', indDF_repertoir_source='ANY', observable='ANY', mod='character', net='ANY'), definition = function(setup, patDF, indDF, indDF_repertoir, indDF_repertoir_source, observable, mod, net){

	### Consolidate repertoir due to INNOVATION
	patPayoff <- setPayoffs(setup, patDF)
  learn <- indDF[,'learn']
  ID <- which(learn == 1 & !is.na(learn)) #'innovate'

  if(length(ID) > 0){
    patches <- c(1:setup$nPat)
		learnAboutPatch <- sample(x = patches, size = length(ID), replace = TRUE)
	  if(sum(learnAboutPatch,na.rm = T) > 0){
	  	k <- cbind(ID, (learnAboutPatch+1))
	  	indDF_repertoir[k] <- patPayoff[ (learnAboutPatch+1) ]
	  }
  }

  ### Consolidate repertoir due to IMITATION
  ID2 <- indDF[,'id'][ learn == 0 & !is.na(learn)] #'imitate'

	# for standard model
	  if(is.matrix(observable)) token<-sum(observable[,2], na.rm = T) > 0 else token<-observable[2] > 0 #either matrix or vector
	  if(length(ID2) > 0 && token){ # if someone observes and something is observable
		  	observed <- unlist(lapply(ID2, function(who){
					if(!is.matrix(observable)){ #there is only one line of results and it is transformed to a vector
						ifelse(who==observable[1], return(0), return(observable[2]))
					} else { #there is more than one line, a matrix is evaluated
			  		x <- observable[matrix(observable,ncol=2)[,1]!=who,2] # with this piece I avoid picking information from the pool that the observer added in the last round (as a demonstrator)
			  		mySample(n = x, size = 1, replace = T)
					}
		  	}))
		  	observed[is.na(observed)] <- 0
	  }

if(exists('observed')){ #for network and standard model
	  j <- cbind(ID2, (observed+1))
		indDF_repertoir[j] <- patPayoff[ (observed+1) ] # now individuals learn which payoff the patch has, i.e. if there are already individuals they learn competition payoff
	}


return(indDF_repertoir)

})
