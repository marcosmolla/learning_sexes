### Individual choose one of two behaviors
selectGOTO <- function(k){
	if(sum(!is.na(k)) == 0) {return(0)} else {
	        res <- which(x=k==max(k, na.rm = T))
	        if(length(res)!=1) res <- sample(x = res, size = 1, replace = F)
	        return(res-1) # otherwise i get the position in the data.frame and not the patch ID
	      }
}

setGeneric('whichAction', function(indDF, nInd, indDF_repertoir, patDF, allExploit, setup) standardGeneric('whichAction'))
setMethod(f = 'whichAction', signature = c(indDF='matrix', nInd='numeric', indDF_repertoir='ANY', patDF='ANY', allExploit='ANY', setup='data.frame'), definition = function(indDF, nInd, indDF_repertoir, patDF, allExploit, setup){
  ### First everyone is a learner
	indDF[,'atPatch'] <- 0
	atPatch <- indDF[,'atPatch']

	# Which propensity to exploit
	ifelse(is.null(setup$beta), beta <- .8, beta <- setup$beta)
	# Who is exploiting
  who <- runif(n = nInd,min=0,max=1) <= beta & indDF[, 'nRound']!=1

  ID <- which(who)
	# Any exploiting individuals?
  if(length(ID) > 0){
  	# Only one is exploiting
		if(length(ID)==1){
			goto <- selectGOTO(k=indDF_repertoir[ID,])
		}else{
		# Several are exploiting
	    goto <- apply(indDF_repertoir[ID,],1,selectGOTO)
		}
		# Whiich pates chosen?
		atPatch[ID] <- goto
  }
	# Feedback information into indDF
  indDF[ ,'atPatch'] <- atPatch
  indDF[ indDF[,'nRound']<0 ,'atPatch'] <- 0 # keep dead indivduals from exploiting patches
  return(indDF)
})
