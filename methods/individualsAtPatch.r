setGeneric('individualsAtPatch', function(patDF, indDF, setup) standardGeneric('individualsAtPatch'))
setMethod(f = 'individualsAtPatch', signature = c(patDF='matrix', indDF='matrix', setup='data.frame'), definition = function(patDF, indDF, setup){
	patDF[,'agentsPresent'] <- 0 # reset count
	patDF[,'competitionPayOff'] <- 0

	learn   <- indDF[,'learn']
	atPatch <- indDF[,'atPatch']
	if(sum( is.na(learn) ) != sum( atPatch>0) ) stop('Something goes wrong with the foraging, in or before the individualsAtPatch() function.')
	if(any( is.na(learn) ) && any( atPatch>0 )){
		nAtPatchAndID <- table(atPatch[ atPatch != 0])#  &  indDF$learn == 'not'])
		nAtPatch <- as.numeric(nAtPatchAndID)

    	if(length(nAtPatch) != 0){
      	patchID <- as.numeric(names(nAtPatchAndID))
  			where <- patDF[,'id'] %in% patchID
    		patDF[,'agentsPresent'][ where ] <- nAtPatch
    	}

			  if(setup$competition){
			  	payoff <- patDF[,'payoff']
			  	patDF[,'competitionPayOff'] <- payoff
			  	# for simulations where individuals can only take a limited amount of resources
			  	if(!is.null(setup$share)){
			  		fullPay <- setup$share*4*patDF[,'agentsPresent'][ where ] <= payoff[ where ]
			  		# if individuals take fewer resources than available, they get a full payoff
			  		patDF[,'competitionPayOff'][ where ][ fullPay ] <- setup$share*4
			  		# if individuals need more resources than available, they have to share equally
			  		patDF[,'competitionPayOff'][ where ][!fullPay ] <- payoff[ where ][!fullPay ] / ((patDF[,'agentsPresent'][ where ][!fullPay ]) ^ setup$competitionStrength)
			  	} else {
					# for all other simulations
			  		patDF[,'competitionPayOff'][ where ] <- payoff[ where ] / (patDF[,'agentsPresent'][ where ] ^ setup$competitionStrength) # represents an pi/(n^c) where c is the scaling factor of how severe number of individuals are
			  	}
			  }
	}
  return(patDF)
})
