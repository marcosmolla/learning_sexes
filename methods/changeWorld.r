### Create and Change world payoffs using different distributions
setGeneric('changeWorld', function(patDF, pat, create, payoffDis, mod, setup) standardGeneric('changeWorld'))
setMethod(f = 'changeWorld', signature = c(patDF='matrix', pat='numeric', create='logical', payoffDis='character', mod='character', setup='data.frame'), definition = function(patDF, pat, create, payoffDis, mod, setup){
	change <- is.na(patDF[,'payoff'])
	if(sum(change)==0 & mod==''){
	  base <- T
	  change <- runif(n = pat,min = 0,max = 1) <= setup$payoffChangeRate
	  if(create) change <- rep(x = TRUE, pat) # This is for initialisation of the world
	}

  if(payoffDis == 'gamma')  patDF[change,'payoff'] <- rgamma(n = sum(change), shape=setup$k, scale=(setup$theta<-(4/setup$k)) )
  if(base) {patDF[1,'payoff'] <- 0}else{ #no payoff when learning happens in base and base_network model (there is a patch "0")
  	patDF[is.na(patDF[,'id']), 'payoff'] <- 0 # set payoff of the unknown patches (ID == NA) to 0
  }
  return(patDF)
})
