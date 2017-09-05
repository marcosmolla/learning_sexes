### Create individuals (pure/mixed strategies); setting innovate proportion
setGeneric('setInnovationProp', function(stratProp, nInd, indDF, net) standardGeneric('setInnovationProp'))
setMethod(f = 'setInnovationProp', signature = c(stratProp='numeric', nInd='numeric', indDF='matrix', net='missing'), definition = function(stratProp, nInd, indDF, net){
  if(stratProp==1) indDF[ ,'innovateProp'] <- 1
  if(stratProp==0) indDF[ ,'innovateProp'] <- 0
  if(all(stratProp!=c(0,1))) {
    if(stratProp==2){ # 2 is reserved as a code for random proportions
      indDF[ ,'innovateProp'] <- runif(n = nInd,min = 0,max = 1)
			indDF[ ,'innovateProp'][indDF[ ,'innovateProp'] < 0] <- 0
			indDF[ ,'innovateProp'][indDF[ ,'innovateProp'] > 1] <- 1
    } else { # any proportion of 0 and 1
    	id <- sample(c(1:nInd), size = nInd*stratProp, replace = F)
    	indDF[  , 'innovateProp'] <- 0
    	indDF[id, 'innovateProp'] <- 1
    }
  }
  return(indDF)
})
