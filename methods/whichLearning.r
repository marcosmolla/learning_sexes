### Individuals learn
setGeneric('whichLearning', function(nInd, indDF) standardGeneric('whichLearning'))
setMethod(f = 'whichLearning', signature = c(nInd='numeric', indDF='matrix'), definition = function(nInd, indDF){
#   browser()
	aP <- indDF[ ,'atPatch']
	indDF[ ,'learn'][aP != 0] <- NA ## == 'not' #(!) # Individuals that go to a patch, do not learn
  indDF[ ,'learn'][aP == 0] <- 0 ## =='imitate' # all other individuals learn
  innovate <- runif(nInd,0,1) <= indDF[ ,'innovateProp']
  imitator <- indDF[ ,'learn']==0 ## 'imitate'
    indDF[ ,'learn'][innovate & imitator] <- 1 ## =='innovate'
	return(indDF)
})