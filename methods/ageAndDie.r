setGeneric('ageAndDie', function(setup, indDF, mod) standardGeneric('ageAndDie'))
setMethod(f = 'ageAndDie', signature = c(setup='data.frame', indDF='matrix', mod='ANY'), definition = function(setup, indDF, mod){
  # Who is exploiting
  exploiter <- indDF[,'atPatch']!=0
  # Who is going to die
	die <- runif(n = setup$nInd, min = 0, max = 1) <= setup$deathRate | indDF[,'nRound']==setup$maxAge
  # Increase age
  whoDies <- indDF[,'nRound']+1
  # Mark dead individual with -1 flag
  whoDies[die] <- -1
  # Store age values
  indDF[,'nRound'] <- whoDies
  return(indDF)
})
