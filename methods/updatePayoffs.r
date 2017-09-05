## Collected payoffs in a given round are now added to the amount of payoffs collected over the whole life time ('collectedPayoff'), and subsequently the fitness is calculted

temporalDiscount <- function(x, windowSize, meanFit, setup){
  payoffs <- x[!is.na(x)]
  r = setup$temporalDiscRate
  gamma = 1
  if((j<-length(payoffs))>0){
    if(j<windowSize){
      J <- 0:(j-1)
      m <- mean( exp(-(r*J)) * payoffs[j:1]^gamma )
    }else{
      y<-(j-(windowSize-1))
      J <- 0:(windowSize-1)
      m <- mean(exp(-(r*J)) * payoffs[j:y]^gamma)
    }
    return(c(fit_sum=m, fit_var=NA))
  } else {
    return(c(fit_sum=NA, fit_var=NA))
  }
}

setGeneric('updatePayoffs', function(indDF, indDF_collectedPayoff, setup) standardGeneric('updatePayoffs'))
# New version to calculate fitness, including a temporal discount
setMethod(f = 'updatePayoffs', signature = c(indDF='matrix', indDF_collectedPayoff='matrix', setup="data.frame"), definition = function(indDF, indDF_collectedPayoff, setup){
  # Increment exploit turns for exploiters
  indDF[,'nExploit'] <- indDF[,'nExploit'] + as.numeric(indDF[,'atPatch']!=0)
  # Who is alive? (exclude 0 for fitness calculation to work)
  rueq0<-indDF[,'nRound']>0
  popMeanFit <- mean(indDF[,"fitness"], na.rm=T)
  tmp <- do.call(rbind,lapply(which(rueq0), function(i){
    dur <- indDF[i,"nRound"]
    temporalDiscount(x=indDF_collectedPayoff[i,1:dur], windowSize = 5, meanFit=popMeanFit, setup=setup)
  }))#;tmp
  tmp[,"fit_sum"][is.nan(tmp[,"fit_sum"])] <- NA
  indDF[,'fitness'][rueq0]<- tmp[,"fit_sum"]
  indDF[,'fitness_variance'][rueq0]<- tmp[,"fit_var"]
  return(indDF)
})
