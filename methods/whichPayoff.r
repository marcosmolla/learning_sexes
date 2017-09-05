### This is the functino to calculate the competition payoff q, depending on the patch payoff Q, the number of individuals at the patch N and the competition strength c:
q_calc <- function(Q,N,c){return(Q /(N^c))}

p_calc <- function(pi, d, di){
				if(!any(round(d,5)==round(di,5))) stop("Input value for d does not occour again in di. If it is an individual in the same patch it really should.")
				n=length(di)
				c=1
				return( (d/sum(di) * (n*pi)/(n^c)) )
			}

setGeneric('whichPayoff', function(setup, patDF, indDF, indDF_repertoir_source, mod) standardGeneric('whichPayoff'))
setMethod(f = 'whichPayoff', signature = c(setup='data.frame', patDF='matrix', indDF='matrix', indDF_repertoir_source='ANY', mod='character'), definition = function(setup, patDF, indDF, indDF_repertoir_source, mod){
	indDF[,'yield2'] <- indDF[,'yield'] # fill yield2 with values from yield
  # Clean yield from last round
	indDF[,'yield'] <- 0
	# Any exploiter?
	atPatch <- indDF[,'atPatch']
	indID <- atPatch!=0
	# if at least one exploits
  if(sum(indID) > 0){
  	# use the correct payoffs depending on absence or presence of competition
		if(setup$competition){
			patPayoff <- patDF[,'competitionPayOff']
		} else {
			patPayoff <- patDF[,'payoff']
		}
		indDF[, 'yield'] <- patPayoff[( c(atPatch)+1 )] #this only works when patch ID 0 with payoff 0 is at poition 1!
  }
	return(indDF)
})
