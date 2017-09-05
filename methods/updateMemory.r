### An individual that forages at a given patch, receives current information about the offered reward. This can be different from what is stored in memory, it is updated here.
setGeneric('updateMemory', function(indDF_repertoir, indDF, mod) standardGeneric('updateMemory'))
setMethod(f = 'updateMemory', signature = c(indDF_repertoir='matrix', indDF='matrix', mod='character'), definition = function(indDF_repertoir, indDF,mod){
	atPatch <- indDF[,'atPatch']
	if(any(atPatch!=0)){
		ID <- indDF[,'id'][ atPatch!=0 ]
		j <- cbind( ID , atPatch[ (ID) ]+1 )
		indDF_repertoir[ j ] <- indDF[,'yield'][ ID ]
		if(mod=='collectiveKnowledge') indDF_repertoir[ , atPatch[ (ID) ]+1 ] <- rep(indDF[,'yield'][ ID ], each=nrow(indDF_repertoir))
	}
	return( indDF_repertoir )
})
