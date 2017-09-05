### Patche IDs that were observed from exploiters in last round
setGeneric('observation', function(indDF, net, adjL) standardGeneric('observation'))
setMethod(f = 'observation', signature = c(indDF='matrix', net='missing', adjL='missing'), definition = function(indDF, net){
  observable <- indDF[, 'atPatch']
  if(any(observable!=0)){
  	tmp <- matrix(c(1:nrow(indDF),observable), ncol=2)
  	observable <- tmp[tmp[,2]!=0,]
  } else {
  	observable <- matrix(c(0,0), ncol=2)
  }
  return(observable)
})
