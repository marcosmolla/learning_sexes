setGeneric('mySample', function(n, size=1, replace=F, prob=NULL) standardGeneric('mySample'))
setMethod(f = 'mySample', signature = c(n='ANY', size='ANY', replace='ANY', prob='ANY'), definition = function(n, size=1, replace=F, prob=NULL){
  if(length(n)>1)  res <- sample(x = n, size = size, replace = replace, prob = prob)
  if(length(n)==1) res <- rep(n, times=size)
  if(length(n)==0) res <- NA
  return(res)
})
