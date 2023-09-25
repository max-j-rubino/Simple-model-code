Replicates = function(lowsuccess.V, highsuccess.V,qual.lo.V,qual.hi.V) {  
  
  replicates = expand.grid(lowsuccess.V, highsuccess.V,qual.lo.V,qual.hi.V)
  colnames(replicates) = c("lowsuccess", "highsuccess","qual.lo","qual.hi")
  
  return(replicates)
} 

