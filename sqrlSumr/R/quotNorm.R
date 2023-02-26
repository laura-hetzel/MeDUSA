# Implementation of quotient implementation 
# according to Dieterle et al., 2006; https://www.ncbi.nlm.nih.gov/pubmed/16808434
#
# Jan K
# last update: 2017-12-17
#
quotNorm = function(X, vars=1:dim(X)[2], NAerror=F, refsamples=NA) {
  # x:          data frame to be normalized
  # vars:       index vector of variables o be used, default: all
  # NAerrors:   throw error for NA's or just ignore?
  # refsamples: indices of samples to calculate reference sample on (e.g. only on control samples)
  
  # crash if there are any negative values
  if (any(unlist(X)[!is.na(unlist(X))]<0)) stop("Matrix contains negative values. Did you input logged data?")
  
  # check if there are any NAs
  if (sum(is.na(X[,vars]))>0) {
    # throw warning or error?
    if (NAerror) {
      stop('Data matrix contains NAs')
    } else {
      warning('Data matrix contains NAs')
    }
  }
  
  # if reference samples not given -> all samples
  if (is.na(refsamples)) refsamples <- 1:nrow(X)
  
  # median reference sample
  ref = apply(X[refsamples,vars],2,function(x)median(x,na.rm=T))
  # get dilution factors
  d = apply(X[,vars],1,  function(s) median(as.numeric(s/ref),na.rm=T))
  # apply to each sample  (for each row=sample, divide values by median dilution factor)
  Y = t(sapply(1:dim(X)[1], function(i)unlist(X[i,]/d[i])))
  
  #Y = t(apply(X,1,  function(s) s /  d) )
  rownames(Y) = rownames(X)
  
  # return
  list(X=Y,dilution=d)
}
