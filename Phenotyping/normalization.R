
#rescale<-function(raw,label){
#  model <- glm(label~raw, family="binomial")
#  return(model$fitted.values)
#}

## rescale intensity values to [0,1]
rescale<-function(raw,label){
  npos<-sum(label)
  ## identify the human-chosen threshold
  if(npos > 0){
    model <- glm(label~raw, family="binomial")
    #thres <- -model$coefficients[1]/model$coefficients[2]
    pratio<-npos/length(label)
    ## P(X = 1) = pratio (the empirical ratio) 
    thres <-(-log(1/pratio-1)-model$coefficients[1])/model$coefficients[2] 
    return(logist1(raw,thres))
  }
  if(npos == 0){
    thres <- max(raw)
    return(logist1(raw,thres))
  }
}

## scale cell intensity to Prob = (0,1). 
## if x < thres, the function return a value <
## idea is borrowed from the celesta paper
logist1<-function(x,thre,scale = 1){
  y = 1/(1+exp(-(x-thre)/scale))
  return(y)
}

vrescale<-Vectorize(rescale)

