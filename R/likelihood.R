#Define log-likelihood function
#Using the Bernoulli quasi-log-likelihood function
#Logit model
ll_logit<-function(params,x,y){
  f <- NA
  for(i in 1:nrow(x)){
    f[i] <- y[i]*log(exp(x[i,]%*%params)/(1+exp(x[i,]%*%params))) + (1-y[i])*log(1-exp(x[i,]%*%params)/(1+exp(x[i,]%*%params)))
  }
  return(sum(f))
}
x <- cbind(rep(1,nrow(X)), X)
ll_logit(a$p, x, y)
#Probit model
ll_probit<-function(params,x,y){
  f <- NA
  for(i in 1:nrow(x)){
    f[i] <- y[i]*log(pnorm(x[i,]%*%params)) + (1-y[i])*log(1-pnorm(x[i,]%*%params))
  }
  return(sum(f))
}
#Loglog model
ll_loglog<-function(params,x,y){
  f <- NA
  for(i in 1:nrow(x)){
    f[i] <- y[i]*log(exp(-exp(x[i,]%*%params))) + (1-y[i])*log(1-exp(-exp(x[i,]%*%params)))
  }
  return(sum(f))
}
#Cloglog model
ll_cloglog<-function(params,x,y){
  f <- NA
  for(i in 1:nrow(x)){
    f[i] <- y[i]*log(1-exp(-exp(x[i,]%*%params))) + (1-y[i])*log(exp(-exp(x[i,]%*%params)))
  }
  return(sum(f))
}
#Cauchit model
ll_cauchit<-function(params,x,y){
  f <- NA
  for(i in 1:nrow(x)){
    f[i] <- y[i]*log((atan(x[i,]%*%params)+pi/2)/pi) + (1-y[i])*log(1-(atan(x[i,]%*%params)+pi/2)/pi)
  }
  return(sum(f))
}
