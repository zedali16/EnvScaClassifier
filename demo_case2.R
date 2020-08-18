setwd("D:/Dropbox/spectral envolope/EnvSccaClassifier_code")
source('functions.r')
library(astsa)
library(abind)
nrep = 100
n = c(100,200,500) #length of time series
rates_proposed = matrix(0,nrep,length(n))
rates_gamma = matrix(0,nrep,length(n))
rates_beta = matrix(0,nrep,length(n))
rates_seql = matrix(0,nrep,length(n))
nsub = 20 # number of time series per group in training data
ntest = 50
kappa = seq(0,1,0.1)
ind = matrix(0,nrow=nrep, ncol=length(n))
######################################################################
## Case 1
######################################################################
set.seed(20190106)
for(case in 1:length(n)){
  for(rep in 1:nrep){
    cat("repetition",rep,'\n')
    nobs = n[case]
    m = 4   # number of categories
    
    beta1 = cbind(c(0, 1.2, 1, 1),
                  c(0, 1, 0.8, 1),
                  c(0, 1, 1, 0.4))
    beta2 = cbind(c(0, 0.8, 1, 1),
                  c(0, 1, 0.4, 1),
                  c(0, 1, 1, 1.2))
    
    ##################################################################
    # simulate training set
    ##################################################################
    
    # Group 1
    yt1 = array(0,c(nobs,m-1,nsub))
    index1 = c()
    index1[1] = 1
    for(sub in 1:nsub){
      yt1[1,1,sub] = 1                            
      zt = matrix(0,nobs,m)
      for (i in 2:nobs){    
        zt[i-1,] = cbind(1, t(yt1[i-1,,sub]))
        deno = 1 + (exp(zt[i-1,]%*% beta1[,1]) + exp(zt[i-1,]%*%beta1[,2]) +
                      exp(zt[i-1,]%*% beta1[,3]))
        pi1 = exp(zt[i-1,]%*% beta1[,1])/deno
        pi2 = exp(zt[i-1,]%*% beta1[,2])/deno
        pi3 = exp(zt[i-1,]%*% beta1[,3])/deno
        pi4 = 1/deno
        
        index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
        if (index1[i]==4){
          yt1[i,1:(m-1),sub] = rep(0,m-1)
        } else {
          yt1[i,index1[i],sub]=1
        }
      }
    }
    
    # Group 2
    yt2 = array(0,c(nobs,m-1,nsub))
    index2 = c()
    index2[1] = 1
    for(sub in 1:nsub){
      yt2[1,1,sub] = 1                            
      zt = matrix(0,nobs,m)
      for (i in 2:nobs){    
        zt[i-1,] = cbind(1, t(yt2[i-1,,sub]))
        deno = 1 + (exp(zt[i-1,]%*% beta2[,1]) + exp(zt[i-1,]%*%beta2[,2]) +
                      exp(zt[i-1,]%*% beta2[,3]))
        pi1 = exp(zt[i-1,]%*% beta2[,1])/deno
        pi2 = exp(zt[i-1,]%*% beta2[,2])/deno
        pi3 = exp(zt[i-1,]%*% beta2[,3])/deno
        pi4 = 1/deno
        
        index2[i] = which(rmultinom(1,1, c(pi1,pi2,pi3,pi4))==1)
        if (index2[i]==4){
          yt2[i,1:(m-1),sub] = rep(0,m-1)
        } else {
          yt2[i,index2[i],sub]=1
        }
      }
    }
    
    
    ##################################################################
    # simulate test set
    ##################################################################
    
    # Group 1
    yt1_test = array(0,c(nobs,m-1,ntest))
    for(sub in 1:ntest){
      yt1_test[1,1,sub] = 1                            
      zt = matrix(0,nobs,m)
      for (i in 2:nobs){    
        zt[i-1,] = cbind(1, t(yt1_test[i-1,,sub]))
        deno = 1 + (exp(zt[i-1,]%*% beta1[,1]) + exp(zt[i-1,]%*%beta1[,2]) +
                      exp(zt[i-1,]%*% beta1[,3]))
        pi1 = exp(zt[i-1,]%*% beta1[,1])/deno
        pi2 = exp(zt[i-1,]%*% beta1[,2])/deno
        pi3 = exp(zt[i-1,]%*% beta1[,3])/deno
        pi4 = 1/deno
        
        index = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
        if (index==4){
          yt1_test[i,1:(m-1),sub] = rep(0,m-1)
        } else {
          yt1_test[i,index,sub]=1
        }
      }
    }
    
    # Group 2
    yt2_test = array(0,c(nobs,m-1,ntest))
    for(sub in 1:ntest){
      yt2_test[1,1,sub] = 1                            
      zt = matrix(0,nobs,m)
      for (i in 2:nobs){    
        zt[i-1,] = cbind(1, t(yt2_test[i-1,,sub]))
        deno = 1 + (exp(zt[i-1,]%*% beta2[,1]) + exp(zt[i-1,]%*%beta2[,2]) +
                      exp(zt[i-1,]%*% beta2[,3]))
        pi1 = exp(zt[i-1,]%*% beta2[,1])/deno
        pi2 = exp(zt[i-1,]%*% beta2[,2])/deno
        pi3 = exp(zt[i-1,]%*% beta2[,3])/deno
        pi4 = 1/deno
        
        index = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
        if (index==4){
          yt2_test[i,1:(m-1),sub] = rep(0,m-1)
        } else {
          yt2_test[i,index,sub]=1
        }
      }
    }
    
    ##################################################################
    # convert test and train datasets into format for sequence learner
    ##################################################################
    
    yt1train.seql <- character(length=nsub);
    yt2train.seql <- character(length=nsub);
    yt1test.seql <- character(length=ntest);
    yt2test.seql <- character(length=ntest);
    
    for (i in 1:nsub){
      tmp1 <- apply(yt1[,,i],1,function(x) ifelse(length(which(x==1))>0,which(x==1),0));
      yt1train.seql[i] <- gsub(", ","",toString(letters[tmp1+1],sep=""));
      
      tmp2 <- apply(yt2[,,i],1,function(x) ifelse(length(which(x==1))>0,which(x==1),0));
      yt2train.seql[i] <- gsub(", ","",toString(letters[tmp2+1],sep=""));
    }
    for (i in 1:ntest){
      tmp3 <- apply(yt1_test[,,i],1,function(x) ifelse(length(which(x==1))>0,which(x==1),0));
      yt1test.seql[i] <- gsub(", ","",toString(letters[tmp3+1],sep=""));
      
      tmp4 <- apply(yt2_test[,,i],1,function(x) ifelse(length(which(x==1))>0,which(x==1),0));
      yt2test.seql[i] <- gsub(", ","",toString(letters[tmp4+1],sep=""));
    }
    
    #create text file with proper format for sequence learner
    ytrain <- rbind(cbind("+1",yt1train.seql),cbind("-1",yt2train.seql));
    ytest <- rbind(cbind("+1",yt1test.seql),cbind("-1",yt2test.seql));
    write.table(ytrain,"seql executables/sim.char.train",row.names=FALSE,col.names=FALSE,sep=" ",quote=FALSE)
    write.table(ytest,"seql executables/sim.char.test",row.names=FALSE,col.names=FALSE,sep=" ",quote=FALSE)
    
    #run sequence learner
    seql.results <- shell("cd seql executables & seql_learn -n 1 -v 2 sim.char.train sim.seql.char.model & seql_mkmodel -i sim.seql.char.model -o sim.seql.char.model.bin -O sim.seql.char.model.predictors & seql_classify -n 1 -v 2 sim.char.test sim.seql.char.model.bin",intern=TRUE)
    
    #get results by search for accuracy
    tmp.acc <- seql.results[grep("Accuracy:",seql.results)]
    rates_seql[rep,case] <- eval(parse(text=sub("\\).*", "", sub(".*\\(", "", tmp.acc)))) 
    
    ##################################################################
    # run classification procedures
    ##################################################################
    yt = abind(yt1,yt2)
    yt_test = abind(yt1_test,yt2_test)
    group = c(rep(1,nsub), rep(2,nsub))
    test_group = c(rep(1,ntest), rep(2,ntest))
    cv = env_classifier_crossv(yt,group,n^(1/3),kappa)
    ind[rep,case] = min(which(cv==max(cv)))
    classes1 = env_classifier(yt,group,n^(1/3),yt_test,kappa[ind[rep,case]])
    rates_proposed[rep,case] = sum(classes1==test_group)/(2*ntest)
    classes2 = gamma_classifier(yt,group,n^(1/3),yt_test)
    rates_gamma[rep,case] = sum(classes2==test_group)/(2*ntest)
    classes3 = beta_classifier(yt,group,n^(1/3),yt_test)
    rates_beta[rep,case] = sum(classes3==test_group)/(2*ntest)
    
  }
}

apply(rates_proposed,2,mean)*100
apply(rates_proposed,2,sd)*100

apply(rates_beta,2,mean)*100
apply(rates_beta,2,sd)*100

apply(rates_gamma,2,mean)*100
apply(rates_gamma,2,sd)*100

apply(rates_seql,2,mean)*100
apply(rates_seql,2,sd)*100

