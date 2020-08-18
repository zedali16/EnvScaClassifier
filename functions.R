#==============================================================
# Convert a categorical ts into multivariate indicator process
#==============================================================
cat_convert <- function(xt){
  stage = sort(unique(xt))
  nobs = length(xt)
  yt = matrix(0,nobs,length(stage)-1)
  for (j in 1:length(stage)-1){
    yt[,j] = (xt==stage[j])*1
  }
  return(yt)
}
#==============================================
# Given a m times p-1 categorical time series,
# compute: 
# (1) the Fourier frequencies
# (2) the spectral envelope 
# (3) the optimal scaling 
#=============================================
env.get <- function(yt,L){
  dimen = dim(yt)[2]
  v = var(yt)
  fyy = mvspec(yt, spans = c(L,L), plot=FALSE,kernel="fejer") # spectral
  fyy_re = Re(fyy$fxx)
  vv = eigen(v)
  
  Q = diag(dimen)
  num   = fyy$n.used	
  nfreq = length(fyy$freq)
  specenv = matrix(0, nfreq,1)                                 
  beta = matrix(0, nfreq, dimen)
  for (k in 1:nfreq){
    ev = eigen(2*Q%*%fyy_re[,,k]%*%Q/num, symmetric = TRUE)
    specenv[k] = ev$values[1]
    b = Q%*%ev$vectors[,1]
    beta[k,] = b/sqrt(sum(b^2))
  }
  freq = fyy$freq
  output = list(freq = freq, envelope = specenv, scale = beta)
  return(output)
}
#=======================================================
# Given a several time series within a group
# Compute 
# (1) the group level spectral envelope
# (2) the group level optimal scaling
#=======================================================
group_env <- function(yt_group,L){
  nsub = dim(yt_group)[3]
  beta = 0
  specenv = 0
  dist = 0
  for (k in 1:nsub){
    output = env.get(yt_group[,,k], L)
    tmp_env = output$envelope
    tmp_scl = output$scale
    tmp_dist = output$dist
    specenv = specenv + tmp_env/nsub
    beta = beta + tmp_scl/nsub
  }
  freq = output$freq
  output = list(freq = freq, envelope = specenv, scale = beta)
  return(output)
}
#========================================================
# The main function
# assign class to a time series
#========================================================
env_classifier <- function(yt, group, L, yt_new, kappa){
  nnew = dim(yt_new)[3]
  if(is.na(nnew)){
    nnew=1
  }else{
    nnew = nnew
  }
  nclass = length(unique(group))
  classes = c()
  env = list()
  scal = list()
  # calculate group level statistics based on training time series
  for (j in 1:nclass){
    env[[j]] = group_env(yt[,,group==j],L)$envelope
    scal[[j]] = group_env(yt[,,group==j],L)$scale
  }
  # for each of testing time series, assign a group to it
  for (k in 1:nnew){
    if(nnew==1){
      new_env = env.get(yt_new,L)$envelope
      new_scal = env.get(yt_new,L)$scale
    }else{
      new_env = env.get(yt_new[,,k],L)$envelope
      new_scal = env.get(yt_new[,,k],L)$scale
    }
    g1 =  new_env
    g2 =  new_scal
    g = c()
    for (j in 1:nclass){
      temp1 = env[[j]]
      temp2 = scal[[j]]
      g[j] = kappa*sum((new_env-temp1)^2)/sum(new_env^2) + 
              (1-kappa)*sum((new_scal-temp2)^2)/sum(new_scal^2)
      }
      classes[k] = which(g==min(g))
    }
  return(classes)
}
#======================================================
# Use cross-validation to do grid search of kappa
#======================================================
env_classifier_crossv <- function(yt, group, L, kappa){
  nclass = length(unique(group))
  ntun = length(kappa)
  classes = rep(0,length(kappa))
  for (jj in 1:ntun){
    for (k in 1:length(group)){
      yt_temp = yt[,,-k]
      group_temp = group[-k]
      yt_test = yt[,,k]
      group_test = group[k]
      env = list()
      scal = list()
      for (j in 1:nclass){
        output1 =  group_env(yt_temp[,,group_temp==j],L)
        env[[j]] = output1$envelope
        scal[[j]] = output1$scale
      }
      output2 = env.get(yt_test,L)
      new_env = output2$envelope
      new_scal = output2$scale
      g = c()
      for (j in 1:nclass){
        temp1 = env[[j]]
        temp2 = scal[[j]]
        g[j] = kappa[jj]*sum((new_env-temp1)^2)/sum(new_env^2) + 
              (1-kappa[jj])*sum((new_scal-temp2)^2)/sum(new_scal^2)
      }
      ig = which(g==min(g))
      classes[jj] = classes[jj] + (group_test==ig)
    }
  }
  return(classes)
}    
#========================================================
# Classification use ONLY optimal scaling
#========================================================
beta_classifier <- function(yt, group, L, yt_new){
  nnew = dim(yt_new)[3]
  if(is.na(nnew)){
    nnew=1
  }else{
    nnew = nnew
  }
  nclass = length(unique(group))
  classes = c()
  dist = list()
  for (j in 1:nclass){
    dist[[j]] = group_env(yt[,,group==j],L)$scale
  }
  for (k in 1:nnew){
    if(nnew==1){
      new_dist = env.get(yt_new,L)$scale
    }else{
      new_dist = env.get(yt_new[,,k],L)$scale
    }
    g = c()
    for (j in 1:nclass){
      temp = dist[[j]]
      g[j] = sum((new_dist-temp)^2)	
    }
    classes[k] = which(g==min(g))
  }
  return(classes)
}
#========================================================
# Classification use ONLY spectral envelope
#========================================================
gamma_classifier <- function(yt, group, L, yt_new){
  nnew = dim(yt_new)[3]
  if(is.na(nnew)){
    nnew=1
  }else{
    nnew = nnew
  }
  nclass = length(unique(group))
  classes = c()
  dist = list()
  for (j in 1:nclass){
    dist[[j]] = group_env(yt[,,group==j],L)$envelope
  }
  for (k in 1:nnew){
    if(nnew==1){
      new_dist = env.get(yt_new,L)$envelope
    }else{
        new_dist = env.get(yt_new[,,k],L)$envelope
    }
    g = c()
    for (j in 1:nclass){
      temp = dist[[j]]
      g[j] = sum((new_dist-temp)^2)	
    }
    classes[k] = which(g==min(g))
  }
  return(classes)
}