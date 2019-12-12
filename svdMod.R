library(svdComp5q0)

# from SVD-Comp.Rmd line 1040-1460

svdMod <- function(ql,Ql,N,S,offset,retAll,adult,q0Fix,smooth,C=4,printS=FALSE) {
  
  # ql is input qs
  # Ql is input summary indicators (child and adult ql)
  # N is number of samples
  # S is sample fraction
  # offset is the SVD offset
  # retAll is switch to return 'all' or just error summaries
  # adult is a switch indicating whether to include adult mortality in the model
  # q0Fix is a switch to execute the q0 fix
  # smooth is a switch to smooth left singular vectors
  # C is number of components, default C=4
  # printS is a switch to turn off printing the sample number at each iteration
  
  ret.ql.samp <- list(0) # sampled qls
  ret.ql.nsamp <- list(0) # out of sample qls
  ret.samp.names <- list(0) # sampled LT names
  ret.svd <- list(0) # svd of sampled qls
  ret.svd.sm <- list(0) # smooth svd of samlped qls
  ret.mods <- list(0) # models
  ret.recon.samp <- list(0) # sample reconstructions
  ret.error.samp <- list(0) # sample errors
  ret.recon.nsamp <- list(0) # out of sample reconstructions
  ret.error.nsamp <- list(0) # out of sample errors
  ret.errsum.samp <- list(0) # summary of sample errors
  ret.errsum.nsamp <- list(0) # summary of out of sample errors
  
  # Ensure C is in reasonable range: 1-4
  if (C < 1 | C > 4) {
    C = 4
    print("Setting C=4")
  }
  
  cat("\n")
  print(paste("Adult mortality is direct input to predictions:",adult))
  print(paste("SVD model is smoothed:",smooth))
  print(paste(N,"iterations"))
  print(paste(round(S*100,0),"% sample fraction",sep=""))
  print(paste(C,"components"))
  
  if (S > 0) {
    
    for (i in 1:N) {
      
      if (printS) {print(paste("  Sample:",i))}
      
      # pick the sample
      if (S == 1) {
        samp <- colnames(ql) # the sample is all LTs
        nsamp <- NA # nothing in the out of sample
      } else {
        # identify sample
        samp <- sample(colnames(ql),ncol(ql)*S) 
        # identify out of sample
        nsamp <- colnames(ql)[-which(colnames(ql) %in% samp)]	
      }
      name <- paste("s",i,sep="") # give this sample a name
      
      # store the sample
      ret.samp.names[[i]] <- samp # store sampled LT names to return list
      names(ret.samp.names)[i] <- name # name the sampled LT names
      
      # store the sampled qls
      ret.ql.samp[[i]] <- ql[,samp] # store sampled qls in return list
      names(ret.ql.samp)[i] <- name	 # name the sampled qls
      
      # store the out of sample qls
      if (S == 1) {
        ret.ql.nsamp[[i]] <- NA # no out of sample LTs
      } else {
        ret.ql.nsamp[[i]] <- ql[,nsamp] # store out of sample qls in return list
      }
      names(ret.ql.nsamp)[i] <- name # name out of sample qls
      
      # calculate the svd of the sampled qls
      svd <- svd(ql[,samp] - offset) # subtract offset before calculating svd
      
      # store the SVD
      ret.svd[[i]] <- svd # store svd in return list
      names(ret.svd)[i] <- name # name the svd
      
      # calculate transformations of *sampled* child mortality: input is logged
      cm <- expit(Ql[1,samp]) # child mortality, natural scale
      cml <- Ql[1,samp] # child mortality, logged
      cmls <- cml^2 # square of logged child mortality
      cmlc <- cml^3 # cube of logged child mortality
      
      # calcualte transformations of *sampled* adult mortality: input is logged
      am <- expit(Ql[2,samp]) # adult mortality, natural scale
      aml <- Ql[2,samp] # adult mortality, logged
      amls <- aml^2 # square of logged adult mortality
      amlc <- aml^3 # cube of logged adult mortality
      
      # calcualte one-way interaction of *sampled* child and adult mortality
      cmlaml <- cml*aml
      
      # model *sampled* adult mortality ~ child mortality
      aml.betas <- lm(aml ~ cm + cml + cmls + cmlc)
      
      # model *sampled* first four vs ~ child mortality and adult mortality
      v1.betas <- lm(svd$v[,1] ~ cm + cml + cmls + cmlc + am + amls + amlc + cmlaml)
      v2.betas <- lm(svd$v[,2] ~ cm + cml + cmls + cmlc + am + amls + amlc + cmlaml)
      v3.betas <- lm(svd$v[,3] ~ cm + cml + cmls + cmlc + am + amls + amlc + cmlaml)
      v4.betas <- lm(svd$v[,4] ~ cm + cml + cmls + cmlc + am + amls + amlc + cmlaml)
      
      # predictions for all LTs, both sampled and out of sample
      # start by transforming child mortality for all LTs
      cml.p <- Ql[1,]
      cm.p <- expit(cml.p)
      cmls.p <- cml.p^2
      cmlc.p <- cml.p^3
      
      # predict the adult mortality that goes with this child mortality
      # data frame of predictors
      predictors.aml <- data.frame(cbind(cm.p,cml.p,cmls.p,cmlc.p))
      # names for predictors that match the variable in the original model
      colnames(predictors.aml) <- c("cm","cml","cmls","cmlc")
      # predictions for adult mortality
      if (adult) {
        aml.p <- Ql[2,]
      } else {
        aml.p <- predict.lm(aml.betas,newdata=predictors.aml)
      }
      
      # predict vs using child mortality and (predicted) adult mortality
      # transform predicted adult mortality
      am.p <- expit(aml.p)
      amls.p <- aml.p^2
      amlc.p <- aml.p^3
      cmlaml.p <- cml.p*aml.p
      # data frame of predictors
      predictors.vs <- data.frame(cbind(
        cm.p,cml.p,cmls.p,cmlc.p,am.p,amls.p,amlc.p,cmlaml.p))
      # names for predictors that match the variables in the original regressions
      colnames(predictors.vs) <- c(
        "cm","cml","cmls","cmlc","am","amls","amlc","cmlaml")
      # predictions for each v
      v1.p <- predict.lm(v1.betas,newdata=predictors.vs)
      v2.p <- predict.lm(v2.betas,newdata=predictors.vs)
      v3.p <- predict.lm(v3.betas,newdata=predictors.vs)
      v4.p <- predict.lm(v4.betas,newdata=predictors.vs)
      
      # smooth left singular vectors
      if (smooth) {
        for (k in 2:6) {
          t <- ksmooth(seq(1,dim(svd$u)[1],1),svd$u[,k],kernel = "normal", bandwidth = (k+1))
          t$y[1:(k-1)] <- svd$u[,k][1:(k-1)]
          svd$u[,k] <- t$y
        }
        # store the smooth SVD
        ret.svd.sm[[i]] <- svd # store the smooth svd in return list
        names(ret.svd.sm)[i] <- name # name the smooth svd
      } else {
        ret.svd.sm[[i]] <- NA
        names(ret.svd.sm)[i] <- name # name the smooth svd
      }
      
      # reconstruct the predicted LTs
      v <- cbind(v1.p,v2.p,v3.p,v4.p)  	# matrix of new predicted vs
      r.p <- ql - ql # data frame for reconstructed values with names
      for (z in 1:C) { # loop over C components
        # svd reconstruction; sum of rank-1 matrices, one for each v
        r.p <- r.p + svd$d[z] * svd$u[,z] %*% t(v[,z]) 
      }
      r.p <- r.p + offset # add the offset back in
      
      if (q0Fix) {
        # fix up q0 prediction
        # child mortality for sample LTs
        cml <- Ql[1,samp] # sample child mortality
        cmls <- cml^2 # square of sample child mortality
        q0.betas <- lm(as.numeric(ql[1,samp]) ~ cml + cmls)	# q0 ~ cml + cmls
        # predictors for all LTs
        cml.p <- Ql[1,] # child mortality for all LTs
        cmls.p <- cml.p^2 # square of child mortality for all LTs
        predictors.q0 <- data.frame(cbind(cml.p,cmls.p))
        colnames(predictors.q0) <- c("cml","cmls")
        q0.p <- predict.lm(q0.betas,newdata=predictors.q0)
        # replace the predicted values for q0 with those from model above
        r.p[1,] <- q0.p 
      }	else {
        q0.betas <- NA
      }
      
      # store the models
      ret.mods[[i]] <- list(
        aml=aml.betas
        ,v1=v1.betas
        ,v2=v2.betas
        ,v3=v3.betas
        ,v4=v4.betas
        ,q0=q0.betas)
      names(ret.mods)[i] <- name
      
      # results: sampled LTs
      # store the reconstructed sampled LTs
      ret.recon.samp[[i]] <- r.p[,samp]
      names(ret.recon.samp)[i] <- name
      
      # store the errors in the reconstructed sampled LTs
      ret.error.samp[[i]] <- expit(ql[,samp]) - expit(r.p[,samp])
      names(ret.error.samp)[i] <- name
      
      # store summaries of errors in sampled LTs
      ret.errsum.samp[[i]] <- 
        summary(as.vector(as.matrix(ret.error.samp[[i]])))
      names(ret.errsum.samp)[i] <- name
      
      if (S == 1) {
        # results: out of sample LTs
        # store the reconstructed out of sample LTs
        ret.recon.nsamp[[i]] <- NA
        names(ret.recon.nsamp)[i] <- name
        
        # store the errors in the reconstructed out of sample LTs
        ret.error.nsamp[[i]] <- NA
        names(ret.error.nsamp)[i] <- name
        
        # store the summaries of errors in out of sample LTs
        ret.errsum.nsamp[[i]] <- NA
        names(ret.errsum.nsamp)[i] <- name
      } else {
        # results: out of sample LTs
        # store the reconstructed out of sample LTs
        ret.recon.nsamp[[i]] <- r.p[,nsamp]
        names(ret.recon.nsamp)[i] <- name
        
        # store the errors in the reconstructed out of sample LTs
        ret.error.nsamp[[i]] <- expit(ql[,nsamp]) - expit(r.p[,nsamp])
        names(ret.error.nsamp)[i] <- name
        
        # store the summaries of errors in out of sample LTs
        ret.errsum.nsamp[[i]] <- 
          summary(as.vector(as.matrix(ret.error.nsamp[[i]])))
        names(ret.errsum.nsamp)[i] <- name
      }
      
    }
    
    # put all the return lists together into one big list
    if (retAll == TRUE) {
      return(list(
        ql.samp = ret.ql.samp # sampled qls
        ,ql.nsamp = ret.ql.nsamp # out of sample qls
        ,names = ret.samp.names # sampled LT names
        ,svd = ret.svd # svd of sampled qls
        ,svd.sm = ret.svd.sm # smooth svd of sampled qls
        ,mods = ret.mods # models
        ,recon.samp = ret.recon.samp # sample reconstructions
        ,error.samp = ret.error.samp # sample errors
        ,recon.nsamp = ret.recon.nsamp # out of sample reconstructions
        ,error.nsamp = ret.error.nsamp # out of sample errors
        ,errsum.samp = ret.errsum.samp # summary of sample errors
        ,errsum.nsamp = ret.errsum.nsamp # summary of out of sample errors
        ,offset = offset # the offset necessary to reconstruct
        ,retAll = retAll
        ,adult = adult
        ,q0fix = q0Fix
        ,smooth = smooth
        ,C=C
      ))
    } else {
      return(list(
        errsum.samp = ret.errsum.samp # summary of sample errors
        ,errsum.nsamp = ret.errsum.nsamp # summary of out of sample errors
      ))
    }
  } else {
    print("S must be larger than 0.0")
    return()
  }
  
  print("Done")
  
}

ltPredict <- function(mods,smooth,cml,aml) {
  
  # mods is an output object from svdMod()
  # cml is a vector of logit child mortality rates
  # aml is a vector of logit adult mortality rates
  # if aml not supplied, then aml predicted from cml
  # smooth is a switch to use smoothed SVD if it's available
  
  cm <- expit(cml)
  cmls <- cml^2
  cmlc <- cml^3
  preds.aml <- data.frame(
    cm = as.numeric(cm),
    cml = as.numeric(cml),
    cmls = as.numeric(cmls),
    cmlc = as.numeric(cmlc)		
  )
  if(missing(aml)) {
    # predict adult mortality
    aml <- predict(mods$mods$s1$aml,newdata=preds.aml)
  }
  
  # predict vs
  am <- expit(aml)
  amls <- aml^2
  amlc <- aml^3
  cmlaml <- cml*aml
  preds.vs <- data.frame(
    cm = as.numeric(cm),
    cml = as.numeric(cml),
    cmls = as.numeric(cmls),
    cmlc = as.numeric(cmlc),
    am = as.numeric(am),
    amls = as.numeric(amls),
    amlc = as.numeric(amlc),
    cmlaml = as.numeric(cmlaml)	
  )
  v1 <- predict(mods$mods$s1$v1,newdata=preds.vs)
  v2 <- predict(mods$mods$s1$v2,newdata=preds.vs)
  v3 <- predict(mods$mods$s1$v3,newdata=preds.vs)
  v4 <- predict(mods$mods$s1$v4,newdata=preds.vs)
  
  # if smoothed SVD available
  if (smooth & mods$smooth) {
    svd <- mods$svd.sm$s1
  } else {
    svd <- mods$svd$s1
  }
  
  # construct LTs
  v <- cbind(v1,v2,v3,v4)
  r.p <- matrix(data=0,ncol=length(cml),nrow=length(svd$u[,1]))						
  for (z in 1:4) {		
    r.p <- r.p + svd$d[z] * svd$u[,z] %*% t(v[,z]) 	
  }
  r.p <- r.p + mods$offset 
  
  # fix q0
  if (mods$q0fix) {
    cmls <- cml^2
    preds.q0 <- data.frame(
      cml = as.numeric(cml),
      cmls = as.numeric(cmls)
    )
    r.p[1,] <- predict(mods$mods$s1$q0,newdata=preds.q0)
  }
  
  # fix up r.p
  r.p <- data.frame(r.p)
  colnames(r.p) <- paste("cml.",cml,sep="")
  rownames(r.p) <- rownames(mods$ql.samp$s1)
  
  # returns the matrix of predicted values
  return(r.p)
  
}

# do comparison 
# line 1615-1746

doComparison <- function(q1logit.f,Qlogit.f,q1logit.m,Qlogit.m
                         ,q5.f,q5.m,N,S,offset,adult,smooth,C) {
  
  # q1logit.f - female logit 1qx
  # Qlogit.f - female child and adult mortality levels
  # q1logit.m - male logit 1qx
  # Qlogit.m - male child and adult mortality levels
  # q5.f - female 5qx values from HMD site
  # q5.m - male 5qx values from HMD site
  # N - number of samples taken
  # S - sample fraction
  # offset - size of offset
  # adult - include adult mortality directly
  # smooth - use smoothing
  # C - number of components to use
  
  # rerun models setting using adult mortality directly
  mod.1_0.m <- svdMod(q1logit.m,Qlogit.m,N,S,10,TRUE,adult,TRUE,smooth,C)
  mod.1_0.f <- svdMod(q1logit.f,Qlogit.f,N,S,10,TRUE,adult,TRUE,smooth,C)
  
  # store the predicted values from the model in five-year age groups
  q5p.f <- convert1qxTo5qxApply(expit(mod.1_0.f$recon.samp$s1))
  q5p.m <- convert1qxTo5qxApply(expit(mod.1_0.m$recon.samp$s1))
  
  # create Log-Quad predictions
  # fit the log-quad using child mortality only
  
  # Source functions file
  source("../R/logQuad/DataProgramsExamples/R/functions.R")
  
  # Create labels for age vectors
  ages.5x1 <- c("0","1-4",paste(seq(5,105,5),seq(9,109,5),sep="-"),"110+")
  sexes <- c("Female","Male","Total")
  
  # Import matrix of model coefficients
  tmp1 <- read.csv("../R/logQuad/DataProgramsExamples/Data/coefs.logquad.HMD719.csv")
  tmp2 <- array(c(as.matrix(tmp1[, 3:6]))
                , dim=c(24, 3, 4)
                , dimnames=list(ages.5x1, sexes, c("ax", "bx", "cx", "vx")))
  coefs <- aperm(tmp2, c(1,3,2))
  
  # female
  q5.lq.f <- q5.f - q5.f
  e5.lq.f <- rbind(q5.f - q5.f,rep(0,ncol(q5.f)))
  a5.lq.f <- rbind(q5.f - q5.f,rep(0,ncol(q5.f)))
  l5.lq.f <- rbind(q5.f - q5.f,rep(0,ncol(q5.f)))
  for (i in 1:ncol(q5.f)) {
    if (adult) {
      lqfit <- lthat.any2.logquad(coefs,"Female",Q5=Q.f[1,i],QQa=Q.f[2,i]) # with adult
    } else {
      lqfit <- lthat.any2.logquad(coefs,"Female",Q5=Q.f[1,i],k=0) # without adult
    }
    q5.lq.f[,i] <- lqfit$lt[1:23,2]
    a5.lq.f[,i] <- lqfit$lt[1:24,3]
    e5.lq.f[,i] <- lqfit$lt[1:24,8]
    l5.lq.f[,i] <- lqfit$lt[1:24,4]
  }
  q5l.lq.f <- log(q5.lq.f)
  q5logit.lq.f <- logit(q5.lq.f)
  
  # male
  q5.lq.m <- q5.m - q5.m
  e5.lq.m <- rbind(q5.m - q5.m,rep(0,ncol(q5.m)))
  a5.lq.m <- rbind(q5.m - q5.m,rep(0,ncol(q5.m)))
  l5.lq.m <- rbind(q5.m - q5.m,rep(0,ncol(q5.m)))
  for (i in 1:ncol(q5.m)) {
    if (adult) {
      lqfit <- lthat.any2.logquad(coefs,"Male",Q5=Q.m[1,i],QQa=Q.m[2,i]) # with adult
    } else {
      lqfit <- lthat.any2.logquad(coefs,"Male",Q5=Q.m[1,i],k=0) # without adult
    }	
    q5.lq.m[,i] <- lqfit$lt[1:23,2]
    a5.lq.m[,i] <- lqfit$lt[1:24,3]
    e5.lq.m[,i] <- lqfit$lt[1:24,8]
    l5.lq.m[,i] <- lqfit$lt[1:24,4]
  }
  q5l.lq.m <- log(q5.lq.m)
  q5logit.lq.m <- logit(q5.lq.m)
  
  # compare the fits using the 5qx values obtained 
  #   directly from the HMD web site, q5.f and q5.m
  # construct a vector of comparison descriptors
  
  # females
  comps.f <- matrix(data = 0, nrow = 2, ncol = 3)
  colnames(comps.f) <- c("total.abs.error","mean.abs.error","max.error")
  rownames(comps.f) <- c("comp","lq")
  comps.f[1,1] <- sum(abs(q5.f - q5p.f))
  comps.f[2,1] <- sum(abs(q5.f - q5.lq.f))
  comps.f[,2] <- comps.f[,1]/(ncol(q5p.f)*nrow(q5p.f))
  comps.f[1,3] <- max(q5.f - q5p.f)
  comps.f[2,3] <- max(q5.f - q5.lq.f)
  
  # males
  comps.m <- matrix(data = 0, nrow = 2, ncol = 3)
  colnames(comps.m) <- c("total.abs.error","mean.abs.error","max.error")
  rownames(comps.m) <- c("comp","lq")
  comps.m[1,1] <- sum(abs(q5.m - q5p.m))
  comps.m[2,1] <- sum(abs(q5.m - q5.lq.m))
  comps.m[,2] <- comps.m[,1]/(ncol(q5p.m)*nrow(q5p.m))
  comps.m[1,3] <- max(q5.m - q5p.m)
  comps.m[2,3] <- max(q5.m - q5.lq.m)
  
  comps <- list(
    female = comps.f
    ,male = comps.m
    
    ,q5p.f = q5p.f
    ,q5p.m = q5p.m
    
    ,q5.lq.f = q5.lq.f
    ,e5.lq.f = e5.lq.f
    ,a5.lq.f = a5.lq.f
    ,l5.lq.f = l5.lq.f
    ,q5l.lq.f = q5l.lq.f
    ,q5logit.lq.f = q5logit.lq.f
    
    ,q5.lq.m = q5.lq.m
    ,e5.lq.m = e5.lq.m
    ,a5.lq.m = a5.lq.m
    ,l5.lq.m = l5.lq.m
    ,q5l.lq.m = q5l.lq.m
    ,q5logit.lq.m = q5logit.lq.m
  )
  
  return(comps)
  
}
