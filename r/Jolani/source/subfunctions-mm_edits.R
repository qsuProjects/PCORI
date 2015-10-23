# Project: IPD
# date: 10 May 2014 only for binary predictor

# data generation (full data)
datagen <- function(setup = 0, stnum = "large",...){
  num <- switch(stnum, large = length(setup$studysize), small = 6)
  data <- matrix(NA,1,2+length(dimnames(setup$prop)[[2]]))
  colnames(data) <- c(dimnames(setup$prop)[[2]], "DV", "STID")
  for (i in 1:num){
    tempx <- rmvbin(setup$studysize[i],margprob=setup$prop[,,i],bincorr=setup$corr[,,i])
    tempbeta <- rmvnorm(1,setup$fixef,setup$tau)
    logit <- cbind(1,tempx)%*%t(tempbeta)
    tempdv <- rbinom(setup$studysize[i], 1, 1/(1 + exp(-logit))) #logit link
    data <- rbind(data, cbind(tempx, tempdv, i))
  }
  data <- data[-1,]
  n.class <- length(unique(data[,"STID"]))
  data[,"STID"] <- factor(data[,"STID"], labels=1:n.class)
  return(data)
}

# missing data generation
misgen <- function(pred = "one", perc = "low",...){
  mispred <- switch(pred, one = 1, two = 2)
  misperc <- switch(perc, low = .2, high = .5)
  data <- fulldata <- datagen(...) 
  tomis <- round(max(data[,"STID"])*misperc + 0.01)
  stid <- sample(1:max(data[,"STID"]), tomis)
  if (mispred == 2){
      loc <- round(length(stid)/2)
      data[data[,"STID"] %in% stid[1:loc],1] <- NA
      data[data[,"STID"] %in% stid[loc+1:length(stid)],2] <- NA
  }else{data[data[,"STID"] %in% stid,1] <- NA}
  list(fulldata = fulldata, misdata = data, sysmis = mispred, tstmis=tomis)
}

# pooling function
poolres <- function(fit, m){
  dim <- length(fixef(fit$analyses[[1]]))
  dimr <- length(attr(VarCorr(fit$analyses[[1]])$STID, "stddev"))
  qfix <- matrix(NA, m, dim)
  ufix <- matrix(0,  dim, dim)
  qran <- matrix(0,  dimr, dimr)
  for (i in 1:m){
    qfix[i,] <- fixef(fit$analyses[[i]])
    ufix <- ufix + vcov(fit$analyses[[i]])
    qran <- qran + VarCorr(fit$analyses[[i]])$STID[1:dimr,1:dimr]
  }
  Qbarfix <- apply(qfix, 2, mean)
  Ubarfix <- ufix/m
  Bfix <- cov(qfix)
  Tfix <- Ubarfix + Bfix*(1 + 1/m)
  resfix <- t(rbind(Qbarfix, sqrt(diag(Tfix)), Qbarfix/sqrt(diag(Tfix)), 
                    dnorm(abs(Qbarfix/sqrt(diag(Tfix))))))
  rownames(resfix) <- rownames(summary(fit$analyses[[1]])$coefficients)
  colnames(resfix) <- colnames(summary(fit$analyses[[1]])$coefficients)
  Qbarran <- qran/m
  riv <- (1 + 1/m)*diag(Bfix/Ubarfix)
  lambda <- (1 + 1/m)*diag(Bfix/Tfix)
  df <- mice:::mice.df(m, lambda, df.residual(fit$analyses[[1]]), method = "smallsample")
  fmi <- (riv + 2/(df + 3))/(riv + 1)
  list (fix.pool = resfix, ran.pool = Qbarran, riv = riv, fmi = fmi)
}

# Multilevel imputation for binary data
# standard prior for psi 
mice.impute.2l.bin <- function(y, ry, x,type){
  # the main code
  x <- data.frame(cbind(1, as.matrix(x)))
  names(x) <- paste("V",1:ncol(x),sep="")
  
  type <- c(2, type)
  
  clust <- names(x)[type==(-2)]
  rande <- names(x)[type==2]
  fixe <- names(x)[type>0]
  
  n.class <- length(unique(x[,clust]))
  x[,clust] <- factor(x[,clust], labels=1:n.class)
  lev<-levels(x[,clust])
  
  X<-as.matrix(x[,fixe])
  Z<-as.matrix(x[,rande])
  xobs<-x[ry,]
  yobs<-y[ry]
  Xobs<-as.matrix(X[ry,])
  Zobs<-as.matrix(Z[ry,])
  
  randmodel <- paste("yobs ~ ", paste(fixe[-1], collapse="+"), "+ ( 1 +", 
                     paste(rande[-1],collapse="+"), "|", clust, ")") # [-1] to remove intercept
  
  
  ################
  # MM added to catch case with only random intercept
  if (length(rande) < 2) {
    randmodel <- paste("yobs ~ ", paste(fixe[-1], collapse="+"),
              "+ ( 1 ", paste(rande[-1],collapse="+"), "|", clust, ")") # [-1] to remove intercept
  } 
  print(randmodel)
  ################
  
  
  suppressWarnings(fit <- try(glmer(formula(randmodel), data = data.frame(yobs,xobs), 
                                    family = binomial),silent=T))
  if(!is.null(attr(fit,"class"))){
    if(attr(fit,"class")=="try-error"){
      warning("glmer cannot be run, sorry!")
      return(y[!ry])
    }
  }  
  
  # draw beta*
  beta <- fixef(fit)
  rv <- t(chol(vcov(fit)))
  beta.star <- beta + rv %*% rnorm(ncol(rv))
  
  # calculate psi*
  rancoef <- as.matrix(ranef(fit)[[1]]) 
  lambda <- t(rancoef)%*%rancoef
  temp <- lambda
  temp <- ginv(temp)
  ev <- eigen(temp)
  if(sum(ev$values<0)>0)
  {
    ev$values[ev$values<0]<-0
    temp <- ev$vectors%*%diag(ev$values)%*%t(ev$vectors)
  }
  deco <- (ev$vectors)%*%sqrt(diag(ev$values))
  temp.psi.star <- rWishart(1, nrow(rancoef), diag(nrow(lambda)))[,,1]
  psi.star <- ginv(deco%*%temp.psi.star%*%t(deco)) 
  
  #### psi.star positive definite?
  if (!isSymmetric(psi.star)) psi.star <- (psi.star + t(psi.star))/2
  valprop<-eigen(psi.star)
  if(sum(valprop$values<0)>0)
  {
    valprop$values[valprop$values<0]<-0
    psi.star<-valprop$vectors%*%diag(valprop$values)%*%t(valprop$vectors)
  }
  
  # the main imputation task
  misindicator <- which((unique(x[,clust]) %in% unique(xobs[,clust])) == F)
  for (i in misindicator){
    suppressWarnings(bi.star <- t(rmvnorm(1,mean = rep(0,nrow(psi.star)), sigma = psi.star, meth="chol"))) # draw bi
    logit <- X[!ry & x[,clust]==i,] %*% beta.star + Z[!ry & x[,clust]==i,]%*% bi.star
    y[!ry & x[,clust]==i]<- rbinom(nrow(logit), 1, as.vector(1/(1 + exp(-logit))))
  }
  return(y[!ry])
}

#traditional multiple imputation 
tmi <- function(data, rm, MAX = 1, M = 1, ...){
  index <- which((apply(!is.na(data), 2, mean) == 1) == F)
  for (i in index) data[,i] <- as.factor(data[,i])
  data <- as.data.frame(data)
  imp0 <- mice(data, print = F, maxit = 0)
  pred <- imp0$pred
  pred[,"STID"] <- 0
  imp <- mice(data, print = F, pred = pred, meth = "logreg", maxit = MAX, m = M)
  suppressWarnings(repfit <- with(imp, glmer(formula(rm), family = binomial)))
  res <- poolres(repfit, m = imp$m)
  list(fix = res$fix.pool, rand = res$ran.pool, riv = res$riv, fmi = res$fmi)
}  

#multiple imputation with the fixed effect appraoch (random intercept) 
fixmi <- function(data, rm, MAX = 1, M = 1, ...){
  index <- which((apply(!is.na(data), 2, mean) == 1) == F)
  for (i in index) data[,i] <- as.factor(data[,i])
  data[,"STID"] <- as.factor(data[,"STID"])
  data <- as.data.frame(data)
  imp <- mice(data, print = F, meth = "logreg", maxit = MAX, m = M)
  suppressWarnings(repfit <- with(imp, glmer(formula(rm), family = binomial)))
  res <- poolres(repfit, m = imp$m)
  list(fix = res$fix.pool, rand = res$ran.pool, riv = res$riv, fmi = res$fmi)
} 

# multilevel multiple imputation
mlmi <- function(data, rm, MAX = 1, M = 1, ...){
  imp0 <- mice(data, print = F, maxit=0)
  pred <- imp0$pred  
  pred[pred == 1] <- 2
  pred[pred[,"STID"] != 0,"STID"] <- -2
  imp <- mice(data, print = F, pred=pred, method="2l.bin", maxit = MAX, m = M)
  suppressWarnings(repfit <- with(imp, glmer(formula(rm), family = binomial)))
  res <- poolres(repfit, m = imp$m)
  list(fix = res$fix.pool, rand = res$ran.pool, iter = MAX, impu = M, 
       riv = res$riv, fmi = res$fmi)  
}

# CCA: complete case analysis
fca <- function(data, rm = 1, fm = 1, ...){
  dim <- ncol(data) - 1
  fcadata <- na.omit(data)
  if (nrow(fcadata) > 0){
    if (!all(fcadata[,"STID"] == fcadata[1,"STID"])){
      suppressWarnings(fitfca <- glmer(formula(rm), data = fcadata, family = binomial))
      fixfca <- summary(fitfca)$coefficients
      ranfca <- VarCorr(fitfca)$STID[1:dim,1:dim]
    }else{
      fitfca <- glm(formula(fm), data = fcadata, family = binomial)
      fixfca <- summary(fitfca)$coefficients
      ranfca <- matrix(0,dim,dim)
    }
  }else 
    fixfca <- ranfca <- matrix(0,dim,dim)
  list(fix = fixfca, rand = ranfca)
}    

# analysis of the data
annals <- function(...){
  gen <- misgen(...)
  fdata <- as.data.frame(gen$fulldata)
  mdata <- as.data.frame(gen$misdata)
  randmodel <- paste("DV ~ ",paste(colnames(fdata)[1:(ncol(fdata)-2)],collapse="+"),
                             "+ ( 1 +", paste(colnames(fdata)[1:(ncol(fdata)-2)],collapse="+"), 
                             "|", "STID )")
  fixmodel <- paste("DV ~ ",paste(colnames(fdata)[1:(ncol(fdata)-2)],collapse="+"))
  # the number of parameters (this version fix and random are the same)
  dim <- ncol(fdata) - 1

  # full data (no missing)
  suppressWarnings(fitfull <- glmer(formula(randmodel), data = fdata, family = binomial))
  fixfull <- summary(fitfull)$coefficients
  ranfull <- VarCorr(fitfull)$STID[1:dim,1:dim]

  # CCA: complete case analysis (delition)
  fitfca <- fca(mdata, rm = randmodel, fm = fixmodel, ...)
  fixeffca <- fitfca$fix
  raneffca <- fitfca$rand
  
  # TMI: traditional multiple imputation (ignoring heterogeneity)
  fittmi <- tmi(mdata, randmodel, ...)
  fixeftmi <- fittmi$fix
  raneftmi <- fittmi$rand
  rivtmi <- fittmi$riv
  fmitmi <- fittmi$fmi
    
  # MIFIX: multiple imputation with fixed effect (random intercept only)
  fitfix <- fixmi(mdata, randmodel, ...)
  fixeffix <- fitfix$fix
  raneffix <- fitfix$rand
  rivfix <- fitfix$riv
  fmifix <- fitfix$fmi
    
  # MLMI: multilevel multiple imputation
  fitmlmi <- mlmi(mdata, randmodel, ...)
  fixefmlmi <- fitmlmi$fix
  ranefmlmi <- fitmlmi$rand
  rivmlmi <- fitmlmi$riv
  fmimlmi <- fitmlmi$fmi  
  timemlmi <- fitmlmi$time
  itermlmi <- fitmlmi$iter
  impumlmi <- fitmlmi$impu
  
  list(f.full = fixfull, f.fca = fixeffca, f.tmi = fixeftmi, f.mifix = fixefmifix, f.mlmi = fixefmlmi,
       r.full = ranfull, r.fca = raneffca, r.tmi = raneftmi, r.mifix = ranefmifix, r.mlmi = ranefmlmi,
       sysmis = gen$sysmis, tstmis = gen$tstmis ,time.mlmi = timemlmi, iter = itermlmi, impu = impumlmi,
       fdata = fdata, mdata=mdata)
}


