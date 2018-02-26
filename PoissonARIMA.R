#Script to run the AR(1) model on the data
#on woodrat oxalate feeding trial
library(glmnet)

#set the directory and read in the data
setwd("~/Dropbox/Shared/CMCI Projects/cmci_oxalate")

#the file OxalateDataProcessing performs the task of 
#taking the raw data and creating acceptable X, Y
#matrices for the analysis. X has data from time t
#while Y has data from time t+1, thus allowing us 
#to perform the AR(1) analysis. Also provides oxy.taxa
#that contains taxanomic designations for the OTUs.
source("OxalateDataProcessing.R")

#findBetas(...)
#workhorse for performing glmnet and grid searching for 
#optimal alpha, lambda combination
#ARGUMENTS:
#X the matrix of data from time t
#Y the matrix of data from time t+1 (or to be regressed against X)
#index is the column index of the Y variable to be analyzed
#seq.alpha defines the alpha values to loop over
#nCV is the number of times to run cv.glmnet for a particular alpha values (to determine lambda)
#pSuccess is the percent of the time an alpha value needed to successfully return a model
#nlambda is as defined in cv.glmnet
#grouped is as defined in cv.glmnet
#... other flags to pass to glmnet
findBetas <- function (X, Y, index, seq.alpha = seq(0.5,1,0.1), nCV = 500, pSuccess = 0.2, nlambda = 100, grouped = TRUE, ...){
  require(doParallel)
  nAlpha <- length(seq.alpha) #number of alphas being tested
  maxIter <- nAlpha*nCV #perform nCV runs per alpha level, maxIter is the total number of loops needed
  
  #loop across alpha values to find the CVM which is either the MSE 
  #or deviance depending on the error family (see cv.glmnet)
  #use parallelization
  CVM <- foreach(i = 1:maxIter) %dopar% tryCatch(
    cv.glmnet(as.matrix(X),Y[,index], alpha=seq.alpha[(i-1) %/% nCV + 1], nlambda=nlambda, grouped = grouped, ...)$cvm, 
    error = function(e){
      print(paste("Index",i,"produced no output.")); 
      rep(NA,nlambda)}
  )
  
  #convert the CVM list into a matrix
  CVM.mat <- sapply(CVM,"[",1:nlambda)
  CVM.long <- CVM.mat[,1:nCV]
  for(i in 1:(nAlpha-1)) CVM.long <- rbind(CVM.long,CVM.mat[,i*nCV+1:nCV])
  
  #had to have at least pSuccess to keep data from that alpha value
  failures <- matrix(rowSums(is.na(CVM.long)) > (1-pSuccess)*nCV,nlambda,nAlpha)
  #get mean CVM values across lambda values (rows) then reshape
  mCVM <- matrix(rowMeans(CVM.long, na.rm=T),nlambda,nAlpha)  
  #get the sd of the CVM
  sdCVM <- matrix(apply(CVM.long,1,var, na.rm=T)^0.5,nlambda,nAlpha)
  mCVM[failures] <- NA
  sdCVM[failures] <- NA
  indices.min <- apply(mCVM+sdCVM,2,which.min) #find the index of minimum for each alpha
  
  if(length(indices.min) == 0) {print("No valid models were found."); return(data.frame())}
  
  #calculate the associated AIC values to decide which alpha, lambda combination to use
  findAIC <- function(glmnet.mod, nY,lambda){
    nP <- nrow(summary(coef(glmnet.mod,lambda)))
    theOffset <- eval(as.list(glmnet.mod$call)$offset)
    #this isn't thorough, only going to check for poisson family
    #otherwise use normal family
    prob <- ifelse(sum(class(glmnet.mod) == "fishnet"),
                   sum(dpois(Y[,nY],predict(glmnet.mod,as.matrix(X),s = lambda,type = "response",offset = theOffset),log=T)),
                   sum(dnorm(Y[,nY],predict(glmnet.mod,as.matrix(X),s = lambda,type = "response",offset = theOffset),log=T))
    )
    2*nP - 2*prob
  }
  
  AIC <- c()
  for(i in 1:nAlpha){
    if(1-length(indices.min[[i]])){AIC[i] <- Inf; next} #need this for cases where no index for the minimum exists
    cur.mod <- glmnet(as.matrix(X),Y[,index], alpha=seq.alpha[i], nlambda = nlambda, ...)
    AIC[i] <- findAIC(cur.mod,index,cur.mod$lambda[indices.min[[i]]])
  }
  
  if(min(AIC) == Inf) return(data.frame()) #return nothing if every alpha level failed
  
  mod.index <- which.min(AIC)
  final.model <- glmnet(as.matrix(X),Y[,index], alpha=seq.alpha[mod.index], nlambda = nlambda, ...)
  final.lambda <- final.model$lambda[indices.min[[mod.index]]]
  
  #return the coefficients of the chosen model
  out <- summary(coef(final.model,final.lambda))
  #as long as some coefficients were found, add
  #some additional data to the returned value for 
  #later analysis
  if(nrow(out) > 0){
    out$j <- index
    out$alpha <- seq.alpha[mod.index]
    out$lambda <- final.lambda
    out$dev.ratio <- final.model$dev.ratio[which(final.lambda == final.model$lambda)]
    out$AIC <- AIC[mod.index]
  }
  as.data.frame(out)
  
}


#####
#####
# Actually running the analysis
#####
#####

library(doParallel) 
registerDoParallel(cores=40)

#pf <- rep(1,ncol(X)-1) #do not force oxalate.consumed, last column X (Total.Reads) will not be used as a predictor, so substract 1 element
#Used for forcing oxalate
pf <- (names(X[,-1*ncol(X)]) != "Oxalate.Consumed")*1 #penalty factor vector that forces inclusion of oxalate.consumed

results <- c()

#some of the variables (e.g. oxalate consumed) are not count data
#so we don't want to use a Poisson error for those
#usePoisson tracks which Y variables we want to use a Poisson error for 
usePoisson <- !names(Y) %in% c("Oxalate.Consumed","Oxalate.Excreted","Oxalate.Degraded","Total.Reads") 

#loop across all Y values of interest
#for convenience the loop can run over whichever
#Y variable are of interest.
runIndices <- seq(1,ncol(Y))[ (! names(Y) %in% c("Total.Reads")) ]
for(i in runIndices){
  print(i)
  try(ifelse(usePoisson[i],
             results <- rbind(results, findBetas(X[,-632],Y,i, offset = log(X[,which(names(Y)[i]==names(X))]*X$Total.Reads + 1e-16) + log(X$Total.Reads), family = "poisson", grouped = FALSE, penalty = pf)),
             results <- rbind(results, findBetas(X[,-632],Y,i, grouped = FALSE, penalty = pf))
  ))
} 

write.csv(results,"PoissonARIMA-forcedOxalate.csv",row.names = FALSE) 