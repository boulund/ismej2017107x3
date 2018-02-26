library(igraph)
library(Matrix)
library(glmnet)
library(doParallel)
library(ggplot2)
library(gridExtra)
library(grid)
setwd("~/Dropbox/Shared/CMCI Projects/cmci_oxalate")

# generateCommunityAdjacency
# a function to return a random, weighted, directed, small world network
# nSpecies: defines the number of species in the community
# max.weight: defines the maximum weight for an edge; the minimum is -1 * maximum
# returns a nSpecies x nSpecies sparse matrix
generateCommunityAdjacency <- function(nSpecies,max.weight=10){
  #create a random, directed, small world network
  n500 <- as.directed(sample_smallworld(1,nSpecies,5,0.1,loops=TRUE), mode="arbitrary")
  n500 <- permute(n500, sample(1:nSpecies))
  
  E(n500)$weight <- runif(length(E(n500)),-1*max.weight,max.weight) # randomly assign weights to edges of the network
  
  #get the adjacency matrix for the network
  A0 <- as_adj(n500) 
  attr(A0,"x") <- E(n500)$weight
  
  A0
}


# generateSampleData
# a function to create samples of time series (t, t+1)
# community: an adjacency matrix defining interactions in the community
# nSample: the number of sample to generate for analysis (i.e. rows of data)
# meanPopSize: mean populations size (read counts) of the bacterial species present
# dispersion: defines the variance in population size where the var = meanPopSize + meanPopSize^2 / dispersion
# returns a list of with an x data matrix (time t) and a y data matrix (time t+1)
generateSampleData <- function(community, nSample, meanPopSize = 100, dispersion = 2){
  nSpp <- ncol(community)
  xData <- matrix(NA,nrow = nSample, ncol = nSpp)
  yData <- xData 
  
  #perform sampling
  for(i in 1:nSample){
    xData[i,]<-rnbinom(nSpp,size=dispersion,mu=meanPopSize) #generate a random state of the community at time t
    yData[i,] <- round(as.vector(xData[i,]*exp(community %*% (xData[i,]/sum(xData[i,]))))) #calculate state at t+1
  }
  
  list(X = xData, Y = yData)
  
}

#estimateCommunity
#a function to run the glmnet analyses on time series data
#timeSeries: the data to analyze, should be a list with elements X, Y which are both matrices
#returns a nSpecies x nSpecies sparseMatrix of estimated interactions
estimateCommunity<-function(timeSeries){
  X <- timeSeries$X
  Y <- timeSeries$Y
  nSpp <- ncol(Y)
  community <- sparseMatrix(NULL, NULL, x = 0, dims = c(nSpp,nSpp))
  for(i in 1:nSpp){
    model <- c()
    # the next line will fail if no non-zero parameters are estimated, the error will be thrown
    # but the loop should continue
    try(model <- cv.glmnet(X/rowSums(X),Y[,i],offset = log(X[,i]+1e-16), family = "poisson", intercept = F))
    if(is.null(model)) next
    if(model$nzero[which(model$lambda == model$lambda.min)] == 0) next
    parms <- t(t(summary(coef(model,s=model$lambda.min)))-c(1,0,0))
    community[i,parms[,1]]<-parms[,3]
  }
  community
}

#calculateSS
#a function to calculate the sensitivity and specificity 
#actual: true community matrix
#estimated: estimated community matrix
#threshold: quantile of small values to discard, must be between 0,1
#returns a names vector of the sensitivity, specificity, and MSE of the estimates
calculateSS <- function(actual, estimated, threshold = 0){
  if(threshold != 0 ){
    nzVals <- summary(estimated)[,3] #non-zero parameter estimates
    nzVals <- nzVals * (abs(nzVals) > quantile(abs(nzVals),threshold)) #make values less than threshold 0
    attr(estimated,"x") <- nzVals #pass the filtered values back to the estimated set
  }
  actual <- as.vector(actual)
  estimated <- as.vector(estimated) 
  res <- table( factor((actual != 0) + 2*(estimated != 0), levels = paste(0:3)) ) #0 = true neg, 1 = false neg, 2 = false pos, 3 = true pos
  out <- c(threshold, res["3"]/(res["1"]+res["3"]), res["0"]/(res["0"]+res["2"]), mean((actual-estimated)^2))
  names(out) <- c("threshold","sensitivity","specificity","MSE")
  out
}


###okay run across various values of community size and sampling effort
###
simVals <- data.frame(nSpp=sort(rep(c(50,250,500),4)),nSamp=rep(c(25,50,100,200),3))
subReplicates <- 100 #number of replicates per species + sample number combination
nReplicates <- nrow(simVals)*subReplicates #total number of replicates to run
registerDoParallel(cores=12)
simData <- foreach(j=0:(nReplicates-1),.combine = rbind, .errorhandling = 'remove') %dopar%{
  index <- j %/% 100 + 1
  test <- generateCommunityAdjacency(simVals[index,1])
  test.data <- generateSampleData(test,simVals[index,2])
  test.comm <- estimateCommunity(test.data)
  ROC <- matrix(NA,ncol=4,nrow=100)
  tVals <- seq(0,.99,length=100)
  for(i in 1:100) ROC[i,] <- calculateSS(test,test.comm,tVals[i])
  cbind(ROC,simVals[index,1], simVals[index,2])
}

save.image("Benchmarks.RData")

###saved

simData <- data.frame(simData)
names(simData) <- c("Threshold","Sensitivity","Specificity","MSE","CommunitySize","Samples")
head(simData)
simData$d.perf <- with(simData,((1-Sensitivity)^2 + (Specificity-1)^2)^0.5)

thePlots <- list()
index <- 1

for(i in unique(simVals$nSpp)){
  for(j in unique(simVals$nSamp)){
    p0 <- textGrob(bquote(paste(n[samples]==.(j))), gp = gpar(fontsize=20), rot = 90)
    
    p1 <- ggplot(subset(simData, CommunitySize == i & Samples == j), aes(x = 1 - Specificity, y = Sensitivity, col = Threshold))
    p1 <- p1 + geom_point(size = 0.5) + scale_colour_gradientn(colours = terrain.colors(10)) + theme_bw() + ylim(0,1)
    
    p2 <- ggplot(subset(simData, CommunitySize == i & Samples == j), aes(x = Threshold, y = MSE, color = d.perf))
    p2 <- p2 + geom_point(size = 0.5) + scale_color_gradient2(midpoint = 0.5, name = "Distance from\nPerfect") + theme_bw()
    
    p3 <- ggplot(subset(simData, CommunitySize == i & Samples == j & Threshold < 0.4), aes(x=Threshold, y=d.perf))
    p3 <- p3 + geom_smooth() + theme_bw() + ylab("Distance from Perfect")
    
    thePlots[[index]] <- arrangeGrob(p0,p1,p2,p3, ncol = 4, widths = c(1,10,10,10))
    index <- index + 1
    
  }
}

#I can't figure out how to get text size changed on "top" without using grid.text
#ggsave("benchmarks.pdf",marrangeGrob(thePlots, ncol=1, nrow=4, top = quote(paste("Community Size = ", unique(simVals$nSpp)[g]))), width=18)

#using grid.text I can't do it in one command
ggsave("benchmarks-1.pdf",marrangeGrob(thePlots[1:4], ncol=1, nrow=4, top = grid.text("Community Size = 50",  gp = gpar(fontsize = 20))), width=18)
ggsave("benchmarks-2.pdf",marrangeGrob(thePlots[5:8], ncol=1, nrow=4, top = grid.text("Community Size = 250",  gp = gpar(fontsize = 20))), width=18)
ggsave("benchmarks-3.pdf",marrangeGrob(thePlots[9:12], ncol=1, nrow=4, top = grid.text("Community Size = 500",  gp = gpar(fontsize = 20))), width=18)


with(subset(simData, CommunitySize == nSpp & Samples == nSamp),min(((1-Sensitivity)^2 + (Specificity-1)^2)^0.5))
with(subset(simData, CommunitySize == nSpp & Samples == nSamp),aggregate(((1-Sensitivity)^2 + (Specificity-1)^2)^0.5, data.frame(Threshold),mean))


