#Analysis of the results of the AR(1) model
#for the oxalate feeding trials.
library(xtable)
library(igraph)
library(gridExtra)
library(ggplot2)
library(glmnet)

setwd("~/Dropbox/Shared/CMCI Projects/cmci_oxalate")

source("OxalateDataProcessing.R") #provides useful matrices
#results <- read.csv("PoissonARIMA-unforcedOxalate.csv") 
results <- read.csv("PoissonARIMA-forcedOxalate.csv") 

indexOC <- which(c(1,names(X)) == "Oxalate.Consumed") #need the index of oxalate consumed repeatedly
#note that an extra element was added to the head of the vector to account for the intercept

#####
# Basic facts about the resulting models and some data cleaning
#####

goodDep <- c(1:624,626,627) #indices of the dependent parameters; 625 is oxalate consumed, avoid it
noModel <- ! goodDep %in% results$j #vector showing which y-variable for which model fitting failed
sum(noModel)
View(Y[,noModel])
#the number models fit
nMod<-sum(unique(results$j) %in% goodDep)
nMod * (ncol(X) - 1) #number of possible parameters (less the intercept)

nMod.good <- with(results, length(unique(j[j %in% goodDep & dev.ratio > 0.02]))) #number of models having dev.ratio > 0.02
nMod - nMod.good #number of low quality models

results <- results[results$i !=1, ]#drop intercept parameters
nrow(results) #number of interactions in the model
length(unique(results$j[results$i != indexOC])) #number of models that had more than the intercept and oxalate consumed

####Look at the unforced model
uf <- read.csv("PoissonARIMA-unforcedOxalate.csv")
uf <- uf[uf$i !=1, ]
nrow(uf) #number of non intercept parameters in the unforced model
length(unique(uf$j[uf$i != indexOC]))
uf <- uf[uf$dev.ratio > 0.02,]
uf.oxy <- uf$j[uf$i == indexOC] #models with oxalate consumed effects

results <- results[results$dev.ratio > 0.02,] #drop models will low deviance ratios
length(unique(results$j[results$i != indexOC])) #see how many "non-trivial" models are left after dropping the above
meanReads <- colMeans(X)
meanReads[1:6] <- 1 #animal effects need to have a 1
meanReads <- c(1,meanReads) #add a 1 for the intercept
results$x.scaled <- results$x * meanReads[results$i] #create scaled coefficients

#########
#Figure S4
#########
par(mfrow=c(2,1))
temp <- results[results$i != indexOC,c("i","j","x.scaled")]
temp$j <- temp$j + 7 # make indices i,j match
temp <- temp[temp$i > 7,] #drop the subject effects
names(temp)[3]<-"weight"
out.graph <- graph.data.frame(temp)
rm(temp)
in.dd <- hist(degree(out.graph,mode="in"), 40, right=F, xlab ="In-Degree", main=NULL, xlim=c(0,35), col="gray")
hist(degree(out.graph,mode="out"), 40, right=F, xlab ="Out-Degree", main=NULL,xlim=c(0,35),col = "gray")
#exported at 450 x 700 EPS - "Histograms-DD.eps"

#############################
#Calculation of Network Stats
#############################

#transitivity values for weight and unweighted
cc.weighted <- transitivity(out.graph,type="weighted")
mean(cc.weighted[! is.nan(cc.weighted) & ! is.infinite(cc.weighted)]) #mean weighted
transitivity(out.graph) #global clustering
mean(transitivity(out.graph,type="local"),na.rm=T) #mean unweighted

#same thing but for the random graph
rg <- erdos.renyi.game(length(V(out.graph)),length(E(out.graph)), type="gnm", directed=T)
rg <- set.edge.attribute(rg,"weight",value=get.edge.attribute(out.graph,"weight"))
rg.weighted <- transitivity(rg,type="weighted")
mean(rg.weighted, na.rm=T)
transitivity(rg)
mean(transitivity(rg,type="local"), na.rm=T)

#look at path lengths
mean_distance(out.graph) #unweighted
mean_distance(rg)

out.dMat <- distances(out.graph, weights = abs(get.edge.attribute(out.graph,"weight")), mode = "out")
diag(out.dMat) <- NA
out.dMat[is.infinite(out.dMat)] <- NA
mean(out.dMat, na.rm = T)

rg.dMat <- distances(rg,weights = abs(get.edge.attribute(rg,"weight")), mode = "out")
diag(rg.dMat) <- NA
rg.dMat[is.infinite(rg.dMat)] <- NA
mean(rg.dMat, na.rm = T)
rm(rg.dMat, rg, rg.weighted, out.dMat)
### end of network statistics

#########
#Figure 2
#########

par(mfrow=c(1,2))

hist(unique(results$dev.ratio),30,main=NULL,xlab = expression(paste("pseudo-",italic(R)^2,sep=""))) #pretty even spread
alphaBetter <- unique(results[,c("alpha","dev.ratio")])
plot(alphaBetter,xlab=expression(alpha),ylab=expression(paste("pseudo-",italic(R)^2,sep="")))
lines(supsmu(alphaBetter$alpha,alphaBetter$dev.ratio),col="red",lty=2,lwd=2)
 
#do a little non-linear model
anova(lm(dev.ratio ~ poly(alpha,3), alphaBetter),lm(dev.ratio ~ poly(alpha,4), alphaBetter))#stop at 3-rd degree
theModel <-lm(dev.ratio ~ poly(alpha,3), alphaBetter)
lines(seq(0.5,1,length=30),predict(theModel,list(alpha=seq(0.5,1,length=30))),col="blue",lwd=2)
summary(theModel)

rm(alphaBetter)
#exported at 700 x 350 EPS - "Deviance.eps"

###############
#Table 1 values
###############
Oxy.Taxa <- t(mapply(as.character,oxy.taxa))
oxyCoefs <- results[results$i == indexOC,c("j","x.scaled")]
good.j <- oxyCoefs$j %in% grep("denovo",names(Y))# indices of coefficients that correspond to OTUs
oxyCoefs <- oxyCoefs[good.j,] #drop unwanted coefficients
oxyCoefs$theTaxa <- names(Y)[oxyCoefs$j]
oxyCoefs <- oxyCoefs[order(oxyCoefs$x.scaled),]
oxyCoefs$Taxa <- Oxy.Taxa[oxyCoefs$theTaxa,1]
oxyCoefs$Division <- Oxy.Taxa[oxyCoefs$theTaxa,2]

TaxaEffects <- data.frame(taxaNeg = oxyCoefs$Taxa[1:15], negEff = oxyCoefs$x.scaled[1:15], taxaPos = rev(oxyCoefs$Taxa)[1:15], posEff = rev(oxyCoefs$x.scaled)[1:15])

# get the table for LaTeX
tab.TaxaEffects<-xtable(TaxaEffects)
print(tab.TaxaEffects,include.rownames=F)
# write csv for Excel tables
#write.csv(TaxaEffects,"OxyTaxaEffects-v2.csv",row.names = F)

#here are the taxa that had oxalate consumed effects in the unforced model
good.uf.oxy <-uf.oxy[grep("denovo",names(Y)[uf.oxy])]
View(data.frame(Taxa = Oxy.Taxa[names(Y)[good.uf.oxy],1], Effect=uf$x[uf$j %in% good.uf.oxy & uf$i == indexOC]))

###############
#Table 2 values
###############

meanEffects <- aggregate(data.frame(N=1,x=oxyCoefs$x.scaled),by=oxyCoefs[,c("Taxa","Division")],sum)
meanEffects$x <- meanEffects$x/meanEffects$N
meanEffects <- meanEffects[order(meanEffects$x),]
names(meanEffects)[4] <- "Mean" 

#write files for LaTeX and Excel
print(xtable(cbind(head(meanEffects,10),tail(meanEffects,10)[10:1,])),include.rownames=F)
#write.csv(cbind(head(meanEffects,10),tail(meanEffects,10)[10:1,]),"MeanOxyEffects-v2.csv",row.names = F)



#########
#Figure 4
#########

par(mfrow=c(3,1))
h1 <- hist(oxyCoefs$x.scaled, 100, main = "All OTUs", xlim =c(-5,5), xlab = expression(hat(beta)[C[2]*O[4]^{"2-"}]))
hist(oxyCoefs$x.scaled[oxyCoefs$Taxa == "Ruminococcaceae"],breaks=h1$breaks, xlim =c(-5,5), xlab = expression(hat(beta)[C[2]*O[4]^{"2-"}]),main=expression(paste(bold("Ruminococcaceae OTUs"))))
abline(v=mean(oxyCoefs$x.scaled[oxyCoefs$Taxa == "Ruminococcaceae"]),col="red",lty=2,lwd=2)
hist(meanEffects$Mean, xlim =c(-5,5), breaks=h1$breaks,xlab = expression(E*bgroup("[",hat(beta)[C[2]*O[4]^{"2-"}],"]")),main="Mean Effects by Taxonomic ID")
#exported at 450 x 700 EPS - "OxyEffects.eps"


#########
#Figure 3
#########

par(mfrow=c(1,1))
timepoint <- c(2:4,2:6,2:4,2:6,3:4,2:6) #these are the actual timepoints for the rows of the Y data
subject <- LETTERS[as.matrix(X[,1:6]) %*% 1:6] #give subjects a letter designation
tVals <- c("t[0]","t[1]","t[2]","t[3]","t[4]","t[5]")
xlabs <- c(bquote(italic(t)[1]),bquote(italic(t)[2]),bquote(italic(t)[3]),bquote(italic(t)[4]),bquote(italic(t)[5])) #pretty labels for plotting

#set the penalty factor appropriately for the model that was used
pf <- rep(1,ncol(X)-1)
#pf <- (names(X[,-1*ncol(X)]) != "Oxalate.Consumed")*1

#Grab the most abundant species within the "good" fit category (r^2 > 0.8)
mostAbund <- which.max(colSums(Y[,unique(results[results$dev.ratio > 0.8 ,"j"])]))
testIndex <- unique(results[results$dev.ratio > 0.8 ,"j"])[mostAbund]
with(results, unique(dev.ratio[j == testIndex])) #see the dev.ratio

test <- glmnet(as.matrix(X[,-632]), Y[,testIndex], alpha = results[results$j==testIndex,"alpha"][1], family = "poisson", offset = log(X[,names(Y)[testIndex]]*X$Total.Reads + 1e-16), penalty.factor = pf)

preds <- predict(test, as.matrix(X[,-632]), offset = log(X[,names(Y)[testIndex]]*X$Total.Reads + 1e-16), s = results[results$j==testIndex,"lambda"][1], type = "response") 

actual <- Y[,testIndex]
predicted <- preds

otuID <- Oxy.Taxa[names(Y)[testIndex],1]
otuID <- paste(otuID," (",names(otuID),")", sep="")
df.plot <- data.frame(timepoint=tVals[timepoint],subject,actual,predicted=unname(predicted))

g1 <- ggplot(data=df.plot, aes(x=timepoint,y=actual, color = subject))
g1 <- g1 + geom_point(size = 3) + scale_y_log10(limits=c(10,10000), breaks=10^(1:4)) + theme_bw() + theme(legend.position="none") + geom_line(aes(x = as.numeric(timepoint), y = predicted), size = 1) + xlab("Observation") + ylab(bquote("Log"[10]*" Reads")) + ggtitle(bquote(.(otuID))) + scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'), name = "Subject") + scale_x_discrete(labels = xlabs)
g1

#just take Oxalobacter for the mid-level fit because we are interested in it.
#obacter <- rownames(Oxy.Taxa)[Oxy.Taxa[,1] == "Oxalobacteraceae"]
#testIndex  <- which(names(Y) == obacter)

mostAbund <- which.max(colSums(Y[,unique(results[results$dev.ratio < 0.8 & results$dev.ratio > 0.5 ,"j"])]))
testIndex <- unique(results[results$dev.ratio < 0.8 & results$dev.ratio > 0.5,"j"])[mostAbund]
with(results, unique(dev.ratio[j == testIndex])) #see the dev.ratio

test <- glmnet(as.matrix(X[,-632]), Y[,testIndex], alpha = results[results$j==testIndex,"alpha"][1], family = "poisson", offset = log(X[,names(Y)[testIndex]]*X$Total.Reads + 1e-16), penalty.factor = pf)

preds <- predict(test, as.matrix(X[,-632]), offset = log(X[,names(Y)[testIndex]]*X$Total.Reads + 1e-16), s = results[results$j==testIndex,"lambda"][1], type = "response") 

actual <- Y[,testIndex]
predicted <- preds

otuID <- Oxy.Taxa[names(Y)[testIndex],1]
otuID <- paste(otuID," (",names(otuID),")", sep="")
df.plot <- data.frame(timepoint=tVals[timepoint],subject,actual,predicted=unname(predicted))

g2 <- ggplot(data=df.plot, aes(x=timepoint,y=actual, color = subject))
g2 <- g2 + geom_point(size = 3) + scale_y_log10(limits=c(10,10000), breaks=10^(1:4)) + theme_bw() + theme(legend.position="none") + geom_line(aes(x = as.numeric(timepoint), y = predicted), size = 1) + xlab("Observation") + ylab(bquote("Log"[10]*" Reads")) + ggtitle(bquote(.(otuID))) + scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'), name = "Subject") + scale_x_discrete(labels = xlabs )
g2

#### grab a low quality one
mostAbund <- which.max(colSums(Y[,unique(results[results$dev.ratio < 0.5 &  results$dev.ratio > 0.2,"j"])]))
testIndex <- unique(results[results$dev.ratio < 0.5 &  results$dev.ratio > 0.2,"j"])[mostAbund]
with(results, unique(dev.ratio[j == testIndex])) #see the dev.ratio

test <- glmnet(as.matrix(X[,-632]), Y[,testIndex], alpha = results[results$j==testIndex,"alpha"][1], family = "poisson", offset = log(X[,names(Y)[testIndex]]*X$Total.Reads + 1e-16), penalty.factor = pf)

preds <- predict(test, as.matrix(X[,-632]), offset = log(X[,names(Y)[testIndex]]*X$Total.Reads + 1e-16), s = results[results$j==testIndex,"lambda"][1], type = "response") 

actual <- Y[,testIndex]
predicted <- preds

otuID <- Oxy.Taxa[names(Y)[testIndex],1]
otuID <- paste(otuID," (",names(otuID),")", sep="")
df.plot <- data.frame(timepoint=tVals[timepoint],subject,actual,predicted=unname(predicted))

g3 <- ggplot(data=df.plot, aes(x=timepoint,y=actual, color = subject))
g3 <- g3 + geom_point(size = 3) + scale_y_log10(limits=c(10,10000), breaks=10^(1:4))  + theme_bw() + theme(legend.position="none") + geom_line(aes(x = as.numeric(timepoint), y = predicted), size = 1) + xlab("Observation") + ylab(bquote("Log"[10]*" Reads")) + ggtitle(bquote(.(otuID))) + scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'), name = "Subject") + scale_x_discrete(labels = xlabs)
g3

grid.arrange(g1,g2,g3,nrow=3)
#Exported at 450 x 700 EPS, "Dynamics.eps"
