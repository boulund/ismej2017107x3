library(ggplot2)
setwd("C://Users/JT/Dropbox/cmci_oxalate/")
source("OxalateDataProcessing.R") #provides useful matrices
results <- read.csv("PoissonARIMA-forcedOxalate.csv") 
indexOC <- which(c(1,names(X)) == "Oxalate.Consumed") #need the index of oxalate consumed repeatedly
#note that an extra element was added to the head of the vector to account for the intercept
results <- results[results$dev.ratio > 0.02,] #drop models will low deviance ratios
results <- results[results$i !=1, ]#drop intercept parameters
meanReads <- colMeans(X)
meanReads[1:6] <- 1 #animal effects need to have a 1
meanReads <- c(1,meanReads) #add a 1 for the intercept
results$x.scaled <- results$x * meanReads[results$i] #create scaled coefficients

#########
#Figure 6
#########

#generating Figure 6 takes some time because it find the taxonomy information from greengenes

library(reshape2)
library(stringr)
#cut.off <- 0.0
wide <- dcast(results, j ~ i, value.var = "x.scaled")
names(wide) <- c("Y.OTU",names(X)[as.numeric(names(wide)[-1])-1]) #get OTU names from "Y"
wide$Y.OTU <- names(Y)[wide$Y.OTU]
out.hist <- results
out.hist$i.denovo <- as.factor(c("intercept", names(X))[out.hist$i])
out.hist$j.denovo <- as.factor(names(Y)[out.hist$j])
out.hist$j.taxa <- as.character(out.hist$j.denovo)
out.hist$i.taxa <- as.character(out.hist$i.denovo)
#fill in some taxa names
out.hist$j.taxa <- sapply(oxy.taxa[1,],as.character)[out.hist$j.taxa]
out.hist$i.taxa <- sapply(oxy.taxa[1,],as.character)[out.hist$i.taxa]
out.hist$j.taxa[is.na(out.hist$j.taxa)] <- as.character(out.hist$j.denovo[is.na(out.hist$j.taxa)])
out.hist$i.taxa[is.na(out.hist$i.taxa)] <- as.character(out.hist$i.denovo[is.na(out.hist$i.taxa)])
out.hist$i.taxa <- as.character(str_trim(out.hist$i.taxa))
out.hist$j.taxa <- as.character(str_trim(out.hist$j.taxa))
out.hist <- out.hist[!out.hist$i.taxa %in% c("Bacteria"), ] #get rid of the high level groups (for now just Bacteria?)

#print off basic heatmap
#if you have subject effects, might want to split into two heatmaps because x.scaled values are on based on read counts for subjects (or treatments for that matter)
heat <- aggregate(x.scaled ~ j.taxa + i.taxa, data=out.hist, "mean")
#simple heatmap
ggplot(data = heat, aes(x=i.taxa,y=j.taxa,fill=x.scaled)) + geom_tile(colour="white") + theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5, size=5),axis.text.y=element_text(size=5)) + scale_fill_gradient2(na.value="grey", name="Interaction Strength",low="blue",mid="white",high="red")

#add taxonomy table
library(phytools)
library(phyloseq)
library(data.table)
library(ggtree)
library(ggdendro)
library(grid)
library(gridExtra)
#need to download and load greengenes phylogeny from
#http://greengenes.secondgenome.com/downloads/database/13_5
grntree <- read_tree_greengenes("gg_13_5_otus_99_annotated.tree")
grntax <- read.delim(file = "gg_13_5_taxonomy.txt", header=FALSE, stringsAsFactors=FALSE)
grntax[,2] <- gsub(";|\\[|\\]| ","",grntax[,2]) #get rid of bad characters
grntax[,2] <- gsub("Enterobacteriaceae","Enterobacteraceae",grntax[,2]) #fix difference in spelling between dbcamplicons and greengenes.
df <- matrix(,nrow=nrow(grntax),ncol=8)
#parse out individual names from concatendated names. Takes a few minuites
for(i in 1:nrow(grntax)){
  df[i,1:length(df2)] <- df2 <- t(as.matrix(unlist(strsplit(grntax[i,2],".__"))))
  if( !is.na(df[i,8])){
    df[i,8] <- paste(df[i,7],df[i,8],sep=" ")
  }
}
df<-df[,-1]
grntax <- cbind(grntax[,1],grntax[,2],df)
colnames(grntax) <- c("idnum","string","k","p","c","o","f","g","s")
grntax[grntax==""] <- NA   #fill in empty cells with NAs
grntax[,1] <- as.numeric(as.character(grntax[,1]))
tax <- as.matrix(unique(out.hist$j.taxa))
colnames(tax) <- "name"
species <- merge(grntax,tax,by.x="s",by.y="name") #just keep tips that showed up in our analysis (by name of course)
species <- subset(species, select = c(idnum,k,p,c,o,f,g,s))
species$idnum <- as.numeric(as.character(species$idnum))
species <- species[sample(nrow(species)),] #rearrange the order of the taxa so they aren't all in blocks
short_s <- aggregate(idnum ~ s, FUN = max, data=species) #condense the species list. grab a greengenes # at random
genre <- merge(grntax,tax,by.x="g",by.y="name")
genre <- genre[which(is.na(genre$s)),]
genre <- subset(genre, select = c(idnum,k,p,c,o,f,g,s))
genre$idnum <- as.numeric(as.character(genre$idnum))
genre <- genre[sample(nrow(genre)),]
short_g <- aggregate(idnum ~ g, FUN = max, data=genre)
family <- merge(grntax,tax,by.x="f",by.y="name")            
family <- family[which(is.na(family$g)),]
family <- subset(family, select = c(idnum,k,p,c,o,f,g,s))
family$idnum <- as.numeric(as.character(family$idnum))
family <- family[sample(nrow(family)),]
short_f <- aggregate(idnum ~ f, FUN = max, data=family)
order <- merge(grntax,tax,by.x="o",by.y="name")
order <- order[which(is.na(order$f)),]
order <- subset(order, select = c(idnum,k,p,c,o,f,g,s))
order$idnum <- as.numeric(as.character(order$idnum))
order <- order[sample(nrow(order)),]
short_o <- aggregate(idnum ~ o, FUN = max, data=order)
class <- merge(grntax,tax,by.x="c",by.y="name")
class <- order[which(is.na(order$o)),]
class <- subset(class, select = c(idnum,k,p,c,o,f,g,s))  #######code is incomplete if you want to go past class!! (none in our case)
class <- class[sample(nrow(class)),]
tips_to_plot <- rbind(species,genre,family,order,class)
short_tips_to_plot <- rbindlist(list(short_s,short_g,short_f,short_o))
tips_to_plot$idnum <- as.character(tips_to_plot$idnum)


#write.table(tips_to_plot,file="tips_to_plot_all_forced.txt")
#tips_to_plot <- read.csv("tips_to_plot_all.txt", sep="")
#write.table(short_tips_to_plot,file="shor_tips_to_plot_all_forced.txt")
#short_tips_to_plot <- read.csv("short_tips_to_plot_all.txt", sep="")

#need to go list off taxa that need to be added to tree by hand
print("Some taxa that you supplied were not found in greengenes;")
tax.found <- list()
tax.unfound <- list()
for(i in 1:nrow(tax)){
  if(tax[i,1] %in% tips_to_plot$s){
    tax.found <- c(tax.found,tax[i,1])
  }else if(tax[i,1] %in% tips_to_plot$g){
    tax.found <- c(tax.found,tax[i,1])
  }else if(tax[i,1] %in% tips_to_plot$f){
    tax.found <- c(tax.found,tax[i,1])
  }else if(tax[i,1] %in% tips_to_plot$o){
    tax.found <- c(tax.found,tax[i,1])
  }else if(tax[i,1] %in% tips_to_plot$c){
    tax.found <- c(tax.found,tax[i,1])
  }else{
    tax.unfound <- c(tax.unfound,tax[i,1])
    print(tax[i,1])
  }
}
tax.found <- as.matrix(tax.found)
tax.unfound <- as.matrix(tax.unfound)

#plot complete greengenes taxonomy for all the matches that you provided (big tree)
pruned.tree <- drop.tip(grntree,setdiff(grntree$tip.label,tips_to_plot$idnum))
#plot(pruned.tree)

#plot the randomly selected tips
pruned.tree.short <- drop.tip(grntree,setdiff(grntree$tip.label,short_tips_to_plot$idnum))
plot(pruned.tree.short)

#there are way more taxa in the greengenes taxonomy file than in
#the greengenes tree file. Need to pick ids that exist in greengenes tree file
df <- tips_to_plot
new.tree <- pruned.tree
new.tree2 <- pruned.tree
rownames(df) <- df$idnum
df <- df[,-1]
tax.found2 <- tax.found
row.names(tax.found2) <- tax.found2[,1]


#loops through phylogeny and replaces tip numbers with names
#takes a while
for(i in 1:length(new.tree$tip.label)){
  count <- as.matrix(colSums(!is.na(df[pruned.tree$tip.label[i],])))
  name <- as.character(df[pruned.tree$tip.label[i],sum(count[,1])])
  new.tree$tip.label[i] <- as.character(tax.found2[name,1])
  tax.found2[name,1] <- paste(name,"1",sep=" ")
}
#picked a Clostridiales that is not clustered with the majority of clostridiales.
#I looked at the phylogeny plotted above and picked one of the taxids that is correct.
new.tree$tip.label[new.tree$tip.label=="Clostridiales"] <- "Clostridiales 1"
new.tree$tip.label[[7469]] <- "Clostridiales"
#get rid of all the tips except for 1
new.tree.pick <- drop.tip(new.tree,as.character(tax.found2[,1]))

#plot tree
par(mfrow = c(1,1))
plot(new.tree.pick,use.edge.length=FALSE,cex=0.75)
#exported at 700 x 800 EPS - "PhyloTree.eps"


#need to plot the heatmap now.
#pull out labels that you don't want
#check this out if you are running different data though here;
rm.from <- c("Oxalate.Excreted","Oxalate.Degraded")
rm.to <- c("Oxalate.Consumed","theAnimalNA66","theAnimalNA67","theAnimalNA58","theAnimalNA70","theAnimalNA71","theAnimalNA68","denovo804699")
for(i in 1:length(rm.from)){
  heat <- subset(heat, i.taxa!=rm.from[i]) #drop oxalate from x-axis on plot
}
for(i in 1:length(rm.to)){
  heat <- subset(heat, j.taxa!=rm.to[i])   #drop animal effects from y-axis on heatmap
}
#scaling factor for oxalate.consumed
heat[heat$i.taxa=="Oxalate.Consumed",]$x.scaled <- heat[heat$i.taxa=="Oxalate.Consumed",]$x.scaled * 0.45  #change this value
heat$i.taxa <- factor(heat$i.taxa, levels=unique(heat$i.taxa))
heat$j.taxa <- as.character(heat$j.taxa)
heat$j.taxa <- factor(heat$j.taxa, levels=unique(heat$j.taxa))
heat$j.taxa <- factor(heat$j.taxa, levels = c(tax.unfound,new.tree.pick$tip.label))
wide2 <- dcast(heat, j.taxa ~ i.taxa, value.var = "x.scaled")
heat2 <- melt(wide2)
names(heat2) <- c("j.taxa","i.taxa","x.scaled")
p2 <- ggplot(heat2, aes(x=i.taxa, y=j.taxa,fill=x.scaled)) + 
  geom_tile(colour = "white") + 
  theme(axis.text.y=element_blank(), 
        axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title.y=element_blank(), 
        panel.grid.major=element_blank(), 
        axis.ticks=element_blank()) + 
  scale_fill_gradient2(na.value="grey", 
        low="darkblue", mid='white', high='red4', 
        midpoint = 0,limits=c(-1,1))
p2

p2 <- ggplot(heat2, aes(x=i.taxa, y=j.taxa,fill=x.scaled)) + 
  geom_tile(colour = "white") + 
  theme( 
        axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title.y=element_blank(), 
        panel.grid.major=element_blank(), 
        axis.ticks=element_blank()) + 
  scale_fill_gradient2(na.value="grey", 
                       low="darkblue", mid='white', high='red4', 
                       midpoint = 0,limits=c(-1,1))
p2

#exported at 800 x 800 EPS - "PhyloTree.eps"

#use to index tree on heatmap
ggplot(heat, aes(x=i.taxa, y=j.taxa)) + geom_tile() + theme(axis.text.x=element_blank())
