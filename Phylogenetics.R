setwd("/Users/jordanbone/Documents/Coding_Ref/Phylogenetics/")
source("MegaLibrary.R")
# install_github("rgriff23/btw")
library("btw")
library(tidyr)
library(stringr)
library(magrittr)
# ape::ace function could also be used to do ancestral reconstruction.

m1_root <- read.nexus("iqtree/alpha iqtree/M1.nexus")

### For discrete characters:
# x <- ifelse(m1_root$tip.label %>% str_split_i("_",3)=="Avian",0,1) %>% as.factor
y <- m1_root$tip.label %>% str_split_i("_",3) %>% as.factor
ans <- ace(y, phy=m1_root,type="discrete")

### Some random data...
data(bird.orders)
x <- rnorm(23)
### Compare the three methods for continuous characters:
ace(x, bird.orders)
ace(x, bird.orders, method = "pic")
ace(x, bird.orders, method = "GLS", corStruct = corBrownian(1, bird.orders))
### For discrete characters:
x <- factor(c(rep(0, 5), rep(1, 18)))
ans <- ace(x, bird.orders, type = "d")
#### Showing the likelihoods on each node:
plot(bird.orders, type = "c", FALSE, label.offset = 1)
co <- c("blue", "yellow")
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 2, adj = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.75)

library(phytools)
data(anoletree)
## this is just to pull out the tip states from the
## data object - normally we would read this from file
x<-getStates(anoletree,"tips")
tree<-anoletree
rm(anoletree)
tree
plotTree(tree,type="fan",fsize=0.8,ftype="i")
cols<-setNames(palette()[1:length(unique(x))],sort(unique(x)))
tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree)),fsize=0.8)
fitER<-ace(x,tree,model="ER",type="discrete")
fitER
round(fitER$lik.anc,3)

plotTree(tree,type="fan",fsize=0.8,ftype="i")
nodelabels(node=1:tree$Nnode+Ntip(tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.5)
tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree)),fsize=0.8)
