setwd("/Users/jordanbone/Documents/Coding_Ref/Phylogenetics/")
source("MegaLibrary.R")
# install_github("rgriff23/btw")
library("btw")
library(tidyr)
library(stringr)
library(magrittr)

uid_ref <- read.csv("/Users/jordanbone/Documents/GitHub/IAV-Mach-Learning/USED UIDS.csv")
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

##########################################

tree <- read.nexus("iqtree/IAV/HA_1.nexus")
anno <- read.csv("iqtree/IAV/HA_anno.tsv",sep="\t")

# groupInfo <- split(tree$tip.label, gsub("_\\w+", "", tree$tip.label))
# groupInfo <- split(tree$tip.label, str_split_i(tree$tip.label,"\\_",3))
chiroptera <- groupOTU(tree, anno)

ggtree(chiroptera,aes(colour=group))+geom_nodelab(aes(label = node), hjust = -0.5)
ggtree(extract.clade(tree, 2727))+geom_nodelab(aes(label = node), hjust = -0.5)

##########################################
library(PhylogeneticEM)
library(TreeSim)
set.seed(17920902)
ntaxa = 80
tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0,
                                 age = 1, mrca = TRUE)[[1]]

params <- params_process("OU",                             ## Process
                         p = 2,                            ## Dimension
                         variance = diag(0.5, 2, 2) + 0.5, ## Rate matrix
                         selection.strength = 3,           ## Selection Strength
                         random = TRUE,                    ## Root is random
                         stationary.root = TRUE,           ## Root is stationary
                         edges = c(16, 81, 124),           ## Positions of the shifts
                         values = cbind(c(5, 4),           ## Values of the shifts
                                        c(-4, -5),
                                        c(5, -3)))
plot(params, phylo = tree, traits = 1, value_in_box = TRUE, shifts_bg = "white")
plot(params, phylo = tree, traits = 2, value_in_box = TRUE, shifts_bg = "white")

sim <- simul_process(params, tree)

data <- extract(sim,             ## The simul_process object
                what = "states", ## We want the actual values
                where = "tips")  ## Only at the tips of the tree

rownames(data) <- c("A", "B")

plot(params, phylo = tree, data = data)

## Inference

res <- PhyloEM(phylo = tree,
               Y_data = data,
               process = "scOU",                   ## scalar OU model
               random.root = TRUE,                 ## Root is stationary (true model)
               stationary.root = TRUE,
               alpha = alpha_grid,                 ## On a grid of alpha
               K_max = 10,                         ## Maximal number of shifts
               parallel_alpha = TRUE,              ## This can be set to TRUE for
               Ncores = 2)                         ## parallel computations
res
plot(res)
plot_criterion(res)
plot(res, params = params_process(res, method.selection = "DDSE"))
plot(res, params = params_process(res, K = 3, alpha = 1))

params_process(res, K = 3, alpha = 1)$shifts
params_process(res, K = 3, alpha = 3)$shifts
