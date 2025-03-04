setwd("/Users/jordanbone/Documents/Coding_Ref/Phylogenetics/")
source("MegaLibrary.R")
# install_github("rgriff23/btw")
library("btw")
# ape::ace function could also be used to do ancestral reconstruction.

m1_root <- read.nexus("iqtree/alpha iqtree/M1.nexus")
### For discrete characters:
x <- ifelse(m1_root$tip.label %>% str_split_i("_",3)=="Avian",0,1) %>% as.factor
ans <- ace(x, phy=m1_root, type="discrete", use.expm=TRUE, use.eigen=FALSE)

#### Showing the likelihoods on each node:
plot(bird.orders, type = "c", FALSE, label.offset = 1)
co <- c("blue", "yellow")
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 2, adj = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.75)


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
