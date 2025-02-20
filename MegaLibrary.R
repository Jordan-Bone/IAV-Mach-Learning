# update.packages(ask=FALSE)

### A
suppressMessages(library(ape))
suppressMessages(library(ade4))
suppressMessages(library(aplot))
### B
suppressMessages(library(bigmemory))
suppressMessages(library(BiocManager))
suppressMessages(library(brms))
# suppressMessages(library(Bios2cor)) 
suppressMessages(library(bio3d))
suppressMessages(library(beepr))
suppressMessages(library(bayesplot))
### C
suppressMessages(library(circlize))
suppressMessages(library(ComplexHeatmap)) #BiocManager
suppressMessages(library(colorspace))
suppressMessages(library(caret))
suppressMessages(library(caretEnsemble))
suppressMessages(library(coRdon)) #BiocManager
### D
suppressMessages(library(dplyr))
suppressMessages(library(DiagrammeR))
suppressMessages(library(devtools))
suppressMessages(library(deepSNV)) #BiocManager
suppressMessages(library(doParallel))
### E
suppressMessages(library(entropy))
suppressMessages(library(e1071))
suppressMessages(library(Ecfun))
### F
suppressMessages(library(factoextra))
suppressMessages(library(formattable))
suppressMessages(library(forcats))
suppressMessages(library(foreach))
# suppressMessages(library(formal))
### G
suppressMessages(library(gapminder))
# suppressMessages(library(GenomicRanges))
suppressMessages(library(gganimate))
suppressMessages(library(ggplot2))
suppressMessages(library(ggraph))
suppressMessages(library(ggalt))
suppressMessages(library(ggfortify))
suppressMessages(library(grid))
suppressMessages(library(ggtree)) #BiocManager
suppressMessages(library(ggpubr))
suppressMessages(library(ggmsa)) #BiocManager
suppressMessages(library(ggeffects))
suppressMessages(library(ggpattern))
# suppressMessages(library(Gviz))
suppressMessages(library(gifski))
# suppressMessages(library(GenomicFeatures))
suppressMessages(library(ggvenn))
suppressMessages(library(ggtext))
suppressMessages(library(ggrepel))
suppressMessages(library(ggh4x))
### H
suppressMessages(library(hillR))
### I
suppressMessages(library(igraph))
suppressMessages(library(IRanges)) #BiocManager
### J
suppressMessages(library(janitor))
### K
suppressMessages(library(kableExtra))
formals(kable_styling)$bootstrap_options <- c("striped","bordered","condensed")
options(knitr.kable.NA = '')
# suppressMessages(library(KEGGgraph))
# suppressMessages(library(klaR))
suppressMessages(library(kernlab))
### L
suppressMessages(library(latex2exp))
suppressMessages(library(lattice))
suppressMessages(library(lme4))
suppressMessages(library(lmerTest))
suppressMessages(library(lubridate))
suppressMessages(library(loo))
### M
suppressMessages(library(magrittr))
# suppressMessages(library(mgcv))
suppressMessages(library(mice))
suppressMessages(library(microbenchmark))
### N
### O
suppressMessages(library(officer))
### P
suppressMessages(library(paletteer))
suppressMessages(library(pROC))
suppressMessages(library(plyr))
# suppressMessages(library(PopGenome)) 
suppressMessages(library(phangorn))
suppressMessages(library(parallel))
suppressMessages(library(pwr))
# suppressMessages(library(picante))
### Q
suppressMessages(library(qqman))
suppressMessages(library(QuasR)) #BiocManager
suppressMessages(library(qrcode))
### R
suppressMessages(library(rstan))
suppressMessages(library(rstanarm))
suppressMessages(library(RColorBrewer))
suppressMessages(library(readr))
suppressMessages(library(readxl))
suppressMessages(library(reshape2))
suppressMessages(library(rmutil))
# suppressMessages(library(RBGL))
# suppressMessages(library(Rgraphviz))
suppressMessages(library(reticulate))
# suppressMessages(library(rtracklayer))
suppressMessages(library(rdiversity))
suppressMessages(library(rasterdiv))
suppressMessages(library(ranger))
suppressMessages(library(rpart))
suppressMessages(library(rpart.plot))
### S
suppressMessages(library(stringr))
suppressMessages(library(scales))
suppressMessages(library(scatterpie))
suppressMessages(library(seqinr))
suppressMessages(library(sjmisc))
suppressMessages(library(slider))
suppressMessages(library(svMisc))
# suppressMessages(library(sjPlot))
# suppressMessages(library(strengejacke))
### T 
suppressMessages(library(tidygraph))
suppressMessages(library(tidyr))
suppressMessages(library(tidyverse))
suppressMessages(library(treeio)) #BiocManager
suppressMessages(library(tree))
suppressMessages(library(tidytree))
# suppressMessages(library(trackViewer))
suppressMessages(library(taxize))
# suppressMessages(library(taxizedb))
suppressMessages(library(tictoc))
### U
### V
suppressMessages(library(VIM))
suppressMessages(library(vivaldi))
### W
suppressMessages(library(wesanderson))
suppressMessages(library(webshot))
### X
suppressMessages(library(xgboost))

#################

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
per.done <- function(j){ifelse(is.wholenumber(j/10000),print(j),"")}
rsq <- function (x, y) cor(x, y) ^ 2
sl <- data.frame(Segment=c("08NS","04HA","05NP","01PB2","06NA","03PA","02PB1","07MP"),len=c(890,1762,1565,2341,1460,2233,2341,1027))

# Mul	Naive [N~M~]{style="color:#F8766D;"}
# Mul	Vacc [V~M~]{style="color:#7CAE01;"}
# Sin	Naive [N~S~]{style="color:#53CCCD;"}
# Sin Vacc [V~S~]{style="color:#C77CFF;"}

# <span style="background-color: #F8766D">N~M~</span>
# <span style="background-color: #7CAE01">V~M~</span>
# <span style="background-color: #53CCCD">N~S~</span>
# <span style="background-color: #C77CFF">V~S~</span>
