setwd("/Users/jordanbone/Documents/Coding_Ref/Phylogenetics/BayesTraits/")
source("MegaLibrary.R")
# install_github("rgriff23/btw")
library("btw")
# ape::ace function could also be used to do ancestral reconstruction.

HA_tree <- read.nexus("IAV/HA_tree.nexus")

# ace(x, HA_tree, type = "discrete", 
#     method = if (type == "continuous")
#   "REML" else "ML", CI = TRUE,
#   model = if (type == "continuous") "BM" else "ER",
#   scaled = TRUE, kappa = 1, corStruct = NULL, ip = 0.1,
#   use.expm = FALSE, use.eigen = TRUE, marginal = FALSE)


ace(HA_tree$tip.label, HA_tree, type = "discrete", method = "ML", CI = TRUE,
  model = "ER", scaled = TRUE, kappa = 1, corStruct = NULL, ip = 0.1,
  use.expm = FALSE, use.eigen = TRUE, marginal = FALSE)
