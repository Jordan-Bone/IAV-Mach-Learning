setwd("/Users/jordanbone/Documents/GitHub/IAV-Mach-Learning")
source("MegaLibrary.R")
library(MLmetrics)

uid_ref <- read.csv("USED UIDS.csv")
ft_data <- read.csv("Comp/Validata.csv") %>% left_join(uid_ref,by="UID",keep = F)
# ft_data %>% filter(!complete.cases(.)) %>% write.csv("Bad cases.csv")
# for(i in 1:length(ft_data)){if(!is.numeric(ft_data[,i])){print(i)}}
master2 <- readRDS("Comp/ConfMatrix 2class Master.rds")
masterM <- readRDS("Comp/ConfMatrix Multiclass Master.rds")
prds <- c()
vlis <- c()
pred_vals <- c()

mods <- list.files("Comp",pattern="mclass.rds",full.names = T)
for(j in c(1:length(mods))){
  print(j)
  print(mods[j])
  
  MOD <- readRDS(mods[j])
  STY <- mods[j] %>% str_split_i("_",2)
  PRT <- mods[j] %>% str_split_i("_",1) %>% str_split_i("\\/",2)
  validate <- ft_data %>% subset(UID %in% subset(uid_ref,Subtype==STY)$UID) %>%
    select("Classification",ends_with(PRT)&!starts_with("ptAcc")) %>% filter(complete.cases(.))
  if(mods[j] %>% str_contains("2class")){
    validate$Classification <- ifelse(validate$Classification=="Avian","Avian","Mammal")
  }
  
  predict_class <- predict(MOD, newdata=validate)
  prvals <- predict_class %>% unlist
  predict_class$Prediction <- apply(predict_class, 1, function(x) colnames(predict_class)[which.max(x)])
  
  # try({
  #   cm <- confusionMatrix(predict_class$Prediction %>% as.factor %>% droplevels,
  #                         validate$Classification %>% as.factor %>% droplevels)
  #   saveRDS(cm,file=paste0(mods[j] %>% str_replace(".rds","_CoMa"),".rds"))
  # },silent = T)
  # 
  # predict_class <- c()
  # validate <- c()
  prds <- append(prds,predict_class$Prediction)
  vlis <- append(vlis,validate$Classification)
  pred_vals <- append(pred_vals,prvals)
}

# confusionMatrix(prds %>% as.factor, vlis %>% as.factor) %>% saveRDS("ConfMatrix Master 2class.rds")

# matrix_one_vs_all <- vector("list", length(levels(MOD$models$CTDc.ranger$trainingData$.outcome)))
# for (i in seq_along(matrix_one_vs_all)) {
# positive.class <- levels(MOD$models$CTDc.ranger$trainingData$.outcome)[i]
# # print(positive.class)
# matrix_one_vs_all[[i]] <- confusionMatrix(predict_class$Prediction %>% droplevels, 
#                                           validate %>% pull(Classification) %>% as.factor %>% droplevels, 
#                                           positive = positive.class)
# t1_df <- cm$overall %>% round(., 3) %>% t() %>%
#   cbind(.,
#         AUC = multiclass.roc(response = validate %>% pull(Classification),
#                              predictor = predict_class) %>%
#           .$auc %>% as.numeric() %>% round(3),
#         micro_F1 = matrix_one_vs_all %>% get.micro.f1() %>% round(3),
#         macro_F1 = matrix_one_vs_all %>% get.macro.f1() %>% round(3))
# if (output == "classwise"){return(matrix_test$byClass %>% reshape2::melt())}
# else {return(as.data.frame(t1_df))}
# }

F1_Score_micro(vlis,prds) %>% percent(accuracy=0.01) # 83.57%
F1_Score_macro(vlis,prds) %>% percent(accuracy=0.01) # 47.94%

# AUC = multiclass.roc(response = vlis %>% as.factor, predictor = prds %>% as.factor,
                     # levels=levels(as.factor(prds)))
# %>% .$auc %>% as.numeric() %>% round(3)