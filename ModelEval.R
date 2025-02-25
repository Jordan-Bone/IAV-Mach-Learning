setwd("/Users/jordanbone/Documents/GitHub/IAV-Mach-Learning")
# source("MegaLibrary.R")
library(dplyr)
library(stringr)
library(tidyr)
library(caret)
library(caretEnsemble)

uid_ref <- read.csv("USED UIDS.csv")
ft_data <- read.csv("Comp/Validata.csv") %>% left_join(uid_ref,by="UID",keep = F)
# ft_data %>% filter(!complete.cases(.)) %>% write.csv("Bad cases.csv")
# for(i in 1:length(ft_data)){
#   if(!is.numeric(ft_data[,i])){
#     # print(ft_data[,i])
#     print(i)
#   }
# }


mods <- list.files("Comp",pattern="mclass.rds",full.names = T)
predict_class <- c()
predict_prob <- data.frame()
VALD <- data.frame()

for(j in c(1:3)){
  MOD <- readRDS(mods[j])
  STY <- mods[j] %>% str_split_i("_",2)
  # FEA <- mods[j] %>% str_split_i("_",2)
  PRT <- mods[j] %>% str_split_i("_",1) %>% str_split_i("\\/",2)
  validate <- ft_data %>% subset(UID %in% subset(uid_ref,Subtype==STY)$UID) %>%
    select("UID",ends_with(PRT)&!starts_with("ptAcc"))
  validate <- validate %>% filter(complete.cases(.)) %>% select(-"UID") %>% mutate_if(is.character,as.numeric)
  
  predict_class_test <- predict(MOD, newdata=validate, type="raw")
  # %>% unlist %>% as.factor
  predict_prob_test <- predict(MOD, newdata=validate, type="prob") %>% bind_rows
  
  predict_class <- c(predict_class,predict_class_test)
  predict_prob <- bind_rows(predict_prob,predict_prob_test)
  VALD <- bind_rows(VALD,validate %>% pull(Classification))
}
predict_class <- predict_class %>% as.character() %>%
  str_replace_all("1","Avian") %>% str_replace_all("2","Canidae") %>% str_replace_all("3","Equidae") %>%
  str_replace_all("4","Hominidae") %>% str_replace_all("5","Phyllostomidae") %>% str_replace_all("6","Suidae")

matrix_test <- confusionMatrix(predict_class %>% as.factor %>% droplevels,
                              VALD %>% pull(label) %>% as.factor %>% droplevels)
  
  matrix_one_vs_all <- vector("list", length(levels(MOD$trainingData$.outcome)))
  for (i in seq_along(matrix_one_vs_all)) {
    positive.class <- levels(MOD$trainingData$.outcome)[i]
    print(positive.class)
    matrix_one_vs_all[[i]] <- confusionMatrix(predict_class_test %>% droplevels, 
                                              validate %>% pull(label) %>% as.factor, 
                                              positive = positive.class)}
  
  t1_df <- matrix_test$overall %>% round(., 3) %>% t() %>%
    cbind(., 
          AUC = multiclass.roc(response = validate %>% pull(label), predictor = predict_prob_test) %>% 
            .$auc %>% as.numeric() %>% round(3),
          micro_F1 = matrix_one_vs_all %>% get.micro.f1() %>% round(3),
          macro_F1 = matrix_one_vs_all %>% get.macro.f1() %>% round(3)) 
  
  if (output == "classwise"){
    return(matrix_test$byClass %>% reshape2::melt())
  } else {
    return(as.data.frame(t1_df))
    
  }
  
  
  temp_summary <- data.frame(STY,FEA,PRT,MOD_TYPE,matrix_test$overall)
  temp_split <- data.frame(STY,FEA,PRT,MOD_TYPE,
                           Class=rownames(matrix_test$byClass) %>% str_split_i(":",2) %>% str_trim,
                           matrix_test$byClass)
#####
  

MOD3 <- readRDS(mods[3])
MOD4 <- readRDS(mods[4])
MOD5 <- readRDS(mods[5])

validate3 <- ft_data %>% subset(UID %in% subset(uid_ref,Subtype=="H17N10")$UID) %>%
  select(ends_with(".ha")&!starts_with("ptAcc")) %>% filter(complete.cases(.))
validate4 <- ft_data %>% subset(UID %in% subset(uid_ref,Subtype=="H18N11")$UID) %>%
  select(ends_with(".ha")&!starts_with("ptAcc")) %>% filter(complete.cases(.))
validate5 <- ft_data %>% subset(UID %in% subset(uid_ref,Subtype=="H1N1")$UID) %>%
    select("Classification",ends_with(".ha")&!starts_with("ptAcc")) %>% filter(complete.cases(.))

# predict_class <- predict(object=MOD5, newdata=validate5) %>% unlist
predict_class <- predict(object=MOD5, newdata=validate5)

predict_class$Prediction <- apply(predict_class, 1, function(x) colnames(predict_class)[which.max(x)])


confusionMatrix(predict_class$Prediction %>% as.factor %>% droplevels,
                validate5$Classification %>% as.factor %>% droplevels)
