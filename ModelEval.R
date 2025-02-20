setwd("/Users/jordanbone/Documents/GitHub/IAV-Mach-Learning")
source("MegaLibrary.R")

uid_ref <- read.csv("USED UIDS.csv")
names(uid_ref)[3] <- "label"
mods <- list.files("Models",pattern="Mclass.rds",full.names = T)
predict_class <- c()
predict_prob <- data.frame()
VALD <- data.frame()

# for(j in 1:length(mods)){
for(j in c(1:25)){
  MOD <- readRDS(mods[j])
  STY <- mods[j] %>% str_split_i("_",3)
  FEA <- mods[j] %>% str_split_i("_",2)
  PRT <- mods[j] %>% str_split_i("_",1) %>% str_split_i("\\/",2)
  og_data <- read.csv(paste0("feats/",FEA,".csv"))
  validate <- og_data %>% subset(UID %in% subset(uid_ref,Subtype==STY)$UID) %>% select("UID",ends_with(PRT))
  validate <- left_join(validate,uid_ref)
  
  predict_class_test <- predict(MOD, newdata=validate, type="raw") %>% unlist %>% as.factor
  predict_prob_test <- predict(MOD, newdata=validate, type="prob") %>% bind_rows
  
  predict_class <- c(predict_class,predict_class_test)
  predict_class <- predict_class %>% replace(1,"Avian")
  predict_prob <- bind_rows(predict_prob,predict_prob_test)
  VALD <- bind_rows(VALD,validate)
}
  matrix_test <- confusionMatrix(predict_class %>% droplevels,
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
