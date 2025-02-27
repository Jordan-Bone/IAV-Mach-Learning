setwd("/Users/jordanbone/Documents/GitHub/IAV-Mach-Learning")
source("MegaLibrary.R")

uid_ref <- read.csv("USED UIDS.csv")
ft_data <- read.csv("Comp/Validata.csv") %>% left_join(uid_ref,by="UID",keep = F)
# ft_data %>% filter(!complete.cases(.)) %>% write.csv("Bad cases.csv")
# for(i in 1:length(ft_data)){if(!is.numeric(ft_data[,i])){print(i)}}

mods <- list.files("Comp",pattern="mclass.rds",full.names = T)
for(j in c(5:7)){
  MOD <- readRDS(mods[j])
  STY <- mods[j] %>% str_split_i("_",2)
  PRT <- mods[j] %>% str_split_i("_",1) %>% str_split_i("\\/",2)
  validate <- ft_data %>% subset(UID %in% subset(uid_ref,Subtype==STY)$UID) %>%
    select("Classification",ends_with(PRT)&!starts_with("ptAcc")) %>% filter(complete.cases(.))
  
  predict_class <- predict(MOD, newdata=validate)
  predict_class$Prediction <- apply(predict_class, 1, function(x) colnames(predict_class)[which.max(x)])
  cm <- confusionMatrix(predict_class$Prediction %>% as.factor %>% droplevels,
                        validate$Classification %>% as.factor %>% droplevels)
  saveRDS(cm,file=paste0(mods[j] %>% str_replace(".rds","_CoMa"),".rds"))
  
  matrix_one_vs_all <- vector("list", length(levels(MOD$models$CTDc.ranger$trainingData$.outcome)))
  for (i in seq_along(matrix_one_vs_all)) {
    positive.class <- levels(MOD$models$CTDc.ranger$trainingData$.outcome)[i]
    print(positive.class)
    matrix_one_vs_all[[i]] <- confusionMatrix(predict_class$Prediction %>% as.factor %>% droplevels, 
                                              validate %>% pull(Classification) %>% as.factor %>% droplevels, 
                                              positive = positive.class)
  t1_df <- cm$overall %>% round(., 3) %>% t() %>%
    cbind(., 
          AUC = multiclass.roc(response = validate %>% pull(Classification), predictor = predict_class) %>% 
            .$auc %>% as.numeric() %>% round(3),
          micro_F1 = matrix_one_vs_all %>% get.micro.f1() %>% round(3),
          macro_F1 = matrix_one_vs_all %>% get.macro.f1() %>% round(3)) 
  }
  if (output == "classwise"){return(matrix_test$byClass %>% reshape2::melt())}
  else {return(as.data.frame(t1_df))}
}  
  
  temp_summary <- data.frame(STY,FEA,PRT,MOD_TYPE,matrix_test$overall)
  temp_split <- data.frame(STY,FEA,PRT,MOD_TYPE,
                           Class=rownames(matrix_test$byClass) %>% str_split_i(":",2) %>% str_trim,
                           matrix_test$byClass)
#####
res <- rbind(cm3$overall,rbind(cm4$overall,cm5$overall)) %>% as.data.frame
res$Protein <- c(rep("HA",3))
res$Holdout <- c("H3N2","H1N2","H1N1")
ggplot(res,aes(Holdout,y=Accuracy,ymin=AccuracyLower,ymax=AccuracyUpper,colour=Holdout))+geom_pointrange()
