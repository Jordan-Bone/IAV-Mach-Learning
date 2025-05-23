setwd("/Users/jordanbone/Documents/GitHub/IAV-Mach-Learning")
source("MegaLibrary.R")
library(MLmetrics)
library(foreach) # enables use of the %do% operator. This is another way of applying a function over and over to different items of a list or vector. But it has a parallelisable version to run these tasks in parallel (best done on a computing cluster, we can do this with %dopar%)

# Run variable importances? NB will take 6hrs per protein = 48 hours for all proteins on single LB desktop
# run_varimp <- TRUE
run_varimp <- FALSE

# Select 2class or mclass models
MODTYPE <- "mclass"

# Set random seed for reproducibility
set.seed(1303)

uid_ref <- read.csv("USED UIDS.csv")
ft_data <- read.csv("Comp/Validata.csv") %>% left_join(uid_ref,by="UID",keep = F)
# ft_data %>% filter(!complete.cases(.)) %>% write.csv("Bad cases.csv")
# for(i in 1:length(ft_data)){if(!is.numeric(ft_data[,i])){print(i)}}
master2 <- readRDS("Comp/ConfMatrix 2class Master.rds")
masterM <- readRDS("Comp/ConfMatrix Multiclass Master.rds")

# Select protein
for (PRT in c("PB2","PB1","PA","HA","NP","NA","M1","NS1")){
  
  prps <- c()
  prds <- c()
  vlis <- c()
  pred_vals <- c()
  mismatches <- c()
  
  mods <- list.files("Comp",pattern=paste0(MODTYPE,".rds"),full.names = T) %>% .[grepl(PRT, .)]     # select protein
  for(j in c(1:length(mods))){
    # print(mods[j])
    
    MOD <- readRDS(mods[j])
    STY <- mods[j] %>% str_split_i("_",2)
    validate <- ft_data %>% subset(UID %in% subset(uid_ref,Subtype==STY)$UID) %>%
      select("Classification","UID",ends_with(PRT)&!starts_with("ptAcc")) %>% filter(complete.cases(.))
    if(MODTYPE == "2class"){
      validate$Classification <- ifelse(validate$Classification=="Avian","Avian","Mammal")
    }
    
    predict_class <- predict(MOD, newdata=validate)
    prps <- bind_rows(prps, predict_class)
    prvals <- predict_class %>% unlist
    predict_class$Prediction <- apply(predict_class, 1, function(x) colnames(predict_class)[which.max(x)])
    
    prds <- append(prds,predict_class$Prediction)
    vlis <- append(vlis,validate$Classification)
    pred_vals <- append(pred_vals,prvals)
    
    mismat <- cbind(validate$UID,PRT,STY,validate$Classification,predict_class) %>% data.frame()
    mismatches <- bind_rows(mismatches,mismat)
  }
  write.csv(mismatches,paste0("Comp/",PRT,"_predictions.csv"),row.names = F,quote = F)
  
  # confusionMatrix(prds %>% as.factor, vlis %>% as.factor) %>% saveRDS(paste0("Comp/",PRT," 2class ConfMatrix.rds"))
  # multiclass.roc(response = validate %>% pull(Classification),predictor = predict_class) %>% print
  
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
  
  # predictor_values <- matrix(ncol=length(levels(as.factor(vlis))),data=pred_vals)
  # colnames(predictor_values) <- levels(as.factor(vlis))
  # AUC = multiclass.roc(response = as.factor(vlis), predictor = predictor_values,
  # levels=levels(as.factor(prds))) %>% .$auc %>% as.numeric() %>% round(3)
  # 
  # f1mi <- F1_Score_micro(vlis,prds)
  # f1ma <- F1_Score_macro(vlis,prds)
  # data.frame(PRT,f1mi,f1ma,AUC) %>% print
  
  ###############
  # VARIABLE IMPORTANCE
  ###############
  
  if(run_varimp == TRUE){
    
    AUC_base = multiclass.roc(response = vlis,
                              predictor = prps) %>%
      .$auc %>%
      as.numeric()
    
    
    varnames <- ft_data %>% select(ends_with(PRT)&!starts_with("ptAcc")) %>% filter(complete.cases(.)) %>% names
    
    varimp_perm <- foreach(varname = varnames) %do% {
      
      start <- Sys.time()
      
      valid_shuffle <- ft_data %>%
        select("UID","Classification",ends_with(PRT)&!starts_with("ptAcc")) %>%
        filter(complete.cases(.)) %>%
        mutate_at(vars(varname), sample) # permute the single given column
      
      # initalise empty vectors
      prps <- c()
      prds <- c()
      vlis <- c()
      
      for(j in c(1:length(mods))){
        MOD <- readRDS(mods[j])
        STY <- mods[j] %>% str_split_i("_",2)
        validate <- valid_shuffle %>% subset(UID %in% subset(uid_ref,Subtype==STY)$UID) %>% filter(complete.cases(.)) %>% select(-UID)
        if(mods[j] %>% str_contains("2class")){
          validate$Classification <- ifelse(validate$Classification=="Avian","Avian","Mammal")
        }
        
        predict_class <- predict(MOD, newdata=validate)
        
        prps <- bind_rows(prps, predict_class)
        
        predict_class$Prediction <- apply(predict_class, 1, function(x) colnames(predict_class)[which.max(x)])
        
        prds <- append(prds,predict_class$Prediction)
        vlis <- append(vlis,validate$Classification)
      }  
      
      AUC_perm = multiclass.roc(response = vlis,
                                predictor = prps) %>%
        .$auc %>%
        as.numeric()
      
      end <- Sys.time()
      
      print(paste("Variable importance:", varname))
      print(end - start)
      
      return(data.frame(var = varname, AUC_loss = AUC_base - AUC_perm))
      
    } # END loop over features
    
    # Takes around 20 seconds per feature = ~6hrs for all 1330 feats for one protein and modeltype
    varimp_perm %>% bind_rows %>% write.csv(paste0("Varimp/varimp_", PRT, "_", MODTYPE, ".csv"))
    
  } # END varimp calculation
  
} # END loop over proteins
