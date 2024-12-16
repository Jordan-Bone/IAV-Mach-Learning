source("MegaLibrary.R")

reference_table <- read.csv("RefTable.csv")
prez_subs <- c("H10N1","H10N2","H10N3","H10N4","H10N5","H10N6","H10N7","H10N8","H10N9","H11N1","H11N2","H11N3","H11N4","H11N6","H11N7","H11N8","H11N9","H12N3","H12N5","H12N7","H12N9","H13N2","H13N6","H13N8","H13N9","H14N5","H14N7","H15N5","H15N9","H16N3","H16N9","H17N10","H18N11","H1N1","H1N2","H1N3","H1N5","H1N8","H1N9","H2N1","H2N2","H2N3","H2N4","H2N5","H2N6","H2N7","H2N9","H3N1","H3N2","H3N3","H3N6","H3N7","H3N8","H4N1","H4N2","H4N3","H4N4","H4N5","H4N6","H4N8","H4N9","H5N1","H5N2","H5N3","H5N4","H5N5","H5N6","H5N7","H5N8","H5N9","H6N1","H6N2","H6N3","H6N4","H6N5","H6N6","H6N8","H6N9","H7N1","H7N2","H7N3","H7N4","H7N5","H7N6","H7N7","H7N8","H7N9","H8N3","H8N4","H9N1","H9N2","H9N3","H9N5","H9N6","H9N7","H9N8","H9N9")

fea <- list.files("feats",pattern = ".csv",full.names = T,recursive = T)
# fea <- fea[4:48]
mod_stats <- data.frame()
for(i in 1:length(fea)){
  print(fea[i])
  for(j in 1:length(prez_subs)){
    mo_nm <- fea[i] %>% str_split_i("\\/",2) %>% str_split("\\_")
    tmp_ft <- read.csv(fea[i])
    tmp_ft$ptAcc <- tmp_ft$X %>% str_split_i(":",1)
    tmp_ft <- right_join(reference_table,tmp_ft)
    tmp_ft <- subset(tmp_ft, str_to_lower(str_split_i(tmp_ft$Sample,"\\/",2))!="reassortant")
    tmp_ft <- subset(tmp_ft,Subtype!=prez_subs[j])
    tmp_mosta <- data.frame(Protein=mo_nm[[1]][1],FeatSet=str_replace_all(mo_nm[[1]][8], ".csv",""),Holdout=prez_subs[j],Avian=table(tmp_ft$Class)[1],Mammal=table(tmp_ft$Class)[2])
    fold_indices <- createMultiFolds(tmp_ft$UID, k = 5, times = 1)
    preds <- tmp_ft[,c(-1:-15)] %>% remove_constant %>% names
    # tmp_ft[,c(-2:-16)] %>% remove_constant %>% pivot_longer(cols = -1, names_to = "Ft",values_to = "Val") %>% subset(Val==0) %>% print()
    # subset(tmp_ft,is.na(Class))[,1:15] %>% write.csv(paste0(substr(fea[i],1,10),"Lost.csv"))
    
    mo1 <- suppressWarnings(train(x = tmp_ft %>% select(all_of(preds)),y = tmp_ft %>% pull(Class),method="ranger",metric="ROC",preProc=c("center","scale"),num.trees=1000,importance="impurity",weights=ifelse(tmp_ft$Class =="Aves",(1/table(tmp_ft$Class)[1]) * 0.5,(1/table(tmp_ft$Class)[2]) * 0.5),trControl = trainControl(method = "repeatedcv",number = 5,repeats = 1,index = fold_indices,classProbs = TRUE,savePredictions = TRUE,summaryFunction=twoClassSummary),tuneGrid = expand.grid(.splitrule = "gini",.min.node.size = seq(from = 5, to = 45, length = 3),.mtry = round(sqrt(length(preds))))))
    
    saveRDS(mo1, file=paste0("Models/",mo_nm[[1]][1],"_ex",prez_subs[j],"_",
                             str_replace_all(mo_nm[[1]][8], ".csv",""),".rds"))
    mod_stats <- rbind(mod_stats,tmp_mosta)
  }
}
write.csv(mod_stats,"Models/Mod Stats.csv",row.names = F)

subset(mod_stats,FeatSet=="ctdd") %>% pivot_longer(cols=4:5,names_to = "Class",values_to = "N") %>% ggplot(aes(Holdout,N,colour=Class))+geom_point(alpha=0.6)
subset(mod_stats,Avian+Mammal!=2325)

# # Model Stacks

mods <- list.files("Models",pattern = ".rds",full.names = T,recursive = F)
mods <- unique(str_split_i(mods,"_",1))
# for(i in 1:length(mods)){
for(i in 1:24){
  print(str_split_i(mods[i],"_",1))
  for(j in 1:length(prez_subs)){
    # print(paste0(str_split_i(mods[i],"_",1),"_ex",prez_subs[j],"_2mer.rds"))
    tmer_mod <- readRDS(paste0(str_split_i(mods[i],"_",1),"_ex",prez_subs[j],"_2mer.rds"))
    ctdc_mod <- readRDS(paste0(str_split_i(mods[i],"_",1),"_ex",prez_subs[j],"_ctdc.rds"))
    ctdd_mod <- readRDS(paste0(str_split_i(mods[i],"_",1),"_ex",prez_subs[j],"_ctdd.rds"))
    ctdt_mod <- readRDS(paste0(str_split_i(mods[i],"_",1),"_ex",prez_subs[j],"_ctdt.rds"))
    ctri_mod <- readRDS(paste0(str_split_i(mods[i],"_",1),"_ex",prez_subs[j],"_ctriad.rds"))
    psea_mod <- readRDS(paste0(str_split_i(mods[i],"_",1),"_ex",prez_subs[j],"_pseaac.rds"))
    
    caretStack(c("two_mer"=tmer_mod,"CTDc"=ctdc_mod,"CTDd"=ctdd_mod,
                 "CTDt"=ctdt_mod,"triad"=ctri_mod,"pseaac"=psea_mod),
               # weights=ifelse(tmp_ft$Class =="Aves",
               #                (1/table(tmp_ft$Class)[1]) * 0.5,
               #                (1/table(tmp_ft$Class)[2]) * 0.5),
               preProc=c("center","scale"), metric = "ROC", method="glmnet") %>%
      saveRDS(paste0("Comp/Stack/",str_split_i(mods[i],"_",1),"_",prez_subs[j],".rds"))
  }
}

modc <- list.files("Comp/Stack/Models",pattern = ".rds",full.names = T,recursive = F)
moli <- c()
for(i in 1:length(modc)){
  a <- readRDS("Comp/Stack/Models/ha_H10N1.rds") 
  b <- readRDS("Comp/Stack/Models/ha_H10N2.rds") 
  caretStack(c(a,b),
             # weights=ifelse(tmp_ft$Class =="Aves",
             #                (1/table(tmp_ft$Class)[1]) * 0.5,
             #                (1/table(tmp_ft$Class)[2]) * 0.5),
             preProc=c("center","scale"), metric = "ROC", method="glmnet") %>%
    saveRDS(paste0("Comp/Stack/",str_split_i(mods[i],"_",1),".rds"))
}



modc <- list.files("Comp/Stack/Models",pattern = ".rds",full.names = T,recursive = F)
for(i in 1:length(modc)){
  # ModMeth <- str_split_i(modc[i],"\\/",2)
  ModParam <- str_split_i(modc[i],"\\/",4) %>% str_replace_all(".rds","")
  MoDo <- paste0(str_to_upper(ModParam),"_complete")
  assign(MoDo, readRDS(modc[i]))
}

#######

ha_exH11N3_ctdd_mod <- readRDS("Models/ha_exH11N3_ctdd.rds")
ha_exH10N5_ctdd_mod <- readRDS("Models/ha_exH10N5_ctdd.rds")
ha_exH10N8_ctdd_mod <- readRDS("Models/ha_exH10N8_ctdd.rds")
ha_exH13N8_ctdd_mod <- readRDS("Models/ha_exH13N8_ctdd.rds")
ha_exH3N1_ctdd_mod <- readRDS("Models/ha_exH3N1_ctdd.rds")

ha_exH11N3_ctdc_mod <- readRDS("Models/ha_exH11N3_ctdc.rds")
ha_exH10N5_ctdc_mod <- readRDS("Models/ha_exH10N5_ctdc.rds")
ha_exH10N8_ctdc_mod <- readRDS("Models/ha_exH10N8_ctdc.rds")
ha_exH13N8_ctdc_mod <- readRDS("Models/ha_exH13N8_ctdc.rds")
ha_exH3N1_ctdc_mod <- readRDS("Models/ha_exH3N1_ctdc.rds")

ha_exH11N3_2mer_mod <- readRDS("Models/ha_exH11N3_2mer.rds")
ha_exH10N5_2mer_mod <- readRDS("Models/ha_exH10N5_2mer.rds")
ha_exH10N8_2mer_mod <- readRDS("Models/ha_exH10N8_2mer.rds")
ha_exH13N8_2mer_mod <- readRDS("Models/ha_exH13N8_2mer.rds")
ha_exH3N1_2mer_mod <- readRDS("Models/ha_exH3N1_2mer.rds")

ha_exH3N1_2mer_mod$finalModel$num.samples

tstmo <- caretStack(c(ha_exH11N3_ctdd_mod,ha_exH10N5_ctdd_mod,ha_exH10N8_ctdd_mod,
                      ha_exH13N8_ctdd_mod,ha_exH3N1_ctdd_mod,
                      ha_exH11N3_ctdc_mod,ha_exH10N5_ctdc_mod,ha_exH10N8_ctdc_mod,
                      ha_exH13N8_ctdc_mod,ha_exH3N1_ctdc_mod,
                      ha_exH11N3_2mer_mod,ha_exH10N5_2mer_mod,ha_exH10N8_2mer_mod,
                      ha_exH13N8_2mer_mod,ha_exH3N1_2mer_mod),
                    weights=ifelse(tmp_ft$Class =="Aves",(1/table(tmp_ft$Class)[1]) * 0.5,(1/table(tmp_ft$Class)[2]) * 0.5),
                    preProc=c("center","scale"), metric = "ROC", method="glmnet")

tstmoS <- data.frame(protein="ha",holdout=c(rep(c("H11N3","H10N5","H10N8","H13N8","H3N1"),3)),feature=c(rep("CTDd",5),rep("CTDc",5),rep("2mer",5)),imp=summary(tstmo)$imp,summary(tstmo)$results[-1],row.names = NULL)

ggplot(tstmoS,aes(feature,ymin=value-sd,y=value,ymax=value+sd,colour=imp))+geom_pointrange()
ggplot(tstmoS,aes(feature,value))+geom_violin(scale = "area")+geom_jitter(aes(colour=holdout),height = 0, width=0.1,alpha=0.8)


