setwd("/Users/jordanbone/Documents/GitHub/IAV-Mach-Learning")
source("MegaLibrary.R")

# Model Stacks

ft_mods <- list.files("Models", pattern=".rds",full.names = T)
stacks_ls <- unique(paste0(str_split_i(ft_mods,"_",1) %>% str_split_i("/",2),"_",str_split_i(ft_mods,"_",3)))
for(i in 1:length(stacks_ls)){
# for(i in 2:10){
  tic(paste("Stacking",stacks_ls[i]))
  print(percent(i/length(stacks_ls)))
  print(stacks_ls[i])
  
  tmer_2c_mod <- readRDS(paste0("Models/",str_split_i(stacks_ls[i],"_",1),"_2mer_",str_split_i(stacks_ls[i],"_",2),"_holdout_2class.rds"))
  ctdc_2c_mod <- readRDS(paste0("Models/",str_split_i(stacks_ls[i],"_",1),"_CTDc_",str_split_i(stacks_ls[i],"_",2),"_holdout_2class.rds"))
  ctdd_2c_mod <- readRDS(paste0("Models/",str_split_i(stacks_ls[i],"_",1),"_CTDd_",str_split_i(stacks_ls[i],"_",2),"_holdout_2class.rds"))
  ctdt_2c_mod <- readRDS(paste0("Models/",str_split_i(stacks_ls[i],"_",1),"_CTDt_",str_split_i(stacks_ls[i],"_",2),"_holdout_2class.rds"))
  ctri_2c_mod <- readRDS(paste0("Models/",str_split_i(stacks_ls[i],"_",1),"_triad_",str_split_i(stacks_ls[i],"_",2),"_holdout_2class.rds"))
  psea_2c_mod <- readRDS(paste0("Models/",str_split_i(stacks_ls[i],"_",1),"_paac_",str_split_i(stacks_ls[i],"_",2),"_holdout_2class.rds"))
  #
  caretStack(c("two_mer"=tmer_2c_mod,"CTDc"=ctdc_2c_mod,"CTDd"=ctdd_2c_mod,
               "CTDt"=ctdt_2c_mod,"triad"=ctri_2c_mod,"pseaac"=psea_2c_mod),weights=ifelse(psea_2c_mod$trainingData$.outcome =="Avian",1/(table(psea_2c_mod$trainingData$.outcome)[1]/length(psea_2c_mod$trainingData$.outcome)),1/(table(psea_2c_mod$trainingData$.outcome)[2]/length(psea_2c_mod$trainingData$.outcome))),preProc=c("center","scale"), metric = "ROC", method="glmnet") %>% saveRDS(paste0("Comp/",stacks_ls[i],"_2class.rds"))
  
  tmer_mc_mod <- readRDS(paste0("Models/",str_split_i(stacks_ls[i],"_",1),"_2mer_",str_split_i(stacks_ls[i],"_",2),"_holdout_mclass.rds"))
  ctdc_mc_mod <- readRDS(paste0("Models/",str_split_i(stacks_ls[i],"_",1),"_CTDc_",str_split_i(stacks_ls[i],"_",2),"_holdout_mclass.rds"))
  ctdd_mc_mod <- readRDS(paste0("Models/",str_split_i(stacks_ls[i],"_",1),"_CTDd_",str_split_i(stacks_ls[i],"_",2),"_holdout_mclass.rds"))
  ctdt_mc_mod <- readRDS(paste0("Models/",str_split_i(stacks_ls[i],"_",1),"_CTDt_",str_split_i(stacks_ls[i],"_",2),"_holdout_mclass.rds"))
  ctri_mc_mod <- readRDS(paste0("Models/",str_split_i(stacks_ls[i],"_",1),"_triad_",str_split_i(stacks_ls[i],"_",2),"_holdout_mclass.rds"))
  psea_mc_mod <- readRDS(paste0("Models/",str_split_i(stacks_ls[i],"_",1),"_paac_",str_split_i(stacks_ls[i],"_",2),"_holdout_mclass.rds"))
  
  w2 <- psea_mc_mod$trainingData %>% group_by(.outcome) %>% 
    summarise(N=length(.outcome)/nrow(psea_mc_mod$trainingData))
  
  caretStack(c("two_mer"=tmer_mc_mod,"CTDc"=ctdc_mc_mod,
               "CTDd"=ctdd_mc_mod,"CTDt"=ctdt_mc_mod,
               "triad"=ctri_mc_mod,"pseaac"=psea_mc_mod),
weights=ifelse(psea_mc_mod$trainingData$.outcome =="Avian",1/subset(w2,.outcome=="Avian")$N,
ifelse(psea_mc_mod$trainingData$.outcome =="Canidae",1/subset(w2,.outcome=="Canidae")$N,
ifelse(psea_mc_mod$trainingData$.outcome =="Equidae",1/subset(w2,.outcome=="Equidae")$N,
ifelse(psea_mc_mod$trainingData$.outcome =="Hominidae",1/subset(w2,.outcome=="Hominidae")$N,
ifelse(psea_mc_mod$trainingData$.outcome =="Suidae",1/subset(w2,.outcome=="Suidae")$N,
1/subset(w2,.outcome=="Phyllostomidae")$N))))),
             preProc=c("center","scale"),method="ranger",num.trees=10) %>% saveRDS(paste0("Comp/",stacks_ls[i],"_mclass.rds"))
  toc()
}
#####

#####
# tc_raw <- list.files("Comp",pattern="2class.rds",full.names = T)
# model_compare <- data.frame()
# tic("Model Outputs")
# for(i in 1:length(tc_raw)){
#   print(i)
#   temp_mo <- readRDS(tc_raw[i])
#   temp_compare <- data.frame(PROT=str_split_i(tc_raw[i],"_",1) %>% str_split_i("/",2),HOLD=str_split_i(tc_raw[i],"_",2),FEAT=str_split_i(names(summary(temp_mo)$imp),"\\.",1),IMP=summary(temp_mo)$imp,summary(temp_mo)$results[-1][,3],summary(temp_mo)$results[-1][,4],row.names = NULL)
#   model_compare <- rbind(model_compare,temp_compare)
# }
# toc()
# write.csv(model_compare,"Comp/2 class model stats.csv",row.names = F)
#####
tc_stats <- read.csv("Comp/2 class model stats.csv")
ggplot(tc_stats,aes(FEAT,IMP,colour=HOLD))+geom_point(position=position_dodge(width=0.5),alpha=0.6)+facet_wrap(~PROT,ncol=4)+theme(legend.position="bottom",axis.text.x = element_text(angle = 90))+ guides(col = guide_legend(nrow = 3))

ggplot(tc_stats,aes(FEAT,ymin=value-sd,y=value,ymax=value+sd,colour=HOLD))+ geom_pointrange(position=position_dodge(width=0.4),alpha=0.6)+facet_wrap(~PROT,ncol=4)+theme(legend.position="bottom")+ guides(col = guide_legend(nrow = 3))


#####
# mc_raw <- list.files("Comp",pattern="mclass.rds",full.names = T)
# model_compare <- data.frame()
# tic("Model Outputs")
# for(i in 1:length(mc_raw)){
#   print(i)
#   temp_mo <- readRDS(mc_raw[i])
#   temp_compare <- data.frame(PROT=str_split_i(mc_raw[i],"_",1) %>% str_split_i("/",2),HOLD=str_split_i(mc_raw[i],"_",2),FEAT=str_split_i(names(summary(temp_mo)$imp),"\\.",1),CLASS=str_split_i(names(summary(temp_mo)$imp),"\\.",2) %>% str_split_i("_",2),IMP=summary(temp_mo)$imp,summary(temp_mo)$results[-1][,3],summary(temp_mo)$results[-1][,4],row.names = NULL)
#   model_compare <- rbind(model_compare,temp_compare)
# }
# toc()
# write.csv(model_compare,"Comp/Multiclass model stats.csv",row.names = F)
#####

mc_stats <- read.csv("Comp/Multiclass model stats.csv")
ggplot(mc_stats,aes(FEAT,IMP,colour=HOLD))+geom_point(position=position_dodge(width=0.5),alpha=0.6)+facet_wrap(~PROT,ncol=4)+theme(legend.position="bottom",axis.text.x = element_text(angle = 90))+ guides(col = guide_legend(nrow = 3))
ggplot(mc_stats,aes(FEAT,ymin=value-sd,y=value,ymax=value+sd,colour=HOLD))+ geom_pointrange(position=position_dodge(width=0.4),alpha=0.6)+facet_wrap(~PROT,ncol=4)+theme(legend.position="bottom")+ guides(col = guide_legend(nrow = 3))

# model_compare$CLASS <- ifelse(model_compare$CLASS=="Canidae","Dog",ifelse(model_compare$CLASS=="Equidae","Horse",ifelse(model_compare$CLASS=="Hominidae","Human",ifelse(model_compare$CLASS=="Phyllostomidae","Bat","Swine"))))
# 
# new <- c(emojifont::emoji("dog"),emojifont::emoji("horse"),emojifont::emoji("man"),emojifont::emoji("bat"),emojifont::emoji("pig")) 
# names(new) <- c("Dog", "Human", "Horse", "Bat","Swine") 
# 
# ggplot(model_compare,aes(FEAT,IMP,colour=HOLD))+geom_point(position=position_dodge(width=0.4))+facet_wrap(~CLASS,labeller = labeller(CLASS = new))+theme(legend.position = "bottom",axis.text.x = element_text(angle = 90))+ guides(col = guide_legend(nrow = 3))
# ggplot(model_compare,aes(FEAT,ymin=value-sd,y=value,ymax=value+sd,colour=HOLD))+ geom_pointrange(alpha=0.8,position=position_dodge(width=0.4))+facet_wrap(~CLASS)+theme(legend.position = "bottom")+ guides(col = guide_legend(nrow = 3))

##### Test #####
# tmer_2c_mod <- readRDS("Models/HA_2mer_H16N3_holdout_2class.rds")
# ctdc_2c_mod <- readRDS("Models/HA_CTDc_H16N3_holdout_2class.rds")
# ctdd_2c_mod <- readRDS("Models/HA_CTDd_H16N3_holdout_2class.rds")
# ctdt_2c_mod <- readRDS("Models/HA_CTDt_H16N3_holdout_2class.rds")
# ctri_2c_mod <- readRDS("Models/HA_triad_H16N3_holdout_2class.rds")
# psea_2c_mod <- readRDS("Models/HA_paac_H16N3_holdout_2class.rds")
#  
# c2 <- caretStack(c("two_mer"=tmer_2c_mod,"CTDc"=ctdc_2c_mod,"CTDd"=ctdd_2c_mod,
#              "CTDt"=ctdt_2c_mod,"triad"=ctri_2c_mod,"pseaac"=psea_2c_mod),
#            weights=ifelse(psea_2c_mod$trainingData$.outcome =="Avian",1/(table(psea_2c_mod$trainingData$.outcome)[1]/length(psea_2c_mod$trainingData$.outcome)),1/(table(psea_2c_mod$trainingData$.outcome)[2]/length(psea_2c_mod$trainingData$.outcome))), preProc=c("center","scale"), metric = "ROC", method="glmnet")
# summary(c2)
# 
tmer_mc_mod <- readRDS("Models/HA_2mer_H10N7_holdout_Mclass.rds")
ctdc_mc_mod <- readRDS("Models/HA_CTDc_H10N7_holdout_Mclass.rds")
ctdd_mc_mod <- readRDS("Models/HA_CTDd_H10N7_holdout_Mclass.rds") #
ctdt_mc_mod <- readRDS("Models/HA_CTDt_H10N7_holdout_Mclass.rds") #
ctri_mc_mod <- readRDS("Models/HA_triad_H10N7_holdout_Mclass.rds")
psea_mc_mod <- readRDS("Models/HA_paac_H10N7_holdout_Mclass.rds")

w2 <- psea_mc_mod$trainingData %>% group_by(.outcome) %>% summarise(N=length(.outcome)/nrow(psea_mc_mod$trainingData))
fold_indices <- createMultiFolds(psea_mc_mod$trainingData$.outcome, k = 5, times = 1)

mc <- caretStack(as.caretList(list("CTDc"=ctdc_mc_mod,"CTDd"=ctdd_mc_mod,"CTDt"=ctdt_mc_mod,
                                   "triad"=ctri_mc_mod,"pseaac"=psea_mc_mod,"twomer"=tmer_mc_mod)),
weights=ifelse(psea_mc_mod$trainingData$.outcome =="Avian",1/subset(w2,.outcome=="Avian")$N,
ifelse(psea_mc_mod$trainingData$.outcome =="Canidae",1/subset(w2,.outcome=="Canidae")$N,
ifelse(psea_mc_mod$trainingData$.outcome =="Equidae",1/subset(w2,.outcome=="Equidae")$N,
ifelse(psea_mc_mod$trainingData$.outcome =="Hominidae",1/subset(w2,.outcome=="Hominidae")$N,
ifelse(psea_mc_mod$trainingData$.outcome =="Suidae",1/subset(w2,.outcome=="Suidae")$N,
       1/subset(w2,.outcome=="Phyllostomidae")$N))))),
                 preProc=c("center","scale"), method="ranger")
# ,method="glmnet", family = "multinomial", type.multinomial = "grouped"

# ctdd_mc_mod$trainingData %>% filter(!complete.cases(.))
summary(mc)
beep(sound=3)
