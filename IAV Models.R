# source("functions.R")
source("MegaLibrary.R")

reference_table <- read.csv("RefTable.csv")

# for(i in 1:length(reference_table$UID)){
#   if(reference_table[i,1] %in% sty$UID){
#     subs <- subset(sty,UID==reference_table[i,1])
# 
#     # reference_table[i,12] <- ifelse(reference_table[i,12]!=subs$Protein|is.na(reference_table[i,12]),subs$Protein,reference_table[i,12])
#     # reference_table[i,13] <- ifelse(reference_table[i,13]!=subs$Length|is.na(reference_table[i,13]),subs$Length,reference_table[i,13])
#     reference_table[i,14] <- ifelse(reference_table[i,14]!=subs$Subtype|is.na(reference_table[i,14]),subs$Subtype,reference_table[i,14])
#   }
# }
# write.csv(reference_table,"RefTable.csv",row.names = F)

SEG <- c("01PB2","02PB1","03PA","04HA","05NP","06NA","07MP","08NS")
PRT <- c("01PB2","02PB1","03PA","04HA","05NP","06NA","07M1","08NS1","07M2","08NEP")

######

mam <- subset(reference_table,Class=="Mammalia") %>%
  group_by(Order,Family,Subtype) %>% summarise(N=length(Subtype)) %>%
  subset(Family %in% c("Hominidae","Suidae","Canidae","Equidae","Phyllostomidae"))

ggplot(subset(mam,N>19),aes(Family,log10(N),fill=Subtype))+geom_col(position="dodge")+
  geom_text(aes(label=Subtype,group=Subtype),position = position_dodge(width = .9),size=3)
# ggsave("Mammal Subtypes.pdf",width=13,height=8)

######
# ssnl <- read_xlsx("Human Seasonals.xlsx",sheet="Sheet5")
# ggplot(ssnl,aes(Release_Date,WGS,colour=Genotype))+geom_point(alpha=0.7)

######
mod_cls <- read.csv("Sequences/Clusters.csv")
ggplot(mod_cls,aes(Protein,Clusters,fill=Group))+geom_col(position="dodge")+
  facet_wrap(~Identity,scales = "free")+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90))

mod_cls %>% pivot_wider(names_from = c("Group","Identity"),values_from = "Clusters") 
# %>% kbl(col.names = c("Protein","100%","95%","85%","75%","100%","95%","85%","75%"),align='c',caption="Number of representative sequences from each group") %>% kable_styling(full_width = F) %>% add_header_above(c(" "=1,"Avian"=4,"Mammalian"=4))

# a minimum coverage (option -c [0,1], which is defined by the number of aligned residue pairs divided by either the maximum of the length of query/centre and target/non-centre sequences alnRes/max(qLen,tLen) (default mode, --cov-mode 0), or by the length of the target/non-centre sequence alnRes/tLen (--cov-mode 1), or by the length of the query/centre alnRes/qLen (--cov-mode 2).

# a minimum sequence identity (--min-seq-id [0,1]) with option --alignment-mode 3 defined as the number of identical aligned residues divided by the number of aligned columns including internal gap columns, or, by default, defined by a highly correlated measure, the equivalent similarity score of the local alignment (including gap penalties) divided by the maximum of the lengths of the two locally aligned sequence segments. The score per residue equivalent to a certain sequence identity is obtained by a linear regression using thousands of local alignments as training set.
######

### Sequence Diversity

haBifas <- seqinr::read.fasta("Sequences/Avian/04HA/HA_cluster_95.fasta")
haBi <- data.frame("Base"=1:567,"Entropy"=Bios2cor::entropy(haBifas),row.names = NULL)
haMafas <- seqinr::read.fasta("Sequences/Mammals/04HA/HA_cluster_95.fasta")
haMa <- data.frame("Base"=1:566,"Entropy"=Bios2cor::entropy(haMafas),row.names = NULL)

ggplot(haBi,aes(Base,Entropy))+geom_line(colour="indianred",alpha=0.8)+geom_line(data=haMa,colour="slateblue",alpha=0.8)+
  geom_label(aes(300,0.58,label="Avian\nH = 93.13"),colour="indianred")+
  geom_label(aes(400,0.4,label="Mammalian\nH = 90.30"),colour="slateblue")+labs(title="Sitewise Diversity in Avian and Mammalian Haemagglutinin Clusters")

### Lolliplots
HASNPcon <- c(144,467)
sample.HAcon <- GRanges("HA", IRanges(HASNPcon, width=1, names=c("Gly 144 Asp","Arg 467 Ser")))
features.HA <- GRanges("HA", IRanges(c(1,1,17,52,276,330,347,514),width=c(565,16,51,223,53,16,166,20),names=c("HA","Signal Peptide","HA1 - Stalk","HA1 - Head","HA1 - Stalk","Fusion Peptide","HA2 - Stalk","Transmembrane Domain")))
features.HA$fill <- c("ivory","tan","plum3","purple3","plum3","firebrick2","lightslateblue","sienna1")
features.HA$height <- 0.08
sample.HAcon$color <- rep("dodgerblue2",2)
sample.HAcon$shape <- rep("circle",2)
trackViewer::lolliplot(sample.HAcon,features.HA,legend = legend)
grid.text("HA Segment", x=.5, y=.98, just="top",gp=gpar(cex=1.5, fontface="bold"))

######

num_clus <- read.csv("Sequences/Clustering Tests.csv")
# subset(num_clus,Identity!=100 & Group=="Avian")[,-1] %>% pivot_wider(names_from = c("Identity","Coverage"), values_from = "Count") %>% cbind(Total=c(15507,15507,15447,15447,15441,15441,15380,15380,14682,14682,15442,15442,15335,15335,15124,15124,15228,15228,15293,15293)) %>%  kbl(col.names = c("Protein","Mode",rep(c("c 50","c 70","c 80"),5),"Total"), caption="Proteins of Avian Viruses") %>% kable_styling(full_width = F) %>% add_header_above(c(" "=1," "=1,"75% Identity"=3,"85% Identity"=3,"90% Identity"=3,"95% Identity"=3,"99% Identity"=3," "=1)) %>% column_spec(12:14,background = "palegreen") %>% collapse_rows(c(1:18))

# subset(num_clus,Identity!=100 & Group=="Mammal")[,-1] %>% pivot_wider(names_from = c("Identity","Coverage"), values_from = "Count") %>% cbind(Total=c(21880,21880,21829,21829,21856,21856,22106,22106,21865,21865,21845,21845,21791,21791,21770,21770,21760,21760,21737,21737)) %>% arrange(Protein) %>%  kbl(col.names = c("Protein","Mode",rep(c("c 50","c 70","c 80"),5),"Total"), caption="Proteins of Mammalian Viruses") %>% kable_styling(full_width = F) %>% add_header_above(c(" "=1," "=1,"75% Identity"=3,"85% Identity"=3,"90% Identity"=3,"95% Identity"=3,"99% Identity"=3," "=1)) %>% column_spec(12:14,background = "palegreen") %>% collapse_rows(c(1:18))

clu <- list.files("feats",pattern = "pseaac.csv",full.names = T,recursive = T)
cls <- data.frame()
for(i in 1:length(clu)){
  clts <- data.frame(ptAcc=read.csv(clu[i])[,1] %>% str_split_i(":",1),
                     prt=clu[i] %>% str_split_i("\\/",2) %>% str_split_i("_",1))
  cls <- rbind(cls,clts)
}

cls2 <- left_join(cls,reference_table)
subset(cls2,Protein!="02PB1-F2"&Protein!="03PA-X"&Protein!="07M2"&Protein!="08NEP") %>% group_by(Protein,Subtype,Class) %>% summarise(N=length(UID)) %>% pivot_wider(names_from = Protein,values_from = N) %>% 
  # write.csv("Cluster Selection.csv",row.names = F)
  print(n=130) %>% kbl() %>% kable_styling(full_width = F)

cls3 <- subset(cls2,Protein!="02PB1-F2"&Protein!="03PA-X"&Protein!="07M2"&Protein!="08NEP") %>% group_by(Subtype,Class,Family,Superorder,Infraclass) %>% summarise(N=length(UID))

ggplot(subset(cls3,Class=="Mammalia"),aes(Family,N,fill=Subtype))+
  geom_col(position="dodge")+geom_text(aes(y=N+20,label=Subtype),check_overlap = T,position = position_dodge(width = .9),size=2.5)+
  theme(legend.position = "none")

ggplot(subset(cls3,Class=="Aves"),aes(Family,N,fill=Subtype))+geom_col(position="dodge")+geom_text(aes(y=N*.5,label=Subtype),check_overlap = T,position = position_dodge(width = .9),size=2.5)+
  theme(legend.position = "none")+facet_wrap(Infraclass~Superorder, scales = "free")

ggplot(subset(cls3,Superorder=="Galloanserae"),aes(Family,N,fill=Subtype))+geom_col(position="dodge")+geom_text(aes(y=N*.5,label=Subtype),check_overlap = T,position = position_dodge(width = .9),size=2.5)+
  theme(legend.position = "none")
ggplot(subset(cls3,Superorder=="Neoaves"),aes(Family,N,fill=Subtype))+geom_col(position="dodge")+geom_text(aes(y=N*.5,label=Subtype),check_overlap = T,position = position_dodge(width = .9),size=2.5)+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))

###### Feature Models 

fea <- list.files("feats",pattern = ".csv",full.names = T,recursive = T)
# fea <- fea[1:18]
for(i in 1:length(fea)){
  print(fea[i])
  tmp_ft <- read.csv(fea[i])
  tmp_ft$ptAcc <- tmp_ft$X %>% str_split_i(":",1)
  tmp_ft <- right_join(reference_table,tmp_ft)
  tmp_ft <- subset(tmp_ft, str_to_lower(str_split_i(tmp_ft$Sample,"\\/",2))!="reassortant")
  #print(table(tmp_ft$Class))
  fold_indices <- createMultiFolds(tmp_ft$UID, k = 5, times = 1)
  preds <- tmp_ft[,c(-1:-15)] %>% remove_constant %>% names
  # tmp_ft[,c(-2:-16)] %>% remove_constant %>% pivot_longer(cols = -1, names_to = "Ft",values_to = "Val") %>% subset(Val==0) %>% print()
  # subset(tmp_ft,is.na(Class))[,1:15] %>% write.csv(paste0(substr(fea[i],1,10),"Lost.csv"))

  mo1 <- suppressWarnings(train(x = tmp_ft %>% select(all_of(preds)),y = tmp_ft %>% pull(Class),
                         method="ranger",metric="ROC",
                         preProc=c("center","scale"),num.trees=1000,importance="impurity",
                         weights=ifelse(tmp_ft$Class =="Aves",
                                        (1/table(tmp_ft$Class)[1]) * 0.5,
                                        (1/table(tmp_ft$Class)[2]) * 0.5),
                         trControl = trainControl(method = "repeatedcv",number = 5,
                                                  repeats = 1,index = fold_indices,
                                                  classProbs = TRUE,savePredictions = TRUE,
                         summaryFunction=twoClassSummary),
                         tuneGrid = expand.grid(.splitrule = "gini",
                                                .min.node.size = seq(from = 5, to = 45, length = 3),
                                                .mtry = round(sqrt(length(preds))))))

  mo_nm <- fea[i] %>% str_split_i("\\/",2) %>% str_split("\\_")
  saveRDS(mo1, file=paste0("Models/",mo_nm[[1]][1],"_",mo_nm[[1]][3],"_",
                           mo_nm[[1]][4],"_",mo_nm[[1]][5],"_",str_replace_all(mo_nm[[1]][8], ".csv",""),".rds"))
}

# Model Stacks

mods <- list.files("Models",pattern = ".rds",full.names = T,recursive = F)
mods <- unique(str_split_i(mods,"_",1))
for(i in 1:length(mods)){
  print(str_to_upper(str_split_i(mods[i],"\\/",2)))
  tmer_mod <- readRDS(paste0(mods[i],"_2mer.rds"))
  ctdc_mod <- readRDS(paste0(mods[i],"_ctdc.rds"))
  ctdd_mod <- readRDS(paste0(mods[i],"_ctdd.rds"))
  ctdt_mod <- readRDS(paste0(mods[i],"_ctdt.rds"))
  ctri_mod <- readRDS(paste0(mods[i],"_ctriad.rds"))
  psea_mod <- readRDS(paste0(mods[i],"_pseaac.rds"))

  # caretEnsemble(c("two_mer"=tmer_mod,"CTDc"=ctdc_mod,"CTDd"=ctdd_mod,"CTDt"=ctdt_mod,"triad"=ctri_mod,"pseaac"=psea_mod),preProc = c("scale"),metric = "ROC") %>% saveRDS(paste0("Comp/Ensemble/",mods[i],".rds"))
  caretStack(c("two_mer"=tmer_mod,"CTDc"=ctdc_mod,"CTDd"=ctdd_mod,"CTDt"=ctdt_mod,"triad"=ctri_mod,"pseaac"=psea_mod),preProc=c("center","scale"),metric = "ROC",method="glmnet") %>% saveRDS(paste0("Comp/Stack/",mods[i],".rds"))
}

modc <- list.files("Comp/Stack/Models",pattern = ".rds",full.names = T,recursive = F)
for(i in 1:length(modc)){
  # ModMeth <- str_split_i(modc[i],"\\/",2)
  ModParam <- str_split_i(modc[i],"\\/",4) %>% str_replace_all(".rds","")
  MoDo <- paste0(str_to_upper(ModParam),"_complete")
  assign(MoDo, readRDS(modc[i]))
}

Stack_ha_95_c80_cov0
summary(Stack_ha_95_c80_cov0)

###### Cluster Tests
stack_res <- read_xlsx("Comp/HA test results.xlsx",sheet="Stacks")
ensemb_res <- read_xlsx("Comp/HA test results.xlsx",sheet="Ensemble")

ggplot(subset(stack_res,Model!="ensemble"),aes(paste(Coverage,Mode),Importance,colour=Model))+
  geom_point()
ggplot(stack_res,aes(paste(Coverage,Mode),ymin=ROC-sdROC,y=ROC,ymax=ROC+sdROC,fill=Model))+
  geom_crossbar(position="dodge")

stack_res %>% kbl() %>% kable_styling(full_width = F) %>% collapse_rows(4:5)
ensemb_res %>% kbl() %>% kable_styling(full_width = F) %>% collapse_rows(4:5)

###### Model Outputs
# protes <- c("pb2","pb1","pa","ha","np","na","m1","ns1")

mos_out <- rbind(data.frame(protein="07M1",feature=str_split_i(names(summary(M1_complete)$imp),"\\.",1),imp=summary(M1_complete)$imp,summary(M1_complete)$results[2:7],row.names = NULL),rbind(data.frame(protein="04HA",feature=str_split_i(names(summary(HA_complete)$imp),"\\.",1),imp=summary(HA_complete)$imp,summary(HA_complete)$results[2:7],row.names = NULL),rbind(data.frame(protein="06NA",feature=str_split_i(names(summary(NA_complete)$imp),"\\.",1),imp=summary(NA_complete)$imp,summary(NA_complete)$results[2:7],row.names = NULL),rbind(data.frame(protein="03PA",feature=str_split_i(names(summary(PA_complete)$imp),"\\.",1),imp=summary(PA_complete)$imp,summary(PA_complete)$results[2:7],row.names = NULL),rbind(data.frame(protein="08NS1",feature=str_split_i(names(summary(NS1_complete)$imp),"\\.",1),imp=summary(NS1_complete)$imp,summary(NS1_complete)$results[2:7],row.names = NULL),rbind(data.frame(protein="02PB1",feature=str_split_i(names(summary(PB1_complete)$imp),"\\.",1),imp=summary(PB1_complete)$imp,summary(PB1_complete)$results[2:7],row.names = NULL),rbind(data.frame(protein="05NP",feature=str_split_i(names(summary(NP_complete)$imp),"\\.",1),imp=summary(NP_complete)$imp,summary(NP_complete)$results[2:7],row.names = NULL),data.frame(protein="01PB2",feature=str_split_i(names(summary(PB2_complete)$imp),"\\.",1),imp=summary(PB2_complete)$imp,summary(PB2_complete)$results[2:7],row.names = NULL))))))))

ggplot(mos_out,aes(feature,ymin=value-sd,y=value,ymax=value+sd,colour=imp))+
  geom_pointrange()+facet_wrap(~protein,ncol=4)+scale_colour_gradient(low="grey",high="red")+theme(legend.position = "bottom",axis.text.x = element_text(angle = 45))
