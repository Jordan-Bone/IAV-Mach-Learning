setwd("/Users/jordanbone/Documents/GitHub/IAV-Mach-Learning")
source("MegaLibrary.R")

reference_table <- read.csv("RefTable.csv")
rep_uids <- read.csv("USED UIDS.txt")

clu <- list.files("feats/combo",pattern = ".csv",full.names = T,recursive = T)
cls <- data.frame()
for(i in 1:length(clu)){
  clts <- data.frame(ptAcc=read.csv(clu[i])[,1] %>% str_split_i(":",1),
                     prt=clu[i] %>% str_split_i("\\/",2) %>% str_split_i("_",1))
  cls <- rbind(cls,clts)
}

cls2 <- left_join(cls,reference_table, relationship = "many-to-many")

####
usrs <- read_xlsx("/Users/jordanbone/Documents/IAV.AI/Used Seqs.xlsx")[,c(2:10)]
us2 <- usrs %>% pivot_longer(2:9, names_to = "Protein", values_to = "ptAcc2")
reference_table <- left_join(reference_table,us2,by=c("UID","Protein"))

for(i in 1:length(reference_table$UID)){
#   if(reference_table[i,12]!="08NEP" & reference_table[i,12]!="07M2" & reference_table[i,12]!="02PB1-F2" & reference_table[i,12]!="03PA-X"){
#     if(reference_table[i,1] %in% us2$UID){
#       subs <- subset(us2,UID %in% reference_table[i,1] & Protein == reference_table[i,12])
#       reference_table[i,11] <- ifelse(reference_table[i,11]!=subs$ptAcc,subs$ptAcc,reference_table[i,11])
      #     # reference_table[i,9] <- ifelse(reference_table[i,9]!=subs$ntAcc|is.na(reference_table[i,9]),subs$ntAcc,reference_table[i,9])
      #     # reference_table[i,14] <- ifelse(reference_table[i,14]!=subs$Subtype|is.na(reference_table[i,14]),subs$Subtype,reference_table[i,14])
#     }}
  reference_table[i,11] <- ifelse(reference_table[i,11]!=reference_table[i,16] &!is.na(reference_table[i,16]),reference_table[i,16],reference_table[i,11])
  per.done(i)
}
beep(2)
write.csv(reference_table[,-16],"RefTable.csv",row.names = F)

SEG <- c("01PB2","02PB1","03PA","04HA","05NP","06NA","07MP","08NS")
PRT <- c("01PB2","02PB1","03PA","04HA","05NP","06NA","07M1","08NS1","07M2","08NEP")
FET <- c("pseaac","2mer","ctdc","ctdd","ctdt","ctriad")
######

mam <- subset(reference_table,Class=="Mammalia") %>%
  group_by(Order,Family,Subtype) %>% summarise(N=length(Subtype)) %>%
  subset(Family %in% c("Hominidae","Suidae","Canidae","Equidae","Phyllostomidae"))

ggplot(subset(mam,N>19),aes(Family,log10(N),fill=Subtype))+geom_col(position="dodge")+
  geom_text(aes(label=Subtype,group=Subtype),position = position_dodge(width = .9),size=3)
# ggsave("Mammal Subtypes.pdf",width=13,height=8)
######

fetc <- list.files("feats",pattern = ".csv",full.names = T,recursive = F)
for(i in 1:length(fetc)){
  FeatParam <- str_split_i(fetc[i],"\\/",2) %>% str_replace_all(".csv","")
  SegFeat <- paste0(str_split_i(FeatParam,"_",1),"_",str_split_i(FeatParam,"_",8))
  assign(SegFeat, read.csv(fetc[i]))
}

ha_ctdc$Seg <- "HA"
ha_ctdc$Feat <- "CTDc"

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

# haBifas <- seqinr::read.fasta("Sequences/Avian/ha_cluster_95_c80_cov0_rep_seq.fasta")
# haBi <- data.frame("Base"=1:567,"Entropy"=Bios2cor::entropy(haBifas),row.names = NULL)
# haMafas <- seqinr::read.fasta("Sequences/Mammals/04HA/HA_cluster_95.fasta")
# haMa <- data.frame("Base"=1:566,"Entropy"=Bios2cor::entropy(haMafas),row.names = NULL)
# 
# ggplot(haBi,aes(Base,Entropy))+geom_line(colour="indianred",alpha=0.8)+geom_line(data=haMa,colour="slateblue",alpha=0.8)+
#   geom_label(aes(300,0.58,label="Avian\nH = 93.13"),colour="indianred")+
#   geom_label(aes(400,0.4,label="Mammalian\nH = 90.30"),colour="slateblue")+labs(title="Sitewise Diversity in Avian and Mammalian Haemagglutinin Clusters")

### Lolliplots
# HASNPcon <- c(144,467)
# sample.HAcon <- GRanges("HA", IRanges(HASNPcon, width=1, names=c("Gly 144 Asp","Arg 467 Ser")))
# features.HA <- GRanges("HA", IRanges(c(1,1,17,52,276,330,347,514),width=c(565,16,51,223,53,16,166,20),names=c("HA","Signal Peptide","HA1 - Stalk","HA1 - Head","HA1 - Stalk","Fusion Peptide","HA2 - Stalk","Transmembrane Domain")))
# features.HA$fill <- c("ivory","tan","plum3","purple3","plum3","firebrick2","lightslateblue","sienna1")
# features.HA$height <- 0.08
# sample.HAcon$color <- rep("dodgerblue2",2)
# sample.HAcon$shape <- rep("circle",2)
# trackViewer::lolliplot(sample.HAcon,features.HA,legend = legend)
# grid.text("HA Segment", x=.5, y=.98, just="top",gp=gpar(cex=1.5, fontface="bold"))

######

# num_clus <- read.csv("Sequences/Clustering Tests.csv")
# subset(num_clus,Identity!=100 & Group=="Avian")[,-1] %>% pivot_wider(names_from = c("Identity","Coverage"), values_from = "Count") %>% cbind(Total=c(15507,15507,15447,15447,15441,15441,15380,15380,14682,14682,15442,15442,15335,15335,15124,15124,15228,15228,15293,15293)) %>%  kbl(col.names = c("Protein","Mode",rep(c("c 50","c 70","c 80"),5),"Total"), caption="Proteins of Avian Viruses") %>% kable_styling(full_width = F) %>% add_header_above(c(" "=1," "=1,"75% Identity"=3,"85% Identity"=3,"90% Identity"=3,"95% Identity"=3,"99% Identity"=3," "=1)) %>% column_spec(12:14,background = "palegreen") %>% collapse_rows(c(1:18))

# subset(num_clus,Identity!=100 & Group=="Mammal")[,-1] %>% pivot_wider(names_from = c("Identity","Coverage"), values_from = "Count") %>% cbind(Total=c(21880,21880,21829,21829,21856,21856,22106,22106,21865,21865,21845,21845,21791,21791,21770,21770,21760,21760,21737,21737)) %>% arrange(Protein) %>%  kbl(col.names = c("Protein","Mode",rep(c("c 50","c 70","c 80"),5),"Total"), caption="Proteins of Mammalian Viruses") %>% kable_styling(full_width = F) %>% add_header_above(c(" "=1," "=1,"75% Identity"=3,"85% Identity"=3,"90% Identity"=3,"95% Identity"=3,"99% Identity"=3," "=1)) %>% column_spec(12:14,background = "palegreen") %>% collapse_rows(c(1:18))

# subset(cls2,Protein!="02PB1-F2"&Protein!="03PA-X"&Protein!="07M2"&Protein!="08NEP") %>% group_by(Protein,Subtype,Class) %>% summarise(N=length(UID)) %>% pivot_wider(names_from = Protein,values_from = N) %>% 
  # write.csv("Cluster Selection.csv",row.names = F)
  # kbl() %>% kable_styling(full_width = F)

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

cls3 %>% group_by(Subtype) %>% summarise(sn=sum(N)) %>% subset(sn>25)

ggplot(cls3,aes(x=reorder(Subtype,-N),y=N,fill=ifelse(Class=="Aves",Class,Family)))+geom_col(position="dodge")+theme(legend.position = "bottom",axis.text.x = element_text(angle = 90))+labs(x="Subtype",fill="Group")+geom_vline(xintercept = 58.5, alpha=0.4, linetype="longdash")
ggsave("Subtype Frequencies.pdf",width=33,height=19,units = "cm")

ggplot(subset(cls3,N>=50),aes(x=reorder(Subtype,-N),y=N,fill=ifelse(Class=="Aves",Class,Family)))+geom_col(position="dodge")+theme(legend.position = "bottom",axis.text.x = element_text(angle = 90))+labs(x="Subtype",fill="Group")+geom_vline(xintercept = 58.5, alpha=0.4, linetype="longdash")
ggsave("Subtype Frequencies zoom.pdf",width=33,height=19,units = "cm")

cls3$Classification <- ifelse(cls3$Class=="Mammalia",cls3$Family,str_replace_all(cls3$Superorder,"Aves","Galloanserae"))
cls4 <- cls3 %>% group_by(Subtype,Classification) %>% summarise(N=sum(N))
ggplot(cls4,aes(x=reorder(Subtype,-N),y=N,fill=Classification))+geom_col()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 90))

###### Cluster Tests
stack_res <- read_xlsx("Comp/HA test results.xlsx",sheet="Stacks")
ensemb_res <- read_xlsx("Comp/HA test results.xlsx",sheet="Ensemble")

ggplot(subset(stack_res,Model!="ensemble"),aes(paste(Coverage,Mode),Importance,colour=Model))+
  geom_point()
ggplot(stack_res,aes(paste(Coverage,Mode),ymin=ROC-sdROC,y=ROC,ymax=ROC+sdROC,fill=Model))+
  geom_crossbar(position="dodge")

# stack_res %>% kbl() %>% kable_styling(full_width = F) %>% collapse_rows(4:5)
# ensemb_res %>% kbl() %>% kable_styling(full_width = F) %>% collapse_rows(4:5)

###### Model Outputs

# mos_out <- rbind(data.frame(protein="07M1",feature=str_split_i(names(summary(M1_complete)$imp),"\\.",1),imp=summary(M1_complete)$imp,summary(M1_complete)$results[2:7],row.names = NULL),rbind(data.frame(protein="04HA",feature=str_split_i(names(summary(HA_complete)$imp),"\\.",1),imp=summary(HA_complete)$imp,summary(HA_complete)$results[2:7],row.names = NULL),rbind(data.frame(protein="06NA",feature=str_split_i(names(summary(NA_complete)$imp),"\\.",1),imp=summary(NA_complete)$imp,summary(NA_complete)$results[2:7],row.names = NULL),rbind(data.frame(protein="03PA",feature=str_split_i(names(summary(PA_complete)$imp),"\\.",1),imp=summary(PA_complete)$imp,summary(PA_complete)$results[2:7],row.names = NULL),rbind(data.frame(protein="08NS1",feature=str_split_i(names(summary(NS1_complete)$imp),"\\.",1),imp=summary(NS1_complete)$imp,summary(NS1_complete)$results[2:7],row.names = NULL),rbind(data.frame(protein="02PB1",feature=str_split_i(names(summary(PB1_complete)$imp),"\\.",1),imp=summary(PB1_complete)$imp,summary(PB1_complete)$results[2:7],row.names = NULL),rbind(data.frame(protein="05NP",feature=str_split_i(names(summary(NP_complete)$imp),"\\.",1),imp=summary(NP_complete)$imp,summary(NP_complete)$results[2:7],row.names = NULL),data.frame(protein="01PB2",feature=str_split_i(names(summary(PB2_complete)$imp),"\\.",1),imp=summary(PB2_complete)$imp,summary(PB2_complete)$results[2:7],row.names = NULL))))))))

# data.frame(protein="01PB2",feature=str_split_i(names(summary(PB2_complete)$imp),"\\.",1),imp=summary(PB2_complete)$imp,summary(PB2_complete)$results[2:7],row.names = NULL)

mos_out <- read.csv("Comp/Model Outputs.csv")
ggplot(mos_out,aes(feature,ymin=value-sd,y=value,ymax=value+sd,colour=imp))+
  geom_pointrange()+facet_wrap(~protein,ncol=4)+scale_colour_gradient(low="grey",high="red")+theme(legend.position = "bottom",axis.text.x = element_text(angle = 45))
#####

# Subtype Distribution
sbd <- cls2 %>% group_by(Order,Class,Superorder,Subtype) %>% summarise(N=length(UID))
sbd$Classification <- ifelse(sbd$Class=="Mammalia",sbd$Order,sbd$Superorder)
# sbd %>% group_by(Classification,Subtype) %>% summarise(N=sum(N)) %>% pivot_wider(names_from = "Classification",values_from = N) %>% write.csv("Subtype Distribution.csv",row.names = F)

reference_table[,c(1,14,15,11,12)] %>% subset(Protein!="02PB1-F2"&Protein!="03PA-X") %>% group_by(UID,Subtype,Classification,Protein) %>% summarise(N=length(ptAcc)) %>% pivot_wider(names_from = Protein, values_from = N)

# dupes <- reference_table[,c(1,14,15,11,12)] %>% subset(Protein!="02PB1-F2"&Protein!="03PA-X") |> dplyr::summarise(n = dplyr::n(), .by = c(UID, Subtype,Classification, Protein)) |> dplyr::filter(n > 1L)

ax <- reference_table[,c(1,14,15,11,12)] %>% subset(Protein!="02PB1-F2"&Protein!="03PA-X"&Protein!="08NEP"&Protein!="07M2") 
ax$Rep <- ifelse(ax$ptAcc %in% cls$ptAcc,"Rep","Not")

badax <- ax |> dplyr::summarise(n = dplyr::n(), .by = c(UID, Subtype, Classification, Rep, Protein)) |> dplyr::filter(n > 1L)

for(i in 1:length(ax$UID)){
  ax[i,4] <- ifelse(ax[i,6]=="Rep",cell_spec(ax[i,4], background = "forestgreen"),ax[i,4])
}
# ax$ptAcc <- cell_spec(ax$ptAcc, background = ifelse(ax$Rep=="Rep","forestgreen","white"))

subset(ax,!UID %in% badax$UID)[,-6] %>% group_by(UID,Subtype,Classification) %>% pivot_wider(names_from = Protein, values_from = ptAcc) %>% kbl(escape=F) %>% kable_styling() %>% save_kable(file = "table1.html")
# %>% write.csv("Full data.csv",row.names = F)

mod_stats <- read.csv("Models/Mod Stats.csv")
subset(mod_stats,FeatSet=="ctdd") %>% pivot_longer(cols=4:5,names_to = "Class",values_to = "N") %>% ggplot(aes(Holdout,N,colour=Class))+geom_point(alpha=0.6)
subset(mod_stats,Avian+Mammal!=2325)

mostat <- read.csv("MoStat.csv")
ggplot(subset(mostat,F=="ctdd"),aes(str_replace_all(X,"ex",""),N))+geom_point(alpha=0.6)+facet_wrap(~P,scales="free",nrow=2)+theme(legend.position = "bottom",axis.text.x = element_text(angle = 90))+labs(x="Subtype")
ggsave("Num Holdouts.pdf",width=50,height=19,units = "cm")

ggplot(subset(mostat,F=="ctdd"&P=="04ha"),aes(str_replace_all(X,"ex",""),N,colour=N==max(N)))+geom_point(alpha=0.9)+theme(legend.position = "bottom",axis.text.x = element_text(angle = 90))+labs(x="Subtype")
ggsave("HA Holdouts.pdf",width=30,height=19,units = "cm")
ha2324 <- subset(mostat,P=="ha"&N==2324&F=="2mer")$X

ggplot(subset(mostat,F=="ctdd"&P=="ha"&N<2300),aes(str_replace_all(X,"ex",""),y=max(mostat$N)-N))+geom_point(alpha=0.9)+labs(x="Excluded Subtype",y="Number of Sequences")
