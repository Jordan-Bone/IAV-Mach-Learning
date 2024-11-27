knitr::opts_chunk$set(echo=FALSE,warning=FALSE,message=FALSE,root.dir= "/Users/jordanbone/Documents/IAV.AI")
source("functions.R")
# source("MegaLibrary.R")
# taxize_options(ncbi_sleep = 1)
reference_table <- read.csv("Reference Table 2.csv")

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

mod_cls <- read_excel("Sequences/Clusters.xlsx")
ggplot(mod_cls,aes(Protein,Clusters,fill=Group))+geom_col(position="dodge")+
  facet_wrap(~Identity,scales = "free")+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90))

mod_cls %>% pivot_wider(names_from = c("Group","Identity"),values_from = "Clusters") %>% kbl(col.names = c("Protein","100%","95%","85%","75%","100%","95%","85%","75%"),align='c',caption="Number of representative sequences from each group") %>% kable_styling(full_width = F) %>% add_header_above(c(" "=1,"Avian"=4,"Mammalian"=4))

######

### Sequence Diversity

haBifas <- read.fasta("Sequences/Avian/ha_cluster_95_rep_seq.fasta")
haMafas <- read.fasta("Sequences/Mammals/ha_cluster_rep_95.fasta")
haBi <- data.frame("Base"=1:567,"Entropy"=Bios2cor::entropy(haBifas),row.names = NULL)
haMa <- data.frame("Base"=1:566,"Entropy"=Bios2cor::entropy(haMafas),row.names = NULL)

ggplot(haBi,aes(Base,Entropy))+geom_line(colour="indianred",alpha=0.8)+geom_line(data=haMa,colour="slateblue",alpha=0.8)+
  geom_label(aes(350,0.955,label="Avian\nH = 93.13"),colour="indianred")+
  geom_label(aes(300,0.85,label="Mammalian\nH = 90.30"),colour="slateblue")+labs(title="Sitewise Diversity in Avian and Mammalian Haemagglutinin Clusters")

### Lolliplots
library(trackViewer)
HASNPcon <- c(144,467)
sample.HAcon <- GRanges("HA", IRanges(HASNPcon, width=1, names=c("Gly 144 Asp","Arg 467 Ser")))
features.HA <- GRanges("HA", IRanges(c(1,1,17,52,276,330,347,514),width=c(565,16,51,223,53,16,166,20),names=c("HA","Signal Peptide","HA1 - Stalk","HA1 - Head","HA1 - Stalk","Fusion Peptide","HA2 - Stalk","Transmembrane Domain")))
features.HA$fill <- c("ivory","tan","plum3","purple3","plum3","firebrick2","lightslateblue","sienna1")
features.HA$height <- 0.08
sample.HAcon$color <- rep("dodgerblue2",2)
sample.HAcon$shape <- rep("circle",2)
trackViewer::lolliplot(sample.HAcon,features.HA,legend = legend)
grid.text("HA Segment", x=.5, y=.98, just="top",gp=gpar(cex=1.5, fontface="bold"))
