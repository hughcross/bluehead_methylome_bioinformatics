############################################################
# Methylome data analysis - Coupling                       #
############################################################

# Load libraries
library(cowplot)
library(Gviz)
library(GenomicRanges)

# Change folder

setwd("~/Syncplicity Folders/BHW_ToddOrtegaLiu2018/Figures/Figure_6_epigenome/PCA methylome")

############################################################
# PCA Methylome                                            #
############################################################

probes <- read.csv("probes_10k_more100.txt",sep= "\t")
probes[1:13] <- c()

probes_t <- as.data.frame(t(probes))

samples_name <- c("TBK14_31G","TBK14_36G","TBK14_63G","TBK14_72G","TBK14_7G","TBK14_82G",
                  "TBK14_11G","TBK14_14G","TBZ13_17G","TBK14_74G","TBZ13_6G","TBK14_50G",
                  "TBK14_17G","TBK14_20G","TPZ16_112G","TPZ16_B109G","TPZ16_B110G")

samples_stage <- c(rep("CF",3),rep("S2",3),rep("S3",3),rep("S4",2),rep("S5a",1),
                   rep("S5b",2),rep("TP",3))


samples_df <- data.frame("name" = samples_name, "stage" = samples_stage)
rownames(probes_t) <- samples_df$name
pca <- prcomp(probes_t)

summary(pca)

PCA_values <- as.data.frame(pca$x)
PCA_values$group <- c(rep("CF",3),rep("S2",3),rep("S3",3), rep("S4",2),
                      rep("S5a",1),rep("S5b",2), rep("TP",3))

cbPalette <- c("#8f0000","#ff8000","#ffd940","#b3ebff",
               "#4dd9ff","#1c79b2","#00048b",3)

plot_PCA <- ggplot(PCA_values, aes(x = PC1, y = PC2))+
  geom_point(aes(color = PCA_values$group),size=5,alpha=0.9) +
  scale_color_manual(values=cbPalette) +
  labs(x = "PC1: 64.08% variance", y = "PC2: 10.72% variance") +
  theme(axis.text = element_text(size=14),axis.title = element_text(size=16))

plot_PCA

# Clean environment
rm(list=ls())

############################################################
# Coupling                                                 #
############################################################

setwd("~/Syncplicity Folders/BHW_ToddOrtegaLiu2018/Figures/Figure_6_epigenome/Coupling")

meth_data <- read.csv2("meth_probes_no3-4-5a.txt", sep = "\t")
RNA_data <- read.csv2("RNA_probes.txt", sep = "\t")

# CF
meth_data_CF <- data.frame("ID"= sub('_upstream', '', meth_data$Probe),
                           "meth" = as.numeric(as.character(meth_data$CF.met_mod,as.numeric)))

RNA_data_CF <- data.frame("ID"= RNA_data$Probe,
                          "rna" = as.numeric(as.character(RNA_data$CF_RNA)),stringsAsFactors = FALSE)

data_CF <- merge(meth_data_CF,RNA_data_CF, by = "ID")

quantile_CF <- quantile(data_CF$rna, seq(0, 1, 0.25))
data_CF$quantile <- findInterval(data_CF$rna, quantile_CF, all.inside = TRUE)

table(data_CF$quantile)

data_CF$quantile <- factor(data_CF$quantile)

colfunc_CF <- colorRampPalette(c("white", "#8f0000"))
col_CF <- colfunc_CF(8)[2:5]

plot_CF <- ggplot(data_CF, aes(x=quantile, y=meth,fill=quantile)) +
  geom_violin(scale = "width")+
  stat_summary(fun.y=median, geom="point", size=4, color="black") +
  theme(axis.text = element_text(size=20)) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_blank())+
  scale_fill_manual(values=col_CF) + 
  guides(fill=FALSE)

plot_CF

# S3
meth_data_S3 <- data.frame("ID"= sub('_upstream', '', meth_data$Probe),
                           "meth" = as.numeric(as.character(meth_data$S3.met,as.numeric)))

RNA_data_S3 <- data.frame("ID"= RNA_data$Probe,
                          "rna" = as.numeric(as.character(RNA_data$S3_RNA)),stringsAsFactors = FALSE)

data_S3 <- merge(meth_data_S3,RNA_data_S3, by = "ID")

quantile_S3 <- quantile(data_S3$rna, seq(0, 1, 0.25))
data_S3$quantile <- findInterval(data_S3$rna, quantile_S3, all.inside = TRUE)

table(data_S3$quantile)

data_S3$quantile <- factor(data_S3$quantile)
str(data_S3)

colfunc_S3 <- colorRampPalette(c("white", "#ffd940"))
col_S3 <- colfunc_S3(8)[2:5]

plot_S3 <- ggplot(data_S3, aes(x=quantile, y=meth,fill=quantile)) +
  geom_violin(scale = "width")+
  stat_summary(fun.y=median, geom="point", size=4, color="black") +
  theme(axis.text = element_text(size=20)) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_blank())+
  scale_fill_manual(values=col_S3) + 
  guides(fill=FALSE)

plot_S3

# S5b
meth_data_S5b <- data.frame("ID"= sub('_upstream', '', meth_data$Probe),
                            "meth" = as.numeric(as.character(meth_data$S5b.met,as.numeric)))

RNA_data_S5b <- data.frame("ID"= RNA_data$Probe,
                           "rna" = as.numeric(as.character(RNA_data$S5b_RNA)),stringsAsFactors = FALSE)

data_S5b <- merge(meth_data_S5b,RNA_data_S5b, by = "ID")

quantile_S5b <- quantile(data_S5b$rna, seq(0, 1, 0.25))
data_S5b$quantile <- findInterval(data_S5b$rna, quantile_S5b, all.inside = TRUE)

table(data_S5b$quantile)

data_S5b$quantile <- factor(data_S5b$quantile)

colfunc_S5b <- colorRampPalette(c("white", "#1D79B3"))
col_S5b <- colfunc_S5b(8)[2:5]

plot_S5b <- ggplot(data_S5b, aes(x=quantile, y=meth,fill=quantile)) +
  geom_violin(scale = "width")+
  stat_summary(fun.y=median, geom="point", size=4, color="black") +
  theme(axis.text = element_text(size=20)) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_blank())+
  scale_fill_manual(values=col_S5b) + 
  guides(fill=FALSE)

plot_S5b

# TP
meth_data_TP <- data.frame("ID"= sub('_upstream', '', meth_data$Probe),
                           "meth" = as.numeric(as.character(meth_data$TP.met,as.numeric)))

RNA_data_TP <- data.frame("ID"= RNA_data$Probe,
                          "rna" = as.numeric(as.character(RNA_data$TP_RNA)),stringsAsFactors = FALSE)

data_TP <- merge(meth_data_TP,RNA_data_TP, by = "ID")

quantile_TP <- quantile(data_TP$rna, seq(0, 1, 0.25))
data_TP$quantile <- findInterval(data_TP$rna, quantile_TP, all.inside = TRUE)

table(data_TP$quantile)

data_TP$quantile <- factor(data_TP$quantile)

colfunc_TP <- colorRampPalette(c("white", "#00048b"))
col_TP <- colfunc_TP(8)[2:5]

plot_TP <- ggplot(data_TP, aes(x=quantile, y=meth,fill=quantile)) +
  geom_violin(scale = "width")+
  stat_summary(fun.y=median, geom="point", size=4, color="black") +
  theme(axis.text = element_text(size=20)) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_blank())+
  scale_fill_manual(values=col_TP) + 
  guides(fill=FALSE)

plot_TP

plot_coupling <- plot_grid(plot_CF,plot_S3,plot_S5b,plot_TP,nrow = 1)
plot_coupling

############################################################
# Annotation tracks                                                 #
############################################################

# Change folders

setwd("~/Syncplicity Folders/BHW_ToddOrtegaLiu2018/Figures/Figure_6_epigenome/gene_tracks")

# Methylation cyp19a1a
data_cyp19a1_meth <- read.csv("cyp19a1_meth.txt",sep= "\t")

df_cyp19a1_meth <- data.frame("chr"=rep(5,nrow(data_cyp19a1_meth)),"start"=data_cyp19a1_meth$Start,
                              "end"=data_cyp19a1_meth$End,"strand"=rep("*",nrow(data_cyp19a1_meth)))

Gr_cyp19a1_meth_F <- makeGRangesFromDataFrame(df_cyp19a1_meth)
values(Gr_cyp19a1_meth_F) <- DataFrame(female = data_cyp19a1_meth$CF.met_mod)

Gr_cyp19a1_meth_M <- makeGRangesFromDataFrame(df_cyp19a1_meth)
values(Gr_cyp19a1_meth_M) <- DataFrame(male = data_cyp19a1_meth$TP.met)

cyp19a1_meth_track_F <- AnnotationTrack(Gr_cyp19a1_meth_F, name = "Female")
cyp19a1_meth_track_M <- AnnotationTrack(Gr_cyp19a1_meth_M, name = "Male")

plotTracks(list(DataTrack(Gr_cyp19a1_meth_F, name ="Female",ylim= c(0,100),col = "black"),
                DataTrack(Gr_cyp19a1_meth_M, name ="Male",ylim = c(0,100),col = "black")),
           type = c("p","histogram","g"), background.title = "darkblue",
           from = 4168500, to = 4181900)

# Expression cyp19a1a
data_cyp19a1_RNA <- read.csv("cyp19a1_RNA.txt",sep= "\t")

df_cyp19a1_RNA <- data.frame("chr"=rep(5,nrow(data_cyp19a1_RNA)),"start"=data_cyp19a1_RNA$Start,
                             "end"=data_cyp19a1_RNA$End,"strand"=rep("*",nrow(data_cyp19a1_RNA)))

Gr_cyp19a1_RNA_F <- makeGRangesFromDataFrame(df_cyp19a1_RNA)  
values(Gr_cyp19a1_RNA_F) <- DataFrame(female = data_cyp19a1_RNA$CF_RNA)

Gr_cyp19a1_RNA_M <- makeGRangesFromDataFrame(df_cyp19a1_RNA)  
values(Gr_cyp19a1_RNA_M) <- DataFrame(male = data_cyp19a1_RNA$TP_RNA)

cyp19a1_RNA_track_F <- AnnotationTrack(Gr_cyp19a1_RNA_F, name = "Female")
cyp19a1_RNA_track_M <- AnnotationTrack(Gr_cyp19a1_RNA_M, name = "Male")

plotTracks(list(DataTrack(Gr_cyp19a1_RNA_F, name ="Female",ylim= c(0,3)),
                DataTrack(Gr_cyp19a1_RNA_M, name ="Male",ylim = c(0,3))),
           type = c("g","histogram"), background.title = "darkblue",
           from = 4168500, to = 4181900)

# CpG islands and gene model cyp19a1a
geneModel <- read.csv("cyp19a1_gene_A.csv")
grtrack <- GeneRegionTrack(geneModel,name = "cyp19a")

cyp19a1_scale <- data.frame("chr"=5,"start"=4170000,"end"=4172000,"strand"="*")
Gr_cyp19a1_scale <- makeGRangesFromDataFrame(cyp19a1_scale)

cyp19a1_CpG <- data.frame("chr"=c(5,5,5),"start"=c(4171450,4173380,4176420),
                          "end"=c(4171700,4173584,4176728),
                          "strand"="*")
Gr_cyp19a1_CpG <- makeGRangesFromDataFrame(cyp19a1_CpG)

plotTracks(list(grtrack,AnnotationTrack(Gr_cyp19a1_CpG, name ="CpG")),
           type = c("g","histogram"), background.title = "darkblue",
           from = 4168500, to = 4181900) 

# Methylation dmrt1
data_dmrt1_meth <- read.csv("dmrt1_meth_def.txt",sep= "\t")
df_dmrt1_meth <- data.frame("chr"=rep(5,nrow(data_dmrt1_meth)),"start"=data_dmrt1_meth$Start,
                            "end"=data_dmrt1_meth$End,"strand"=rep("*",nrow(data_dmrt1_meth)))

Gr_dmrt1_meth_F <- makeGRangesFromDataFrame(df_dmrt1_meth)
values(Gr_dmrt1_meth_F) <- DataFrame(female = data_dmrt1_meth$CF.met_mod)

Gr_dmrt1_meth_M <- makeGRangesFromDataFrame(df_dmrt1_meth)
values(Gr_dmrt1_meth_M) <- DataFrame(male = data_dmrt1_meth$TP.met)

dmrt1_meth_track_F <- AnnotationTrack(Gr_dmrt1_meth_F, name = "Female")
dmrt1_meth_track_M <- AnnotationTrack(Gr_dmrt1_meth_M, name = "Male")

plotTracks(list(DataTrack(Gr_dmrt1_meth_F, name ="Female",ylim= c(0,100),col = "black"),
                DataTrack(Gr_dmrt1_meth_M, name ="Male",ylim = c(0,100),col = "black")),
           type = c("p","histogram","g"), background.title = "darkblue",
           from = 31110000, to = 31146500) # B, histogram


# Expression dmrt1
data_dmrt1_RNA <- read.csv("dmrt1_RNA.txt",sep= "\t")

df_dmrt1_RNA <- data.frame("chr"=rep(5,nrow(data_dmrt1_RNA)),"start"=data_dmrt1_RNA$Start,
                           "end"=data_dmrt1_RNA$End,"strand"=rep("*",nrow(data_dmrt1_RNA)))

Gr_dmrt1_RNA_F <- makeGRangesFromDataFrame(df_dmrt1_RNA)  
values(Gr_dmrt1_RNA_F) <- DataFrame(female = data_dmrt1_RNA$CF_RNA)

Gr_dmrt1_RNA_M <- makeGRangesFromDataFrame(df_dmrt1_RNA)  
values(Gr_dmrt1_RNA_M) <- DataFrame(male = data_dmrt1_RNA$TP_RNA)

dmrt1_RNA_track_F <- AnnotationTrack(Gr_dmrt1_RNA_F, name = "Female")
dmrt1_RNA_track_M <- AnnotationTrack(Gr_dmrt1_RNA_M, name = "Male")

plotTracks(list(DataTrack(Gr_dmrt1_RNA_F, name ="Female",ylim= c(0,3)),
                DataTrack(Gr_dmrt1_RNA_M, name ="Male",ylim = c(0,3))),
           type = c("g","histogram"), background.title = "darkblue",
           from = 31110000, to = 31146500)

# CpG islands and gene model dmrt1
geneModel <- read.csv("dmrt1_gene.csv")
grtrack <- GeneRegionTrack(geneModel,name = "dmrt1")

dmrt1_scale <- data.frame("chr"=5,"start"=31120000,"end"=31122000,"strand"="*")
Gr_dmrt1_scale <- makeGRangesFromDataFrame(dmrt1_scale)

dmrt1_CpG <- data.frame("chr"=c(5,5,5,5,5,5,5,5),
                        "start"=c(31110863,31118894,31128992,31129847,31131311,31135982,31137743,31143718),
                        "end"=c(31111062,31119258,31129191,31130089,31131824,31136331,31137993,31144087),"strand"="*")
Gr_dmrt1_CpG <- makeGRangesFromDataFrame(dmrt1_CpG)

plotTracks(list(grtrack,AnnotationTrack(Gr_dmrt1_CpG, name ="CpG")),
           type = c("g","histogram"), background.title = "darkblue",
           from = 31110000, to = 31146500) 
