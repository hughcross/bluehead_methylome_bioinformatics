############################################################
# RNA-seq data analysis                                    #
############################################################

# Load libraries
library(DESeq2)
library(ggpubr)
library(ggplot2)
library(eulerr)
library(cowplot)
library(org.Dr.eg.db)
library(dplyr)
library(topGO)

# Change folder
setwd("~/Syncplicity Folders/BHW_ToddOrtegaLiu2018/Figures/Figure_2_PCA_analysis")

############################################################
# Analysis of gonad and forebrain samples using DESeq2     #
############################################################

# Load data
Gonad_expression <- read.csv('~/Syncplicity Folders/BHW_ToddOrtegaLiu2018/Figures/Data/gonad_cmd_matrix_2018-06-05.isoform.counts.matrix', header=TRUE, sep="\t", row.names = "X")
Forebrain_expression <- read.csv('~/Syncplicity Folders/BHW_ToddOrtegaLiu2018/Figures/Data/forebrain_cmd_matrix_2018-06-05.isoform.counts.matrix', header=TRUE, sep="\t", row.names = "X")

# Bind tables
all_expression <- cbind(Gonad_expression,Forebrain_expression)

# Transform to integer
all_expression <- as.data.frame(sapply(all_expression, as.integer),row.names = row.names(all_expression))

# Remove all rows with less than 10 counts across all samples
low_count_mask <- rowSums(all_expression) < 10
sprintf("Removing %d low-count genes (%d remaining).", sum(low_count_mask), 
        sum(!low_count_mask))

all_expression_clean <- all_expression[!low_count_mask,]

order <- c(sort(colnames(G_expression)),sort(colnames(F_expression)))

data_sc <- all_expression_clean[c(sort(colnames(all_expression_clean)))]

samples_sc <- data.frame(row.names=c(colnames(data_sc)),
                         stage=as.factor(c(rep("G_Ctrl_female",6), rep("G_Stage_1",3),rep("G_Stage_2",7),
                                           rep("G_Stage_3",3), rep("G_Stage_4",3), rep("G_Stage_5a",3),
                                           rep("G_Stage_5b",5), rep("G_Stage_6",3), rep("G_TP_male",8),
                                           rep("F_Ctrl_female",6), rep("F_Stage_1",3),rep("G_Stage_2",7),
                                           rep("G_Stage_3",3), rep("G_Stage_4",3), rep("G_Stage_5a",3),
                                           rep("G_Stage_5b",5), rep("G_Stage_6",3), rep("G_TP_male",8))))

# Process data DESeq2
dds_sc <- DESeqDataSetFromMatrix(countData = data_sc, colData=samples_sc, design= ~stage)
dds_sc <- estimateSizeFactors(dds_sc)
vsd_sc <- varianceStabilizingTransformation(dds_sc)

# PCA forebrain and gonads analysis and plot
rv <- rowVars(assay(vsd_sc))
select <- order(rv, decreasing=TRUE)[seq_len(min(10000, length(rv)))]
pca <- prcomp(t(assay(vsd_sc)[select,]))
summary(pca)

PCA_values <- as.data.frame(pca$x)
PCA_values$shape <- c(rep("17",41),rep("16",41))
PCA_values$group <- c(rep("CF",6), rep("S1",3),rep("S2",7),
                      rep("S3",3), rep("S4",3), rep("S5a",3),
                      rep("S5b",5), rep("S6",3), rep("TP",8),
                      rep("CF",6), rep("S1",3),rep("S2",7),
                      rep("S3",3), rep("S4",3), rep("S5a",3),
                      rep("S5b",5), rep("S6",3), rep("TP",8))


cbPalette <- c("#8f0000", "#ff0000", "#ff8000", "#ffd940",
               "#b3ebff", "#4dd9ff", "#1c79b2", "#2463ff", "#00048b")

plot_figure_2_a <- ggplot(PCA_values, aes(x = PC1, y = PC2))+
  geom_point(aes(shape = PCA_values$shape,color = PCA_values$group),size=5,alpha=0.8) +
  scale_color_manual(values=cbPalette) + 
  labs(x = "PC1: 98% variance", y = "PC2: 1% variance") +
  theme(axis.text = element_text(size=14),axis.title = element_text(size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

plot_figure_2_a

# Clean environment
rm(list=ls())

############################################################
# Analysis gonad samples using DESeq2                      #
############################################################

# Load data
Gonad_expression <- read.csv('~/Syncplicity Folders/BHW_ToddOrtegaLiu2018/Figures/Data/gonad_cmd_matrix_2018-06-05.isoform.counts.matrix', header=TRUE, sep="\t", row.names = "X")

# Transform to integer
Gonad_expression_data <- as.data.frame(sapply(Gonad_expression, as.integer),row.names = row.names(Gonad_expression))

# Remove all rows with less than 10 counts across all samples.
low_count_mask <- rowSums(Gonad_expression_data) < 10
sprintf("Removing %d low-count genes (%d remaining).", sum(low_count_mask), 
        sum(!low_count_mask))

Gonad_expression_clean <- Gonad_expression_data[!low_count_mask,]

data_sc <- Gonad_expression_clean[c(sort(colnames(Gonad_expression_clean)))]

samples_sc <- data.frame(row.names=c(colnames(data_sc)),
                         stage=as.factor(c(rep("G_Ctrl_female",6), rep("G_Stage_1",3),rep("G_Stage_2",7),
                                           rep("G_Stage_3",3), rep("G_Stage_4",3), rep("G_Stage_5a",3),
                                           rep("G_Stage_5b",5), rep("G_Stage_6",3), rep("G_TP_male",8))))

# Process data DESeq2
dds_sc <- DESeqDataSetFromMatrix(countData = data_sc, colData=samples_sc, design= ~stage)
dds_sc <- estimateSizeFactors(dds_sc)
vsd_sc <- varianceStabilizingTransformation(dds_sc)

# PCA forebrain and gonads analysis and plot
rv <- rowVars(assay(vsd_sc))
select <- order(rv, decreasing=TRUE)[seq_len(min(10000, length(rv)))]
pca <- prcomp(t(assay(vsd_sc)[select,]))
summary(pca)

PCA_values <- as.data.frame(pca$x)
PCA_values$group <- c(rep("CF",6), rep("S1",3),rep("S2",7),
                      rep("S3",3), rep("S4",3), rep("S5a",3),
                      rep("S5b",5), rep("S6",3), rep("TP",8))

cbPalette <- c("#8f0000", "#ff0000", "#ff8000", "#ffd940",
               "#b3ebff", "#4dd9ff", "#1c79b2", "#2463ff", "#00048b")

plot_figure_2_b <- ggplot(PCA_values, aes(x = PC1, y = PC2))+
  geom_point(aes(color = PCA_values$group),size=5,alpha=0.9) +
  scale_color_manual(values=cbPalette) +
  stat_conf_ellipse(aes(color = PCA_values$group,fill=PCA_values$group),alpha=0.2,geom = "polygon",bary = TRUE) +
  labs(x = "PC1: 58% variance", y = "PC2: 17% variance") +
  theme(axis.text = element_text(size=14),axis.title = element_text(size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

plot_figure_2_b

# Clean environment
rm(list=ls())

############################################################
# Analysis of PCA data for gonad samples                   #
############################################################

# Load data
Gonad_expression <- read.csv('~/Syncplicity Folders/BHW_ToddOrtegaLiu2018/Figures/Data/gonad_cmd_matrix_2018-06-05.isoform.counts.matrix', header=TRUE, sep="\t", row.names = "X")

# Transform to integer
Gonad_expression_data <- as.data.frame(sapply(Gonad_expression, as.integer),row.names = row.names(Gonad_expression))

# Remove all rows with less than 10 counts across all samples.
low_count_mask <- rowSums(Gonad_expression_data) < 10
sprintf("Removing %d low-count genes (%d remaining).", sum(low_count_mask), 
        sum(!low_count_mask))

Gonad_expression_clean <- Gonad_expression_data[!low_count_mask,]

data_sc <- Gonad_expression_clean[c(sort(colnames(Gonad_expression_clean)))]

samples_sc <- data.frame(row.names=c(colnames(data_sc)),
                         stage=as.factor(c(rep("G_Ctrl_female",6), rep("G_Stage_1",3),rep("G_Stage_2",7),
                                           rep("G_Stage_3",3), rep("G_Stage_4",3), rep("G_Stage_5a",3),
                                           rep("G_Stage_5b",5), rep("G_Stage_6",3), rep("G_TP_male",8))))

# Process data DESeq2
dds_sc <- DESeqDataSetFromMatrix(countData = data_sc, colData=samples_sc, design= ~stage)
dds_sc <- estimateSizeFactors(dds_sc)
vsd_sc <- varianceStabilizingTransformation(dds_sc)

# PCA gonad analysis
rv <- rowVars(assay(vsd_sc))
select <- order(rv, decreasing=TRUE)[seq_len(min(10000, length(rv)))]
pca <- prcomp(t(assay(vsd_sc)[select,]))

pca_loadings <- as.data.frame(pca$rotation,row.names = row.names(pca$rotation))

# Plot loadings
plot(pca$rotation[,1], pca$rotation[,2],
     xlim=c(-0.05,0.05), ylim=c(-0.05,0.05),
     main='Loadings for PC1 vs. PC2', type = "p")

# Thresholds
ranges_perc = c(0.05,0.95) # 500 values

PC1_quantiles <- quantile(pca_loadings$PC1,ranges_perc)
PC1_min = as.numeric(PC1_quantiles[1])
PC1_max = as.numeric(PC1_quantiles[2])

PC2_quantiles <- quantile(pca_loadings$PC2,ranges_perc)
PC2_min = as.numeric(PC2_quantiles[1])
PC2_max = as.numeric(PC2_quantiles[2])

# Components
comp_1_left = pca_loadings[pca_loadings$PC1<=PC1_min,]
comp_1_right = pca_loadings[pca_loadings$PC1>=PC1_max,]
comp_2_up = pca_loadings[pca_loadings$PC2>=PC2_max,]
comp_2_down = pca_loadings[pca_loadings$PC2<=PC2_min,]

# Euler diagram
left_names = rownames(comp_1_left)
right_names = rownames(comp_1_right)
up_names = rownames(comp_2_up)
down_names = rownames(comp_2_down)

# Intersect
intersect_left_up = intersect(left_names,up_names)
intersect_left_down = intersect(left_names,down_names)
intersect_right_up = intersect(right_names,up_names)
intersect_right_down = intersect(right_names,down_names)

LU = length(intersect_left_up)
LD = length(intersect_left_down)
RU = length(intersect_right_up)
RD = length(intersect_right_down)

strict_left = setdiff(left_names,union(intersect_left_up,intersect_left_down))
strict_right = setdiff(right_names,union(intersect_right_up,intersect_right_down))
strict_up = setdiff(up_names,union(intersect_left_up,intersect_right_up))
strict_down = setdiff(down_names,union(intersect_left_down,intersect_right_down))

L = length(strict_left)
R = length(strict_right)
U = length(strict_up)
D = length(strict_down)

length(left_names)==L+LU+LD
length(right_names)==R+RU+RD
length(up_names)==U+LU+RU
length(down_names)==D+LD+RD

euler_data <- c(Not_committed = U, Committed = D,Male = R, Female = L,"Male&Committed" = RD,
                "Male&Not_committed" = RU,"Female&Committed" = LD,"Female&Not_committed" = LU)

euler_fit <- euler(euler_data)

plot(euler_fit,quantities = TRUE,labels = list(cex = 2))

############################################################
# Component trends                                         #
############################################################

## Extract data DESeq2 analysis
normalizedCounts = as.data.frame(counts(dds_sc, normalized = TRUE))
colnames(normalizedCounts) = c(rep("CF",6), rep("S1",3),rep("S2",7), rep("S3",3), rep("S4",3),
                               rep("S5a",3), rep("S5b",5), rep("S6",3), rep("TP",8))

mean_CF = rowMeans(normalizedCounts[,1:6])
mean_S1 = rowMeans(normalizedCounts[,7:9])
mean_S2 = rowMeans(normalizedCounts[,10:16])
mean_S3 = rowMeans(normalizedCounts[,17:19])
mean_S4 = rowMeans(normalizedCounts[,20:22])
mean_S5a = rowMeans(normalizedCounts[,23:25])
mean_S5b = rowMeans(normalizedCounts[,26:30])
mean_S6 = rowMeans(normalizedCounts[,31:33])
mean_TP = rowMeans(normalizedCounts[,34:41])

mean_df <- data.frame(cbind(mean_CF,mean_S1,mean_S2,mean_S3,mean_S4,mean_S5a,mean_S5b,mean_S6,mean_TP))
mean_df_log <- log2(mean_df)

female <- mean_df_log[strict_left,]
male <- mean_df_log[strict_right,]
transitionary <- mean_df_log[strict_up,]
differentiated <- mean_df_log[strict_down,]
male_differentiated <- mean_df_log[intersect_right_down,]

# Female component
subset_female <- data.frame(t(female),stringsAsFactors = FALSE)
subset_female$stage <- as.numeric(seq(1,9))

subset_female_long <- melt(subset_female, id="stage")
female_trend_plot <- ggplot(data=subset_female_long,aes(x=stage, y=value, colour=variable)) +  geom_line(size=1) + 
  theme(legend.position="none",plot.title = element_text(size=32)) + ggtitle("PCA Female component") +
  labs(x = "Stage", y = "log2 normalized expression") +
  theme(axis.text.x = element_blank(),axis.title = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

female_trend_plot

# Male component
subset_male <- data.frame(t(male),stringsAsFactors = FALSE)
subset_male$stage <- as.numeric(seq(1,9))

subset_male_long <- melt(subset_male, id="stage")
male_trend_plot <- ggplot(data=subset_male_long,aes(x=stage, y=value, colour=variable)) +  geom_line(size=1) + 
  theme(legend.position="none",plot.title = element_text(size=32)) + ggtitle("PCA Male component") +
  labs(x = "Stage", y = "log2 normalized expression") +
  theme(axis.text.x = element_blank(),axis.title = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

male_trend_plot

# Differentiated component
subset_differentiated <- data.frame(t(differentiated),stringsAsFactors = FALSE)
subset_differentiated$stage <- as.numeric(seq(1,9))

subset_differentiated_long <- melt(subset_differentiated, id="stage")
differentiated_trend_plot <- ggplot(data=subset_differentiated_long,aes(x=stage, y=value, colour=variable)) +  geom_line(size=1) + 
  theme(legend.position="none",plot.title = element_text(size=32)) + ggtitle("PCA differentiated component") +
  labs(x = "Stage", y = "log2 normalized expression") +
  theme(axis.text.x = element_blank(),axis.title = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

differentiated_trend_plot

# Transitionary component
subset_transitionary <- data.frame(t(transitionary),stringsAsFactors = FALSE)
subset_transitionary$stage <- as.numeric(seq(1,9))

subset_transitionary_long <- melt(subset_transitionary, id="stage")
transitionary_trend_plot <- ggplot(data=subset_transitionary_long,aes(x=stage, y=value, colour=variable)) +  geom_line(size=1) + 
  theme(legend.position="none",plot.title = element_text(size=32)) + ggtitle("PCA Transitionary component") +
  labs(x = "Stage", y = "log2 normalized expression") +
  theme(axis.text.x = element_blank(),axis.title = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

transitionary_trend_plot

# Male differentiated component
subset_male_differentiated <- data.frame(t(male_differentiated),stringsAsFactors = FALSE)
subset_male_differentiated$stage <- as.numeric(seq(1,9))

subset_male_differentiated_long <- melt(subset_male_differentiated, id="stage")
male_differentiated_trend_plot <- ggplot(data=subset_male_differentiated_long,aes(x=stage, y=value, colour=variable)) +  geom_line(size=1) + 
  theme(legend.position="none",plot.title = element_text(size=32)) + ggtitle("PCA Male Committed component") +
  labs(x = "Stage", y = "log2 normalized expression") +
  theme(axis.text.x = element_blank(),axis.title = element_text(size=20))

male_differentiated_trend_plot

plot_trends <- plot_grid(female_trend_plot,male_trend_plot,
                         transitionary_trend_plot,differentiated_trend_plot,ncol= 1)
plot_trends

############################################################
# Gene ontology analysis                                   #
############################################################

## Gene universe
annotations_raw <- read.csv('~/Syncplicity Folders/BHW_ToddOrtegaLiu2018/Figures/Data/Gonad_isoform_lognorm_counts_annotated.csv', header=TRUE,stringsAsFactors = FALSE)
annotations_all_genes <- annotations_raw[,c(1,43:47)]
annotations_all_genes <- annotations_all_genes[!is.na(annotations_all_genes$Zfish_name),]

gene_ontology_analysis <- function(gene_list,ontology){
  # Create GO gene lists
  gene_list = rownames(female)
  gene_annotations <- data.frame("contigs" = gene_list,stringsAsFactors = FALSE)
  gene_annotations <- left_join(gene_annotations,annotations_all_genes, by = "contigs" )
  
  inSelection = gene_annotations$Zfish_name
  inUniverse = annotations_all_genes$Zfish_name
  
  geneList <- factor(as.integer(inUniverse %in% inSelection))
  names(geneList) <- inUniverse
  
  tgd <- new("topGOdata", ontology=ontology, allGenes = geneList, nodeSize=5,
             annot=annFUN.org, mapping="org.Dr.eg.db", ID = "ensembl" )
  
  ## Run tests
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
  resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
  resultTopGO.weight01 <- runTest(tgd, algorithm = "weight01", statistic = "Fisher" )

  ## Create table and save results
  GO_results <- GenTable(tgd, Fisher.elim = resultTopGO.elim, 
                               Fisher.classic = resultTopGO.classic,
                               weight01 = resultTopGO.weight01,
                               orderBy = "weight01" , topNodes = 200)
  
  return(GO_results)
  }

male_GO_BP <- gene_ontology_analysis(rownames(male),"BP")
male_GO_MF <- gene_ontology_analysis(rownames(male),"MF")
female_GO_BP <- gene_ontology_analysis(rownames(female),"BP")
female_GO_MF <- gene_ontology_analysis(rownames(female),"MF")
transitionary_GO_BP <- gene_ontology_analysis(rownames(transitionary),"BP")
transitionary_GO_MF <- gene_ontology_analysis(rownames(transitionary),"MF")
differentiated_GO_BP <- gene_ontology_analysis(rownames(differentiated),"BP")
differentiated_GO_MF <- gene_ontology_analysis(rownames(differentiated),"MF")
