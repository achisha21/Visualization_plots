#############DESEQ2 for info.csv control vs tested##############################

setwd("/Users/achisha_saikia/Documents/Salivary glands/Mouse salivary glands/Par/CD1_vs_C57/All_Females")
library("DESeq2")
library(ggplot2)
library(ggrepel)
library("RColorBrewer")
library("pheatmap")
library(cowplot)
library('genefilter')
library('gplots')
library("limma")
library(apeglm)
library(EnhancedVolcano)
library(dplyr)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(genefilter)
library(PoiClaClu)

#Create the dds object with the original featurecounts and include batch in your design formula. Your design formula may look like ~ batch + condition 
#Conduct a differential expression analysis in the normal way for your factor of interest (e.g., condition) - the effect of batch will be 'adjusted'
#when the statistical inferences are made due to the fact that batch is included in the design formula
#Run the commands with limma::removeBatchEffect() if you want to use your data matrix downstream for clustering, PCA, machine learning
#########################################################################################################################################################
# Read in counts table 'cts' from featureCount output
getwd()
#setwd("/Users/achisha_saikia/Documents/Salivary glands/Mouse salivary glands/Par/All_Females")
dat <- read.csv("readCounts.tsv",
                header = TRUE, 
                sep = '\t',
                row.names = 'Geneid',
                skip=1)

dat.subset<- dat[,c(12,14,16,13,15,17)]
dat<-dat.subset
# Drop extra columns from featureCount output
# Drop extra columns from featureCount output
#drops <- c("Chr","Start","End","Strand","Length","X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.Par_1A_Mal.bam", "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.Par_1B_Mal.bam", "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.Par_2A_Mal.bam", "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.Par_2B_Mal.bam", "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.Par_3A_Mal.bam", "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.Par_3B_Mal.bam",
#           "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SL_1A_Mal.bam","X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SL_1B_Mal.bam", "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SL_2A_Mal.bam","X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SL_2B_Mal.bam", "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SL_3A_Mal.bam","X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SL_3B_Mal.bam", "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SL_4A_Fem.bam","X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SL_4B_Fem.bam", "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SL_5A_Fem.bam", "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SL_5B_Fem.bam", "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SL_6A_Fem.bam","X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SL_6B_Fem.bam",
#           "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SM_1A_Mal.bam", "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SM_1B_Mal.bam","X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SM_2A_Mal.bam", "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SM_2B_Mal.bam", "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SM_3A_Mal.bam", "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SM_3B_Mal.bam","X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SM_4A_Fem.bam","X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SM_4B_Fem.bam", "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SM_5A_Fem.bam", "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SM_5B_Fem.bam","X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SM_6A_Fem.bam", "X.projects.academic.omergokc.Achisha.Workflow.outputs.mapping.SM_6B_Fem.bam" )
#head(drops)
#dat <- dat[ , !(names(dat) %in% drops)]


colnames(dat)
colnames(dat) <- c("Par_4A_CD1","Par_5A_CD1","Par_6A_CD1","Par_4B_C57","Par_5B_C57","Par_6B_C57")

# Read in sample information table 'coldata'
info <- read.csv("info.csv", header = TRUE, 
                 sep = ',',
                 stringsAsFactors = TRUE)

dds <- DESeqDataSetFromMatrix(countData = dat, colData = info, design = ~ Condition)

#remove lowly expressed genes
keep<-rowSums(counts(dds))>=10
dds<-dds[keep,]

#main DESeq
dds<-DESeq(dds)

#Compare and contrast
result<- results(dds, contrast = c("Condition", "C57", "CD1"))
resultsNames(dds)

#res1_B6_vs_CD1 <- results(dds, contrast = list("Condition_B6_vs_CD1"))
#res1_B6_vs_C57 <- results(dds, contrast = list("Condition_B6_vs_C57"))
#res1_CD1_vs_C57 <- results(dds, contrast = list("Condition_B6_vs_CD1", "Condition_B6_vs_C57"))
#resultsNames(ddsDE)

#export normalized read counts
normCounts <- counts(dds, normalized=T)
write.csv(normCounts, "normalized_All_Par_Fem_counts.csv")

#DESeq results
res<-results(dds, alpha = 0.05)

resOrdered <- res[order(res$padj),]
write.csv(resOrdered, "deseq_All_Par_Fem_padj_counts.csv")

summary(res)
#resultsNames(dds)


####################################plotting#####################################

normCount <- read.csv("normalized_All_Par_Fem_counts.csv", row.names = 1)
deSeqRes <- read.csv("deseq_All_Par_Fem_padj_counts.csv", row.names = 1)
View(deSeqRes)

#deSeqRes$sig <- ifelse((deSeqRes$log2FoldChange > 1.5 & deSeqRes$padj <0.05), "yes", "no")
deSeqRes$sig <- ifelse(deSeqRes$padj <= 0.05, "yes", "no")
deSeqRes <- na.omit(deSeqRes)

write.csv(deSeqRes, "Deseq2_significant_res_All_Par_Fem.csv")

View(deSeqRes)
#plotMA
ggplot(deSeqRes, aes(x= log10(baseMean), y= log2FoldChange, color = sig)) + geom_point()  + theme_minimal()

#volcana plot
ggplot(deSeqRes, aes(x= log2FoldChange, y = -log10(padj), color= sig)) + geom_point() + scale_color_manual(values=c("orange","blue")) + theme_minimal() 

#pheatmap
signi <- subset(deSeqRes,padj <= 0.05)
#View(signi)
allsig <- merge(normCount, signi, by = 0)
#View(allsig)

sigCounts <- allsig[,2:7]
row.names(sigCounts) <- allsig$Row.names
#View(sigCounts)
pheatmap(log2(sigCounts +1), scale = 'row', show_rownames = F, treeheight_row = 0, treeheight_col = 0, cluster_cols = FALSE )
heatmap(sigCounts, cluster_cols = T)
#PCA
vsdata<-vst(dds, blind = FALSE)
plotPCA(vsdata, intgroup = "Condition") 
plotDispEsts(dds)


######counts data

d <- plotCounts(dds, gene="Pla2g2a", intgroup=c("Sample","Condition", "Design"), 
                returnData=TRUE)
ggplot(d, aes(x=Sample, y=count,color=Condition,shape=Condition)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + theme_minimal() + ggtitle("Pla2g2a")



#####################################ENHANCED VOLCANO PLOT#########################################################################

library(EnhancedVolcano)
View(df)

EnhancedVolcano(df, x="log2FoldChange", y="padj",, lab=df$X, pCutoff = 0.05, FCcutoff = 1)

#################################################################################################################################

deSeqRes <- read.csv("deseq_All_Par_Fem_padj_counts.csv", row.names = 1)

deSeqRes$sig <- ifelse(deSeqRes$padj <= 0.05, "yes", "no")
deSeqRes <- na.omit(deSeqRes)
#plotMA
ggplot(deSeqRes, aes(x= log10(baseMean), y= log2FoldChange, color = sig)) + geom_point() + ggtitle("MA Plot") + theme_bw() 

#volcana plot
ggplot(deSeqRes, aes(x= log2FoldChange, y = -log10(padj), color= sig)) + geom_point() + scale_color_manual(values=c("pink","blue")) + theme_minimal()  + ggtitle("Volcano Plot")

#pheatmap
signi <- subset(deSeqRes,padj <= 0.05)
#View(signi)
allsig <- merge(normCount, signi, by = 0)
#View(allsig)

sigCounts <- allsig[,2:7]
row.names(sigCounts) <- allsig$Row.names
#View(sigCounts)
pheatmap(log2(sigCounts +1), scale = 'row', show_rownames = F, treeheight_row = 0, treeheight_col = 0,cluster_cols = FALSE) 


#distance-to-distance
sampleDists <- dist(t(assay(sigCounts)))
head(assay(deSeqRes), 3)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
ntd <- normTransform(dds)



library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Condition","Design")])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Condition, vsd$Design, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#PCAplot
plotPCA(vsd, intgroup=c("Condition", "Design"))

#####################Upreg and downreg volcano###############################

df <- read.csv("deseq_All_Par_Fem_padj_counts.csv")

ggplot(data=df, aes(x=log2FoldChange, y=-log(padj))) + geom_point()

# Convert directly in the aes()
p <- ggplot(data=df, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point()
# Add more simple "theme"
p <- ggplot(data=df, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point() + theme_minimal()
# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

p2

# add a column of NAs
df$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df$diffexpressed[df$log2FoldChange > 1.0 & df$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"

df$diffexpressed[df$log2FoldChange < -1.0 & df$pvalue < 0.05] <- "DOWN"

View(df)
# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=df, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)


# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
df$delabel <- NA
df$delabel[df$diffexpressed != "NO"] <- df$gene_symbol[df$diffexpressed != "NO"]

View(df)
ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()


ggplot(data = df, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "#FFDB6D"), # to set the colours of our variable<br /><br /><br />
                     labels = c("Downregulated", "Not significant", "Upregulated")) # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)</p><br /><br />


# Edit axis labels and limits
ggplot(data = df, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 50), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Severe', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) # to customise the breaks in the x axis

# Note. with coord_cartesian() even if we have genes with p-values or log2FC ourside our limits, they will still be plotted.

# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
df$delabel <- ifelse(df$X %in% head(df[order(df$padj), "gene_symbol"], 30), df$gene_symbol, NA)


ggplot(data = df, aes(x = log2FoldChange , y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 250), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Severe', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  ggtitle('Thf-like cells in severe COVID vs healthy patients') + # Plot title 
  geom_text_repel(max.overlaps = Inf) # To show all labels 


