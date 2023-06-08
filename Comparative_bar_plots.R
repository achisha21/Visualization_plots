
setwd("~/Documents/Salivary glands/Mouse salivary glands/BarPlots/")

SalExp <- read.csv("Human_SG.csv")

library(ggplot2)
library(cowplot)

#For Sjorgren's syndrome

Genes = c("STAT4", "IL12A", "IRF5", "TNPO3", "FAM167A", "BLK", "DDX6","CXCR5", "TNIP1","TNFAIP3","PTTG1", "PRDM1","DGKQ","IRAK1BP1","GTF2I")

c = SalExp[SalExp$hgnc_symbol %in% Genes,] #Substets rows by column value, 
                                          #in this case, it substets the rows with the names in Genes
#-------------------------Human Salivary Glands -----------------------

# Retrieving the Expression values for each of the Glands for Adults
SM <- data.frame(c$hgnc_symbol, log(c[,27],10), "SM") 
colnames(SM)=c("Genes", "log10_Expression", "Gland")
SM[SM == "-Inf"] <- c(0)

PAR <- data.frame(c$hgnc_symbol, log(c[,28],10), "PAR") 
colnames(PAR)=c("Genes", "log10_Expression", "Gland")
PAR[PAR == "-Inf"] <- c(0)

SL <- data.frame(c$hgnc_symbol, log(c[,29],10), "SL") 
colnames(SL)=c("Genes", "log10_Expression", "Gland")
SL[SL == "-Inf"] <- c(0)

Human = rbind(SM, PAR, SL) #paste one below the other
colnames(Human)=c("Genes", "log10_Expression", "Gland")
#Plots
#Submandibular
SM_gene <- ggplot(Human, aes(fill=Gland, x=factor(Genes), y=log10_Expression)) + 
  geom_bar(position="dodge",stat="identity", width=0.5, color = "black") + 
  theme_minimal() + 
  geom_text(position = position_dodge2(width = .9, preserve = "single"), #positioning the number labels above the bars
            aes(label=round(Human[,2], 2), #I want only two algorithms after the comma
                hjust=+0.5, vjust=-0.3,), size = 2.5) + #To add the values to each column 
  labs(y="Log10(Expression)", x="Genes", title="Gene Expression in the Adult Salivary Glands") +
  theme(plot.title = element_text(hjust = 0.5)) + #to centralize the title of the graph
  scale_x_discrete(limit = Genes) + # Ordering the genes in the x axis as the same as the variable
  theme(axis.text.x = element_text(face = "bold",colour = c("blue", "blue", "green", "green", "green", "blue", "blue","red", "red", "red", "blue","red","green","green","green","green","purple", "purple", "purple")))

#-------------------------Mouse Salivary Glands -----------------------

Genes_m = c("STAT4", "IL12A", "IRF5", "TNPO3", "FAM167A", "BLK", "DDX6","CXCR5", "TNIP1","TNFAIP3","PTTG1", "PRDM1","DGKQ","IRAK1BP1","GTF2I")

Mice <- read.csv("Mice_major_SG_Gene_Exp.csv")
Mice$Geneid = toupper(Mice$Geneid) # genes in upper case
d = Mice[Mice$Geneid %in% Genes_m,] #Substets rows by column value, 

#Subset each gland
Mice_Par = d[,2:13]
Par_means = rowMeans(Mice_Par)
Par_m <- data.frame(d[,1], log(Par_means, 10), "PAR")
Par_m[Par_m == "-Inf"] <- c(0)
colnames(Par_m)=c("Genes", "log10_Expression", "Gland")

Mice_SL = d[,14:25]
SL_means = rowMeans(Mice_SL)
SL_m <- data.frame(d[,1], log(SL_means, 10), "SL") 
SL_m[SL_m == "-Inf"] <- c(0)
colnames(SL_m)=c("Genes", "log10_Expression", "Gland")

Mice_SM = d[,26:37]
SM_means = rowMeans(Mice_SM)
SM_m <- data.frame(d[,1], log(SM_means, 10), "SM") 
SM_m[SM_m == "-Inf"] <- c(0)
colnames(SM_m)=c("Genes", "log10_Expression", "Gland")

Mouse = rbind(Par_m, SM_m, SL_m)

#Plots
#Submandibular
Gland_mouse <- ggplot(Mouse, aes(fill=Gland, x=factor(Genes), y=log10_Expression)) +
  geom_bar(position="dodge",stat="identity", width=0.5, color = "black") + 
  theme_minimal() + 
  geom_text(position = position_dodge2(width = .9, preserve = "single"), #positioning the number labels above the bars
            aes(label=round(Mouse[,2], 2), #I want only two algorithms after the comma
                hjust=+0.5, vjust=-0.3,), size = 2.5) + #To add the values to each column 
  labs(y="Log10(Expression)", x="Genes", title="Gene Expression in Mouse Salivary Glands") +
  theme(plot.title = element_text(hjust = 0.5)) + #to centralize the title of the graph
  scale_x_discrete(limit = Genes_m) + # Ordering the genes in the x axis as the same as the variable
  theme(axis.text.x = element_text(face = "bold", colour = c("blue", "blue", "green", "green", "green", "blue", "blue","red", "red", "red", "blue", "red","green","green","green","green","purple", "purple", "purple")))
                                             

#-----------------------------------Mice liver and Pancreas --------------------------------------#

Mice_LP <- read.csv("Mice_Liv_Panc_Gene_Exp.csv")

Genes_m = c("STAT4", "IL12A", "IRF5", "TNPO3", "FAM167A", "BLK", "DDX6","CXCR5", "TNIP1","TNFAIP3","PTTG1", "PRDM1","DGKQ","IRAK1BP1","GTF2I")

Mice_LP$Geneid = toupper(Mice_LP$Geneid) # genes in upper case
e = Mice_LP[Mice_LP$Geneid %in% Genes_m,] #Substets rows by column value, 

#Subset each gland
Mice_Liver = e[,2:13]
Liver_means = rowMeans(Mice_Liver)
Liver_m <- data.frame(e[,1], log(Liver_means, 10), "Liver")
Liver_m[Liver_m == "-Inf"] <- c(0)
colnames(Liver_m)=c("Genes", "log10_Expression", "Gland")

Mice_Pan = e[,14:25]
Pan_means = rowMeans(Mice_Pan)
Pan_m <- data.frame(e[,1], log(Pan_means, 10), "Pancreas") 
Pan_m[Pan_m == "-Inf"] <- c(0)
colnames(Pan_m)=c("Genes", "log10_Expression", "Gland")

Mice_LP_2 = rbind(Liver_m, Pan_m)

#Plots
#Submandibular
LP_mouse <- ggplot(Mice_LP_2, aes(fill=Gland, x=factor(Genes), y=log10_Expression)) +
  geom_bar(position="dodge",stat="identity", width=0.5, color = "black") + 
  theme_minimal() + 
  geom_text(position = position_dodge2(width = .9, preserve = "single"), #positioning the number labels above the bars
            aes(label=round(Mice_LP_2[,2], 2), #I want only two algorithms after the dot
                hjust=+0.5, vjust=-0.3,), size = 2.5) + #To add the values to each column 
  labs(y="Log10(Expression)", x="Genes", title="Gene Expression in Mouse Liver and Pancreas") +
  theme(plot.title = element_text(hjust = 0.5)) + #to centralize the title of the graph
  scale_x_discrete(limit = Genes_m) + # Ordering the genes in the x axis as the same as the variable
  theme(axis.text.x = element_text(face = "bold", colour = c("blue", "blue", "green", "green", "green", "blue", "blue","red", "red", "red", "blue", "red","green","green","green","green","purple", "purple", "purple")))

pdf("Comparative_bar_plots.pdf", width=20, height=12)

plot_grid(SM_gene,Gland_mouse,LP_mouse, ncol = 1)

dev.off()

