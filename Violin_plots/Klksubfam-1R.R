setwd("~/Documents/Salivary glands/Mouse salivary glands/violin_plots")

#Rcode for violin plots

library(RColorBrewer)
library(ggplot2)
genes <- read.csv("klksubfam-1.csv")
head(genes)

mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))

display.brewer.pal(11, "Paired")
brewer.pal(11, "Paired")

#tiff('test.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')
p=ggplot(genes, aes(x=Sex, y=log(Expression+1,10)))
p+geom_violin(aes(col=Gene), alpha=0.3)+
  facet_grid(Gene~Strain+Glands, scales="free") +
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired")+
  geom_jitter(width=0.1)+
  theme_light()
#dev.off()