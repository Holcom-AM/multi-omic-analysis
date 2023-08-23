library(tidyverse)
library(readxl)
library(ggrepel)
library(paletteer)
library(ggthemes)
library(ggfortify)


#RNA volcano plots
rna_degs1 <- read.csv("./input/angelina_code/corrected DGE day 1 CK10_vs_WT.csv")
rna_degs2 <- read.csv("./input/angelina_code/corrected DGE day 1 CL2355_vs_WT.csv")
rna_degs3 <- read.csv("./input/angelina_code/corrected DGE day 1 GL405_vs_WT.csv")


ggplot(rna_degs2, aes(log2FoldChange, -log(padj) ) ) + 
  geom_point( data=subset(rna_degs2, padj > 0.05  | abs(log2FoldChange) < 0.85), color="#E2DEDD", size=0.7) +   
  geom_point( data=subset(rna_degs2, padj < 0.05 & log2FoldChange > 0.85), color="#cf2d21", size=0.7) + 
  geom_point( data=subset(rna_degs2, padj < 0.05 & log2FoldChange < -0.85), color="#177fb5", size=0.7) + 
  geom_text_repel(data=subset(rna_degs2, padj < 0.05 & log2FoldChange > 4 ), aes(log2FoldChange,-log(padj),label=X), size=3, xlim = c(2.46,10) ) + 
  geom_text_repel(data=subset(rna_degs2, padj < 0.001 & log2FoldChange < -4.7 ), aes(log2FoldChange,-log(padj),label=X), size=3, xlim = c(-10,-1.5) ) +
  theme_linedraw() +
  theme(panel.grid = element_line("gray", 0.5), text = element_text(size = 10)) +
  xlab("log2(Aβ/WT)") +
  xlim(-9,9) +
  ylim(0,66)


ggplot(rna_degs1, aes(log2FoldChange, -log(padj) ) ) + 
  geom_point( data=subset(rna_degs1, padj > 0.05  | abs(log2FoldChange) < 0.85), color="#E2DEDD", size=0.7) +   
  geom_point( data=subset(rna_degs1, padj < 0.05 & log2FoldChange > 0.85), color="#cf2d21", size=0.7) + 
  geom_point( data=subset(rna_degs1, padj < 0.05 & log2FoldChange < -0.85), color="#177fb5", size=0.7) + 
  geom_text_repel(data=subset(rna_degs1, padj < 0.05 & log2FoldChange > 4 ), aes(log2FoldChange,-log(padj),label=X), size=3, xlim = c(2.46,10) ) + 
  geom_text_repel(data=subset(rna_degs1, padj < 0.001 & log2FoldChange < -4.7 ), aes(log2FoldChange,-log(padj),label=X), size=3, xlim = c(-10,-1.5) ) +
  theme_linedraw() +
  theme(panel.grid = element_line("gray", 0.5), text = element_text(size = 10)) +
  xlab("log2(Tau/WT)") +
  xlim(-9,9) +
  ylim(0,66)

ggplot(rna_degs3, aes(log2FoldChange, -log(padj) ) ) + 
  geom_point( data=subset(rna_degs3, padj > 0.05  | abs(log2FoldChange) < 0.85), color="#E2DEDD", size=0.7) +   
  geom_point( data=subset(rna_degs3, padj < 0.05 & log2FoldChange > 0.85), color="#cf2d21", size=0.7) + 
  geom_point( data=subset(rna_degs3, padj < 0.05 & log2FoldChange < -0.85), color="#177fb5", size=0.7) + 
  geom_text_repel(data=subset(rna_degs3, padj < 0.05 & log2FoldChange > 4 ), aes(log2FoldChange,-log(padj),label=X), size=3, xlim = c(2.46,10) ) + 
  geom_text_repel(data=subset(rna_degs3, padj < 0.001 & log2FoldChange < -4.7 ), aes(log2FoldChange,-log(padj),label=X), size=3, xlim = c(-10,-1.5) ) +
  theme_linedraw() +
  theme(panel.grid = element_line("gray", 0.5), text = element_text(size = 10)) +
  xlab("log2(Tau;Aβ/WT)") +
  xlim(-9,9) +
  ylim(0,66)