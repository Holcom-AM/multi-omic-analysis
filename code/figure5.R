library(UpSetR)
library(ggplot2)
library(grid)
library(plyr)
library(CellPlot)
library(tidyverse)

# Read the CSV file into a data frame
day3vsday1_data <- read.csv("./input/angelina_code/day 3 vs day 1 upset plot data.csv")
colnames(day3vsday1_data) <- c("gene","WT","Tau","AB","Tau:AB")

# Create the upset plot
pdf(file = "./output/figure5_a1.pdf",width=5.38,height=5.13)
upset(day3vsday1_data,
      mainbar.y.label = "Intersections",
      sets.x.label = "Insoluble proteins",
      sets.bar.color = c("grey60","#6699CC","#EE99AA","#EECC66"),
      text.scale = 2,
      point.size = 3,
      order.by = "freq")
dev.off()

##############################################################################
#############day 5 vs day 1 #################################################
# Read the CSV file into a data frame
day5vsday1_data <- read.csv("./input/angelina_code/day 5 vs day 1 upset plot data.csv")
colnames(day5vsday1_data) <- c("gene","WT","Tau","AB","Tau:AB")

# Create the upset plot
pdf(file = "./output/figure5_a2.pdf",width=5.38,height=5.13)
upset(day5vsday1_data,
      mainbar.y.label = "Intersections",
      sets.x.label = "Insoluble proteins",
      sets.bar.color = c("grey60","#6699CC","#EE99AA","#EECC66"),
      text.scale = 2,
      point.size = 3,
      order.by = "freq")
dev.off()

##############################################################################
#############day 8 vs day 1 #################################################
# Read the CSV file into a data frame
day8vsday1_data <- read.csv("./input/angelina_code/day 8 vs day 1 upset plot data.csv")
colnames(day8vsday1_data) <- c("gene","WT","Tau","AB","Tau:AB")

pdf(file = "./output/figure5_a3.pdf",width=5.38,height=5.13)
upset(day8vsday1_data,
      mainbar.y.label = "Intersections",
      sets.x.label = "Insoluble proteins",
      sets.bar.color = c("grey60","#6699CC","#EE99AA","#EECC66"),
      text.scale = 2,
      point.size = 3,
      order.by = "freq")
dev.off()

#plot enrichment results
D8_vs_D1_WT_protein_GO_MF <- read.csv("./input/angelina_code/day 8 vs day 1 WT protein GO_MF.csv")
pdf(file = "./output/figure5_b.pdf",width=5.38,height=5.13)
cell.plot(x = setNames(D8_vs_D1_WT_protein_GO_MF$X.log10.padj., D8_vs_D1_WT_protein_GO_MF$Term),
          cells = lapply(D8_vs_D1_WT_protein_GO_MF$Gene.Count, function(x) list(rep(1,x)) %>% unlist) ,
          main ="WT\n(Day 5 vs Day 1)",
          xlab = "-log10(FDR)",
          x.mar = c(0.7, 0),
          cell.col	= c("#2166AC","white","#44AA99"),
          cell.bounds = c(-0.05,0.05),
          xlab.ticks	= 2,
          y.mar = c(0, 0),
          cex = 1.6,
          key = FALSE,
          sym = TRUE,
          key.lab	= "Coefficient",
          cell.limit = 1,
          cell.outer = 1,
          bar.scale = 0.1,
          space = 0.1)
dev.off()

D5_vs_D1_CK_protein_GO_MF <- read.csv("./input/angelina_code/day 5 vs day 1 CK10 protein GO_MF.csv")
pdf(file = "./output/figure5_c.pdf",width=5.38,height=5.13)
cell.plot(x = setNames(D5_vs_D1_CK_protein_GO_MF$X.log10.padj., D5_vs_D1_CK_protein_GO_MF$Term),
          cells = lapply(D5_vs_D1_CK_protein_GO_MF$Gene.Count, function(x) list(rep(1,x)) %>% unlist) ,
          main ="Tau\n(Day 5 vs Day 1)",
          xlab = "-log10(FDR)",
          x.mar = c(0.7, 0),
          cell.col	= c("#2166AC","white","#44AA99"),
          cell.bounds = c(-0.05,0.05),
          xlab.ticks	= 2,
          y.mar = c(0, 0),
          cex = 1.6,
          key = FALSE,
          sym = TRUE,
          key.lab	= "Coefficient",
          cell.limit = 1,
          cell.outer = 1,
          bar.scale = 0.1,
          space = 0.1)
dev.off()

D8_vs_D1_CL_protein_GO_MF <- read.csv("./input/angelina_code/Day 8 vs Day 1 CL2355 protein GSEA.csv")
pdf(file = "./output/figure5_d.pdf",width=5.38,height=5.13)
cell.plot(x = setNames(D8_vs_D1_CL_protein_GO_MF$X.log10.qvalue., D8_vs_D1_CL_protein_GO_MF$Term),
          cells = lapply(D8_vs_D1_CL_protein_GO_MF$Gene.Count, function(x) list(rep(1,x)) %>% unlist) ,
          main ="AB\n(Day 8 vs Day 1)",
          xlab = "-log10(FDR)",
          x.mar = c(0.7, 0),
          cell.col	= c("#2166AC","white","#44AA99"),
          cell.bounds = c(-0.05,0.05),
          xlab.ticks	= 2,
          y.mar = c(0, 0),
          cex = 1.6,
          key = FALSE,
          sym = TRUE,
          key.lab	= "Coefficient",
          cell.limit = 1,
          cell.outer = 1,
          bar.scale = 0.1,
          space = 0.1)
dev.off()

D8_vs_D1_GL_protein_GO_MF <- read.csv("./input/angelina_code/day 8 vs day 1 GL protein GO_MF.csv")
pdf(file = "./output/figure5_e.pdf",width=5.38,height=5.13)
cell.plot(x = setNames(D8_vs_D1_GL_protein_GO_MF$X.log10.padj., D8_vs_D1_GL_protein_GO_MF$Term),
          cells = lapply(D8_vs_D1_GL_protein_GO_MF$Gene.Count, function(x) list(rep(1,x)) %>% unlist) ,
          main ="Tau:AB\n(Day 8 vs Day 1)",
          xlab = "-log10(FDR)",
          x.mar = c(0.7, 0),
          cell.col	= c("#2166AC","white","#44AA99"),
          cell.bounds = c(-0.05,0.05),
          xlab.ticks	= 2,
          y.mar = c(0, 0),
          cex = 1.6,
          key = FALSE,
          sym = TRUE,
          key.lab	= "Coefficient",
          cell.limit = 1,
          cell.outer = 1,
          bar.scale = 0.1,
          space = 0.1)
dev.off()

