library(tidyverse)
library(ComplexHeatmap)
library(ggpubr)
library(circlize)
library(RColorBrewer)

#load differentially expressed genes from metabolomics analysis day 1
metabolite_degs <- read.csv("./input/metabolomics/metabolites_degs.csv")
sig_metabolites <- filter(metabolite_degs, denominator == "N2")
sig_metabolites$numerator <- gsub("CK10","Tau",sig_metabolites$numerator)
sig_metabolites$numerator <- gsub("CL2355","AB",sig_metabolites$numerator)
sig_metabolites$numerator <- gsub("GL405","Tau:AB",sig_metabolites$numerator)
sig_metabolites$numerator <- factor(sig_metabolites$numerator, levels = c("Tau","AB","Tau:AB"))

#plot fold changes for the comparing of tau/ab expressing worms vs wild-type
p0 <- ggplot(sig_metabolites, aes(x = numerator, y = abs(log(FC, 2)), color = numerator)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha = 0.3, width = 0.2)+
  scale_color_manual(values = c("#EECC66","#6699CC","#EE99AA"))+
  theme_pubr()+
  theme(legend.position = "none")+
  labs(x = "", y = "Absolute log2 fold change")

#transform data into a matrix
matrix <- sig_metabolites %>% filter(numerator=="Tau:AB") %>% dplyr::select(metabolite,FC) %>% mutate(FC = log2(FC)) %>% set_names("metabolite","Tau:AB") %>% 
  left_join(sig_metabolites %>% filter(numerator=="Tau") %>% dplyr::select(metabolite,FC) %>% mutate(FC = log2(FC)) %>% set_names("metabolite","Tau")) %>%
  left_join(sig_metabolites %>% filter(numerator=="AB") %>% dplyr::select(metabolite,FC) %>% mutate(FC = log2(FC)) %>% set_names("metabolite","AB"))

#perform correlation analysis between the genes
cors <- cor(matrix[,2:4])

p1 <- ggplot(matrix, aes(AB,Tau))+
  geom_point(shape = 21, fill = "#44AA99", color = "black", size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", color = "#009988")+
  theme_pubr()+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  geom_label(label = paste0("R = ",round(cors["AB","Tau"],2)), stat = "unique", x = 3, y = -5)+
  labs(x = "log2(AB/WT)", y = "log2(Tau/WT)")

p2 <- ggplot(matrix, aes(`Tau:AB`,Tau))+
  geom_point(shape = 21, fill = "#44AA99", color = "black", size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", color = "#009988")+
  theme_pubr()+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  geom_label(label = paste0("R = ",round(cors["Tau:AB","Tau"],2)), stat = "unique", x = 3, y = -5)+
  labs(x = "log2(Tau:AB/WT)", y = "log2(Tau/WT)")

p3 <- ggplot(matrix, aes(`Tau:AB`,AB))+
  geom_point(shape = 21, fill = "#44AA99", color = "black", size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", color = "#009988")+
  theme_pubr()+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  geom_label(label = paste0("R = ",round(cors["Tau:AB","AB"],2)), stat = "unique", x = 3, y = -5)+
  labs(x = "log2(Tau:AB/WT)", y = "log2(AB/WT)")

#generate a single figure with all the plots
pdf(file = "./output/figure2_ab.pdf",width=10.73,height=2.84)
ggarrange(p0,p1,p2,p3, nrow = 1)
dev.off()

#read metabolites information
features <- read.csv("./input/metabolomics/metabolite_features.csv")[,c(1,5,12)] %>% set_names("variable","type","name") %>% mutate_all(as.character)

#read sample metadata
metabolite_metadata <- read.csv("./input/metabolomics/metabolite_metadata.csv")[,c(1,21,23)] %>% set_names("PARENT_SAMPLE_NAME","strain","time")

#combine count data with metadata
metabolite_counts <- reshape2::melt(read_csv("./input/metabolomics/metabolite_counts.csv")) %>% 
  left_join(features) %>% left_join(metabolite_metadata) %>% 
  filter(time=="D1") %>% group_by(strain,name,type) %>% summarise(value = mean(value))
metabolite_counts$strain <- gsub("CK10","Tau",metabolite_counts$strain)
metabolite_counts$strain <- gsub("CL2355","AB",metabolite_counts$strain)
metabolite_counts$strain <- gsub("GL405","Tau:AB",metabolite_counts$strain)
metabolite_counts$strain <- gsub("N2","WT",metabolite_counts$strain)

#log transform fold changes
sig_metabolites$logFC <- log2(sig_metabolites$FC)

#define up- and down-regulated metabolites
up_metabolites <- sig_metabolites %>% filter(numerator=="Tau:AB"&q.value<0.05&logFC>2&(!grepl("X-",metabolite))) %>% pull(metabolite)
down_metabolites <- sig_metabolites %>% filter(numerator=="Tau:AB"&q.value<0.05&logFC<(-2)&(!grepl("X-",metabolite))) %>% pull(metabolite)

#generate independent matrices of metabolites with increasing and decreasing abudance
down_matrix <- reshape2::acast(metabolite_counts %>% filter(name%in%down_metabolites), name~strain, value.var = "value")
down_matrix <- down_matrix[,c("WT","Tau","AB","Tau:AB")]
up_matrix <- reshape2::acast(metabolite_counts %>% filter(name%in%up_metabolites), name~strain, value.var = "value")
up_matrix <- up_matrix[,c("WT","Tau","AB","Tau:AB")]

#define color palette
col_fun <- colorRamp2(seq(-4,4,1), rev(brewer.pal(9, "RdBu")))

#generate heatmap of the differentially abundant metabolites
pdf(file = "./output/figure2_d.pdf",width=4.4,height=5.76)
Heatmap(as.matrix(down_matrix),
        border = TRUE,
        rect_gp = gpar(col = "black", lwd = 1, lty = 3),
        row_split = sapply(rownames(down_matrix), function(x) unique(metabolite_counts$type[metabolite_counts$name==x])),
        row_title_side = "right",
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        show_row_names = TRUE,
        show_row_dend = FALSE,
        show_column_names = TRUE,
        row_names_side = "left",
        col = col_fun)
dev.off()

pdf(file = "./output/figure2_c.pdf",width=4.4,height=11.54)
Heatmap(as.matrix(up_matrix),
        border = TRUE,
        rect_gp = gpar(col = "black", lwd = 1, lty = 3),
        row_split = sapply(rownames(up_matrix), function(x) unique(metabolite_counts$type[metabolite_counts$name==x])),
        row_title_side = "right",
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        show_row_names = TRUE,
        show_row_dend = FALSE,
        show_column_names = TRUE,
        row_names_side = "left",
        col = col_fun)
dev.off()

