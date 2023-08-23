library(tidyverse)
library(eulerr)
library(ggpubr)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

#read proteomics analysis results
data <- read_rds("./input/insolublome/pairwise.rds")
aging <- data %>% filter(control=="N2_D1"&mutant=="N2_D8"&q<0.05) %>% pull(id)
mutant <- data %>% filter(control=="N2_D1"&mutant=="GL405_D1"&q<0.05) %>% pull(id)

#generate list with the proteins
list1 <- euler(c(A = length(aging), B = length(mutant),  
                 "A&B" = length(intersect(aging,mutant))))

#plot venn diagram
p1 <- plot(list1,
     fills = c("grey", "pink"),
     edges = TRUE,
     fontsize = 60,
     quantities = list(fontsize = 20))

#read metabolomics analysis results
data <- read_rds("./input/metabolomics/pairwise_comparisons.rds")
aging <- data %>% filter(control=="N2_D1"&mutant=="N2_D8"&q<0.05) %>% pull(metabolite)
mutant <- data %>% filter(control=="N2_D1"&mutant=="GL405_D1"&q<0.05) %>% pull(metabolite)

#generate list with metabolites
list2 <- euler(c(A = length(aging), B = length(mutant),  
                "A&B" = length(intersect(aging,mutant))))

#plot venn diagram
p2 <- plot(list2,
     fills = c("grey", "pink"),
     edges = TRUE,
     fontsize = 60,
     quantities = list(fontsize = 20))

#merge previous plot in a single figure
pdf(file = "./output/figure6_ab.pdf",width=5.38,height=5.13)
ggarrange(p1,NULL,p2, ncol = 1, heights = c(1,0.2,1))
dev.off()

#read metabolite abundance
features <- read.csv("./input/angelina_code/metabolite_features.csv")[,c(1,5,12)] %>% set_names("variable","type","name") %>% mutate_all(as.character)
metabolite_metadata <- read.csv("./input/angelina_code/metabolite_metadata.csv")[,c(1,21,23)] %>% set_names("PARENT_SAMPLE_NAME","strain","time")
overlap <- intersect(aging,mutant)
metabolite_counts <- reshape2::melt(read_csv("./input/angelina_code/metabolite_counts.csv")) %>% 
  left_join(features) %>% left_join(metabolite_metadata) %>% 
  filter(name%in%overlap) %>% group_by(strain,time,type,name) %>% summarise(mean = mean(value))
metabolite_counts$strain <- gsub("GL405","2Tau:AB",metabolite_counts$strain)
metabolite_counts$strain <- gsub("N2","1WT",metabolite_counts$strain)
metabolite_counts$time <- sapply(metabolite_counts$time, function(x) gsub("D","",x) %>% as.numeric)
metabolite_counts <- metabolite_counts %>% filter((time==1&strain=="1WT")|(time==8&strain=="1WT")|(time==1&strain=="2Tau:AB"))

matrix <- reshape2::acast(metabolite_counts,name~interaction(time,strain), value.var = "mean")
col_fun <- colorRamp2(seq(-2,2,0.5), rev(brewer.pal(9, "RdBu")))

#plot metabolite abundance
pdf(file = "./output/figure6_c.pdf",width=4.25,height=9.64)
Heatmap(as.matrix(matrix),
        border = TRUE,
        rect_gp = gpar(col = "black", lwd = 1, lty = 3),
        row_split = sapply(rownames(matrix), function(x) unique(metabolite_counts$type[metabolite_counts$name==x])),
        column_split = sapply(colnames(matrix), function(x) strsplit(x, split = "\\.")[[1]][2]),
        row_title_side = "right",
        column_labels = sapply(colnames(matrix), function(x) strsplit(x, split = "\\.")[[1]][1]), 
        column_names_rot = 0,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        show_row_names = TRUE,
        show_row_dend = FALSE,
        show_column_names = TRUE,
        row_names_side = "left",
        col = col_fun)
dev.off()




