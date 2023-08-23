library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)

#read metabolomics analysis
data <- read_rds("./input/metabolomics/pairwise_comparisons.rds")
data <- rbind(data[grepl("N2_D1",data$control)&grepl("N2",data$mutant),] %>% filter(q<0.05) %>% mutate(strain = "WT"),
      data[grepl("CK10_D1",data$control)&grepl("CK10",data$mutant),] %>% filter(q<0.05) %>% mutate(strain = "Tau"),
      data[grepl("CL2355_D1",data$control)&grepl("CL2355",data$mutant),] %>% filter(q<0.05) %>% mutate(strain = "AB"),
      data[grepl("GL405_D1",data$control)&grepl("GL405",data$mutant),] %>% filter(q<0.05) %>% mutate(strain = "Tau:AB"))
data$day <- sapply(data$mutant, function(x) gsub("D","",strsplit(x, split = "_")[[1]][2]) %>% as.numeric)
data$dir <- ifelse(data$logfc>0,"Increase","Decrease")
data <- data %>% group_by(strain,day,dir) %>% summarise(n = n()) %>% mutate(day = as.character(day))
data$strain <- factor(data$strain, levels = c("WT","Tau","AB","Tau:AB"))

#plot number of significant metabolites
p1 <- ggplot(data, aes(fill=dir, y=n, x=day)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x="Group", y="Count") +
  facet_wrap(.~strain, nrow = 1)+
  scale_fill_manual(values=c(brewer.pal(9, "RdBu")[8], brewer.pal(9, "RdBu")[2]))+
  theme_pubr()+
  theme(plot.background=element_blank(), strip.background = element_blank(), legend.position = "top") +
  labs(x = "Compared to day 1", y = "Number of metabolites\n(q-value < 0.05)", fill = "")

#read metadata and modules
features <- read.csv("./input/angelina_code/metabolite_features.csv")[,c(1,5,12)] %>% set_names("variable","type","name") %>% mutate_all(as.character)
metabolite_metadata <- read.csv("./input/angelina_code/metabolite_metadata.csv")[,c(1,21,23)] %>% set_names("PARENT_SAMPLE_NAME","strain","time")
green_module <- read_csv("./input/angelina_code/new_wgcna/green/green module.csv")

#merge abundance with metadata
metabolite_counts <- reshape2::melt(read_csv("./input/angelina_code/metabolite_counts.csv")) %>% 
  left_join(features) %>% left_join(metabolite_metadata) %>% 
  filter(name%in%green_module$PLOT_NAME) %>% group_by(strain,time,PARENT_SAMPLE_NAME) %>% summarise(mean = mean(value))
metabolite_counts$strain <- gsub("CK10","Tau",metabolite_counts$strain)
metabolite_counts$strain <- gsub("CL2355","AB",metabolite_counts$strain)
metabolite_counts$strain <- gsub("GL405","Tau:AB",metabolite_counts$strain)
metabolite_counts$strain <- gsub("N2","WT",metabolite_counts$strain)
metabolite_counts$time <- sapply(metabolite_counts$time, function(x) gsub("D","",x) %>% as.numeric)
metabolite_counts$strain <- factor(metabolite_counts$strain, levels = c("WT","Tau","AB","Tau:AB"))

#plot average abundance levels for genes in green module
p2 <- ggplot(metabolite_counts, aes(time, mean))+
  geom_point(shape = 21, fill = "#44AA99", color = "black", size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", color = "#009988")+
  facet_wrap(.~strain, nrow = 1)+
  theme_pubr()+
  theme(plot.background=element_blank(), strip.background = element_blank(), legend.position = "top") +
  scale_x_continuous(breaks = c(1,3,5,8))+
  labs(x = "Day", y = "Log transformed abundance")

#merge previous figures
pdf(file = "./output/figure3_ab.pdf",width=6.59,height=5.53)
ggarrange(p1,p2, nrow = 2, heights = c(1,1))
dev.off()

#read abundance levels for metabolites
metabolite_counts <- reshape2::melt(read_csv("./input/angelina_code/metabolite_counts.csv")) %>% 
  left_join(features) %>% left_join(metabolite_metadata) %>% 
  filter(name%in%green_module$PLOT_NAME) %>% group_by(strain,type,time,name) %>% summarise(value = mean(value))
metabolite_counts$strain <- gsub("CK10","2Tau",metabolite_counts$strain)
metabolite_counts$strain <- gsub("CL2355","3AB",metabolite_counts$strain)
metabolite_counts$strain <- gsub("GL405","4Tau:AB",metabolite_counts$strain)
metabolite_counts$strain <- gsub("N2","1WT",metabolite_counts$strain)
metabolite_counts$time <- sapply(metabolite_counts$time, function(x) gsub("D","",x) %>% as.character)
metabolite_counts <- metabolite_counts[!grepl("X-",metabolite_counts$name),]

matrix <- reshape2::acast(metabolite_counts,name~interaction(time,strain), value.var = "value")
col_fun <- colorRamp2(seq(-2,2,0.5), rev(brewer.pal(9, "RdBu")))

#plot metabolite abundance heatmap
pdf(file = "./output/figure3_c.pdf",width=6.72,height=5.89)
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
