library(gplots)
library(ggplot2)

downreg_GL_Tissue<-read.csv("./input/angelina_code/downregulated genes tissue enrichment day 1 GL vs WT.csv")
downreg_GL_GO<-read.csv("./input/angelina_code/downregulated genes GO enrichment day 1 GL vs WT.csv")
upeg_GL_Tissue<-read.csv("./input/angelina_code/upregulated genes tissue enrichment day 1 GL vs WT.csv")
upreg_GL_GO<-read.csv("./input/angelina_code/upregulated genes GO enrichment day 1 GL vs WT.csv")


downreg_GL_Tissue <- downreg_GL_Tissue[order(downreg_GL_Tissue$Gene.Count, decreasing = TRUE), ]

# Create color palette
my_palette <- colorRampPalette(c("white", "black"))(max(downreg_GL_Tissue$Gene.Count))

# Create barplot
ggplot(downreg_GL_Tissue, aes(x = X.log10.qvalue., y = Term, fill = Gene.Count)) +
  geom_bar(stat = "identity", width = 0.3, position = "identity") +
  scale_fill_gradientn(colours = my_palette, 
                       na.value = "white",
                       limits = c(0, max(downreg_GL_Tissue$Gene.Count))) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "-log10(qvalue)", y = "Terms",
       title = "") +
  theme_classic()

###################################################

# Create reordered factor levels for the Term variable based on X.log10.qvalue.
downreg_GL_GO$Term <- reorder(downreg_GL_GO$Term, downreg_GL_GO$X.log10.qvalue., FUN = function(x) -log10(x))

# Create color palette
my_palette2 <- colorRampPalette(c("white", "black"))(max(downreg_GL_GO$Gene.Count))

# Create barplot
ggplot(downreg_GL_GO, aes(x = X.log10.qvalue., y = Term, fill = Gene.Count)) +
  geom_bar(stat = "identity", width = 0.5, position = "identity") +
  scale_fill_gradientn(colours = my_palette2, 
                       na.value = "white",
                       limits = c(0, max(downreg_GL_GO$Gene.Count))) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "-log10(qvalue)", y = "Terms",
       title = "") +
  theme_classic()


###################################################
# Create reordered factor levels for the Term variable based on X.log10.qvalue.
upeg_GL_Tissue$Term <- reorder(upeg_GL_Tissue$Term, upeg_GL_Tissue$X.log10.qvalue., FUN = function(x) -log10(x))

# Create color palette
my_palette3 <- colorRampPalette(c("white", "black"))(max(upeg_GL_Tissue$Gene.Count))

# Create barplot
ggplot(upeg_GL_Tissue, aes(x = X.log10.qvalue., y = Term, fill = Gene.Count)) +
  geom_bar(stat = "identity", width = 0.5, position = "identity") +
  scale_fill_gradientn(colours = my_palette3, 
                       na.value = "white",
                       limits = c(0, max(upeg_GL_Tissue$Gene.Count))) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "-log10(qvalue)", y = "Terms",
       title = "") +
  theme_classic()



###################################################
# Create reordered factor levels for the Term variable based on X.log10.qvalue.
upreg_GL_GO$Term <- reorder(upreg_GL_GO$Term, upreg_GL_GO$X.log10.qvalue., FUN = function(x) -log10(x))

# Create color palette
my_palette4 <- colorRampPalette(c("white", "black"))(max(upreg_GL_GO$Gene.Count))

# Create barplot
ggplot(upreg_GL_GO, aes(x = X.log10.qvalue., y = Term, fill = Gene.Count)) +
  geom_bar(stat = "identity", width = 0.5, position = "identity") +
  scale_fill_gradientn(colours = my_palette4, 
                       na.value = "white",
                       limits = c(0, max(upreg_GL_GO$Gene.Count))) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "-log10(qvalue)", y = "Terms",
       title = "") +
  theme_classic()