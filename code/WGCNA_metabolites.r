
install.packages("WGCNA")

library(WGCNA)


metabolite_counts <- read.csv("./input/metabolomics/metabolites_degs.csv")


counts <- metabolite_counts
rownames(counts) <- counts$PARENT_SAMPLE_NAME
counts$PARENT_SAMPLE_NAME <- NULL



gsg = goodSamplesGenes(counts, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(counts)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(counts)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  counts = counts[gsg$goodSamples, gsg$goodGenes]
}





sampleTree = hclust(dist(counts), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


nGenes = ncol(counts)
nSamples = nrow(counts)



# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=50, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(counts, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")







###############################################COMPUTE MODULES

counts <- sapply(counts, as.numeric)

allowWGCNAThreads(nThreads = 8)


net = blockwiseModules(counts, power = 9,
                       TOMType = "unsigned", minModuleSize = 7,
                       reassignThreshold = 0, mergeCutHeight = 0.35,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       verbose = 3)











###############################################PLOT DENDOGRAM


sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)










###############################################PLOT ASSOCIATION HEATMAP
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];


metabolite_metadata <- read.csv("./input/angelina_code/metabolite_metadata_wgcna.csv")

metabolite_metadata <- data.frame(lapply(metabolite_metadata, as.factor))

# Convert factor columns to numeric by assigning an integer value to each unique factor level
metabolite_metadata[] <- lapply(metabolite_metadata, function(col) {
  if (is.factor(col)) {
    return(as.numeric(col))
  } else {
    return(col)
  }
})


MEs0 = moduleEigengenes(counts, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, metabolite_metadata, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

mcor = cor(metabolite_metadata, metabolite_metadata, use = "p");



sizeGrWindow(10,10)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(metabolite_metadata),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))



ggsave("metaheatmap.pdf", device = "pdf", width = 20, height = 20, units = "in")
















# Define variable weight containing the weight column of datTrait
weight = as.data.frame(metabolite_metadata$TIME);
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(counts, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");

names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(counts, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");











###############################################PLOT GENE MEMBERSHIP

first <- T
for(i in colnames(MEs) ){
  
  module = gsub("ME", "", i)
  print(module)
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  if(first){
    df <- data.frame(moduleMembership = abs(geneModuleMembership[moduleGenes, column]),
                    module = module)
    first <- F
  }else{
    df2 <- data.frame(moduleMembership = abs(geneModuleMembership[moduleGenes, column]),
                     module = module)
    df <- bind_rows(df,df2)
  }

}

ggplot(df, aes(x=module, y=moduleMembership)) + geom_boxplot()



  
  
  
  ###############################################PLOT scatter plots for modules!!!


meta <- read.csv("./input/metabolomics/metabolite_metadata.csv")

for(i in colnames(MEs) ){
  dt <- data.frame(module = MEs[,i], time = as.numeric(as.factor(meta$TIME)), strain = meta$STRAIN_CN)
  
  a <- ggplot(dt, aes(x=time, y=module, group=time, color=strain)) + facet_grid(. ~ strain) + geom_boxplot() + ggtitle(i) + geom_smooth(method = "lm")
  
  
  
  calc_cor <- function(data){
    cor_test <- cor.test(data$time, data$module)
    cor <- cor_test$estimate
    p.value <- cor_test$p.value
    return(data.frame(cor = cor, p.value = p.value, stringsAsFactors = FALSE))
  }
  
  # Calculate correlation and p-values for each subset of data
  dt_cor <- dt %>% group_by(strain) %>% do(calc_cor(.))
  
  # Add a new column to the data frame that contains the correlation coefficient and p-value as a string
  dt_cor$label <- paste0("r = ", round(dt_cor$cor, 2), "\n",
                         "P = ", formatC(dt_cor$p.value, digits = 3, format = "f"))
  
  # Merge the original data with the new data frame
  dt <- merge(dt, dt_cor, by = "strain")
  i <- gsub("ME", "", i)
  dt$time <- case_match(dt$time, 1 ~ 1, 2 ~ 3, 3 ~ 5, 4 ~ 8)
  # Plot
  b <- ggplot(dt, aes(x=time, y=module)) + 
    facet_grid(. ~ factor(strain, levels = c("WT", "Tau", "Aβ","Tau/Aβ") ) )+ 
    geom_point() +
    geom_smooth(method = "lm") +
    geom_text(aes(label = label), x = Inf, y = Inf, hjust = 1, vjust = 1) + ggtitle(i)
  print(b)
  
  ggsave( paste0("scatterplots/", i, ".pdf"), device = cairo_pdf, width = 8, height = 4 )
  
}













###########################################################EXTRACT GENES FROM MODULES

modules <- data.frame(moduleLabels, moduleColors)
modules$CHEM_ID <- rownames(modules)
modules$CHEM_ID <- gsub("X", "", modules$CHEM_ID)
modules$CHEM_ID <- as.integer(modules$CHEM_ID)

metabolite_features <- read.csv("./input/metabolomics/metabolite_features.csv")


modules <- inner_join(modules, metabolite_features, by = "CHEM_ID")


write_csv(modules, "./output/angelina_code/metabolite_modules_V5.csv")




















