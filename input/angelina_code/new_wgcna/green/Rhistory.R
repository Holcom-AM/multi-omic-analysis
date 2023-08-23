# PID of current job: 84062
mSet<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-c("HMDB0000755","HMDB0000134","HMDB0031518","HMDB0002166","HMDB0002088","HMDB0004080","HMDB0004666","HMDB0000176","HMDB0011567","HMDB0011537","HMDB0011538","HMDB0011568","HMDB11578","HMDB0011569","HMDB0011561","HMDB0000912","HMDB0011577","HMDB0001851","HMDB0011565","HMDB0011534","HMDB11580","HMDB0011550")
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "hmdb");
mSet<-CreateMappingResultTable(mSet)
mSet<-Setup.HMDBReferenceMetabolome(mSet, "background list based on pubchem ids.txt");
mSet<-SetMetabolomeFilter(mSet, T);
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_10_", "net", "png", 72, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_10_", "png", 72, width=NA)
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_11_", "net", "png", 72, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_11_", "png", 72, width=NA)
mSet<-SaveTransformedData(mSet)
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_12_", "net", "png", 72, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_12_", "png", 72, width=NA)
mSet<-CreateMappingResultTable(mSet)
mSet<-Setup.HMDBReferenceMetabolome(mSet, "background list based on pubchem ids.txt");
mSet<-SetMetabolomeFilter(mSet, T);
mSet<-SetCurrentMsetLib(mSet, "kegg_pathway", 2);
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_13_", "net", "png", 72, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_13_", "png", 72, width=NA)
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_14_", "net", "png", 72, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_14_", "png", 72, width=NA)
mSet<-SaveTransformedData(mSet)
