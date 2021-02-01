library(SummarizedExperiment)
library(SingleCellExperiment)
library(clusterExperiment)
library(slingshot)
library(destiny, quietly = TRUE)
library(mclust, quietly = TRUE)
library(RColorBrewer)


rd <- read.table('../data_sc_regulation/tiklova/rd.tab', sep = '\t', header = TRUE, row.names = 1)
tpm_norm_counts <- read.table('../data_sc_regulation/tiklova/scone_norm_tpm.tab', sep = '\t', header = TRUE, row.names = 1)
sc_meta <- read.table('../data_sc_regulation/tiklova/SraRunTable.txt', sep = ',', header = TRUE, row.names = 1)
mrna <- read.table('../data_sc_regulation/tiklova/mrna_per_event_adj.tab', sep = '\t', header = TRUE, row.names = 1)

rd <- rd[colnames(mrna),]
tpm_norm_counts <- tpm_norm_counts[,colnames(mrna)]
sc_meta <- sc_meta[colnames(mrna),]

neurogenesis <- SingleCellExperiment(assays=list(normalized=tpm_norm_counts),colData=sc_meta)
#reducedDims(neurogenesis) <- SimpleList(PCA = as.matrix(rd[,c("PC_1", "PC_2", "PC_3", "PC_4", "PC_5")]), 
#                                        tSNE = as.matrix(rd[,c("tSNE_1", "tSNE_2"),]),
#                                        UMAP = as.matrix(rd[,c("UMAP_1", "UMAP_2", "UMAP_3", "UMAP_4", "UMAP_5")]))

reducedDims(neurogenesis) <- SimpleList(PCA = as.matrix(rd[,c("PC_1", "PC_2")]), 
                                        tSNE = as.matrix(rd[,c("tSNE_1", "tSNE_2"),]),
                                        UMAP = as.matrix(rd[,c("UMAP_1", "UMAP_2", "UMAP_3", "UMAP_4", "UMAP_5")]))

#cl1 <- Mclust(rd[,c("PC_1", "PC_2", "PC_3", "PC_4", "PC_5")])$classification
cl1 <- Mclust(rd[,c("PC_1", "PC_2")])$classification
colData(neurogenesis)$GMM <- cl1

plot(rd[,c("PC_1", "PC_2")], col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)




#sce <- slingshot(neurogenesis,  reducedDim = 'PCA',  start.clus = 'embryonic day 13.5', 
#                 end.clus = 'postnatal day 90', clusterLabels = 'Age')


#sce <- slingshot(neurogenesis,  reducedDim = 'PCA',  start.clus = '1', end.clust='9', clusterLabels = 'GMM')
sce <- slingshot(neurogenesis,  reducedDim = 'PCA',  start.clus = '1', clusterLabels = 'GMM')


colors <- colorRampPalette(brewer.pal(9,'Spectral')[-6])(100)
plot(reducedDims(sce)$PCA[,c('PC_1', 'PC_2')], col = colors[cut(sce$slingPseudotime_1,breaks=100)], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2)



pseudo_curve <- data.frame(attributes(SlingshotDataSet(sce))$curves$curve1$s)
pseudo_curve$lineage_1_order <- attributes(SlingshotDataSet(sce))$curves$curve1$ord
pseudo_curve$lineage_1_pseudotime <- sce$slingPseudotime_1

colnames(pseudo_curve) <- c('lineage_1_PC_1', 'lineage_1_PC_2', 'lineage_1_order', 'lineage_1_pseudotime')

pseudo_curve$lineage_2_PC_1 <- attributes(SlingshotDataSet(sce))$curves$curve2$s[,'PC_1']
pseudo_curve$lineage_2_PC_2 <- attributes(SlingshotDataSet(sce))$curves$curve2$s[,'PC_2']
pseudo_curve$lineage_2_order <- attributes(SlingshotDataSet(sce))$curves$curve1$ord
pseudo_curve$lineage_2_pseudotime <- sce$slingPseudotime_2


write.table(pseudo_curve, '../data_sc_regulation/tiklova/pseudotime.tab', quote = FALSE, sep='\t')
