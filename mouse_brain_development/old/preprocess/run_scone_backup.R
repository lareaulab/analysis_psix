library(ggplot2)
library(gplots)
library(scone)
library(biomaRt)
library(RColorBrewer)
library(pheatmap)
library(destiny, quietly = TRUE)
library(mclust, quietly = TRUE)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(clusterExperiment)
library(slingshot)
library(ggplot2)
library(gplots)
library(scone)
library(biomaRt)
library(pheatmap)

cc <- c(brewer.pal(9, "Set1"))
par(mar=c(1,1,1,1))

# 1. Get the data and metadata.

# tpm_counts <- read.table('Tiklova_tables_spike/rsem_gene_tpm.tab', sep = '\t', header = TRUE, row.names = 1)
tpm_counts <- read.table('/mnt/lareaulab/cfbuenabadn/data_sc_regulation/tiklova/rsem_gene_tpm.tab', sep = '\t', 
                         header = TRUE, row.names = 1)

#star_counts <- read.table('Tiklova_tables_spike/star_counts.tab', sep = '\t', header = TRUE, row.names = 1)
star_counts <- read.table('/mnt/lareaulab/cfbuenabadn/data_sc_regulation/tiklova/star_counts.tab', sep = '\t', 
                          header = TRUE, row.names = 1)

meta <- read.table('/mnt/lareaulab/cfbuenabadn/data_sc_regulation/tiklova/SraRunTable.txt', sep = ',', header = TRUE, row.names = 1)
star_meta <- as.data.frame(t(read.table('/mnt/lareaulab/cfbuenabadn/data_sc_regulation/tiklova/star_meta.tab', sep = '\t', header = TRUE, row.names = 1)))
#meta <- read.table('Tiklova_tables_spike/SraRunTable.txt', sep = ',', header = TRUE, row.names = 1)
#star_meta <- as.data.frame(t(read.table('Tiklova_tables_spike/star_meta.tab', sep = '\t', header = TRUE, row.names = 1)))

cells_meta <- rownames(star_meta[star_meta$input_reads > 1000000,])
cells_meta <- cells_meta[cells_meta %in% rownames(star_meta[star_meta$unique_percent > 25,])]
cells_meta <- cells_meta[cells_meta %in% rownames(star_meta[star_meta$map_len/51 > 0.9,])]

tpm_counts <- tpm_counts[,cells_meta]
star_counts <- star_counts[,cells_meta]
meta <- meta[cells_meta,]
star_meta <- star_meta[cells_meta,]

neurogenesis_types <- rownames(table(meta$Age))

filtered_genes <- rownames(star_counts)[rowMeans(star_counts >= 20) > 0.1]
filtered_genes_tpm <- rownames(tpm_counts)[rowMeans(tpm_counts >= 20) > 0.1]

star_filtered <- star_counts[filtered_genes,]
tpm_filtered <- tpm_counts[filtered_genes_tpm,]

log_counts <- log10(star_counts[filtered_genes,]+1)
log_tpm <- log10(tpm_counts[filtered_genes_tpm,]+1)

color_cell = c()
for (i in meta$Age){
  if (i=='embryonic day 13.5'){color_cell<-c(color_cell, 'red')}
  else if(i=='embryonic day 15.5'){color_cell<-c(color_cell, 'chocolate')}
  else if(i=='embryonic day 18.5'){color_cell<-c(color_cell, 'orange')}
  else if(i=='postnatal day 1'){color_cell<-c(color_cell, 'forestgreen')}
  else if(i=='postnatal day 7'){color_cell<-c(color_cell, 'deepskyblue3')}
  else {color_cell<-c(color_cell, 'darkblue')}
}

pca <- prcomp(t(log_counts))
plot(pca$x[,1], pca$x[,2], col = color_cell, type = 'p',
     pch=16, xlab='PC1', ylab='PC2')

pca <- prcomp(t(log_tpm))
plot(pca$x[,1], pca$x[,2], col = color_cell, type = 'p',
     pch=16, xlab='PC1', ylab='PC2')

#
cells_e13 <- rownames(meta[meta$Age == 'embryonic day 13.5',])
cells_e18 <- rownames(meta[meta$Age == 'embryonic day 18.5',])
cells_d1 <- rownames(meta[meta$Age == 'postnatal day 1',])
cells_d90 <- rownames(meta[meta$Age == 'postnatal day 90',])

spike_cells = c(cells_e13, cells_e18, cells_d1, cells_d90)
#

spike_genes <- c()
for (i in rownames(star_filtered)){if (substr(i, 1, 4)=='ERCC'){spike_genes <- c(spike_genes, i)}}

mouse_pos <- c('Cd24a', 'Nrep', 'Sbk1', 'Dcx', 'Carhsp1', 'Mkrn3', 'Cdkn1a', 'S100a11', 
               'Mlkl', 'Fanci', 'Impact', 'Nrsn2', 'Eno2', 'Faim2', 'Tuba4a', 'Mkrn3', 
               'Pitx3', 'Slc6a3', 'Th', 'Nxph4', 'Gad2', 'Aldh1a1', 'Vip')

mpos <- c('Dat', 'Slc18a2', 'Ddc', 'Lmx1a', 'Lmx1b', 'Nr4a2', 'Foxa1', 'Foxa2', 'Foxp2', 'En1',
  'Nfia', 'C1ql1', 'Fam19a2', 'Lypd1', 'Car10', 'Slc32a1', 'Wnt7b', 'Ebf3', 'Crhbp', 'Plcl1',
  'Dlk1', 'Slc10a4', 'Cck', 'Zfhx3', 'Drd2', 'Grik3', 'Sncg', 'Npw', 'Unc5d', 'Vwc2l', 'Rprm', 'Gucy2c',
  'Nrp1', 'Kcnab1', 'Grp', 'Lpl', 'Mob3b', 'Id4', 'Gipr', 'Pou2f2', 'Tdh', 'Zfp804a')


positives <- c(mouse_pos, mpos)
positives <- positives[positives %in% rownames(star_filtered)]

mouse_hk <- c('Actb', 'Actg1', 'Anxa2', 'B2m', 'Calm1', 'Calm2', 'Cox6b1', 'Cox7c', 'Cstb',
              'Ctsd', 'Dynll1', 'Eef1a1', 'Eef2', 'Eif4g2', 'Ezr', 'Fau', 'Fth1', 'Gnb2l1',
              'H3f3a', 'Hsp90ab1', 'Hspa8', 'Naca', 'Nfkbia', 'Nme2', 'Pabpc1', 'Prdx1', 'Ptma',
              'Rpl13', 'Rpl13a', 'Rpl29', 'Rpl32', 'Rpl35', 'Rpl36al', 'Rpl37', 'Rpl5', 'Rpl8',
              'Rplp1', 'Rps11', 'Rps12', 'Rps13', 'Rps14', 'Rps16', 'Rps24', 'Rps25', 'Rps27a',
              'Rps5', 'Rps9', 'Tuba1b', 'Ubc','Uqcrh','Ywhaz')

#negatives <- c(mouse_hk, spike_genes)
#negatives <- negatives[positives %in% rownames(star_filtered)]
negatives <- spike_genes


star_meta$log_input_reads <- log10(star_meta$input_reads)
qcNames = c("log_input_reads", "unique_percent")
goodQC <- as.data.frame(star_meta[spike_cells,qcNames])
goodQC$mapping <- star_meta[spike_cells,]$map_len/51

#goodQC$genotype <- c(meta[spike_cells,]$Genotype)

batch = factor(meta[spike_cells,]$Genotype)



filtered <- as.matrix(star_filtered[,spike_cells])

#num_reads = quantile(filtered[filtered > 0])[4]
#num_cells = ncol(filtered)
#is_common = rowSums(filtered >= num_reads ) >= num_cells

bio = factor(meta[spike_cells,]$Age)




ppq = scale(goodQC[,apply(goodQC,2,sd) > 0],center = TRUE,scale = TRUE)

star_spike <- star_filtered[,spike_cells]

my_scone <- SconeExperiment(filtered,
                            qc=ppq, bio = bio,
                            negcon_ruv = rownames(star_spike) %in% negatives,
                            poscon = rownames(star_spike) %in% positives
)




EFF_FN = function (ei)
{
  sums = colSums(ei > 0)
  eo = t(t(ei)*sums/mean(sums))
  return(eo)
}

## ----- Scaling Argument -----

scaling=list(none=identity, # Identity - do nothing
             
             eff = EFF_FN, # User-defined function
             
             sum = SUM_FN, # SCONE library wrappers...
             tmm = TMM_FN, 
             uq = UQ_FN,
             fq = FQT_FN,
             deseq = DESEQ_FN)

BiocParallel::register(
  BiocParallel::SerialParam()
) # Register BiocParallel Serial Execution

my_scone <- scone(my_scone,
                  scaling=scaling,
                  #imputation = imputation, impute_args = impute_args,
                  run=TRUE,
                  k_qc=3, k_ruv = 3,
                  adjust_bio="no",
                  eval_kclust = 2:6,
                  return_norm = "in_memory",
                  zero = "postadjust")



color_cell_ercc = c()
for (i in meta[spike_cells,]$Age){
  if (i=='embryonic day 13.5'){color_cell_ercc<-c(color_cell_ercc, 'red')}
  #else if(i=='embryonic day 15.5'){color_cell<-c(color_cell, 'chocolate')}
  else if(i=='embryonic day 18.5'){color_cell_ercc<-c(color_cell_ercc, 'orange')}
  else if(i=='postnatal day 1'){color_cell_ercc<-c(color_cell_ercc, 'forestgreen')}
  #else if(i=='postnatal day 7'){color_cell<-c(color_cell, 'deepskyblue3')}
  else {color_cell_ercc<-c(color_cell_ercc, 'darkblue')}
}



# View Metric Scores
head(get_scores(my_scone))

# View Mean Score Rank
head(get_score_ranks(my_scone))

# Extract normalized data from top method
out_norm = get_normalized(my_scone,
                          method = rownames(get_params(my_scone))[1], log=TRUE)



pc_obj = prcomp(apply(t(get_scores(my_scone)),1,rank),
                center = TRUE,scale = FALSE)
bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)



bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)

points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1)
points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1.5)

points(t(bp_obj[dim(bp_obj)[1],]), pch = 1, col = "blue", cex = 1)
points(t(bp_obj[dim(bp_obj)[1],]), pch = 1, col = "blue", cex = 1.5)

#points(t(bp_obj[rownames(bp_obj) == "none,uq,qc_k=2,no_bio,no_batch",]),
#       pch = 1, col = "blue", cex = 1)
#points(t(bp_obj[rownames(bp_obj) == "none,uq,qc_k=2,no_bio,no_batch",]),
#       pch = 1, col = "blue", cex = 1.5)

arrows(bp_obj[dim(bp_obj)[1],][1],
       bp_obj[dim(bp_obj)[1],][2],
       bp_obj[1,][1],
       bp_obj[1,][2],
       lty = 2, lwd = 2)



good_cells_norm <- get_normalized(my_scone,1, log=TRUE)

best_norm <- good_cells_norm

nrows <- dim(best_norm)[1]
ncols <- dim(best_norm)[2]
tpm <- filtered
normalized <- best_norm

expressed <- rownames(best_norm[rowMedians(best_norm) > 2, ])

fano_factor <- apply(best_norm[expressed,],1,function(x){var(x)/mean(x)})
fano_order <- names(fano_factor[order(fano_factor, decreasing=TRUE)])
top_fano <- fano_order[1:500]


norm_top_fano <- best_norm[top_fano,]



pca <- prcomp(t(as.matrix(norm_top_fano)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = color_cell_ercc, pch=16, asp = 1)



################################
################################



qcNames = c("log_input_reads", "unique_percent")
goodQC <- as.data.frame(star_meta[,qcNames])
goodQC$mapping <- star_meta$map_len/51

batch = factor(meta$Genotype)

mouse_hk <- c('Actb', 'Actg1', 'Anxa2', 'B2m', 'Calm1', 'Calm2', 'Cox6b1', 'Cox7c', 'Cstb',
              'Ctsd', 'Dynll1', 'Eef1a1', 'Eef2', 'Eif4g2', 'Ezr', 'Fau', 'Fth1', 'Gnb2l1',
              'H3f3a', 'Hsp90ab1', 'Hspa8', 'Naca', 'Nfkbia', 'Nme2', 'Pabpc1', 'Prdx1', 'Ptma',
              'Rpl13', 'Rpl13a', 'Rpl29', 'Rpl32', 'Rpl35', 'Rpl36al', 'Rpl37', 'Rpl5', 'Rpl8',
              'Rplp1', 'Rps11', 'Rps12', 'Rps13', 'Rps14', 'Rps16', 'Rps24', 'Rps25', 'Rps27a',
              'Rps5', 'Rps9', 'Tuba1b', 'Ubc','Uqcrh','Ywhaz')

negatives <- mouse_hk

spike_genes <- c(spike_genes, 'EGFP')
star_filtered <- star_filtered[rownames(star_filtered)[! rownames(star_filtered) %in% spike_genes],]

filtered <- as.matrix(star_filtered)

bio = factor(meta$Age)

ppq = scale(goodQC[,apply(goodQC,2,sd) > 0],center = TRUE,scale = TRUE)

negatives <- negatives[negatives %in% rownames(star_filtered)]
positives <- positives[positives %in% rownames(star_filtered)]

my_scone <- SconeExperiment(filtered,
                            qc=ppq, bio = bio,
                            negcon_ruv = rownames(star_filtered) %in% negatives,
                            poscon = rownames(star_filtered) %in% positives
)

EFF_FN = function (ei)
{
  sums = colSums(ei > 0)
  eo = t(t(ei)*sums/mean(sums))
  return(eo)
}

## ----- Scaling Argument -----

scaling=list(none=identity, # Identity - do nothing
             
             eff = EFF_FN, # User-defined function
             
             sum = SUM_FN, # SCONE library wrappers...
             tmm = TMM_FN, 
             uq = UQ_FN,
             fq = FQT_FN,
             deseq = DESEQ_FN)

BiocParallel::register(
  BiocParallel::SerialParam()
) # Register BiocParallel Serial Execution

my_scone <- scone(my_scone,
                  scaling=scaling,
                  #imputation = imputation, impute_args = impute_args,
                  run=TRUE,
                  k_qc=3, k_ruv = 3,
                  adjust_bio="no",
                  eval_kclust = 2:6,
                  return_norm = "in_memory",
                  zero = "postadjust")







# View Metric Scores
head(get_scores(my_scone))

# View Mean Score Rank
head(get_score_ranks(my_scone))

# Extract normalized data from top method
out_norm = get_normalized(my_scone,
                          method = rownames(get_params(my_scone))[1], log=TRUE)



pc_obj = prcomp(apply(t(get_scores(my_scone)),1,rank),
                center = TRUE,scale = FALSE)
bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)



bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)

points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1)
points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1.5)

points(t(bp_obj[dim(bp_obj)[1],]), pch = 1, col = "blue", cex = 1)
points(t(bp_obj[dim(bp_obj)[1],]), pch = 1, col = "blue", cex = 1.5)

#points(t(bp_obj[rownames(bp_obj) == "none,uq,qc_k=2,no_bio,no_batch",]),
#       pch = 1, col = "blue", cex = 1)
#points(t(bp_obj[rownames(bp_obj) == "none,uq,qc_k=2,no_bio,no_batch",]),
#       pch = 1, col = "blue", cex = 1.5)

arrows(bp_obj[dim(bp_obj)[1],][1],
       bp_obj[dim(bp_obj)[1],][2],
       bp_obj[1,][1],
       bp_obj[1,][2],
       lty = 2, lwd = 2)



good_cells_norm <- get_normalized(my_scone,1, log=TRUE)

best_norm <- good_cells_norm

nrows <- dim(best_norm)[1]
ncols <- dim(best_norm)[2]
tpm <- filtered
normalized <- best_norm


##
expressed <- rownames(best_norm[rowMedians(best_norm) > 0.5, ])

fano_factor <- apply(best_norm[expressed,],1,function(x){var(x)/mean(x)})
fano_order <- names(fano_factor[order(fano_factor, decreasing=TRUE)])
top_fano <- fano_order[1:500]

norm_top_fano <- best_norm[top_fano,]

pca <- prcomp(t(as.matrix(norm_top_fano)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = color_cell, pch=16, asp = 1)

write.table(best_norm, 'Tiklova_tables_spike/normalized_star_counts.tab', quote = FALSE, sep='\t')




##### TPM normalization #####


tpm_filtered <- tpm_filtered[rownames(tpm_filtered)[! rownames(tpm_filtered) %in% spike_genes],]

filtered <- as.matrix(tpm_filtered)

bio = factor(meta$Age)

ppq = scale(goodQC[,apply(goodQC,2,sd) > 0],center = TRUE,scale = TRUE)

negatives <- mouse_hk
positives <- c(mouse_pos, mpos)


negatives <- negatives[negatives %in% rownames(tpm_filtered)]
positives <- positives[positives %in% rownames(tpm_filtered)]

my_scone <- SconeExperiment(filtered,
                            qc=ppq, bio = bio,
                            negcon_ruv = rownames(tpm_filtered) %in% negatives,
                            poscon = rownames(tpm_filtered) %in% positives
)

EFF_FN = function (ei)
{
  sums = colSums(ei > 0)
  eo = t(t(ei)*sums/mean(sums))
  return(eo)
}

## ----- Scaling Argument -----

scaling=list(none=identity, # Identity - do nothing
             
             eff = EFF_FN, # User-defined function
             
             sum = SUM_FN, # SCONE library wrappers...
             tmm = TMM_FN, 
             uq = UQ_FN,
             fq = FQT_FN,
             deseq = DESEQ_FN)

BiocParallel::register(
  BiocParallel::SerialParam()
) # Register BiocParallel Serial Execution

my_scone <- scone(my_scone,
                  scaling=scaling,
                  #imputation = imputation, impute_args = impute_args,
                  run=TRUE,
                  k_qc=3, k_ruv = 3,
                  adjust_bio="no",
                  eval_kclust = 2:6,
                  return_norm = "in_memory",
                  zero = "postadjust")







# View Metric Scores
head(get_scores(my_scone))

# View Mean Score Rank
head(get_score_ranks(my_scone))

# Extract normalized data from top method
out_norm = get_normalized(my_scone,
                          method = rownames(get_params(my_scone))[1], log=TRUE)



pc_obj = prcomp(apply(t(get_scores(my_scone)),1,rank),
                center = TRUE,scale = FALSE)
bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)



bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)

points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1)
points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1.5)

points(t(bp_obj[dim(bp_obj)[1],]), pch = 1, col = "blue", cex = 1)
points(t(bp_obj[dim(bp_obj)[1],]), pch = 1, col = "blue", cex = 1.5)

#points(t(bp_obj[rownames(bp_obj) == "none,uq,qc_k=2,no_bio,no_batch",]),
#       pch = 1, col = "blue", cex = 1)
#points(t(bp_obj[rownames(bp_obj) == "none,uq,qc_k=2,no_bio,no_batch",]),
#       pch = 1, col = "blue", cex = 1.5)

arrows(bp_obj[dim(bp_obj)[1],][1],
       bp_obj[dim(bp_obj)[1],][2],
       bp_obj[1,][1],
       bp_obj[1,][2],
       lty = 2, lwd = 2)



good_cells_norm <- get_normalized(my_scone,1, log=TRUE)

best_norm <- good_cells_norm

nrows <- dim(best_norm)[1]
ncols <- dim(best_norm)[2]
tpm <- filtered
normalized <- best_norm


##
expressed <- rownames(best_norm[rowMedians(best_norm) > 0.5, ])

fano_factor <- apply(best_norm[expressed,],1,function(x){var(x)/mean(x)})
fano_order <- names(fano_factor[order(fano_factor, decreasing=TRUE)])
top_fano <- fano_order[1:500]

norm_top_fano <- best_norm[top_fano,]

pca <- prcomp(t(as.matrix(norm_top_fano)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = color_cell, pch=16, asp = 1)

# write.table(best_norm, 'Tiklova_tables_spike/normalized_tpm_counts.tab', quote = FALSE, sep='\t') # original name
write.table(best_norm, 'Tiklova_tables_spike/scone_norm_tpm.tab', quote = FALSE, sep='\t')

#####

library(Rtsne)

norm_tsne <- Rtsne::Rtsne(best_norm[top_fano,])
