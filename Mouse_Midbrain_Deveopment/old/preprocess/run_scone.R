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

# Make sure to decompress diles first

tpm_counts <- read.table('../rsem_gene_tpm.tab', sep = '\t', 
                         header = TRUE, row.names = 1)

meta <- read.table('../SraRunTable.txt', sep = ',', header = TRUE, row.names = 1)
star_meta <- as.data.frame(t(read.table('../star_meta.tab', sep = '\t', header = TRUE, row.names = 1)))

cells_meta <- rownames(star_meta[star_meta$input_reads > 1000000,])
cells_meta <- cells_meta[cells_meta %in% rownames(star_meta[star_meta$unique_percent > 25,])]
cells_meta <- cells_meta[cells_meta %in% rownames(star_meta[star_meta$map_len/51 > 0.9,])]

tpm_counts <- tpm_counts[,cells_meta]
meta <- meta[cells_meta,]
star_meta <- star_meta[cells_meta,]

neurogenesis_types <- rownames(table(meta$Age))

filtered_genes_tpm <- rownames(tpm_counts)[rowMeans(tpm_counts >= 20) > 0.1]

tpm_filtered <- tpm_counts[filtered_genes_tpm,]

#
cells_e13 <- rownames(meta[meta$Age == 'embryonic day 13.5',])
cells_e18 <- rownames(meta[meta$Age == 'embryonic day 18.5',])
cells_d1 <- rownames(meta[meta$Age == 'postnatal day 1',])
cells_d90 <- rownames(meta[meta$Age == 'postnatal day 90',])

spike_cells = c(cells_e13, cells_e18, cells_d1, cells_d90)
#

spike_genes <- c('EGFP')
for (i in rownames(tpm_filtered)){if (substr(i, 1, 4)=='ERCC'){spike_genes <- c(spike_genes, i)}}

mouse_pos <- c('Cd24a', 'Nrep', 'Sbk1', 'Dcx', 'Carhsp1', 'Mkrn3', 'Cdkn1a', 'S100a11', 
               'Mlkl', 'Fanci', 'Impact', 'Nrsn2', 'Eno2', 'Faim2', 'Tuba4a', 'Mkrn3', 
               'Pitx3', 'Slc6a3', 'Th', 'Nxph4', 'Gad2', 'Aldh1a1', 'Vip', 
               'Dat', 'Slc18a2', 'Ddc', 'Lmx1a', 'Lmx1b', 'Nr4a2', 'Foxa1', 'Foxa2', 'Foxp2', 'En1',
  'Nfia', 'C1ql1', 'Fam19a2', 'Lypd1', 'Car10', 'Slc32a1', 'Wnt7b', 'Ebf3', 'Crhbp', 'Plcl1',
  'Dlk1', 'Slc10a4', 'Cck', 'Zfhx3', 'Drd2', 'Grik3', 'Sncg', 'Npw', 'Unc5d', 'Vwc2l', 'Rprm', 'Gucy2c',
  'Nrp1', 'Kcnab1', 'Grp', 'Lpl', 'Mob3b', 'Id4', 'Gipr', 'Pou2f2', 'Tdh', 'Zfp804a')


positives <- mouse_pos#c(mouse_pos, mpos)
positives <- positives[positives %in% rownames(tpm_filtered)]

mouse_hk <- c('Actb', 'Actg1', 'Anxa2', 'B2m', 'Calm1', 'Calm2', 'Cox6b1', 'Cox7c', 'Cstb',
              'Ctsd', 'Dynll1', 'Eef1a1', 'Eef2', 'Eif4g2', 'Ezr', 'Fau', 'Fth1', 'Gnb2l1',
              'H3f3a', 'Hsp90ab1', 'Hspa8', 'Naca', 'Nfkbia', 'Nme2', 'Pabpc1', 'Prdx1', 'Ptma',
              'Rpl13', 'Rpl13a', 'Rpl29', 'Rpl32', 'Rpl35', 'Rpl36al', 'Rpl37', 'Rpl5', 'Rpl8',
              'Rplp1', 'Rps11', 'Rps12', 'Rps13', 'Rps14', 'Rps16', 'Rps24', 'Rps25', 'Rps27a',
              'Rps5', 'Rps9', 'Tuba1b', 'Ubc','Uqcrh','Ywhaz')


star_meta$log_input_reads <- log10(star_meta$input_reads)
qcNames = c("log_input_reads", "unique_percent")
goodQC <- as.data.frame(star_meta[,qcNames])
goodQC$mapping <- star_meta$map_len/51

batch = factor(meta$Genotype)
##### TPM normalization #####


tpm_filtered <- tpm_filtered[rownames(tpm_filtered)[! rownames(tpm_filtered) %in% spike_genes],]

filtered <- as.matrix(tpm_filtered)

bio = factor(meta$Age)

ppq = scale(goodQC[,apply(goodQC,2,sd) > 0],center = TRUE,scale = TRUE)

negatives <- mouse_hk
positives <- c(mouse_pos)


negatives <- negatives[negatives %in% rownames(tpm_filtered)]
positives <- positives[positives %in% rownames(tpm_filtered)]

my_scone <- SconeExperiment(filtered,
                            qc=ppq, bio = bio, batch=batch,
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


best_norm <- get_normalized(my_scone,1, log=TRUE)

write.table(best_norm, 'scone_norm_tpm.tab', quote = FALSE, sep='\t')

#####
