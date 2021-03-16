# This script is now unsupported. Keeping it for references.

# Expression data for Song et al
library(ggplot2)
library(gplots)
library(scone)
library(biomaRt)
library(RColorBrewer)
library(pheatmap)
library(slingshot)

cc <- c(brewer.pal(9, "Set1"))
par(mar=c(1,1,1,1))

#dir <- '/data/yosef2/users/cfbuenabadn/Neurogenesis/Botvinnik'
counts <- read.table('rsemTPM.tab', header = TRUE, row.names = 1)
meta <- read.table('metaData.tab', sep = '\t', 
                   header = TRUE, row.names = 1)

# This selects only single cells
counts <- counts[,rownames(meta)]

meta$Percent_unique <-  meta$Percent_unique * 100
meta$Percent_unique_annotated <- meta$Percent_unique_annotated * 100
meta$Percent_annotated_v_unique <- meta$Percent_annotated_v_unique * 100

filtered <- counts[rowSums(counts > 20) > 20,]
log_counts <- log1p(filtered)
pca <- prcomp(t(log_counts))


#human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
#hg38 <- getBM(c("ensembl_gene_id", "hgnc_symbol"),mart=human)

hg38 <- readRDS('hg38_translate.rds')
   
SRSF <- c('SRSF1', 'SRSF2', 'SRSF3', 'SRSF4', 'SRSF5', 'SRSF6', 
          'SRSF8', 'SRSF7', 'SRSF9', 'SRSF10', 'SRSF11', 'SRSF12')
ensembl_SRSF <- SRSF
for (i in 1:12) {ensembl_SRSF[i] <- hg38$ensembl_gene_id[hg38$hgnc_symbol == SRSF[i]][1]}



################ Trying SCONE

ralign_human = meta$Percent_unique
batch = factor(meta$cell_type_s)

o = order(ralign_human)[order(batch[order(ralign_human)])] # Order by batch, then value

barplot(ralign_human[o], col=cc[batch][o], 
        border=cc[batch][o], main="Percentage of reads mapped")
legend("bottomleft", legend=levels(batch), fill=cc,cex=0.4)

#

ralign_human = meta$Percent_unique
batch = factor(meta$LibraryLayout_s)

o = order(ralign_human)[order(batch[order(ralign_human)])] # Order by batch, then value

barplot(ralign_human[o], col=cc[batch][o], 
        border=cc[batch][o], main="Percentage of reads mapped")
legend("bottomleft", legend=levels(batch), fill=cc,cex=0.4)

#

ralign_human = meta$Percent_unique
batch = factor(paste(meta$cell_type_s,meta$LibraryLayout_s, sep='-'))

o = order(ralign_human)[order(batch[order(ralign_human)])] # Order by batch, then value

barplot(ralign_human[o], col=cc[batch][o], 
        border=cc[batch][o], main="Percentage of reads mapped")
legend("bottomleft", legend=levels(batch), fill=cc,cex=0.4)

#############################


# Barplot of total read number
nreads = meta$Total_reads
batch = factor(paste(meta$cell_type_s,meta$LibraryLayout_s, sep='-'))

o = order(nreads)[order(batch[order(nreads)])] # Order by batch, then value

barplot(nreads[o], col=cc[batch][o], 
        border=cc[batch][o], main="Total number of reads")
legend("topright", legend=levels(batch), fill=cc, cex=0.4)

#############################

qcNames = c("MBases_l", "Total_reads", "Uniquely_mapped_reads", "unique_annotated_reads",    
            "Percent_unique", "Percent_unique_annotated", "Percent_annotated_v_unique")

qc <- as.matrix(meta[,qcNames])

## ----- PCA of QC matrix -----
qpc = prcomp(qc,center = TRUE,scale. = TRUE)
barplot((qpc$sdev^2)/sum(qpc$sdev^2), border="gray", 
        xlab="PC", ylab="Proportion of Variance", main="Quality PCA")



# Barplot of PC1 of the QC matrix
qc1 = as.vector(qpc$x[,1])
o = order(qc1)[order(batch[order(qc1)])]

barplot(qc1[o], col=cc[batch][o], 
        border=cc[batch][o], main="Quality PC1")
legend("bottomright", legend=levels(batch), 
       fill=cc, cex=0.8)


##############################

# Extract Housekeeping Genes
library(scRNAseq)
data(housekeeping)
hk = housekeeping$V1

hk <- intersect(hg38$ensembl_gene_id[hg38$hgnc_symbol %in% hk], rownames(filtered))

########################################


hg_obs = rowMeans(log10(as.matrix(filtered[hk,])+1))
drop_outs = as.matrix(filtered[hk,]) == 0


ref.glms = list()
for (si in 1:dim(drop_outs)[2]){
  fit = glm(cbind(drop_outs[,si],1 - drop_outs[,si]) ~ hg_obs,
            family=binomial(logit))
  ref.glms[[si]] = fit$coefficients
}

###

par(mfrow=c(1,2))

# Plot Failure Curves and Calculate AUC
plot(NULL, main = "False Negative Rate Curves",
     ylim = c(0,1),xlim = c(0,6), 
     ylab = "Failure Probability", xlab = "Mean log10 Expression")
x = (0:60)/10
AUC = NULL
for(si in 1:ncol(filtered)){
  y = 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1)
  AUC[si] = sum(y)/10
  lines(x, 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1),
        type = 'l', lwd = 2, col = cc[batch][si])
}


# Barplot of FNR AUC
o = order(AUC)[order(batch[order(AUC)])]

barplot(AUC[o], col=cc[batch][o], border=cc[batch][o], main="FNR AUC")
legend("topright", legend=levels(batch), fill=cc, cex=0.4)

##########

yeo_counts <- as.matrix(filtered)

num_reads = quantile(yeo_counts[yeo_counts > 0])[4]
num_cells = ncol(yeo_counts)
is_common = rowSums(yeo_counts >= num_reads ) >= num_cells

# Metric-based Filtering
mfilt = metric_sample_filter(yeo_counts,
                             nreads = meta$Total_reads,
                             ralign = meta$Percent_unique,
                             gene_filter = is_common,
                             pos_controls = rownames(yeo_counts) %in% hk,
                             
                             zcut = 3, mixture = FALSE,
                             plot = TRUE)

# Simplify to a single logical
mfilt = !apply(simplify2array(mfilt[!is.na(mfilt)]),1,any)

logCounts = log1p(counts[,mfilt])

par(mfrow=c(1,1))
layout(c(1,1))


######################################


hist(meta$Percent_unique, breaks=0:100)
# Hard threshold
abline(v = 15, col = "yellow", lwd = 2)
# 3 (zcut) standard deviations below the mean ralign value
abline(v = mean(meta$Percent_unique) - 3*sd(meta$Percent_unique), col = "green", lwd = 2)
# 3 (zcut) MADs below the median ralign value
abline(v = median(meta$Percent_unique) - 3*mad(meta$Percent_unique), col = "red", lwd = 2)
# Sufficient threshold
abline(v = NULL, col = "grey", lwd = 2)

# Final threshold is the minimum of 
# 1) the sufficient threshold and 
# 2) the max of all others
thresh = min(NULL,
             max(c(15,mean(meta$Percent_unique) - 3*sd(meta$Percent_unique),
                   median(meta$Percent_unique) - 3*mad(meta$Percent_unique))))
abline(v = thresh, col = "blue", lwd = 2, lty = 2)

legend("topleft",legend = c("Hard","SD","MAD","Sufficient","Final"),
       lwd = 2, col = c("yellow","green","red","grey","blue"),
       lty = c(1,1,1,1,2), cex = .5)



#goodCounts = counts_male[,mfilt]

# This is a modification I made...
mfilt <- !(meta$Percent_unique < thresh)

goodCounts = yeo_counts

goodMeta = meta
goodQC = qc


# Final Gene Filtering: Highly expressed in at least 5 cells
num_reads = quantile(yeo_counts[yeo_counts > 0])[4]
num_cells = 10
is_quality = rowSums(yeo_counts >= num_reads ) >= num_cells

###

### Check for SRSF

##### Creating SCONE experiment

# Expression Data (Required)
expr = goodCounts[is_quality,]
log_expr = logCounts[is_quality,]

# Biological Origin - Variation to be preserved (Optional)
bio = factor(goodMeta$cell_type_s)

# Processed Alignment Metrics - Variation to be removed (Optional)
ppq = scale(goodQC[,apply(goodQC,2,sd) > 0],center = TRUE,scale = TRUE)




hk <- scan("hk_filtered.txt", what="", sep="\n")
poscon <- scan("pos_controls_filtered.txt", what="", sep="\n")



##############################

# Positive Control Genes - Prior knowledge of DE (Optional)
poscon = intersect(rownames(goodCounts), poscon)

# Negative Control Genes - Uniformly expressed transcripts (Optional)
negcon = intersect(rownames(goodCounts), hk)

#fnr_out = estimate_ziber(x = goodCounts, bulk_model = TRUE,
#                         pos_controls = rownames(goodCounts) %in% poscon,
#                         maxiter = 100)

#imputation=list(none=impute_null, # No imputation
#                expect=impute_expectation) # Replace zeroes

## ----- Imputation Function Arguments -----
# accessible by functions in imputation list argument
#impute_args = list(p_nodrop = fnr_out$p_nodrop, mu = exp(fnr_out$Alpha[1,]))



# Creating a SconeExperiment Object
my_scone <- SconeExperiment(goodCounts,
                            qc=ppq, bio = bio,
                            negcon_ruv = rownames(goodCounts) %in% negcon,
                            poscon = rownames(goodCounts) %in% poscon
)



####################################

## ----- User-defined function: Dividing by number of detected genes -----

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

########################################

my_scone <- scone(my_scone,
                  scaling=scaling,
                  #imputation = imputation, impute_args = impute_args,
                  k_qc=3, k_ruv = 3,
                  adjust_bio="no",
                  run=FALSE)

head(get_params(my_scone))

apply(get_params(my_scone),2,unique)


is_screened = ((get_params(my_scone)$imputation_method == "expect") &
                 (get_params(my_scone)$scaling_method %in% c("none",
                                                             "eff")))

my_scone = select_methods(my_scone,
                          rownames(get_params(my_scone))[!is_screened ])


BiocParallel::register(
  BiocParallel::SerialParam()
) # Register BiocParallel Serial Execution

my_scone <- scone(my_scone,
                  scaling=scaling,
                  #imputation = imputation, impute_args = impute_args,
                  run=TRUE,
                  eval_kclust = 2:6,
                  return_norm = "in_memory",
                  zero = "postadjust")





# View Metric Scores
head(get_scores(my_scone))

# View Mean Score Rank
head(get_score_ranks(my_scone))

# Extract normalized data from top method
out_norm = get_normalized(my_scone,
                          method = rownames(get_params(my_scone))[1])



pc_obj = prcomp(apply(t(get_scores(my_scone)),1,rank),
                center = TRUE,scale = FALSE)
bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)



bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)

points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1)
points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1.5)

points(t(bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",]),
       pch = 1, col = "blue", cex = 1)
points(t(bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",]),
       pch = 1, col = "blue", cex = 1.5)

arrows(bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",][1],
       bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",][2],
       bp_obj[1,][1],
       bp_obj[1,][2],
       lty = 2, lwd = 2)




best_norm <- get_normalized(my_scone,1, log=TRUE)

pca_norm <- prcomp(t(as.matrix(best_norm)))
cells <- goodMeta[rownames(pca_norm$x),]$cell_type
colors <- c()
for (i in 1:length(cells)) {
  if (cells[i]=='iPSC') {
    colors <- c(colors,'red')
  } else if (cells[i]=='NPC') {
    colors <- c(colors,'forestgreen')
  } else if (cells[i]=='MN') {
    colors <- c(colors,'blue')
  }
}


plot(pca_norm$x[,1], pca_norm$x[,2], col = colors, type = 'p', 
     pch=16, xlab='PC1', ylab='PC2')
title('Best normalization(counts)')
legend(-45000,-27000, c('mES2i', 'mES', 'Epi', 'MN'), pch=16, 
       col=c('red', 'orange','forestgreen','blue'))




library(slingshot)

pca <- prcomp(t(best_norm), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = colors, pch=16, asp = 1)

cl1 <- kmeans(rd1, centers = 3)$cluster
#cl1 <- as.numeric(bio)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

lin1 <- getLineages(rd1, cl1)
lin1
plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(lin1, lwd = 3)

crv1 <- getCurves(lin1)
crv1
plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(crv1, lwd = 3)

plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
title('iPCS to Motor neuron lineage (human)')
legend("bottomright",  c('iPSC', 'NPC', 'MN'), pch=16, col=c("#E41A1C", "#377EB8","#4DAF4A"))
lines(crv1, lwd = 3)

'''
require(gam)
t <- pseudotime(crv1)
y <- best_norm
gam.pval <- apply(y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})

'''

####################

# ordered by fano

fano_order <- names(apply(goodCounts,1, function(x){var(x)/mean(x)})[order(apply(goodCounts,1, 
                                                                                 function(x){var(x)/mean(x)}), decreasing=TRUE)])

top_fano <- fano_order[1:500]

fano_log_order <- names(apply(log1p(goodCounts),1, function(x){var(x)/mean(x)})[order(apply(log1p(goodCounts),1, 
                                                                                 function(x){var(x)/mean(x)}), decreasing=TRUE)])


top_log_fano <- fano_log_order[1:500]



fano_norm <- names(apply(best_norm,1, function(x){var(x)/mean(x)})[order(apply(best_norm,1, 
                                                                                function(x){var(x)/mean(x)}), decreasing=TRUE)])

top_norm_fano <- fano_norm[1:500]

#################

# Slingshot with top fano

pca <- prcomp(t(as.matrix(best_norm[top_norm_fano,])))
cells <- goodMeta[rownames(pca$x),]$cell_type
colors <- c()
for (i in 1:length(cells)) {
  if (cells[i]=='iPSC') {
    colors <- c(colors,'red')
  } else if (cells[i]=='NPC') {
    colors <- c(colors,'forestgreen')
  } else if (cells[i]=='MN') {
    colors <- c(colors,'blue')
  }
}

rd1 <- pca$x[,1:2]

plot(rd1, col = colors, pch=16, asp = 1)

#cl1 <- kmeans(rd1, centers = 3)$cluster
cl1 <- as.numeric(bio)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)


lin1 <- getLineages(rd1, cl1)
lin1
plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(lin1, lwd = 3)

crv1 <- getCurves(lin1)
crv1
plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(crv1, lwd = 3)

plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
title('iPCS to Motor neuron lineage (human)')
legend("bottomright",  c('iPSC', 'NPC', 'MN'), pch=16, 
       col=c("#E41A1C", "#377EB8","#4DAF4A"))
lines(crv1, lwd = 3)

##

#meta$pseudotime = t
#ordered_meta <- meta[order(meta$pseudotime),]
#ordered <- best_norm[,rownames(ordered_meta)]


###################

rd1 <- pca$x[,1:2]

plot(rd1, col = colors, pch=16, asp = 1)

cl1 <- kmeans(rd1, centers = 3)$cluster
#cl1 <- as.numeric(bio)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

lin1 <- getLineages(rd1, cl1)
lin1
plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(lin1, lwd = 3)

crv1 <- getCurves(lin1)
crv1
plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(crv1, lwd = 3)
'''
require(gam)
t <- pseudotime(crv1)
y <- best_norm
gam.pval <- apply(y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})
'''

##########################

# known states

rd1 <- pca$x[,1:2]

plot(rd1, col = colors, pch=16, asp = 1)

#cl1 <- kmeans(rd1, centers = 3)$cluster
cl1 <- bio
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

lin1 <- getLineages(rd1, cl1, end.clus = 'MN')
lin1
plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(lin1, lwd = 3)

crv1 <- getCurves(lin1)
crv1
plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(crv1, lwd = 3)
'''
require(gam)
t <- pseudotime(crv1)
y <- best_norm
gam.pval <- apply(y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})
'''

############ SRSF

SRSF_psi <- read.table('SRSF.tab', header = TRUE, row.names = 1)
#SRSF_asc <- read.table('ase_counts.tab', header = TRUE, row.names = 1)

SRSF10 <- SRSF_psi['SRSF10',]

#pseudotime <- max(t) - t
'''
plot(pseudotime[SRSF10 != -1], SRSF10[SRSF10 != -1], col=bio[SRSF10 != -1], xlab='Pseudotime', 
     ylab= expression(paste("SRSF10 ", Psi)))
title(expression(paste("Change in SRSF10 ", Psi)))
legend("bottomright",  c('iPSC', 'NPC', 'MN'), pch=1, col=c("black", "green","red"))

cor(pseudotime[SRSF10 != -1], SRSF10[SRSF10 != -1])

SRSF10_expression <- best_norm["ENSG00000188529",]

plot(pseudotime[SRSF10 != -1],SRSF10_expression[SRSF10 != -1], col=bio[SRSF10 != -1])



plot(pseudotime[SRSF10 != -1], SRSF10[SRSF10 != -1]/SRSF10_expression[SRSF10 != -1], col=bio[SRSF10 != -1])

plot(pseudotime[SRSF10 != -1], SRSF10[SRSF10 != -1]/goodCounts["ENSG00000188529",][SRSF10 != -1], col=bio[SRSF10 != -1])

#####################
plot(pseudotime[SRSF10 != -1], SRSF10[SRSF10 != -1], col=bio[SRSF10 != -1])



srsf6 <- "chr20:43458361:43458509:+@chr20:43459153:43459420:+@chr20:43459771:43459895:+"
srsf6 <- SRSF_psi[srsf6,rownames(t)]

plot(pseudotime[srsf6 >= 0], srsf6[srsf6 >= 0], col=bio[srsf6 >= 0], xlab='Pseudotime', 
     ylab= expression(paste("SRSF6 ", Psi)))
title(expression(paste("Change in SRSF6 ", Psi)))
legend("bottomright",  c('iPSC', 'NPC', 'MN'), pch=1, col=c("black", "green","red"))

##################################

srsf3 <- "chr6:36598849:36598983:+@chr6:36599821:36600276:+@chr6:36601152:36601190:+"
srsf3 <- SRSF_psi[srsf3,rownames(t)]

plot(pseudotime[srsf3 >= 0], srsf3[srsf3 >= 0], col=bio[srsf3 >= 0], xlab='Pseudotime', 
     ylab= expression(paste("SRSF3 ", Psi)))
title(expression(paste("Change in SRSF3 ", Psi)))
legend("topleft",  c('iPSC', 'NPC', 'MN'), pch=1, col=c("black", "green","red"))


plot(pseudotime[srsf3 > 0], srsf3[srsf3 > 0], col=bio[srsf3 > 0])


plot(pseudotime[srsf3 >= 0], best_norm["ENSG00000112081",][srsf3 >= 0], col=bio[srsf3 >= 0])



######################
######################

for (srsf in rownames(SRSF_psi)) {
  srsf_dat <- SRSF_psi[srsf,rownames(t)]
  plot(pseudotime[srsf_dat >= 0], srsf_dat[srsf_dat >= 0], col=bio[srsf_dat >= 0])
  title(srsf)
}

srsf7 <- "chr2:38749529:38749705:-@chr2:38748898:38749173:-@chr2:38748579:38748653:-"
srsf7 <- SRSF_psi[srsf7,rownames(t)]
plot(pseudotime[srsf7 >= 0], srsf7[srsf7 >= 0], col=bio[srsf7 >= 0], xlab='Pseudotime', 
     ylab= expression(paste("SRSF7 (SE) ", Psi)))
title(expression(paste("Change in SRSF7 SE ", Psi)))
legend("topleft",  c('iPSC', 'NPC', 'MN'), pch=1, col=c("black", "green","red"))


srsf7_1 <- "chr2:38749529:38749705:-@chr2:38748898:38749346:-@chr2:38748579:38748653:-"

srsf7_1 <- SRSF_psi[srsf7_1,rownames(t)]
plot(pseudotime[srsf7_1 >= 0], srsf7_1[srsf7_1 >= 0], col=bio[srsf7_1 >= 0], xlab='Pseudotime', 
     ylab= expression(paste("SRSF7 (SE) ", Psi)))
title(expression(paste("Change in SRSF7 SE2 ", Psi)))
legend("topleft",  c('iPSC', 'NPC', 'MN'), pch=1, col=c("black", "green","red"))

srsf7_2 <- "chr2:38749705-38749529:-@chr2:38748653-38748579:-" 
srsf7_2 <- SRSF_psi[srsf7_2,rownames(t)]
plot(pseudotime[srsf7_2 >= 0], srsf7_2[srsf7_2 >= 0], col=bio[srsf7_2 >= 0], xlab='Pseudotime', 
     ylab= expression(paste("SRSF7 (RI) ", Psi)))
title(expression(paste("Change in SRSF7 RI ", Psi)))
legend("topleft",  c('iPSC', 'NPC', 'MN'), pch=1, col=c("black", "green","red"))

srsf1 <- "chr17:58005973-58004923:-@chr17:58004001-58000919:-"

srsf1 <- SRSF_psi[srsf1,rownames(t)]
plot(pseudotime[srsf1 >= 0], srsf1[srsf1 >= 0], col=bio[srsf1 >= 0], xlab='Pseudotime', 
     ylab= expression(paste("SRSF1 (RI) ", Psi)))
title(expression(paste("Change in SRSF1 RI ", Psi)))
legend("bottomleft",  c('iPSC', 'NPC', 'MN'), pch=1, col=c("black", "green","red"))


'''


################# BULK DATA ###################

library(clusterExperiment)

srsf1 <- "chr17:58005973-58004923:-@chr17:58004001-58000919:-"
srsf3 <- "chr6:36598849:36598983:+@chr6:36599821:36600276:+@chr6:36601152:36601190:+"
srsf6 <- "chr20:43458361:43458509:+@chr20:43459153:43459420:+@chr20:43459771:43459895:+"
srsf7 <- "chr2:38749529:38749705:-@chr2:38748898:38749173:-@chr2:38748579:38748653:-"
srsf7_1 <- "chr2:38749529:38749705:-@chr2:38748898:38749346:-@chr2:38748579:38748653:-"
srsf7_2 <- "chr2:38749705-38749529:-@chr2:38748653-38748579:-" 

pooled_samples <- c('SRR4047331', 'SRR4047332', 'SRR4047260', 'SRR4047278', 'SRR4047390', 'SRR4047295', 'SRR4047315', 'SRR4047323')

pooled_psi <- SRSF_psi[c(srsf1, srsf3, srsf6, srsf7, srsf7_1, srsf7_2, 'SRSF10'),pooled_samples]

colnames(pooled_psi) <- c('iPSC_1', 'iPSC_2', 'NPC_1', 'NPC_2', 'NPC_3', 'MN_1', 'MN_2', 'MN_3')
rownames(pooled_psi) <- c('SRSF1', 'SRSF3', 'SRSF6', 'SRSF7_1', 'SRSF7_2', 'SRSF7_3', 'SRSF10')
#write.table(pooled_psi, file='pooled.tsv', quote=FALSE, sep='\t')

par(mar=c(1,1,1,1))

heatmap.2(as.matrix(pooled_psi), trace='none', dendrogram='none', 
          Colv=FALSE, Rowv = FALSE)


bk2 = unique(c(-1, seq(-0.1, 1, length=11)))

#set different color vectors for each interval
col1 ="gray40" #set the order of greys
col2 = colorRampPalette(c("blue", 'white', "red"))(10)
colors2 <- c(col1, col2)

normalized <- pooled_psi
normalized['SRSF6', c('iPSC_1', 'MN_1', 'MN_2', 'MN_3')] <- 0.03
normalized <- t(apply(normalized, 1, function(x)(x-min(x))/(max(x)-min(x))))
normalized['SRSF6', c('iPSC_1', 'MN_1', 'MN_2', 'MN_3')] <- -1



#draw heatmap
hm2 <- pheatmap(as.matrix(normalized), color=colors2, breaks=bk2, scale="none", cluster_rows=F, cluster_cols=F)


#############################
'''
proba <- c()
for (i in 1:5) {
  n <- i*20
  proba <- c(proba, mean(as.numeric(SRSF10[,(pseudotime[,1] <= n) & (pseudotime[,1] > n-10) & (SRSF10[1,] != -1)])))
}

plot(pseudotime[SRSF10 != -1], SRSF10[SRSF10 != -1], col=bio[SRSF10 != -1], xlab='Pseudotime', 
     ylab= expression(paste("SRSF10 ", Psi)))
#lines((0:4)*20+10, proba, type='l')
title(expression(paste("Change in SRSF10 ", Psi)))
legend("bottomright",  c('iPSC', 'NPC', 'MN'), pch=1, col=c("black", "green","red"))


############################################

bonferroni <- p.adjust(gam.pval, method='bonferroni')

#for visuals, only top 1000
topgenes1 <- names(sort(gam.pval[bonferroni <= 0.05], decreasing = FALSE))[1:1000]
heatdata1 <- best_norm[rownames(best_norm) %in% topgenes1, order(pseudotime, na.last = NA)]
heatclus1 <- bio[order(pseudotime, na.last = NA)]


ce1 <- clusterExperiment(heatdata1, heatclus1, transformation=identity)
plotHeatmap(ce1, clusterSamplesData="orderSamplesValue")

SF_TF <- read.table("HanList.txt", sep="\t", row.names = 1, col.names = c('ensembl_id', 'annotation'))

SF <- rownames(SF_TF)[SF_TF$annotation %in% c('Splicing factor/RBP', 'Overlap')]

convertMouseENSid <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = x , 
                   mart = mouse, attributesL = c("ensembl_gene_id"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}


hSF <- convertMouseENSid(SF)

sigSF <- hSF[hSF %in% topgenes1]

heatdata1_SF <- best_norm[rownames(best_norm) %in% sigSF, order(pseudotime, na.last = NA)]
heatclus1_SF <- bio[order(pseudotime, na.last = NA)]

ce1_SF <- clusterExperiment(heatdata1_SF, heatclus1_SF, transformation=identity)
plotHeatmap(ce1_SF, clusterSamplesData="orderSamplesValue")



#######
topgenes <- names(sort(gam.pval[bonferroni <= 0.05], decreasing = FALSE))
heatdata <- best_norm[rownames(best_norm) %in% topgenes, order(pseudotime, na.last = NA)]
heatclus <- bio[order(pseudotime, na.last = NA)]

ce <- clusterExperiment(heatdata, heatclus, transformation=identity)
plotHeatmap(ce, clusterSamplesData="orderSamplesValue")
'''
############################################


sigSF <- hSF[hSF %in% topgenes]

heatdata_SF <- best_norm[rownames(best_norm) %in% sigSF, order(pseudotime, na.last = NA)]
heatclus_SF <- bio[order(pseudotime, na.last = NA)]

ce_SF <- clusterExperiment(heatdata_SF, heatclus_SF, transformation=identity)
plotHeatmap(ce_SF, clusterSamplesData="orderSamplesValue")


pca_SF <- prcomp(t(best_norm[sigSF,]))
plot(pca_SF$x[SRSF10 != -1,'PC1'], pca_SF$x[SRSF10 != -1,'PC2'], col=colores, asp = 1, pch = 16)

colfunc<-colorRampPalette(c("blue","grey","red"))
colores <- colfunc(100)[as.numeric(cut(SRSF10[SRSF10 != -1], breaks = 100))]


plot(pca_SF[SRSF10 != -1,], col = colores, asp = 1, pch = 16)

##############

sigSRSF <- hg38$ensembl_gene_id[hg38$hgnc_symbol %in% SRSF]
sigSRSF <- sigSRSF[sigSRSF %in% sigSF]

heatdata_SRSF <- best_norm[rownames(best_norm) %in% sigSRSF, order(pseudotime, na.last = NA)]
heatclus_SRSF <- bio[order(pseudotime, na.last = NA)]

ce_SRSF <- clusterExperiment(heatdata_SRSF, heatclus_SRSF, transformation=identity)
plotHeatmap(ce_SRSF, clusterSamplesData="orderSamplesValue")


############################################
# Correlations


correlations <- sapply(sigSF, function(x) cor(best_norm[x, SRSF10 != -1], 
                                              SRSF10[SRSF10 != -1]))

corSF <- sigSF[abs(correlations) >= 0.2]



pca_SF <- prcomp(t(best_norm[corSF,]))
plot(pca_SF$x[SRSF10 != -1,'PC1'], pca_SF$x[SRSF10 != -1,'PC2'], col=colores, asp = 1, pch = 16)

colfunc<-colorRampPalette(c("blue","grey","red"))
colores <- colfunc(100)[as.numeric(cut(SRSF10[SRSF10 != -1], breaks = 100))]


plot(pca_SF$x[SRSF10 != -1,], col = colores, asp = 1, pch = 16)


heatdata_cor <- best_norm[rownames(best_norm) %in% corSF[corSF %in% topgenes1], order(pseudotime, na.last = NA)]
heatclus_cor <- bio[order(pseudotime, na.last = NA)]

ce_cor <- clusterExperiment(heatdata_cor, heatclus_cor, transformation=identity)
plotHeatmap(ce_cor, clusterSamplesData="orderSamplesValue")


##################################

SRSF10_bin <- as.numeric(SRSF10[SRSF10 != -1] >= 0.4)


heatdata_cor <- best_norm[rownames(best_norm[,SRSF10 != -1]) %in% corSF[corSF %in% topgenes1], 
                          order(pseudotime[SRSF10 != -1,], na.last = NA)]
heatclus_cor <- SRSF10_bin[order(pseudotime[SRSF10 != -1,], na.last = NA)]

ce_cor <- clusterExperiment(heatdata_cor, heatclus_cor, transformation=identity)
plotHeatmap(ce_cor, clusterSamplesData="orderSamplesValue")



pca_SRSF10 <- prcomp(t(best_norm[corSF,SRSF10 != -1]))
plot(pca_SRSF10$x[,'PC1'], pca_SRSF10$x[,'PC2'], col=SRSF10_bin+1, asp = 1, pch = 16)

#####################

bin_logit <- glm(SRSF10_bin ~ t(best_norm[corSF,SRSF10 != -1]), family = "binomial")
summary(bin_logit)

logit <- glm(SRSF10[SRSF10 != -1] ~ t(best_norm[corSF,SRSF10 != -1]), family = "binomial")
summary(logit)

#correlations <- sapply(sigSF, function(x) cor.test(best_norm[x, SRSF10 != -1], 
#                                              SRSF10[SRSF10 != -1], method = 'pearson'))


SRSF_id <- hg38$ensembl_gene_id[hg38$hgnc_symbol %in% SRSF]
heatdata_SRSF <- best_norm[rownames(best_norm) %in% SRSF_id, order(pseudotime, na.last = NA)]
heatclus_SRSF <- bio[order(pseudotime, na.last = NA)]

ce_SRSF <- clusterExperiment(heatdata_SRSF, heatclus_SRSF, transformation=identity)
plotHeatmap(ce_SRSF, clusterSamplesData="orderSamplesValue")

ens2sym <- function(x) {
  symbol <- hg38$hgnc_symbol[hg38$ensembl_gene_id %in% x]
  return(symbol)
}

sym2ens <- function(x) {
  eid <- hg38$ensembl_gene_id[hg38$hgnc_symbol %in% x]
  return(eid)
}


############################################
# Checking correlations of SRSF genes
# All SF correlated to SRSF10, not only sig changes


expSF <- hSF[hSF %in% rownames(best_norm)]

cor_all <- sapply(expSF, function(x) cor(best_norm[x, SRSF10 != -1], 
                                              SRSF10[SRSF10 != -1]))

corSF_all <- expSF[abs(cor_all) >= 0.3]





pca_SF_all <- prcomp(t(best_norm[corSF_all,]))
plot(pca_SF_all$x[SRSF10 != -1,'PC1'], pca_SF_all$x[SRSF10 != -1,'PC2'], col=colores, asp = 1, pch = 16)

colfunc<-colorRampPalette(c("blue","grey","red"))
colores <- colfunc(100)[as.numeric(cut(SRSF10[SRSF10 != -1], breaks = 100))]

plot(pca_SF_all$x[SRSF10 != -1,], col = colores, asp = 1, pch = 16)

heatdata_cor <- best_norm[rownames(best_norm) %in% corSF[corSF %in% topgenes1], order(pseudotime, na.last = NA)]
heatclus_cor <- bio[order(pseudotime, na.last = NA)]

ce_cor <- clusterExperiment(heatdata_cor, heatclus_cor, transformation=identity)
plotHeatmap(ce_cor, clusterSamplesData="orderSamplesValue")



logit <- glm(SRSF10[SRSF10 != -1] ~ t(best_norm[corSF,SRSF10 != -1]), family = "binomial")
summary(logit)
###############################################
# Inside individual types of cells
# This is the plot I sent to Liana

SRSF10 <- SRSF10[,rownames(meta)]

npc <- best_norm[,bio == 'NPC' & SRSF10 != -1]
mn <- best_norm[,bio == 'MN' & SRSF10 != -1]

expSF <- hSF[hSF %in% rownames(mn)]

cor_all <- sapply(expSF, function(x) cor(mn[x,], 
                                         SRSF10[bio == 'MN' & SRSF10 != -1]))

corSF_all <- expSF[abs(cor_all) >= 0.25 & ! is.na(cor_all)]


colfunc<-colorRampPalette(c("blue","grey","red"))
colores <- colfunc(100)[as.numeric(cut(SRSF10[bio == 'MN' & SRSF10 != -1], breaks = 100))]



pca_SF_all <- prcomp(t(mn[corSF_all,]))
plot(pca_SF_all$x[,'PC1'], pca_SF_all$x[,'PC2'], col=colores, asp = 1, pch = 16, cex=2)

col.labels <- round(seq(0,1,length=3) ,digits=1)
color.legend( xl =20 , yb = -7, xr = 23, yt = 0 , # the coordinates
              legend = col.labels , gradient="y", 
              rect.col=colfunc(100), align="rwb")



ordered <- names(SRSF10[,bio == 'MN' & SRSF10 != -1][order(SRSF10[,bio == 'MN' & SRSF10 != -1])])



bind_genes <- c('AARS',
                'AGO2',
                'AGO3',
                'ATXN2',
                'CAPRIN1',
                'CPSF1',
                'CPSF2',
                'CPSF3',
                'CPSF4',
                'CPSF6',
                'CPSF7',
                'CSTF2',
                'CSTF2T',
                'DDX24',
                'DDX55',
                'DGCR8',
                'EIF3H',
                'EIF4A3',
                'ELAVL1',
                'FBL',
                'FIP1L1',
                'FMR1',
                'FUBP3',
                'FUS',
                'FXR1',
                'FXR2',
                'GPKOW',
                'GRWD1',
                'HLTF',
                'HNRNPA1',
                'HNRNPA2B1',
                'HNRNPC',
                'HNRNPD',
                'HNRNPF',
                'HNRNPH',
                'HNRNPM',
                'HNRNPU',
                'HNRNPUL1',
                'IGF2BP1',
                'IGF2BP2',
                'IGF2BP3',
                'ILF3',
                'KHDRBS1',
                'KHSRP',
                'LARP4',
                'LIN28A',
                'LIN28B',
                'MOV10',
                'NONO',
                'NOP56',
                'NOP58',
                'NUDT21',
                'PPIG',
                'PRPF8',
                'PTBP1',
                'PTBP1PTBP2',
                'PUM2',
                'QKI',
                'RBFOX2',
                'RBM22',
                'RBM47',
                'RPS3',
                'RPS5',
                'RTCB',
                'SAFB2',
                'SF3A3',
                'SF3B4',
                'SFPQ',
                'SMNDC1',
                'SND1',
                'SRRM4',
                'SRSF1',
                'SRSF10',
                'SRSF3',
                'SRSF7',
                'SRSF9',
                'SUB1',
                'SUGP2',
                'TAF15',
                'TARDBP',
                'TBRG4',
                'TIA1',
                'TRA2A',
                'U2AF1',
                'U2AF2',
                'UCHL5',
                'UPF1',
                'YBX3',
                'YTHDC1',
                'YTHDF2',
                'YWHAG',
                'ZC3H7B',
                'ZNF622')

bind_ens <- hg38$ensembl_gene_id[hg38$hgnc_symbol %in% bind_genes]

pca_SF_all <- prcomp(t(mn[rownames(mn)[(rownames(mn) %in% bind_ens)],]))
plot(pca_SF_all$x[,'PC1'], pca_SF_all$x[,'PC2'], col=colores, asp = 1, pch = 16, cex=2)


###############


heatdata_cor <- mn[corSF_all, ordered]
heatclus_cor <- as.numeric(SRSF10[,ordered] >= 0.5) + 1

ce_cor <- clusterExperiment(heatdata_cor, heatclus_cor, transformation=identity)
plotHeatmap(ce_cor, clusterSamplesData="orderSamplesValue", colorScale='rainbow')


heatdata_cor <- t(scale(t(mn[corSF_all, ordered])))
heatclus_cor <- as.numeric(SRSF10[,ordered] >= 0.5) + 1

ce_cor <- clusterExperiment(heatdata_cor, heatclus_cor, transformation=identity)
plotHeatmap(ce_cor, clusterSamplesData="orderSamplesValue")


heatmap.2(heatdata_cor, col=bluered(149), scale="row",key=TRUE, symkey=TRUE, Colv = FALSE,
          density.info="none", trace="none", cexRow=0.5)


my.breaks <- c(seq(-4.4, -1.6, length.out=30),seq(-1.5, 1.5, length.out=90),seq(1.6,4.4, length.out=30))

heatmap.2(heatdata_cor, col=bluered(149),key=TRUE, symkey=TRUE, Colv = FALSE,
          density.info="none", trace="none", cexRow=0.5, dendrogram = 'row',breaks = my.breaks)



my.breaks <- c(seq(-4.4, -0.6, length.out=30),seq(-0.5, 1, length.out=90),seq(1.1,4.4, length.out=30))

colors <- colorpanel(149,"midnightblue","mediumseagreen","yellow") 
heatmap.2(heatdata_cor, col=colors,key=TRUE, symkey=FALSE, density.info="none", 
          trace="none", Colv=FALSE, dendrogram = 'row', cexRow=0.5, breaks=my.breaks)


#########################
# all SF


pca_SF_all <- prcomp(t(mn[expSF,]))
plot(pca_SF_all$x[,'PC1'], pca_SF_all$x[,'PC2'], col=colores, asp = 1, pch = 16)

ordered <- names(SRSF10[,bio == 'MN' & SRSF10 != -1][order(SRSF10[,bio == 'MN' & SRSF10 != -1])])


heatdata_cor <- mn[expSF, ordered]
heatclus_cor <- as.numeric(SRSF10[,ordered] >= 0.5) + 1

ce_cor <- clusterExperiment(heatdata_cor, heatclus_cor, transformation=identity)
plotHeatmap(ce_cor, clusterSamplesData="orderSamplesValue", colorScale='rainbow')


heatdata_cor <- t(scale(t(mn[expSF, ordered])))
heatclus_cor <- as.numeric(SRSF10[,ordered] >= 0.5)

ce_cor <- clusterExperiment(heatdata_cor, heatclus_cor, transformation=identity)
plotHeatmap(ce_cor, clusterSamplesData="orderSamplesValue", colorScale='rainbow')


#########################
### colores:
#
#'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral', 'Accent', 
#'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3', 'Blues', 'BuGn', 'BuPu', 
#'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 
#'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd'.



####

### Slingshot again
### K means
### Separate by SRSF10
### Same genes as in NPC?
### In imputed data


npc <- best_norm[,bio == 'NPC' & SRSF10 != -1]

expSFnpc <- hSF[hSF %in% rownames(npc)]

cor_allnpc <- sapply(expSFnpc, function(x) cor(npc[x,], 
                                         SRSF10[bio == 'NPC' & SRSF10 != -1]))

corSF_allnpc <- expSF[abs(cor_allnpc) >= 0.25 & ! is.na(cor_allnpc)]


colfunc<-colorRampPalette(c("blue","grey","red"))
colores <- colfunc(100)[as.numeric(cut(SRSF10[bio == 'NPC' & SRSF10 != -1], breaks = 100))]



pca_SF_all_npc <- prcomp(t(npc[corSF_all,]))
plot(pca_SF_all_npc$x[,'PC1'], pca_SF_all_npc$x[,'PC2'], col=colores, asp = 1, pch = 16)




###

heatdata_cor <- best_norm[rownames(best_norm) %in% corSF_all[corSF_all %in% topgenes1], 
                          order(pseudotime, na.last = NA)]
heatclus_cor <- bio[order(pseudotime, na.last = NA)]

ce_cor <- clusterExperiment(heatdata_cor, heatclus_cor, transformation=identity)
plotHeatmap(ce_cor, clusterSamplesData="orderSamplesValue")



############################################

# Un-comment to write tables
#write.table(rd1, file='rsem_PCA.tsv', quote=FALSE, 
#            sep='\t', col.names = NA)

#write.table(pseudotime, file='rsem_pseudotime.tsv', quote=FALSE, 
#            sep='\t', col.names = NA)

#write.table(best_norm, file='rsem_normalized_scone.tsv', quote=FALSE, 
#            sep='\t', col.names = NA)





pca_SRSF <- prcomp(t(best_norm[sym2ens(SRSF),]))
plot(pca_SRSF$x[SRSF10 != -1,'PC1'], pca_SRSF$x[SRSF10 != -1,'PC2'], col=colores, asp = 1, pch = 16)

colfunc<-colorRampPalette(c("blue","grey","red"))
colores <- colfunc(100)[as.numeric(cut(SRSF10[SRSF10 != -1], breaks = 100))]

plot(pca_SRSF$x[SRSF10 != -1,], col = colores, asp = 1, pch = 16)

heatdata_cor <- best_norm[rownames(best_norm) %in% sym2ens(SRSF), order(pseudotime, na.last = NA)]
heatclus_cor <- bio[order(pseudotime, na.last = NA)]

ce_cor <- clusterExperiment(heatdata_cor, heatclus_cor, transformation=identity)
plotHeatmap(ce_cor, clusterSamplesData="orderSamplesValue")



logit <- glm(SRSF10[SRSF10 != -1] ~ t(best_norm[sym2ens(SRSF),SRSF10 != -1]), family = "binomial")
summary(logit)



