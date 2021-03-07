
## OBTAINING NECESSARY LIBRARIES

library('SymSim')


## SETTING RANDOM SEED FOR CONSISTENCY

seed <- 1
set.seed(seed)


## CREATING A 2-BRANCH PHYLOGENETIC TREE AND SETTING UP PARAMETERS

#phylo_2tips <- read.tree("../data/newick_AB.txt")
ngenes = 5000
ncells_total = 1000
#plot(phylo_2tips)


Phyla2 <- function(plotting=F){
  # par(mfrow=c(2,2))
  phyla<-rtree(2)
  phyla <- compute.brlen(phyla,1)
  edges <- cbind(phyla$edge,phyla$edge.length)
  edges <- cbind(c(1:length(edges[,1])),edges)
  connections <- table(c(edges[,2],edges[,3]))
  root <- as.numeric(names(connections)[connections==2])
  tips <- as.numeric(names(connections)[connections==1])
  phyla$tip.label <- as.character(tips)
  
  if(plotting==T){
    plot(phyla,show.tip.label = F,lwd=2)
    tiplabels(cex=2)
    nodelabels(cex=2)
  }
  return(phyla)
}

phylo_2tips <- Phyla2()


## GETTING ISOFORM LENGTHS

# get gene length pool, for simulation of sequencing
# gene_len_pool: array of exon lengths
data(gene_len_pool)

# selecting lengths between 0.1-0.9 quantiles
gene_len_pool_selected <- gene_len_pool[((gene_len_pool >= quantile(gene_len_pool, probs=0.1)) & (gene_len_pool <= quantile(gene_len_pool, probs=0.9)))]

# get exon length list from a pool of exon lengths from the human genome (the distribution of exon lengths is similar in mouse).
# exon_lengths: array of exon lengths
exon_lengths <- scan('../data/exon_lengths.txt')

# This is needed for the fragmentation simulation.
load("../data/len2nfrag.RData")

# Here we are creating two arrays of isoform lenghts: one for the included, and one for the excluded
# We do this by selecting randomly and independently a gene length X (from gene_len_pool_selected), 
# and an exon length Y from (exon_lengths).
# We then set the the length of the excluded isoform to be X, and the included isoform to be X + Y.

# This for loop is a bit cumbersome because we have to make sure that the lengths are included in SymSim's pool.
# SymSim has an internal pool of gene lengths that for some reason are fixed; probably for optimizing the fragmentation simulation. 
# For this reason, we have to make sure that X + Y is also in the pool of gene lengths, so that it doesn't throw an error. 
# There's probably a way to fix this, but I don't think it really matters, since it is only a handful of cases.

# len_excluded: array of selected excluded isoform length
# len_included: array of included isoforms (excluded + exon length)

counter <- 0
len_excluded <- c()
len_included <- c()

while (counter < ngenes) {
    
    excl <- sample(gene_len_pool_selected, 1)
    exon_len <- sample(exon_lengths, 1)
    incl <- excl + exon_len
    
    if (as.character(incl) %in% rownames(len2nfrag)) {
        
        len_excluded <- c(len_excluded, excl)
        len_included <- c(len_included, incl)
        counter <- counter + 1
        
    }
    
}    

# merge included and excluded lengths into a single array. Maybe isoform_len would be a more appropriate name?
gene_len <- c(len_included, len_excluded)


## SIMULATING TRUE COUNTS USING SYMSIM

true_counts_res <- SimulateTrueCounts(ncells_total=ncells_total, ngenes=ngenes, nevf=20, evf_type="continuous",
                                     n_de_evf=10, vary="s", Sigma=0.6, phyla=phylo_2tips, randseed=1, prop_hge = 0.0)


#true_counts_res <- SimulateTrueCounts(ncells_total=ncells_total, ngenes=ngenes, nevf=20, evf_type="continuous",
#                                      n_de_evf=8, vary="s", Sigma=0.4, phyla=phylo_2tips, randseed=1, prop_hge = 0.0)

write.table(true_counts_res$counts, 'sim_output/gene_counts.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')

print('saved true counts')
#print(true_counts_res$cell_meta)


## LINEAR TRANSFORMATION TO CHANGE THE DEPTHS OF CELLS AND OBTAIN A SINGLE LINEAGE

branch = true_counts_res$cell_meta['pop'][[1]]
#print(true_counts_res$cell_meta['depth'][[1]])
#print(branch)
#print(length(branch))
#print(c(1:length(branch)))


lineage1 <- rownames(true_counts_res$cell_meta)[true_counts_res$cell_meta$pop == '3_1']
lineage2 <- rownames(true_counts_res$cell_meta)[true_counts_res$cell_meta$pop == '3_2']

depth1 <- true_counts_res$cell_meta[lineage1, 'depth']*0.5 + 0.5
depth2 <- abs(1-true_counts_res$cell_meta[lineage2, 'depth'])*0.5

depth <- c(depth1, depth2)*100




# for (i in 1:(length(branch)-1) )
#   if(branch[i] != branch[i+1]){
#     bridge = i+1
#   }
# depth = true_counts_res$cell_meta['depth'][[1]]
# for (i in bridge:length(depth)){
#   depth[i] = depth[i]+50
# }
#print(depth)


## CHOOSING PARAMETERS OF PLATONIC EXONIC IMPULSES

Impulse_Parameters = matrix(nrow = ngenes, ncol=7)
for (i in c( 1:ngenes ) )  {
expit <- function(z) exp(z)/(1+exp(z))
logit <- function(z) log(z/(1-z))
hx = runif(1, min=-3, max=logit(0.7))
hy = runif(1, min=logit(expit(hx)+0.2), max=3) 
hz = runif(1, min=hx, max=hy)
h_vector = c(hx, hy, hz)
h_vector = sample(h_vector)
if (i <= 2500){
  hx = h_vector[1]
  hy = h_vector[2]
  hz = h_vector[3]
} else {
  hx = h_vector[1]
  hy = h_vector[1]
  hz = h_vector[1]
}
t1 = runif(1, min=1, max=70)
t2 = runif(1, min=t1+20, max=99)
beta = runif(1, min=0.05, max=0.25)
sigma = runif(1, min=0.5, max=1)
Impulse_Parameters[i, 1:7] = c(hx, hy, hz, t1, t2, beta, sigma)
}



## IMPULSE EVALUATION FUNCTION

impulse_evaluation <- function (imp_parameters, t) {
hx = imp_parameters[1]
hy = imp_parameters[2]
hz = imp_parameters[3]
t1 = imp_parameters[4]
t2 = imp_parameters[5]
beta = imp_parameters[6]
impulse_value = (1/hy)*( hx + (hy-hx)*(1/(1+exp(-beta*(t - t1)))) )*(hz + (hy-hz)*(1/(1+exp(beta*(t - t2)))))
return(impulse_value)
}


## ADDING BIOLOGICAL NOISE TO SPLICING

Splicing_Z_Actual = matrix(nrow = ngenes, ncol = ncells_total)
Splicing_Z_Sampled = matrix(nrow = ngenes, ncol = ncells_total)
for (i in c( 1:ngenes ) )  {
  parameters = Impulse_Parameters[i,]
  for (j in c(1:ncells_total) ) {
    pos = depth[j]
    splicing_mean = impulse_evaluation(parameters, pos)
    Splicing_Z_Actual[i, j] = splicing_mean
    Splicing_Z_Sampled[i, j] = rnorm(1, mean = splicing_mean, sd = Impulse_Parameters[i,7])
  }   
}
#print(Splicing_Z_Actual)
#print(Splicing_Z_Sampled)


## SIMULATIING mRNA SPLICING

Spliced_In = matrix(nrow = ngenes, ncol = ncells_total)
Spliced_Out = matrix(nrow = ngenes, ncol = ncells_total)
for (i in c(1:ngenes)){
  for (j in c(1:ncells_total)){
    Spliced_In[i,j] =  rbinom(1, true_counts_res$counts[i,j], expit(Splicing_Z_Sampled[i,j]))
    Spliced_Out[i,j] = true_counts_res$counts[i,j] - Spliced_In[i,j]
  }
}
mrna_matrix <- rbind(Spliced_In, Spliced_Out)
Splicing_Z_True <- Spliced_In/(true_counts_res$counts)

#print(dim(mrna_matrix))


## SIMULATING SEQUENCING

#observed_reads <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=true_counts_res[[3]], protocol="nonUMI", alpha_mean=0.1, 
#                  alpha_sd=0.01, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)

true_counts_res$cell_meta$lineage_depth = depth
## SAVING OUTPUTS

write.table(Splicing_Z_Actual, 'sim_output/psi_platonic.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(Splicing_Z_Sampled, 'sim_output/psi_underlying.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(Splicing_Z_True, 'sim_output/psi_true.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(true_counts_res$cell_meta, 'sim_output/meta.tab', quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
write.table(mrna_matrix, 'sim_output/isoform_counts.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')

lapply(gene_len, write, "sim_output/isoform_length.txt", append=TRUE, ncolumns=1000)



observed_reads <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=true_counts_res[[3]], protocol="nonUMI", alpha_mean=0.1, 
                  alpha_sd=0.05, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)

write.table(observed_reads[[1]], 'sim_output/observed_counts_0.1.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')

write.table(observed_reads$cell_meta, 'sim_output/observed_meta_0.1.tab', quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')

##

observed_reads <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=true_counts_res[[3]], protocol="nonUMI", alpha_mean=0.05, 
                  alpha_sd=0.02, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)

write.table(observed_reads[[1]], 'sim_output/observed_counts_0.05.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads$cell_meta, 'sim_output/observed_meta_0.05.tab', quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
##

observed_reads <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=true_counts_res[[3]], protocol="nonUMI", alpha_mean=0.01, 
                  alpha_sd=0.01, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)

write.table(observed_reads[[1]], 'sim_output/observed_counts_0.01.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads$cell_meta, 'sim_output/observed_meta_0.01.tab', quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
