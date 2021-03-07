# ---
# title: "R Notebook"
# output: html_notebook
# ---
## OBTAINING NECESSARY LIBRARIES

library('SymSim')


## SETTING RANDOM SEED FOR CONSISTENCY

seed <- 1
set.seed(seed)


## CREATING A 2-BRANCH PHYLOGENETIC TREE AND SETTING UP PARAMETERS

# phylo_3tips <- Phyla3()
phylo_3tips <- read.tree("tree.txt")
ngenes = 5000
ncells_total = 2000
#plot(phylo_3tips)
#print(phylo_3tips$edge)
#print(phylo_3tips$edge.length)


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

true_counts_res <- SimulateTrueCounts(ncells_total=ncells_total, ngenes=ngenes, nevf=20, evf_type="continuous", n_de_evf=18, vary="s", Sigma=0.8, phyla=phylo_3tips, randseed=1, prop_hge=0.0)


write.table(true_counts_res$counts, 'sim_output/gene_counts.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(true_counts_res$cell_meta, 'sim_output/meta.tab', quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')


## SPLITTING CELLS INTO THREE DIFFERENT LINEAGES

lineages = true_counts_res$cell_meta['pop'][[1]]
counter = 1
lineage_switch = vector()
for (i in 1:(length(lineages)-1) )
  if(lineages[i] != lineages[i+1]){
    lineage_switch[counter] = i+1
    counter = counter+1
  }
print(lineage_switch)
### Lineage 1 is on branch 4_1
lineage1_meta = true_counts_res$cell_meta[c(lineage_switch[3]:ncells_total),]
#print(lineage1_meta)
lineage1_gene_expression = true_counts_res$counts[,c(lineage_switch[3]:ncells_total)]
#print(lineage1_gene_expression)


### Lineage 2 is on branches 4_5, 5_2
#print(lineage2_meta)
lineage2_meta = true_counts_res$cell_meta[c(1:(lineage_switch[2]-1)),]
lineage2_gene_expression = true_counts_res$counts[,c(1: (lineage_switch[2]-1) )]
#print(lineage2_gene_expression)
lineage3_meta = true_counts_res$cell_meta[c(lineage_switch[2]: (lineage_switch[3]-1) ),]
#print(lineage3_meta)
lineage3_gene_expression = true_counts_res$counts[,c(lineage_switch[2]: (lineage_switch[3]-1) )]
#print(lineage3_gene_expression)
#print(dim(true_counts_res$counts))
lineage1_cell_count = dim(lineage1_meta)[1]
lineage2_cell_count = dim(lineage2_meta)[1]
lineage3_cell_count = dim(lineage3_meta)[1]
#print(lineage1_cell_count)
#print(lineage2_cell_count)
#print(lineage3_cell_count)
#print(depth)


## CHOOSING PARAMETERS OF PLATONIC EXONIC SIGMOIDS

Impulse_Parameters_L1 = matrix(nrow = ngenes, ncol=5)
Impulse_Parameters_L2 = matrix(nrow = ngenes, ncol=5)
Impulse_Parameters_L3 = matrix(nrow = ngenes, ncol=5)
for (l in c(1:3)) {
for (i in c( 1:ngenes ) )  {
expit <- function(z) exp(z)/(1+exp(z))
logit <- function(z) log(z/(1-z))
hx = runif(1, min=-4, max=logit(0.7))
hy = runif(1, min=logit(expit(hx)+0.2), max=4) 
#hz = runif(1, min=hx, max=hy)
#h_vector = c(hx, hy, hz)
h_vector = c(hx, hy)
#print(h_vector)
h_vector = sample(h_vector)
#print(h_vector)
hx = h_vector[1]
hy = h_vector[2]
#hz = h_vector[3]
#t1 = runif(1, min=1, max=70)
#t2 = runif(1, min=t1+20, max=99)
####t1 = t1 = runif(1, min=0.1, max=0.9)
#print(t1)
####beta = runif(1, min=0.05, max=0.25)
beta = runif(1, min=5, max=15)
sigma = runif(1, min=0.25, max=1)
if(l==1){
  size_lineage = max(lineage1_meta$depth) - min(lineage1_meta$depth)
  min_l = min(lineage1_meta$depth) + (0.1*size_lineage)
  max_l = min(lineage1_meta$depth) + (0.7*size_lineage)
  t1 = runif(1, min=min_l, max=max_l)
  Impulse_Parameters_L1[i, 1:5] = c(hx, hy, t1, beta, sigma)
}
else if(l==2){
  hx = Impulse_Parameters_L1[i, 1]
  
  h_vector <= c()
  if (hx < logit(0.7)) {h_vector = c(h_vector, runif(1, min=logit(expit(hx)+0.2), max=4))}
  if (hx > logit(0.3)) {h_vector = c(h_vector, runif(1, min=-4, max=logit(expit(hx)-0.2)))}
  
  h_vector = sample(h_vector)
  hy = h_vector[1]

  size_lineage = max(lineage2_meta$depth) - min(lineage2_meta$depth)
  min_l = min(lineage2_meta$depth) + (0.1*size_lineage)
  max_l = min(lineage2_meta$depth) + (0.7*size_lineage)
  t1 = runif(1, min=min_l, max=max_l)
  Impulse_Parameters_L2[i, 1:5] = c(hx, hy, t1, beta, sigma)
}
else{
  size_lineage = max(lineage3_meta$depth) - min(lineage3_meta$depth)
  min_l = min(lineage3_meta$depth) + (0.1*size_lineage)
  max_l = min(lineage3_meta$depth) + (0.7*size_lineage)
  t1 = runif(1, min=min_l, max=max_l)
  Impulse_Parameters_L3[i, 1:5] = c(hx, hy, t1, beta, sigma)
}
}
}


## IMPULSE EVALUATION FUNCTION

impulse_evaluation <- function (offset, imp_parameters, t) {
if (offset== -1000){
    hx = imp_parameters[1]}
else {hx = offset}
    
hy = imp_parameters[2]
#hz = imp_parameters[3]
t1 = imp_parameters[3]
#t2 = imp_parameters[5]
beta = imp_parameters[4]
impulse_value = 0 + ( hx + (hy-hx)*(1/(1+exp(-beta*(t - t1)))) )
return(impulse_value)
}


## ADDING BIOLOGICAL NOISE TO SPLICING


#print(lineage1_meta['depth'][[1]])
Splicing_L1_Z_Actual = matrix(nrow = ngenes, ncol = lineage1_cell_count)
Splicing_L1_Z_Sampled = matrix(nrow = ngenes, ncol = lineage1_cell_count)

Splicing_L2_Z_Actual = matrix(nrow = ngenes, ncol = lineage2_cell_count)
Splicing_L2_Z_Sampled = matrix(nrow = ngenes, ncol = lineage2_cell_count)

Splicing_L3_Z_Actual = matrix(nrow = ngenes, ncol = lineage3_cell_count)
Splicing_L3_Z_Sampled = matrix(nrow = ngenes, ncol = lineage3_cell_count)

l1_diff = rbinom(ngenes, 1, 0.25)
l2_diff = rbinom(ngenes, 1, 0.25)
l3_diff = rbinom(ngenes, 1, 0.25)

for (i in c( 1:ngenes ) )  {
  for (j in c(1:lineage1_cell_count) ) {
    pos = lineage1_meta['depth'][[1]][j]
    splicing_mean = impulse_evaluation(-1000, Impulse_Parameters_L1[i,], pos)
    
    if (j == 1) {x = splicing_mean}
    
    if (l1_diff[i] == 1){
      Splicing_L1_Z_Actual[i, j] = splicing_mean
      Splicing_L1_Z_Sampled[i, j] = rnorm(1, mean = splicing_mean, sd = Impulse_Parameters_L1[i,5])
    }
    else{
      Splicing_L1_Z_Actual[i, j] =x
      Splicing_L1_Z_Sampled[i, j] = rnorm(1, mean = x, sd = Impulse_Parameters_L1[i,5])
    }
  }   


  for (j in c(1:lineage2_cell_count) ) {
    pos = lineage2_meta['depth'][[1]][j]
    splicing_offset = impulse_evaluation(-1000, Impulse_Parameters_L1[i,], 1)
    splicing_mean = impulse_evaluation(splicing_offset, Impulse_Parameters_L2[i,], pos)
    
    if (j == 1) {
      y = splicing_mean
      diff_x_y = y - x
    }
    
    if (l2_diff[i] == 1){
      Splicing_L2_Z_Actual[i, j] = splicing_mean - diff_x_y
      Splicing_L2_Z_Sampled[i, j] = rnorm(1, mean = splicing_mean - diff_x_y, sd = Impulse_Parameters_L2[i,5])
      
      if(j == (lineage2_cell_count - lineage3_cell_count)){
      inflection_point = splicing_mean - diff_x_y
    }
    }
    else{
      Splicing_L2_Z_Actual[i, j] =x
      Splicing_L2_Z_Sampled[i, j] = rnorm(1, mean = x, sd = Impulse_Parameters_L2[i,5])
      
      if(j == (lineage2_cell_count - lineage3_cell_count)){
      inflection_point = x
    }
    }
    
    #Splicing_L2_Z_Actual[i, j] = splicing_mean - diff_x_y
    
    
    
    #Splicing_L2_Z_Sampled[i, j] = rnorm(1, mean = splicing_mean - diff_x_y, sd = Impulse_Parameters_L2[i,5])
  }   

  bridge_position = true_counts_res$cell_meta['depth'][[1]][lineage_switch[2]]

  for (j in c(1:lineage3_cell_count) ) {
    pos = lineage3_meta['depth'][[1]][j]
    splicing_offset = impulse_evaluation(-1000, Impulse_Parameters_L2[i,], bridge_position)
    splicing_mean = impulse_evaluation(splicing_offset, Impulse_Parameters_L3[i,], pos)
    
    if (j == 1) {
      z = splicing_mean
      diff_y_z = z - inflection_point
      #diff_y_z = z - splicing_offset - diff_x_y
    }
    
    if (l3_diff[i] == 1){
      Splicing_L3_Z_Actual[i, j] = splicing_mean - diff_y_z
      Splicing_L3_Z_Sampled[i, j] = rnorm(1, mean = splicing_mean - diff_y_z, sd = Impulse_Parameters_L3[i,5])
    }
    else{
      Splicing_L3_Z_Actual[i, j] =inflection_point
      Splicing_L3_Z_Sampled[i, j] = rnorm(1, mean = inflection_point, sd = Impulse_Parameters_L3[i,5])
    }
    
    
    #Splicing_L3_Z_Actual[i, j] = splicing_mean - diff_y_z
    
    #Splicing_L3_Z_Sampled[i, j] = rnorm(1, mean = splicing_mean - diff_y_z, sd = Impulse_Parameters_L3[i,5])
  }   
}

Splicing_Z_Actual = cbind(Splicing_L2_Z_Actual, Splicing_L3_Z_Actual, Splicing_L1_Z_Actual)
Splicing_Z_Sampled = cbind(Splicing_L2_Z_Sampled, Splicing_L3_Z_Sampled, Splicing_L1_Z_Sampled)
#print(dim(Splicing_Z_Actual))
#print(dim(Splicing_Z_Sampled))
# plot(true_counts_res$cell_meta$depth, Splicing_Z_Actual[1, ])


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

## SAVING mRNA

write.table(Splicing_Z_Actual, 'sim_output/psi_platonic.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(Splicing_Z_Sampled, 'sim_output/psi_underlying.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(Splicing_Z_True, 'sim_output/psi_true.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
#write.table(true_counts_res$cell_meta, 'sim_output/meta.tab', quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
write.table(mrna_matrix, 'sim_output/isoform_counts.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')

lapply(l1_diff, write, "sim_output/l1_diff.txt", append=TRUE, ncolumns=1000)
lapply(l2_diff, write, "sim_output/l2_diff.txt", append=TRUE, ncolumns=1000)
lapply(l3_diff, write, "sim_output/l3_diff.txt", append=TRUE, ncolumns=1000)

lapply(gene_len, write, "sim_output/isoform_length.txt", append=TRUE, ncolumns=1000)

## SIMULATING SEQUENCING & SAVING OUTPUTS

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
