library(zinbwave)
library(matrixStats)
library(magrittr)
library(ggplot2)
library(biomaRt)


counts <- as.matrix(read.table('/mnt/lareaulab/cfbuenabadn/RNASeq/Mouse/Tiklova/star_counts.tab.gz', sep = '\t', 
                         header = TRUE, row.names = 1))

meta <- read.table('/mnt/lareaulab/cfbuenabadn/psix_project/analysis_psix/Mouse_Midbrain_Deveopment/data/SraRunTable.txt.gz', sep = ',', header = TRUE, row.names = 1)

neurogenesis <- SingleCellExperiment(assays=list(counts=counts),colData=meta)


filter <- rowSums(assay(neurogenesis)>5)>5
table(filter)

neurogenesis <- neurogenesis[filter,]

assay(neurogenesis) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(neurogenesis)
vars <- sort(vars, decreasing = TRUE)
head(vars)


assayNames(neurogenesis)[1] <- "counts"

neurogenesis_zinb <- zinbwave(neurogenesis, K = 2, epsilon=1000)

W <- reducedDim(neurogenesis_zinb)

data.frame(W, bio=colData(neurogenesis)$Biological_Condition,
           coverage=colData(neurogenesis)$Coverage_Type) %>%
    ggplot(aes(W1, W2, colour=bio, shape=coverage)) + geom_point() + 
    scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()