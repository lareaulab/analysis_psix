import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os
import sys
from tqdm import tqdm
sys.path.insert(0, '/mnt/lareaulab/cfbuenabadn/psix_project/psix/psix/')
import psix

psix_three_lineage = psix.Psix()
psix_three_lineage.process_rnaseq(
        'three_lineages/processed_tables/SE_counts_0.1.tab.gz',
        'three_lineages/processed_tables/constitutive_introns_0.1.tab.gz',
        'three_lineages/processed_tables/tpm_0.1.tab.gz',
        minJR = 1,
        minCell=1,
        min_observed = 0.25)

psix_three_lineage.compute_psix_scores(latent='three_lineages/processed_tables/pc2_rd.tab.gz', n_jobs=5, 
                                n_neighbors=100, pvals_bins=20, 
                                turbo='/mnt/lareaulab/cfbuenabadn/psix_project/psix/psix/psix_turbo/',
                                           n_random_exons = 1000000)

psix_three_lineage.psix_results.to_csv('psix_three_lineages.tab.gz', sep='\t', index = True, header=True)
