#import numpy as np
import pandas as pd
import psix
from numpy import float16
import os

# psix_single_lineage = psix.Psix()
# psix_single_lineage.junctions2psi(
#        'processed_tables/splice_junction_counts_0.1.tab.gz',
#        '',
#        'processed_tables/tpm_0.1.tab.gz',
#        save_files_in='psix_output_10000/', dtype=float16)

# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=1000)

# psix_results = psix_single_lineage.psix_results
# psix_single_lineage.save_psix_object(psix_dir = 'psix_output_10000/psix_object/')



# os.mkdir('psix_output_10000/k_number')


# #################################

# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=10)
# psix_single_lineage.psix_results.to_csv('psix_output_10000/k_number/k_10.tab.gz', sep='\t', index=True, header=True)

# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=25)
# psix_single_lineage.psix_results.to_csv('psix_output_10000/k_number/k_25.tab.gz', sep='\t', index=True, header=True)


# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=50)
# psix_single_lineage.psix_results.to_csv('psix_output_10000/k_number/k_50.tab.gz', sep='\t', index=True, header=True)

# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=100)
# psix_single_lineage.psix_results.to_csv('psix_output_10000/k_number/k_100.tab.gz', sep='\t', index=True, header=True)


# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=200)
# psix_single_lineage.psix_results.to_csv('psix_output_10000/k_number/k_200.tab.gz', sep='\t', index=True, header=True)

# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=500)
# psix_single_lineage.psix_results.to_csv('psix_output_10000/k_number/k_500.tab.gz', sep='\t', index=True, header=True)

# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=1000)
# psix_single_lineage.psix_results.to_csv('psix_output_10000/k_number/k_1000.tab.gz', sep='\t', index=True, header=True)


# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=2500)
# psix_single_lineage.psix_results.to_csv('psix_output_10000/k_number/k_2500.tab.gz', sep='\t', index=True, header=True)

# del psix_single_lineage



os.mkdir('psix_output_5000/k_number')
os.mkdir('psix_output_1000/k_number')
os.mkdir('psix_output_500/k_number')
os.mkdir('psix_output_100/k_number')


psix_single_lineage = psix.Psix()
psix_single_lineage.junctions2psi(
        'processed_tables/splice_junction_counts_0.1_5000.tab.gz',
        '',
        'processed_tables/tpm_0.1_5000.tab.gz',
        save_files_in='psix_output_5000/', dtype=float16)

#psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=500)
#
#psix_results = psix_single_lineage.psix_results
#psix_single_lineage.save_psix_object(psix_dir = 'psix_output_5000/psix_object/')


psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=10)
psix_single_lineage.psix_results.to_csv('psix_output_5000/k_number/k_10.tab.gz', sep='\t', index=True, header=True)

psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=25)
psix_single_lineage.psix_results.to_csv('psix_output_5000/k_number/k_25.tab.gz', sep='\t', index=True, header=True)

psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=50)
psix_single_lineage.psix_results.to_csv('psix_output_5000/k_number/k_50.tab.gz', sep='\t', index=True, header=True)

psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=100)
psix_single_lineage.psix_results.to_csv('psix_output_5000/k_number/k_100.tab.gz', sep='\t', index=True, header=True)


psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=200)
psix_single_lineage.psix_results.to_csv('psix_output_5000/k_number/k_200.tab.gz', sep='\t', index=True, header=True)


psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=500)
psix_single_lineage.psix_results.to_csv('psix_output_5000/k_number/k_500.tab.gz', sep='\t', index=True, header=True)


psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=1000)
psix_single_lineage.psix_results.to_csv('psix_output_5000/k_number/k_1000.tab.gz', sep='\t', index=True, header=True)


psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=2500)
psix_single_lineage.psix_results.to_csv('psix_output_5000/k_number/k_2500.tab.gz', sep='\t', index=True, header=True)

#psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=2500)
#psix_single_lineage.psix_results.to_csv('psix_output_5000/k_number/k_0.5.tab.gz', sep='\t', index=True, header=True)




del psix_single_lineage







psix_single_lineage = psix.Psix()
psix_single_lineage.junctions2psi(
        'processed_tables/splice_junction_counts_0.1_1000.tab.gz',
        '',
        'processed_tables/tpm_0.1_1000.tab.gz',
        save_files_in='psix_output_1000/', dtype=float16)

#psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_1000.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=100)
#
#psix_results = psix_single_lineage.psix_results
#psix_single_lineage.save_psix_object(psix_dir = 'psix_output_1000/psix_object/')




psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_1000.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=10)
psix_single_lineage.psix_results.to_csv('psix_output_1000/k_number/k_10.tab.gz', sep='\t', index=True, header=True)

psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_1000.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=25)
psix_single_lineage.psix_results.to_csv('psix_output_1000/k_number/k_25.tab.gz', sep='\t', index=True, header=True)

psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_1000.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=50)
psix_single_lineage.psix_results.to_csv('psix_output_1000/k_number/k_50.tab.gz', sep='\t', index=True, header=True)


psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_1000.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=100)
psix_single_lineage.psix_results.to_csv('psix_output_1000/k_number/k_100.tab.gz', sep='\t', index=True, header=True)

psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_1000.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=200)
psix_single_lineage.psix_results.to_csv('psix_output_1000/k_number/k_200.tab.gz', sep='\t', index=True, header=True)

psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_1000.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=500)
psix_single_lineage.psix_results.to_csv('psix_output_1000/k_number/k_500.tab.gz', sep='\t', index=True, header=True)


del psix_single_lineage








psix_single_lineage = psix.Psix()
psix_single_lineage.junctions2psi(
        'processed_tables/splice_junction_counts_0.1_500.tab.gz',
        '',
        'processed_tables/tpm_0.1_500.tab.gz',
        save_files_in='psix_output_500/', dtype=float16)

#psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_500.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=50)
#
#psix_results = psix_single_lineage.psix_results
#psix_single_lineage.save_psix_object(psix_dir = 'psix_output_500/psix_object/')


psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_500.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=10)
psix_single_lineage.psix_results.to_csv('psix_output_500/k_number/k_10.tab.gz', sep='\t', index=True, header=True)

psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_500.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=25)
psix_single_lineage.psix_results.to_csv('psix_output_500/k_number/k_25.tab.gz', sep='\t', index=True, header=True)


psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_500.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=50)
psix_single_lineage.psix_results.to_csv('psix_output_500/k_number/k_50.tab.gz', sep='\t', index=True, header=True)


psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_500.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=100)
psix_single_lineage.psix_results.to_csv('psix_output_500/k_number/k_100.tab.gz', sep='\t', index=True, header=True)


psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_500.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=100)
psix_single_lineage.psix_results.to_csv('psix_output_500/k_number/k_0.2.tab.gz', sep='\t', index=True, header=True)


psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_500.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=200)
psix_single_lineage.psix_results.to_csv('psix_output_500/k_number/k_200.tab.gz', sep='\t', index=True, header=True)

#psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_500.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=250)
#psix_single_lineage.psix_results.to_csv('psix_output_500/k_number/k_0.5.tab.gz', sep='\t', index=True, header=True)




del psix_single_lineage






psix_single_lineage = psix.Psix()
psix_single_lineage.junctions2psi(
        'processed_tables/splice_junction_counts_0.1_100.tab.gz',
        '',
        'processed_tables/tpm_0.1_100.tab.gz',
        save_files_in='psix_output_100/', dtype=float16)

#psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_100.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=10)
#
#psix_results = psix_single_lineage.psix_results
#psix_single_lineage.save_psix_object(psix_dir = 'psix_output_100/psix_object/')





psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_100.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=10)
psix_single_lineage.psix_results.to_csv('psix_output_100/k_number/k_10.tab.gz', sep='\t', index=True, header=True)

psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_100.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=25)
psix_single_lineage.psix_results.to_csv('psix_output_100/k_number/k_25.tab.gz', sep='\t', index=True, header=True)

psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_100.tab.gz', n_jobs=25,
                                n_random_exons=2000, n_neighbors=50)
psix_single_lineage.psix_results.to_csv('psix_output_100/k_number/k_50.tab.gz', sep='\t', index=True, header=True)

