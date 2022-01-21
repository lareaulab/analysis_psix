#import numpy as np
import pandas as pd
import psix
from numpy import float16


psix_single_lineage = psix.Psix()
psix_single_lineage.junctions2psi(
       'processed_tables/splice_junction_counts_0.1.tab.gz',
       '',
       'processed_tables/tpm_0.1.tab.gz',
       save_files_in='psix_output_10000/', dtype=float16)

#psix_single_lineage.run_psix(latent='processed_tables/pc2_rd.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=1000)
#
#psix_results = psix_single_lineage.psix_results
#psix_single_lineage.save_psix_object(psix_dir = 'psix_output_10000/psix_object/')



# os.mkdir('psix_output_10000/k_sensitivity')


##################################


#psix_single_lineage.run_psix(latent='processed_tables/pc2_rd.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=10)
#psix_single_lineage.psix_results.to_csv('psix_output_10000/k_sensitivity/k_0.001.tab.gz', sep='\t', index=True, header=True)


#psix_single_lineage.run_psix(latent='processed_tables/pc2_rd.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=100)
#psix_single_lineage.psix_results.to_csv('psix_output_10000/k_sensitivity/k_0.01.tab.gz', sep='\t', index=True, header=True)
#
#psix_single_lineage.run_psix(latent='processed_tables/pc2_rd.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=500)
#psix_single_lineage.psix_results.to_csv('psix_output_10000/k_sensitivity/k_0.05.tab.gz', sep='\t', index=True, header=True)


#psix_single_lineage.run_psix(latent='processed_tables/pc2_rd.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=1000)
#psix_single_lineage.psix_results.to_csv('psix_output_10000/k_sensitivity/k_0.1.tab.gz', sep='\t', index=True, header=True)

# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=1500)
# psix_single_lineage.psix_results.to_csv('psix_output_10000/k_sensitivity/k_0.15.tab.gz', sep='\t', index=True, header=True)

#psix_single_lineage.run_psix(latent='processed_tables/pc2_rd.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=2000)
#psix_single_lineage.psix_results.to_csv('psix_output_10000/k_sensitivity/k_0.2.tab.gz', sep='\t', index=True, header=True)


#psix_single_lineage.run_psix(latent='processed_tables/pc2_rd.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=2500)
#psix_single_lineage.psix_results.to_csv('psix_output_10000/k_sensitivity/k_0.25.tab.gz', sep='\t', index=True, header=True)

#psix_single_lineage.run_psix(latent='processed_tables/pc2_rd.tab.gz', n_jobs=25,
#                                n_random_exons=2000, n_neighbors=5000)
#psix_single_lineage.psix_results.to_csv('psix_output_10000/k_sensitivity/k_0.5.tab.gz', sep='\t', index=True, header=True)





#del psix_single_lineage





# psix_single_lineage = psix.Psix()
# psix_single_lineage.junctions2psi(
#         'processed_tables/splice_junction_counts_0.1_5000.tab.gz',
#         '',
#         'processed_tables/tpm_0.1_5000.tab.gz',
#         save_files_in='psix_output_5000/', dtype=float16)

# #psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
# #                                n_random_exons=2000, n_neighbors=500)
# #
# #psix_results = psix_single_lineage.psix_results
# #psix_single_lineage.save_psix_object(psix_dir = 'psix_output_5000/psix_object/')




# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=5)
# psix_single_lineage.psix_results.to_csv('psix_output_5000/k_sensitivity/k_0.001.tab.gz', sep='\t', index=True, header=True)


# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=50)
# psix_single_lineage.psix_results.to_csv('psix_output_5000/k_sensitivity/k_0.01.tab.gz', sep='\t', index=True, header=True)

# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=250)
# psix_single_lineage.psix_results.to_csv('psix_output_5000/k_sensitivity/k_0.05.tab.gz', sep='\t', index=True, header=True)


# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=500)
# psix_single_lineage.psix_results.to_csv('psix_output_5000/k_sensitivity/k_0.1.tab.gz', sep='\t', index=True, header=True)


# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=750)
# psix_single_lineage.psix_results.to_csv('psix_output_5000/k_sensitivity/k_0.15.tab.gz', sep='\t', index=True, header=True)


# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=1000)
# psix_single_lineage.psix_results.to_csv('psix_output_5000/k_sensitivity/k_0.2.tab.gz', sep='\t', index=True, header=True)


# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=1250)
# psix_single_lineage.psix_results.to_csv('psix_output_5000/k_sensitivity/k_0.25.tab.gz', sep='\t', index=True, header=True)

# #psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_5000.tab.gz', n_jobs=25,
# #                                n_random_exons=2000, n_neighbors=2500)
# #psix_single_lineage.psix_results.to_csv('psix_output_5000/k_sensitivity/k_0.5.tab.gz', sep='\t', index=True, header=True)




# del psix_single_lineage







# psix_single_lineage = psix.Psix()
# psix_single_lineage.junctions2psi(
#         'processed_tables/splice_junction_counts_0.1_1000.tab.gz',
#         '',
#         'processed_tables/tpm_0.1_1000.tab.gz',
#         save_files_in='psix_output_1000/', dtype=float16)

# #psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_1000.tab.gz', n_jobs=25,
# #                                n_random_exons=2000, n_neighbors=100)
# #
# #psix_results = psix_single_lineage.psix_results
# #psix_single_lineage.save_psix_object(psix_dir = 'psix_output_1000/psix_object/')


# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_1000.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=10)
# psix_single_lineage.psix_results.to_csv('psix_output_1000/k_sensitivity/k_0.01.tab.gz', sep='\t', index=True, header=True)

# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_1000.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=50)
# psix_single_lineage.psix_results.to_csv('psix_output_1000/k_sensitivity/k_0.05.tab.gz', sep='\t', index=True, header=True)


# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_1000.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=100)
# psix_single_lineage.psix_results.to_csv('psix_output_1000/k_sensitivity/k_0.1.tab.gz', sep='\t', index=True, header=True)

# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_1000.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=150)
# psix_single_lineage.psix_results.to_csv('psix_output_1000/k_sensitivity/k_0.15.tab.gz', sep='\t', index=True, header=True)

# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_1000.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=200)
# psix_single_lineage.psix_results.to_csv('psix_output_1000/k_sensitivity/k_0.2.tab.gz', sep='\t', index=True, header=True)


# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_1000.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=250)
# psix_single_lineage.psix_results.to_csv('psix_output_1000/k_sensitivity/k_0.25.tab.gz', sep='\t', index=True, header=True)

# #psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_1000.tab.gz', n_jobs=25,
# #                                n_random_exons=2000, n_neighbors=500)
# #psix_single_lineage.psix_results.to_csv('psix_output_1000/k_sensitivity/k_0.5.tab.gz', sep='\t', index=True, header=True)




# del psix_single_lineage








# psix_single_lineage = psix.Psix()
# psix_single_lineage.junctions2psi(
#         'processed_tables/splice_junction_counts_0.1_500.tab.gz',
#         '',
#         'processed_tables/tpm_0.1_500.tab.gz',
#         save_files_in='psix_output_500/', dtype=float16)

# #psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_500.tab.gz', n_jobs=25,
# #                                n_random_exons=2000, n_neighbors=50)
# #
# #psix_results = psix_single_lineage.psix_results
# #psix_single_lineage.save_psix_object(psix_dir = 'psix_output_500/psix_object/')





# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_500.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=5)
# psix_single_lineage.psix_results.to_csv('psix_output_500/k_sensitivity/k_0.01.tab.gz', sep='\t', index=True, header=True)

# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_500.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=25)
# psix_single_lineage.psix_results.to_csv('psix_output_500/k_sensitivity/k_0.05.tab.gz', sep='\t', index=True, header=True)


# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_500.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=50)
# psix_single_lineage.psix_results.to_csv('psix_output_500/k_sensitivity/k_0.1.tab.gz', sep='\t', index=True, header=True)


# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_500.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=75)
# psix_single_lineage.psix_results.to_csv('psix_output_500/k_sensitivity/k_0.15.tab.gz', sep='\t', index=True, header=True)


# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_500.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=100)
# psix_single_lineage.psix_results.to_csv('psix_output_500/k_sensitivity/k_0.2.tab.gz', sep='\t', index=True, header=True)


# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_500.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=125)
# psix_single_lineage.psix_results.to_csv('psix_output_500/k_sensitivity/k_0.25.tab.gz', sep='\t', index=True, header=True)

# #psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_500.tab.gz', n_jobs=25,
# #                                n_random_exons=2000, n_neighbors=250)
# #psix_single_lineage.psix_results.to_csv('psix_output_500/k_sensitivity/k_0.5.tab.gz', sep='\t', index=True, header=True)




# del psix_single_lineage






# psix_single_lineage = psix.Psix()
# psix_single_lineage.junctions2psi(
#         'processed_tables/splice_junction_counts_0.1_100.tab.gz',
#         '',
#         'processed_tables/tpm_0.1_100.tab.gz',
#         save_files_in='psix_output_100/', dtype=float16)

# #psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_100.tab.gz', n_jobs=25,
# #                                n_random_exons=2000, n_neighbors=10)
# #
# #psix_results = psix_single_lineage.psix_results
# #psix_single_lineage.save_psix_object(psix_dir = 'psix_output_100/psix_object/')





# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_100.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=5)
# psix_single_lineage.psix_results.to_csv('psix_output_100/k_sensitivity/k_0.05.tab.gz', sep='\t', index=True, header=True)


# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_100.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=10)
# psix_single_lineage.psix_results.to_csv('psix_output_100/k_sensitivity/k_0.1.tab.gz', sep='\t', index=True, header=True)

# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_100.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=15)
# psix_single_lineage.psix_results.to_csv('psix_output_100/k_sensitivity/k_0.15.tab.gz', sep='\t', index=True, header=True)

# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_100.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=20)
# psix_single_lineage.psix_results.to_csv('psix_output_100/k_sensitivity/k_0.2.tab.gz', sep='\t', index=True, header=True)


# psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_100.tab.gz', n_jobs=25,
#                                 n_random_exons=2000, n_neighbors=25)
# psix_single_lineage.psix_results.to_csv('psix_output_100/k_sensitivity/k_0.25.tab.gz', sep='\t', index=True, header=True)

# #psix_single_lineage.run_psix(latent='processed_tables/pc2_rd_100.tab.gz', n_jobs=25,
# #                                n_random_exons=2000, n_neighbors=50)
# #psix_single_lineage.psix_results.to_csv('psix_output_100/k_sensitivity/k_0.5.tab.gz', sep='\t', index=True, header=True)
