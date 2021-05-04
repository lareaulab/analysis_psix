DATA='/mnt/lareaulab/cfbuenabadn/data_sc_regulation/data_autocorrelation'

#python ~/psix/utils/psix.py -psi $DATA/tiklova_neurogenesis/skipped_exons_psi.tab -mrna $DATA/tiklova_neurogenesis/mrna_per_event.tab -rd $DATA/tiklova_neurogenesis/rd_pc2.tab -o tiklova_neurogenesis -k 100

#python ~/psix/utils/psix.py -psi $DATA/tiklova_neurogenesis/skipped_exons_psi.tab -mrna $DATA/tiklova_neurogenesis/mrna_per_event.tab -rd ~/analysis_psix/scvi_runs/tiklova_scvi3_rd.tab -o scvi3 -k 100

#python ~/psix/utils/psix.py -psi $DATA/tiklova_neurogenesis/skipped_exons_psi.tab -mrna $DATA/tiklova_neurogenesis/mrna_per_event.tab -rd ~/analysis_psix/scvi_runs/tiklova_scvi5_batch_rd.tab -o scvi5_batch -k 100

python ~/psix/utils/psix.py -psi $DATA/tiklova_neurogenesis/skipped_exons_psi.tab -mrna $DATA/tiklova_neurogenesis/mrna_per_event.tab -rd ~/analysis_psix/scvi_runs/tiklova_scvi5_rd.tab -o scvi5 -k 100
