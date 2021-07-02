#brie-count -a ~/brie_gencode.vM10/mm10_tiklova_ase.gtf -S brie2_files/bam_E13_v_E15.tsv -o brie_output_E13_v_E15 -p 20
#brie-count -a ~/brie_gencode.vM10/mm10_tiklova_ase.gtf -S brie2_files/bam_E13_v_E18.tsv -o brie_output_E13_v_E18 -p 20
#brie-count -a ~/brie_gencode.vM10/mm10_tiklova_ase.gtf -S brie2_files/bam_E13_v_P1.tsv -o brie_output_E13_v_P1 -p 20
#brie-count -a ~/brie_gencode.vM10/mm10_tiklova_ase.gtf -S brie2_files/bam_E13_v_P7.tsv -o brie_output_E13_v_P7 -p 20
#brie-count -a ~/brie_gencode.vM10/mm10_tiklova_ase.gtf -S brie2_files/bam_E13_v_P90.tsv -o brie_output_E13_v_P90 -p 20

brie-quant -i brie_output_E13_v_E15/brie_count.h5ad -o brie_output_E13_v_E15/brie_quant_cell.h5ad -c brie2_files/cells_E13_v_E15.tsv --interceptMode gene --LRTindex=All
brie-quant -i brie_output_E13_v_E18/brie_count.h5ad -o brie_output_E13_v_E18/brie_quant_cell.h5ad -c brie2_files/cells_E13_v_E18.tsv --interceptMode gene --LRTindex=All
brie-quant -i brie_output_E13_v_P1/brie_count.h5ad -o brie_output_E13_v_P1/brie_quant_cell.h5ad -c brie2_files/cells_E13_v_P1.tsv --interceptMode gene --LRTindex=All
brie-quant -i brie_output_E13_v_P7/brie_count.h5ad -o brie_output_E13_v_P7/brie_quant_cell.h5ad -c brie2_files/cells_E13_v_P7.tsv --interceptMode gene --LRTindex=All
brie-quant -i brie_output_E13_v_P90/brie_count.h5ad -o brie_output_E13_v_P90/brie_quant_cell.h5ad -c brie2_files/cells_E13_v_P90.tsv --interceptMode gene --LRTindex=All
