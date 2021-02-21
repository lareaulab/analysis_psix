import numpy as np
import pandas as pd
import os
import subprocess as sp
from tqdm import tqdm 

ase_bed = pd.read_csv('/mnt/lareaulab/cfbuenabadn/Genomes/pipeline_files/mm10_skipped_exons_extended.bed',
                     names = ['chrom', 'start', 'end', 'intron', 'event', 'strand'], sep='\t')


# Exons in Tiklova neurogenesis lineage only


# event_list = [x.rstrip() for x in open('exons_list.txt', 'r').readlines()]

# fh = open('binding_sites.bed', 'w')

# for event in tqdm(event_list, position=0, leave=True):
#     i1 = ase_bed.loc[ase_bed.intron == event + '_I1']
#     i2 = ase_bed.loc[ase_bed.intron == event + '_I2']
    
#     e1 = (int(i1.start) -101, 
#           int(np.min([int(i1.start) + 100, 
#                       np.mean([int(i1.start), int(i1.end)])])))
    
#     e2 = (int(np.max([int(i2.end) - 100, 
#                       np.mean([int(i2.start), int(i2.end)])])), int(i2.end)+101)
#     se = (int(np.max([int(i1.end)-100, np.mean([int(i1.start), int(i1.end)])])), 
#           int(np.min([int(i2.start)+100, np.mean([int(i2.start), int(i2.end)])]))) # changed max to min
    
#     e1_row = '\t'.join([str(i1.loc[i1.index[0], 'chrom']), str(e1[0]), str(e1[1]), 
#                         event + '_e1', event, str(i1.loc[i1.index[0], 'strand'])]) + '\n'
#     se_row = '\t'.join([str(i1.loc[i1.index[0], 'chrom']), str(se[0]), str(se[1]), event + '_se', 
#                         event, str(i1.loc[i1.index[0], 'strand'])]) + '\n'
#     e2_row = '\t'.join([str(i1.loc[i1.index[0], 'chrom']), str(e2[0]), str(e2[1]), event + '_e2', 
#                         event, str(i1.loc[i1.index[0], 'strand'])]) + '\n'
    
#     fh.write(e1_row + se_row + e2_row)
    
# fh.close()


####### 5'-3' and separate SE into junctions

event_list = [x.rstrip() for x in open('exons_list.txt', 'r').readlines()]

fh = open('binding_sites.bed', 'w')

for event in tqdm(event_list, position=0, leave=True):
    i1 = ase_bed.loc[ase_bed.intron == event + '_I1']
    i2 = ase_bed.loc[ase_bed.intron == event + '_I2']
    
    e1 = (int(i1.start) -101, 
          int(np.min([int(i1.start) + 100, 
                      np.mean([int(i1.start), int(i1.end)])])))
    
    e2 = (int(np.max([int(i2.end) - 100, 
                      np.mean([int(i2.start), int(i2.end)])])), int(i2.end)+101)
    
    se1 = (int(np.max([int(i1.end)-100, np.mean([int(i1.start), int(i1.end)])])), 
          int(np.min([int(i1.end)+100, np.mean([int(i1.end), int(i2.start)])])))
    
    se2 = (int(np.max([int(i2.start)-100, np.mean([int(i1.end), int(i2.start)])])), 
          int(np.min([int(i2.start)+100, np.mean([int(i2.start), int(i2.end)])])))
    
    if str(i1.loc[i1.index[0], 'strand']) == '+':
        
    
        e1_row = '\t'.join([str(i1.loc[i1.index[0], 'chrom']), str(e1[0]), str(e1[1]), event + '_e1', 
                            event, str(i1.loc[i1.index[0], 'strand'])]) + '\n'
        se1_row = '\t'.join([str(i1.loc[i1.index[0], 'chrom']), str(se1[0]), str(se1[1]), event + '_s1', 
                            event, str(i1.loc[i1.index[0], 'strand'])]) + '\n'
        se2_row = '\t'.join([str(i1.loc[i1.index[0], 'chrom']), str(se2[0]), str(se2[1]), event + '_s2', 
                            event, str(i1.loc[i1.index[0], 'strand'])]) + '\n'
        e2_row = '\t'.join([str(i1.loc[i1.index[0], 'chrom']), str(e2[0]), str(e2[1]), event + '_e2', 
                            event, str(i1.loc[i1.index[0], 'strand'])]) + '\n'
        
    elif str(i1.loc[i1.index[0], 'strand']) == '-':
        e1_row = '\t'.join([str(i1.loc[i1.index[0], 'chrom']), str(e2[0]), str(e2[1]), event + '_e1', 
                            event, str(i1.loc[i1.index[0], 'strand'])]) + '\n'
        se1_row = '\t'.join([str(i1.loc[i1.index[0], 'chrom']), str(se2[0]), str(se2[1]), event + '_s1', 
                            event, str(i1.loc[i1.index[0], 'strand'])]) + '\n'
        se2_row = '\t'.join([str(i1.loc[i1.index[0], 'chrom']), str(se1[0]), str(se1[1]), event + '_s2', 
                            event, str(i1.loc[i1.index[0], 'strand'])]) + '\n'
        e2_row = '\t'.join([str(i1.loc[i1.index[0], 'chrom']), str(e1[0]), str(e1[1]), event + '_e2', 
                            event, str(i1.loc[i1.index[0], 'strand'])]) + '\n'
    else:
        raise Exception('strand error')

    fh.write(e1_row + se1_row + se2_row + e2_row)
    
fh.close()
