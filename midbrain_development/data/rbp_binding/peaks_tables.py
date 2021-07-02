import pandas as pd
import numpy as np
import subprocess as sp
import os
from matplotlib import pyplot as plt
from tqdm import tqdm

def get_clip_peaks(exons_bed, rbp_bed, min_overlap=4):
    
    intersect_cmd = 'bedtools intersect -s -a {exons_bed} -b {rbp_bed} > temp.bed'
    intersect_cmd = intersect_cmd.format(exons_bed=exons_bed, rbp_bed=rbp_bed)
    sp.run(intersect_cmd, shell=True)
    
    peaks = pd.read_csv('temp.bed', sep='\t', 
                        names=['chrom', 'start', 'end', 'position', 'exon', 'strand']).sort_values(['position', 'start'])
    
    peaks_df = pd.DataFrame()
    chrom = []
    start = []
    end = []
    position = []
    exon = []
    strand = []
    
    current_exon = ''
    for idx, rows in peaks.iterrows():
        overlaps = False
        if current_exon == rows.position:
            
            if (int(rows.start) >= current_start) and (int(rows.start) <= current_end):
                overlaps = True
            if (int(rows.end) >= current_start) and (int(rows.end) <= current_end):
                overlaps = True
            if (current_start <= int(rows.start)) and (current_end >= int(rows.end)):
                overlaps = True
                
            if (current_start >= int(rows.start)) and (current_end <= int(rows.end)):
                overlaps = True
        if overlaps:

            start[-1] = np.min([int(rows.start), current_start])
            end[-1] = np.max([int(rows.end), current_end])
            
        else:
            chrom.append(rows.chrom)
            start.append(int(rows.start))
            end.append(int(rows.end))
            position.append(rows.position)
            exon.append(rows.exon)
            strand.append(rows.strand)
            
        current_exon = position[-1]
        current_start = start[-1]
        current_end = end[-1]
                
    peaks_df['chrom'] = chrom
    peaks_df['start'] = start
    peaks_df['end'] = end
    peaks_df['position'] = position
    peaks_df['exon'] = exon
    peaks_df['strand'] = strand
        
    peaks_df = peaks_df.loc[(peaks_df.end - peaks_df.start) >= min_overlap]
    
    overlap_len = peaks_df.end - peaks_df.start
    peaks_df['overlap_len'] = overlap_len
    
    
    sp.run('rm temp.bed', shell=True)
    
    return peaks_df


def get_clip_table(exon_list, rbp_list, rbp_dir, exon_bed, min_overlap=4, report_len=False):
    clip_df = pd.DataFrame(np.zeros((len(exon_list), len(rbp_list))))
    clip_df.index = exon_list
    clip_df.columns = rbp_list
    
    bed_df = pd.read_csv(exon_bed, sep='\t', names=['chrom', 'start', 'end', 'exon', 'ase', 'strand'])
    bed_df.index = bed_df.exon
    
    for rbp in tqdm(rbp_list, position=0, leave=True):
        rbp_bed = '/'.join([rbp_dir, rbp, rbp + '.bed'])
        peaks = get_clip_peaks(exon_bed, rbp_bed, min_overlap=min_overlap).drop_duplicates(keep='first')
        
        if len(peaks) > 0:
            for idx, row in peaks.iterrows():
                if report_len:
                    exon_len = int(bed_df.loc[row.position].end) - int(bed_df.loc[row.position].start)
                    clip_df.loc[row.position, rbp] += row.overlap_len/exon_len
                else:
                    clip_df.loc[row.position, rbp] += 1
    return clip_df


with open('exons_list.txt', 'r') as fh:
    exons_list = [x.rstrip() for x in fh.readlines()]
    
event_list = [x+'_e1' for x in exons_list]
event_list += [x+'_s1' for x in exons_list]
event_list += [x+'_s2' for x in exons_list]
event_list += [x+'_e2' for x in exons_list]
    
import os
rbp_list = os.listdir('/mnt/lareaulab/cfbuenabadn/Network/Mouse/CLIPSeq/merged_datasets/peaks_mm10/')

clip_table = get_clip_table(event_list, rbp_list, 
               '~/Network/Mouse/CLIPSeq/merged_datasets/peaks_mm10/', 'binding_sites.bed', min_overlap=4)
clip_table.to_csv('peaks_clip_tags.tab', sep='\t', index=True, header=True)
            
    
clip_table_len = get_clip_table(event_list, rbp_list,
               '~/Network/Mouse/CLIPSeq/merged_datasets/peaks_mm10/', 'binding_sites.bed', min_overlap=4,report_len=True)
clip_table_len.to_csv('peaks_clip_tags_overlap_len.tab', sep='\t', index=True, header=True)
