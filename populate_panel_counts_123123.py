#### MJL 12/31/23, added recombination 
#### SCRIPT TO ADJUST PANEL FILES WITH ADDED OFFSPRING COUNT COLUMNS
#### THIS SCRIPT IS A VARIANT OF "populate_snp_matrix.pl"
#### USES TOTAL COUNTS IN ADDED COLUMNS INSTEAD OF 1/0 PRESENCE/ABSENCE

import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description="new populate snp matrix script")
parser.add_argument('--p', help='Input Panel File', required=True, type=str)
parser.add_argument('--s', help='Input Sync File', required=True, type=str)
parser.add_argument('--c', help='Chromosome Arm', required=True, type=str)
args = parser.parse_args()

### Input to dataframes
df_panel = pd.read_csv(args.p, sep='\t', header=None)
df_sync = pd.read_csv(args.s, sep='\t', header=None)

### Add headers for simple combine
#A1 - allele1, A2 - allele2, 
#P1A1 - Population1/allele1, etc.., P2A2 - population2/allele2
df_panel.columns = ['chr', 'pos', 'A1', 'A2', 'P1A1', 'P1A2', 'P2A1', 'P2A2', 'recomb']
df_sync.columns = ['chr', 'pos', 'refA', 'count']

### Merge Sync file and Panel files
df_combined = pd.merge(df_panel, df_sync[['pos','count']], on='pos')
df_combined[['A', 'T', 'C', 'G', 'countN', 'countDel']] = df_combined['count'].str.split(':', expand=True)
df_combined.drop(df_combined.columns[[9,14,15]], axis=1, inplace=True) #drop unneeded columns(total counts, Indel and Del counts)

### Match offspring allele count with panel allele count
def allele_match(row):
	#match allele to count
	alelle_1 = str(row['A1'])
	allele_2 = str(row['A2'])
	count_a1 = int(row[alelle_1])
	count_a2 = int(row[allele_2])
	return(count_a1, count_a2)

### Generate new columns
df_combined[["A1f", "A2f"]] = df_combined.apply(lambda row: allele_match(row), axis='columns', result_type='expand')

### Drop rows with no allele counts in input sync (no information rows) and remove intermediate columns for output
index_group = df_combined[(df_combined['A1f'] == 0) & (df_combined['A2f'] == 0)].index
df_combined.drop(index_group, inplace=True)
df_combined.drop(df_combined.columns[[9,10,11,12]], axis=1, inplace=True)

### Print panel output
outchrom = args.c + '.ahmm_in.panel'
df_combined.to_csv(outchrom, sep='\t', header=False, index=False)
