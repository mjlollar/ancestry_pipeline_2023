#### updated: MJL 01/29/24
#### SCRIPT TO ADJUST PANEL FILES WITH ADDED OFFSPRING COUNT COLUMNS
#### THIS SCRIPT IS A VARIANT OF "populate_snp_matrix.pl"
#### USES TOTAL COUNTS IN ADDED COLUMNS INSTEAD OF 1/0 PRESENCE/ABSENCE

import argparse
import numpy as np
import pandas as pd
import operator

parser = argparse.ArgumentParser(description="new populate snp matrix script")
parser.add_argument('--p', help='Input Panel File', required=True, type=str)
parser.add_argument('--s', help='Input Sync File', required=True, type=str)
parser.add_argument('--c', help='Chromosome Arm', required=True, type=str)
parser.add_argument('--rc', help='Minimum Read Count at sites used in ahmm_in.panel', required=True, type=int)
parser.add_argument('--off', help='Offset value of sync position coordinates for 0/1 index checking (default -1 as sync usually 1 indexed)', required=True, default=-1, type=int)
args = parser.parse_args()

print("Generating ahmm panel inputs with " + str(args.rc) + " read threshold")
print("Sync positions offset by: " + str(args.off))

### Input to dataframes
df_panel = pd.read_csv(args.p, sep='\t', header=None)
df_sync = pd.read_csv(args.s, sep='\t', header=None)

### Add headers for simple combine
#A1 - allele1, A2 - allele2, 
#P1A1 - Population1/allele1, etc.., P2A2 - population2/allele2
df_panel.columns = ['chr', 'pos', 'A1', 'A2', 'P1A1', 'P1A2', 'P2A1', 'P2A2', 'recomb']
df_sync.columns = ['chr', 'pos', 'refA', 'count']

df_panel['pos'] = df_panel['pos'] + args.off #set position by offset
#df_sync['pos'] = df_sync['pos'] + args.off #set position by offset

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

### Drop rows with no allele counts in input sync (no information rows)
#index_group = df_combined[(df_combined['A1f'] == 0) & (df_combined['A2f'] == 0)].index
#df_combined.drop(index_group, inplace=True)

### Drop rows with <# counts specified with input
index_group = df_combined[((df_combined['A'].astype(int) + df_combined['T'].astype(int) + df_combined['C'].astype(int) + df_combined['G'].astype(int)) < args.rc)].index
df_combined.drop(index_group, inplace=True)

#Drop columns not used by ancestry_hmm
df_combined.drop(df_combined.columns[[2,3,9,10,11,12]], axis=1, inplace=True)

### Print panel output
outchrom = args.c + '.ahmm_in.panel.txt'
df_combined.to_csv(outchrom, sep='\t', header=False, index=False)
