### MJL 02/08/2024
### Takes two parental sync files and generates a panel
### Optional flag to adjust minimum read count in parental sync for consideration of a site

import sys
import argparse
import numpy as np
import pandas as pd
import operator

parser = argparse.ArgumentParser(description="Takes parental and sample sync files and generates ahmm panel input")
parser.add_argument("--p1", help="Parental Sync File 1", required=True, type=str)
parser.add_argument("--p2", help="Parental Sync File 2", required=True, type=str)
parser.add_argument("--prc", help="Per Parent Minimum Total Read Count at Panel Site (Default=1)", required=True, default=1, type=int)
parser.add_argument("--c", help="Chromosome Arm", required=True, type=str)
args = parser.parse_args()

### Load in data
print("Loading Syncs....")
df_par1 = pd.read_csv(args.p1, sep='\t', header=None)
df_par2 = pd.read_csv(args.p2, sep='\t', header=None)
df_par1.columns = ['chr', 'pos', 'refA', 'p1_count']
df_par2.columns = ['chr', 'pos', 'refA', 'p2_count']

### Expand sync count columns
df_par1[['p1A', 'p1T', 'p1C', 'p1G', 'p1countN', 'p1countDel']] = df_par1['p1_count'].str.split(':', expand=True)
df_par2[['p2A', 'p2T', 'p2C', 'p2G', 'p2countN', 'p2countDel']] = df_par2['p2_count'].str.split(':', expand=True)

### Find overlapping read sites and marge syncs to common dataframe
print("Merging overlapping read sites in parental syncs...")
df_panel = pd.merge(df_par1, df_par2, on='pos')
df_panel.drop(df_panel.columns[[2,3,8,9,10,11,12,17,18]], axis=1, inplace=True)

### Drop rows with <# counts specified with input (per parent)
print("Dropping Reads with Threshold: >=" + str(args.prc))
df_panel = df_panel[(df_panel['p1A'].astype(int) + df_panel['p1T'].astype(int) + df_panel['p1C'].astype(int) + df_panel['p1G'].astype(int)).ge(args.prc)]
df_panel = df_panel[(df_panel['p2A'].astype(int) + df_panel['p2T'].astype(int) + df_panel['p2C'].astype(int) + df_panel['p2G'].astype(int)).ge(args.prc)]

allele_dict = {'0':'A','1':'T','2':'C','3':'G'}
### Function to determine Major and Minor Allele and add columns to panel
def majorminor_grouped(row): 
	A_count = int(row['p1A']) + int(row['p2A'])
	T_count = int(row['p1T']) + int(row['p2T'])
	C_count = int(row['p1C']) + int(row['p2C'])
	G_count = int(row['p1G']) + int(row['p2G'])

	counts = np.array([A_count, T_count, C_count, G_count])
	count_index = np.argsort(counts)
	major = allele_dict.get(str(count_index[3]))
	#If minor count is less than 10% of major count sample, do not consider minor allele
	if counts[count_index[2]] <= round(0.1 * counts[count_index[3]]):
		minor = "N"
	else:
		minor = allele_dict.get(str(count_index[2]))
	return major, minor

### Function to Count Parental alleles and add four new columns for panel output
def count_parental_alleles(row):
	par1_alleles = np.array([int(row['p1A']),int(row['p1T']),int(row['p1C']),int(row['p1G'])])
	par2_alleles = np.array([int(row['p2A']),int(row['p2T']),int(row['p2C']),int(row['p2G'])])
	p1_count_index = np.argsort(par1_alleles)
	p2_count_index = np.argsort(par2_alleles)

	### Get allele ID of parental allele for A1 and A2 allele count
	p1_A1 = allele_dict.get(str(p1_count_index[3]))
	p2_A1 = allele_dict.get(str(p2_count_index[3]))
	p1_A2 = allele_dict.get(str(p1_count_index[2]))
	p2_A2 = allele_dict.get(str(p2_count_index[2]))

	### If A2 count represents less than 10% of reads within parent do not consider as minor read
	### Catches zero counts in minor to exclude as well (zero never 10% of any positive value)
	if par1_alleles[p1_count_index[2]] <= round(0.1 * par1_alleles[p1_count_index[3]]):
		p1_A2 = "N"
	else:
		pass
	if par2_alleles[p2_count_index[2]] <= round(0.1 * par2_alleles[p2_count_index[3]]):
		p2_A2 = "N"
	else:
		pass

	### Store counts of alleles for each population (Input Parent 1 will represent Pop1, or first two count columns)
	## Parent 1
	## If A1/A2 parental allele matches A1(major) allele, count
	if (p1_A1 == row['A1']) | (p1_A2 == row['A1']):
		p1A1 = 1
	else:
		p1A1 = 0
	## If A1/A2 parental allele matches A2(minor) allele, count
	if (p1_A1 == row['A2']) or (p1_A2 == row['A2']):
		p1A2 = 1
	else:
		p1A2 = 0
	## Parent2
	if (p2_A1 == row['A1']) or (p2_A2 == row['A1']):
		p2A1 = 1
	else:
		p2A1 = 0
	if (p2_A1 == row['A2']) or (p2_A2 == row['A2']):
		p2A2 = 1
	else:
		p2A2 = 0

	return p1A1, p1A2, p2A1, p2A2

### Get Major and Minor alleles at merged sites
try:
	print("Determining Major/Minor Alleles and Counts...")
	df_panel[["A1", "A2"]] = df_panel.apply(lambda row: majorminor_grouped(row), axis='columns', result_type='expand')
except ValueError: #catch instances where empty dataframe
	print("Emtpy Dataframe Detected. Are your read cutoffs too strict? (--prc flag)")
	print("Exiting....")
	sys.exit()

### Drop rows with no minor allele (Represents same allele in both parents, not informative)
df_panel = df_panel[df_panel["A2"]!="N"]

### Get Parental Counts
try:
	print("Gettting Parental Counts....")
	df_func = df_panel.apply(lambda row: count_parental_alleles(row), axis='columns', result_type='expand')
	df_func.columns = ['P1A1', 'P1A2', 'P2A1', 'P2A2']
	df_panel = pd.concat([df_panel, df_func], axis='columns')
except ValueError: #catch instances where empty dataframe, unlikely here
	print("Emtpy Dataframe Detected. Insufficient minor counts in panel")
	print("Exiting....")
	sys.exit()

### Drop rows with both Major/Minor counted (all 1's in added four columns)
df_panel = df_panel[(df_panel['P1A1'].astype(int) + df_panel['P1A2'].astype(int) + df_panel['P2A1'].astype(int) + df_panel['P2A2'].astype(int)).ne(4)]
### Drop rows not needed in panel output
df_panel.drop(df_panel.columns[[2,3,4,5,6,7,8,9]], axis=1, inplace=True)

### Add dummy recombination columns
df_panel['recomb'] = 0.00000001
### Replace chromosome code with corrected arm ID
df_panel['chr_x'] = args.c

### Output to Dataframe
try:
	print("Creating panel file output...")
	outname = args.c + ".panel"
	df_panel.to_csv(outname, sep='\t', index=False, header=None)
except ValueError: #catch instances where empty dataframe
	print("Emtpy Dataframe Detected. Insufficient informative sites to make panel")
	print("Exiting....")
	sys.exit()

