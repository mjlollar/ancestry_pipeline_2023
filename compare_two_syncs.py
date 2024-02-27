#### updated: MJL 01/29/24
#### Hey, do you need a script that uses really slow apply methods to loop through sync files and tell you % matches, mismatches, and minor allele matches?
#### Well do I have something for you

import argparse
import numpy as np
import pandas as pd
import sys

parser = argparse.ArgumentParser(description="Compare to sync files")
parser.add_argument('--a', help='Sync File One', required=True, type=str)
parser.add_argument('--b', help='Sync File Two', required=True, type=str)
args = parser.parse_args()

df_s1 = pd.read_csv(args.a, sep='\t', header=None)
df_s2 = pd.read_csv(args.b, sep='\t', header=None)

df_s1.columns = ['chr', 'pos', 'refA', 'countA']
df_s2.columns = ['chr', 'pos', 'refA', 'countB']


### Merge Sync files at common sites
df = pd.merge(df_s1[['pos', 'countA']], df_s2[['pos', 'countB']], on='pos')
df[['A1', 'T1', 'C1', 'G1', 'countN1', 'countDel1']] = df['countA'].str.split(':', expand=True)
df[['A2', 'T2', 'C2', 'G2', 'countN2', 'countDel2']] = df['countB'].str.split(':', expand=True)
df.drop(df.columns[[1,2,7,8,13,14]], axis=1, inplace=True) #drop unneeded columns(total counts, Indel and Del counts)

###Drop no read columns
df = df[(df['A1'].astype(int) + df['T1'].astype(int) + df['C1'].astype(int) + df['G1'].astype(int)).ge(1)]
df = df[(df['A2'].astype(int) + df['T2'].astype(int) + df['C2'].astype(int) + df['G2'].astype(int)).ge(1)]


allele_dict = {'0':'A','1':'T','2':'C','3':'G'}
### Function to determine Major and Minor Allele and add columns to panel
def majorminor_grouped(row, set_num): 
	if set_num == 1:
		A_count = int(row['A1'])
		T_count = int(row['T1'])
		C_count = int(row['C1'])
		G_count = int(row['G1'])
	elif set_num == 2:
		A_count = int(row['A2'])
		T_count = int(row['T2'])
		C_count = int(row['C2'])
		G_count = int(row['G2'])
	else:
		print("Sanity Error Print (function majorminor_grouped)")
	counts = np.array([A_count, T_count, C_count, G_count])
	count_index = np.argsort(counts)
	major = allele_dict.get(str(count_index[3]))
	#If minor count is less than 10% of major count sample, do not consider minor allele
	if counts[count_index[2]] <= round(0.1 * counts[count_index[3]]):
		minor = "N"
	else:
		minor = allele_dict.get(str(count_index[2]))
	return major, minor

### Get Major and Minor alleles at merged sites
try:
	#print("Major 1 calculate")
	df[["1A1", "1A2"]] = df.apply(lambda row: majorminor_grouped(row, 1), axis='columns', result_type='expand')
except ValueError: #catch instances where empty dataframe
	print("Empty Dataframe, Exiting....")
	sys.exit()
try:
	#print("Major 2 calculate")
	df[["2A1", "2A2"]] = df.apply(lambda row: majorminor_grouped(row, 2), axis='columns', result_type='expand')
except ValueError: #catch instances where empty dataframe
	print("Empty Dataframe, Exiting....")
	sys.exit()

### Compare Alleles
def allele_match(row):
	A1 = str(row['1A1'])
	a1 = str(row['1A2'])
	A2 = str(row['2A1'])
	a2 = str(row['2A2'])
	status = "N"

	if A1 == A2: #If major alleles match
		status = "MA"
		return status
	
	elif A1 != A2: #If major alleles don't match
		#Should be no N read in major column so don't need to check N match to minor as false count
		if A1 == a2: #If major 1 is equal to minor 2
			status = "MI"
			return status
		elif A1 != a2: #If major 1 doesn't match minor 2
			if A2 == a1: #If major 2 matches minor 1
				status = "MI"
				return status
			elif A2 != a1: #No matches to either major
				if a1 == a2: #If minor alleles match
					status == "MII"
					return status
				else: # NO matches return default N status
					return status
			else:
				print("Sanity print Error: allele match function loop 3")
		else:
			print("Sanity print Error: allele match function loop 2")	 
	else:
		print("Sanity print Error: allele match function loop 1")

try:
	#print("Comparing alleles")
	df["status"] = df.apply(lambda row: allele_match(row), axis='columns', result_type='expand')
except ValueError: #catch instances where empty dataframe
	print("Empty Dataframe, Exiting....")
	sys.exit()

total_sites = int(df.shape[0])
try:
	major_matches = int(df['status'].value_counts()['MA'])
except KeyError: #catch no counts of this value
	major_matches = 0
try:
	minor_matches = int(df['status'].value_counts()['MI'])
except KeyError: #catch no counts of this value
	minor_matches = 0
try:
	minorpair_matches = int(df['status'].value_counts()['MII'])
except KeyError: #catch no counts of this value
	minorpair_matches = 0

major_per = major_matches / total_sites * 100
minor_per = minor_matches / total_sites * 100
minorpair_per = minorpair_matches / total_sites * 100
mismatches = total_sites - major_matches - minor_matches - minorpair_matches
mismatch_per =  mismatches / total_sites * 100

#print("=========================================")
print("Major Allele Match Percentage: " + str(major_per) + " (" + str(major_matches) + " sites)")
print("Minor Allele Match Percentage: " + str(minor_per) + " (" + str(minor_matches) + " sites)")
print("Minor/Minor Allele Match Percentage: " + str(minorpair_per) + " (" + str(minorpair_matches) + " sites)")
print("Mismatch Percentage: " + str(mismatch_per) + " (" + str(mismatches) + " sites)")
