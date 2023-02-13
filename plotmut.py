#import sys,os,subprocess
import argparse
import allel
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('vcf_file')
parser.add_argument('outstem')
args = parser.parse_args()

vcf = allel.read_vcf(args.vcf_file, numbers={'GT': 2, 'ALT':1})

outstem = args.outstem

#print(vcf['samples'])
gt = allel.GenotypeArray(vcf['calldata/GT'])
#print(gt)

# obtain calls for the second sample in all variants
#print(gt[:, 0])

#print(str(gt[0, 1]))

print(f'working on file: {args.vcf_file}\n')

samples = vcf['samples'] 
print(f'list of samples: {samples}\n')

ref = vcf['variants/REF']
#print(ref)

alt = vcf['variants/ALT']
#print(alt)

l = len(gt)
print(f'number of SNPs = {l}\n')


homref = '[0 0]'
het = '[0 1]'
homalt = '[1 1]'

df = pd.DataFrame(samples)

dfr = pd.DataFrame(samples)

homreflist = np.empty(len(samples), dtype=object)
hetlist = np.empty(len(samples), dtype=object)
homaltlist = np.empty(len(samples), dtype=object)

tis = 'AG', 'GA', 'CT', 'TC'

tvs = 'AC', 'CA', 'GT', 'TG', 'AT','TA', 'GC', 'CG'

for ti in tis:
	exec(f'het{ti}list = np.empty(len(samples), dtype=object)')
	exec(f'hom{ti}list = np.empty(len(samples), dtype=object)')

for tv in tvs:
	exec(f'het{tv}list = np.empty(len(samples), dtype=object)')
	exec(f'hom{tv}list = np.empty(len(samples), dtype=object)')

hettilist = np.empty(len(samples), dtype=object)
hettvlist = np.empty(len(samples), dtype=object)

homtilist = np.empty(len(samples), dtype=object)
homtvlist = np.empty(len(samples), dtype=object)

possumlist = np.empty(len(samples), dtype=object)

for i in range(len(samples)):

	print("\nsample",f'{i+1} /',len(samples))
	print(samples[i])

	homrefcount = 0

	#Ti
	hetAGcount, hetGAcount, hetCTcount, hetTCcount = 0, 0, 0, 0
	homAGcount, homGAcount, homCTcount, homTCcount = 0, 0, 0, 0

	#Tv
	hetACcount, hetCAcount, hetGTcount, hetTGcount, hetATcount, hetTAcount, hetGCcount, hetCGcount = 0, 0, 0, 0, 0, 0, 0, 0
	homACcount, homCAcount, homGTcount, homTGcount, homATcount, homTAcount, homGCcount, homCGcount = 0, 0, 0, 0, 0, 0, 0, 0

	for j in range(len(ref)):

		if str(gt[j, i]) == homref:
			if ref[j] == 'A' or ref[j] == 'G' or ref[j] == 'C' or ref[j] == 'T':
				if alt[j] == 'A' or alt[j] == 'G' or alt[j] == 'C' or alt[j] == 'T':

					homrefcount += 1

		elif str(gt[j, i]) == het:

			#Ti
		    if ref[j] == 'A' and alt[j] == 'G':
		        hetAGcount += 1
		    elif ref[j] == 'G' and alt[j] == 'A':
		        hetGAcount += 1
		    elif ref[j] == 'C' and alt[j] == 'T':
			    hetCTcount += 1
		    elif ref[j] == 'T' and alt[j] == 'C':
		    	hetTCcount += 1

			#Tv
		    elif ref[j] == 'A' and alt[j] == 'C':
		    	hetACcount += 1
		    elif ref[j] == 'C' and alt[j] == 'A':
		    	hetCAcount += 1
		    elif ref[j] == 'G' and alt[j] == 'T':
		    	hetGTcount += 1
		    elif ref[j] == 'T' and alt[j] == 'G':
	    		hetTGcount += 1
		    elif ref[j] == 'A' and alt[j] == 'T':
		        hetATcount += 1
		    elif ref[j] == 'T' and alt[j] == 'A':
			    hetTAcount += 1
		    elif ref[j] == 'G' and alt[j] == 'C':
		    	hetGCcount += 1
		    elif ref[j] == 'C' and alt[j] == 'G':
		    	hetCGcount += 1

		elif str(gt[j, i]) == homalt:

			#Ti
			if ref[j] == 'A' and alt[j] == 'G':
				homAGcount += 1
			elif ref[j] == 'G' and alt[j] == 'A':
				homGAcount += 1
			elif ref[j] == 'C' and alt[j] == 'T':
				homCTcount += 1
			elif ref[j] == 'T' and alt[j] == 'C':
				homTCcount += 1

			#Tv
			elif ref[j] == 'A' and alt[j] == 'C':
				homACcount += 1
			elif ref[j] == 'C' and alt[j] == 'A':
				homCAcount += 1
			elif ref[j] == 'G' and alt[j] == 'T':
				homGTcount += 1
			elif ref[j] == 'T' and alt[j] == 'G':
				homTGcount += 1
			elif ref[j] == 'A' and alt[j] == 'T':
				homATcount += 1
			elif ref[j] == 'T' and alt[j] == 'A':
				homTAcount += 1
			elif ref[j] == 'G' and alt[j] == 'C':
				homGCcount += 1
			elif ref[j] == 'C' and alt[j] == 'G':
				homCGcount += 1

		hetticount = hetAGcount + hetGAcount + hetCTcount + hetTCcount
		hettvcount = hetACcount + hetCAcount + hetGTcount + hetTGcount + hetATcount + hetTAcount + hetGCcount + hetCGcount
		hetcount = hetticount + hettvcount

		homticount = homAGcount + homGAcount + homCTcount + homTCcount
		homtvcount = homACcount + homCAcount + homGTcount + homTGcount + homATcount + homTAcount + homGCcount + homCGcount
		homaltcount = homticount + homtvcount

		possum = homrefcount + hetcount + homaltcount

    #general
	possumlist[i] = possum

	homreflist[i] = homrefcount
	hetlist[i] = hetcount
	homaltlist[i] = homaltcount

	hettilist[i] = hetticount
	hettvlist[i] = hettvcount

	homtilist[i] = homticount
	homtvlist[i] = homtvcount
	
	#Ti
	for ti in tis:
		exec(f'het{ti}list[i] = het{ti}count')
		exec(f'hom{ti}list[i] = hom{ti}count')
	#Tv
	for tv in tvs:
		exec(f'het{tv}list[i] = het{tv}count')
		exec(f'hom{tv}list[i] = hom{tv}count')

	print(f'hom ref count={homrefcount}')
	print(f'het count={hetcount}')
	print(f'\thet Ti count={hetticount}')
	print(f'\thet Tv count={hettvcount}')
	print(f'hom alt count={homaltcount}')
	print(f'\thom Ti count={homticount}')
	print(f'\thom Tv count={homtvcount}')

df["sum"] = possumlist
df["homref"] = homreflist
df["het"] = hetlist
df["homalt"] = homaltlist
df["hetTi"] = hettilist
df["hetTv"] = hettvlist
df["homTi"] = homtilist
df["homTv"] = homtvlist

#Ti
for ti in tis:
	exec(f'df["het{ti}"] = het{ti}list')
#Tv
for tv in tvs:
	exec(f'df["het{tv}"] = het{tv}list')

for ti in tis:
	exec(f'df["hom{ti}"] = hom{ti}list')

for tv in tvs:
	exec(f'df["hom{tv}"] = hom{tv}list')

#	df = pd.DataFrame({homaltcount, hetcount, homrefcount})
#print(homreflist)
#print(hetlist)
#print(homaltlist)

#print(df)

mutsum = df["hetTi"] + df["hetTv"] + df["homTi"] + df["homTv"]

dfr["homref"] = df["homref"] / df["sum"]
dfr["het"] = df["het"] / df["sum"]
dfr["homalt"] = df["homalt"] / df["sum"]
dfr["hetTi"] = df["hetTi"] / mutsum
dfr["hetTv"] = df["hetTv"] / mutsum
dfr["homTi"] = df["homTi"] / mutsum
dfr["homTv"] = df["homTv"] / mutsum

#Ti
for ti in tis:
	exec(f'dfr["het{ti}"] =  df["het{ti}"] / df["sum"]')
#Tv
for tv in tvs:
	exec(f'dfr["het{tv}"] =  df["het{tv}"] / df["sum"]')

for ti in tis:
	exec(f'dfr["hom{ti}"] =  df["hom{ti}"] / df["sum"]')

for tv in tvs:
	exec(f'dfr["hom{tv}"] =  df["hom{tv}"] / df["sum"]')

#save tables
print("Saving counts to table")
df.to_csv(f'{outstem}_count.tsv', index=False, sep="\t")
dfr.to_csv(f'{outstem}_ratio.tsv', index=False, sep="\t")

print("Plotting")
df.plot(x=0, y=["sum", "homref", "het", "homalt"], kind="bar", width = 0.9, figsize=(len(samples)/3, 8))
plt.legend(loc= 'upper left')
plt.savefig(f'{outstem}_counts.pdf', dpi=300, bbox_inches = "tight")

dfr.plot(x=0, y=["homref", "het", "homalt"], kind="bar", width = 0.9, figsize=(len(samples)/3, 8))
plt.legend(loc= 'upper left')
plt.savefig(f'{outstem}_hom_het_ratio.pdf', dpi=300, bbox_inches = "tight")

dfr.plot(x=0, y=["hetTi", "homTi", "hetTv", "homTv"], kind="bar", width = 0.9, figsize=(len(samples)/3, 8))
plt.legend(loc= 'upper left')
plt.savefig(f'{outstem}_Ti_Tv_ratio.pdf', dpi=300, bbox_inches = "tight")

dfr.plot(x=0, y=["hetAG", "hetGA", "hetCT", "hetTC", "homAG", "homGA", "homCT", "homTC"], kind="bar", width = 0.9, figsize=(len(samples)/2, 8))
plt.legend(loc= 'upper left')
plt.savefig(f'{outstem}_Tis_ratio.pdf', dpi=300, bbox_inches = "tight")

dfr.plot(x=0, y=["hetAC", "hetCA", "hetGT", "hetTG", "hetAT", "hetTA", "hetGC", "hetCG", "homAC", "homCA", "homGT", "homTG", "homAT", "homTA", "homGC", "homCG"], kind="bar", width = 0.9, figsize=(len(samples), 8))
plt.legend(loc= 'upper left')
plt.savefig(f'{outstem}_Tvs_ratio.pdf', dpi=300, bbox_inches = "tight")
