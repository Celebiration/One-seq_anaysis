#!/usr/local/anaconda/bin/python
import sys
import pandas as pd
import numpy as np
import subprocess
import os
import multiprocessing
from Bio import SeqIO

#inputs
file1="cas_offinder_6_0_0.out"
file2="cas_offinder_4_2_2.out"

def merge_hits(file_list):
	tmp=[]
	merged=pd.DataFrame()
	for ff in file_list:
		dd=pd.read_csv(ff,sep="\t")
		for i in range(len(dd)):
			ind=0
			for j in range(-5,6):
				if (dd["Chromosome"][i],dd["Direction"][i],dd["Position"][i]+j) in tmp:
					ind=1
					index=tmp.index((dd["Chromosome"][i],dd["Direction"][i],dd["Position"][i]+j))
					#比较，保留bulge size低的那一个，若相同，则保留任意一个
					if merged["Bulge Size"][index] > dd["Bulge Size"][i]:
						merged.iloc[index]=dd.iloc[i]
					break
			if ind==0:
				merged=pd.concat([merged,dd.iloc[i:(i+1)]],ignore_index=True)
				tmp.append((dd["Chromosome"][i],dd["Direction"][i],dd["Position"][i]))
	return(merged)

merged=merge_hits([file1,file2])
merged.to_csv("merged_6.tsv",sep="\t",index=False)

#sort
chr_mapping={"chr1":1,"chr2":2,"chr3":3,"chr4":4,"chr5":5,"chr6":6,"chr7":7,"chr8":8,"chr9":9,"chr10":10,"chr11":11,"chr12":12,"chr13":13,"chr14":14,"chr15":15,"chr16":16,"chr17":17,"chr18":18,"chr19":19,"chr20":20,"chr21":21,"chr22":22,"chrX":23,"chrY":24}
def custom_sort(value):
	return chr_mapping.get(value.split(' ')[0], float('inf'))
result=merged
result['custom_sorted'] = result['Chromosome'].map(custom_sort)
result.sort_values(by=['custom_sorted',"Chromosome","Position"], ascending=True,inplace=True)
result.drop(columns=["custom_sorted"],inplace=True)
result.to_csv("merged_6.sorted.tsv",sep="\t",index=False)
