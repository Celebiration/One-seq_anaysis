#!/usr/local/anaconda/bin/python
import pandas as pd
import numpy as np
import os
import random
from collections import Counter
import sys
from Bio import SeqIO

# 读取基因组序列文件
genome_file = "/Yol_Data/resources/Gencode_human_ref/GRCh38.p13.genome.fa"

# 创建一个空字典用于存储基因组序列
GRCh38 = {}

# 使用 BioPython 读取基因组序列文件，并存储到字典中
with open(genome_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        chromosome_name = record.id
        sequence = str(record.seq)  # 将序列转换为字符串形式
        GRCh38[chromosome_name] = sequence

def extract_genome_seq(chr,start,end):
	return(GRCh38.get(chr)[start:end])

comp_dict={'A':'T','T':'A','C':'G','G':'C','N':'N',
          'a':'t','t':'a','c':'g','g':'c','n':'n'}

def recomp(seq):
    return(''.join(reversed(list(map(lambda x:comp_dict.get(x,x),seq)))))

#加入barcode
bases = ['A', 'T', 'C', 'G']
def hamming_distance(s1, s2):
    # 确保字符串长度相等
    if len(s1) != len(s2):
        raise ValueError("字符串长度不相等")
    # 计算不同字符的数量
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def is_valid_barcode(barcode,ref_seqs):
    return(all(hamming_distance(barcode,i.upper())>=2 for i in ref_seqs))

def generate_barcode(length, ref_seqs):
    barcode = ''.join(random.choice(bases) for _ in range(length))
    ind = 1
    while not is_valid_barcode(barcode, ref_seqs):
        barcode = ''.join(random.choice(bases) for _ in range(length))
        ind += 1
        if ind >= 1000:
            raise ValueError("too strict to generate barcode!")
    return(barcode)

barcodes=[]
barcodes_f=[]
barcodes_r=[]
data = pd.read_csv('ttnAGGACAGAGGGTCAGCATGCCA.calitas_searchreference.out.filtered.tsv',sep="\t")
seqs = []
for i in range(len(data)):
	if data["strand"][i]=="+":
		seqs.append(extract_genome_seq(data["chromosome"][i],data["coordinate_start"][i]-3-9,data["coordinate_start"][i]+30))
	elif data["strand"][i]=="-":
		seqs.append(recomp(extract_genome_seq(data["chromosome"][i],data["coordinate_end"][i]-30,data["coordinate_end"][i]+3+9)))
	else:
		raise ValueError(str(data["strand"][i]))

new_seqs=[]
i=1
for seq in seqs:
    barcode_f = generate_barcode(14, barcodes)
    barcodes.append(barcode_f)
    barcodes_f.append(barcode_f)
    barcode_r = generate_barcode(14, barcodes)
    barcodes.append(barcode_r)
    barcodes_r.append(barcode_r)
    new_seqs.append('AGACGTTCTCACAGCAATTCGTACAGTCGACGTCGATTCGTGT'+barcode_f+'TTGACATTCTGCAATTA'+seq+'AGTATGTATGCTTCGCGCAGTGCGACTTCGCAGCGCATCACTTCA'+barcode_r+'AGAGCTGCGAGTCTTACAGCATTGCA')
    print(i)
    i+=1
data['barcode_f']=barcodes_f
data['barcode_r']=barcodes_r
data['synthesize']=new_seqs

data.to_csv("calitas.Y7_oneseq_lib.tsv",sep="\t",index=False)
