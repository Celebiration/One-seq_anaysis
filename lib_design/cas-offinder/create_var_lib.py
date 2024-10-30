#!/usr/local/anaconda/bin/python
import sys
import pandas as pd
import numpy as np
import subprocess
import os
import multiprocessing
from Bio import SeqIO

max_processes = 50
seq1="NNNNNNNNNNNNNNNNNNNNNGG"
seq2="CCCGCACCTTGGCGCAGCGGNNN"

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


#merged=merge_hits([file1,file2])
#merged.to_csv("merged_6.tsv",sep="\t",index=False)

#读取
merged=pd.read_csv("merged_6.tsv",sep="\t")

#tool functions
def extract_var(chr,start,end,type="snps"):
	vcf="ALL."+chr+".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
	if not os.path.exists("/home/fengchr/staff/zhk_11_24_one_seq/new/"+vcf):
		return([])
	snps=[]
	command=["bcftools", "view", "-r", chr[3:]+":"+str(start)+"-"+str(end),"/home/fengchr/staff/zhk_11_24_one_seq/new/"+vcf, "-v", type, "-H"]
	process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
	for line in process.stdout:
		l=line.strip().split("\t")
		snps.append((int(l[1]),l[3],l[4]))
	process.wait()
	return(snps)

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

# def extract_genome_seq(chr,start,end):
# 	command=["seqkit","grep","-p",chr,"/Yol_Data/resources/Gencode_human_ref/GRCh38.p13.genome.fa","|","seqkit","subseq","--quiet","-r",str(start)+":"+str(end),"-w","0","|","tail","-1"]
# 	return(str(subprocess.check_output(" ".join(command),shell=True).decode('utf-8').strip()))

def extract_genome_seq(chr,start,end):
	return(GRCh38.get(chr)[start - 1:end])

comp_dict={'A':'T','T':'A','C':'G','G':'C','N':'N',
          'a':'t','t':'a','c':'g','g':'c','n':'n'}
def recomp(seq):
    return(''.join(reversed(list(map(lambda x:comp_dict.get(x,x),seq)))))

def process_one(lock,var_lib,i):
	strand=merged["Direction"][i]
	chr=merged["Chromosome"][i].split(" ")[0]
	start=merged["Position"][i]+1
	target=merged["DNA"][i]
	target_len=len(target)-merged["DNA"][i].count("-")
	target_indices=[k for k in range(len(target)) if target[k]!="-"]
	no_gap_target=target.replace("-","")
	end=start+target_len-1

	snps=extract_var(chr,start,end,"snps")
	indels=extract_var(chr,start,end,"indels")

	#original
	if strand=="+":
		core_seq=extract_genome_seq(chr,start-10,start+32)
	else:
		core_seq=recomp(extract_genome_seq(chr,end-32,end+10))
	if core_seq[10:(10+target_len)].upper()!=no_gap_target.upper():
		raise ValueError("%s第%d位不为%s！" %(core_seq, 11, no_gap_target))
	flank_l=core_seq[:10]
	flank_r=core_seq[(10+target_len):]
	var_lib.append(("ori","",merged["#Bulge type"][i],merged["crRNA"][i],target,target,core_seq,chr,start,end,strand,merged["Mismatches"][i],merged["Bulge Size"][i]))

	for snp in snps:
		if strand=="+":
			relative_target_pos=target_indices[snp[0]-start]
			if target[relative_target_pos].upper()!=snp[1].upper():
				raise ValueError("%s第%d位不为%s！" %(target, relative_target_pos+1, snp[1]))
			new_target=target[:relative_target_pos]+snp[2]+target[(relative_target_pos+1):]
		else:
			relative_target_pos=target_indices[target_len-(snp[0]-start)-1]
			if target[relative_target_pos].upper()!=recomp(snp[1].upper()):
				raise ValueError("%s第%d位不为%s！" %(target, relative_target_pos+1, recomp(snp[1])))
			new_target=target[:relative_target_pos]+recomp(snp[2])+target[(relative_target_pos+1):]
		var_lib.append(("snp",str(snp),merged["#Bulge type"][i],merged["crRNA"][i],target,new_target,flank_l+new_target.replace("-","")+flank_r,chr,start,end,strand,merged["Mismatches"][i],merged["Bulge Size"][i]))
	for indel in indels:#直接纳入，固定PAM的位置
		old_region=extract_genome_seq(chr,min(indel[0],start)-25,max(indel[0]+len(indel[1]),end)+25)
		relative_old_region_indel_pos=indel[0]-(min(indel[0],start)-25)+1
		if (old_region[(relative_old_region_indel_pos-1):(relative_old_region_indel_pos-1+len(indel[1]))]).upper()!=indel[1].upper():
			raise ValueError("%s第%d位处的indel突变ref不为%s！" %(old_region, relative_old_region_indel_pos, indel[1]))
		new_region=old_region[:(relative_old_region_indel_pos-1)]+indel[2]+old_region[(relative_old_region_indel_pos+len(indel[1])-1):]
		if strand=="+":
			with lock:
				with open("tmp.fa","w") as f:
					f.write(">tmp\n"+new_region)
				with open("tmp.config","w") as f:
					f.write("tmp.fa\n"+seq1+" 2 2\n"+seq2+" 5")
				with open("tmp.config2","w") as f:
					f.write("tmp.fa\n"+seq1+" 0 0\n"+seq2+" 7")
				subprocess.run(["cas-offinder-bulge","tmp.config","G","tmp.out1"])
				subprocess.run(["cas-offinder-bulge","tmp.config2","G","tmp.out2"])
				tmp=merge_hits(["tmp.out1","tmp.out2"]).iloc[0:1]
			if len(tmp)==0:
				continue
			try:
				tmp_strand=tmp["Direction"][0]
			except:
				print(tmp)
				print(str(len(tmp)))
			tmp_start=tmp["Position"][0]+1
			tmp_target=tmp["DNA"][0]
			tmp_target_len=len(tmp_target)-tmp_target.count("-")
			tmp_target_indices=[k for k in range(len(tmp_target)) if tmp_target[k]!="-"]
			tmp_no_gap_target=tmp_target.replace("-","")
			tmp_end=tmp_start+tmp_target_len-1
			if tmp_strand!="+":
				continue
			if tmp_start<11 or tmp_start+32 > len(new_region):
				continue
			core_seq=new_region[(tmp_start-11):(tmp_start+32)]
			# if tmp_end<33 or tmp_end+10 > len(new_region):
			# 	continue
			# core_seq=recomp(new_region[(tmp_end-33):(tmp_end+10)])
			if core_seq[10:(10+tmp_target_len)].upper()!=tmp_no_gap_target.upper():
				raise ValueError("%s第%d位不为%s！" %(core_seq, 11, tmp_no_gap_target))
			var_lib.append(("indel",str(indel),tmp["#Bulge type"][0],tmp["crRNA"][0],target,tmp_target,core_seq,chr,start,end,strand,tmp["Mismatches"][0],tmp["Bulge Size"][0]))
		else:
			new_region=recomp(new_region)
			with lock:
				with open("tmp.fa","w") as f:
					f.write(">tmp\n"+new_region)
				with open("tmp.config","w") as f:
					f.write("tmp.fa\n"+seq1+" 2 2\n"+seq2+" 5")
				with open("tmp.config2","w") as f:
					f.write("tmp.fa\n"+seq1+" 0 0\n"+seq2+" 7")
				subprocess.run(["cas-offinder-bulge","tmp.config","G","tmp.out1"])
				subprocess.run(["cas-offinder-bulge","tmp.config2","G","tmp.out2"])
				tmp=merge_hits(["tmp.out1","tmp.out2"]).iloc[0:1]
			if len(tmp)==0:
				continue
			tmp_strand=tmp["Direction"][0]
			tmp_start=tmp["Position"][0]+1
			tmp_target=tmp["DNA"][0]
			tmp_target_len=len(tmp_target)-tmp_target.count("-")
			tmp_target_indices=[k for k in range(len(tmp_target)) if tmp_target[k]!="-"]
			tmp_no_gap_target=tmp_target.replace("-","")
			tmp_end=tmp_start+tmp_target_len-1
			if tmp_strand!="+":
				continue
			if tmp_start<11 or tmp_start+32 > len(new_region):
				continue
			core_seq=new_region[(tmp_start-11):(tmp_start+32)]
			# if tmp_end<33 or tmp_end+10 > len(new_region):
			# 	continue
			# core_seq=recomp(new_region[(tmp_end-33):(tmp_end+10)])
			if core_seq[10:(10+tmp_target_len)].upper()!=tmp_no_gap_target.upper():
				raise ValueError("%s第%d位不为%s！" %(core_seq, 11, tmp_no_gap_target))
			var_lib.append(("indel",str(indel),tmp["#Bulge type"][0],tmp["crRNA"][0],target,tmp_target,core_seq,chr,start,end,strand,tmp["Mismatches"][0],tmp["Bulge Size"][0]))

#生成文库
with multiprocessing.Manager() as manager:
	var_lib = manager.list()  # 创建共享列表
	locker=manager.Lock()
	with multiprocessing.Pool(max_processes) as pool:
		pool.starmap(process_one, [(locker,var_lib, i) for i in range(len(merged))])
	varlib = list(var_lib)

result=pd.DataFrame(varlib,columns=["type","variant","Bulge_type","crRNA","original_target","target","core_seq","chr","start","end","strand","Mismatches","Bulge_Size"])
result.to_csv("varlib.tsv",sep="\t",index=False)


##读取文库
result=pd.read_csv("varlib.tsv",sep="\t")
#sort
chr_mapping={"chr1":1,"chr2":2,"chr3":3,"chr4":4,"chr5":5,"chr6":6,"chr7":7,"chr8":8,"chr9":9,"chr10":10,"chr11":11,"chr12":12,"chr13":13,"chr14":14,"chr15":15,"chr16":16,"chr17":17,"chr18":18,"chr19":19,"chr20":20,"chr21":21,"chr22":22,"chrX":23,"chrY":24}
def custom_sort(value):
	return chr_mapping.get(value, float('inf'))

result['custom_sorted'] = result['chr'].map(custom_sort)
result.sort_values(by=['custom_sorted',"chr","start"], ascending=True,inplace=True)
result.drop(columns=["custom_sorted"],inplace=True)
result.to_csv("varlib.sorted.tsv",sep="\t",index=False)
#删除core_seq重复的行
result['core_seq'] = result['core_seq'].str.upper()
result.drop_duplicates(subset=['core_seq'], inplace=True)
result.to_csv("varlib.sorted.core_seq_uniq.tsv",sep="\t",index=False)
