import os
import sys
import numpy
import re
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import numpy as np
import pandas as pd

fq = sys.argv[1]
barcode_mapping_file = sys.argv[2]#barcode, id, targetseq
output_file = sys.argv[3]
constant1 = 'GCTGACTAGACACTGCTATCACACTCTCTCA'
constant2 = 'AGACGTTCTCACAGCAATTCGTACAGTCGACGTCGATTCGTGT'
constant3 = 'TTGACATTCTGCAATTA'
constant4 = 'AGTATGTATGCTTCGCGCAGTGCGACTTCGCAGCGCATCACTTCA'
constant5 = 'AGAGCTGCGAGTCTTACAGCATTGCA'
constant6 = 'TGTCTCATTTACGAAGCATCAGTGATTACGC'
fwd_cut_site = 27
fwd_tolerate = 2
rev_cut_site = 17
rev_tolerate = 2
barcode_len = 14
umi_len = 11

mapping = pd.read_csv(barcode_mapping_file, header=None)
barcode_to_id = mapping.set_index(0)[1].to_dict()
barcode_to_target = mapping.set_index(0)[2].to_dict()

#functions
comp_dict={'A':'T','T':'A','C':'G','G':'C','N':'N',
		  'a':'t','t':'a','c':'g','g':'c','n':'n'}

def recomp(seq):
	return(''.join(reversed(list(map(lambda x:comp_dict.get(x,x),seq)))))

aligner = PairwiseAligner(mode="local", match_score=1.0, mismatch_score=-1.0, open_gap_score=-2, extend_gap_score=-1)

def seq_is_in_ref(seq, ref, max_penalty = 5):
	if seq in ref:
		return (True, next(re.finditer(seq,ref)).span())
	alignment = aligner.align(ref, seq)[0]
	unaligned_len = len(seq) - (alignment.coordinates[1][-1] - alignment.coordinates[1][0])
	if len(seq) - alignment.score + unaligned_len > 5:
		return (False, (None, None))
	else:
		return(True, (alignment.coordinates[0][0], alignment.coordinates[0][-1]))
	
class seq_obj:
	def __init__(self, sequence):
		self.seq = sequence.upper()
		self.cate = None
		self.constant1_ = None
		self.constant2_ = None
		self.constant3_ = None
		self.constant4_ = None
		self.constant5_ = None
		self.constant6_ = None
		self.umi_ = None
		self.barcode_ = None
		self.target_ = None
		self.target = None
		self.targetid = None
		self.target_coverage = (None, None)
		self.target_matched = None


res = []
for record in SeqIO.parse(fq, "fastq"):
	seq = str(record.seq)
	revseq = recomp(seq)
	constant1_ = seq_is_in_ref(constant1, seq)
	constant2_ = seq_is_in_ref(constant2, seq)
	constant3_ = seq_is_in_ref(constant3, seq)
	constant4_ = seq_is_in_ref(constant4, seq)
	constant5_ = seq_is_in_ref(constant5, seq)
	constant6_ = seq_is_in_ref(constant6, seq)

	rev_constant1_ = seq_is_in_ref(constant1, revseq)
	rev_constant2_ = seq_is_in_ref(constant2, revseq)
	rev_constant3_ = seq_is_in_ref(constant3, revseq)
	rev_constant4_ = seq_is_in_ref(constant4, revseq)
	rev_constant5_ = seq_is_in_ref(constant5, revseq)
	rev_constant6_ = seq_is_in_ref(constant6, revseq)

	if sum([i[0] for i in [rev_constant1_, rev_constant2_, rev_constant3_, rev_constant4_, rev_constant5_, rev_constant6_]]) > sum([i[0] for i in [constant1_, constant2_, constant3_, constant4_, constant5_, constant6_]]):
		#说明是反向互补的
		tmp = seq_obj(revseq)
		tmp.constant1_ = rev_constant1_
		tmp.constant2_ = rev_constant2_
		tmp.constant3_ = rev_constant3_
		tmp.constant4_ = rev_constant4_
		tmp.constant5_ = rev_constant5_
		tmp.constant6_ = rev_constant6_
	else:
		tmp = seq_obj(seq)
		tmp.constant1_ = constant1_
		tmp.constant2_ = constant2_
		tmp.constant3_ = constant3_
		tmp.constant4_ = constant4_
		tmp.constant5_ = constant5_
		tmp.constant6_ = constant6_
	if sum([i[0] for i in [tmp.constant1_, tmp.constant2_, tmp.constant3_, tmp.constant4_, tmp.constant5_, tmp.constant6_]]) == 0:
		res.append((tmp.cate, tmp.constant1_[0], tmp.constant2_[0], tmp.constant3_[0], tmp.constant4_[0], tmp.constant5_[0], tmp.constant6_[0], tmp.umi_, tmp.barcode_, tmp.target_, tmp.targetid, tmp.target, tmp.target_coverage, tmp.target_matched))
		continue
	if sum([i[0] for i in [tmp.constant1_, tmp.constant2_, tmp.constant3_]]) > 0 and sum([i[0] for i in [tmp.constant4_, tmp.constant5_, tmp.constant6_]]) > 0:
		tmp.cate = "uncut"
		res.append((tmp.cate, tmp.constant1_[0], tmp.constant2_[0], tmp.constant3_[0], tmp.constant4_[0], tmp.constant5_[0], tmp.constant6_[0], tmp.umi_, tmp.barcode_, tmp.target_, tmp.targetid, tmp.target, tmp.target_coverage, tmp.target_matched))
		continue
	if sum([i[0] for i in [tmp.constant1_, tmp.constant2_, tmp.constant3_]]) > sum([i[0] for i in [tmp.constant4_, tmp.constant5_, tmp.constant6_]]):
		#说明是proto side而非PAM side
		tmp.cate = "proto_side"
		if tmp.constant1_[0] and tmp.constant2_[0]:
			tmp.umi_ = tmp.seq[tmp.constant1_[1][1]:tmp.constant2_[1][0]]
		else:
			tmp.umi_ = None
		if tmp.constant2_[0] and tmp.constant3_[0]:
			tmp.barcode_ = tmp.seq[tmp.constant2_[1][1]:tmp.constant3_[1][0]]
		else:
			tmp.barcode_ = None
		if tmp.constant3_[0]:
			tmp.target_ = tmp.seq[tmp.constant3_[1][1]:]
		else:
			tmp.target_ = None
	else:
		tmp.cate = "PAM_side"
		if tmp.constant5_[0] and tmp.constant6_[0]:
			tmp.umi_ = tmp.seq[tmp.constant5_[1][1]:tmp.constant6_[1][0]]
		else:
			tmp.umi_ = None
		if tmp.constant4_[0] and tmp.constant5_[0]:
			tmp.barcode_ = tmp.seq[tmp.constant4_[1][1]:tmp.constant5_[1][0]]
		else:
			tmp.barcode_ = None
		if tmp.constant4_[0]:
			tmp.target_ = tmp.seq[:tmp.constant4_[1][0]]
		else:
			tmp.target_ = None
	#匹配barcode
	if tmp.barcode_ in barcode_to_id.keys():
		tmp.targetid = barcode_to_id[tmp.barcode_]
		tmp.target = barcode_to_target[tmp.barcode_]
	#比对target
	if tmp.target_ and tmp.targetid:
		tmpp = seq_is_in_ref(tmp.target_, tmp.target, 5)
		if tmpp[0]:
			tmp.target_coverage = tmpp[1]
			#检查比对坐标是否match
			if tmp.cate == "proto_side":
				if tmp.target_coverage[0] == 0 and abs(tmp.target_coverage[1] - fwd_cut_site) <= fwd_tolerate:
					tmp.target_matched = True
				else:
					tmp.target_matched = False
			else:
				if tmp.target_coverage[1] == len(tmp.target) and abs(tmp.target_coverage[0] - rev_cut_site) <= rev_tolerate:
					tmp.target_matched = True
				else:
					tmp.target_matched = False
	#检查是否有A-to-G
	#...
	res.append((tmp.cate, tmp.constant1_[0], tmp.constant2_[0], tmp.constant3_[0], tmp.constant4_[0], tmp.constant5_[0], tmp.constant6_[0], tmp.umi_, tmp.barcode_, tmp.target_, tmp.targetid, tmp.target, tmp.target_coverage, tmp.target_matched))
	
res = pd.DataFrame(res, columns=["type","constant1","constant2","constant3","constant4","constant5","constant6","umi","barcode","target_seq","targetid","target","target_coverage","target_matched"])
res.to_csv(output_file, index=False)
