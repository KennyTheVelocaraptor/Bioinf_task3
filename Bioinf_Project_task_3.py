from Bio import pairwise2
from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
from Bio.Emboss.Applications import NeedleCommandline
import sys


rna = "RNAs/schmallenberg_nucleprotein_RNA_complex"

def pairwise():
	## simple sequences` aligning variables
	seq1 = SeqIO.read("proteins/protein_c.faa", "fasta")
	seq2 = SeqIO.read("proteins/protein_z.faa", "fasta")

	## Calculating best global alignments using blosum62 matrix of matches
	alignments_glob = pairwise2.align.globalds(seq1.seq, seq2.seq, blosum62, -8, -0.3)
	print ("#####Best Global aligments:")
	for a in alignments_glob:
		print(pairwise2.format_alignment(*a))
	## Calculating best local alignments using extended penalties system for match/not match
	alignments_loc = pairwise2.align.localms(seq1.seq, seq2.seq, 15, -12, -8, -0.3)
	print ("#####Best Local aligments:")
	for a in alignments_loc:
		print(pairwise2.format_alignment(*a))



def clustal():
	##RNA clustalw
	cline = ClustalwCommandline("clustalw", infile=rna + ".fasta")
	cline()
	align = AlignIO.read(rna + ".aln", "clustal")
	print (align)	
	tree = Phylo.read(rna +	".dnd", "newick")
	Phylo.draw_ascii(tree)


def muscle():
	##Muscle

	muscle_cline = MuscleCommandline(input=rna + ".fasta")
	stdout, stderr = muscle_cline()
	muscle_align = AlignIO.read(StringIO(stdout), "fasta")
	print(muscle_align)

def mboss():
	##Water and Needle
	needle_cline = NeedleCommandline(asequence=rna + ".fasta", bsequence="proteins/protein_c.faa", gapopen=10, gapextend=0.5, outfile="water_needle.txt")
	needle_cline()
	needle_align = AlignIO.read("water_needle.txt", "emboss")
	print(needle_align)

switcher = {1 : pairwise,
           2 : clustal,
           3 : muscle,
           4 : mboss,
}

switcher[int(sys.argv[1])]()
##For RNA coded sequences 
## Variables
# rna1 = open("RNAs/ebola_glycoprotein_RNA.pdb", "rU")
# rna2 = open("RNAs/schmallenberg_nucleprotein_RNA_complex.pdb", "rU")
# rna1_parsed = list(SeqIO.parse(rna1, "pdb-seqres"))
# rna2_parsed = list(SeqIO.parse(rna2, "pdb-seqres"))
# #print (rna2_parsed)
# #for record in rna2_parsed:
# #    print (record.description)
# #Function, which globally aligns all sequences inside the RNA using blosum62 matching matrix and taking file_name,gap_open_penalty and gap_extension_penalty as input arguments
# #=def AlignGlobalRNA(seq_list,gap_open_penalty,gap_extension_penalty):
# 	ultimate_record = seq_list[0]
# 	count = len(seq_list)
# 	for record in seq_list[:count-1:]:
# 		alignments_glob_in_rna = pairwise2.align.globalds(ultimate_record, record, blosum62, gap_open_penalty, gap_extension_penalty)
# 		ultimate_record = alignments_glob_in_rna[0]
# 	return (pairwise2.format_alignment(*alignments_glob_in_rna[0]))	

# #Function, which locally aligns all sequences inside the RNA using extended penalties system and taking file_name,match_penalty,not_match_penalty,gap_open_penalty,gap_extension_penalty as input arguments
# def AlignLocalRNA(seq_list,match_penalty,not_match_penalty,gap_open_penalty,gap_extension_penalty):
# 	ultimate_record = seq_list[0]
# 	count = len(seq_list)
# 	for record in seq_list[:count-1:]:
# 		alignments_loc_in_rna = pairwise2.align.localms(ultimate_record, record, match_penalty, not_match_penalty, gap_open_penalty, gap_extension_penalty)
# 		ultimate_record = alignments_loc_in_rna[0]
# 	return (pairwise2.format_alignment(*alignments_loc_in_rna[0]))	

# #Function to exclude sequences of different types from list, cause only sequences of same type can be aligned
# def SortSeq(sequences_list,seq_description):
# 	sorted_list = list()
# 	for i in sequences_list:
		
# 		if seq_description == i.description[:3]:
# 			print (i.description[:3] + i.seq)
# 			sorted_list.append("".join(i.seq))
# 	return (sorted_list)


# #na1_glob_aligment = DecodeGlobalRNA(rna2,-8,-0.3)
# test = SortSeq(rna2_parsed,"UNP")

# #rna2_loc_aligment = AlignGlobalRNA(test,-8,-0.3)
# alignments_loc_in_rna = pairwise2.align.globalds(test[0],test[1], blosum62, -8, -0.3)
# alignments_loc_in_rna2 = pairwise2.align.globalds(alignments_loc_in_rna,test[2], blosum62, -8, -0.3)
# print (alignments_loc_in_rna2[0])



