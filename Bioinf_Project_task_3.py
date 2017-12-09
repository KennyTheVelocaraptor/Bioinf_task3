from Bio import pairwise2
from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.Align import MultipleSeqAlignment


seq1 = SeqIO.read("proteins/protein_c.faa", "fasta")
seq2 = SeqIO.read("proteins/protein_z.faa", "fasta")
alignments_glob = pairwise2.align.globalds(seq1.seq, seq2.seq, blosum62, -8, -0.3)
print ("#####Best Global aligments:")
for a in alignments_glob:
	print(pairwise2.format_alignment(*a))
alignments_loc = pairwise2.align.localms(seq1.seq, seq2.seq, 15, -12, -8, -0.3)
print ("#####Best Local aligments:")
for a in alignments_loc:
	print(pairwise2.format_alignment(*a))

