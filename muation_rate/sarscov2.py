"""
@Author:weidong wu
@Author_Email:weidongwu404@gmail.com
@Time:2022/10/11 下午2:59
"""
import os
import re
from Bio import SeqIO
from collections import OrderedDict, defaultdict
from typing import Dict, List, Tuple

from aa import nuc2aa

# from .align import NEXTCLADE_DATA
REFERENCE_SEQ = None  # loaded lazily
REFERENCE_PATH = '/media/wvdon/sdata/covid_evolution/sequence/'
# Adapted from https://github.com/nextstrain/ncov/blob/50ceffa/defaults/annotation.gff
# Note these are 1-based positions
annotation_tsv = """\
seqname	source	feature	start	end	score	strand	frame	attribute
.	.	gene	26245	26472	.	+	.	 gene_name "E"
.	.	gene	26523	27191	.	+	.	 gene_name "M"
.	.	gene	28274	29533	.	+	.	 gene_name "N"
.	.	gene	29558	29674	.	+	.	 gene_name "ORF10"
.	.	gene	28734	28955	.	+	.	 gene_name "ORF14"
.	.	gene	266	13468	.	+	.	 gene_name "ORF1a"
.	.	gene	13468	21555	.	+	.	 gene_name "ORF1b"
.	.	gene	25393	26220	.	+	.	 gene_name "ORF3a"
.	.	gene	27202	27387	.	+	.	 gene_name "ORF6"
.	.	gene	27394	27759	.	+	.	 gene_name "ORF7a"
.	.	gene	27756	27887	.	+	.	 gene_name "ORF7b"
.	.	gene	27894	28259	.	+	.	 gene_name "ORF8"
.	.	gene	28284	28577	.	+	.	 gene_name "ORF9b"
.	.	gene	21563	25384	.	+	.	 gene_name "S"
"""

annotation_tsv_SMEN = """\
seqname	source	feature	start	end	score	strand	frame	attribute
.	.	gene	26245	26472	.	+	.	 gene_name "E"
.	.	gene	26523	27191	.	+	.	 gene_name "M"
.	.	gene	28274	29533	.	+	.	 gene_name "N"
.	.	gene	21563	25384	.	+	.	 gene_name "S"
"""

def sars_position():
    genes = []
    rows = annotation_tsv.split("\n")
    header, rows = rows[0].split("\t"), rows[1:]
    for row in rows:
        if row:
            row = dict(zip(header, row.split("\t")))
            gene_name = row["attribute"].split('"')[1]
            start = int(row["start"])
            end = int(row["end"])
            genes.append(((start, end), gene_name))
    genes.sort()
    return OrderedDict((gene_name, pos) for pos, gene_name in genes)


# This maps gene name to the nucleotide position in the genome,
# as measured in the original Wuhan virus.
GENE_TO_POSITION: Dict[str, Tuple[int, int]] = sars_position()
#print(GENE_TO_POSITION)

def load_reference_sequence():
    # with open(os.path.join(REFERENCE_PATH, "EPI_ISL_402124.fasta")) as f:
    #     ref = "".join(line.strip() for line in f if not line.startswith(">"))
    # #assert len(ref) == 29891, len(ref)
    RECORD = list(SeqIO.parse(REFERENCE_PATH+"EPI_ISL_402124.fasta", "fasta"))[0]
    return RECORD.seq

"""\
seqname	source	feature	start	end	score	strand	frame	attribute
.	.	gene	26245	26472	.	+	.	 gene_name "E"
.	.	gene	26523	27191	.	+	.	 gene_name "M"
.	.	gene	28274	29533	.	+	.	 gene_name "N"
.	.	gene	21563	25384	.	+	.	 gene_name "S"
"""
def NSEM_protein(sequence):
    #global GENE_TO_POSITION
    #for GENE,(begin,end) in GENE_TO_POSITION.items():
        #print(GENE)
        #if GENE=='S':
    #26245	26472 E

    E = nuc2aa(sequence[26244:26471])
    M = nuc2aa(sequence[26522:27190])
    N = nuc2aa(sequence[28273:29532])
    S = nuc2aa(sequence[21562:25383])
    return E,M,N,S

def test():
    REFERENCE_SEQ = load_reference_sequence()
    print(NSEM_protein(REFERENCE_SEQ))
    # GENE_TO_POSITION: Dict[str, Tuple[int, int]] = sars_position()
    # E, M, N, S = '','','',''
    #
    # #GENE_list = ['E', 'M', 'N', 'S']
    # for GENE,(begin,end) in GENE_TO_POSITION.items():
    #     aa = nuc2aa(REFERENCE_SEQ[begin-1:end-1])
