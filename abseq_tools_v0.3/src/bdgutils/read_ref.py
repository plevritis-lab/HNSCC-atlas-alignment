"""'
read_ref.py
Read GTF/Fasta files for AbSeq counting.
"""
from collections import namedtuple
import re

def read_gtf(gtf_file):
    # read gtf file
    gtf_data = []
    with open(gtf_file,'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                gtf_data.append(line.strip().split('\t'))

    Gtf = namedtuple("GTF", "contig start end gene_id")
    for rec in gtf_data:
        yield Gtf(contig=rec[0],
                  start=int(rec[3]),
                  end=int(rec[4]),
                  gene_id=re.search('gene_id "(.*?)"', rec[8]).group(1))


def allTags(gtf_file):
    return [g.gene_id for g in read_gtf(gtf_file)]
