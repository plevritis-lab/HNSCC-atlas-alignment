#!/usr/bin/env python
"""
Create AbSeq fasta and gtf files
"""

import sys
import argparse
import os

# structure of the gtf data:
# seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below."
# source - name of the program that generated this feature, or the data source (database or project name)"
# feature - feature type name, e.g. Gene, Variation, Similarity"
# start - Start position of the feature, with sequence numbering starting at 1."
# end - End position of the feature, with sequence numbering starting at 1."
# score - A floating point value.
# strand - defined as + (forward) or - (reverse).
# frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on.."
# attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature."

def generate_gtf(fasta_file):
    out_gtf = fasta_file.split('.')[0]+ '.gtf'
    with open(fasta_file) as fr, open(out_gtf,'w+') as fw:
        for line in fr:
            if line.startswith('>'):
                gene = line[1:].split('\n')[0]
            else:
                n = len(line.split('\n')[0])
                fw.write(gene + '\t' + 'AbSeq\texon\t1\t' + str(n) + '\t' + '.\t+\t0\t' + 'gene_id "'+gene+'"; '+
                         'gene_name "' + gene + '"; ' + 'transcript_id "' + gene + '"\n')
    fw.close()
    fr.close()

def generate_ref(reflist, filename):
    # read list of barcode number and generate fasta and gtf file
    bc_names = []

    # read list of barcodes supplied by
    with open(reflist) as f:
        for line in f:
            bc_names.append(line.strip('\n').strip().strip(','))
    out_fasta = filename + '.fasta'
    script_dir = os.path.dirname(__file__)
    rel_path = 'AbSeqV2Final.fasta'

    # keep track of number of barcodes found and make sure it is consistent with list provided
    num_bc = 0
    with open(os.path.join(script_dir, rel_path)) as fr, open(out_fasta,'w+') as fw:
        for line in fr:
            if line.startswith('>'):
                gene = line[1:].split('\n')[0]
                bc_id = gene.split('|')[2]

            if bc_id in bc_names:
                    num_bc += 1
                    fw.write(line)
    generate_gtf(out_fasta)
    if num_bc/2 < len(bc_names):
        print("Warning: Only {} out of {} barcodes were found in the reference file.".format(num_bc/2, len(bc_names) ))


def dumpGtf(filename, species):
    out_fasta = filename + '.gtf'
    script_dir = os.path.dirname(__file__)
    if species=='h':
        rel_path = 'BDSampleTags.gtf'
    else:
        rel_path = 'BDSampleTagsMM.gtf'
    with open(os.path.join(script_dir, rel_path)) as fr, open(out_fasta,'w+') as fw:
        for line in fr:
            fw.write(line)


def dumpBDG(filename, species):
    out_fasta = filename + '.fasta'
    script_dir = os.path.dirname(__file__)
    if species=='h':
        rel_path = 'BDSampleTags.fasta'
    else:
        rel_path = 'BDSampleTagsMM.fasta'
    with open(os.path.join(script_dir, rel_path)) as fr, open(out_fasta,'w+') as fw:
        for line in fr:
            fw.write(line)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
        parser = argparse.ArgumentParser()

    # if len(argv) > 0 and argv[0] == 'bdgenome':
    #
        parser.add_argument('--bdgenome', help='generate reference files for sample multiplexing (human)',
                            action='store_true')
        parser.add_argument('--bdgenome_mm', help='generate reference file for sample multiplexing (mus musculus)',
                            action='store_true')
        parser.add_argument('--bdabseq', help='generate reference files for BD AbSeq', action='store_true')
        parser.add_argument('--genome', help='name of output fasta and gtf files', default='BDreference')
        parser.add_argument('--fasta', help='AbSeq FASTA reference file')
        parser.add_argument('--reflist', help='list of BD AbSeq barcode IDs in a txt file, one barcode Id per line')
        args = parser.parse_args()
        if args.bdgenome:
            if args.genome=='BDreference':
                args.genome='BDSampleTags'
            dumpBDG(args.genome,'h')
            dumpGtf(args.genome,'h')
        elif args.bdgenome_mm:
            if args.genome=='BDreference':
                args.genome='BDSampleTagsMM'
            dumpBDG(args.genome, 'm')
            dumpGtf(args.genome, 'm')
        elif args.bdabseq:
            if args.genome=='BDreference':
                args.genome='BDAbSeq'
            if args.fasta is not None:
                if not os.path.isfile(args.fasta):
                    sys.stderr.write("Invalid fasta file\n")
                    sys(-1)
                elif os.path.isfile(args.fasta.split('.')[0] + '.gtf'):
                    sys.stderr.write("gtf file with the same name already exists\n")
                    exit(-1)
                generate_gtf(args.fasta)
                print("{} is generated.".format(args.fasta.split('.')[0] + '.gtf'))
            elif args.reflist is not None:
                if not os.path.isfile(args.reflist):
                    sys.stderr.write("Invalid reference list file\n")
                    exit(-1)
                else:
                    generate_ref(args.reflist, args.genome)
                    print("{} and {} are generated.".format(args.genome + '.gtf', args.genome + '.fasta' ))
            else:
                sys.stderr.write("Please enter valid --fasta or --reflist input")
                exit(-1)


if __name__ == '__main__':
    sys.exit(main(None))
