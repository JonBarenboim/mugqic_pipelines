#!/usr/bin/env python

import sys
import os
from Bio import SeqIO

# if len(sys.argv) != 3:
#     raise Exception('Need an input directory and output directory')
#
# input_dir = sys.argv[1]
# output_dir = sys.argv[2]
input_dir = 'remove_duplicates'
output_dir = 'remove_rrna'

for i in (1, 2):
    id_filename = os.path.join(input_dir, 'cow{}_qual_all_unique_IDs.txt'.format(i))
    fasta_filename = os.path.join(input_dir, 'cow{}_qual_all_unique.fasta'.format(i))

    out_filename = os.path.join(output_dir, 'cow{}_IDs_length.txt'.format(i))

    reads = dict()
    with open(id_filename) as id_file:
        for line in id_file:
            id, value = line.rstrip().split('\t')
            reads[id] = value

    print('Number of unique reads: ' + str(len(reads)))

    print('Writing to ' + out_filename)
    with open(out_filename, 'w+') as out_file:
        for record in SeqIO.parse(fasta_filename, 'fasta'):
            out_file.write('\t'.join([record.id, str(reads[record.id]), str(len(record.seq))]) + '\n')
