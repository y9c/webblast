#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-06-29 00:18


import os
import gzip
import pysam

def read_file(file_path, record_limit):
    """Reads files depending on their format, handling gzipped files and returning sequences in FASTA format without annotations."""
    _, file_ext = os.path.splitext(file_path)
    sequences = []
    
    # Decide on the method to open the file based on its compression
    def open_file(file_path):
        if file_path.endswith('.gz'):
            return gzip.open(file_path, 'rt')
        else:
            return open(file_path, 'rt')
    
    # Processing the file based on its extension
    with open_file(file_path) as file:
        if file_ext in ['.fasta', '.fa', '.fasta.gz', '.fa.gz']:
            name, seq = None, []
            for line in file:
                if line.startswith('>'):
                    if name:
                        sequences.append(f">{name}\n{''.join(seq)}")
                        if len(sequences) >= record_limit:
                            break
                    name = line[1:].split(None, 1)[0].strip()  # Take the first part of the header
                    seq = []
                else:
                    seq.append(line.strip())
            if name and len(sequences) < record_limit:
                sequences.append(f">{name}\n{''.join(seq)}")
        elif file_ext in ['.fastq', '.fq', '.fastq.gz', '.fq.gz']:
            counter = 0
            for line in file:
                if counter % 4 == 0:  # Name line
                    name = line[1:].strip().split()[0]  # Remove '@' and take the first part
                elif counter % 4 == 1:  # Sequence line
                    sequences.append(f">{name}\n{line.strip()}")
                    if len(sequences) >= record_limit:
                        break
                counter += 1
        elif file_ext in ['.sam', '.bam']:
            with pysam.AlignmentFile(file_path, "r") as samfile:
                for aln in samfile:
                    if aln.query_sequence:
                        sequences.append(f">{aln.query_name}\n{aln.query_sequence}")
                        if len(sequences) >= record_limit:
                            break

    return sequences

