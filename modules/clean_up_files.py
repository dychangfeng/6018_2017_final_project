
"""modules Description
Copyright (c) 2017 Yun Ding <u0674686@utah.edu>
This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).
@status:  experimental
@author:  Yun Ding
@contact: u0674686@utah.edu
This module clean up the G4_bed files, find tss from GFF3 files,
and and the overlap of G4 with coding regions and TSS regions
fasta_parser is a function to parse fasta files into dict"""
import os
import gzip
def get_tss_list(file_path):
    """take a gff3 file object
    return a list of tuples with all the genes from gff3 file
"""
    with gzip.open(file_path) as f0:
        lines = f0.readlines()
        genes = []
        for line in lines:
            if not line.startswith('#'):
                line = line.split()
                if line[2] == 'gene': ## only pick out 'genes'from all the features
                    if line[6] == '+':
                        genes.append((line[0], int(line[3]), \
                        int(line[3]), '.', '.', line[6]))
                        ## append seq_id, tss_start(twice), strand
                    else:
                        genes.append((line[0], int(line[4]), \
                        int(line[4]), '.', '.', line[6]))
                        # on '-' strand, append the end, instead of start
    if len(genes) > 0:
        return sorted(genes, key=lambda x: x[0:2])
        ## sort by gene_seq, start, and end

    else: return "no genes in this file"

def create_tss_from_gff3(base_dir):
    """input: base directionary
    output: save tss as bed6 format
    """
    #base_dir='/Users/Yun/Documents/bacteria_G4/D_thermus/'
    base_dir_gff = base_dir + 'gff_files/'
    if not os.path.isdir(base_dir + 'tss_files/'):
        os.mkdir(base_dir + 'tss_files/')
    for i in os.listdir(base_dir_gff):
        if os.path.isfile(base_dir + 'tss_files' + i):
            print "file {} exist".format(i)
            continue
        else:
            gff = get_tss_list(os.path.join(base_dir_gff+i))
            ## get only genes out of gff3 files
            if not type(gff) == str:
                ## when there is no genes in gff3, read_gff3 return a string
                with open(os.path.join(base_dir + 'tss_files/', \
                i.split('.gff')[0] + '.bed'), 'w') as f0:
                    f0.write('\n'.join('{}\t{}\t{}\t{}\t{}\t{}'.format(\
                    x[0], x[1], x[2], x[3], x[4], x[5]) for x in gff))
            else:
                print "no genes in file: %s"%(i)

def trim_g4_chr(base_dir):
    """trim off the extra text to match chromsome name in
    the GFF files to use bedfiles,
    save the new file to the folder all_G4_clean"""
    #base_dir='/Users/Yun/Documents/bacteria_G4/D_thermus/'
    G4_dir = base_dir + "all_G4/"
    if not os.path.isdir(base_dir + 'all_G4_clean/'):
        os.mkdir(base_dir + 'all_G4_clean/')
    for i in os.listdir(G4_dir):
        with open(G4_dir+i, 'r') as fp:
            lines = fp.readlines()
            newlines = []
            for line in lines:
                line = line.split('\t')
                seq_name = line[0].split(' ')[0]
                newlines.append((seq_name, line[1], line[2], '.', \
                line[4], line[5]))
                ## save as bed6 format later
        if len(newlines) > 0:
            with open(base_dir+'all_G4_clean/' + i, 'w') as f0:
                ## substitude GCF with GCA to match GFF files
                f0.write('\n'.join('{}\t{}\t{}\t{}\t{}\t{}'.format(\
                x[0], x[1], x[2], x[3], x[4], x[5]) for x in newlines))
        else:
            continue

def fasta_parser(fasta_path):
    """parse a fasta file, keep all the seq name as the key
    and the sequence as the value in the seq_dict,
    return the dictionary"""
    g_dict={}
    with open(fasta_path, 'r') as f0:
        lines = f0.readlines()
        seq=''
        for line in lines:
            if line.startswith('>'):
                try:
                    g_dict[key]=seq
                except:
                    pass
                key=line.strip()[1:]
                seq=''
            else:
                seq+=line.strip()
    g_dict[key]=seq ## assign the last sequence to the key
    return g_dict
