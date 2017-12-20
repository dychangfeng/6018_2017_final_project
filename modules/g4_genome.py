"""modules Description
Copyright (c) 2017 Yun Ding <u0674686@utah.edu>
This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).
@status:  experimental
@author:  Yun Ding
@contact: u0674686@utah.edu
This modules deal with a genome and its GFF3 file,
to find the GC content, the genome length, G4 denstiy,
and popular G4s"""
from collections import Counter, defaultdict

def reverse_comp(seq):
    """take a sequence and return the reverse complementary strand of this seq
    input: a DNA sequence
    output: the reverse complementary of this seq"""
    bases_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    ##get the complementary if key is in the dict, otherwise return the base
    return ''.join(bases_dict.get(base, base) for base in seq[::-1])

class Genome(defaultdict):
    """genome class, each genome has one G4 class"""
    def __init__(self, ID, *args, **kwargs):
        super(Genome, self).__init__(*args, **kwargs)
        self.ID = ID
    @property
    def seq(self):
        """add all the sequence together"""
        seq = ''
        for val in self.itervalues():
            seq += val
        return seq
    @property
    def gc(self):
        """calcuate GC percentage of a genome"""
        try:
            c_seq = Counter(self.seq.upper())
            return "%.2f"%(float((c_seq['G']+c_seq['C']))/sum(c_seq.viewvalues()))
        except ZeroDivisionError:
            return "Not provided"
    @property
    def genome_len(self):
        """calculate genome length in million base pairs"""
        return "%.2f"%(float(len(self.seq))/1000000)

    def __str__(self):
        return "{} has a genome of {} Mb, the GC content is {}".format(\
    self.ID, self.genome_len, self.gc)



class G4_genome(Genome):
    """G4 info in a specific genome
    G4_path: file path of G4_bed file"""
    """G4 info in a specific genome
    G4_path: file path of G4_bed file"""
    def __init__(self, G4_path, *args, **kwargs):
        self.G4_path = G4_path
        super(G4_genome, self).__init__(*args, **kwargs)
    @property
    def G4_density(self):
        """return the number of G4 per Mb"""
        try:
            with open(self.G4_path, 'r') as f0:
                g4_line = f0.readlines()
            return "%d"%(len(g4_line)/float(self.genome_len))
        except ZeroDivisionError:
            return "please provide the genome first"
    @property
    def G4_list(self):
        """keep all G4 in the genome in one list"""
        with open(self.G4_path, 'r') as f0:
            lines = f0.readlines()
            g4_ls = []
            for line in lines:
                line = line.strip()
                line = line.split('\t')
                if line[6][0] == 'G':
                    g4_ls.append(line[6])
                else:
                    g4_ls.append(reverse_comp(line[6]))
        return g4_ls
    def most_freq_G4(self, n=10):
        """count the most frequent G4 from G4_list, defaut number is top 10 """
        count_g4 = Counter(self.G4_list)
        return count_g4.most_common(n)

    def __str__(self):
        try:
            float(self.G4_density)
            return "The G4_density of {} is {} per Mb".format(self.ID, self.G4_density)
        except ValueError:
            return "The {} has {} G4s".format(self.ID, len(self.G4_list))
