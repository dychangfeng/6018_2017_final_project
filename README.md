# intro_to_programming_with_python_FINAL_PROJECT_YUN

Youtube video that describe this project:
https://youtu.be/Uf4yLfONwjo

This package is desinged to look for G-quadruplexes (G4s) in a genome or a list of genomes in terms where are they and whether they are enriched in the gene regulation region.

Three types of files are needed to use this package:
1. The genome of interest in fasta format. All genomes are in the genomes folder under the root directory.
2. The genome annotation file in GFF3 format. All GFF3 are in the gff_files folder under the root directory.
3. All G4s found in the genome by Quadparser in bed format. All G4 bed files are in the all_G4 foler under the root directory.

Environment and modules required for this package:
1. This package is developed under PYTHON 2.7. But it should also work with PYTHON 3.7 with minor changes.
2. All the package required: pandas, numpy, OS, gzip, pybedtools, re, matplotlib, holoviews, jupyter, and notebook.

Modules in this package:
There are three modules in this package:
clean_up_files.py has functions to prepare and clean up files required in this package.
gff_g4_functions.py has most of the functions to look for G4s in different genomic elements.
g4_genome.py has two classes designed to deal with genome file and G4 bed file to find the most popular G4s and G4 density in the genome of interest.

Instructions:
the D_thermus notebook contains a step by step instructions to use this package.
base_dir which is the directory to the root directory (where all the files are stored) is needed for this module.

Due to the sizes of all the files, original files required for the D_thermus notebook are not uploaded. But all the results are uploaded to the results folder.

Questions asked in this package:

1. Does different genus of bacteria have different G4 density (number of G4s per million base pairs)?
    This question can be answered by loading both genome sequence (fasta file) and list of G4s in the genome in to G4_genome class. The G4_genome class has properties such as GC_percentage, G4_density, and most popular G4s.
    
    
2. Is G4s enriched in any genomic element such as around TSS (transcription start site) region?
    This question was approached by different functions in gff_G4_functions to calculate the percentage of G4s in the coding regions and around TSS regions. 


