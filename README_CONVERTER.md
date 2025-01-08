# bam2csv converter 
The bam2csv converter is intended to convert bam files into table format (CSV) as part of the next-generation sequencing (NGS) pipeline.
# Introduction
The NGS is a method that is widely used for basic research tasks as well as for translational projects. It is also actively used in industry as well as for medical diagnosis and outcome prediction. 
NGS produces raw data that should be aligned to reference sequences to find similarities/differences, mutations, etc. Usually, different program aligners are used depending on the task, preferences, etc. But most of them (if not all) create files in bam format that can be used for downstream analysis.
# The proposal
Development of a tool that converts the bam files into a data table containing information about:
  1. The coordinates of nucleotides, 
  2. Genes present in each position, 
  3. Reference sequence, 
  4. Total coverage (in each position), 
  5. Separate coverage for each of nucleotides (A, G, T, C) in each position
  6. Some other information (e.g. percentage of each nucleotide). 
  Together with the information about reference sequence, this information will be a starting point for downstream analysis - various genes and changes can be analyzed (substitution, deletions, insertions), data can be used for graph creation, etc.

# Technical implementation
The bam2csv converter will be implemented based on Python script. After importing the necessary information to extract data from the bam file, it will create data tables for all the bam files in the folder and will create a CSV file, that can be processed further. 
Input parameters will include:
input folder (containing .bam and .bai (indexed .bam ) files)
output folder
fasta file (with reference sequence)
bed file (with data about genes)
The output file will contain columns with the following:
  1. The coordinates of nucleotides, 
  2. Genes present in each position, 
  3. Reference sequence, 
  4. Total coverage, 
  5. Separate coverage for each of nucleotides (A, G, T, C) in each position
  6. Some other information (e.g. percentage of each nucleotide). 
# Conclusion 
The bam2csv converter will be a convenient intermediate data file that can be used for further processing of data (depending on scientist demands) to analyze or graphically present data/analysis results.

# How to run bam2csv converter:
python "path to bam2csv.py" --input_folder "path to input folder" --output_folder "path to output folder" --reference_fasta "path to fasta file"  --bed_file "path to bed file"
