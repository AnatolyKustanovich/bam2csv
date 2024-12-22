## bam2csv converter 
Program for conversion of bam files to table format (csv) after alignment
# Introduction
The next-generation sequencing (NGS) is a method that widely used for the basic research tasks as well as for translational projects. It is also actively used in industry and for medical diagnosis and outcome prediction. NGS produce raw data that should be aligned to reference sequence to find similarities/differences, mutations, etc. Usually, different programs-aligners are used depending on the task, preferences, etc. But most of them (if not all) create files in bam format that can be used for downstream analysis.
# The proposal
Development a tool that converts the bam file into a data table with information about:
  the coordinates of nucleotides, 
  genes present in each position, 
  reference sequence, 
  total coverage, 
  separate coverage for each of nucleotides (A, G, T, C) in each position
  +some other information (e.g. percentage of each nucleotide). 
  Compared with the information about reference sequence, this information will be a starting point for downstream analysis - various genes and changes can be analyzed (substitution, deletions, insertions), data can be used for graph creation, etc.
  
