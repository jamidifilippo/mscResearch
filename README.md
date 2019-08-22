# mscResearch

Rscripts and python notebooks for my master's research project at King's College London.

The organisms used in this study were _Hydra vulgaris, Drosophila melanogaster, Homo sapiens.
[] = hydra, dmel, hsap

Promoter and protein sequence fasta files were generated by running
  []_promoter.R <in> <feature> <out> ; where in was a gff file, feature was what defined a gene
  []_translation.R <in> <out> ; where in was a protein sequence file in fasta format
  
The "exact match" pattern search was run with 
  prediction/patternSearch.R <in> <out>
  
Analysis of the generated files was performed in python jupyter notebooks listed under the results directory.
