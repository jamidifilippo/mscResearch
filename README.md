# mscResearch

Rscripts and python notebooks for my master's research project at King's College London.

The organisms used in this study were Hydra vulgaris, Drosophila melanogaster, Homo sapiens.
[] = hydra, dmel, hsap

Promoter and protein sequence fasta files were generated by running
  []_promoter.R <in> <feature> <out> ; where in was a gff file, feature was what defined a gene
  []_translation.R <in> <out> ; where in was a protein sequence file in fasta format
  
The "exact match" pattern search was run with 
  prediction/patternSearch.R <in> <out>
  
Analysis of the generated files was performed in python jupyter notebooks listed under the results directory.



Input files can be downloaded from:
Hydra

DNA: https://research.nhgri.nih.gov/hydra/download/assembly/Hm105_Dovetail_Assembly_1.0.fa.gz

protein: https://research.nhgri.nih.gov/hydra/download/genemodels_proteins/hydra2.0_genemodels.aa.gz

gff: https://research.nhgri.nih.gov/hydra/download/genemodels_gff3/hydra2.0_genemodels.gff3.gz

Dmel

DNA: ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.29.fasta.gz

protein: ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-translation-r6.29.fasta.gz



Hsap

DNA: ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

protein: ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.abinitio.fa.gz

gff: ftp://ftp.ensembl.org/pub/release-96/gff3/homo_sapiens/Homo_sapiens.GRCh38.96.abinitio.gff3.gz
