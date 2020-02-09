# Assembly of Genome

Separate assemblies were performed for each TrueSeq library using the DISCOVAR de novo assembler (Weisenfeld et al. 2014) with default parameters. Based on the higher number of input reads and post-assembly statistics, the 350 bp insert assembly was used as a substrate for scaffolding. 

```
DiscovarDeNovo MAX_MEM_GB=900 READS="/data/wrasse/nzgl01779/Raw/H3WYYBCXX-1779-01-26-1_L001_R1.fastq.gz,/data/wrasse/nzgl01779/Raw/H3WYYBCXX-1779-01-26-1_L001_R2.fastq.gz,/data/wrasse/nzgl01779/Raw/H3WYYBCXX-1779-01-26-1_L002_R1.fastq.gz,/data/wrasse/nzgl01779/Raw/H3WYYBCXX-1779-01-26-1_L002_R2.fastq.gz" OUT_DIR=350bp
```