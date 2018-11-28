# Annotation of Bluehead genome

Annotation was performed following the Trinotate version 3.1.1 and PASA version 2.2 pipelines. Other scripts were then used to modify gff3 files for viewing in SeqMonk.

## Mapping transcripts to the genome

Following PASA (Program to Assemble Spliced Alignments) guidelines, the transcripts from the published transcriptome were mapped to the genome using both GMAP and blat. Due to the size of the transcriptome, mapping was done separately from the PASA pipeline, on the New Zealand national cluster (NeSI).

`blat -t=dna -q=dna bluehead_genome_v1_27_02_2018.fa ann-trinity-beta-jaccard_clip-v3-SS-RF-2015-01-09-fix.fasta blat_mapA.pl`

```
$GMAP_DIR/gmap -d blue_genomeV1 \
  -D ../genomes/gmapDB --nthreads=10 -B 5 \ 
  -f 3 -n 0 -x 50 --intronlength=100000 \
  ann-trinity-beta-jaccard_clip-v3-SS-RF-2015-01-09-fix.fasta \
   > blue_transA_to_genomeV1_pasa.gff3
```

## PASA: Assembling transcripts to genome

The mapping results were incorporated into the PASA pipeline using the *Launch_PASA_pipeline.pl* script. As per [guidelines](https://github.com/PASApipeline/PASApipeline/wiki), multiple rounds of this script were used, incorporating results from each successive round to build a comprehensive transcriptome database. 

## Trinotate Steps

