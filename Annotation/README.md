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

The mapping results were incorporated into the PASA pipeline using the *Launch_PASA_pipeline.pl* script. As per [guidelines](https://github.com/PASApipeline/PASApipeline/wiki), multiple rounds of this script were used, incorporating results from each successive round to build a comprehensive transcriptome database. The details are included in the file *run_pasa_on_transcriptome.txt*.

## Trinotate: Transcriptome Functional Annotation

The Trinotate pipeline was used to annotate the transcripts mapped to the genome. The methods followed the guidelines on their [webpage](https://github.com/Trinotate/Trinotate.github.io/wiki). 

First, [TransDecoder](https://github.com/TransDecoder) version 5.0.2  was used to predict coding regions from the assemblies produced by the PASA analysis. Predicted peptide sequences and transcripts were then used as queries to search multiple protein and nucleotide databases. The protein database SwissProt was queried using BLAST, using blastp for peptide sequences and blastx for transcripts. Additionally, protein databases of the Zebrafish and Tilapia genomes were searched. The program [HMMER 3.1b2](http://hmmer.org/) was used to identify protein domains in the peptide sequences using the Pfam database. Search results were consolidated into a Trinotate sqlite database to produce an annotation report. Details of the Trinotate steps are included in the file *trinotate_on_transcriptome.txt*

## Visualising Annotations

A custom Python script (included here as the [Jupyter](http://jupyter.org/) notebook *Trinotate_annotations_to_SeqMonk-final.ipynb*) was used to convert the gff3 annotation files so the gene names and information could be visualised easily in the program [SeqMonk](https://www.bioinformatics.babraham.ac.uk/projects/seqmonk/). This step allowed for easy tracking of which genes were undergoing differential methylation and/or gene expression. 

Additionally, tables were produced of all annotations that included Ensembl gene descriptions in order to catalog and better explore all functional gene information. These tables were produced using Python and R scripts. Separate tables were made for sprot, and for zebrafish and tilapia genome matches, as well as one that combined all information. The scripts used to produce the zebrafish table are provided (*Produce_table_of_zebrafish_annots_w_ensembl_info.ipynb* and *add_annots_biomart.R*).






