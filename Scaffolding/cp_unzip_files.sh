#!/usr/bin/env bash

# copy files to working folder
cp /scratch/bluehead_wrasse/trimmomatic/trailing26/H3WYYBCXX-1779-01-26-1_L002_R1.fastq.gz.trimmed.xz /data/bluehead_scaffold/trim_data/
cp /scratch/bluehead_wrasse/trimmomatic/trailing26/H3WYYBCXX-1779-01-27-1_L001_R1.fastq.gz.trimmed.xz /data/bluehead_scaffold/trim_data/
cp /scratch/bluehead_wrasse/trimmomatic/trailing26/H3WYYBCXX-1779-01-27-1_L001_R2.fastq.gz.trimmed.xz /data/bluehead_scaffold/trim_data/
cp /scratch/bluehead_wrasse/trimmomatic/trailing26/H3WYYBCXX-1779-01-27-1_L002_R1.fastq.gz.trimmed.xz /data/bluehead_scaffold/trim_data/
cp /scratch/bluehead_wrasse/trimmomatic/trailing26/H3WYYBCXX-1779-01-27-1_L002_R2.fastq.gz.trimmed.xz /data/bluehead_scaffold/trim_data/

cp /scratch/bluehead_wrasse/trimmomatic/trailing26/H3WYYBCXX-1779-01-26-1_L002_R1.fastq.gz.unpaired.xz /data/bluehead_scaffold/trim_data/
cp /scratch/bluehead_wrasse/trimmomatic/trailing26/H3WYYBCXX-1779-01-26-1_L002_R2.fastq.gz.unpaired.xz /data/bluehead_scaffold/trim_data/
cp /scratch/bluehead_wrasse/trimmomatic/trailing26/H3WYYBCXX-1779-01-27-1_L001_R1.fastq.gz.unpaired.xz /data/bluehead_scaffold/trim_data/
cp /scratch/bluehead_wrasse/trimmomatic/trailing26/H3WYYBCXX-1779-01-27-1_L001_R2.fastq.gz.unpaired.xz /data/bluehead_scaffold/trim_data/
cp /scratch/bluehead_wrasse/trimmomatic/trailing26/H3WYYBCXX-1779-01-27-1_L002_R1.fastq.gz.unpaired.xz /data/bluehead_scaffold/trim_data/
cp /scratch/bluehead_wrasse/trimmomatic/trailing26/H3WYYBCXX-1779-01-27-1_L002_R2.fastq.gz.unpaired.xz /data/bluehead_scaffold/trim_data/

# unzip files
unxz *.xz


