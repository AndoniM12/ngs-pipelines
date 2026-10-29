#!/usr/bin/env bash
# Descarga las lecturas públicas del aislado SAL4 de Salmonella enterica serovar
# Typhimurium (queso, Michoacán; BioProject PRJNA802840, run SRR18072959) y el
# genoma de referencia S. Typhimurium LT2 (RefSeq) en GenBank.
set -euo pipefail
mkdir -p data && cd data
fasterq-dump --split-files --threads 4 SRR18072959
gzip -f SRR18072959_1.fastq SRR18072959_2.fastq
datasets download genome accession GCF_000006945.2 --include gbff --filename LT2.zip
unzip -o LT2.zip -d LT2
cp LT2/ncbi_dataset/data/GCF_000006945.2/genomic.gbff LT2_genomic.gbk
