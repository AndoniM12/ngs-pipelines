#!/usr/bin/env bash
# Datos públicos del estudio GSE60450 (Fu et al., Nat Cell Biol 2015).
set -euo pipefail
mkdir -p data && cd data
# Lecturas de una muestra (single-end 100 pb, HiSeq 2000)
fasterq-dump --threads 4 SRR1552444 && gzip -f SRR1552444.fastq
# Índice de HISAT2 para mm10 y anotación GENCODE vM10 (misma nomenclatura chrN)
wget -nc https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz && tar -xzf mm10_genome.tar.gz
wget -nc https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.annotation.gtf.gz
gunzip -kf gencode.vM10.annotation.gtf.gz
# Matriz de recuentos de las 12 muestras (para el análisis de expresión diferencial)
wget -nc https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60450/suppl/GSE60450_Lactation-GenewiseCounts.txt.gz
