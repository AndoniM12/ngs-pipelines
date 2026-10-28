#!/usr/bin/env bash
# Descarga el cromosoma 12 de GRCh38 (Ensembl release 115).
# Las lecturas de NA12880 y el BED del panel TruSight One proceden del material
# del máster y no se incluyen en el repositorio.
set -euo pipefail
mkdir -p ref && cd ref
wget -nc https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.12.fa.gz
gunzip -kf Homo_sapiens.GRCh38.dna.chromosome.12.fa.gz
