#!/usr/bin/env bash
# Datos del tutorial "Parkinson's Mouse" de QIIME 2 (versión 2021.4): trasplante
# de microbiota fecal de donantes sanos o con Parkinson a ratones wild type y
# susceptibles (ASO). 16S rRNA V4, Illumina single-end 150 pb.
set -euo pipefail
V=2021.4
mkdir -p data && cd data
wget -nc -O sample_metadata.tsv  "https://data.qiime2.org/${V}/tutorials/pd-mice/sample_metadata.tsv"
wget -nc -O demultiplexed_seqs.zip "https://data.qiime2.org/${V}/tutorials/pd-mice/demultiplexed_seqs.zip"
wget -nc -O manifest              "https://data.qiime2.org/${V}/tutorials/pd-mice/manifest"
unzip -o -q demultiplexed_seqs.zip
# El manifest usa rutas $PWD/...: se sustituyen por la ruta absoluta actual
sed "s|\$PWD|${PWD}|g" manifest > manifest.tsv
wget -nc "https://data.qiime2.org/${V}/common/gg-13-8-99-515-806-nb-classifier.qza"
wget -nc "https://data.qiime2.org/${V}/common/sepp-refs-gg-13-8.qza"
