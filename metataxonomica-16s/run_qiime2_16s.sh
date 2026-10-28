#!/usr/bin/env bash
# =============================================================================
# Metataxonómica 16S rRNA (región V4, single-end) con QIIME 2
#
#   importación (manifest) -> resumen de calidad -> DADA2 -> tabla de ASVs
#   -> clasificación taxonómica (Naive Bayes, Greengenes 13_8) -> barplot
#   -> árbol filogenético (SEPP) -> rarefacción -> diversidad alfa y beta
#   -> exportación para R (phyloseq)
#
# Autor: Andoni Moreno Lanceta
# Uso:   ./run_qiime2_16s.sh -m manifest.tsv -d metadata.tsv -c classifier.qza \
#                            -r sepp-refs.qza [-o resultados] [-t 4]
# Entorno: qiime2-2021.4 (ver README)
# =============================================================================
set -euo pipefail

OUTDIR="resultados"
THREADS=4
TRIM_LEFT=0          # DADA2: bases a recortar al inicio
TRUNC_LEN=150        # DADA2: longitud de truncado (revisar en demux.qzv)
SAMPLING_DEPTH=2000  # profundidad de rarefacción para core-metrics
MAX_DEPTH=4000       # profundidad máxima de la curva de rarefacción
GROUP_COL="donor_status"

usage() {
  grep '^# Uso' -A1 "$0" | sed 's/^# //'
  cat <<EOF

Opciones:
  -m  manifest (SingleEndFastqManifestPhred33V2)   -d  metadatos de las muestras
  -c  clasificador taxonómico (.qza)               -r  referencias SEPP (.qza)
  -o  directorio de salida [${OUTDIR}]             -t  hilos [${THREADS}]
  -L  truncado DADA2 [${TRUNC_LEN}]                -D  profundidad de rarefacción [${SAMPLING_DEPTH}]
  -M  profundidad máxima de rarefacción [${MAX_DEPTH}]
  -g  columna de metadatos para comparar grupos [${GROUP_COL}]
EOF
  exit 1
}

while getopts "m:d:c:r:o:t:L:D:M:g:h" opt; do
  case $opt in
    m) MANIFEST=$OPTARG ;;   d) METADATA=$OPTARG ;;
    c) CLASSIFIER=$OPTARG ;; r) SEPP_REFS=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;     t) THREADS=$OPTARG ;;
    L) TRUNC_LEN=$OPTARG ;;  D) SAMPLING_DEPTH=$OPTARG ;;
    M) MAX_DEPTH=$OPTARG ;;  g) GROUP_COL=$OPTARG ;;
    *) usage ;;
  esac
done
[[ -z "${MANIFEST:-}" || -z "${METADATA:-}" || -z "${CLASSIFIER:-}" || -z "${SEPP_REFS:-}" ]] && usage
for f in "$MANIFEST" "$METADATA" "$CLASSIFIER" "$SEPP_REFS"; do [[ -s "$f" ]] || { echo "No existe: $f" >&2; exit 1; }; done

log() { echo "[$(date '+%F %T')] $*" | tee -a "${OUTDIR}/pipeline.log"; }
mkdir -p "${OUTDIR}"/{01_import,02_dada2,03_taxonomy,04_phylogeny,05_diversity,06_export}
: > "${OUTDIR}/pipeline.log"   # registro nuevo en cada ejecución, sin restos de las anteriores
log "Manifest: ${MANIFEST} | Metadatos: ${METADATA} | Hilos: ${THREADS}"

# ---------------------------- 0. Metadatos -----------------------------------
qiime metadata tabulate --m-input-file "$METADATA" \
  --o-visualization "${OUTDIR}/01_import/metadata.qzv"

# ---------------------------- 1. Importación ---------------------------------
# Las lecturas ya vienen demultiplexadas (una FASTQ por muestra). Con lecturas
# multiplexadas EMP se usaría: qiime demux emp-single --i-seqs ... \
#   --m-barcodes-file metadata.tsv --m-barcodes-column BarcodeSequence ...
log "1. Importación de las lecturas y resumen de calidad"
qiime tools import --type 'SampleData[SequencesWithQuality]' \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path "$MANIFEST" --output-path "${OUTDIR}/01_import/demux_seqs.qza"
qiime demux summarize --i-data "${OUTDIR}/01_import/demux_seqs.qza" \
  --o-visualization "${OUTDIR}/01_import/demux_seqs.qzv"

# ---------------------------- 2. DADA2 ---------------------------------------
log "2. DADA2: filtrado, corrección de errores, ASVs y eliminación de quimeras (truncado ${TRUNC_LEN})"
qiime dada2 denoise-single --i-demultiplexed-seqs "${OUTDIR}/01_import/demux_seqs.qza" \
  --p-trim-left "$TRIM_LEFT" --p-trunc-len "$TRUNC_LEN" --p-n-threads "$THREADS" \
  --o-table "${OUTDIR}/02_dada2/table.qza" \
  --o-representative-sequences "${OUTDIR}/02_dada2/rep_seqs.qza" \
  --o-denoising-stats "${OUTDIR}/02_dada2/denoising_stats.qza"
qiime metadata tabulate --m-input-file "${OUTDIR}/02_dada2/denoising_stats.qza" \
  --o-visualization "${OUTDIR}/02_dada2/denoising_stats.qzv"
qiime feature-table summarize --i-table "${OUTDIR}/02_dada2/table.qza" \
  --m-sample-metadata-file "$METADATA" --o-visualization "${OUTDIR}/02_dada2/table.qzv"
qiime feature-table tabulate-seqs --i-data "${OUTDIR}/02_dada2/rep_seqs.qza" \
  --o-visualization "${OUTDIR}/02_dada2/rep_seqs.qzv"

# ---------------------------- 3. Taxonomía -----------------------------------
log "3. Clasificación taxonómica (Naive Bayes) y gráfico de barras"
qiime feature-classifier classify-sklearn --i-classifier "$CLASSIFIER" \
  --i-reads "${OUTDIR}/02_dada2/rep_seqs.qza" --p-n-jobs "$THREADS" \
  --o-classification "${OUTDIR}/03_taxonomy/taxonomy.qza"
qiime metadata tabulate --m-input-file "${OUTDIR}/03_taxonomy/taxonomy.qza" \
  --o-visualization "${OUTDIR}/03_taxonomy/taxonomy.qzv"
qiime taxa barplot --i-table "${OUTDIR}/02_dada2/table.qza" \
  --i-taxonomy "${OUTDIR}/03_taxonomy/taxonomy.qza" --m-metadata-file "$METADATA" \
  --o-visualization "${OUTDIR}/03_taxonomy/taxa_barplot.qzv"

# ---------------------------- 4. Filogenia -----------------------------------
log "4. Inserción de fragmentos en el árbol de referencia (SEPP)"
qiime fragment-insertion sepp --i-representative-sequences "${OUTDIR}/02_dada2/rep_seqs.qza" \
  --i-reference-database "$SEPP_REFS" --p-threads "$THREADS" \
  --o-tree "${OUTDIR}/04_phylogeny/tree.qza" \
  --o-placements "${OUTDIR}/04_phylogeny/tree_placements.qza"

# ---------------------------- 5. Diversidad ----------------------------------
log "5. Rarefacción (hasta ${MAX_DEPTH}) y métricas de diversidad (profundidad ${SAMPLING_DEPTH})"
qiime diversity alpha-rarefaction --i-table "${OUTDIR}/02_dada2/table.qza" \
  --i-phylogeny "${OUTDIR}/04_phylogeny/tree.qza" --m-metadata-file "$METADATA" \
  --p-max-depth "$MAX_DEPTH" --o-visualization "${OUTDIR}/05_diversity/alpha_rarefaction.qzv"

# Las muestras con menos lecturas que SAMPLING_DEPTH se excluyen de este paso.
qiime diversity core-metrics-phylogenetic --i-table "${OUTDIR}/02_dada2/table.qza" \
  --i-phylogeny "${OUTDIR}/04_phylogeny/tree.qza" --m-metadata-file "$METADATA" \
  --p-sampling-depth "$SAMPLING_DEPTH" --p-n-jobs-or-threads "$THREADS" \
  --output-dir "${OUTDIR}/05_diversity/core_metrics"

CM="${OUTDIR}/05_diversity/core_metrics"
for a in shannon_vector faith_pd_vector; do
  qiime diversity alpha-group-significance --i-alpha-diversity "${CM}/${a}.qza" \
    --m-metadata-file "$METADATA" --o-visualization "${CM}/${a}_group_significance.qzv"
done
for b in unweighted_unifrac weighted_unifrac bray_curtis; do
  qiime diversity beta-group-significance --i-distance-matrix "${CM}/${b}_distance_matrix.qza" \
    --m-metadata-file "$METADATA" --m-metadata-column "$GROUP_COL" --p-pairwise \
    --o-visualization "${CM}/${b}_${GROUP_COL}_permanova.qzv"
done

# ---------------------------- 6. Exportación ---------------------------------
log "6. Exportación de tabla, taxonomía y árbol para R"
E="${OUTDIR}/06_export"
qiime tools export --input-path "${OUTDIR}/02_dada2/table.qza" --output-path "$E"
biom convert -i "${E}/feature-table.biom" -o "${E}/ASV.tsv" --to-tsv
qiime tools export --input-path "${OUTDIR}/03_taxonomy/taxonomy.qza" --output-path "$E"
qiime tools export --input-path "${OUTDIR}/04_phylogeny/tree.qza" --output-path "$E"

log "Fin. Visualizaciones (.qzv) en https://view.qiime2.org; datos para R en ${E}"
