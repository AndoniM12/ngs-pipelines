#!/usr/bin/env bash
# =============================================================================
# Pipeline de llamada y anotación de variantes germinales en panel de captura
# Illumina paired-end · Referencia GRCh38 (Ensembl)
#
#   FastQC -> fastp -> FastQC -> BWA-MEM -> SAMtools -> flagstat/Qualimap
#   -> Picard MarkDuplicates -> FreeBayes -> RTG vcfstats -> VEP
#
# Autor: Andoni Moreno Lanceta
# Uso:   ./run_variantes_panel.sh -1 R1.fastq.gz -2 R2.fastq.gz -r ref.fa \
#                                 -b panel.bed -s NA12880 [-o resultados] [-t 4]
# =============================================================================
set -euo pipefail

# ---------------------------- Parámetros por defecto -------------------------
OUTDIR="resultados"
THREADS=4
MIN_LEN=100          # fastp: longitud mínima tras el recorte
MIN_QUAL=28          # fastp: calidad media de la ventana (inicio, final y media)
MIN_ALT_COUNT=30     # freebayes -C: lecturas mínimas con el alelo alternativo
                     # (30 = valor de la práctica; ver README, es muy restrictivo)
VEP_CACHE="${HOME}/.vep"

usage() {
  grep '^# Uso' -A1 "$0" | sed 's/^# //'
  cat <<EOF

Opciones:
  -1  FASTQ R1 (gz)                 -2  FASTQ R2 (gz)
  -r  genoma de referencia (FASTA)  -b  regiones de captura (BED)
  -s  nombre de la muestra          -o  directorio de salida [${OUTDIR}]
  -t  hilos [${THREADS}]            -C  min. lecturas alelo alternativo [${MIN_ALT_COUNT}]
  -c  directorio de la caché de VEP [${VEP_CACHE}]
EOF
  exit 1
}

while getopts "1:2:r:b:s:o:t:C:c:h" opt; do
  case $opt in
    1) R1=$OPTARG ;;        2) R2=$OPTARG ;;
    r) REF=$OPTARG ;;       b) BED=$OPTARG ;;
    s) SAMPLE=$OPTARG ;;    o) OUTDIR=$OPTARG ;;
    t) THREADS=$OPTARG ;;   C) MIN_ALT_COUNT=$OPTARG ;;
    c) VEP_CACHE=$OPTARG ;; *) usage ;;
  esac
done
[[ -z "${R1:-}" || -z "${R2:-}" || -z "${REF:-}" || -z "${BED:-}" || -z "${SAMPLE:-}" ]] && usage
for f in "$R1" "$R2" "$REF" "$BED"; do [[ -s "$f" ]] || { echo "No existe: $f" >&2; exit 1; }; done

log() { echo "[$(date '+%F %T')] $*" | tee -a "${OUTDIR}/pipeline.log"; }

mkdir -p "${OUTDIR}"/{01_qc_raw,02_trimming,03_mapping,04_dedup,05_variants,06_vep}
log "Muestra: ${SAMPLE} | Referencia: ${REF} | Panel: ${BED} | Hilos: ${THREADS}"

# ---------------------------- 1. Calidad inicial -----------------------------
log "1. FastQC sobre las lecturas crudas"
fastqc -t "$THREADS" -o "${OUTDIR}/01_qc_raw" "$R1" "$R2"

# ---------------------------- 2. Limpieza ------------------------------------
# --cut_front/--cut_tail: ventana deslizante en los extremos; --cut_mean_quality
# fija el umbral de calidad media de la ventana.
log "2. Limpieza con fastp (Q${MIN_QUAL}, longitud mínima ${MIN_LEN} pb)"
T1="${OUTDIR}/02_trimming/${SAMPLE}_R1.trim.fastq.gz"
T2="${OUTDIR}/02_trimming/${SAMPLE}_R2.trim.fastq.gz"
fastp -i "$R1" -I "$R2" -o "$T1" -O "$T2" \
      --length_required "$MIN_LEN" \
      --cut_front --cut_tail --cut_mean_quality "$MIN_QUAL" \
      --detect_adapter_for_pe --trim_poly_g --trim_poly_x \
      --thread "$THREADS" \
      --html "${OUTDIR}/02_trimming/${SAMPLE}.fastp.html" \
      --json "${OUTDIR}/02_trimming/${SAMPLE}.fastp.json"
fastqc -t "$THREADS" -o "${OUTDIR}/02_trimming" "$T1" "$T2"

# ---------------------------- 3. Referencia y mapeo --------------------------
log "3. Indexado de la referencia (solo si falta)"
[[ -s "${REF}.bwt" ]] || bwa index "$REF"
[[ -s "${REF}.fai" ]] || samtools faidx "$REF"

# Read group a partir de la cabecera del FASTQ (Illumina):
#   @instrumento:run:flowcell:carril:...  ->  ID = flowcell.carril
# "|| true": con pipefail, el SIGPIPE de zcat al cortar head abortaría el script
HEADER=$(zcat "$R1" | head -n 1 || true)
FLOWCELL=$(echo "$HEADER" | cut -d: -f3)
LANE=$(echo "$HEADER" | cut -d: -f4)
RG="@RG\tID:${FLOWCELL}.${LANE}\tSM:${SAMPLE}\tLB:${SAMPLE}_lib1\tPL:ILLUMINA\tPU:${FLOWCELL}.${LANE}.${SAMPLE}"
log "   Read group: ${RG}"

log "3. Mapeo con BWA-MEM, ordenación e indexado del BAM"
BAM="${OUTDIR}/03_mapping/${SAMPLE}.sorted.bam"
bwa mem -t "$THREADS" -R "$RG" "$REF" "$T1" "$T2" 2> "${OUTDIR}/03_mapping/bwa.log" \
  | samtools sort -@ "$THREADS" -o "$BAM" -
samtools index "$BAM"

# ---------------------------- 4. Estadísticas de mapeo -----------------------
log "4. samtools flagstat y Qualimap sobre las regiones del panel"
samtools flagstat "$BAM" > "${OUTDIR}/03_mapping/${SAMPLE}.flagstat.txt"

# Armonizar la nomenclatura de cromosomas entre BED (UCSC: chr12) y
# referencia (Ensembl: 12) sin modificar el BED original.
BED_USE="${OUTDIR}/03_mapping/panel.harmonized.bed"
if awk 'NR==1{exit !($1 ~ /^chr/)}' "${REF}.fai"; then
  cp "$BED" "$BED_USE"
else
  sed 's/^chr//' "$BED" > "$BED_USE"
fi
# Restringir el BED a los contigs presentes en la referencia (p. ej. solo chr12)
cut -f1 "${REF}.fai" > "${OUTDIR}/03_mapping/contigs.txt"
awk 'NR==FNR{c[$1]; next} ($1 in c)' "${OUTDIR}/03_mapping/contigs.txt" "$BED_USE" \
  > "${BED_USE}.tmp" && mv "${BED_USE}.tmp" "$BED_USE"

qualimap bamqc -bam "$BAM" -gff "$BED_USE" -c -nt "$THREADS" -gd HUMAN \
         -outdir "${OUTDIR}/03_mapping/qualimap" --java-mem-size=4G

# ---------------------------- 5. Duplicados ----------------------------------
log "5. Marcado de duplicados con Picard"
DEDUP="${OUTDIR}/04_dedup/${SAMPLE}.dedup.bam"
picard MarkDuplicates -I "$BAM" -O "$DEDUP" \
       -M "${OUTDIR}/04_dedup/${SAMPLE}.markdup_metrics.txt" \
       --ASSUME_SORT_ORDER coordinate
samtools index "$DEDUP"

# ---------------------------- 6. Llamada de variantes ------------------------
log "6. Llamada de variantes con FreeBayes en las regiones del panel (-C ${MIN_ALT_COUNT})"
VCF="${OUTDIR}/05_variants/${SAMPLE}.vcf"
freebayes -f "$REF" -t "$BED_USE" -C "$MIN_ALT_COUNT" "$DEDUP" > "$VCF"
bgzip -f "$VCF" && tabix -p vcf "${VCF}.gz"
rtg vcfstats "${VCF}.gz" > "${OUTDIR}/05_variants/${SAMPLE}.vcfstats.txt"

# ---------------------------- 7. Anotación con VEP ---------------------------
log "7. Anotación funcional con VEP (GRCh38, caché offline)"
VEP_OUT="${OUTDIR}/06_vep/${SAMPLE}.vep.vcf"
vep -i "${VCF}.gz" -o "$VEP_OUT" --vcf \
    --species homo_sapiens --assembly GRCh38 \
    --cache --offline --dir_cache "$VEP_CACHE" \
    --everything --fork "$THREADS" \
    --stats_file "${OUTDIR}/06_vep/${SAMPLE}.vep_summary.html"

# Ejemplo de filtrado (pregunta de la práctica): variantes stop_gained en VWF
filter_vep -i "$VEP_OUT" -o "${OUTDIR}/06_vep/${SAMPLE}.VWF_stop_gained.vcf" \
           --filter "SYMBOL is VWF and Consequence is stop_gained" --force_overwrite

log "Fin. Para IGV: ${DEDUP}, ${VCF}.gz y ${BED_USE}"
