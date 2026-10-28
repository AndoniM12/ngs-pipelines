#!/usr/bin/env bash
# =============================================================================
# RNA-seq: de lecturas a recuentos por gen (single-end, Illumina)
#
#   FastQC -> Trim Galore -> FastQC -> HISAT2 -> SAMtools -> Qualimap rnaseq
#   -> htseq-count -> resumen de asignación
#
# Autor: Andoni Moreno Lanceta
# Uso:   ./run_rnaseq_counts.sh -r SRR1552444.fastq.gz -x mm10/genome \
#                               -g gencode.vM10.annotation.gtf -s SRR1552444 [-o resultados] [-t 4]
# =============================================================================
set -euo pipefail

OUTDIR="resultados"
THREADS=4
HISAT_K=1          # alineamientos por lectura que reporta HISAT2 (1 = como en la práctica)
STRANDED="no"      # TruSeq RNA v2 (poli-A, sin hebra) -> "no"; TruSeq Stranded -> "reverse"
MINAQUAL=10        # htseq-count: MAPQ mínimo (HISAT2 da MAPQ 0/1 a multimapeo)

usage() {
  grep '^# Uso' -A1 "$0" | sed 's/^# //'
  cat <<EOF

Opciones:
  -r  FASTQ single-end (gz)            -x  prefijo del índice de HISAT2
  -g  anotación GTF                     -s  nombre de la muestra
  -o  directorio de salida [${OUTDIR}]  -t  hilos [${THREADS}]
  -k  alineamientos por lectura en HISAT2 [${HISAT_K}]
  -S  hebra para htseq-count: no | yes | reverse [${STRANDED}]
EOF
  exit 1
}

while getopts "r:x:g:s:o:t:k:S:h" opt; do
  case $opt in
    r) READS=$OPTARG ;;  x) INDEX=$OPTARG ;;
    g) GTF=$OPTARG ;;    s) SAMPLE=$OPTARG ;;
    o) OUTDIR=$OPTARG ;; t) THREADS=$OPTARG ;;
    k) HISAT_K=$OPTARG ;; S) STRANDED=$OPTARG ;;
    *) usage ;;
  esac
done
[[ -z "${READS:-}" || -z "${INDEX:-}" || -z "${GTF:-}" || -z "${SAMPLE:-}" ]] && usage
for f in "$READS" "$GTF" "${INDEX}.1.ht2"; do [[ -s "$f" ]] || { echo "No existe: $f" >&2; exit 1; }; done

log() { echo "[$(date '+%F %T')] $*" | tee -a "${OUTDIR}/pipeline.log"; }
mkdir -p "${OUTDIR}"/{01_qc_raw,02_trimming,03_alignment,04_counts}
: > "${OUTDIR}/pipeline.log"   # registro nuevo en cada ejecución, sin restos de las anteriores
log "Muestra: ${SAMPLE} | Índice: ${INDEX} | GTF: ${GTF}"

# ---------------------------- 1. Calidad inicial -----------------------------
log "1. FastQC sobre las lecturas crudas"
fastqc -t "$THREADS" --outdir "${OUTDIR}/01_qc_raw" "$READS"

# ---------------------------- 2. Recorte -------------------------------------
# Parámetros por defecto: detección automática del adaptador, Q20 en el extremo
# 3', longitud mínima 20 pb.
log "2. Trim Galore (parámetros por defecto) + FastQC"
trim_galore --cores "$THREADS" --fastqc -o "${OUTDIR}/02_trimming" "$READS"
TRIMMED=$(ls "${OUTDIR}"/02_trimming/*_trimmed.fq.gz)

# ---------------------------- 3. Alineamiento --------------------------------
log "3. HISAT2 (-k ${HISAT_K}), ordenación e indexado del BAM"
BAM="${OUTDIR}/03_alignment/${SAMPLE}.sorted.bam"
hisat2 -p "$THREADS" -k "$HISAT_K" -x "$INDEX" -U "$TRIMMED" \
       2> "${OUTDIR}/03_alignment/${SAMPLE}.hisat2.log" \
  | samtools sort -@ "$THREADS" -o "$BAM" -
samtools index "$BAM"
cat "${OUTDIR}/03_alignment/${SAMPLE}.hisat2.log" | tee -a "${OUTDIR}/pipeline.log"

# Distribución de FLAG y MAPQ (preguntas de la práctica)
samtools view "$BAM" | cut -f2 | sort | uniq -c | sort -k1,1nr \
  > "${OUTDIR}/03_alignment/${SAMPLE}.flags.txt"
samtools view "$BAM" | cut -f5 | sort -n | uniq -c \
  > "${OUTDIR}/03_alignment/${SAMPLE}.mapq.txt"
samtools flagstat "$BAM" > "${OUTDIR}/03_alignment/${SAMPLE}.flagstat.txt"

# ---------------------------- 4. Calidad del alineamiento --------------------
log "4. Qualimap rnaseq (origen genómico de las lecturas y sesgo 5'-3')"
qualimap rnaseq -bam "$BAM" -gtf "$GTF" -outdir "${OUTDIR}/03_alignment/qualimap" \
         --java-mem-size=4G

# ---------------------------- 5. Recuento por gen ----------------------------
log "5. htseq-count (exones agrupados por gene_id, hebra: ${STRANDED})"
COUNTS="${OUTDIR}/04_counts/${SAMPLE}.counts.tsv"
htseq-count -f bam -r pos -s "$STRANDED" -t exon -i gene_id -a "$MINAQUAL" \
            "$BAM" "$GTF" > "$COUNTS"

# Resumen: lecturas asignadas a genes frente a las categorías especiales (__)
awk -F'\t' '
  /^__/ { cat[$1]=$2; total+=$2; next }
        { assigned+=$2; total+=$2 }
  END {
    printf "categoria\tlecturas\tporcentaje\n"
    printf "assigned\t%d\t%.2f\n", assigned, 100*assigned/total
    for (c in cat) printf "%s\t%d\t%.2f\n", substr(c,3), cat[c], 100*cat[c]/total
  }' "$COUNTS" > "${OUTDIR}/04_counts/${SAMPLE}.assignment_summary.tsv"
tee -a "${OUTDIR}/pipeline.log" < "${OUTDIR}/04_counts/${SAMPLE}.assignment_summary.tsv"

log "Fin. Recuentos: ${COUNTS}"
