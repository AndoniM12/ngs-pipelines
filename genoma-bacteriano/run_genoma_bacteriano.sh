#!/usr/bin/env bash
# =============================================================================
# Pipeline de ensamblaje de novo y caracterización de genomas bacterianos
# Illumina paired-end
#
#   FastQC -> fastp -> FastQC -> SPAdes -> (filtrado de contigs) -> QUAST
#   -> Prokka (proteoma de referencia) -> AMRFinderPlus
#
# Autor: Andoni Moreno Lanceta
# Uso:   ./run_genoma_bacteriano.sh -1 R1.fastq.gz -2 R2.fastq.gz \
#             -g referencia.gbk -s SRR18072959 [-O Salmonella] [-o resultados] [-t 4]
# =============================================================================
set -euo pipefail

# ---------------------------- Parámetros por defecto -------------------------
OUTDIR="resultados"
THREADS=4
MIN_LEN=151          # fastp: longitud mínima tras el recorte
MIN_QUAL=25          # fastp: calidad media de la ventana (inicio, final y media)
MIN_CONTIG=0         # contigs más cortos se descartan antes de anotar
                     # (0 = como en la práctica; recomendado 500, ver README)
ORGANISM="Salmonella"
GENUS="Salmonella"
SPECIES="enterica"
AMR_IDENT=0.9        # AMRFinderPlus: identidad mínima
AMR_COV=0.9          # AMRFinderPlus: cobertura mínima

usage() {
  grep '^# Uso' -A1 "$0" | sed 's/^# //'
  cat <<EOF

Opciones:
  -1  FASTQ R1 (gz)                    -2  FASTQ R2 (gz)
  -g  genoma de referencia anotado (GenBank) para Prokka
  -s  nombre de la muestra             -O  organismo para AMRFinderPlus [${ORGANISM}]
  -o  directorio de salida [${OUTDIR}] -t  hilos [${THREADS}]
  -m  longitud mínima de contig [${MIN_CONTIG}]
EOF
  exit 1
}

while getopts "1:2:g:s:O:o:t:m:h" opt; do
  case $opt in
    1) R1=$OPTARG ;;       2) R2=$OPTARG ;;
    g) REF_GBK=$OPTARG ;;  s) SAMPLE=$OPTARG ;;
    O) ORGANISM=$OPTARG ;; o) OUTDIR=$OPTARG ;;
    t) THREADS=$OPTARG ;;  m) MIN_CONTIG=$OPTARG ;;
    *) usage ;;
  esac
done
[[ -z "${R1:-}" || -z "${R2:-}" || -z "${REF_GBK:-}" || -z "${SAMPLE:-}" ]] && usage
for f in "$R1" "$R2" "$REF_GBK"; do [[ -s "$f" ]] || { echo "No existe: $f" >&2; exit 1; }; done

log() { echo "[$(date '+%F %T')] $*" | tee -a "${OUTDIR}/pipeline.log"; }

mkdir -p "${OUTDIR}"/{01_qc_raw,02_trimming,03_assembly,04_quast,05_annotation,06_amr}
log "Muestra: ${SAMPLE} | Referencia: ${REF_GBK} | Hilos: ${THREADS}"

# ---------------------------- 1. Calidad inicial -----------------------------
log "1. FastQC sobre las lecturas crudas"
fastqc -t "$THREADS" -o "${OUTDIR}/01_qc_raw" "$R1" "$R2"

# ---------------------------- 2. Limpieza ------------------------------------
log "2. Limpieza con fastp (Q${MIN_QUAL}, longitud mínima ${MIN_LEN} pb)"
T1="${OUTDIR}/02_trimming/${SAMPLE}_R1.trim.fastq.gz"
T2="${OUTDIR}/02_trimming/${SAMPLE}_R2.trim.fastq.gz"
fastp -i "$R1" -I "$R2" -o "$T1" -O "$T2" \
      --length_required "$MIN_LEN" \
      --cut_front --cut_tail --cut_mean_quality "$MIN_QUAL" \
      --detect_adapter_for_pe --thread "$THREADS" \
      --html "${OUTDIR}/02_trimming/${SAMPLE}.fastp.html" \
      --json "${OUTDIR}/02_trimming/${SAMPLE}.fastp.json"
fastqc -t "$THREADS" -o "${OUTDIR}/02_trimming" "$T1" "$T2"

# ---------------------------- 3. Ensamblaje de novo --------------------------
# --careful: corrección de errores (BayesHammer) + ensamblaje + corrección de
# mismatches/indels (MismatchCorrector). k-meros automáticos según la longitud.
log "3. Ensamblaje con SPAdes (--careful, k automático)"
ASM_DIR="${OUTDIR}/03_assembly/spades"
spades.py -1 "$T1" -2 "$T2" --careful -k auto -t "$THREADS" -o "$ASM_DIR"

CONTIGS="${OUTDIR}/03_assembly/${SAMPLE}.contigs.fasta"
if (( MIN_CONTIG > 0 )); then
  log "   Filtrando contigs < ${MIN_CONTIG} pb"
  seqkit seq -m "$MIN_CONTIG" "${ASM_DIR}/contigs.fasta" > "$CONTIGS"
else
  cp "${ASM_DIR}/contigs.fasta" "$CONTIGS"
fi
log "   Contigs: $(grep -c '>' "$CONTIGS")"

# ---------------------------- 4. Evaluación del ensamblaje -------------------
log "4. QUAST (todos los contigs, con mapeo de las lecturas)"
quast.py -o "${OUTDIR}/04_quast" -m 0 -t "$THREADS" \
         --pe1 "$T1" --pe2 "$T2" "$CONTIGS" "${ASM_DIR}/scaffolds.fasta"

# ---------------------------- 5. Anotación -----------------------------------
log "5. Proteoma de referencia y anotación con Prokka"
REF_PROTS="${OUTDIR}/05_annotation/ref_prots.faa"
prokka-genbank_to_fasta_db "$REF_GBK" > "$REF_PROTS"
log "   Proteínas de referencia: $(grep -c '>' "$REF_PROTS")"

prokka --outdir "${OUTDIR}/05_annotation/prokka" --prefix "$SAMPLE" --force \
       --kingdom Bacteria --genus "$GENUS" --species "$SPECIES" --usegenus \
       --proteins "$REF_PROTS" --addgenes --addmrna --compliant \
       --mincontiglen "$MIN_CONTIG" --cpus "$THREADS" "$CONTIGS"

# ---------------------------- 6. Resistencias y virulencia -------------------
# Modo combinado (proteínas + nucleótidos + GFF): detecta genes adquiridos y
# mutaciones puntuales, y evita duplicados entre ambos tipos de búsqueda.
log "6. AMRFinderPlus (--plus, identidad ${AMR_IDENT}, cobertura ${AMR_COV})"
P="${OUTDIR}/05_annotation/prokka/${SAMPLE}"
amrfinder --update || log "   Aviso: no se pudo actualizar la base de datos de AMRFinderPlus"
amrfinder -p "${P}.faa" -n "${P}.fna" -g "${P}.gff" --annotation_format prokka \
          --organism "$ORGANISM" --plus -i "$AMR_IDENT" -c "$AMR_COV" \
          --threads "$THREADS" -o "${OUTDIR}/06_amr/${SAMPLE}.amrfinder.tsv"
amrfinder --database_version > "${OUTDIR}/06_amr/amrfinder_version.txt" 2>&1 || true

log "Fin. Resumen QUAST: ${OUTDIR}/04_quast/report.txt"
