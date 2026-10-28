#!/usr/bin/env Rscript
# =============================================================================
# RNA-seq: de recuentos a genes diferencialmente expresados (DEG) y rutas
# Datos: GSE60450 (Fu et al., Nat Cell Biol 2015) - glándula mamaria de ratón,
#        células basales y luminales en estado virgen, gestante y lactante.
#
#   anotación -> DGEList -> filterByExpr -> TMM -> MDS -> estimateDisp
#   -> glmQLFit / glmQLFTest -> volcano -> Venn/UpSet -> heatmap -> GO (goana)
#
# Autor: Andoni Moreno Lanceta
# Uso:   Rscript de_analysis.R [GSE60450_Lactation-GenewiseCounts.txt.gz] [resultados_de]
# =============================================================================
suppressPackageStartupMessages({
  library(edgeR)          # incluye limma
  library(org.Mm.eg.db)
  library(ggplot2)
  library(pheatmap)
})

args        <- commandArgs(trailingOnly = TRUE)
counts_file <- if (length(args) >= 1) args[1] else "data/GSE60450_Lactation-GenewiseCounts.txt.gz"
outdir      <- if (length(args) >= 2) args[2] else "resultados_de"
FDR_CUT <- 0.05
LFC_CUT <- 1
fig_dir <- file.path(outdir, "figuras"); tab_dir <- file.path(outdir, "tablas")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

# ----------------------------------------------------------------------------
# 1. Matriz de recuentos y diseño experimental
# ----------------------------------------------------------------------------
raw    <- read.delim(counts_file, row.names = 1, check.names = FALSE)   # filas = ENTREZID
counts <- as.matrix(raw[, colnames(raw) != "Length"])
colnames(counts) <- substr(colnames(counts), 1, 7)                       # MCL1.DG, ...
msg("Genes en la matriz inicial: ", nrow(counts))

# Orden de columnas en GEO: DG DH DI DJ DK DL (basales) LA LB LC LD LE LF (luminales)
samples <- data.frame(
  sample   = colnames(counts),
  celltype = rep(c("basal", "luminal"), each = 6),
  status   = rep(rep(c("virgin", "pregnant", "lactate"), each = 2), 2)
)
group <- factor(paste(samples$celltype, samples$status, sep = "_"))
samples$group <- group
write.table(samples, file.path(tab_dir, "muestras.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ENTREZID -> SYMBOL
symbols <- mapIds(org.Mm.eg.db, keys = rownames(counts), keytype = "ENTREZID",
                  column = "SYMBOL", multiVals = "first")
msg("Genes sin SYMBOL (NA): ", sum(is.na(symbols)))

# Gen con más recuentos en cada muestra
top_gene <- data.frame(
  sample = colnames(counts),
  group  = group,
  gene   = symbols[apply(counts, 2, which.max)],
  counts = apply(counts, 2, max)
)
write.table(top_gene, file.path(tab_dir, "gen_maximo_por_muestra.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
p <- ggplot(top_gene, aes(sample, counts, fill = gene)) + geom_col() +
  geom_text(aes(label = counts), angle = 90, hjust = -0.1, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Gen con más recuentos por muestra", x = "Muestra", y = "Recuentos", fill = "Gen") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "01_gen_maximo_por_muestra.png"), p, width = 7, height = 5, dpi = 300)

# Porcentaje de genes detectados (> 0) por muestra
detected <- data.frame(sample = colnames(counts), pct = 100 * colMeans(counts > 0))
p <- ggplot(detected, aes(sample, pct)) + geom_col(fill = "#4DAF9C") +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.4, size = 3) + ylim(0, 100) +
  labs(title = "Genes con recuentos > 0", x = "Muestra", y = "% de genes") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "02_genes_detectados.png"), p, width = 7, height = 5, dpi = 300)

# ----------------------------------------------------------------------------
# 2. DGEList y filtrado de genes poco expresados
# ----------------------------------------------------------------------------
y <- DGEList(counts = counts, group = group,
             genes = data.frame(ENTREZID = rownames(counts), SYMBOL = symbols))
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]
msg(sprintf("Genes tras filterByExpr: %d (%.2f%%)", nrow(y), 100 * nrow(y) / nrow(counts)))

lab <- paste(group, ave(seq_along(group), group, FUN = seq_along), sep = ".")
box_logcpm <- function(obj, file, title) {
  png(file, width = 2400, height = 1600, res = 300)
  par(mar = c(9, 4, 3, 1))
  boxplot(cpm(obj, log = TRUE), las = 2, col = as.integer(group) + 1, names = lab,
          ylab = "log2 CPM", main = title, cex.axis = 0.7)
  abline(h = median(cpm(obj, log = TRUE)), lty = 2, col = "grey40")
  dev.off()
}
box_logcpm(y, file.path(fig_dir, "03_logcpm_sin_normalizar.png"), "logCPM antes de la normalización TMM")

# ----------------------------------------------------------------------------
# 3. Normalización TMM y exploración (MDS)
# ----------------------------------------------------------------------------
y <- calcNormFactors(y, method = "TMM")
box_logcpm(y, file.path(fig_dir, "04_logcpm_TMM.png"), "logCPM tras la normalización TMM")
write.table(y$samples, file.path(tab_dir, "factores_normalizacion.tsv"), sep = "\t", quote = FALSE)
msg("Muestras con norm.factors < 1: ", paste(rownames(y$samples)[y$samples$norm.factors < 1], collapse = ", "))

cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")
png(file.path(fig_dir, "05_MDS.png"), width = 2000, height = 1800, res = 300)
plotMDS(y, col = cols[as.integer(group)], pch = 16, cex = 1.4, main = "MDS (top 500 genes)")
legend("topleft", legend = levels(group), col = cols, pch = 16, cex = 0.7, bty = "n")
dev.off()

# ----------------------------------------------------------------------------
# 4. Dispersión y modelo GLM quasi-likelihood
# ----------------------------------------------------------------------------
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
y <- estimateDisp(y, design, robust = TRUE)
msg("Dispersión común: ", signif(y$common.dispersion, 4), " (BCV = ", signif(sqrt(y$common.dispersion), 3), ")")
png(file.path(fig_dir, "06_BCV.png"), width = 2000, height = 1800, res = 300); plotBCV(y); dev.off()

fit <- glmQLFit(y, design, robust = TRUE)
png(file.path(fig_dir, "07_QL_dispersion.png"), width = 2000, height = 1800, res = 300); plotQLDisp(fit); dev.off()

genes_interes <- c("Imp4", "Xkr4", "Sox17")
idx <- match(genes_interes, y$genes$SYMBOL)
disp_tab <- data.frame(gene = genes_interes, common = y$common.dispersion,
                       trended = y$trended.dispersion[idx], tagwise = y$tagwise.dispersion[idx])
write.table(disp_tab, file.path(tab_dir, "dispersion_genes_interes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
coef_tab <- data.frame(gene = genes_interes, fit$coefficients[idx, , drop = FALSE], check.names = FALSE)
write.table(coef_tab, file.path(tab_dir, "coeficientes_genes_interes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ----------------------------------------------------------------------------
# 5. Expresión diferencial: lactante vs gestante en cada tipo celular
# ----------------------------------------------------------------------------
contr <- makeContrasts(
  luminal = luminal_lactate - luminal_pregnant,
  basal   = basal_lactate - basal_pregnant,
  levels  = design
)

de <- list()
for (cn in colnames(contr)) {
  qlf <- glmQLFTest(fit, contrast = contr[, cn])
  tt  <- topTags(qlf, n = Inf)$table
  tt$Diff_Expr <- ifelse(tt$FDR <= FDR_CUT & tt$logFC >=  LFC_CUT, "up",
                  ifelse(tt$FDR <= FDR_CUT & tt$logFC <= -LFC_CUT, "down", "ns"))
  de[[cn]] <- tt
  write.table(tt, file.path(tab_dir, paste0("DE_", cn, "_lactate_vs_pregnant.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  msg(sprintf("%s lactate vs pregnant: %d DEG (up %d, down %d) con FDR <= %.2f y |logFC| >= %g",
              cn, sum(tt$Diff_Expr != "ns"), sum(tt$Diff_Expr == "up"), sum(tt$Diff_Expr == "down"),
              FDR_CUT, LFC_CUT))

  top3 <- head(tt[order(-abs(tt$logFC)), ], 3)
  p <- ggplot(tt, aes(logFC, -log10(FDR), colour = Diff_Expr)) +
    geom_point(size = 0.6, alpha = 0.6) +
    scale_colour_manual(values = c(down = "#E41A1C", ns = "grey60", up = "#1F9FBF")) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = 2) +
    geom_hline(yintercept = -log10(FDR_CUT), linetype = 2) +
    geom_text(data = top3, aes(label = SYMBOL), colour = "black", vjust = -0.8, size = 3.2) +
    labs(title = sprintf("Volcano plot: %s lactante vs gestante", cn),
         x = "log2 fold change", y = "-log10(FDR)", colour = "") +
    theme_minimal()
  ggsave(file.path(fig_dir, sprintf("08_volcano_%s.png", cn)), p, width = 7, height = 5.5, dpi = 300)
}

# Genes comunes y únicos entre ambas comparaciones (dirección incluida)
res_mat <- sapply(de, function(tt) {
  v <- setNames(c(up = 1, down = -1, ns = 0)[tt$Diff_Expr], tt$ENTREZID)
  v[rownames(y)]
})
png(file.path(fig_dir, "09_venn_luminal_basal.png"), width = 2000, height = 1800, res = 300)
vennDiagram(res_mat, include = c("up", "down"), counts.col = c("#1F9FBF", "#E41A1C"),
            names = c("Luminal", "Basal"), main = "DEG lactante vs gestante (azul = up, rojo = down)")
dev.off()

if (requireNamespace("UpSetR", quietly = TRUE)) {
  sets <- list(
    Luminal_up   = rownames(res_mat)[res_mat[, "luminal"] ==  1],
    Luminal_down = rownames(res_mat)[res_mat[, "luminal"] == -1],
    Basal_up     = rownames(res_mat)[res_mat[, "basal"]   ==  1],
    Basal_down   = rownames(res_mat)[res_mat[, "basal"]   == -1]
  )
  png(file.path(fig_dir, "10_upset.png"), width = 2400, height = 1600, res = 300)
  print(UpSetR::upset(UpSetR::fromList(sets), order.by = "freq"))
  dev.off()
}

# ----------------------------------------------------------------------------
# 6. Genes sobreexpresados en ambos tipos celulares: heatmap y GO
# ----------------------------------------------------------------------------
shared_up <- rownames(res_mat)[res_mat[, "luminal"] == 1 & res_mat[, "basal"] == 1]
msg("Genes sobreexpresados en luminal y basal: ", length(shared_up))

logcpm <- cpm(y, log = TRUE)
colnames(logcpm) <- lab
mat <- logcpm[shared_up, , drop = FALSE]
sym <- y$genes$SYMBOL[match(shared_up, y$genes$ENTREZID)]
rownames(mat) <- make.unique(ifelse(is.na(sym), shared_up, sym))
pheatmap(mat, scale = "row", show_rownames = FALSE,
         annotation_col = data.frame(row.names = lab, celltype = samples$celltype, status = samples$status),
         main = "Genes sobreexpresados en lactancia (luminal y basal)",
         filename = file.path(fig_dir, "11_heatmap_shared_up.png"), width = 8, height = 7)

go <- goana(shared_up, universe = rownames(y), species = "Mm")
for (ont in c("BP", "MF")) {
  res_go <- topGO(go, ontology = ont, number = 10)
  write.table(res_go, file.path(tab_dir, paste0("GO_", ont, "_shared_up.tsv")), sep = "\t", quote = FALSE)
  msg("Top 10 GO ", ont, ":"); print(res_go[, c("Term", "N", "DE", paste0("P.DE"))])
}

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
msg("Fin. Resultados en ", outdir)
