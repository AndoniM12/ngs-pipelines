#!/usr/bin/env Rscript
# =============================================================================
# Metataxonómica 16S: análisis estadístico en R a partir de la salida de QIIME 2
#
#   phyloseq -> control de profundidad -> composición por familia
#   -> diversidad alfa (Shannon) y contraste entre grupos -> diversidad beta
#   (PCoA + PERMANOVA) -> árbol filogenético de un subconjunto taxonómico
#
# Autor: Andoni Moreno Lanceta
# Uso:   Rscript analisis_phyloseq.R <dir_resultados_qiime2> <metadata.tsv> [resultados_R]
# =============================================================================
suppressPackageStartupMessages({
  library(qiime2R)
  library(phyloseq)
  library(microbiome)
  library(vegan)
  library(ggplot2)
})

args    <- commandArgs(trailingOnly = TRUE)
q2dir   <- if (length(args) >= 1) args[1] else "resultados"
meta    <- if (length(args) >= 2) args[2] else "data/sample_metadata.tsv"
outdir  <- if (length(args) >= 3) args[3] else "resultados_R"
MIN_READS <- 1000                 # muestras por debajo se excluyen de los contrastes
GROUP     <- "donor_status"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
save_plot <- function(p, f, w = 9, h = 6) ggsave(file.path(outdir, f), p, width = w, height = h, dpi = 300)

# ----------------------------------------------------------------------------
# 1. Objeto phyloseq
# ----------------------------------------------------------------------------
physeq <- qza_to_phyloseq(
  features = file.path(q2dir, "02_dada2/table.qza"),
  tree     = file.path(q2dir, "04_phylogeny/tree.qza"),
  taxonomy = file.path(q2dir, "03_taxonomy/taxonomy.qza"),
  metadata = meta
)
msg("Muestras: ", nsamples(physeq), " | ASVs: ", ntaxa(physeq))

# Profundidad por muestra: identificar muestras con pocas lecturas tras DADA2
depth <- data.frame(sample = sample_names(physeq), reads = sample_sums(physeq))
depth <- depth[order(depth$reads), ]
write.table(depth, file.path(outdir, "profundidad_por_muestra.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
low <- depth$sample[depth$reads < MIN_READS]
msg("Muestras con < ", MIN_READS, " lecturas (excluidas de los contrastes): ",
    if (length(low)) paste(low, collapse = ", ") else "ninguna")
ps <- prune_samples(sample_sums(physeq) >= MIN_READS, physeq)
# Columna "grupo" con la variable de interés, para no depender de su nombre
sample_data(ps)$grupo <- factor(sample_data(ps)[[GROUP]])

# ----------------------------------------------------------------------------
# 2. Subconjunto taxonómico (ejemplo: orden Clostridiales)
# ----------------------------------------------------------------------------
clos <- subset_taxa(physeq, Order == "Clostridiales")
msg("Taxones totales: ", ntaxa(physeq), " | Clostridiales: ", ntaxa(clos))

# ----------------------------------------------------------------------------
# 3. Composición por familia (abundancia relativa)
# ----------------------------------------------------------------------------
fam <- tax_glom(ps, taxrank = "Family", NArm = FALSE)
fam_rel <- transform_sample_counts(fam, function(x) x / sum(x))
top_fams <- names(sort(taxa_sums(fam_rel), decreasing = TRUE))[1:12]
df <- psmelt(fam_rel)
df$Family <- ifelse(df$OTU %in% top_fams & !is.na(df$Family), as.character(df$Family), "Otras")
p <- ggplot(df, aes(Sample, Abundance, fill = Family)) +
  geom_col(width = 0.95) +
  facet_grid(~ grupo, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Composición por familia (12 más abundantes)", y = "Abundancia relativa", x = NULL) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, size = 5, vjust = 0.5))
save_plot(p, "01_familias_abundancia_relativa.png", w = 12)

# ----------------------------------------------------------------------------
# 4. Diversidad alfa (Shannon) y contraste entre grupos
# ----------------------------------------------------------------------------
alpha_df <- data.frame(sample_data(ps), check.names = FALSE)
alpha_df$shannon <- microbiome::alpha(ps, index = "shannon")$diversity_shannon
write.table(alpha_df, file.path(outdir, "diversidad_alfa.tsv"), sep = "\t", quote = FALSE)
msg("Muestra con menor Shannon: ", rownames(alpha_df)[which.min(alpha_df$shannon)],
    " (", round(min(alpha_df$shannon), 3), ")")
print(aggregate(shannon ~ grupo, data = alpha_df, FUN = mean))

# (a) Contraste por muestra: ignora que hay varias muestras del mismo ratón
w_samples <- wilcox.test(shannon ~ grupo, data = alpha_df)
# (b) Contraste por ratón: una media por animal, unidades independientes
per_mouse <- aggregate(shannon ~ mouse_id + grupo, data = alpha_df, FUN = mean)
w_mouse <- wilcox.test(shannon ~ grupo, data = per_mouse)
msg(sprintf("Wilcoxon Shannon ~ %s | por muestra: p = %.3g | por ratón (n = %d): p = %.3g",
            GROUP, w_samples$p.value, nrow(per_mouse), w_mouse$p.value))

p <- ggplot(alpha_df, aes(grupo, shannon, fill = grupo)) +
  geom_violin(trim = FALSE, alpha = 0.6) + geom_jitter(width = 0.08, size = 1) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "red") +
  labs(title = "Diversidad alfa (Shannon)",
       subtitle = sprintf("Wilcoxon por ratón: p = %.3g", w_mouse$p.value),
       x = GROUP, y = "Índice de Shannon") +
  theme_minimal() + theme(legend.position = "none")
save_plot(p, "02_shannon_por_grupo.png", w = 6, h = 5)

# ----------------------------------------------------------------------------
# 5. Diversidad beta: PCoA (Bray-Curtis y UniFrac ponderado) + PERMANOVA
# ----------------------------------------------------------------------------
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
for (m in c("bray", "wunifrac")) {
  d   <- phyloseq::distance(ps_rel, method = m)
  ord <- ordinate(ps_rel, method = "PCoA", distance = d)
  p <- plot_ordination(ps_rel, ord, color = "grupo", shape = "genotype") +
    geom_point(size = 2.5) + labs(title = paste("PCoA -", m)) + theme_minimal()
  save_plot(p, sprintf("03_PCoA_%s.png", m), w = 7, h = 5.5)
  sd <- data.frame(sample_data(ps_rel), check.names = FALSE)
  # PERMANOVA con donante y genotipo. Las muestras del mismo ratón no son
  # independientes, así que el p-valor es orientativo (ver README).
  perm <- adonis2(d ~ grupo + genotype, data = sd, permutations = 999, by = "margin")
  capture.output(perm, file = file.path(outdir, sprintf("PERMANOVA_%s.txt", m)))
  msg("PERMANOVA (", m, "):"); print(perm)
}

# ----------------------------------------------------------------------------
# 6. Árbol filogenético de la clase Betaproteobacteria
# ----------------------------------------------------------------------------
beta <- subset_taxa(physeq, Class == "Betaproteobacteria")
p <- plot_tree(beta, color = "genotype_and_donor_status", shape = "Family",
               label.tips = "Genus", size = "abundance", plot.margin = 0.5, ladderize = TRUE) +
  labs(title = "Betaproteobacteria: géneros por familia")
save_plot(p, "04_arbol_Betaproteobacteria.png", w = 10, h = 7)

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
msg("Fin. Resultados en ", outdir)
