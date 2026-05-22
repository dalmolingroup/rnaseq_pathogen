# =============================================================================
# Enriquecimento funcional via g:Profiler (gprofiler2::gost)
# -----------------------------------------------------------------------------
# Script didático que faz ORA (Over-Representation Analysis) usando apenas o
# pacote gprofiler2. g:Profiler suporta C. albicans SC5314 via Ensembl Fungi
# (código de organismo: "calbicans"). A versão do banco fica registrada em
# gost_res$meta$version (citável: "g:Profiler version eXXX_egYY").
#
# Para rodar:
#   Rscript analysis/enriquecimento_gprofiler.R
# =============================================================================

suppressPackageStartupMessages({
  library(gprofiler2)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(here)
})

# -----------------------------------------------------------------------------
# 1. Carregar os DEGs
# -----------------------------------------------------------------------------
res_dge <- readRDS(here("data/DGE/res_dge.rds"))

genes_up   <- rownames(subset(res_dge, signif == "Upregulated"))
genes_down <- rownames(subset(res_dge, signif == "Downregulated"))
genes_all  <- rownames(subset(res_dge, signif != "Not significant"))

cat("Genes upregulados  :", length(genes_up),   "\n")
cat("Genes downregulados:", length(genes_down), "\n")
cat("Total DEGs         :", length(genes_all),  "\n\n")

# -----------------------------------------------------------------------------
# 2. Rodar gost() para cada conjunto
# -----------------------------------------------------------------------------
# Argumentos principais:
#   organism          = "calbicans"       # C. albicans SC5314 (Ensembl Fungi)
#   sources           = bancos consultados (GO BP/MF/CC)
#   user_threshold    = 0.5               # cutoff didático (relaxado)
#   correction_method = "fdr"             # Benjamini-Hochberg

gost_up <- gost(
  query             = genes_up,
  organism          = "calbicans",
  sources           = c("GO:BP"),
  user_threshold    = 0.5,
  correction_method = "fdr"
)

gost_down <- gost(
  query             = genes_down,
  organism          = "calbicans",
  sources           = c("GO:BP"),
  user_threshold    = 0.5,
  correction_method = "fdr"
)

gost_all <- gost(
  query             = genes_all,
  organism          = "calbicans",
  sources           = c("GO:BP"),
  user_threshold    = 0.5,
  correction_method = "fdr"
)

# Versão do banco usado (para citar nos métodos)
cat("g:Profiler version:", gost_all$meta$version, "\n\n")

cat("=== Termos significativos (FDR <= 0.5) ===\n")
cat("UP   :", if (is.null(gost_up))   0 else nrow(gost_up$result),   "\n")
cat("DOWN :", if (is.null(gost_down)) 0 else nrow(gost_down$result), "\n")
cat("ALL  :", if (is.null(gost_all))  0 else nrow(gost_all$result),  "\n\n")

# -----------------------------------------------------------------------------
# 3. Diretório de saída
# -----------------------------------------------------------------------------
fig_dir <- here("analysis", "gprofiler_figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Salvar os objetos para reuso (não precisa rodar a API outra vez)
#saveRDS(gost_up,   here("data/DGE/gost_up.rds"))
#saveRDS(gost_down, here("data/DGE/gost_down.rds"))
#saveRDS(gost_all,  here("data/DGE/gost_all.rds"))

# -----------------------------------------------------------------------------
# 4. Manhattan plot (nativo do gprofiler2)
# -----------------------------------------------------------------------------
# O gostplot mostra -log10(p) por fonte. Útil para visão geral de quais
# bancos (GO BP, GO MF, GO CC) contribuíram com os enriquecimentos.

if (!is.null(gost_up)) {
  p1 <- gostplot(gost_up, interactive = FALSE)
  ggsave(file.path(fig_dir, "manhattan_up.png"), p1,
         width = 10, height = 6, dpi = 150)
  cat("Salvo:", file.path(fig_dir, "manhattan_up.png"), "\n")
}

if (!is.null(gost_down)) {
  p2 <- gostplot(gost_down, interactive = FALSE)
  ggsave(file.path(fig_dir, "manhattan_down.png"), p2,
         width = 10, height = 6, dpi = 150)
  cat("Salvo:", file.path(fig_dir, "manhattan_down.png"), "\n")
}

if (!is.null(gost_all)) {
  p3 <- gostplot(gost_all, interactive = FALSE)
  ggsave(file.path(fig_dir, "manhattan_all.png"), p3,
         width = 10, height = 6, dpi = 150)
  cat("Salvo:", file.path(fig_dir, "manhattan_all.png"), "\n")
}

# -----------------------------------------------------------------------------
# 5. Dotplots — estilo do relatório
# -----------------------------------------------------------------------------
# Convenções:
#   Eixo X  : GeneRatio = intersection_size / query_size
#   Cor     : -log10(p_value)
#   Tamanho : intersection_size (Count)
#   Facets  : sub-ontologia GO (BP / MF / CC)

dotplot_gost <- function(gost_res, titulo, top_n_por_fonte = 5) {
  df <- gost_res$result %>%
    mutate(
      GeneRatio = intersection_size / query_size,
      neglog10p = -log10(p_value),
      term_lbl  = str_wrap(term_name, 45)
    ) %>%
    group_by(source) %>%
    arrange(p_value, .by_group = TRUE) %>%
    slice_head(n = top_n_por_fonte) %>%
    ungroup() %>%
    # term_id é único; usamos como chave de y e exibimos term_lbl
    mutate(term_id = factor(term_id, levels = unique(term_id[order(GeneRatio)])))

  ggplot(df, aes(x = GeneRatio, y = term_id,
                 color = neglog10p, size = intersection_size)) +
    geom_point() +
    facet_wrap(~ source, scales = "free_y", ncol = 1) +
    scale_y_discrete(labels = setNames(df$term_lbl, df$term_id)) +
    scale_color_gradient(low = "grey70", high = "#E31A1C",
                         name = expression(-log[10](p))) +
    scale_size_continuous(name = "Count", range = c(3, 8)) +
    labs(title = titulo, x = "GeneRatio", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title  = element_text(face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 9),
          strip.text  = element_text(face = "bold"))
}

if (!is.null(gost_up) && nrow(gost_up$result) > 0) {
  p <- dotplot_gost(gost_up, "g:Profiler — Genes Upregulados")
  ggsave(file.path(fig_dir, "dotplot_up.png"), p,
         width = 10, height = 7, dpi = 150)
  cat("Salvo:", file.path(fig_dir, "dotplot_up.png"), "\n")
} else {
  cat("Nenhum termo enriquecido para genes upregulados.\n")
}

if (!is.null(gost_down) && nrow(gost_down$result) > 0) {
  p <- dotplot_gost(gost_down, "g:Profiler — Genes Downregulados")
  ggsave(file.path(fig_dir, "dotplot_down.png"), p,
         width = 10, height = 7, dpi = 150)
  cat("Salvo:", file.path(fig_dir, "dotplot_down.png"), "\n")
} else {
  cat("Nenhum termo enriquecido para genes downregulados.\n")
}

if (!is.null(gost_all) && nrow(gost_all$result) > 0) {
  p <- dotplot_gost(gost_all, "g:Profiler — Todos os DEGs")
  ggsave(file.path(fig_dir, "dotplot_all.png"), p,
         width = 10, height = 7, dpi = 150)
  cat("Salvo:", file.path(fig_dir, "dotplot_all.png"), "\n")
} else {
  cat("Nenhum termo enriquecido para todos os DEGs.\n")
}

# -----------------------------------------------------------------------------
# 6. Barplot bidirecional — UP vs. DOWN
# -----------------------------------------------------------------------------
# Barras UP (vermelho) crescem para a direita; DOWN (azul) crescem para a
# esquerda. Ambas compartilham o mesmo eixo y (termos GO).

df_up <- if (!is.null(gost_up) && nrow(gost_up$result) > 0) {
  gost_up$result %>%
    mutate(Direction = "Upregulated", Count_dir = intersection_size)
} else { data.frame() }

df_down <- if (!is.null(gost_down) && nrow(gost_down$result) > 0) {
  gost_down$result %>%
    mutate(Direction = "Downregulated", Count_dir = -intersection_size)
} else { data.frame() }

if (nrow(df_up) > 0 || nrow(df_down) > 0) {
  df_bidir <- bind_rows(
    df_up   %>% arrange(p_value) %>% head(10),
    df_down %>% arrange(p_value) %>% head(10)
  ) %>%
    mutate(term_lbl = str_wrap(term_name, 45))

  p_bidir <- ggplot(df_bidir,
                    aes(x = Count_dir,
                        y = reorder(term_lbl, Count_dir),
                        fill = Direction)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("Upregulated"   = "#E31A1C",
                                 "Downregulated" = "#1F78B4")) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    labs(title = "g:Profiler — UP vs. DOWN",
         x = "Número de Genes", y = NULL, fill = "Regulação") +
    theme_minimal(base_size = 12) +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 9),
          legend.position = "bottom")

  ggsave(file.path(fig_dir, "barplot_bidirecional.png"), p_bidir,
         width = 11, height = 8, dpi = 150)
  cat("Salvo:", file.path(fig_dir, "barplot_bidirecional.png"), "\n")
}

# -----------------------------------------------------------------------------
# 7. Barplot — Todos os DEGs (roxo)
# -----------------------------------------------------------------------------

if (!is.null(gost_all) && nrow(gost_all$result) > 0) {
  df_all <- gost_all$result %>%
    arrange(p_value) %>%
    head(15) %>%
    mutate(term_lbl = str_wrap(term_name, 45))

  p_all <- ggplot(df_all,
                  aes(x = intersection_size,
                      y = reorder(term_lbl, intersection_size))) +
    geom_col(fill = "#6A3D9A", width = 0.7) +
    labs(title = "g:Profiler — Todos os DEGs",
         x = "Número de Genes", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 9))

  ggsave(file.path(fig_dir, "barplot_all.png"), p_all,
         width = 10, height = 7, dpi = 150)
  cat("Salvo:", file.path(fig_dir, "barplot_all.png"), "\n")
}
