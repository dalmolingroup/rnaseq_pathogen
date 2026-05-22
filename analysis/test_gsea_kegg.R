# =============================================================================
# Enriquecimento KEGG via GSEA — Gene Set Enrichment Analysis
# -----------------------------------------------------------------------------
# Paradigma alternativo ao ORA usado no relatório (seção 5 com enrichKEGG).
#
# ORA  : pega só os DEGs (lista discreta, com cutoff) e testa hipergeométrico.
# GSEA : pega TODOS os genes ranqueados por log2FoldChange e procura por vias
#        cujos membros se concentram nas extremidades do ranking (cima = UP,
#        baixo = DOWN). Não precisa de cutoff arbitrário e detecta efeitos
#        coordenados pequenos que o ORA descarta.
#
# Para rodar:
#   Rscript analysis/test_gsea_kegg.R
# =============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(KEGGREST)
  library(enrichplot)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(here)
})

setwd(here())

# -----------------------------------------------------------------------------
# 1. Carregar resultados de expressão diferencial
# -----------------------------------------------------------------------------
res_dge <- readRDS("data/DGE/res_dge.rds")
cat("Total de genes em res_dge:", nrow(res_dge), "\n")

# -----------------------------------------------------------------------------
# 2. Reconstruir o dicionário Ensembl Fungi → KEGG (CAALFM_) — mesma lógica
#    da seção 4 do relatório.
# -----------------------------------------------------------------------------
kegg_genes <- keggList("cal")
kegg_mapping <- data.frame(
  kegg_id = sub("cal:", "", names(kegg_genes)),
  stringsAsFactors = FALSE
) %>%
  mutate(
    core       = sub("^CAALFM_", "", kegg_id),
    ensembl_id = gsub("^(C[0-9R]+)(\\d{5})(\\w)(\\w)$", "\\1_\\2\\3_\\4", core)
  )
ensembl_to_kegg <- setNames(kegg_mapping$kegg_id, kegg_mapping$ensembl_id)

# -----------------------------------------------------------------------------
# 3. Construir o vetor ranqueado para GSEA
# -----------------------------------------------------------------------------
# A métrica de ranqueamento define qual sinal a GSEA vai detectar. Usaremos
# o log2FoldChange (padrão para detectar "magnitude" de mudança).
# Alternativas: stat (DESeq2), -log10(padj) * sign(log2FC), etc.

ranks <- res_dge$log2FoldChange
names(ranks) <- rownames(res_dge)
ranks <- ranks[!is.na(ranks)]

# Converter Ensembl Fungi → KEGG ID
ranks_kegg <- ranks[names(ranks) %in% names(ensembl_to_kegg)]
names(ranks_kegg) <- ensembl_to_kegg[names(ranks_kegg)]

# GSEA exige vetor ordenado decrescente e sem duplicatas
ranks_kegg <- sort(ranks_kegg[!duplicated(names(ranks_kegg))], decreasing = TRUE)

cat("Genes no ranking (KEGG IDs únicos):", length(ranks_kegg), "\n")
cat("log2FC — máximo (top do ranking) :", round(max(ranks_kegg), 2), "\n")
cat("log2FC — mínimo (base do ranking):", round(min(ranks_kegg), 2), "\n\n")

# -----------------------------------------------------------------------------
# 4. Rodar GSEA contra o KEGG
# -----------------------------------------------------------------------------
# Argumentos:
#   organism      : "cal" (Candida albicans no KEGG)
#   keyType       : "kegg" (IDs no formato CAALFM_)
#   minGSSize     : 5 (vias com menos de 5 genes anotados são ignoradas)
#   maxGSSize     : 500
#   pvalueCutoff  : 0.5 (cutoff didático, como no relatório)
#   eps           : 0 (permite p-valores muito pequenos via algoritmo adaptativo)

set.seed(42)  # GSEA é estocástica (permutações)
gsea_res <- gseKEGG(
  geneList      = ranks_kegg,
  organism      = "cal",
  keyType       = "kegg",
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 0.5,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = FALSE
)

cat("=== GSEA-KEGG ===\n")
cat("Vias testadas         :", nrow(gsea_res@result), "\n")
cat("Vias significativas   :", nrow(as.data.frame(gsea_res)), "\n\n")

if (nrow(as.data.frame(gsea_res)) > 0) {
  cat("--- Top 10 vias por |NES| ---\n")
  print(
    as.data.frame(gsea_res) %>%
      arrange(desc(abs(NES))) %>%
      dplyr::select(ID, Description, setSize, NES, pvalue, p.adjust) %>%
      head(10),
    row.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# 5. Persistir resultado
# -----------------------------------------------------------------------------
saveRDS(gsea_res, "data/DGE/gsea_kegg.rds")
cat("\nObjeto salvo em data/DGE/gsea_kegg.rds\n\n")

# -----------------------------------------------------------------------------
# 6. Plots
# -----------------------------------------------------------------------------
# Diferente do ORA, GSEA tem visualizações próprias:
#   - dotplot   : termos por p.adjust e setSize, separados por direção (UP/DOWN)
#   - ridgeplot : distribuição dos valores do ranking dentro de cada via
#   - gseaplot2 : curva de enriquecimento clássica (running ES) por via

fig_dir <- "analysis/test_gsea_figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

if (nrow(as.data.frame(gsea_res)) > 0) {

  # 6a. Dotplot — separar UP (NES > 0) e DOWN (NES < 0)
  p_dot <- dotplot(gsea_res, showCategory = 15, split = ".sign",
                   title = "GSEA-KEGG — Vias enriquecidas") +
    facet_grid(. ~ .sign) +
    scale_color_gradient(low = "#E31A1C", high = "#FCBBA1",
                         name = "p.adjust")
  ggsave(file.path(fig_dir, "gsea_dotplot.png"), p_dot,
         width = 11, height = 7, dpi = 150)
  cat("Salvo:", file.path(fig_dir, "gsea_dotplot.png"), "\n")

  # 6b. Ridgeplot — distribuição de log2FC dentro de cada via
  if (requireNamespace("ggridges", quietly = TRUE)) {
    p_ridge <- ridgeplot(gsea_res, showCategory = 15) +
      labs(title = "GSEA-KEGG — Distribuição de log2FC por via",
           x = "log2 Fold Change") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    ggsave(file.path(fig_dir, "gsea_ridgeplot.png"), p_ridge,
           width = 11, height = 8, dpi = 150)
    cat("Salvo:", file.path(fig_dir, "gsea_ridgeplot.png"), "\n")
  } else {
    cat("ggridges não instalado — pulando ridgeplot.\n")
  }

  # 6c. gseaplot2 — curva de enriquecimento das 3 vias mais significativas
  top3 <- as.data.frame(gsea_res) %>%
    arrange(p.adjust) %>%
    head(3) %>%
    pull(ID)
  if (length(top3) > 0) {
    p_curve <- gseaplot2(gsea_res, geneSetID = top3,
                         title = "GSEA-KEGG — Top 3 vias (running ES)")
    ggsave(file.path(fig_dir, "gsea_runES.png"), p_curve,
           width = 10, height = 7, dpi = 150)
    cat("Salvo:", file.path(fig_dir, "gsea_runES.png"), "\n")
  }

  # 6d. Barplot bidirecional — NES como UP/DOWN (similar ao estilo do .qmd)
  df_bidir <- as.data.frame(gsea_res) %>%
    mutate(Direction = ifelse(NES > 0, "Upregulated", "Downregulated"),
           term_lbl  = str_wrap(Description, 45)) %>%
    arrange(desc(abs(NES))) %>%
    head(15)
  p_bar <- ggplot(df_bidir,
                  aes(x = NES, y = reorder(term_lbl, NES), fill = Direction)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("Upregulated"   = "#E31A1C",
                                 "Downregulated" = "#1F78B4")) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    labs(title = "GSEA-KEGG — NES por via (UP vs. DOWN)",
         x = "Normalized Enrichment Score (NES)",
         y = NULL, fill = "Direção") +
    theme_minimal(base_size = 12) +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 9),
          legend.position = "bottom")
  ggsave(file.path(fig_dir, "gsea_barplot_NES.png"), p_bar,
         width = 11, height = 8, dpi = 150)
  cat("Salvo:", file.path(fig_dir, "gsea_barplot_NES.png"), "\n")

} else {
  cat("Nenhuma via enriquecida — sem plots para gerar.\n")
}

# -----------------------------------------------------------------------------
# 7. Comparação com o ORA da seção 5 do relatório (snapshot kegg_*.rds)
# -----------------------------------------------------------------------------
cat("\n=== ORA vs. GSEA — comparação ===\n")
kegg_ora_all <- readRDS("data/DGE/kegg_all.rds")
ora_ids  <- as.data.frame(kegg_ora_all)$ID
gsea_ids <- as.data.frame(gsea_res)$ID
inter    <- intersect(ora_ids, gsea_ids)

cat(sprintf("ORA  (enrichKEGG, snapshot) : %3d vias\n", length(ora_ids)))
cat(sprintf("GSEA (gseKEGG, este script) : %3d vias\n", length(gsea_ids)))
cat(sprintf("Em comum                    : %3d vias  ->  %s\n",
            length(inter),
            if (length(inter) == 0) "(nenhuma)" else paste(inter, collapse = ", ")))
