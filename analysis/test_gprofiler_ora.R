# =============================================================================
# Teste: ORA via g:Profiler usando IDs Ensembl Fungi (C1_00210C_A)
# -----------------------------------------------------------------------------
# Compara o enriquecimento GO obtido pelo g:Profiler com o resultado do CGD-GAF
# já presente no relatório. g:Profiler suporta C. albicans SC5314 via Ensembl
# Fungi e versiona cada release (citável: "g:Profiler version eXXX_egYY").
# -----------------------------------------------------------------------------
# Para rodar:
#   Rscript analysis/test_gprofiler_ora.R
# Requer: gprofiler2 (CRAN), conexão com biit.cs.ut.ee
# =============================================================================

suppressPackageStartupMessages({
  library(gprofiler2)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(httr)
  library(readr)
  library(clusterProfiler)
  library(enrichplot)
  library(GO.db)
})

setwd(here::here())

# -----------------------------------------------------------------------------
# 1. Carregar DEGs no formato Ensembl Fungi (não converter para KEGG)
# -----------------------------------------------------------------------------
res_dge <- readRDS("data/DGE/res_dge.rds")

genes_up   <- rownames(subset(res_dge, signif == "Upregulated"))
genes_down <- rownames(subset(res_dge, signif == "Downregulated"))
genes_all  <- rownames(subset(res_dge, signif != "Not significant"))
universe   <- rownames(res_dge)

cat("UP   :", length(genes_up),   "genes\n")
cat("DOWN :", length(genes_down), "genes\n")
cat("ALL  :", length(genes_all),  "genes\n")
cat("Universo (background):", length(universe), "genes\n\n")

# -----------------------------------------------------------------------------
# 2. Descobrir o código g:Profiler para C. albicans SC5314
# -----------------------------------------------------------------------------
# A lista oficial está em: https://biit.cs.ut.ee/gprofiler/page/organism-list
# Para C. albicans SC5314 o código costuma ser "calbicans" ou
# "calbicans_sc5314". Vamos confirmar via a função de busca.

cat("Procurando organismo C. albicans no g:Profiler...\n")
orgs <- gprofiler2::get_version_info(organism = "calbicans")
cat("Organismo encontrado:\n")
cat("  Nome científico :", orgs$organism, "\n")
cat("  Versão Ensembl  :", orgs$source_versions$ENSG, "\n")
cat("  Fontes p/ GO    :", paste(grep("^GO", names(orgs$source_versions), value = TRUE), collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# 3. Rodar ORA com g:Profiler para UP / DOWN / ALL
# -----------------------------------------------------------------------------
# Usaremos:
#   - sources = c("GO:BP","GO:MF","GO:CC","KEGG","REAC","WP")  apenas o que estiver disponível
#   - correction_method = "g_SCS" (padrão; mais conservador que BH)
#   - significant = FALSE para devolver TODOS os termos (decidimos cutoff depois)

run_gost <- function(query, label) {
  cat("--- g:Profiler:", label, "---\n")
  res <- gost(
    query             = query,
    organism          = "calbicans",
    sources           = c("GO:BP", "KEGG"),
    user_threshold    = 0.5,            # liberal, mesmo cutoff didático do qmd
    correction_method = "fdr",
    significant       = FALSE,          # devolve todos para inspeção
    #custom_bg         = universe        # usa nosso universo (DEG table)
  )
  if (is.null(res)) {
    cat("Nenhum termo retornado.\n\n")
    return(NULL)
  }
  cat("Termos retornados :", nrow(res$result), "\n")
  cat("Significativos    :", sum(res$result$significant), "\n")
  print(head(res$result %>%
               filter(significant) %>%
               dplyr::select(source, term_id, term_name, p_value, intersection_size, term_size),
             10))
  cat("\n")
  return(res)
}

gp_up   <- run_gost(genes_up,   "UP")
gp_down <- run_gost(genes_down, "DOWN")
gp_all  <- run_gost(genes_all,  "ALL")




# -----------------------------------------------------------------------------
# 4. Salvar os resultados para inspeção posterior
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 5. Comparação rápida com o GO-CGD do relatório (apenas BP, para alinhar)
# -----------------------------------------------------------------------------
# Esta seção é informativa, não obrigatória para o teste.
compare_bp <- function(gp_res, label) {
  if (is.null(gp_res)) return(invisible(NULL))
  bp <- gp_res$result %>% dplyr::filter(source == "GO:BP", significant)
  cat(sprintf("[%s] GO:BP enriquecidos no g:Profiler: %d termos\n",
              label, nrow(bp)))
  if (nrow(bp) > 0) {
    print(head(bp %>% dplyr::select(term_id, term_name, p_value, intersection_size), 5))
  }
}
cat("\n=== Sumário GO:BP (g:Profiler) ===\n")
compare_bp(gp_up,   "UP")
compare_bp(gp_down, "DOWN")
compare_bp(gp_all,  "ALL")

# -----------------------------------------------------------------------------
# 6. Enriquecimento KEGG via g:Profiler
# -----------------------------------------------------------------------------
# g:Profiler já consultou a fonte "KEGG" no passo 3 (faz parte do mesmo gost).
# Aqui apenas isolamos os termos KEGG para inspeção dedicada.

show_kegg <- function(gp_res, label) {
  if (is.null(gp_res)) {
    cat(sprintf("[%s] Nenhum resultado g:Profiler.\n\n", label))
    return(invisible(NULL))
  }
  kegg <- gp_res$result %>% dplyr::filter(source == "KEGG")
  cat(sprintf("--- KEGG via g:Profiler: %s ---\n", label))
  cat("Termos KEGG retornados :", nrow(kegg), "\n")
  cat("Significativos         :", sum(kegg$significant), "\n")
  if (nrow(kegg) > 0) {
    top <- kegg %>%
      dplyr::arrange(p_value) %>%
      head(10) %>%
      dplyr::select(term_id, term_name, p_value, significant,
                    intersection_size, term_size)
    print(top)
  }
  cat("\n")
  invisible(kegg)
}

cat("\n=== Enriquecimento KEGG (g:Profiler) ===\n")
kegg_gp_up   <- show_kegg(gp_up,   "UP")
kegg_gp_down <- show_kegg(gp_down, "DOWN")
kegg_gp_all  <- show_kegg(gp_all,  "ALL")

# -----------------------------------------------------------------------------
# 8. Dotplots dos termos GO enriquecidos (g:Profiler)
# -----------------------------------------------------------------------------
# Constrói dotplots no estilo do relatório (clusterProfiler::dotplot), facetando
# por sub-ontologia GO (BP / MF / CC).
#   Eixo X  : GeneRatio = intersection_size / query_size
#   Cor     : -log10(p_value)
#   Tamanho : intersection_size (Count)

wrap_lbl <- function(x, w = 45) stringr::str_wrap(x, width = w)

dotplot_gprofiler_go <- function(gp_res, label, top_n = 15) {
  if (is.null(gp_res)) return(invisible(NULL))
  go <- gp_res$result %>%
    dplyr::filter(grepl("^GO", source), significant)
  if (nrow(go) == 0) {
    message(sprintf("[%s] Nenhum termo GO significativo para plotar.", label))
    return(invisible(NULL))
  }
  go <- go %>%
    dplyr::mutate(
      GeneRatio    = intersection_size / query_size,
      neglog10p    = -log10(p_value),
      term_wrapped = wrap_lbl(term_name)
    ) %>%
    dplyr::arrange(p_value) %>%
    dplyr::group_by(source) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(term_wrapped = factor(term_wrapped,
                  levels = unique(term_wrapped[order(GeneRatio)])))

  size_breaks <- unique(pretty(go$intersection_size, n = 4))
  size_breaks <- size_breaks[size_breaks >= 1 & size_breaks == as.integer(size_breaks)]

  ggplot(go, aes(x = GeneRatio, y = term_wrapped,
                 color = neglog10p, size = intersection_size)) +
    geom_point() +
    facet_wrap(~ source, scales = "free_y", ncol = 1) +
    scale_color_gradient(low = "grey70", high = "#E31A1C",
                         name = expression(-log[10](p))) +
    scale_size_continuous(name = "Count", range = c(3, 8), breaks = size_breaks) +
    labs(title = sprintf("GO Enriquecido — %s (g:Profiler)", label),
         x = "GeneRatio", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          axis.text.y   = element_text(size = 9),
          strip.text    = element_text(face = "bold"),
          legend.position = "right",
          plot.margin   = margin(10, 14, 10, 14))
}

fig_dir <- "analysis/test_gprofiler_figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n=== Dotplots GO (g:Profiler) ===\n")
for (lbl in c("UP", "DOWN", "ALL")) {
  obj <- get(paste0("gp_", tolower(lbl)))
  p   <- dotplot_gprofiler_go(obj, lbl)
  if (!is.null(p)) {
    out <- file.path(fig_dir, sprintf("gprofiler_go_%s.png", tolower(lbl)))
    #ggsave(out, p, width = 10, height = 6, dpi = 150)
    #cat("Salvo:", out, "\n")
  }
}

# -----------------------------------------------------------------------------
# 9. Enriquecimento GO via APIs REST do GO Consortium / EBI
# -----------------------------------------------------------------------------
# Estratégia alternativa ao g:Profiler — consultamos as APIs canônicas e
# montamos a anotação localmente:
#   (a) QuickGO @ EBI : todas as anotações GO BP para taxonId=237561
#       https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch
#   (b) UniProt REST  : mapping UniProt → Ensembl Fungi (mesmo ID de res_dge)
#       https://rest.uniprot.org/uniprotkb/stream
# Em seguida construímos TERM2GENE / TERM2NAME e rodamos clusterProfiler::enricher()
# (mesmo estimador hipergeométrico + BH usado no relatório).

cat("\n=== Enriquecimento GO via API REST (QuickGO + UniProt) ===\n")

# (a) Anotações QuickGO para C. albicans SC5314, aspect = BP
cat("Baixando anotações QuickGO (BP, taxonId=237561)...\n")
resp_qgo <- httr::GET(
  "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch",
  query = list(
    taxonId         = "237561",
    taxonUsage      = "descendants",
    aspect          = "biological_process",
    geneProductType = "protein",
    downloadLimit   = "50000"
  ),
  httr::accept("text/tsv"),
  httr::user_agent("rnaseq_pathogen-minicurso/1.0")
)
stopifnot(httr::status_code(resp_qgo) == 200)
qgo_df <- readr::read_tsv(
  I(httr::content(resp_qgo, as = "text", encoding = "UTF-8")),
  show_col_types = FALSE
)
cat("Anotações QuickGO retornadas :", nrow(qgo_df), "\n")

# (b) Mapping UniProt → Ensembl Fungi (mesma espécie)
cat("Baixando cross-refs UniProt → Ensembl Fungi...\n")
resp_up <- httr::GET(
  "https://rest.uniprot.org/uniprotkb/stream",
  query = list(
    query  = "organism_id:237561",
    fields = "accession,xref_ensemblfungi",
    format = "tsv"
  ),
  httr::user_agent("rnaseq_pathogen-minicurso/1.0")
)
stopifnot(httr::status_code(resp_up) == 200)
up_df <- readr::read_tsv(
  I(httr::content(resp_up, as = "text", encoding = "UTF-8")),
  show_col_types = FALSE
)
# Coluna EnsemblFungi vem como "C1_06680W_A-T;" — limpamos para "C1_06680W_A"
up_df <- up_df %>%
  dplyr::mutate(ensembl_id = sub("-T;.*$", "", EnsemblFungi)) %>%
  dplyr::filter(!is.na(ensembl_id), ensembl_id != "")
cat("UniProt entries com EnsemblFungi:", nrow(up_df), "\n")

# (c) Construir TERM2GENE / TERM2NAME
qgo_mapped <- qgo_df %>%
  dplyr::transmute(
    uniprot = `GENE PRODUCT ID`,
    GO_ID   = `GO TERM`
  ) %>%
  dplyr::left_join(up_df %>% dplyr::select(Entry, ensembl_id),
                   by = c("uniprot" = "Entry")) %>%
  dplyr::filter(!is.na(ensembl_id)) %>%
  dplyr::distinct(GO_ID, ensembl_id)

term2gene_api <- qgo_mapped %>% dplyr::select(GO_ID, ensembl_id)

go_names <- AnnotationDbi::select(
  GO.db,
  keys     = unique(term2gene_api$GO_ID),
  columns  = c("GOID", "TERM", "ONTOLOGY"),
  keytype  = "GOID"
)
go_names <- go_names[!is.na(go_names$TERM), ]
term2name_api <- data.frame(
  GO_ID = go_names$GOID,
  Name  = go_names$TERM
)

cat("Mapeamentos gene-GO (API)  :", nrow(term2gene_api), "\n")
cat("Termos GO únicos (API)     :", length(unique(term2gene_api$GO_ID)), "\n")
cat("Genes únicos com GO (API)  :", length(unique(term2gene_api$ensembl_id)), "\n\n")

# (d) ORA — mesma assinatura que a do relatório (.qmd)
# Observação: o universo via API é menor que o CGD-GAF (cross-ref UniProt
# → EnsemblFungi só está populado para entries revisados), então usamos
# minGSSize=5 para incluir termos GO mais granulares.
run_enricher_api <- function(genes, label) {
  res <- enricher(
    gene = genes,
    TERM2GENE = term2gene_api, TERM2NAME = term2name_api,
    minGSSize = 5, maxGSSize = 500,
    pvalueCutoff = 0.5, qvalueCutoff = 0.5, pAdjustMethod = "BH"
  )
  cat(sprintf("[%s] Termos significativos (p.adj <= 0.5): %d\n",
              label, if (!is.null(res)) nrow(as.data.frame(res)) else 0))
  res
}

cat("=== Resultados GO via API REST ===\n")
go_api_up   <- run_enricher_api(genes_up,   "UP")
go_api_down <- run_enricher_api(genes_down, "DOWN")
go_api_all  <- run_enricher_api(genes_all,  "ALL")

# (d.1) Bloco diagnóstico: top termos por p-valor bruto, mesmo sem passar BH.
# Útil porque o universo mapeável via API é pequeno (~1148 genes) e a correção
# BH sobre muitos testes engole sinais biologicamente coerentes.
cat("\n--- Top termos por p-valor bruto (diagnóstico, ignora BH) ---\n")
diag_top <- function(res, label, n = 10) {
  if (is.null(res) || nrow(res@result) == 0) {
    cat(sprintf("[%s] sem testes válidos.\n", label))
    return(invisible(NULL))
  }
  top <- res@result %>%
    dplyr::arrange(pvalue) %>%
    head(n) %>%
    dplyr::select(ID, Description, pvalue, p.adjust, Count, GeneRatio, BgRatio)
  cat(sprintf("\n[%s] top %d por p-valor bruto:\n", label, n))
  print(top, row.names = FALSE)
}
diag_top(go_api_up,   "UP")
diag_top(go_api_down, "DOWN")
diag_top(go_api_all,  "ALL")

# (e) Persistir resultados para inspeção
saveRDS(go_api_up,    "data/DGE/go_api_up.rds")
saveRDS(go_api_down,  "data/DGE/go_api_down.rds")
saveRDS(go_api_all,   "data/DGE/go_api_all.rds")
saveRDS(term2gene_api,"data/DGE/go_api_term2gene.rds")
saveRDS(term2name_api,"data/DGE/go_api_term2name.rds")

# (f) Dotplots — se houver termos significativos, usa clusterProfiler::dotplot.
#     Caso contrário, plota top-15 por p-valor bruto (diagnóstico).
dotplot_diag_topp <- function(obj, titulo, top_n = 15) {
  if (is.null(obj) || nrow(obj@result) == 0) return(invisible(NULL))
  df <- obj@result %>%
    dplyr::arrange(pvalue) %>%
    head(top_n) %>%
    dplyr::mutate(
      # GeneRatio vem como "k/n" — extrair n (denominador) corretamente
      query_size    = as.numeric(sub(".*/", "", GeneRatio)),
      GeneRatio_num = Count / query_size,
      neglog10p     = -log10(pvalue),
      Description_w = wrap_lbl(Description),
      Description_w = factor(Description_w,
                             levels = unique(Description_w[order(GeneRatio_num)]))
    )
  size_breaks <- unique(pretty(df$Count, n = 4))
  size_breaks <- size_breaks[size_breaks >= 1 & size_breaks == as.integer(size_breaks)]

  ggplot(df, aes(x = GeneRatio_num, y = Description_w,
                 color = neglog10p, size = Count)) +
    geom_point() +
    scale_color_gradient(low = "grey70", high = "#E31A1C",
                         name = expression(-log[10](p[raw]))) +
    scale_size_continuous(name = "Count", range = c(3, 8), breaks = size_breaks) +
    labs(title = titulo,
         subtitle = "Top termos por p bruto\n(nenhum passa BH p.adjust ≤ 0.5)",
         x = "GeneRatio", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
          axis.text.y   = element_text(size = 9),
          plot.margin   = margin(10, 14, 10, 14))
}

cat("\n=== Dotplots GO via API ===\n")
for (lbl in c("UP", "DOWN", "ALL")) {
  obj <- get(paste0("go_api_", tolower(lbl)))
  has_sig <- !is.null(obj) && nrow(as.data.frame(obj)) > 0
  p <- if (has_sig) {
    dotplot(obj, showCategory = 15,
            title = sprintf("GO BP — %s (API REST)", lbl),
            color = "p.adjust") +
      theme(axis.text.y = element_text(size = 9)) +
      scale_color_gradient(low = "#E31A1C", high = "#FCBBA1",
                           name = "p.adjust")
  } else {
    dotplot_diag_topp(obj, sprintf("GO BP — %s (API REST)", lbl))
  }
  if (!is.null(p)) {
    out <- file.path(fig_dir, sprintf("api_go_%s.png", tolower(lbl)))
    ggsave(out, p, width = 11, height = 7, dpi = 150)
    cat("Salvo:", out, "\n")
  } else {
    cat(sprintf("[%s] sem termos sequer testados.\n", lbl))
  }
}

# -----------------------------------------------------------------------------
# 10. Enriquecimento GO via GAF do CGD (gene_association.cgd → submetido ao GO)
# -----------------------------------------------------------------------------
# O arquivo gene_association.cgd é o GAF que o CGD submete periodicamente ao GO
# Consortium. Diferente do cgd_C_albicans_SC5314.gaf (SC5314-only), este contém
# anotações para múltiplos organismos do gênero Candida. Filtramos para SC5314
# (taxon:237561), construímos TERM2GENE/TERM2NAME e rodamos enricher() —
# mesma estatística que o relatório aplica ao GAF SC5314-only.
# Fonte: https://www.candidagenome.org/download/go/
#
# OBSERVAÇÃO: a versão atual do GO.db (3.23+) tem MAIS termos válidos do que
# tinha quando o relatório foi renderizado. Isso multiplica os testes feitos
# pelo enricher e penaliza o BH — fenômeno análogo ao do KEGG na seção 5.
# Por isso o bloco diagnóstico abaixo também mostra top-N por p-valor bruto.

cat("\n=== Enriquecimento GO via GAF do CGD (gene_association.cgd) ===\n")

gaf_raw <- read.delim(
  "data/cgd/gene_association.cgd",
  header       = FALSE,
  comment.char = "!",
  quote        = "",
  sep          = "\t",
  stringsAsFactors = FALSE
)

# Nomear colunas relevantes do GAF 2.0
colnames(gaf_raw)[c(1, 2, 3, 4, 5, 7, 9, 11, 13)] <- c(
  "DB", "DB_Object_ID", "Gene_Symbol", "Qualifier", "GO_ID",
  "Evidence_Code", "Aspect", "Synonyms", "Taxon"
)

cat("Total de anotações no GAF (todos organismos):", nrow(gaf_raw), "\n")

# Filtrar para SC5314
gaf_sc <- gaf_raw %>% dplyr::filter(Taxon == "taxon:237561")
cat("Anotações SC5314 (taxon:237561)              :", nrow(gaf_sc), "\n")

# Extrair systematic name (formato C1_00210C_A) da coluna Synonyms
extract_systematic <- function(synonyms) {
  syns <- unlist(strsplit(synonyms, "\\|"))
  hits <- grep("^C[0-9R]+_\\d{5}[WC]_A$", syns, value = TRUE)
  if (length(hits) > 0) return(hits[1])
  return(NA_character_)
}
gaf_sc$ensembl_id <- sapply(gaf_sc$Synonyms, extract_systematic)

cat("Anotações com systematic name mapeado        :",
    sum(!is.na(gaf_sc$ensembl_id)), "\n")

# Filtrar gaf válido e montar TERM2GENE preliminar
gaf_mapped <- gaf_sc %>%
  dplyr::filter(!is.na(ensembl_id)) %>%
  dplyr::distinct(GO_ID, ensembl_id)

# Pegar nomes + ontologia direto do GO.db (mais autoritativo que o Aspect do GAF)
go_names_gaf <- AnnotationDbi::select(
  GO.db,
  keys     = unique(gaf_mapped$GO_ID),
  columns  = c("GOID", "TERM", "ONTOLOGY"),
  keytype  = "GOID"
)
go_names_gaf <- go_names_gaf[!is.na(go_names_gaf$TERM), ]

# Filtrar para BP usando ONTOLOGY do GO.db (mesma estratégia do relatório)
bp_ids <- go_names_gaf$GOID[go_names_gaf$ONTOLOGY == "BP"]
term2gene_gaf <- gaf_mapped %>% dplyr::filter(GO_ID %in% bp_ids)
term2name_gaf <- data.frame(
  GO_ID = go_names_gaf$GOID[go_names_gaf$ONTOLOGY == "BP"],
  Name  = go_names_gaf$TERM [go_names_gaf$ONTOLOGY == "BP"]
)

cat("Mapeamentos gene-GO (BP)                     :", nrow(term2gene_gaf), "\n")
cat("Termos GO únicos (BP)                        :", length(unique(term2gene_gaf$GO_ID)), "\n")
cat("Genes únicos com GO (BP)                     :", length(unique(term2gene_gaf$ensembl_id)), "\n\n")

# ORA — mesmas configs do relatório (pvalueCutoff/qvalueCutoff = 0.5)
go_gaf_up <- enricher(
  gene = genes_up,
  TERM2GENE = term2gene_gaf, TERM2NAME = term2name_gaf,
  pvalueCutoff = 0.5, qvalueCutoff = 0.5, pAdjustMethod = "BH"
)
go_gaf_down <- enricher(
  gene = genes_down,
  TERM2GENE = term2gene_gaf, TERM2NAME = term2name_gaf,
  pvalueCutoff = 0.5, qvalueCutoff = 0.5, pAdjustMethod = "BH"
)
go_gaf_all <- enricher(
  gene = genes_all,
  TERM2GENE = term2gene_gaf, TERM2NAME = term2name_gaf,
  pvalueCutoff = 0.5, qvalueCutoff = 0.5, pAdjustMethod = "BH"
)

cat("=== Resultados GO via GAF do CGD ===\n")
cat("UP   :", if (!is.null(go_gaf_up))   nrow(as.data.frame(go_gaf_up))   else 0, "termos\n")
cat("DOWN :", if (!is.null(go_gaf_down)) nrow(as.data.frame(go_gaf_down)) else 0, "termos\n")
cat("ALL  :", if (!is.null(go_gaf_all))  nrow(as.data.frame(go_gaf_all))  else 0, "termos\n\n")

# Bloco diagnóstico: top termos por p-valor bruto, mesmo sem passar BH.
cat("--- Top termos por p-valor bruto (diagnóstico, ignora BH) ---\n")
diag_top(go_gaf_up,   "UP")
diag_top(go_gaf_down, "DOWN")
diag_top(go_gaf_all,  "ALL")

# Persistir os objetos
saveRDS(go_gaf_up,    "data/DGE/go_gaf_up.rds")
saveRDS(go_gaf_down,  "data/DGE/go_gaf_down.rds")
saveRDS(go_gaf_all,   "data/DGE/go_gaf_all.rds")

# Comparar: g:Profiler (sec. 3) vs. GAF do CGD (esta seção)
cat("=== Comparação: g:Profiler vs. GAF do CGD ===\n")
compare_gaf_gp <- function(label, gp_res, gaf_obj) {
  gp_ids <- if (!is.null(gp_res)) {
    gp_res$result$term_id[grepl("^GO:", gp_res$result$term_id) &
                          gp_res$result$significant]
  } else character()
  gaf_ids <- if (!is.null(gaf_obj)) as.data.frame(gaf_obj)$ID else character()
  inter <- intersect(gp_ids, gaf_ids)
  cat(sprintf("[%s] g:Profiler sig.: %3d | GAF sig.: %3d | em comum: %3d\n",
              label, length(gp_ids), length(gaf_ids), length(inter)))
}
compare_gaf_gp("UP",   gp_up,   go_gaf_up)
compare_gaf_gp("DOWN", gp_down, go_gaf_down)
compare_gaf_gp("ALL",  gp_all,  go_gaf_all)

# Dotplots — quando há termos significativos usa clusterProfiler::dotplot;
# caso contrário, plota top-15 por p bruto (mesma função do bloco da API).
cat("\n=== Dotplots GO via GAF do CGD ===\n")
for (lbl in c("UP", "DOWN", "ALL")) {
  obj <- get(paste0("go_gaf_", tolower(lbl)))
  has_sig <- !is.null(obj) && nrow(as.data.frame(obj)) > 0
  p <- if (has_sig) {
    dotplot(obj, showCategory = 15,
            title = sprintf("GO BP — %s (GAF do CGD)", lbl),
            color = "p.adjust") +
      theme(axis.text.y = element_text(size = 9)) +
      scale_color_gradient(low = "#E31A1C", high = "#FCBBA1",
                           name = "p.adjust")
  } else {
    dotplot_diag_topp(obj, sprintf("GO BP — %s (GAF do CGD)", lbl))
  }
  if (!is.null(p)) {
    out <- file.path(fig_dir, sprintf("gaf_go_%s.png", tolower(lbl)))
    ggsave(out, p, width = 11, height = 7, dpi = 150)
    cat("Salvo:", out, "\n")
  } else {
    cat(sprintf("[%s] sem termos sequer testados.\n", lbl))
  }
}
