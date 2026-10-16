# ______________________________________________________________________________
# 08_04  Figura 7 - GSEA running-score del eje de receptores nucleares
#
# Que hace: reconstruye el ranking de prostata y re-ejecuta la GSEA con la misma
#           semilla para dibujar el running-score (gseaplot2) del eje de
#           receptores nucleares, con un panel Reactome complementario.
# Entra:    Meta_Analisis_Objetos/contraste_grupos.rds
# Sale:     Meta_Analisis_Graficos/Figura7_running_score_receptores_nucleares.*,
#           Figura7b_*.*
# ______________________________________________________________________________


# 01. SETUP ____________________________________________________________________

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)

select  <- dplyr::select
filter  <- dplyr::filter
arrange <- dplyr::arrange
mutate  <- dplyr::mutate

# Mismos parametros que 08_01 (reproducibilidad exacta).
parametros_e_lista <- list(
  pvalue_cutoff_enriq = 0.05,
  min_gs_size         = 10,
  max_gs_size         = 500,
  semilla_aleatoria   = 1234
)


# 02. RUTAS Y CARGA ____________________________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_objetos_chr  <- file.path(ruta_base_chr, "Meta_Analisis_Objetos")
ruta_graficos_chr <- file.path(ruta_base_chr, "Meta_Analisis_Graficos")
if (!dir.exists(ruta_graficos_chr)) dir.create(ruta_graficos_chr, recursive = TRUE)

obj_contraste <- readRDS(file.path(ruta_objetos_chr, "contraste_grupos.rds"))
firma_prostata_tbl <- obj_contraste$firmas[["prostata_4ds"]]


# 03. RECONSTRUCCION DEL RANKING (identica a construir_ranking de 08_01) ________

# enriquecimiento_completo.rds guarda solo las tablas resumen (as.data.frame),
# no los objetos gseaResult. El running-score (gseaplot2) necesita el objeto
# completo (geneList + geneSets + curva de ES), asi que se reconstruye el ranking
# exacto (Z combinado firmado de prostata) y se re-ejecuta gseGO / gsePathway con
# los mismos parametros y semilla que 08_01, para que NES/FDR coincidan.

mapear_a_entrez <- function(simbolos_chr) {
  suppressWarnings(suppressMessages(
    AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(simbolos_chr),
                          keytype = "SYMBOL", column = "ENTREZID",
                          multiVals = "first")
  )) %>% purrr::discard(is.na)
}

construir_ranking <- function(firma_tbl) {
  ranking_tbl <- firma_tbl %>%
    select(gen_simbolo, z_combinado) %>%
    filter(!is.na(z_combinado), is.finite(z_combinado))
  entrez_chr <- mapear_a_entrez(ranking_tbl$gen_simbolo)
  ranking_tbl <- ranking_tbl %>%
    mutate(entrez = entrez_chr[gen_simbolo]) %>%
    filter(!is.na(entrez)) %>%
    arrange(desc(abs(z_combinado))) %>%
    distinct(entrez, .keep_all = TRUE)
  ranking_num <- ranking_tbl$z_combinado
  names(ranking_num) <- ranking_tbl$entrez
  sort(ranking_num, decreasing = TRUE)
}

set.seed(parametros_e_lista$semilla_aleatoria)
ranking_prostata_num <- construir_ranking(firma_prostata_tbl)
message("Ranking de próstata: ", length(ranking_prostata_num), " genes")


# 04. RE-EJECUCION DEL GSEA ____________________________________________________

set.seed(parametros_e_lista$semilla_aleatoria)
gsea_go_prostata <- gseGO(
  geneList     = ranking_prostata_num,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  minGSSize    = parametros_e_lista$min_gs_size,
  maxGSSize    = parametros_e_lista$max_gs_size,
  pvalueCutoff = parametros_e_lista$pvalue_cutoff_enriq,
  verbose      = FALSE
)

set.seed(parametros_e_lista$semilla_aleatoria)
gsea_react_prostata <- gsePathway(
  geneList     = ranking_prostata_num,
  minGSSize    = parametros_e_lista$min_gs_size,
  maxGSSize    = parametros_e_lista$max_gs_size,
  pvalueCutoff = parametros_e_lista$pvalue_cutoff_enriq,
  verbose      = FALSE
)


# 05. FIGURA 7: RUNNING-SCORE DEL EJE DE RECEPTORES NUCLEARES ___________________

# Cuatro vias GO del eje hormonal, superpuestas en un solo running-score. Todas
# pican hacia el extremo infraexpresado (ES negativo): la firma del programa
# coordinado de receptores nucleares en prostata.
ejes_go_chr <- c(
  "GO:0048545",  # response to steroid hormone                 (NES ~ -2.19)
  "GO:0009755",  # hormone-mediated signaling pathway          (NES ~ -2.05)
  "GO:0141193",  # nuclear receptor-mediated signaling pathway (NES ~ -2.01)
  "GO:0043401"   # steroid hormone receptor signaling pathway  (NES ~ -1.84)
)

# Comprobacion: que las cuatro vias esten en el resultado (si el cutoff las
# hubiera filtrado, se re-ejecuta reteniendo todo).
faltantes_chr <- setdiff(ejes_go_chr, gsea_go_prostata@result$ID)
if (length(faltantes_chr) > 0) {
  message("Reejecutando gseGO con pvalueCutoff = 1 para retener: ",
          paste(faltantes_chr, collapse = ", "))
  set.seed(parametros_e_lista$semilla_aleatoria)
  gsea_go_prostata <- gseGO(
    geneList = ranking_prostata_num, OrgDb = org.Hs.eg.db, ont = "BP",
    minGSSize = parametros_e_lista$min_gs_size,
    maxGSSize = parametros_e_lista$max_gs_size,
    pvalueCutoff = 1, verbose = FALSE
  )
}

etiquetas_ejes_chr <- c(
  "response to steroid hormone",
  "hormone-mediated signaling",
  "nuclear receptor-mediated signaling",
  "steroid hormone receptor signaling"
)

colores_ejes_chr <- c("#8A3B34", "#C1666B", "#4A7C8C", "#E8A33D")

fig7_go <- gseaplot2(
  gsea_go_prostata,
  geneSetID    = ejes_go_chr,
  title        = "Figura 7. Running-score del eje de receptores nucleares (próstata)",
  color        = colores_ejes_chr,
  base_size    = 12,
  rel_heights  = c(1.5, 0.4, 0.9),
  subplots     = 1:3,
  pvalue_table = FALSE,
  ES_geom      = "line"
)

# Leyenda legible con las etiquetas de via (gseaplot2 usa los ID por defecto).
# Se re-etiqueta la leyenda del primer subplot.
fig7_go[[1]] <- fig7_go[[1]] +
  ggplot2::scale_color_manual(values = colores_ejes_chr,
                              labels = etiquetas_ejes_chr, name = NULL) +
  ggplot2::theme(legend.position = c(0.72, 0.30),
                 legend.background = ggplot2::element_rect(fill = "white",
                                                           color = "grey80"),
                 legend.text = ggplot2::element_text(size = 8.5))

# Panel complementario Reactome: "Signaling by Nuclear Receptors" (R-HSA-9006931).
tiene_react <- "R-HSA-9006931" %in% gsea_react_prostata@result$ID
if (tiene_react) {
  fig7_react <- gseaplot2(
    gsea_react_prostata,
    geneSetID    = "R-HSA-9006931",
    title        = "Reactome — Signaling by Nuclear Receptors (próstata)",
    color        = "#4A7C8C",
    base_size    = 12,
    pvalue_table = TRUE,
    ES_geom      = "line"
  )
}


# 06. EXPORTACION ______________________________________________________________

ggplot2::ggsave(
  file.path(ruta_graficos_chr, "Figura7_running_score_receptores_nucleares.png"),
  fig7_go, width = 9, height = 6.5, dpi = 300, bg = "white")
ggplot2::ggsave(
  file.path(ruta_graficos_chr, "Figura7_running_score_receptores_nucleares.pdf"),
  fig7_go, width = 9, height = 6.5, bg = "white")

if (tiene_react) {
  ggplot2::ggsave(
    file.path(ruta_graficos_chr, "Figura7b_running_score_reactome_nuclear_receptors.png"),
    fig7_react, width = 9, height = 6, dpi = 300, bg = "white")
  ggplot2::ggsave(
    file.path(ruta_graficos_chr, "Figura7b_running_score_reactome_nuclear_receptors.pdf"),
    fig7_react, width = 9, height = 6, bg = "white")
}

message("Figura 7 guardada en: ", ruta_graficos_chr)
