# ______________________________________________________________________________
# 08_01  Enriquecimiento funcional de las firmas (ORA y GSEA)
#
# Que hace: calcula el enriquecimiento funcional de las firmas, ORA por firma,
#           direccion y base, y GSEA por grupo sobre el ranking de Z combinado,
#           usando como fondo el universo evaluado en los tres grupos.
# Entra:    Meta_Analisis_Objetos/contraste_grupos.rds
# Sale:     Meta_Analisis_Enriquecimiento/enriquecimiento_ORA.csv,
#           enriquecimiento_GSEA.csv, enriquecimiento_completo.rds
# ______________________________________________________________________________

# __ INDEX ___
# 01. Setup, librerias y parametros
# 02. Rutas y carga de firmas y contraste
# 03. Definicion de las firmas y del universo de fondo
# 04. Mapeo de simbolos a Entrez
# 05. ORA: sobre-representacion por firma, direccion y base
# 06. GSEA: analisis por ranking por grupo y base
# 07. Exportacion de resultados
# 08. Resumen final


# 01. SETUP ____________________________________________________________________

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)      # enriquecimiento Reactome (ORA y GSEA)

# Preferencia de funciones frente al enmascaramiento de Bioconductor.
select   <- dplyr::select
filter   <- dplyr::filter
rename   <- dplyr::rename
mutate   <- dplyr::mutate
arrange  <- dplyr::arrange
map      <- purrr::map
map_dfr  <- purrr::map_dfr
imap_dfr <- purrr::imap_dfr

parametros_e_lista <- list(
  umbral_fdr_firma      = 0.05,   # significancia usada para definir las firmas
  concordancia_minima   = 0.75,
  umbral_logfc_ausencia = 0.30,
  pvalue_cutoff_enriq   = 0.05,   # significancia de los terminos enriquecidos
  qvalue_cutoff_enriq   = 0.10,
  min_gs_size           = 10,
  max_gs_size           = 500,
  semilla_aleatoria     = 1234
)

set.seed(parametros_e_lista$semilla_aleatoria)


# 02. RUTAS Y CARGA ____________________________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_salida_chr   <- file.path(ruta_base_chr, "Meta_Analisis_Resultados")
ruta_objetos_chr  <- file.path(ruta_base_chr, "Meta_Analisis_Objetos")
ruta_enriq_chr    <- file.path(ruta_base_chr, "Meta_Analisis_Enriquecimiento")

if (!dir.exists(ruta_enriq_chr)) {
  dir.create(ruta_enriq_chr, recursive = TRUE)
}

obj_contraste <- readRDS(file.path(ruta_objetos_chr, "contraste_grupos.rds"))

firmas_lista    <- obj_contraste$firmas
contraste_lista <- obj_contraste$contraste

# La firma prostatica principal es la de 4 datasets (comparador normal).
firma_prostata_tbl  <- firmas_lista[["prostata_4ds"]]
firma_testiculo_tbl <- firmas_lista[["testiculo"]]


# 03. DEFINICION DE LAS FIRMAS Y DEL UNIVERSO DE FONDO _________________________

# Firma 1: la interseccion completa prostata y testiculo (287 genes).
firma1_tbl <- contraste_lista[["prostata_4ds"]]

# Firma 2: el subconjunto hormono-especifico (divergente en vejiga).
firma2_tbl <- firma1_tbl %>%
  filter(categoria_vejiga %in%
           c("A_inversion", "B_ausencia", "B_ausencia_debil"))

# Universo de fondo: genes evaluados en los tres grupos. Es la interseccion de
# los espacios genicos, que representa lo que realmente pudo salir significativo.
# Usar el genoma completo como fondo inflaria la significancia de los terminos.
universo_chr <- reduce(
  list(
    firma_prostata_tbl$gen_simbolo,
    firma_testiculo_tbl$gen_simbolo,
    firmas_lista[["vejiga"]]$gen_simbolo
  ),
  intersect
)

message("Tamano del universo de fondo (evaluado en 3 grupos): ",
        length(universo_chr))
message("Firma 1 (compartida): ", nrow(firma1_tbl), " genes")
message("Firma 2 (hormono-especifica): ", nrow(firma2_tbl), " genes")


# 04. MAPEO DE SIMBOLOS A ENTREZ _______________________________________________

# clusterProfiler y ReactomePA trabajan con identificadores Entrez. Se mapea con
# la anotacion local org.Hs.eg.db, que es de version fija (reproducible).

mapear_a_entrez <- function(simbolos_chr) {
  suppressWarnings(suppressMessages(
    AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(simbolos_chr),
                          keytype = "SYMBOL", column = "ENTREZID",
                          multiVals = "first")
  )) %>%
    purrr::discard(is.na)
}

universo_entrez_chr <- mapear_a_entrez(universo_chr) %>% as.character()


# 05. ORA: SOBRE-REPRESENTACION ________________________________________________

# Ejecuta ORA sobre una lista de genes (simbolos) para las tres bases, con el
# universo de fondo comun. Devuelve una tabla larga etiquetada por base.

ora_una_lista <- function(simbolos_chr, etiqueta_chr) {

  entrez_chr <- mapear_a_entrez(simbolos_chr) %>% as.character()

  if (length(entrez_chr) < 5) {
    message("    [", etiqueta_chr, "] muy pocos genes mapeados (",
            length(entrez_chr), "); se omite")
    return(tibble())
  }

  message("    [", etiqueta_chr, "] ", length(entrez_chr), " genes mapeados")

  # GO Biological Process.
  ego <- enrichGO(
    gene          = entrez_chr,
    universe      = universo_entrez_chr,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = parametros_e_lista$pvalue_cutoff_enriq,
    qvalueCutoff  = parametros_e_lista$qvalue_cutoff_enriq,
    minGSSize     = parametros_e_lista$min_gs_size,
    maxGSSize     = parametros_e_lista$max_gs_size,
    readable      = TRUE
  )

  # KEGG.
  ekegg <- enrichKEGG(
    gene          = entrez_chr,
    universe      = universo_entrez_chr,
    organism      = "hsa",
    pvalueCutoff  = parametros_e_lista$pvalue_cutoff_enriq,
    minGSSize     = parametros_e_lista$min_gs_size,
    maxGSSize     = parametros_e_lista$max_gs_size
  )

  # Reactome.
  ereact <- enrichPathway(
    gene          = entrez_chr,
    universe      = universo_entrez_chr,
    pvalueCutoff  = parametros_e_lista$pvalue_cutoff_enriq,
    qvalueCutoff  = parametros_e_lista$qvalue_cutoff_enriq,
    minGSSize     = parametros_e_lista$min_gs_size,
    maxGSSize     = parametros_e_lista$max_gs_size,
    readable      = TRUE
  )

  # Se unifican los tres resultados en una tabla larga.
  a_tibble <- function(resultado_enrich, base_chr) {
    if (is.null(resultado_enrich) || nrow(as.data.frame(resultado_enrich)) == 0) {
      return(tibble())
    }
    as.data.frame(resultado_enrich) %>%
      as_tibble() %>%
      mutate(base = base_chr, lista = etiqueta_chr, .before = 1)
  }

  bind_rows(
    a_tibble(ego,    "GO:BP"),
    a_tibble(ekegg,  "KEGG"),
    a_tibble(ereact, "Reactome")
  )
}


# Se prepara la rejilla de listas: cada firma partida por direccion.

preparar_listas_ora <- function() {

  particiones_lista <- list(
    "Firma1_Up"   = firma1_tbl %>% filter(direccion_hormonal == "Up"),
    "Firma1_Down" = firma1_tbl %>% filter(direccion_hormonal == "Down"),
    "Firma2_Up"   = firma2_tbl %>% filter(direccion_hormonal == "Up"),
    "Firma2_Down" = firma2_tbl %>% filter(direccion_hormonal == "Down")
  )

  particiones_lista %>%
    map(~ .x$gen_simbolo)
}

message("\n=== ORA (sobre-representacion) ===")

listas_ora <- preparar_listas_ora()

ora_tbl <- listas_ora %>%
  imap_dfr(function(simbolos_chr, etiqueta_chr) {
    message("  ", etiqueta_chr, " (", length(simbolos_chr), " genes)")
    ora_una_lista(simbolos_chr, etiqueta_chr)
  })


# 06. GSEA: ANALISIS POR RANKING _______________________________________________

# GSEA no usa listas cortadas sino el ranking completo de genes por su Z
# combinado. Se ejecuta por grupo tumoral (prostata y testiculo) para las tres
# bases. Captura senal distribuida que ORA podria perder con firmas pequenas.

construir_ranking <- function(firma_tbl) {

  # Vector nombrado por Entrez, ordenado de mayor a menor Z combinado.
  ranking_tbl <- firma_tbl %>%
    select(gen_simbolo, z_combinado) %>%
    filter(!is.na(z_combinado), is.finite(z_combinado))

  entrez_chr <- mapear_a_entrez(ranking_tbl$gen_simbolo)

  ranking_tbl <- ranking_tbl %>%
    mutate(entrez = entrez_chr[gen_simbolo]) %>%
    filter(!is.na(entrez)) %>%
    # Si un Entrez aparece por varios simbolos, se conserva el de mayor |Z|.
    arrange(desc(abs(z_combinado))) %>%
    distinct(entrez, .keep_all = TRUE)

  ranking_num <- ranking_tbl$z_combinado
  names(ranking_num) <- ranking_tbl$entrez

  sort(ranking_num, decreasing = TRUE)
}


gsea_un_grupo <- function(firma_tbl, etiqueta_chr) {

  message("  [", etiqueta_chr, "] construyendo ranking")

  ranking_num <- construir_ranking(firma_tbl)

  message("    ", length(ranking_num), " genes en el ranking")

  a_tibble <- function(resultado_gsea, base_chr) {
    if (is.null(resultado_gsea) || nrow(as.data.frame(resultado_gsea)) == 0) {
      return(tibble())
    }
    as.data.frame(resultado_gsea) %>%
      as_tibble() %>%
      mutate(base = base_chr, grupo = etiqueta_chr, .before = 1)
  }

  # GSEA GO:BP.
  gsea_go <- gseGO(
    geneList     = ranking_num,
    OrgDb        = org.Hs.eg.db,
    ont          = "BP",
    minGSSize    = parametros_e_lista$min_gs_size,
    maxGSSize    = parametros_e_lista$max_gs_size,
    pvalueCutoff = parametros_e_lista$pvalue_cutoff_enriq,
    verbose      = FALSE
  )

  # GSEA Reactome.
  gsea_react <- gsePathway(
    geneList     = ranking_num,
    minGSSize    = parametros_e_lista$min_gs_size,
    maxGSSize    = parametros_e_lista$max_gs_size,
    pvalueCutoff = parametros_e_lista$pvalue_cutoff_enriq,
    verbose      = FALSE
  )

  # GSEA KEGG.
  gsea_kegg <- gseKEGG(
    geneList     = ranking_num,
    organism     = "hsa",
    minGSSize    = parametros_e_lista$min_gs_size,
    maxGSSize    = parametros_e_lista$max_gs_size,
    pvalueCutoff = parametros_e_lista$pvalue_cutoff_enriq,
    verbose      = FALSE
  )

  bind_rows(
    a_tibble(gsea_go,    "GO:BP"),
    a_tibble(gsea_kegg,  "KEGG"),
    a_tibble(gsea_react, "Reactome")
  )
}

message("\n=== GSEA (por ranking) ===")

gsea_tbl <- bind_rows(
  gsea_un_grupo(firma_prostata_tbl,  "Prostata"),
  gsea_un_grupo(firma_testiculo_tbl, "Testiculo")
)


# 07. EXPORTACION ______________________________________________________________

write_csv(ora_tbl,  file.path(ruta_enriq_chr, "enriquecimiento_ORA.csv"))
write_csv(gsea_tbl, file.path(ruta_enriq_chr, "enriquecimiento_GSEA.csv"))

saveRDS(list(ora = ora_tbl, gsea = gsea_tbl, universo = universo_chr),
        file.path(ruta_enriq_chr, "enriquecimiento_completo.rds"))


# 08. RESUMEN FINAL ____________________________________________________________

message("\n________________ RESUMEN DEL ENRIQUECIMIENTO ________________\n")

message("Universo de fondo: ", length(universo_chr), " genes")

message("\n--- ORA: terminos significativos por lista y base ---")
if (nrow(ora_tbl) > 0) {
  ora_tbl %>%
    count(lista, base, name = "n_terminos") %>%
    pmap(function(lista, base, n_terminos) {
      message("  ", lista, " | ", base, ": ", n_terminos, " terminos")
    }) %>%
    invisible()
} else {
  message("  Sin terminos significativos en ORA.")
}

message("\n--- GSEA: rutas significativas por grupo y base ---")
if (nrow(gsea_tbl) > 0) {
  gsea_tbl %>%
    count(grupo, base, name = "n_rutas") %>%
    pmap(function(grupo, base, n_rutas) {
      message("  ", grupo, " | ", base, ": ", n_rutas, " rutas")
    }) %>%
    invisible()
} else {
  message("  Sin rutas significativas en GSEA.")
}

# Se destacan las rutas relacionadas con androgenos, por su relevancia directa.
message("\n--- Rutas relacionadas con androgenos (GSEA) ---")
if (nrow(gsea_tbl) > 0) {
  gsea_tbl %>%
    filter(str_detect(tolower(Description), "androgen|hormone|nuclear receptor")) %>%
    select(grupo, base, Description, NES, p.adjust) %>%
    pmap(function(grupo, base, Description, NES, p.adjust) {
      message("  [", grupo, "] ", Description, " (NES=", round(NES, 2),
              ", padj=", signif(p.adjust, 2), ")")
    }) %>%
    invisible()
}

message("\nTablas exportadas en: ", ruta_enriq_chr)
message("____________________________________________________________\n")
