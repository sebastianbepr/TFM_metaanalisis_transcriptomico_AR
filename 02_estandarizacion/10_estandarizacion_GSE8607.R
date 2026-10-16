# ______________________________________________________________________________
# 02_10  Estandarizacion de GSE8607 (testiculo)
#
# Que hace: pasa el ExpressionSet de GSE8607 (microarray Affymetrix GPL8300) al
#           esquema estandar del proyecto (simbolo HGNC), aplica log2(x + 1)
#           porque viene en escala lineal, y lo guarda.
# Entra:    Testiculo/GSE8607*ExpressionSet*.rds
# Sale:     Datasets_Estandarizados/GSE8607_estandarizado.rds
# ______________________________________________________________________________

# __ INDEX ___
# 01. LIBRERIAS
# 02. RUTAS
# 03. IMPORTAR EL EXPRESSIONSET
# 04. EXTRAER EXPRESION, ANOTACION Y FENOTIPO
# 05. ASIGNAR GRUPO DIAGNOSTICO DESDE TEXTO LIBRE
# 06. CONTROL DE CALIDAD Y TRANSFORMACION LOG2
# 07. ANOTACION SONDA A SIMBOLO
# 08. COLAPSAR SONDAS A SIMBOLOS DE GEN
# 09. ENSAMBLAR LA MATRIZ ESTANDARIZADA
# 10. DESCRIPTOR DEL DATASET
# 11. EXPORTAR EL RDS ESTANDARIZADO
# 12. RESUMEN
# ________________________________________________________________________


# 01. LIBRERIAS ________________________________________________________________

library(tidyverse)
library(Biobase)


# 02. RUTAS ____________________________________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
base_datos_directorio_path <- ruta   # el resto del script usa este nombre sin cambios

testiculo_directorio_path <-
  file.path(base_datos_directorio_path, "Testiculo")

datasets_estandarizados_directorio_path <-
  file.path(base_datos_directorio_path, "Datasets_Estandarizados")

if (!dir.exists(datasets_estandarizados_directorio_path)) {
  dir.create(datasets_estandarizados_directorio_path, recursive = TRUE)
}

gse_objetivo_id_char  <- "GSE8607"
organo_objetivo_char  <- "Testis"
plataforma_tipo_char  <- "Microarray (Affymetrix HG-U95Av2, GPL8300)"


# 03. IMPORTAR EL EXPRESSIONSET ________________________________________________

message("Auditing dataset ", gse_objetivo_id_char)

testiculo_eset_archivo_path <-
  list.files(
    path       = testiculo_directorio_path,
    pattern    = paste0(gse_objetivo_id_char, ".*ExpressionSet.*\\.rds$"),
    recursive  = TRUE,
    full.names = TRUE
  )

if (length(testiculo_eset_archivo_path) == 0) {
  stop("Critical error: ExpressionSet not found for ",
       gse_objetivo_id_char)
}

testiculo_eset <- readRDS(testiculo_eset_archivo_path[1])

message("ExpressionSet imported: ",
        nrow(testiculo_eset), " probes x ",
        ncol(testiculo_eset), " samples")


# 04. EXTRAER EXPRESION, ANOTACION Y FENOTIPO __________________________________

testiculo_exprs_matriz_inicial <- exprs(testiculo_eset)

testiculo_feature_tbl_inicial <-
  fData(testiculo_eset) %>%
  rownames_to_column("probe_id") %>%
  as_tibble()

testiculo_pheno_tbl_inicial <-
  pData(testiculo_eset) %>%
  rownames_to_column("gsm_accession") %>%
  as_tibble()

message("Expression matrix: ",
        nrow(testiculo_exprs_matriz_inicial), " probes x ",
        ncol(testiculo_exprs_matriz_inicial), " samples")


# 05. ASIGNAR GRUPO DIAGNOSTICO DESDE TEXTO LIBRE ______________________________

# La etiqueta diagnostica no viene en un campo categorico propio, sino dentro
# del texto descriptivo de source_name_ch1. Los seminomas son el grupo Tumor y
# los testiculos con espermatogenesis normal (de orquiectomia por carcinoma de
# prostata o post-vasectomia) son el grupo Normal. Guardamos la procedencia
# clinica de cada control en el campo de subtipo.

testiculo_pheno_tbl_curado <-
  testiculo_pheno_tbl_inicial %>%
  mutate(
    grupo_diagnostico = case_when(
      str_detect(source_name_ch1,
                 regex("seminoma", ignore_case = TRUE))
      ~ "Tumor",
      str_detect(source_name_ch1,
                 regex("regular spermatogenesis|normal|control|vasectomy",
                       ignore_case = TRUE))
      ~ "Normal",
      TRUE
      ~ NA_character_
    ),
    tejido_control_subtipo = case_when(
      str_detect(source_name_ch1,
                 regex("orchiectomy.*prostate", ignore_case = TRUE))
      ~ "Normal testis (orchiectomy for prostate carcinoma)",
      str_detect(source_name_ch1,
                 regex("vasectomy", ignore_case = TRUE))
      ~ "Normal testis (post-vasectomy)",
      TRUE
      ~ NA_character_
    ),
    histologia_clase = if_else(
      grupo_diagnostico == "Tumor", "Seminoma", NA_character_
    ),
    estadio_tumoral = str_extract(source_name_ch1, "pT\\d[^,]*"),
    edad_paciente = str_extract(source_name_ch1, "(?<=age:\\s)\\d+")
  )

message("Standardized diagnostic distribution:")
print(table(testiculo_pheno_tbl_curado$grupo_diagnostico,
            useNA = "ifany"))

message("Control subtype distribution:")
print(table(testiculo_pheno_tbl_curado$tejido_control_subtipo,
            useNA = "ifany"))


# 06. CONTROL DE CALIDAD Y TRANSFORMACION LOG2 _________________________________

# Las intensidades vienen en escala lineal, al reves que los valores log2 RMA de
# GSE3218. Aplicamos log2(x + 1) para dejarlas en la escala que limma espera para
# el modelado de expresion diferencial.

message("Raw value range (linear scale):")
message("  min = ",
        round(min(testiculo_exprs_matriz_inicial, na.rm = TRUE), 2))
message("  max = ",
        round(max(testiculo_exprs_matriz_inicial, na.rm = TRUE), 2))
message("  median = ",
        round(median(testiculo_exprs_matriz_inicial, na.rm = TRUE), 2))

testiculo_escala_lineal_lgl <-
  max(testiculo_exprs_matriz_inicial, na.rm = TRUE) > 100

# Ponemos a cero los valores negativos o nulos antes del logaritmo y luego
# aplicamos log2(x + 1).

testiculo_exprs_matriz_log2 <-
  testiculo_exprs_matriz_inicial %>%
  pmax(0) %>%
  {log2(. + 1)}

message("Applied log2(x + 1) transformation: ", testiculo_escala_lineal_lgl)
message("Transformed value range:")
message("  min = ",
        round(min(testiculo_exprs_matriz_log2, na.rm = TRUE), 2))
message("  max = ",
        round(max(testiculo_exprs_matriz_log2, na.rm = TRUE), 2))
message("  median = ",
        round(median(testiculo_exprs_matriz_log2, na.rm = TRUE), 2))

# Comprobamos que el receptor de androgenos (AR) esta en el array.

testiculo_ar_presente_lgl <-
  any(testiculo_feature_tbl_inicial$`Gene symbol` == "AR", na.rm = TRUE)

message("Androgen Receptor (AR) represented on array: ",
        testiculo_ar_presente_lgl)


# 07. ANOTACION SONDA A SIMBOLO ________________________________________________

# Pasamos las sondas a simbolo HGNC y descartamos las que apuntan a varios genes
# a la vez (GENE1///GENE2) para no atribuir la senal de forma ambigua. Misma
# regla que en GSE3218.

testiculo_feature_tbl_curado <-
  testiculo_feature_tbl_inicial %>%
  select(probe_id, gene_symbol = `Gene symbol`) %>%
  mutate(gene_symbol = na_if(str_trim(gene_symbol), "")) %>%
  filter(!is.na(gene_symbol)) %>%
  filter(!str_detect(gene_symbol, fixed("///")))

testiculo_probes_multigen_int <-
  testiculo_feature_tbl_inicial %>%
  filter(str_detect(`Gene symbol`, fixed("///"))) %>%
  nrow()

message("Probes with a unique gene symbol: ",
        nrow(testiculo_feature_tbl_curado), " / ",
        nrow(testiculo_feature_tbl_inicial))
message("Multi-gene probes discarded: ", testiculo_probes_multigen_int)


# 08. COLAPSAR SONDAS A SIMBOLOS DE GEN ________________________________________

# Colapsamos a un registro por gen quedandonos con la de mayor expresion media,
# trabajando ya sobre la matriz transformada a log2.

testiculo_exprs_tbl_inicial <-
  testiculo_exprs_matriz_log2 %>%
  as.data.frame() %>%
  rownames_to_column("probe_id") %>%
  as_tibble() %>%
  inner_join(testiculo_feature_tbl_curado, by = "probe_id")

testiculo_exprs_tbl_curado <-
  testiculo_exprs_tbl_inicial %>%
  mutate(
    expresion_media = rowMeans(
      select(., -probe_id, -gene_symbol),
      na.rm = TRUE
    )
  ) %>%
  arrange(gene_symbol, desc(expresion_media)) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  select(-expresion_media, -probe_id)

message("Genes retained after probe collapse: ",
        nrow(testiculo_exprs_tbl_curado))


# 09. ENSAMBLAR LA MATRIZ ESTANDARIZADA ________________________________________

testiculo_genes_anotacion_tbl_final <-
  tibble(
    gene_id_native = testiculo_exprs_tbl_curado$gene_symbol,
    gene_symbol    = testiculo_exprs_tbl_curado$gene_symbol,
    ensembl_id     = NA_character_
  )

testiculo_expresion_matriz_final <-
  testiculo_exprs_tbl_curado %>%
  select(-gene_symbol) %>%
  as.matrix()

rownames(testiculo_expresion_matriz_final) <-
  testiculo_genes_anotacion_tbl_final$gene_symbol

testiculo_metadatos_tbl_final <-
  testiculo_pheno_tbl_curado %>%
  select(sample_id = gsm_accession,
         gsm_accession,
         grupo_diagnostico,
         tejido_control_subtipo,
         histologia_clase,
         estadio_tumoral,
         edad_paciente) %>%
  slice(match(colnames(testiculo_expresion_matriz_final), sample_id))

message("Standardized matrix assembled: ",
        nrow(testiculo_expresion_matriz_final), " genes x ",
        ncol(testiculo_expresion_matriz_final), " samples")


# 10. DESCRIPTOR DEL DATASET ___________________________________________________

testiculo_descriptor_lista <-
  list(
    gse_id                  = gse_objetivo_id_char,
    organo                  = organo_objetivo_char,
    plataforma_tipo         = plataforma_tipo_char,
    gene_id_convention      = "HGNC gene symbol",
    matrix_value_nature     = "log2(x + 1) transformed microarray intensities (originally linear)",
    dge_branch_recommended  = "limma on log2 intensities",
    n_total                 = ncol(testiculo_expresion_matriz_final),
    n_tumor                 = sum(testiculo_metadatos_tbl_final$grupo_diagnostico
                                  == "Tumor", na.rm = TRUE),
    n_normal                = sum(testiculo_metadatos_tbl_final$grupo_diagnostico
                                  == "Normal", na.rm = TRUE),
    n_sin_etiqueta          = sum(is.na(testiculo_metadatos_tbl_final$grupo_diagnostico)),
    comparator_subtype      = "Normal testis (orchiectomy / post-vasectomy)",
    histology               = "Pure seminoma",
    paired_design           = FALSE,
    log2_transform_applied  = TRUE,
    multi_gene_probes_removed = testiculo_probes_multigen_int,
    feature_filtering       = "Unique-symbol probes only; multi-gene probes discarded",
    duplicate_collapse      = "Maximum mean expression per gene symbol",
    power_limitation        = "Only 3 normal controls; limited statistical power",
    fecha_estandarizacion   = format(Sys.Date(), "%Y-%m-%d")
  )


# 11. EXPORTAR EL RDS ESTANDARIZADO ____________________________________________

testiculo_dataset_estandarizado_lista <-
  list(
    conteos_matriz         = testiculo_expresion_matriz_final,
    genes_anotacion_tbl    = testiculo_genes_anotacion_tbl_final,
    muestras_metadatos_tbl = testiculo_metadatos_tbl_final,
    dataset_descriptor     = testiculo_descriptor_lista
  )

testiculo_rds_destino_path <-
  file.path(
    datasets_estandarizados_directorio_path,
    paste0(gse_objetivo_id_char, "_estandarizado.rds")
  )

saveRDS(testiculo_dataset_estandarizado_lista, testiculo_rds_destino_path)

message("Standardized object exported to: ",
        basename(testiculo_rds_destino_path))


# 12. RESUMEN __________________________________________________________________

message("")
message("Pipeline completed for ", gse_objetivo_id_char)
message("Organ: ",
        testiculo_descriptor_lista$organo)
message("Platform: ",
        testiculo_descriptor_lista$plataforma_tipo)
message("Matrix nature: ",
        testiculo_descriptor_lista$matrix_value_nature)
message("DGE branch recommended: ",
        testiculo_descriptor_lista$dge_branch_recommended)
message("Histology: ",
        testiculo_descriptor_lista$histology)
message("Comparator subtype: ",
        testiculo_descriptor_lista$comparator_subtype)
message("Features retained: ",
        nrow(testiculo_dataset_estandarizado_lista$conteos_matriz),
        " unique gene symbols")
message("Samples retained: ",
        testiculo_descriptor_lista$n_total,
        " (Tumor = ", testiculo_descriptor_lista$n_tumor,
        ", Normal = ", testiculo_descriptor_lista$n_normal,
        if (testiculo_descriptor_lista$n_sin_etiqueta > 0)
          paste0(", Unlabeled = ",
                 testiculo_descriptor_lista$n_sin_etiqueta)
        else "",
        ")")
message("Note: ", testiculo_descriptor_lista$power_limitation)
message("Output: ", testiculo_rds_destino_path)
