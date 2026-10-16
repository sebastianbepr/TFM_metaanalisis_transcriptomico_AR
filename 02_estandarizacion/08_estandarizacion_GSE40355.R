# ______________________________________________________________________________
# 02_08  Estandarizacion de GSE40355 (vejiga)
#
# Que hace: pasa el ExpressionSet de mRNA de GSE40355 (microarray Agilent
#           GPL13497) al esquema estandar del proyecto (simbolo HGNC) y lo guarda.
# Entra:    Vejiga/GSE40355*GPL13497*ExpressionSet*.rds
# Sale:     Datasets_Estandarizados/GSE40355_estandarizado.rds
# ______________________________________________________________________________

# __ INDEX ___
# 01. LIBRERIAS
# 02. RUTAS
# 03. IMPORTAR EL EXPRESSIONSET (PLATAFORMA mRNA)
# 04. EXTRAER EXPRESION, ANOTACION Y FENOTIPO
# 05. ASIGNAR GRUPO DIAGNOSTICO Y GRADO TUMORAL
# 06. CONTROL DE CALIDAD DE LAS SONDAS
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

vejiga_metadatos_directorio_path <-
  file.path(base_datos_directorio_path, "Vejiga")

datasets_estandarizados_directorio_path <-
  file.path(base_datos_directorio_path, "Datasets_Estandarizados")

if (!dir.exists(datasets_estandarizados_directorio_path)) {
  dir.create(datasets_estandarizados_directorio_path, recursive = TRUE)
}

gse_objetivo_id_char  <- "GSE40355"
organo_objetivo_char  <- "Bladder"
plataforma_tipo_char  <- "Microarray (Agilent mRNA, GPL13497)"


# 03. IMPORTAR EL EXPRESSIONSET (PLATAFORMA mRNA) ______________________________

message("Auditing dataset ", gse_objetivo_id_char, " (mRNA platform only)")

# Nos quedamos solo con el ExpressionSet de mRNA; el companero de miRNA
# (GPL8227) se deja fuera a proposito porque el metaanalisis mira genes
# codificantes.

vejiga_eset_archivo_path <-
  list.files(
    path       = vejiga_metadatos_directorio_path,
    pattern    = paste0(gse_objetivo_id_char,
                        ".*GPL13497.*ExpressionSet.*\\.rds$"),
    recursive  = TRUE,
    full.names = TRUE
  )

if (length(vejiga_eset_archivo_path) == 0) {
  stop("Critical error: mRNA ExpressionSet not found for ",
       gse_objetivo_id_char)
}

vejiga_eset <- readRDS(vejiga_eset_archivo_path[1])

message("ExpressionSet imported: ",
        nrow(vejiga_eset), " probes x ",
        ncol(vejiga_eset), " samples")
message("Annotation platform: ", annotation(vejiga_eset))


# 04. EXTRAER EXPRESION, ANOTACION Y FENOTIPO __________________________________

vejiga_exprs_matriz_inicial <- exprs(vejiga_eset)

vejiga_feature_tbl_inicial <-
  fData(vejiga_eset) %>%
  rownames_to_column("probe_id") %>%
  as_tibble()

vejiga_pheno_tbl_inicial <-
  pData(vejiga_eset) %>%
  rownames_to_column("gsm_accession") %>%
  as_tibble()

message("Expression matrix: ",
        nrow(vejiga_exprs_matriz_inicial), " probes x ",
        ncol(vejiga_exprs_matriz_inicial), " samples")


# 05. ASIGNAR GRUPO DIAGNOSTICO Y GRADO TUMORAL ________________________________

# La etiqueta diagnostica viene en 'characteristics_ch1' como una de tres
# categorias de tejido. Los dos grados tumorales se juntan en un unico grupo
# Tumor para el contraste principal; guardamos el grado original en
# 'tumor_grade' por si hace falta un analisis por estratos.

vejiga_pheno_tbl_curado <-
  vejiga_pheno_tbl_inicial %>%
  mutate(
    grupo_diagnostico = case_when(
      str_detect(characteristics_ch1,
                 regex("urothelial carcinoma|cancer|tumor|carcinoma",
                       ignore_case = TRUE))
      ~ "Tumor",
      str_detect(characteristics_ch1,
                 regex("nonmaligant|nonmalignant|normal|benign",
                       ignore_case = TRUE))
      ~ "Normal",
      TRUE
      ~ NA_character_
    ),
    tumor_grade = case_when(
      str_detect(characteristics_ch1,
                 regex("high-grade", ignore_case = TRUE))
      ~ "High-grade",
      str_detect(characteristics_ch1,
                 regex("low-grade", ignore_case = TRUE))
      ~ "Low-grade",
      TRUE
      ~ NA_character_
    ),
    tejido_control_subtipo = case_when(
      grupo_diagnostico == "Normal"
      ~ "Nonmalignant bladder tissue",
      TRUE
      ~ NA_character_
    ),
    histologia_subtipo = "Papillary urothelial carcinoma"
  )

message("Standardized diagnostic distribution:")
print(table(vejiga_pheno_tbl_curado$grupo_diagnostico,
            useNA = "ifany"))

message("Tumor grade distribution (within tumor group):")
print(table(vejiga_pheno_tbl_curado$tumor_grade,
            useNA = "ifany"))


# 06. CONTROL DE CALIDAD DE LAS SONDAS _________________________________________

# Dejamos anotado que los valores son log-ratios de Agilent (microarray de dos
# colores, ya normalizado) y comprobamos algunas cosas antes de anotar.

message("Expression value range:")
message("  min = ",
        round(min(vejiga_exprs_matriz_inicial, na.rm = TRUE), 2))
message("  max = ",
        round(max(vejiga_exprs_matriz_inicial, na.rm = TRUE), 2))
message("  median = ",
        round(median(vejiga_exprs_matriz_inicial, na.rm = TRUE), 2))

# Los log-ratios de Agilent estan centrados cerca de cero, al reves que las
# intensidades log2 (que solo son positivas).
vejiga_valores_log_ratio_lgl <-
  min(vejiga_exprs_matriz_inicial, na.rm = TRUE) < 0 &&
  abs(median(vejiga_exprs_matriz_inicial, na.rm = TRUE)) < 1

message("Values consistent with two-color log-ratio scale: ",
        vejiga_valores_log_ratio_lgl)

# Comprobamos que el receptor de androgenos (AR) esta en el array. En una
# plataforma Agilent de genoma completo deberia estar; que la vejiga sea el
# comparador no lo excluye.

vejiga_ar_presente_lgl <-
  any(vejiga_feature_tbl_inicial$GENE_SYMBOL == "AR", na.rm = TRUE)

message("Androgen Receptor (AR) represented on array: ",
        vejiga_ar_presente_lgl)


# 07. ANOTACION SONDA A SIMBOLO ________________________________________________

# Pasamos las sondas Agilent a simbolo HGNC con el campo GENE_SYMBOL. Las sondas
# sin simbolo se descartan porque no se pueden unir con los otros datasets del
# metaanalisis a nivel de gen.

vejiga_feature_tbl_curado <-
  vejiga_feature_tbl_inicial %>%
  select(probe_id,
         gene_symbol = GENE_SYMBOL,
         ensembl_id  = ENSEMBL_ID) %>%
  mutate(
    gene_symbol = na_if(str_trim(gene_symbol), ""),
    ensembl_id  = na_if(str_trim(ensembl_id),  "")
  ) %>%
  filter(!is.na(gene_symbol))

message("Probes with a unique gene symbol: ",
        nrow(vejiga_feature_tbl_curado), " / ",
        nrow(vejiga_feature_tbl_inicial))


# 08. COLAPSAR SONDAS A SIMBOLOS DE GEN ________________________________________

# Cuando varias sondas miden el mismo gen, nos quedamos con la de mayor
# expresion media. Es la misma regla de colapso que en los otros microarrays
# (GSE3218, GSE8607).

vejiga_exprs_tbl_inicial <-
  vejiga_exprs_matriz_inicial %>%
  as.data.frame() %>%
  rownames_to_column("probe_id") %>%
  as_tibble() %>%
  inner_join(vejiga_feature_tbl_curado, by = "probe_id")

vejiga_exprs_tbl_curado <-
  vejiga_exprs_tbl_inicial %>%
  mutate(
    expresion_media = rowMeans(
      select(., -probe_id, -gene_symbol, -ensembl_id),
      na.rm = TRUE
    )
  ) %>%
  arrange(gene_symbol, desc(expresion_media)) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  select(-expresion_media, -probe_id)

message("Genes retained after probe collapse: ",
        nrow(vejiga_exprs_tbl_curado))


# 09. ENSAMBLAR LA MATRIZ ESTANDARIZADA ________________________________________

vejiga_genes_anotacion_tbl_final <-
  tibble(
    gene_id_native = vejiga_exprs_tbl_curado$gene_symbol,
    gene_symbol    = vejiga_exprs_tbl_curado$gene_symbol,
    ensembl_id     = vejiga_exprs_tbl_curado$ensembl_id
  )

vejiga_expresion_matriz_final <-
  vejiga_exprs_tbl_curado %>%
  select(-gene_symbol, -ensembl_id) %>%
  as.matrix()

rownames(vejiga_expresion_matriz_final) <-
  vejiga_genes_anotacion_tbl_final$gene_symbol

vejiga_metadatos_tbl_final <-
  vejiga_pheno_tbl_curado %>%
  select(sample_id = gsm_accession,
         gsm_accession,
         grupo_diagnostico,
         tejido_control_subtipo,
         tumor_grade,
         histologia_subtipo) %>%
  slice(match(colnames(vejiga_expresion_matriz_final), sample_id))

message("Standardized matrix assembled: ",
        nrow(vejiga_expresion_matriz_final), " genes x ",
        ncol(vejiga_expresion_matriz_final), " samples")


# 10. DESCRIPTOR DEL DATASET ___________________________________________________

vejiga_descriptor_lista <-
  list(
    gse_id                  = gse_objetivo_id_char,
    organo                  = organo_objetivo_char,
    plataforma_tipo         = plataforma_tipo_char,
    gene_id_convention      = "HGNC gene symbol (with Ensembl ID preserved)",
    matrix_value_nature     = "Two-color microarray log-ratio (Agilent, pre-normalized)",
    dge_branch_recommended  = "limma on log-ratio values",
    n_total                 = ncol(vejiga_expresion_matriz_final),
    n_tumor                 = sum(vejiga_metadatos_tbl_final$grupo_diagnostico
                                  == "Tumor", na.rm = TRUE),
    n_normal                = sum(vejiga_metadatos_tbl_final$grupo_diagnostico
                                  == "Normal", na.rm = TRUE),
    n_sin_etiqueta          = sum(is.na(vejiga_metadatos_tbl_final$grupo_diagnostico)),
    n_high_grade            = sum(vejiga_metadatos_tbl_final$tumor_grade
                                  == "High-grade", na.rm = TRUE),
    n_low_grade             = sum(vejiga_metadatos_tbl_final$tumor_grade
                                  == "Low-grade", na.rm = TRUE),
    histology               = "Papillary urothelial carcinoma",
    comparator_subtype      = "Nonmalignant bladder tissue (intra-dataset)",
    paired_design           = FALSE,
    role_in_meta_analysis   = "Low-androgen-dependency comparator",
    platform_companion      = "GPL8227 (miRNA) excluded as out-of-scope",
    feature_filtering       = "Probes with valid HGNC symbol only",
    duplicate_collapse      = "Maximum mean expression per gene symbol",
    fecha_estandarizacion   = format(Sys.Date(), "%Y-%m-%d")
  )


# 11. EXPORTAR EL RDS ESTANDARIZADO ____________________________________________

vejiga_dataset_estandarizado_lista <-
  list(
    conteos_matriz         = vejiga_expresion_matriz_final,
    genes_anotacion_tbl    = vejiga_genes_anotacion_tbl_final,
    muestras_metadatos_tbl = vejiga_metadatos_tbl_final,
    dataset_descriptor     = vejiga_descriptor_lista
  )

vejiga_rds_destino_path <-
  file.path(
    datasets_estandarizados_directorio_path,
    paste0(gse_objetivo_id_char, "_estandarizado.rds")
  )

saveRDS(vejiga_dataset_estandarizado_lista, vejiga_rds_destino_path)

message("Standardized object exported to: ",
        basename(vejiga_rds_destino_path))


# 12. RESUMEN __________________________________________________________________

message("")
message("Pipeline completed for ", gse_objetivo_id_char)
message("Organ: ",
        vejiga_descriptor_lista$organo)
message("Role in meta-analysis: ",
        vejiga_descriptor_lista$role_in_meta_analysis)
message("Platform: ",
        vejiga_descriptor_lista$plataforma_tipo)
message("Matrix nature: ",
        vejiga_descriptor_lista$matrix_value_nature)
message("DGE branch recommended: ",
        vejiga_descriptor_lista$dge_branch_recommended)
message("Histology: ",
        vejiga_descriptor_lista$histology)
message("Comparator subtype: ",
        vejiga_descriptor_lista$comparator_subtype)
message("Features retained: ",
        nrow(vejiga_dataset_estandarizado_lista$conteos_matriz),
        " unique gene symbols")
message("Samples retained: ",
        vejiga_descriptor_lista$n_total,
        " (Tumor = ", vejiga_descriptor_lista$n_tumor,
        ", Normal = ", vejiga_descriptor_lista$n_normal,
        if (vejiga_descriptor_lista$n_sin_etiqueta > 0)
          paste0(", Unlabeled = ",
                 vejiga_descriptor_lista$n_sin_etiqueta)
        else "",
        ")")
message("Tumor grade breakdown: High-grade = ",
        vejiga_descriptor_lista$n_high_grade,
        ", Low-grade = ",
        vejiga_descriptor_lista$n_low_grade)
message("Companion platform excluded: ",
        vejiga_descriptor_lista$platform_companion)
message("Output: ", vejiga_rds_destino_path)
