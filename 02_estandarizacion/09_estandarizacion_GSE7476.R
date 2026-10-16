# ______________________________________________________________________________
# 02_09  Estandarizacion de GSE7476 (vejiga)
#
# Que hace: pasa el ExpressionSet de GSE7476 (microarray Affymetrix GPL570) al
#           esquema estandar del proyecto (simbolo HGNC, log2) y lo guarda. Cada
#           muestra es un pool biologico de 4-5 donantes.
# Entra:    Vejiga/GSE7476*ExpressionSet*.rds
# Sale:     Datasets_Estandarizados/GSE7476_estandarizado.rds
# ______________________________________________________________________________

# __ INDEX ___
# 01. LIBRERIAS
# 02. RUTAS
# 03. IMPORTAR EL EXPRESSIONSET
# 04. EXTRAER EXPRESION, ANOTACION Y FENOTIPO
# 05. ASIGNAR GRUPO DIAGNOSTICO Y GRADO/ESTADIO TUMORAL
# 06. DOCUMENTAR LA ARQUITECTURA DE MUESTRAS AGRUPADAS
# 07. CONTROL DE CALIDAD DE LAS SONDAS
# 08. ANOTACION SONDA A SIMBOLO
# 09. COLAPSAR SONDAS A SIMBOLOS DE GEN
# 10. ENSAMBLAR LA MATRIZ ESTANDARIZADA
# 11. DESCRIPTOR DEL DATASET
# 12. EXPORTAR EL RDS ESTANDARIZADO
# 13. RESUMEN
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

gse_objetivo_id_char  <- "GSE7476"
organo_objetivo_char  <- "Bladder"
plataforma_tipo_char  <- "Microarray (Affymetrix HG-U133 Plus 2.0, GPL570)"


# 03. IMPORTAR EL EXPRESSIONSET ________________________________________________

message("Auditing dataset ", gse_objetivo_id_char)

vejiga_eset_archivo_path <-
  list.files(
    path       = vejiga_metadatos_directorio_path,
    pattern    = paste0(gse_objetivo_id_char,
                        ".*ExpressionSet.*\\.rds$"),
    recursive  = TRUE,
    full.names = TRUE
  )

if (length(vejiga_eset_archivo_path) == 0) {
  stop("Critical error: ExpressionSet not found for ",
       gse_objetivo_id_char)
}

vejiga_eset <- readRDS(vejiga_eset_archivo_path[1])

message("ExpressionSet imported: ",
        nrow(vejiga_eset), " probes x ",
        ncol(vejiga_eset), " pools")
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
        ncol(vejiga_exprs_matriz_inicial), " pools")


# 05. ASIGNAR GRUPO DIAGNOSTICO Y GRADO/ESTADIO TUMORAL ________________________

# La etiqueta diagnostica viene en 'source_name_ch1', que junta grado y estadio
# del tumor. Los tres estratos tumorales (Ta bajo grado, T1 alto grado, T2-T4
# alto grado) se juntan en un unico grupo Tumor para el contraste principal;
# guardamos el grado/estadio original por si hace falta un analisis por estratos.

vejiga_pheno_tbl_curado <-
  vejiga_pheno_tbl_inicial %>%
  mutate(
    grupo_diagnostico = case_when(
      str_detect(source_name_ch1,
                 regex("tumoral|tumor|carcinoma|cancer",
                       ignore_case = TRUE))
      ~ "Tumor",
      str_detect(source_name_ch1,
                 regex("normal|nonmalignant|control|benign",
                       ignore_case = TRUE))
      ~ "Normal",
      TRUE
      ~ NA_character_
    ),
    tumor_grade_stage = case_when(
      str_detect(source_name_ch1,
                 regex("T1 high grade", ignore_case = TRUE))
      ~ "T1 high-grade NMIBC",
      str_detect(source_name_ch1,
                 regex("T2.*T3.*T4|T2, T3 or T4", ignore_case = TRUE))
      ~ "T2-T4 high-grade MIBC",
      str_detect(source_name_ch1,
                 regex("Ta low grade", ignore_case = TRUE))
      ~ "Ta low-grade NMIBC",
      TRUE
      ~ NA_character_
    ),
    tumor_invasion_class = case_when(
      str_detect(tumor_grade_stage, "MIBC") ~ "Muscle-invasive (MIBC)",
      str_detect(tumor_grade_stage, "NMIBC") ~ "Non-muscle-invasive (NMIBC)",
      TRUE                                   ~ NA_character_
    ),
    tejido_control_subtipo = case_when(
      grupo_diagnostico == "Normal"
      ~ "Pooled normal bladder tissue",
      TRUE
      ~ NA_character_
    ),
    histologia_subtipo = "Urothelial carcinoma"
  )

message("Standardized diagnostic distribution:")
print(table(vejiga_pheno_tbl_curado$grupo_diagnostico,
            useNA = "ifany"))

message("Tumor grade / stage distribution:")
print(table(vejiga_pheno_tbl_curado$tumor_grade_stage,
            useNA = "ifany"))

message("Tumor invasion class distribution:")
print(table(vejiga_pheno_tbl_curado$tumor_invasion_class,
            useNA = "ifany"))


# 06. DOCUMENTAR LA ARQUITECTURA DE MUESTRAS AGRUPADAS _________________________

# Sacamos el tamano de cada pool (4 o 5 donantes por libreria) para dejar
# trazabilidad. Ese dato esta en 'characteristics_ch1.2'.

vejiga_pheno_tbl_curado <-
  vejiga_pheno_tbl_curado %>%
  mutate(
    pool_size_text = str_extract(
      characteristics_ch1.2,
      regex("Pool of \\d+ individual samples", ignore_case = TRUE)),
    pool_size_n = as.integer(str_extract(pool_size_text, "\\d+"))
  )

vejiga_pool_sizes_summary_tbl <-
  vejiga_pheno_tbl_curado %>%
  count(grupo_diagnostico, pool_size_n)

message("Pool size composition by diagnostic group:")
print(vejiga_pool_sizes_summary_tbl)

vejiga_n_donantes_estimado_int <-
  sum(vejiga_pheno_tbl_curado$pool_size_n, na.rm = TRUE)

message("Total individual donors estimated across all pools: ",
        vejiga_n_donantes_estimado_int)


# 07. CONTROL DE CALIDAD DE LAS SONDAS _________________________________________

# Confirmamos que los valores estan en escala log2 RMA y que el receptor de
# androgenos (AR) esta en el array.

message("Expression value range:")
message("  min = ",
        round(min(vejiga_exprs_matriz_inicial, na.rm = TRUE), 2))
message("  max = ",
        round(max(vejiga_exprs_matriz_inicial, na.rm = TRUE), 2))
message("  median = ",
        round(median(vejiga_exprs_matriz_inicial, na.rm = TRUE), 2))

vejiga_valores_log2_lgl <-
  min(vejiga_exprs_matriz_inicial, na.rm = TRUE) >= 0 &&
  max(vejiga_exprs_matriz_inicial, na.rm = TRUE) < 30

message("Values consistent with log2 RMA scale: ",
        vejiga_valores_log2_lgl)

vejiga_ar_presente_lgl <-
  any(vejiga_feature_tbl_inicial$`Gene Symbol` == "AR", na.rm = TRUE)

message("Androgen Receptor (AR) represented on array: ",
        vejiga_ar_presente_lgl)


# 08. ANOTACION SONDA A SIMBOLO ________________________________________________

# Pasamos las sondas Affymetrix a simbolo HGNC con el campo 'Gene Symbol'. Las
# sondas que apuntan a varios genes a la vez (GENE1///GENE2) se descartan para
# no atribuir la senal de forma ambigua; misma regla que en GSE3218 y GSE8607.

vejiga_feature_tbl_curado <-
  vejiga_feature_tbl_inicial %>%
  select(probe_id, gene_symbol = `Gene Symbol`) %>%
  mutate(gene_symbol = na_if(str_trim(gene_symbol), "")) %>%
  filter(!is.na(gene_symbol)) %>%
  filter(!str_detect(gene_symbol, fixed("///")))

vejiga_probes_multigen_int <-
  vejiga_feature_tbl_inicial %>%
  filter(str_detect(`Gene Symbol`, fixed("///"))) %>%
  nrow()

message("Probes with a unique gene symbol: ",
        nrow(vejiga_feature_tbl_curado), " / ",
        nrow(vejiga_feature_tbl_inicial))
message("Multi-gene probes discarded: ", vejiga_probes_multigen_int)


# 09. COLAPSAR SONDAS A SIMBOLOS DE GEN ________________________________________

# Cuando varias sondas miden el mismo gen, nos quedamos con la de mayor
# expresion media. Es la misma regla de colapso que en todos los microarrays
# del metaanalisis.

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
      select(., -probe_id, -gene_symbol),
      na.rm = TRUE
    )
  ) %>%
  arrange(gene_symbol, desc(expresion_media)) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  select(-expresion_media, -probe_id)

message("Genes retained after probe collapse: ",
        nrow(vejiga_exprs_tbl_curado))


# 10. ENSAMBLAR LA MATRIZ ESTANDARIZADA ________________________________________

vejiga_genes_anotacion_tbl_final <-
  tibble(
    gene_id_native = vejiga_exprs_tbl_curado$gene_symbol,
    gene_symbol    = vejiga_exprs_tbl_curado$gene_symbol,
    ensembl_id     = NA_character_
  )

vejiga_expresion_matriz_final <-
  vejiga_exprs_tbl_curado %>%
  select(-gene_symbol) %>%
  as.matrix()

rownames(vejiga_expresion_matriz_final) <-
  vejiga_genes_anotacion_tbl_final$gene_symbol

vejiga_metadatos_tbl_final <-
  vejiga_pheno_tbl_curado %>%
  select(sample_id = gsm_accession,
         gsm_accession,
         grupo_diagnostico,
         tejido_control_subtipo,
         tumor_grade_stage,
         tumor_invasion_class,
         histologia_subtipo,
         pool_size_n) %>%
  slice(match(colnames(vejiga_expresion_matriz_final), sample_id))

message("Standardized matrix assembled: ",
        nrow(vejiga_expresion_matriz_final), " genes x ",
        ncol(vejiga_expresion_matriz_final), " pools")


# 11. DESCRIPTOR DEL DATASET ___________________________________________________

vejiga_descriptor_lista <-
  list(
    gse_id                    = gse_objetivo_id_char,
    organo                    = organo_objetivo_char,
    plataforma_tipo           = plataforma_tipo_char,
    gene_id_convention        = "HGNC gene symbol",
    matrix_value_nature       = "log2 RMA-normalized microarray intensities",
    dge_branch_recommended    = "limma on log2 intensities",
    n_total                   = ncol(vejiga_expresion_matriz_final),
    n_tumor                   = sum(vejiga_metadatos_tbl_final$grupo_diagnostico
                                    == "Tumor", na.rm = TRUE),
    n_normal                  = sum(vejiga_metadatos_tbl_final$grupo_diagnostico
                                    == "Normal", na.rm = TRUE),
    n_sin_etiqueta            = sum(is.na(vejiga_metadatos_tbl_final$grupo_diagnostico)),
    sample_architecture       = "Biological pools (4-5 individual donors per pool)",
    n_individual_donors_est   = vejiga_n_donantes_estimado_int,
    n_ta_low_grade            = sum(vejiga_metadatos_tbl_final$tumor_grade_stage
                                    == "Ta low-grade NMIBC", na.rm = TRUE),
    n_t1_high_grade           = sum(vejiga_metadatos_tbl_final$tumor_grade_stage
                                    == "T1 high-grade NMIBC", na.rm = TRUE),
    n_t2_4_high_grade         = sum(vejiga_metadatos_tbl_final$tumor_grade_stage
                                    == "T2-T4 high-grade MIBC", na.rm = TRUE),
    histology                 = "Urothelial carcinoma (Ta NMIBC, T1 NMIBC, T2-4 MIBC)",
    comparator_subtype        = "Pooled normal bladder tissue (intra-dataset)",
    paired_design             = FALSE,
    role_in_meta_analysis     = "Low-androgen-dependency comparator",
    methodological_caveat     = "Each library is a biological pool, not an individual; statistical inferences reflect between-pool variability",
    multi_gene_probes_removed = vejiga_probes_multigen_int,
    feature_filtering         = "Unique-symbol probes only; multi-gene probes discarded",
    duplicate_collapse        = "Maximum mean expression per gene symbol",
    fecha_estandarizacion     = format(Sys.Date(), "%Y-%m-%d")
  )


# 12. EXPORTAR EL RDS ESTANDARIZADO ____________________________________________

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


# 13. RESUMEN __________________________________________________________________

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
message("Sample architecture: ",
        vejiga_descriptor_lista$sample_architecture)
message("Estimated individual donors across all pools: ",
        vejiga_descriptor_lista$n_individual_donors_est)
message("Histology: ",
        vejiga_descriptor_lista$histology)
message("Comparator subtype: ",
        vejiga_descriptor_lista$comparator_subtype)
message("Features retained: ",
        nrow(vejiga_dataset_estandarizado_lista$conteos_matriz),
        " unique gene symbols")
message("Pools retained: ",
        vejiga_descriptor_lista$n_total,
        " (Tumor = ", vejiga_descriptor_lista$n_tumor,
        ", Normal = ", vejiga_descriptor_lista$n_normal,
        if (vejiga_descriptor_lista$n_sin_etiqueta > 0)
          paste0(", Unlabeled = ",
                 vejiga_descriptor_lista$n_sin_etiqueta)
        else "",
        ")")
message("Tumor stratification: ",
        "Ta low-grade NMIBC = ", vejiga_descriptor_lista$n_ta_low_grade,
        ", T1 high-grade NMIBC = ", vejiga_descriptor_lista$n_t1_high_grade,
        ", T2-T4 high-grade MIBC = ", vejiga_descriptor_lista$n_t2_4_high_grade)
message("Multi-gene probes removed: ",
        vejiga_descriptor_lista$multi_gene_probes_removed)
message("Methodological caveat: ",
        vejiga_descriptor_lista$methodological_caveat)
message("Output: ", vejiga_rds_destino_path)
