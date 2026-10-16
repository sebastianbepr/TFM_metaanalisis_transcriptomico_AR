# ______________________________________________________________________________
# 02_07  Estandarizacion de GSE3218 (testiculo)
#
# Que hace: pasa el ExpressionSet de GSE3218 (microarray Affymetrix GPL96) al
#           esquema estandar del proyecto (simbolo HGNC, matriz log2) y lo guarda.
# Entra:    Testiculo/GSE3218*ExpressionSet*.rds
# Sale:     Datasets_Estandarizados/GSE3218_estandarizado.rds
# ______________________________________________________________________________

# __ INDEX ___
# 01. LIBRERIAS
# 02. RUTAS
# 03. IMPORTAR EL EXPRESSIONSET
# 04. EXTRAER EXPRESION, ANOTACION Y FENOTIPO
# 05. ASIGNAR GRUPO DIAGNOSTICO Y CONSERVAR HISTOLOGIA
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

testiculo_directorio_path <-
  file.path(base_datos_directorio_path, "Testiculo")

datasets_estandarizados_directorio_path <-
  file.path(base_datos_directorio_path, "Datasets_Estandarizados")

if (!dir.exists(datasets_estandarizados_directorio_path)) {
  dir.create(datasets_estandarizados_directorio_path, recursive = TRUE)
}

gse_objetivo_id_char  <- "GSE3218"
organo_objetivo_char  <- "Testis"
plataforma_tipo_char  <- "Microarray (Affymetrix HG-U133A, GPL96)"


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

# Separamos el contenedor de Bioconductor en sus tres tablas (expresion,
# anotacion y fenotipo) para seguir trabajando con tidyverse.

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
message("Feature annotation columns: ",
        ncol(testiculo_feature_tbl_inicial))
message("Phenotype annotation columns: ",
        ncol(testiculo_pheno_tbl_inicial))


# 05. ASIGNAR GRUPO DIAGNOSTICO Y CONSERVAR HISTOLOGIA _________________________

# Pasamos la etiqueta diagnostica al vocabulario comun Tumor / Normal. Los
# tumores de celulas germinales son el grupo Tumor y los testiculos normales
# (individuales y el agrupado) son el grupo Normal. Marcamos cual control es un
# pool en el campo de subtipo.

testiculo_pheno_tbl_curado <-
  testiculo_pheno_tbl_inicial %>%
  mutate(
    grupo_diagnostico = case_when(
      str_detect(source_name_ch1,
                 regex("germ cell tumor|tumor|carcinoma|gct",
                       ignore_case = TRUE))
      ~ "Tumor",
      str_detect(source_name_ch1,
                 regex("normal testis|normal|control",
                       ignore_case = TRUE))
      ~ "Normal",
      TRUE
      ~ NA_character_
    ),
    tejido_control_subtipo = case_when(
      str_detect(source_name_ch1,
                 regex("pooled", ignore_case = TRUE))
      ~ "Pooled normal testis",
      grupo_diagnostico == "Normal"
      ~ "Individual normal testis",
      TRUE
      ~ NA_character_
    ),
    histologia_subtipo = characteristics_ch1,
    histologia_clase = case_when(
      str_detect(characteristics_ch1,
                 regex("^seminoma|pure seminoma|classic seminoma",
                       ignore_case = TRUE)) &
        !str_detect(characteristics_ch1,
                    regex("embryonal|teratoma|yolk|chorio|mixed",
                          ignore_case = TRUE))
      ~ "Seminoma",
      grupo_diagnostico == "Tumor"
      ~ "Non-seminoma / mixed",
      TRUE
      ~ NA_character_
    )
  )

message("Standardized diagnostic distribution:")
print(table(testiculo_pheno_tbl_curado$grupo_diagnostico,
            useNA = "ifany"))

message("Control subtype distribution:")
print(table(testiculo_pheno_tbl_curado$tejido_control_subtipo,
            useNA = "ifany"))

message("Histological class distribution (tumors):")
print(table(testiculo_pheno_tbl_curado$histologia_clase,
            useNA = "ifany"))


# 06. CONTROL DE CALIDAD DE LAS SONDAS _________________________________________

# Comprobamos que las intensidades estan en escala log2 antes de anotar.

message("Expression value range:")
message("  min = ",
        round(min(testiculo_exprs_matriz_inicial, na.rm = TRUE), 2))
message("  max = ",
        round(max(testiculo_exprs_matriz_inicial, na.rm = TRUE), 2))
message("  median = ",
        round(median(testiculo_exprs_matriz_inicial, na.rm = TRUE), 2))

testiculo_valores_log2_lgl <-
  max(testiculo_exprs_matriz_inicial, na.rm = TRUE) < 30

message("Values consistent with log2 scale: ",
        testiculo_valores_log2_lgl)

# Comprobamos que el receptor de androgenos (AR) esta en el array, que es el
# eje central del metaanalisis.

testiculo_ar_presente_lgl <-
  any(str_detect(testiculo_feature_tbl_inicial$`Gene symbol`,
                 regex("\\bAR\\b", ignore_case = FALSE)),
      na.rm = TRUE)

message("Androgen Receptor (AR) represented on array: ",
        testiculo_ar_presente_lgl)


# 07. ANOTACION SONDA A SIMBOLO ________________________________________________

# Pasamos las sondas Affymetrix a simbolo HGNC con el campo 'Gene symbol'. Las
# sondas que apuntan a varios genes a la vez (GENE1///GENE2) son ambiguas y se
# descartan porque su senal no se puede atribuir a un solo gen.

testiculo_feature_tbl_curado <-
  testiculo_feature_tbl_inicial %>%
  select(probe_id, gene_symbol = `Gene symbol`) %>%
  mutate(
    gene_symbol = na_if(str_trim(gene_symbol), "")
  ) %>%
  filter(!is.na(gene_symbol)) %>%
  filter(!str_detect(gene_symbol, fixed("///")))

message("Probes with a unique gene symbol: ",
        nrow(testiculo_feature_tbl_curado), " / ",
        nrow(testiculo_feature_tbl_inicial))

testiculo_probes_multigen_int <-
  testiculo_feature_tbl_inicial %>%
  filter(str_detect(`Gene symbol`, fixed("///"))) %>%
  nrow()

message("Multi-gene probes discarded: ", testiculo_probes_multigen_int)


# 08. COLAPSAR SONDAS A SIMBOLOS DE GEN ________________________________________

# Cuando varias sondas miden el mismo gen, nos quedamos con la de mayor
# expresion media. Es la regla de colapso que usamos en todo el metaanalisis.

testiculo_exprs_tbl_inicial <-
  testiculo_exprs_matriz_inicial %>%
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

# Construimos la matriz numerica indexada por simbolo de gen, la tabla de
# anotacion y los metadatos en el mismo orden de columnas que la matriz.

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

# Metadatos de muestra en el esquema comun, en el mismo orden que las columnas
# de la matriz.

testiculo_metadatos_tbl_final <-
  testiculo_pheno_tbl_curado %>%
  select(sample_id = gsm_accession,
         gsm_accession,
         grupo_diagnostico,
         tejido_control_subtipo,
         histologia_subtipo,
         histologia_clase) %>%
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
    matrix_value_nature     = "log2 RMA-normalized microarray intensities",
    dge_branch_recommended  = "limma on log2 intensities",
    n_total                 = ncol(testiculo_expresion_matriz_final),
    n_tumor                 = sum(testiculo_metadatos_tbl_final$grupo_diagnostico
                                  == "Tumor", na.rm = TRUE),
    n_normal                = sum(testiculo_metadatos_tbl_final$grupo_diagnostico
                                  == "Normal", na.rm = TRUE),
    n_sin_etiqueta          = sum(is.na(testiculo_metadatos_tbl_final$grupo_diagnostico)),
    comparator_subtype      = "Normal testis (5 individual + 1 pooled)",
    n_seminoma              = sum(testiculo_metadatos_tbl_final$histologia_clase
                                  == "Seminoma", na.rm = TRUE),
    n_no_seminoma           = sum(testiculo_metadatos_tbl_final$histologia_clase
                                  == "Non-seminoma / mixed", na.rm = TRUE),
    paired_design           = FALSE,
    multi_gene_probes_removed = testiculo_probes_multigen_int,
    feature_filtering       = "Unique-symbol probes only; multi-gene probes discarded",
    duplicate_collapse      = "Maximum mean expression per gene symbol",
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
message("Comparator subtype: ",
        testiculo_descriptor_lista$comparator_subtype)
message("Multi-gene probes removed: ",
        testiculo_descriptor_lista$multi_gene_probes_removed)
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
message("Histological breakdown: Seminoma = ",
        testiculo_descriptor_lista$n_seminoma,
        ", Non-seminoma / mixed = ",
        testiculo_descriptor_lista$n_no_seminoma)
message("Gene identifier convention: ",
        testiculo_descriptor_lista$gene_id_convention)
message("Output: ", testiculo_rds_destino_path)
