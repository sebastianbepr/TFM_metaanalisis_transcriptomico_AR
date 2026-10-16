# ______________________________________________________________________________
# 03_03  Expresion diferencial con DESeq2 (cohorte TCGA-TGCT + GTEx)
#
# Que hace: corre DESeq2 sobre la cohorte de testiculo TCGA-TGCT + GTEx
#           (contraste Tumor vs Normal), reconstruyendo los conteos desde la
#           escala log2 antes del ajuste, y guarda los resultados en el esquema
#           comun del meta-analisis.
# Entra:    Datasets_Estandarizados/TCGA_GTEX_Testis_estandarizado.rds
# Sale:     DGE_Resultados/TCGA_GTEX_Testis_dge_DESeq2.rds
# ______________________________________________________________________________

# __INDEX_________________________________________________________________
# 01. LIBRERIAS
# 02. RUTAS REPRODUCIBLES
# 03. PROCESAR TCGA + GTEX TESTICULO CON DESEQ2
# 04. EXPORTAR RESULTADOS DGE
# 05. RESUMEN DE SESION
# ________________________________________________________________________


# __LIBRERIAS_____________________________________________________________

library(DESeq2)
library(tidyverse)


# __RUTAS REPRODUCIBLES___________________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
base_datos_directorio_path <- ruta   # el resto del script usa este nombre sin cambios

datasets_estandarizados_directorio_path <-
  file.path(base_datos_directorio_path, "Datasets_Estandarizados")

dge_resultados_directorio_path <-
  file.path(base_datos_directorio_path, "DGE_Resultados")

if (!dir.exists(dge_resultados_directorio_path)) {
  dir.create(dge_resultados_directorio_path, recursive = TRUE)
}

dataset_id_char       <- "TCGA_GTEX_Testis"
formula_disenio_char  <- "~ grupo_diagnostico"


# __PROCESAR TCGA + GTEX TESTICULO CON DESEQ2_____________________________

message(">>> Processing ", dataset_id_char,
        " with formula: ", formula_disenio_char)

# Carga el contenedor estandarizado .rds
rds_path <- file.path(
  datasets_estandarizados_directorio_path,
  paste0(dataset_id_char, "_estandarizado.rds")
)

if (!file.exists(rds_path)) {
  stop("Standardized .rds not found: ", rds_path)
}

testiculo_dataset_estandarizado_lista <- readRDS(rds_path)

testiculo_conteos_matriz <-
  testiculo_dataset_estandarizado_lista$conteos_matriz

testiculo_metadatos_tbl <-
  testiculo_dataset_estandarizado_lista$muestras_metadatos_tbl

testiculo_descriptor <-
  testiculo_dataset_estandarizado_lista$dataset_descriptor

message("Input matrix: ",
        nrow(testiculo_conteos_matriz), " genes x ",
        ncol(testiculo_conteos_matriz), " samples")

# Revierte la transformacion log2(conteo + 1)
message("Reversing log2(count + 1) transformation")
testiculo_conteos_matriz <- round(2 ^ testiculo_conteos_matriz - 1)
testiculo_conteos_matriz[testiculo_conteos_matriz < 0] <- 0
testiculo_conteos_matriz[is.na(testiculo_conteos_matriz)] <- 0
storage.mode(testiculo_conteos_matriz) <- "integer"

# Reordena los metadatos para que coincidan con las columnas de la matriz
testiculo_metadatos_tbl <-
  testiculo_metadatos_tbl %>%
  dplyr::slice(match(colnames(testiculo_conteos_matriz), sample_id))

# Filtro de baja expresion
n_grupo_menor_int <-
  min(table(testiculo_metadatos_tbl$grupo_diagnostico))

testiculo_genes_pasan_filtro_lgl <-
  rowSums(testiculo_conteos_matriz >= 10, na.rm = TRUE) >= n_grupo_menor_int

n_genes_pre_filtro_int  <- nrow(testiculo_conteos_matriz)
n_genes_post_filtro_int <- sum(testiculo_genes_pasan_filtro_lgl)

testiculo_conteos_matriz <-
  testiculo_conteos_matriz[testiculo_genes_pasan_filtro_lgl, ]

message("Low-expression filter: ",
        n_genes_post_filtro_int, " / ",
        n_genes_pre_filtro_int, " genes retained")

# Convierte los metadatos a data.frame base con rownames para DESeq2
testiculo_coldata_df <- as.data.frame(testiculo_metadatos_tbl)
rownames(testiculo_coldata_df) <- testiculo_coldata_df$sample_id

testiculo_coldata_df$grupo_diagnostico <-
  factor(testiculo_coldata_df$grupo_diagnostico,
         levels = c("Normal", "Tumor"))

if (!identical(colnames(testiculo_conteos_matriz),
               rownames(testiculo_coldata_df))) {
  stop("Matrix column names do not match coldata rownames")
}

# Comprueba que AR sigue en la matriz filtrada
ar_present_lgl <- "ENSG00000169083" %in% rownames(testiculo_conteos_matriz)
message("AR (ENSG00000169083) present after low-expression filter: ",
        ar_present_lgl)

# Ajuste del modelo DESeq2
message("Fitting DESeq2 model (this can take 10-15 minutes for 411 samples)...")

deseq_dataset <-
  DESeqDataSetFromMatrix(
    countData = testiculo_conteos_matriz,
    colData   = testiculo_coldata_df,
    design    = as.formula(formula_disenio_char)
  )

deseq_fit <- DESeq(deseq_dataset, quiet = TRUE)

deseq_resultados <-
  results(deseq_fit,
          contrast = c("grupo_diagnostico", "Tumor", "Normal"),
          alpha    = 0.05)

# Esquema comun del meta-analisis
dge_resultados_tbl <-
  deseq_resultados %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  transmute(
    gene_id_native       = gene_id,
    log_fold_change      = log2FoldChange,
    average_expression   = log2(baseMean + 1),
    statistic            = stat,
    p_value              = pvalue,
    p_value_adjusted     = padj,
    direction_in_tumor   = case_when(
      log2FoldChange > 0 ~ "Up",
      log2FoldChange < 0 ~ "Down",
      TRUE               ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(p_value))

# Descriptor del metodo
n_significativos_int <-
  sum(dge_resultados_tbl$p_value_adjusted < 0.05, na.rm = TRUE)
n_up_int <-
  sum(dge_resultados_tbl$p_value_adjusted < 0.05 &
        dge_resultados_tbl$direction_in_tumor == "Up", na.rm = TRUE)
n_down_int <-
  sum(dge_resultados_tbl$p_value_adjusted < 0.05 &
        dge_resultados_tbl$direction_in_tumor == "Down", na.rm = TRUE)

dge_method_descriptor <-
  list(
    method                   = "DESeq2",
    formula                  = formula_disenio_char,
    contrast                 = "Tumor vs Normal",
    reference_level          = "Normal",
    log2_reversal_applied    = TRUE,
    low_expression_filter    = paste0("rowSums(counts >= 10) >= ",
                                      n_grupo_menor_int),
    n_genes_pre_filter       = n_genes_pre_filtro_int,
    n_genes_tested           = nrow(dge_resultados_tbl),
    n_significant_padj_005   = n_significativos_int,
    n_up_in_tumor            = n_up_int,
    n_down_in_tumor          = n_down_int,
    n_samples_total          = ncol(testiculo_conteos_matriz),
    n_tumor                  = sum(testiculo_coldata_df$grupo_diagnostico
                                   == "Tumor"),
    n_normal                 = sum(testiculo_coldata_df$grupo_diagnostico
                                   == "Normal"),
    methodological_caveat    = "origen_dataset cannot be modelled as covariate due to perfect confounding with grupo_diagnostico",
    fecha_dge                = format(Sys.Date(), "%Y-%m-%d")
  )

message("DGE completed: ",
        nrow(dge_resultados_tbl), " genes tested, ",
        n_significativos_int, " significant at padj<0.05 (",
        n_up_int, " up / ", n_down_int, " down in Tumor)")


# __EXPORTAR RESULTADOS DGE_______________________________________________

testiculo_dge_container_lista <-
  list(
    dge_results_tbl           = dge_resultados_tbl,
    dge_method_descriptor     = dge_method_descriptor,
    parent_dataset_descriptor = testiculo_descriptor
  )

dge_destino_path <- file.path(
  dge_resultados_directorio_path,
  paste0(dataset_id_char, "_dge_DESeq2.rds")
)

saveRDS(testiculo_dge_container_lista, dge_destino_path)

message("Exported: ", basename(dge_destino_path))


# __RESUMEN DE SESION_____________________________________________________

message("")
message("DESeq2 closing run completed for ", dataset_id_char)
message("Organ: ", testiculo_descriptor$organo)
message("Platform: ", testiculo_descriptor$plataforma_tipo)
message("Formula: ", dge_method_descriptor$formula)
message("Samples: ", dge_method_descriptor$n_samples_total,
        " (Tumor = ", dge_method_descriptor$n_tumor,
        ", Normal = ", dge_method_descriptor$n_normal, ")")
message("Genes pre-filter: ", dge_method_descriptor$n_genes_pre_filter)
message("Genes tested: ", dge_method_descriptor$n_genes_tested)
message("Significant (padj<0.05): ", dge_method_descriptor$n_significant_padj_005)
message("  Up in Tumor:   ", dge_method_descriptor$n_up_in_tumor)
message("  Down in Tumor: ", dge_method_descriptor$n_down_in_tumor)
message("AR detected in filtered matrix: ", ar_present_lgl)
message("Limitation: ", dge_method_descriptor$methodological_caveat)
message("Output: ", dge_destino_path)
