# ______________________________________________________________________________
# 03_01  Expresion diferencial con DESeq2 (datasets de conteos crudos)
#
# Que hace: ajusta DESeq2 por dataset (contraste Tumor vs Normal), aplica el
#           filtro de baja expresion y guarda los resultados en el esquema comun
#           del meta-analisis.
# Entra:    Datasets_Estandarizados/<GSE>_estandarizado.rds (GSE221107, GSE134073, GSE217260)
# Sale:     DGE_Resultados/<GSE>_dge_DESeq2.rds
# ______________________________________________________________________________

# Cubre los datasets de conteos crudos ligeros que corren en DESeq2: GSE221107, GSE134073 y GSE217260.

# __INDEX_________________________________________________________________
# 01. LIBRERIAS
# 02. RUTAS REPRODUCIBLES
# 03. DEFINIR LOS TRES DATASETS DE DESEQ2
# 04. MOTOR DESEQ2 POR DATASET
# 05. ITERAR SOBRE DATASETS Y ENSAMBLAR RESULTADOS
# 06. EXPORTAR CONTENEDORES DE RESULTADOS DGE
# 07. RESUMEN DE SESION
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


# __DEFINIR LOS TRES DATASETS DE DESEQ2___________________________________

deseq2_datasets_objetivo_tbl <-
  tribble(
    ~gse_id,       ~formula_disenio,                          ~requiere_reversion_log2,
    "GSE221107",   "~ paciente_id + grupo_diagnostico",       FALSE,
    "GSE134073",   "~ grupo_diagnostico",                     FALSE,
    "GSE217260",   "~ grupo_diagnostico",                     FALSE
  )

message("DESeq2 datasets included in this tactical run:")
print(deseq2_datasets_objetivo_tbl)


# __MOTOR DESEQ2 POR DATASET______________________________________________

procesar_dataset_deseq2_func <-
  function(gse_id_char,
           formula_disenio_char,
           requiere_reversion_log2_lgl) {

    message("\n>>> Processing ", gse_id_char,
            " with formula: ", formula_disenio_char)

    rds_path <- file.path(
      datasets_estandarizados_directorio_path,
      paste0(gse_id_char, "_estandarizado.rds")
    )

    if (!file.exists(rds_path)) {
      message("  WARNING: standardized .rds not found for ",
              gse_id_char, "; skipping")
      return(NULL)
    }

    dataset_estandarizado_lista <- readRDS(rds_path)

    conteos_matriz <- dataset_estandarizado_lista$conteos_matriz
    metadatos_tbl  <- dataset_estandarizado_lista$muestras_metadatos_tbl
    descriptor     <- dataset_estandarizado_lista$dataset_descriptor

    # Revierte la transformacion log2(conteo + 1) si hace falta
    if (requiere_reversion_log2_lgl) {
      message("  Reversing log2(count + 1) transformation")
      conteos_matriz <- round(2 ^ conteos_matriz - 1)
      conteos_matriz[conteos_matriz < 0] <- 0
    }

    # Sustituye los NA residuales por 0 y convierte a entero
    conteos_matriz[is.na(conteos_matriz)] <- 0
    storage.mode(conteos_matriz) <- "integer"

    # Reordena los metadatos para que coincidan con las columnas de la matriz
    # (dplyr:: explicito para evitar el enmascaramiento por IRanges::slice)
    metadatos_tbl <-
      metadatos_tbl %>%
      dplyr::slice(match(colnames(conteos_matriz), sample_id))

    # Quita las muestras sin etiqueta de diagnostico
    samples_validas_lgl <- !is.na(metadatos_tbl$grupo_diagnostico)

    if (sum(!samples_validas_lgl) > 0) {
      message("  Removing ", sum(!samples_validas_lgl),
              " samples with missing diagnostic label")
      conteos_matriz <- conteos_matriz[, samples_validas_lgl,
                                       drop = FALSE]
      metadatos_tbl  <- metadatos_tbl[samples_validas_lgl, ]
    }

    # Filtro de baja expresion
    n_grupo_menor_int <-
      min(table(metadatos_tbl$grupo_diagnostico))

    genes_pasan_filtro_lgl <-
      rowSums(conteos_matriz >= 10, na.rm = TRUE) >= n_grupo_menor_int

    n_genes_pre_filtro_int  <- nrow(conteos_matriz)
    n_genes_post_filtro_int <- sum(genes_pasan_filtro_lgl)

    conteos_matriz <- conteos_matriz[genes_pasan_filtro_lgl, ]

    message("  Low-expression filter: ",
            n_genes_post_filtro_int, " / ",
            n_genes_pre_filtro_int, " genes retained")

    # Convierte los metadatos a data.frame base con rownames para DESeq2
    coldata_df <- as.data.frame(metadatos_tbl)
    rownames(coldata_df) <- coldata_df$sample_id

    coldata_df$grupo_diagnostico <-
      factor(coldata_df$grupo_diagnostico,
             levels = c("Normal", "Tumor"))

    formula_objeto <- as.formula(formula_disenio_char)
    terminos_formula_vector <- all.vars(formula_objeto)

    if ("paciente_id" %in% terminos_formula_vector) {
      if (!"paciente_id" %in% colnames(coldata_df)) {
        message("  ERROR: paciente_id required by formula but not in metadata; skipping")
        return(NULL)
      }
      coldata_df$paciente_id <- factor(coldata_df$paciente_id)
    }

    if (!identical(colnames(conteos_matriz), rownames(coldata_df))) {
      message("  ERROR: matrix column names do not match coldata rownames")
      return(NULL)
    }

    # Ajuste de DESeq2
    deseq_dataset <-
      DESeqDataSetFromMatrix(
        countData = conteos_matriz,
        colData   = coldata_df,
        design    = formula_objeto
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
        log2_reversal_applied    = requiere_reversion_log2_lgl,
        low_expression_filter    = paste0("rowSums(counts >= 10) >= ",
                                          n_grupo_menor_int),
        n_genes_pre_filter       = n_genes_pre_filtro_int,
        n_genes_tested           = nrow(dge_resultados_tbl),
        n_significant_padj_005   = n_significativos_int,
        n_up_in_tumor            = n_up_int,
        n_down_in_tumor          = n_down_int,
        n_samples_total          = ncol(conteos_matriz),
        n_tumor                  = sum(coldata_df$grupo_diagnostico
                                       == "Tumor"),
        n_normal                 = sum(coldata_df$grupo_diagnostico
                                       == "Normal"),
        fecha_dge                = format(Sys.Date(), "%Y-%m-%d")
      )

    message("  DGE completed: ",
            nrow(dge_resultados_tbl), " genes tested, ",
            n_significativos_int, " significant at padj<0.05 (",
            n_up_int, " up / ", n_down_int, " down in Tumor)")

    # EXPORTACION INMEDIATA: guarda el .rds dentro de la funcion, antes de
    # pasar al siguiente dataset. Asi los datasets ya terminados quedan en
    # disco aunque uno posterior haga fallar el pipeline.

    dge_container_lista <-
      list(
        dge_results_tbl           = dge_resultados_tbl,
        dge_method_descriptor     = dge_method_descriptor,
        parent_dataset_descriptor = descriptor
      )

    dge_destino_path <- file.path(
      dge_resultados_directorio_path,
      paste0(gse_id_char, "_dge_DESeq2.rds")
    )

    saveRDS(dge_container_lista, dge_destino_path)

    message("  Exported: ", basename(dge_destino_path))

    return(dge_container_lista)
  }


# __ITERAR SOBRE DATASETS Y ENSAMBLAR RESULTADOS__________________________

deseq2_resultados_lista <-
  deseq2_datasets_objetivo_tbl %>%
  pmap(function(gse_id, formula_disenio, requiere_reversion_log2) {
    procesar_dataset_deseq2_func(
      gse_id_char                  = gse_id,
      formula_disenio_char         = formula_disenio,
      requiere_reversion_log2_lgl  = requiere_reversion_log2
    )
  }) %>%
  set_names(deseq2_datasets_objetivo_tbl$gse_id)


# __RESUMEN DE SESION_____________________________________________________

message("")
message("Tactical DESeq2 run completed (3 lightweight datasets)")
message("Datasets attempted: ", nrow(deseq2_datasets_objetivo_tbl))
message("Datasets successfully processed: ",
        sum(!sapply(deseq2_resultados_lista, is.null)))

dge_resumen_global_tbl <-
  deseq2_resultados_lista %>%
  discard(is.null) %>%
  imap_dfr(function(resultado, gse_id) {
    descriptor <- resultado$dge_method_descriptor
    tibble(
      gse_id                  = gse_id,
      n_samples               = descriptor$n_samples_total,
      n_tumor                 = descriptor$n_tumor,
      n_normal                = descriptor$n_normal,
      n_genes_tested          = descriptor$n_genes_tested,
      n_significant_padj_005  = descriptor$n_significant_padj_005,
      n_up                    = descriptor$n_up_in_tumor,
      n_down                  = descriptor$n_down_in_tumor
    )
  })

message("")
message("DESeq2 tactical run summary:")
print(dge_resumen_global_tbl)

message("")
message("Files written to: ", dge_resultados_directorio_path)
