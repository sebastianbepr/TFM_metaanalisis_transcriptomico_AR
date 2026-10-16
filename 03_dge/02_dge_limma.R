# ______________________________________________________________________________
# 03_02  Expresion diferencial con limma (matrices prenormalizadas y microarrays)
#
# Que hace: ajusta limma por dataset (contraste Tumor vs Normal) despachando entre
#           lmFit directo sobre matrices prenormalizadas y la via voom para conteos
#           crudos, y guarda los resultados en el esquema comun del meta-analisis.
# Entra:    Datasets_Estandarizados/<GSE>_estandarizado.rds (7 datasets)
# Sale:     DGE_Resultados/<GSE>_dge_limma.rds
# ______________________________________________________________________________

# __INDEX_________________________________________________________________
# 01. LIBRERIAS
# 02. RUTAS REPRODUCIBLES
# 03. DATASETS Y FORMULAS DE LIMMA
# 04. MOTOR LIMMA POR DATASET
# 05. ITERAR SOBRE DATASETS CON EXPORTACION INMEDIATA
# 06. RESUMEN DE SESION
# ________________________________________________________________________


# __LIBRERIAS_____________________________________________________________

library(limma)
library(edgeR)
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


# __DATASETS Y FORMULAS DE LIMMA__________________________________________

# Cada fila declara el identificador del dataset, la formula del modelo lineal y
# si primero hay que aplicar voom (TRUE) o la matriz ya viene prenormalizada
# (FALSE). Los datasets que necesitan duplicateCorrelation lo indican con la
# bandera 'usa_duplicate_correlation'.

limma_datasets_objetivo_tbl <-
  tribble(
    ~gse_id,       ~formula_disenio,                          ~usa_voom, ~usa_duplicate_correlation,
    "GSE114740",   "~ paciente_id + grupo_diagnostico",       FALSE,     FALSE,
    "GSE104131",   "~ paciente_id + grupo_diagnostico",       FALSE,     FALSE,
    "GSE3218",     "~ grupo_diagnostico",                     FALSE,     FALSE,
    "GSE8607",     "~ grupo_diagnostico",                     FALSE,     FALSE,
    "GSE40355",    "~ grupo_diagnostico",                     FALSE,     FALSE,
    "GSE7476",     "~ grupo_diagnostico",                     FALSE,     FALSE,
    "GSE214940",   "~ grupo_diagnostico",                     TRUE,      TRUE
  )

message("limma branch datasets to process:")
print(limma_datasets_objetivo_tbl)


# __MOTOR LIMMA POR DATASET_______________________________________________

procesar_dataset_limma_func <-
  function(gse_id_char,
           formula_disenio_char,
           usa_voom_lgl,
           usa_duplicate_correlation_lgl) {

    message("\n>>> Processing ", gse_id_char,
            " with formula: ", formula_disenio_char,
            " | voom: ", usa_voom_lgl,
            " | dupcor: ", usa_duplicate_correlation_lgl)

    # Localiza y carga el contenedor estandarizado .rds
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

    expresion_matriz <- dataset_estandarizado_lista$conteos_matriz
    metadatos_tbl    <- dataset_estandarizado_lista$muestras_metadatos_tbl
    descriptor       <- dataset_estandarizado_lista$dataset_descriptor

    # Reordena los metadatos para que coincidan con las columnas de la matriz
    # (dplyr:: explicito)
    metadatos_tbl <-
      metadatos_tbl %>%
      dplyr::slice(match(colnames(expresion_matriz), sample_id))

    # Quita las muestras sin etiqueta de diagnostico
    samples_validas_lgl <- !is.na(metadatos_tbl$grupo_diagnostico)

    if (sum(!samples_validas_lgl) > 0) {
      message("  Removing ", sum(!samples_validas_lgl),
              " samples with missing diagnostic label")
      expresion_matriz <- expresion_matriz[, samples_validas_lgl,
                                           drop = FALSE]
      metadatos_tbl    <- metadatos_tbl[samples_validas_lgl, ]
    }

    # Sustituye los NA residuales de la matriz por cero (via voom) o por la media
    # de la fila (via prenormalizada); ambos son valores por defecto conservadores.
    n_nas_int <- sum(is.na(expresion_matriz))
    if (n_nas_int > 0) {
      message("  Imputing ", n_nas_int, " NA values")
      if (usa_voom_lgl) {
        expresion_matriz[is.na(expresion_matriz)] <- 0
      } else {
        row_means_vector <- rowMeans(expresion_matriz, na.rm = TRUE)
        na_indices <- which(is.na(expresion_matriz), arr.ind = TRUE)
        expresion_matriz[na_indices] <-
          row_means_vector[na_indices[, "row"]]
      }
    }

    # Define los factores con Normal como nivel de referencia del contraste
    metadatos_tbl <-
      metadatos_tbl %>%
      mutate(
        grupo_diagnostico = factor(grupo_diagnostico,
                                   levels = c("Normal", "Tumor"))
      )

    formula_objeto <- as.formula(formula_disenio_char)
    terminos_formula_vector <- all.vars(formula_objeto)

    if ("paciente_id" %in% terminos_formula_vector) {
      if (!"paciente_id" %in% colnames(metadatos_tbl)) {
        message("  ERROR: paciente_id required by formula but not in metadata; skipping")
        return(NULL)
      }
      metadatos_tbl$paciente_id <- factor(metadatos_tbl$paciente_id)
    }

    # Aplica un filtro de baja expresion / baja varianza segun la naturaleza del dato
    n_genes_pre_filtro_int <- nrow(expresion_matriz)

    if (usa_voom_lgl) {
      # Entrada de conteos enteros: aplica el filtro estilo DESeq2 para que los
      # pesos ajustados por voom se estimen solo sobre transcritos detectables.
      n_grupo_menor_int <-
        min(table(metadatos_tbl$grupo_diagnostico))
      genes_pasan_filtro_lgl <-
        rowSums(expresion_matriz >= 10, na.rm = TRUE) >= n_grupo_menor_int
      criterio_filtro_char <-
        paste0("rowSums(counts >= 10) >= ", n_grupo_menor_int)
    } else {
      # Entrada prenormalizada: filtra por la distribucion conjunta de expresion
      # y varianza; los genes por debajo de cualquiera de los dos umbrales del
      # percentil 10 se quitan por poco informativos.
      gen_expresion_media_vector <-
        rowMeans(expresion_matriz, na.rm = TRUE)
      gen_varianza_vector <-
        apply(expresion_matriz, 1, var, na.rm = TRUE)
      umbral_expresion_dbl <-
        quantile(gen_expresion_media_vector, 0.10, na.rm = TRUE)
      umbral_varianza_dbl <-
        quantile(gen_varianza_vector, 0.10, na.rm = TRUE)
      genes_pasan_filtro_lgl <-
        gen_expresion_media_vector >= umbral_expresion_dbl &
        gen_varianza_vector       >= umbral_varianza_dbl
      criterio_filtro_char <-
        paste0("mean expression >= P10 (",
               round(umbral_expresion_dbl, 3),
               ") AND variance >= P10 (",
               round(umbral_varianza_dbl, 4), ")")
    }

    expresion_matriz <- expresion_matriz[genes_pasan_filtro_lgl, ]

    n_genes_post_filtro_int <- nrow(expresion_matriz)

    message("  Filter: ", n_genes_post_filtro_int, " / ",
            n_genes_pre_filtro_int, " genes retained (",
            criterio_filtro_char, ")")

    # Construye la matriz de diseno del modelo
    coldata_df <- as.data.frame(metadatos_tbl)
    rownames(coldata_df) <- coldata_df$sample_id

    diseno_matriz <- model.matrix(formula_objeto, data = coldata_df)

    # Comprobacion: rango de la matriz de diseno
    rango_diseno_int <- qr(diseno_matriz)$rank

    if (rango_diseno_int < ncol(diseno_matriz)) {
      message("  WARNING: design matrix is rank-deficient (",
              rango_diseno_int, " / ", ncol(diseno_matriz),
              "); some terms may be dropped by lmFit")
    }

    # Identifica la columna de la matriz de diseno que lleva el contraste
    # Tumor - Normal
    coef_contraste_char <-
      grep("grupo_diagnostico", colnames(diseno_matriz),
           value = TRUE)[1]

    if (is.na(coef_contraste_char)) {
      message("  ERROR: cannot locate grupo_diagnostico coefficient in design matrix")
      return(NULL)
    }

    # Ajusta el modelo por la via de limma que corresponda
    if (usa_voom_lgl) {

      # Modo B: conteos crudos -> voom -> lmFit (con dupcor opcional)
      message("  Running voom (mean-variance modelling for counts)")

      # Asegura almacenamiento entero para edgeR / voom
      storage.mode(expresion_matriz) <- "integer"
      dge_list <- DGEList(counts = expresion_matriz)
      dge_list <- calcNormFactors(dge_list)

      voom_objeto <- voom(dge_list, design = diseno_matriz, plot = FALSE)

      if (usa_duplicate_correlation_lgl) {
        if (!"paciente_id" %in% colnames(coldata_df)) {
          message("  ERROR: paciente_id required for duplicateCorrelation but not in metadata; skipping")
          return(NULL)
        }
        message("  Estimating duplicateCorrelation by paciente_id")
        dupcor_resultado <-
          duplicateCorrelation(voom_objeto,
                               design = diseno_matriz,
                               block  = coldata_df$paciente_id)
        message("  Intra-patient correlation: ",
                round(dupcor_resultado$consensus.correlation, 4))
        # Vuelve a correr voom con la correlacion estimada
        voom_objeto <- voom(dge_list,
                            design       = diseno_matriz,
                            block        = coldata_df$paciente_id,
                            correlation  = dupcor_resultado$consensus.correlation,
                            plot         = FALSE)
        # Reestima duplicateCorrelation sobre el segundo voom (buena practica)
        dupcor_resultado <-
          duplicateCorrelation(voom_objeto,
                               design = diseno_matriz,
                               block  = coldata_df$paciente_id)
        message("  Re-estimated correlation: ",
                round(dupcor_resultado$consensus.correlation, 4))

        ajuste_fit <-
          lmFit(voom_objeto,
                design       = diseno_matriz,
                block        = coldata_df$paciente_id,
                correlation  = dupcor_resultado$consensus.correlation)
      } else {
        ajuste_fit <- lmFit(voom_objeto, design = diseno_matriz)
      }

      metodo_limma_char <- "limma-voom"

    } else {
      # Modo A: matriz prenormalizada -> lmFit directo
      ajuste_fit <- lmFit(expresion_matriz, design = diseno_matriz)
      metodo_limma_char <- "limma (pre-normalized)"
    }

    ajuste_eBayes <- eBayes(ajuste_fit)

    # Extrae el contraste Tumor - Normal como tabla de resultados comun
    top_table_df <-
      topTable(ajuste_eBayes,
               coef      = coef_contraste_char,
               number    = Inf,
               sort.by   = "none",
               adjust.method = "BH")

    dge_resultados_tbl <-
      top_table_df %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>%
      as_tibble() %>%
      transmute(
        gene_id_native       = gene_id,
        log_fold_change      = logFC,
        average_expression   = AveExpr,
        statistic            = t,
        p_value              = P.Value,
        p_value_adjusted     = adj.P.Val,
        direction_in_tumor   = case_when(
          logFC > 0 ~ "Up",
          logFC < 0 ~ "Down",
          TRUE      ~ NA_character_
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
        method                       = metodo_limma_char,
        formula                      = formula_disenio_char,
        contrast                     = "Tumor vs Normal",
        contrast_coefficient         = coef_contraste_char,
        reference_level              = "Normal",
        voom_applied                 = usa_voom_lgl,
        duplicate_correlation_used   = usa_duplicate_correlation_lgl,
        treat_used                   = FALSE,
        low_expression_filter        = criterio_filtro_char,
        n_genes_pre_filter           = n_genes_pre_filtro_int,
        n_genes_tested               = nrow(dge_resultados_tbl),
        n_significant_padj_005       = n_significativos_int,
        n_up_in_tumor                = n_up_int,
        n_down_in_tumor              = n_down_int,
        n_samples_total              = ncol(expresion_matriz),
        n_tumor                      = sum(metadatos_tbl$grupo_diagnostico
                                           == "Tumor", na.rm = TRUE),
        n_normal                     = sum(metadatos_tbl$grupo_diagnostico
                                           == "Normal", na.rm = TRUE),
        fecha_dge                    = format(Sys.Date(), "%Y-%m-%d")
      )

    message("  DGE completed: ",
            nrow(dge_resultados_tbl), " genes tested, ",
            n_significativos_int, " significant at padj<0.05 (",
            n_up_int, " up / ", n_down_int, " down in Tumor)")

    # Ensambla el contenedor comun
    dge_container_lista <-
      list(
        dge_results_tbl           = dge_resultados_tbl,
        dge_method_descriptor     = dge_method_descriptor,
        parent_dataset_descriptor = descriptor
      )

    # EXPORTACION INMEDIATA: guarda el .rds antes de pasar al siguiente dataset
    # para no perder el progreso parcial
    dge_destino_path <- file.path(
      dge_resultados_directorio_path,
      paste0(gse_id_char, "_dge_limma.rds")
    )

    saveRDS(dge_container_lista, dge_destino_path)

    message("  Exported: ", basename(dge_destino_path))

    return(dge_container_lista)
  }


# __ITERAR SOBRE DATASETS CON EXPORTACION INMEDIATA_______________________

limma_resultados_lista <-
  limma_datasets_objetivo_tbl %>%
  pmap(function(gse_id, formula_disenio,
                usa_voom, usa_duplicate_correlation) {
    procesar_dataset_limma_func(
      gse_id_char                       = gse_id,
      formula_disenio_char              = formula_disenio,
      usa_voom_lgl                      = usa_voom,
      usa_duplicate_correlation_lgl     = usa_duplicate_correlation
    )
  }) %>%
  set_names(limma_datasets_objetivo_tbl$gse_id)


# __RESUMEN DE SESION_____________________________________________________

message("")
message("limma differential expression branch completed")
message("Datasets attempted: ", nrow(limma_datasets_objetivo_tbl))
message("Datasets successfully processed: ",
        sum(!sapply(limma_resultados_lista, is.null)))

dge_resumen_global_tbl <-
  limma_resultados_lista %>%
  discard(is.null) %>%
  imap_dfr(function(resultado, gse_id) {
    descriptor <- resultado$dge_method_descriptor
    tibble(
      gse_id                  = gse_id,
      method                  = descriptor$method,
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
message("limma branch summary table:")
print(dge_resumen_global_tbl)

message("")
message("Output directory: ", dge_resultados_directorio_path)
message("Each dataset emits <GSE>_dge_limma.rds for downstream meta-analysis")
