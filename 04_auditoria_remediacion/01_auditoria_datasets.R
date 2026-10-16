# ______________________________________________________________________________
# 04_01  Auditoria de datasets estandarizados y de resultados DGE
#
# Que hace: revisa los .rds estandarizados y los de expresion diferencial,
#           resume su estado (escala, anotacion, alineamiento, tamanos) y genera
#           alertas de compatibilidad para revision manual.
# Entra:    Datasets_Estandarizados/*.rds, DGE_Resultados/*.rds
# Sale:     Meta_Analisis_Resultados/auditoria_datasets.csv, auditoria_dge.csv,
#           alertas_compatibilidad.csv
# ______________________________________________________________________________

# LIBRERIAS

library(tidyverse)
library(org.Hs.eg.db)

# _______01. Configuracion revision _______________________________


# Establecemos los parametros que aplicaremos en los datasets

parametros_analisis_lista <- list(
  umbral_presencia_proporcion = 2 / 3,
  umbral_presencia_minimo_abs = 2,
  umbral_fdr                  = 0.05,
  semilla_aleatoria           = 1234
)

set.seed(1234)

# Dejamos la ruta con los archivos para que recupere los datasets estandarizados y aquellos
# de la matriz de expresion, y establecemos la carpeta de resultados

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_datasets_chr <- file.path(ruta_base_chr, "Datasets_Estandarizados")
ruta_dge_chr      <- file.path(ruta_base_chr, "DGE_Resultados")
ruta_salida_chr   <- file.path(ruta_base_chr, "Meta_Analisis_Resultados")


listar_archivos_rds <- function(ruta_carpeta_chr) {

  if (!dir.exists(ruta_carpeta_chr)) {
    stop("La carpeta no existe: ", ruta_carpeta_chr)
  }

  archivos_chr <- list.files(ruta_carpeta_chr,
                             pattern     = "\\.rds$",
                             full.names  = TRUE,
                             recursive   = TRUE,
                             ignore.case = TRUE)

  archivos_chr <- archivos_chr[!file.info(archivos_chr)$isdir]

  if (length(archivos_chr) == 0) {
    stop("No se encontraron archivos .rds en: ", ruta_carpeta_chr)
  }

  archivos_chr
}

archivos_datasets_chr <- listar_archivos_rds(ruta_datasets_chr)
archivos_dge_chr      <- listar_archivos_rds(ruta_dge_chr)

message("Datasets estandarizados encontrados: ", length(archivos_datasets_chr))
message("Archivos de DGE encontrados: ",         length(archivos_dge_chr))


# Verificamos que los campos Nulos esten bien codificados y codificamos cual es

extraer_campo <- function(lista_origen, campo_chr) {

  valor <- lista_origen[[campo_chr]]

  if (is.null(valor) || length(valor) == 0) {
    return(NA_character_)
  }

  as.character(valor[1])
}

contar_grupo <- function(metadatos_tbl, etiqueta_chr) {

  if (!"grupo_diagnostico" %in% names(metadatos_tbl)) {
    return(NA_integer_)
  }

  sum(metadatos_tbl$grupo_diagnostico == etiqueta_chr, na.rm = TRUE)
}


inferir_escala_matriz <- function(expresion_matriz) {

  n_columnas_int <- min(5, ncol(expresion_matriz))

  valores_num <- as.numeric(expresion_matriz[, seq_len(n_columnas_int)])
  valores_num <- valores_num[is.finite(valores_num)]

  if (length(valores_num) == 0) {
    return("Sin valores finitos")
  }

  todos_enteros_lgl <- all(abs(valores_num - round(valores_num)) < 1e-8)
  maximo_num        <- max(valores_num)

  if (todos_enteros_lgl && maximo_num > 100) {
    return("Probable conteos crudos")
  }

  if (maximo_num < 30) {
    return("Probable escala log2")
  }

  "Indeterminada"
}

# Ya con los cirterios determinados vamos a comprobar nuestras matrices y vamos
# a ver nuestro conteo final

auditar_dataset <- function(ruta_archivo_chr) {

  dataset_lista <- readRDS(ruta_archivo_chr)

  descriptor_lista <- dataset_lista$dataset_descriptor
  expresion_matriz <- dataset_lista$conteos_matriz
  metadatos_tbl    <- dataset_lista$muestras_metadatos_tbl
  anotacion_tbl    <- dataset_lista$genes_anotacion_tbl

  n_tumor_int  <- contar_grupo(metadatos_tbl, "Tumor")
  n_normal_int <- contar_grupo(metadatos_tbl, "Normal")

# Vamos a comprobar nuestras proporciones para determinar las varianzas


  n_efectivo_num <- (n_tumor_int * n_normal_int) / (n_tumor_int + n_normal_int)

  tibble(
    archivo             = basename(ruta_archivo_chr),
    gse_id              = extraer_campo(descriptor_lista, "gse_id"),
    organo              = extraer_campo(descriptor_lista, "organo"),
    plataforma          = extraer_campo(descriptor_lista, "plataforma_tipo"),
    convencion_id       = extraer_campo(descriptor_lista, "gene_id_convention"),
    naturaleza_valores  = extraer_campo(descriptor_lista, "matrix_value_nature"),
    rama_dge_sugerida   = extraer_campo(descriptor_lista, "dge_branch_recommended"),
    colapso_duplicados  = extraer_campo(descriptor_lista, "duplicate_collapse"),
    n_genes             = nrow(expresion_matriz),
    n_muestras          = ncol(expresion_matriz),
    n_tumor             = n_tumor_int,
    n_normal            = n_normal_int,
    n_efectivo          = round(n_efectivo_num, 2),
    escala_inferida     = inferir_escala_matriz(expresion_matriz),
    valor_minimo        = round(min(expresion_matriz, na.rm = TRUE), 2),
    valor_maximo        = round(max(expresion_matriz, na.rm = TRUE), 2),
    n_na_matriz         = sum(is.na(expresion_matriz)),
    simbolos_duplicados = sum(duplicated(rownames(expresion_matriz))),
    muestras_alineadas  = identical(colnames(expresion_matriz),
                                    metadatos_tbl$sample_id),
    genes_alineados     = identical(rownames(expresion_matriz),
                                    anotacion_tbl$gene_id_native)
  )
}

auditoria_datasets_tbl <- archivos_datasets_chr %>%
  map_dfr(auditar_dataset) %>%
  arrange(organo, gse_id)

# Repetimos el mismo proceso para los datos de expresion

auditar_dge <- function(ruta_archivo_chr) {

  dge_lista <- readRDS(ruta_archivo_chr)

  resultados_tbl <- dge_lista$dge_results_tbl
  metodo_lista   <- dge_lista$dge_method_descriptor
  padre_lista    <- dge_lista$parent_dataset_descriptor

  tibble(
    archivo_dge        = basename(ruta_archivo_chr),
    gse_id             = extraer_campo(padre_lista, "gse_id"),
    organo             = extraer_campo(padre_lista, "organo"),
    metodo_dge         = extraer_campo(metodo_lista, "method"),
    contraste          = extraer_campo(metodo_lista, "contrast"),
    nivel_referencia   = extraer_campo(metodo_lista, "reference_level"),
    filtro_expresion   = extraer_campo(metodo_lista, "low_expression_filter"),
    n_genes_pre_filtro = as.integer(extraer_campo(metodo_lista,
                                                  "n_genes_pre_filter")),
    n_genes_testados   = as.integer(extraer_campo(metodo_lista,
                                                  "n_genes_tested")),
    n_significativos   = as.integer(extraer_campo(metodo_lista,
                                                  "n_significant_padj_005")),
    columnas_resultado = paste(names(resultados_tbl), collapse = " | ")
  )
}

message("Auditando resultados de expresion diferencial...")

auditoria_dge_tbl <- archivos_dge_chr %>%
  map_dfr(auditar_dge) %>%
  arrange(organo, gse_id)


# Verificamos la integridad de los datos con los genes perdidos por los filtros
# de baja expresion


auditoria_integrada_tbl <- auditoria_datasets_tbl %>%
  left_join(
    auditoria_dge_tbl %>%
      select(gse_id, metodo_dge, n_genes_testados, n_significativos),
    by = "gse_id"
  ) %>%
  mutate(
    n_genes_perdidos_filtro = n_genes - n_genes_testados,
    pct_genes_perdidos      = round(100 * n_genes_perdidos_filtro / n_genes, 1)
  )


# Generamos alertas de verificacion manual

alertas_compatibilidad_tbl <- auditoria_integrada_tbl %>%
  transmute(
    gse_id,
    organo,
    alerta_escala        = escala_inferida == "Probable conteos crudos",
    alerta_anotacion     = convencion_id != "HGNC gene symbol",
    alerta_desalineado   = !muestras_alineadas | !genes_alineados,
    alerta_duplicados    = simbolos_duplicados > 0,
    alerta_valores_na    = n_na_matriz > 0,
    alerta_grupo_vacio   = n_tumor == 0 | n_normal == 0,
    alerta_n_efectivo    = n_efectivo < 5
  ) %>%
  mutate(n_alertas = rowSums(across(starts_with("alerta_")), na.rm = TRUE))


# 07. Resumen final

write_csv(auditoria_integrada_tbl,
          file.path(ruta_salida_chr, "auditoria_datasets.csv"))

write_csv(auditoria_dge_tbl,
          file.path(ruta_salida_chr, "auditoria_dge.csv"))

write_csv(alertas_compatibilidad_tbl,
          file.path(ruta_salida_chr, "alertas_compatibilidad.csv"))


message("\n________________ RESUMEN DE LA AUDITORIA ________________\n")

message("Datasets auditados: ", nrow(auditoria_integrada_tbl))

auditoria_integrada_tbl %>%
  count(organo, name = "n_datasets") %>%
  pmap(function(organo, n_datasets) {
    message("  ", organo, ": ", n_datasets, " datasets")
  }) %>%
  invisible()

message("\nEscalas de valores detectadas:")

auditoria_integrada_tbl %>%
  count(escala_inferida, name = "n") %>%
  pmap(function(escala_inferida, n) {
    message("  ", escala_inferida, ": ", n)
  }) %>%
  invisible()

message("\nConvenciones de identificador genico:")

auditoria_integrada_tbl %>%
  count(convencion_id, name = "n") %>%
  pmap(function(convencion_id, n) {
    message("  ", convencion_id, ": ", n)
  }) %>%
  invisible()

message("\nDatasets con alertas de compatibilidad: ",
        sum(alertas_compatibilidad_tbl$n_alertas > 0))

if (any(alertas_compatibilidad_tbl$n_alertas > 0)) {
  alertas_compatibilidad_tbl %>%
    filter(n_alertas > 0) %>%
    pmap(function(gse_id, organo, ...) {
      message("  Revisar: ", gse_id, " (", organo, ")")
    }) %>%
    invisible()
}

message("\nTablas exportadas en: ", ruta_salida_chr)
message("_________________________________________________________\n")
