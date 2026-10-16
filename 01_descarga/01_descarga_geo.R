# ______________________________________________________________________________
# 01_01  Descarga de cohortes de GEO (prostata y testiculo)
#
# Que hace: descarga de GEO las matrices de expresion y los metadatos de las
#           series de prostata y testiculo, baja los conteos crudos suplementarios
#           y audita los archivos guardados.
# Entra:    accesiones GEO (getGEO / getGEOSuppFiles)
# Sale:     Prostata/ y Testiculo/ (ExpressionSets, metadatos.csv, conteos suplementarios)
# ______________________________________________________________________________

# __ INDEX ___
# 01. LIBRERIAS
# 02. IMPORTAR DATOS CLINICOS Y MOLECULARES
# 03. AUDITORIA DE INTEGRIDAD DEL REPOSITORIO
# 04. RESUMEN DE SESION
# ________________________________________________________________________


# 01. LIBRERIAS ________________________________________________________________

# Aseguramos que BiocManager este disponible para los paquetes de Bioconductor.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(GEOquery)
})

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)


# 02. IMPORTAR DATOS CLINICOS Y MOLECULARES ____________________________________

# Accesiones GEO validadas para cada organo.
testiculo_gse_inicial_vector <- c(
  "GSE8607",
  "GSE3218"
)

prostata_gse_inicial_vector <- c(
  "GSE214940", "GSE134073", "GSE133626", "GSE221107",
  "GSE217260", "GSE104131", "GSE114740", "GSE82071"
)

# Carpetas de destino dentro de la ruta de trabajo.
testiculo_directorio <- file.path(ruta, "Testiculo")
prostata_directorio <- file.path(ruta, "Prostata")
prostata_conteos_directorio <- file.path(ruta, "Prostata_Counts")

dir.create(testiculo_directorio, recursive = TRUE, showWarnings = FALSE)
dir.create(prostata_directorio, recursive = TRUE, showWarnings = FALSE)
dir.create(prostata_conteos_directorio, recursive = TRUE, showWarnings = FALSE)


# Descarga la matriz de expresion y el fenotipo de una serie de GEO.
descargar_y_procesar_geo_gse <- function(gse_id, carpeta_destino) {

  message("----------------------------------------------------------------")
  message("INFO: Initiating acquisition for dataset: ", gse_id)

  tryCatch({
    # Descargamos la matriz de la serie de GEO.
    geo_gse_list <- GEOquery::getGEO(
      GEO = gse_id,
      GSEMatrix = TRUE,
      AnnotGPL = TRUE,
      destdir = carpeta_destino
    )

    if (length(geo_gse_list) > 0) {
      geo_expressionset_obj <- geo_gse_list[[1]]

      # Extraemos la matriz de expresion y los metadatos del fenotipo.
      gse_expresion_matriz <- Biobase::exprs(geo_expressionset_obj)
      gse_clinica_df <- Biobase::pData(geo_expressionset_obj)

      # Rutas de los archivos de salida.
      gse_expresion_ruta <- file.path(carpeta_destino, paste0(gse_id, "_expresion.csv"))
      gse_clinica_ruta <- file.path(carpeta_destino, paste0(gse_id, "_metadatos.csv"))
      gse_native_rds_ruta <- file.path(carpeta_destino, paste0(gse_id, "_ExpressionSet.rds"))

      # Guardamos los objetos en disco.
      write.csv(gse_expresion_matriz, file = gse_expresion_ruta, row.names = TRUE)
      write.csv(gse_clinica_df, file = gse_clinica_ruta, row.names = TRUE)
      saveRDS(geo_expressionset_obj, file = gse_native_rds_ruta)

      message("SUCCESS: Target ", gse_id, " successfully serialized to disk.")
      message("   - Features retained: ", nrow(gse_expresion_matriz))
      message("   - Clinical cases: ", ncol(gse_expresion_matriz))

    } else {
      message("WARNING: No expression matrix elements discovered for ", gse_id)
    }

  }, error = function(e) {
    message("CRITICAL ERROR: Failed download or parsing step for ", gse_id, " - ", e$message)
  })
}


# Recorremos cada cohorte y descargamos sus series.
message("EXECUTION: Compiling testicular cancer expressionsets")
for (gse_id in testiculo_gse_inicial_vector) {
  descargar_y_procesar_geo_gse(gse_id, testiculo_directorio)
}

message("EXECUTION: Compiling prostate cancer expressionsets")
for (gse_id in prostata_gse_inicial_vector) {
  descargar_y_procesar_geo_gse(gse_id, prostata_directorio)
}


# Descargamos los archivos suplementarios con los conteos crudos de RNA-Seq.
message("EXECUTION: Extracting supplementary raw sequencing count files")
for (gse_id in prostata_gse_inicial_vector) {

  message("----------------------------------------------------------------")
  message("INFO: Querying supplementary file space for: ", gse_id)

  tryCatch({
    gse_suplementarios_df <- GEOquery::getGEOSuppFiles(
      GSE = gse_id,
      baseDir = prostata_conteos_directorio
    )

    if (is.null(gse_suplementarios_df)) {
      message("WARNING: Unresolved supplementary files for target: ", gse_id)
    } else {
      message("SUCCESS: Supplementary raw data successfully stored for: ", gse_id)
    }

  }, error = function(e) {
    message("CRITICAL ERROR: Supplementary download aborted for ", gse_id, " - ", e$message)
  })
}


# 03. AUDITORIA DE INTEGRIDAD DEL REPOSITORIO __________________________________

# Extrae las accesiones GSE unicas a partir de los nombres de archivo.
extraer_gse_vector <- function(carpeta_ruta) {
  if (!dir.exists(carpeta_ruta)) return(character(0))

  archivos_vector <- list.files(carpeta_ruta, recursive = TRUE, full.names = FALSE)
  gse_detectados_vector <- unique(regmatches(archivos_vector, regexpr("GSE\\d+", archivos_vector)))
  gse_detectados_vector[gse_detectados_vector != ""]
}

# Lista los archivos de cada GSE con su tamano en MB.
listar_archivos_gse_df <- function(carpeta_ruta, gse_id) {
  if (!dir.exists(carpeta_ruta)) return(data.frame())

  archivos_vector <- list.files(carpeta_ruta, recursive = TRUE, full.names = TRUE)
  gse_especifico_patron <- paste0(gse_id, "(?!\\d)")
  archivos_filtrados_vector <- archivos_vector[grepl(gse_especifico_patron, archivos_vector, perl = TRUE)]

  if (length(archivos_filtrados_vector) == 0) return(data.frame())

  data.frame(
    gse_id = gse_id,
    archivo_nombre = basename(archivos_filtrados_vector),
    tamano_mb = round(file.info(archivos_filtrados_vector)$size / 1e6, 2),
    stringsAsFactors = FALSE
  )
}

# Inventario de los archivos descargados por carpeta.
message("================================================================")
message("AUDIT: Executing multi-cohort physical repository inventory")

subcarpetas_evaluadas_vector <- c("Testiculo", "Prostata", "Prostata_Counts")
cohortes_inventario_list <- list()

for (sub_carpeta in subcarpetas_evaluadas_vector) {
  target_carpeta_ruta <- file.path(ruta, sub_carpeta)
  if (!dir.exists(target_carpeta_ruta)) next

  message("Processing repository partition: ", sub_carpeta)
  gse_encontrados_vector <- extraer_gse_vector(target_carpeta_ruta)

  if (length(gse_encontrados_vector) > 0) {
    auditoria_detallada_df <- do.call(
      rbind,
      lapply(gse_encontrados_vector, function(g) listar_archivos_gse_df(target_carpeta_ruta, g))
    )

    for (g in sort(gse_encontrados_vector)) {
      sub_auditoria_df <- auditoria_detallada_df[auditoria_detallada_df$gse_id == g, ]
      message("  -> Accession ", g, ": ", nrow(sub_auditoria_df), " files detected (", sum(sub_auditoria_df$tamano_mb), " MB)")

      for (i in seq_len(nrow(sub_auditoria_df))) {
        truncacion_marca_char <- if (sub_auditoria_df$tamano_mb[i] < 0.5) " [POTENTIAL TRUNCATION ALERT]" else ""
        message("     · ", sub_auditoria_df$archivo_nombre[i], " (", sub_auditoria_df$tamano_mb[i], " MB)", truncacion_marca_char)
      }
    }
    cohortes_inventario_list[[sub_carpeta]] <- list(gses = gse_encontrados_vector, detalle = auditoria_detallada_df)
  }
}

# Comparamos lo descargado con la lista de accesiones esperadas.
message("================================================================")
message("AUDIT: Reconciling storage against golden manifest reference lists")

gse_descargados_detectados_vector <- unique(unlist(lapply(cohortes_inventario_list, function(x) x$gses)))
gse_totales_esperados_vector <- c(testiculo_gse_inicial_vector, prostata_gse_inicial_vector)

verificar_integridad_cohortes_funcion <- function(esperados_vector, cohorte_etiqueta_char) {
  presentes_vector <- intersect(esperados_vector, gse_descargados_detectados_vector)
  ausentes_vector <- setdiff(esperados_vector, gse_descargados_detectados_vector)

  message("Cohort: ", cohorte_etiqueta_char, " (n = ", length(esperados_vector), " expected)")
  message("  - Synchronized: ", length(presentes_vector), "/", length(esperados_vector))
  if (length(ausentes_vector) > 0) {
    message("  - MISSING ACCESSIONS: ", paste(ausentes_vector, collapse = ", "))
  } else {
    message("  - Manifest validation status: COMPLETE")
  }
}

verificar_integridad_cohortes_funcion(testiculo_gse_inicial_vector, "TESTICULO")
verificar_integridad_cohortes_funcion(prostata_gse_inicial_vector, "PROSTATA")

gse_no_esperados_vector <- setdiff(gse_descargados_detectados_vector, gse_totales_esperados_vector)
if (length(gse_no_esperados_vector) > 0) {
  message("NOTE: Unlisted testing/residual GSE directories detected: ", paste(gse_no_esperados_vector, collapse = ", "))
}


# 04. RESUMEN DE SESION ________________________________________________________

message("================================================================")
message("Analysis completed successfully")
message("Final localized repository: Data_TFM/")
message("Total GSE accessions validated: ", length(gse_descargados_detectados_vector), "/", length(gse_totales_esperados_vector))
message("Biomedical ingestion pipeline successfully executed")
