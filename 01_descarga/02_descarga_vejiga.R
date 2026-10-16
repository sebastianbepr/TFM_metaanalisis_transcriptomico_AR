# ______________________________________________________________________________
# 01_02  Descarga de cohortes de vejiga (comparador)
#
# Que hace: descarga de GEO las cuatro series de cancer de vejiga que forman el
#           brazo comparador de baja dependencia androgenica, guarda los
#           ExpressionSets y metadatos, baja los suplementarios y audita el repositorio.
# Entra:    series GEO de vejiga (GSE40355, GSE7476, GSE23732, GSE56037) via getGEO
# Sale:     Vejiga/ (ExpressionSets, metadatos.csv) y Vejiga_Counts/
# ______________________________________________________________________________

# __ INDEX ___
# 01. LIBRERIAS
# 02. RUTAS REPRODUCIBLES
# 03. DATASETS DE VEJIGA OBJETIVO
# 04. DESCARGA DE METADATOS Y EXPRESSIONSETS
# 05. DESCARGA DE ARCHIVOS SUPLEMENTARIOS
# 06. AUDITORIA DEL REPOSITORIO
# 07. RESUMEN DE SESION
# ________________________________________________________________________


# 01. LIBRERIAS ________________________________________________________________

library(tidyverse)
library(GEOquery)
library(Biobase)


# 02. RUTAS REPRODUCIBLES ______________________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
base_datos_directorio_path <- ruta   # el resto del script usa este nombre sin cambios

vejiga_metadatos_directorio_path <-
  file.path(base_datos_directorio_path, "Vejiga")

vejiga_conteos_directorio_path <-
  file.path(base_datos_directorio_path, "Vejiga_Counts")

# Creamos las carpetas de destino si no existen.

if (!dir.exists(vejiga_metadatos_directorio_path)) {
  dir.create(vejiga_metadatos_directorio_path, recursive = TRUE)
}

if (!dir.exists(vejiga_conteos_directorio_path)) {
  dir.create(vejiga_conteos_directorio_path, recursive = TRUE)
}

message("Vejiga metadata directory: ", vejiga_metadatos_directorio_path)
message("Vejiga counts directory: ", vejiga_conteos_directorio_path)


# 03. DATASETS DE VEJIGA OBJETIVO ______________________________________________

# Cada fila documenta el motivo de inclusion del dataset.

vejiga_datasets_objetivo_tbl <-
  tribble(
    ~gse_id,      ~n_total, ~plataforma_esperada,                     ~rationale,
    "GSE40355",   48,       "Microarray (mRNA expression)",            "Most robust bladder candidate, full mRNA expression profiling",
    "GSE7476",    12,       "Microarray (expression)",                 "Bladder cancer pathological categories",
    "GSE23732",   8,        "Microarray (whole-genome expression)",    "Whole-genome expression of urothelial carcinomas (small but informative)",
    "GSE56037",   6,        "Microarray (expression)",                 "Urothelial cell carcinoma-focused profiling (smallest cohort retained)"
  )

message("Target bladder datasets:")
print(vejiga_datasets_objetivo_tbl)


# 04. DESCARGA DE METADATOS Y EXPRESSIONSETS ___________________________________

# Para cada dataset bajamos de GEO el ExpressionSet y guardamos tres archivos:
#   - <GSE>_ExpressionSet.rds   (el objeto de Bioconductor completo)
#   - <GSE>_metadatos.csv       (tabla de fenotipo, lista para fread)
#   - <GSE>_series_matrix.txt.gz (que GEOquery cachea sola)
# Algunas series traen varias plataformas (cada una es un ExpressionSet) y se
# guardan por separado.

descargar_serie_geo_func <- function(gse_id) {

  message("\n>>> Downloading ", gse_id, " from GEO")

  serie_resultado_lista <-
    tryCatch(
      getGEO(gse_id,
             destdir = vejiga_metadatos_directorio_path,
             GSEMatrix = TRUE,
             getGPL    = TRUE),
      error = function(e) {
        message("  ERROR downloading ", gse_id, ": ",
                conditionMessage(e))
        return(NULL)
      }
    )

  if (is.null(serie_resultado_lista)) return(invisible(NULL))

  # Una serie puede tener varios ExpressionSets, uno por plataforma. Recorremos
  # todos y guardamos cada uno, anadiendo el identificador de plataforma al
  # nombre cuando hay mas de una.

  n_plataformas_int <- length(serie_resultado_lista)

  for (i in seq_len(n_plataformas_int)) {

    eset_individual <- serie_resultado_lista[[i]]

    sufijo_plataforma_char <-
      if (n_plataformas_int > 1) {
        paste0("_", annotation(eset_individual))
      } else {
        ""
      }

    # Guardamos el ExpressionSet.
    rds_path <- file.path(
      vejiga_metadatos_directorio_path,
      paste0(gse_id, sufijo_plataforma_char, "_ExpressionSet.rds")
    )
    saveRDS(eset_individual, rds_path)

    # Guardamos la tabla de fenotipo.
    pheno_tbl <-
      pData(eset_individual) %>%
      rownames_to_column("geo_accession_row") %>%
      as_tibble()

    csv_path <- file.path(
      vejiga_metadatos_directorio_path,
      paste0(gse_id, sufijo_plataforma_char, "_metadatos.csv")
    )
    write_csv(pheno_tbl, csv_path)

    message("  Saved: ", basename(rds_path))
    message("         ", basename(csv_path),
            " (", nrow(pheno_tbl), " samples)")
  }

  return(invisible(serie_resultado_lista))
}

# Recorremos todos los datasets objetivo acumulando los resultados.

vejiga_descarga_resultados_lista <-
  vejiga_datasets_objetivo_tbl$gse_id %>%
  set_names() %>%
  map(descargar_serie_geo_func)


# 05. DESCARGA DE ARCHIVOS SUPLEMENTARIOS ______________________________________

# Los archivos suplementarios traen las matrices de expresion crudas o
# procesadas de los autores. Se bajan a una subcarpeta por dataset dentro de
# Vejiga_Counts, con la misma estructura que prostata y testiculo.

descargar_suplementarios_func <- function(gse_id) {

  message("\n>>> Downloading supplementary files for ", gse_id)

  tryCatch(
    getGEOSuppFiles(
      GEO         = gse_id,
      makeDirectory = TRUE,
      baseDir     = vejiga_conteos_directorio_path,
      fetch_files = TRUE
    ),
    error = function(e) {
      message("  ERROR downloading supplementary for ", gse_id, ": ",
              conditionMessage(e))
      return(NULL)
    }
  )
}

vejiga_suplementarios_resultados_lista <-
  vejiga_datasets_objetivo_tbl$gse_id %>%
  set_names() %>%
  map(descargar_suplementarios_func)


# 06. AUDITORIA DEL REPOSITORIO ________________________________________________

# Comprobamos que cada dataset produjo al menos un ExpressionSet, un CSV de
# metadatos y un archivo suplementario, y reportamos lo que falte.

message("\n>>> Audit of local bladder repository")

vejiga_auditoria_tbl <-
  vejiga_datasets_objetivo_tbl$gse_id %>%
  map_dfr(function(gse_id) {

    eset_archivos_vector <-
      list.files(vejiga_metadatos_directorio_path,
                 pattern = paste0("^", gse_id, ".*ExpressionSet.*\\.rds$"),
                 recursive = TRUE, full.names = TRUE)

    metadatos_archivos_vector <-
      list.files(vejiga_metadatos_directorio_path,
                 pattern = paste0("^", gse_id, ".*metadatos.*\\.csv$"),
                 recursive = TRUE, full.names = TRUE)

    supplementary_archivos_vector <-
      list.files(file.path(vejiga_conteos_directorio_path, gse_id),
                 recursive = TRUE, full.names = TRUE)

    tibble(
      gse_id           = gse_id,
      n_eset_files     = length(eset_archivos_vector),
      n_metadata_files = length(metadatos_archivos_vector),
      n_supp_files     = length(supplementary_archivos_vector),
      supp_total_mb    = round(
        sum(file.info(supplementary_archivos_vector)$size,
            na.rm = TRUE) / 1e6, 2)
    )
  })

print(vejiga_auditoria_tbl)


# 07. RESUMEN DE SESION ________________________________________________________

message("")
message("Bladder dataset download completed")
message("Target datasets attempted: ",
        nrow(vejiga_datasets_objetivo_tbl))
message("Datasets with ExpressionSet retrieved: ",
        sum(vejiga_auditoria_tbl$n_eset_files > 0))
message("Datasets with metadata retrieved: ",
        sum(vejiga_auditoria_tbl$n_metadata_files > 0))
message("Datasets with supplementary files retrieved: ",
        sum(vejiga_auditoria_tbl$n_supp_files > 0))
message("Cumulative supplementary file size: ",
        sum(vejiga_auditoria_tbl$supp_total_mb), " MB")
message("Repository ready for the dataset-by-dataset standardization phase")
