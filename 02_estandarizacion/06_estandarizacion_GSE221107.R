# ______________________________________________________________________________
# 02_06  Estandarizacion de GSE221107 (prostata, pareado tumor/normal)
#
# Que hace: convierte los conteos crudos y los metadatos de GSE221107 al esquema
#           estandar del proyecto (simbolo HGNC, descriptor comun) y lo guarda.
# Entra:    Prostata/GSE221107*metadatos*.csv,
#           Prostata_Counts/GSE221107*raw_counts*.csv.gz
# Sale:     Datasets_Estandarizados/GSE221107_estandarizado.rds
# ______________________________________________________________________________

# __ INDEX ___
# 01. LIBRERIAS
# 02. RUTAS
# 03. IMPORTAR LOS DATOS CRUDOS VERIFICADOS
# 04. ARMONIZAR ANOTACIONES CLINICAS
# 05. ESTRUCTURAR LA ANOTACION DE GENES
# 06. ENSAMBLAR LA MATRIZ ESTANDARIZADA
# 07. DESCRIPTOR DEL DATASET
# 08. EXPORTAR EL RDS ESTANDARIZADO
# 09. RESUMEN
# ________________________________________________________________________


# 01. LIBRERIAS ________________________________________________________________

library(tidyverse)
library(data.table)


# 02. RUTAS ____________________________________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
base_datos_directorio_path <- ruta   # el resto del script usa este nombre sin cambios

prostata_metadatos_directorio_path <-
  file.path(base_datos_directorio_path, "Prostata")

prostata_conteos_directorio_path <-
  file.path(base_datos_directorio_path, "Prostata_Counts")

datasets_estandarizados_directorio_path <-
  file.path(base_datos_directorio_path, "Datasets_Estandarizados")

if (!dir.exists(datasets_estandarizados_directorio_path)) {
  dir.create(datasets_estandarizados_directorio_path, recursive = TRUE)
}

gse_objetivo_id_char <- "GSE221107"
organo_objetivo_char  <- "Prostate"
plataforma_tipo_char  <- "RNA-seq"


# 03. IMPORTAR LOS DATOS CRUDOS VERIFICADOS ____________________________________

# Buscamos de forma recursiva los dos archivos clave del dataset, ya que
# cada serie de GEO cuelga de su propio subdirectorio.

message("Standardizing dataset ", gse_objetivo_id_char)

prostata_metadatos_archivo_path <-
  list.files(
    path       = prostata_metadatos_directorio_path,
    pattern    = paste0(gse_objetivo_id_char, ".*metadatos.*\\.csv$"),
    recursive  = TRUE,
    full.names = TRUE
  )

prostata_conteos_archivo_path <-
  list.files(
    path       = prostata_conteos_directorio_path,
    pattern    = paste0(gse_objetivo_id_char, ".*raw_counts.*\\.csv\\.gz$"),
    recursive  = TRUE,
    full.names = TRUE
  )

prostata_metadatos_tbl_inicial <-
  prostata_metadatos_archivo_path[1] %>%
  fread() %>%
  as_tibble()

prostata_conteos_tbl_inicial <-
  prostata_conteos_archivo_path[1] %>%
  fread() %>%
  as_tibble()

message("Raw inputs imported: ",
        nrow(prostata_metadatos_tbl_inicial), " samples, ",
        nrow(prostata_conteos_tbl_inicial), " genes")


# 04. ARMONIZAR ANOTACIONES CLINICAS ___________________________________________

# Recuperamos el identificador corto de muestra que va dentro del campo
# 'title' de GEO, que codifica el diseno pareado tumor / normal y es la
# clave de cruce con la matriz de conteos.

extraer_codigo_muestra_func <- function(titulo_text) {
  patron_codigo <- regmatches(titulo_text,
                              regexpr("\\(P[CN]\\d+\\)", titulo_text))
  gsub("[()]", "", patron_codigo)
}

# Pasamos las variables clinicas a un vocabulario comun a todos los
# datasets del meta-analisis. El grupo diagnostico es el contraste
# biologico central (Tumor vs Normal); el codigo de paciente conserva el
# diseno pareado tumor-adyacente de esta cohorte.

prostata_metadatos_tbl_curado <-
  prostata_metadatos_tbl_inicial %>%
  mutate(
    sample_id_corto = extraer_codigo_muestra_func(title),
    grupo_diagnostico = case_when(
      str_detect(source_name_ch1, regex("cancer|tumor", ignore_case = TRUE))
      ~ "Tumor",
      str_detect(source_name_ch1, regex("normal|adjacent", ignore_case = TRUE))
      ~ "Normal",
      TRUE
      ~ NA_character_
    ),
    paciente_id = str_remove(sample_id_corto, "^P[CN]")
  ) %>%
  select(
    sample_id      = sample_id_corto,
    gsm_accession  = geo_accession,
    sample_title   = title,
    grupo_diagnostico,
    paciente_id
  )

message("Clinical annotations harmonized")
print(table(prostata_metadatos_tbl_curado$grupo_diagnostico))


# 05. ESTRUCTURAR LA ANOTACION DE GENES ________________________________________

# Localizamos la columna de identificadores de gen y montamos un tibble
# de anotacion minimo. En GSE221107 los identificadores ya son simbolos
# HGNC, asi que dejamos el campo Ensembl a NA a proposito; solo se rellena
# en los datasets cuyo identificador nativo es un accession de Ensembl.

prostata_id_genes_columna_char <-
  names(prostata_conteos_tbl_inicial)[1]

prostata_genes_anotacion_tbl_inicial <-
  tibble(
    gene_id_native = prostata_conteos_tbl_inicial[[prostata_id_genes_columna_char]]
  ) %>%
  mutate(
    gene_symbol   = gene_id_native,
    ensembl_id    = NA_character_
  )

message("Gene annotation space structured: ",
        nrow(prostata_genes_anotacion_tbl_inicial), " features")


# 06. ENSAMBLAR LA MATRIZ ESTANDARIZADA ________________________________________

# Pasamos la tabla de conteos a una matriz numerica con los genes en
# filas y las muestras en columnas, indexada por el identificador corto
# de muestra que cruza con los metadatos ya armonizados.

prostata_conteos_matriz_final <-
  prostata_conteos_tbl_inicial %>%
  select(-all_of(prostata_id_genes_columna_char)) %>%
  as.matrix()

rownames(prostata_conteos_matriz_final) <-
  prostata_genes_anotacion_tbl_inicial$gene_symbol

# Comprobamos que todas las columnas de la matriz tienen anotacion clinica
# antes de exportar. Una columna sin cruzar meteria muestras sin sentido
# biologico en el meta-analisis.

prostata_columnas_no_anotadas_vector <-
  setdiff(
    colnames(prostata_conteos_matriz_final),
    prostata_metadatos_tbl_curado$sample_id
  )

if (length(prostata_columnas_no_anotadas_vector) > 0) {
  message("Warning: unmatched matrix columns: ",
          paste(prostata_columnas_no_anotadas_vector, collapse = ", "))
  stop("Cannot export with unannotated samples")
}

# Reordenamos los metadatos para que sigan el orden de columnas de la
# matriz, manteniendo la concordancia posicional que necesitan las
# matrices de diseno de DESeq2.

prostata_metadatos_tbl_final <-
  prostata_metadatos_tbl_curado %>%
  slice(match(colnames(prostata_conteos_matriz_final), sample_id))

message("Standardized matrix assembled: ",
        nrow(prostata_conteos_matriz_final), " genes x ",
        ncol(prostata_conteos_matriz_final), " samples")


# 07. DESCRIPTOR DEL DATASET ___________________________________________________

# Descriptor compacto que documenta la procedencia biologica y las
# decisiones de preprocesado del dataset, para que cada contenedor .rds
# quede autodocumentado de cara a las auditorias de reproducibilidad.

prostata_descriptor_lista <-
  list(
    gse_id              = gse_objetivo_id_char,
    organo              = organo_objetivo_char,
    plataforma_tipo     = plataforma_tipo_char,
    gene_id_convention  = "Gene symbol",
    matrix_value_nature = "Raw integer counts",
    n_total             = ncol(prostata_conteos_matriz_final),
    n_tumor             = sum(prostata_metadatos_tbl_final$grupo_diagnostico
                              == "Tumor"),
    n_normal            = sum(prostata_metadatos_tbl_final$grupo_diagnostico
                              == "Normal"),
    paired_design       = TRUE,
    fecha_estandarizacion = format(Sys.Date(), "%Y-%m-%d")
  )


# 08. EXPORTAR EL RDS ESTANDARIZADO ____________________________________________

prostata_dataset_estandarizado_lista <-
  list(
    conteos_matriz         = prostata_conteos_matriz_final,
    genes_anotacion_tbl    = prostata_genes_anotacion_tbl_inicial,
    muestras_metadatos_tbl = prostata_metadatos_tbl_final,
    dataset_descriptor     = prostata_descriptor_lista
  )

prostata_rds_destino_path <-
  file.path(
    datasets_estandarizados_directorio_path,
    paste0(gse_objetivo_id_char, "_estandarizado.rds")
  )

saveRDS(prostata_dataset_estandarizado_lista, prostata_rds_destino_path)

message("Standardized object exported to: ",
        basename(prostata_rds_destino_path))


# 09. RESUMEN __________________________________________________________________

message("")
message("Standardization completed for ", gse_objetivo_id_char)
message("Organ: ",            prostata_descriptor_lista$organo)
message("Platform: ",         prostata_descriptor_lista$plataforma_tipo)
message("Features retained: ",
        nrow(prostata_dataset_estandarizado_lista$conteos_matriz))
message("Samples retained: ",
        prostata_descriptor_lista$n_total,
        " (Tumor = ", prostata_descriptor_lista$n_tumor,
        ", Normal = ", prostata_descriptor_lista$n_normal, ")")
message("Paired tumor-adjacent design: ",
        prostata_descriptor_lista$paired_design)
message("Output: ", prostata_rds_destino_path)
