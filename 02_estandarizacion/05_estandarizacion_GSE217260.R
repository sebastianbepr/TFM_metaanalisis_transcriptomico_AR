# ______________________________________________________________________________
# 02_05  Estandarizacion de GSE217260 (prostata)
#
# Que hace: pasa la matriz de conteos crudos y los metadatos de GSE217260 al
#           esquema estandar del proyecto (IDs Ensembl), suma las replicas
#           tecnicas por paciente y lo guarda.
# Entra:    Prostata/GSE217260*metadatos*.csv,
#           Prostata_Counts/GSE217260*matriz_ensamblada.csv
# Sale:     Datasets_Estandarizados/GSE217260_estandarizado.rds
# ______________________________________________________________________________

# __ INDEX ___
# 01. LIBRERIAS
# 02. RUTAS
# 03. LOCALIZAR ARCHIVOS
# 04. IMPORTAR METADATOS CLINICOS
# 05. IMPORTAR MATRIZ DE CONTEOS CRUDOS
# 06. ALINEAR MUESTRAS ENTRE METADATOS Y CONTEOS
# 07. ASIGNAR GRUPO DIAGNOSTICO Y CONSERVAR VARIABLES CLINICAS
# 08. IDENTIFICAR PACIENTES Y REPLICAS TECNICAS
# 09. CONTROL DE CALIDAD DE LOS GENES
# 10. COLAPSAR REPLICAS TECNICAS POR PACIENTE (SUMA)
# 11. ENSAMBLAR LA MATRIZ ESTANDARIZADA
# 12. DESCRIPTOR DEL DATASET
# 13. EXPORTAR EL RDS ESTANDARIZADO
# 14. RESUMEN
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

gse_objetivo_id_char <- "GSE217260"
organo_objetivo_char  <- "Prostate"
plataforma_tipo_char  <- "RNA-seq"


# 03. LOCALIZAR ARCHIVOS _______________________________________________________

message("Auditing dataset ", gse_objetivo_id_char)

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
    pattern    = paste0(gse_objetivo_id_char, ".*matriz_ensamblada\\.csv$"),
    recursive  = TRUE,
    full.names = TRUE
  )

if (length(prostata_metadatos_archivo_path) == 0 ||
    length(prostata_conteos_archivo_path) == 0) {
  message("Metadata files located: ",
          length(prostata_metadatos_archivo_path))
  message("Pre-assembled count files located: ",
          length(prostata_conteos_archivo_path))
  stop("Critical error: incomplete local repository for ",
       gse_objetivo_id_char)
}

message("Metadata file located: ",
        basename(prostata_metadatos_archivo_path[1]))
message("Pre-assembled count matrix located: ",
        basename(prostata_conteos_archivo_path[1]))


# 04. IMPORTAR METADATOS CLINICOS ______________________________________________

prostata_metadatos_tbl_inicial <-
  prostata_metadatos_archivo_path[1] %>%
  fread() %>%
  as_tibble()

message("Clinical metadata imported: ",
        nrow(prostata_metadatos_tbl_inicial), " sequencing libraries, ",
        ncol(prostata_metadatos_tbl_inicial), " annotation columns")

if ("characteristics_ch1" %in% names(prostata_metadatos_tbl_inicial)) {
  message("Native diagnosis distribution (characteristics_ch1):")
  print(table(prostata_metadatos_tbl_inicial$characteristics_ch1))
}


# 05. IMPORTAR MATRIZ DE CONTEOS CRUDOS ________________________________________

prostata_conteos_tbl_inicial <-
  prostata_conteos_archivo_path[1] %>%
  fread() %>%
  as_tibble()

# La primera columna de la matriz ensamblada lleva los identificadores de
# gen pero se exporto con la cabecera vacia; la renombramos de forma
# explicita a 'ensembl_id' para que quede claro mas adelante.

names(prostata_conteos_tbl_inicial)[1] <- "ensembl_id"

message("Count matrix imported: ",
        nrow(prostata_conteos_tbl_inicial), " features x ",
        ncol(prostata_conteos_tbl_inicial) - 1, " sequencing libraries")

prostata_libraries_columnas_vector <-
  setdiff(names(prostata_conteos_tbl_inicial), "ensembl_id")


# 06. ALINEAR MUESTRAS ENTRE METADATOS Y CONTEOS _______________________________

# Esperamos poder cruzar directamente por GSM, porque la matriz ya
# ensamblada conserva los identificadores de GEO en los nombres de columna.

prostata_match_gsm_vector <-
  intersect(
    prostata_libraries_columnas_vector,
    prostata_metadatos_tbl_inicial$geo_accession
  )

message("Direct GSM alignment: ",
        length(prostata_match_gsm_vector), " / ",
        length(prostata_libraries_columnas_vector))

if (length(prostata_match_gsm_vector) ==
    length(prostata_libraries_columnas_vector)) {
  estrategia_alineacion_char <- "GSM accession"
} else {
  estrategia_alineacion_char <- "Partial; manual review required"
  message("Warning: alignment incomplete; manual review required")
}

message("Alignment strategy selected: ", estrategia_alineacion_char)


# 07. ASIGNAR GRUPO DIAGNOSTICO Y CONSERVAR VARIABLES CLINICAS __________________

# Pasamos las etiquetas fenotipicas al vocabulario Tumor / Normal. Como
# en GSE134073, el comparador es hiperplasia benigna de prostata (BPH);
# lo dejamos anotado en 'tejido_control_subtipo' y declarado en el
# descriptor del dataset.

prostata_metadatos_tbl_armonizado <-
  prostata_metadatos_tbl_inicial %>%
  mutate(
    grupo_diagnostico = case_when(
      str_detect(characteristics_ch1,
                 regex("prostate cancer|adenocarcinoma|tumor|carcinoma",
                       ignore_case = TRUE))
      ~ "Tumor",
      str_detect(characteristics_ch1,
                 regex("benign prostat|bph|hyperplasia|normal|adjacent",
                       ignore_case = TRUE))
      ~ "Normal",
      TRUE
      ~ NA_character_
    ),
    tejido_control_subtipo = case_when(
      str_detect(characteristics_ch1,
                 regex("benign prostat|bph|hyperplasia",
                       ignore_case = TRUE))
      ~ "BPH",
      grupo_diagnostico == "Normal"
      ~ "Adjacent normal",
      TRUE
      ~ NA_character_
    ),
    gleason_score = str_remove(characteristics_ch1.2,
                               regex("^gleason score:\\s*",
                                     ignore_case = TRUE)),
    isup_grade_group = str_remove(characteristics_ch1.3,
                                  regex("^isup grade_group:\\s*",
                                        ignore_case = TRUE)),
    proportion_of_tumor = str_remove(characteristics_ch1.4,
                                     regex("^proportion of_tumor:\\s*",
                                           ignore_case = TRUE)),
    edad_paciente = str_remove(characteristics_ch1.1,
                               regex("^age:\\s*",
                                     ignore_case = TRUE))
  )


# 08. IDENTIFICAR PACIENTES Y REPLICAS TECNICAS ________________________________

# Cada title lleva el donante biologico y la replica tecnica en el patron
# '<codigo_paciente><letra_replica> [<diagnostico>]', por ejemplo
# 'A6a [BPH]' o 'E1c [Medium-risk prostate cancer]'. El codigo del
# paciente es el prefijo alfanumerico antes de la letra final de replica.

extraer_paciente_id_func <- function(titulo_text) {
  patron_codigo <- regmatches(
    titulo_text,
    regexpr("^[A-Z]+\\d+", titulo_text)
  )
  ifelse(length(patron_codigo) == 0, NA_character_, patron_codigo)
}

extraer_replica_letra_func <- function(titulo_text) {
  patron_letra <- regmatches(
    titulo_text,
    regexpr("(?<=\\d)[a-z](?=\\s)", titulo_text, perl = TRUE)
  )
  ifelse(length(patron_letra) == 0, NA_character_, patron_letra)
}

prostata_metadatos_tbl_curado <-
  prostata_metadatos_tbl_armonizado %>%
  rowwise() %>%
  mutate(
    paciente_id    = extraer_paciente_id_func(title),
    replica_letra  = extraer_replica_letra_func(title)
  ) %>%
  ungroup()

message("Patient and replicate annotation:")
message("  Unique biological donors: ",
        length(unique(prostata_metadatos_tbl_curado$paciente_id)))
message("  Distribution of technical replicates per patient:")
print(table(prostata_metadatos_tbl_curado$paciente_id,
            prostata_metadatos_tbl_curado$replica_letra))

message("Patient-level diagnostic balance:")
print(prostata_metadatos_tbl_curado %>%
        distinct(paciente_id, grupo_diagnostico) %>%
        count(grupo_diagnostico))


# 09. CONTROL DE CALIDAD DE LOS GENES __________________________________________

# Comprobamos que la matriz trae IDs Ensembl canonicos sin sufijo de
# version y que los valores son conteos crudos enteros compatibles con
# DESeq2.

prostata_ids_genes_vector <- prostata_conteos_tbl_inicial$ensembl_id

prostata_ids_son_ensembl_lgl <-
  mean(grepl("^ENSG\\d+$", head(prostata_ids_genes_vector, 100))) > 0.9

prostata_ids_tienen_version_lgl <-
  any(grepl("\\.\\d+$", head(prostata_ids_genes_vector, 20)))

message("Gene identifier convention: ",
        if (prostata_ids_son_ensembl_lgl) "Canonical Ensembl"
        else "Mixed / non-canonical",
        if (prostata_ids_tienen_version_lgl)
          " (with version suffix)" else " (no version suffix)")

prostata_valores_muestra_matriz <-
  prostata_conteos_tbl_inicial %>%
  select(2:min(6, ncol(.))) %>%
  as.matrix()

prostata_valores_son_enteros_lgl <-
  all(prostata_valores_muestra_matriz ==
        round(prostata_valores_muestra_matriz), na.rm = TRUE)

message("Integer-valued matrix (raw counts): ",
        prostata_valores_son_enteros_lgl)


# 10. COLAPSAR REPLICAS TECNICAS POR PACIENTE (SUMA) ___________________________

# Sumamos los conteos crudos de las replicas tecnicas de cada donante.
# Con datos de conteos esta es la forma correcta de juntar replicas:
# mantiene la estructura de dispersion binomial negativa que asume DESeq2
# y aumenta la profundidad de secuenciacion de la pseudo-libreria final.

prostata_replicas_mapeo_tbl <-
  prostata_metadatos_tbl_curado %>%
  select(geo_accession, paciente_id) %>%
  filter(!is.na(paciente_id))

# Traducimos los nombres de columna de la matriz (GSM) a codigos de paciente
prostata_conteos_largo_tbl <-
  prostata_conteos_tbl_inicial %>%
  pivot_longer(
    cols      = -ensembl_id,
    names_to  = "geo_accession",
    values_to = "raw_count"
  ) %>%
  left_join(prostata_replicas_mapeo_tbl, by = "geo_accession")

prostata_conteos_colapsado_tbl <-
  prostata_conteos_largo_tbl %>%
  group_by(ensembl_id, paciente_id) %>%
  summarise(
    raw_count = sum(raw_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = paciente_id,
    values_from = raw_count
  )

message("Technical replicates collapsed by summation per patient")
message("Resulting matrix: ",
        nrow(prostata_conteos_colapsado_tbl), " genes x ",
        ncol(prostata_conteos_colapsado_tbl) - 1, " patients")


# 11. ENSAMBLAR LA MATRIZ ESTANDARIZADA ________________________________________

# Pasamos la tabla colapsada a una matriz numerica y construimos los
# metadatos por paciente respetando su orden de columnas.

prostata_genes_anotacion_tbl_final <-
  tibble(
    gene_id_native = prostata_conteos_colapsado_tbl$ensembl_id,
    ensembl_id     = prostata_conteos_colapsado_tbl$ensembl_id,
    gene_symbol    = NA_character_
  )

prostata_conteos_matriz_final <-
  prostata_conteos_colapsado_tbl %>%
  select(-ensembl_id) %>%
  as.matrix()

rownames(prostata_conteos_matriz_final) <-
  prostata_genes_anotacion_tbl_final$ensembl_id

# Metadatos por paciente: una fila por donante biologico. Las variables
# clinicas (Gleason, ISUP, proporcion de tumor, edad) se toman de la
# primera replica de cada paciente, porque son atributos del donante que
# no cambian entre replicas.

prostata_metadatos_tbl_final <-
  prostata_metadatos_tbl_curado %>%
  group_by(paciente_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(paciente_id, grupo_diagnostico, tejido_control_subtipo,
         gleason_score, isup_grade_group, proportion_of_tumor,
         edad_paciente) %>%
  rename(sample_id = paciente_id) %>%
  slice(match(colnames(prostata_conteos_matriz_final), sample_id))

message("Standardized matrix assembled: ",
        nrow(prostata_conteos_matriz_final), " genes x ",
        ncol(prostata_conteos_matriz_final), " patients")


# 12. DESCRIPTOR DEL DATASET ___________________________________________________

prostata_descriptor_lista <-
  list(
    gse_id                = gse_objetivo_id_char,
    organo                = organo_objetivo_char,
    plataforma_tipo       = plataforma_tipo_char,
    gene_id_convention    = "Canonical Ensembl (no version suffix)",
    matrix_value_nature   = if (prostata_valores_son_enteros_lgl)
      "Raw integer counts (summed across technical replicates)"
    else "Non-integer values; verify normalization",
    n_total               = ncol(prostata_conteos_matriz_final),
    n_tumor               = sum(prostata_metadatos_tbl_final$grupo_diagnostico
                                == "Tumor", na.rm = TRUE),
    n_normal              = sum(prostata_metadatos_tbl_final$grupo_diagnostico
                                == "Normal", na.rm = TRUE),
    n_sin_etiqueta        = sum(is.na(prostata_metadatos_tbl_final$grupo_diagnostico)),
    comparator_subtype    = "Benign prostatic hyperplasia (BPH)",
    paired_design         = FALSE,
    alignment_strategy    = estrategia_alineacion_char,
    technical_replicates  = "Three replicates per patient summed into one library",
    n_raw_libraries       = length(prostata_libraries_columnas_vector),
    n_patients_final      = ncol(prostata_conteos_matriz_final),
    fecha_estandarizacion = format(Sys.Date(), "%Y-%m-%d")
  )


# 13. EXPORTAR EL RDS ESTANDARIZADO ____________________________________________

prostata_dataset_estandarizado_lista <-
  list(
    conteos_matriz         = prostata_conteos_matriz_final,
    genes_anotacion_tbl    = prostata_genes_anotacion_tbl_final,
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


# 14. RESUMEN __________________________________________________________________

message("")
message("Pipeline completed for ", gse_objetivo_id_char)
message("Organ: ",            prostata_descriptor_lista$organo)
message("Platform: ",         prostata_descriptor_lista$plataforma_tipo)
message("Alignment strategy: ", prostata_descriptor_lista$alignment_strategy)
message("Comparator subtype: ", prostata_descriptor_lista$comparator_subtype)
message("Technical replicates: ",
        prostata_descriptor_lista$technical_replicates)
message("Raw sequencing libraries: ",
        prostata_descriptor_lista$n_raw_libraries,
        " collapsed into ",
        prostata_descriptor_lista$n_patients_final,
        " biological donors")
message("Features retained: ",
        nrow(prostata_dataset_estandarizado_lista$conteos_matriz),
        " canonical Ensembl genes")
message("Patients retained: ",
        prostata_descriptor_lista$n_total,
        " (Tumor = ", prostata_descriptor_lista$n_tumor,
        ", Normal = ", prostata_descriptor_lista$n_normal,
        if (prostata_descriptor_lista$n_sin_etiqueta > 0)
          paste0(", Unlabeled = ", prostata_descriptor_lista$n_sin_etiqueta)
        else "",
        ")")
message("Gene identifier type: ",
        prostata_descriptor_lista$gene_id_convention)
message("Matrix nature: ",
        prostata_descriptor_lista$matrix_value_nature)
message("Output: ", prostata_rds_destino_path)
