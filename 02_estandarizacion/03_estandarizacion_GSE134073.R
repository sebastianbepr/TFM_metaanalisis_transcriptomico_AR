# ______________________________________________________________________________
# 02_03  Estandarizacion de GSE134073 (prostata)
#
# Que hace: pasa la matriz de conteos crudos y los metadatos de GSE134073 al
#           esquema estandar del proyecto (IDs Ensembl) y lo guarda.
# Entra:    Prostata/GSE134073*metadatos*.csv,
#           Prostata_Counts/GSE134073*matriz_ensamblada.csv
# Sale:     Datasets_Estandarizados/GSE134073_estandarizado.rds
# ______________________________________________________________________________

# __ INDEX ___
# 01. LIBRERIAS
# 02. RUTAS
# 03. LOCALIZAR ARCHIVOS
# 04. IMPORTAR METADATOS CLINICOS
# 05. IMPORTAR MATRIZ DE CONTEOS CRUDOS
# 06. ALINEAR MUESTRAS ENTRE METADATOS Y CONTEOS
# 07. ASIGNAR GRUPO DIAGNOSTICO Y CONSERVAR VARIABLES CLINICAS
# 08. CONTROL DE CALIDAD DE LOS GENES
# 09. FILTRAR GENES ENSEMBL CANONICOS Y QUITAR SUFIJO DE VERSION
# 10. COLAPSAR GENES DUPLICADOS POR MAXIMA EXPRESION MEDIA
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

gse_objetivo_id_char <- "GSE134073"
organo_objetivo_char  <- "Prostate"
plataforma_tipo_char  <- "RNA-seq"


# 03. LOCALIZAR ARCHIVOS _______________________________________________________

# Buscamos de forma recursiva el archivo de metadatos y la matriz de
# conteos ya ensamblada. Preferimos la matriz ensamblada a los 64
# archivos por muestra estilo HTSeq porque evita rehacer el ensamblado y
# la informacion biologica es la misma.

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
        nrow(prostata_metadatos_tbl_inicial), " samples, ",
        ncol(prostata_metadatos_tbl_inicial), " annotation columns")

# Miramos el campo de diagnostico (en characteristics_ch1.1 para esta
# serie de GEO), que distingue el adenocarcinoma de prostata del grupo
# comparador de hiperplasia benigna.

if ("characteristics_ch1.1" %in% names(prostata_metadatos_tbl_inicial)) {
  message("Native diagnosis distribution (characteristics_ch1.1):")
  print(table(prostata_metadatos_tbl_inicial$characteristics_ch1.1))
}


# 05. IMPORTAR MATRIZ DE CONTEOS CRUDOS ________________________________________

prostata_conteos_tbl_inicial <-
  prostata_conteos_archivo_path[1] %>%
  fread() %>%
  as_tibble()

message("Count matrix imported: ",
        nrow(prostata_conteos_tbl_inicial), " features x ",
        ncol(prostata_conteos_tbl_inicial), " columns")

prostata_id_genes_columna_char <-
  names(prostata_conteos_tbl_inicial)[1]

prostata_muestras_columnas_vector <-
  names(prostata_conteos_tbl_inicial)[-1]

message("Feature identifier column: '",
        prostata_id_genes_columna_char, "'")
message("Biological sample columns: ",
        length(prostata_muestras_columnas_vector))


# 06. ALINEAR MUESTRAS ENTRE METADATOS Y CONTEOS _______________________________

# Aqui esperamos poder cruzar directamente por GSM, porque la matriz ya
# ensamblada conserva los identificadores de GEO en las columnas.

prostata_match_gsm_vector <-
  intersect(
    prostata_muestras_columnas_vector,
    prostata_metadatos_tbl_inicial$geo_accession
  )

message("Direct GSM alignment: ",
        length(prostata_match_gsm_vector), " / ",
        length(prostata_muestras_columnas_vector))

if (length(prostata_match_gsm_vector) ==
    length(prostata_muestras_columnas_vector)) {
  estrategia_alineacion_char <- "GSM accession"
  prostata_metadatos_tbl_armonizado <-
    prostata_metadatos_tbl_inicial %>%
    mutate(id_alineacion = geo_accession)
} else {
  estrategia_alineacion_char <- "Manual review required"
  prostata_metadatos_tbl_armonizado <-
    prostata_metadatos_tbl_inicial %>%
    mutate(id_alineacion = geo_accession)
  message("Warning: alignment incomplete; manual review required")
}

message("Alignment strategy selected: ", estrategia_alineacion_char)


# 07. ASIGNAR GRUPO DIAGNOSTICO Y CONSERVAR VARIABLES CLINICAS __________________

# Pasamos las etiquetas fenotipicas al vocabulario Tumor / Normal comun
# a todos los datasets. En GSE134073 el comparador es hiperplasia benigna
# de prostata (BPH); lo dejamos anotado en 'tejido_control_subtipo' para
# los analisis de sensibilidad y para declararlo en la memoria.

prostata_metadatos_tbl_curado <-
  prostata_metadatos_tbl_armonizado %>%
  mutate(
    grupo_diagnostico = case_when(
      str_detect(characteristics_ch1.1,
                 regex("prostate cancer|adenocarcinoma|tumor|carcinoma",
                       ignore_case = TRUE))
      ~ "Tumor",
      str_detect(characteristics_ch1.1,
                 regex("benign prostat|bph|hyperplasia|normal|adjacent",
                       ignore_case = TRUE))
      ~ "Normal",
      TRUE
      ~ NA_character_
    ),
    tejido_control_subtipo = case_when(
      str_detect(characteristics_ch1.1,
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
    risk_group = str_remove(characteristics_ch1.3,
                            regex("^risk group:\\s*",
                                  ignore_case = TRUE))
  )

message("Standardized phenotypic distribution:")
print(table(prostata_metadatos_tbl_curado$grupo_diagnostico,
            prostata_metadatos_tbl_curado$tejido_control_subtipo,
            useNA = "ifany"))

# Construimos un tibble que asocia cada columna de la matriz con su
# anotacion clinica, respetando el orden de columnas de la matriz.

prostata_grupos_alineados_tbl <-
  tibble(sample_id_matriz = prostata_muestras_columnas_vector) %>%
  mutate(
    grupo_diagnostico = prostata_metadatos_tbl_curado$grupo_diagnostico[
      match(sample_id_matriz, prostata_metadatos_tbl_curado$id_alineacion)
    ],
    tejido_control_subtipo = prostata_metadatos_tbl_curado$tejido_control_subtipo[
      match(sample_id_matriz, prostata_metadatos_tbl_curado$id_alineacion)
    ],
    gleason_score = prostata_metadatos_tbl_curado$gleason_score[
      match(sample_id_matriz, prostata_metadatos_tbl_curado$id_alineacion)
    ],
    risk_group = prostata_metadatos_tbl_curado$risk_group[
      match(sample_id_matriz, prostata_metadatos_tbl_curado$id_alineacion)
    ],
    sample_title = prostata_metadatos_tbl_curado$title[
      match(sample_id_matriz, prostata_metadatos_tbl_curado$id_alineacion)
    ],
    gsm_accession = prostata_metadatos_tbl_curado$geo_accession[
      match(sample_id_matriz, prostata_metadatos_tbl_curado$id_alineacion)
    ]
  )

message("Final diagnostic balance:")
print(table(prostata_grupos_alineados_tbl$grupo_diagnostico,
            useNA = "ifany"))


# 08. CONTROL DE CALIDAD DE LOS GENES __________________________________________

# Miramos como es la columna de identificadores. GSE134073 se proceso con
# un pipeline TopHat + Cufflinks que introduce identificadores de
# transcritos nuevos (prefijo XLOC) junto a los IDs Ensembl canonicos.

prostata_ids_genes_vector <-
  prostata_conteos_tbl_inicial[[prostata_id_genes_columna_char]]

prostata_ids_ensembl_lgl_vector <-
  grepl("^ENSG\\d+", prostata_ids_genes_vector)

prostata_ids_xloc_lgl_vector <-
  grepl("^XLOC", prostata_ids_genes_vector)

prostata_ids_otros_lgl_vector <-
  !prostata_ids_ensembl_lgl_vector & !prostata_ids_xloc_lgl_vector

message("Feature composition:")
message("  Canonical Ensembl (ENSG): ",
        sum(prostata_ids_ensembl_lgl_vector))
message("  Cufflinks novel transcripts (XLOC): ",
        sum(prostata_ids_xloc_lgl_vector))
message("  Other / non-canonical: ",
        sum(prostata_ids_otros_lgl_vector))

# Comprobamos que son conteos antes de filtrar nada.

prostata_valores_muestra_matriz <-
  prostata_conteos_tbl_inicial %>%
  select(2:min(6, ncol(.))) %>%
  as.matrix()

prostata_valores_son_enteros_lgl <-
  all(prostata_valores_muestra_matriz ==
        round(prostata_valores_muestra_matriz), na.rm = TRUE)

message("Integer-valued matrix (raw counts): ",
        prostata_valores_son_enteros_lgl)


# 09. FILTRAR GENES ENSEMBL CANONICOS Y QUITAR SUFIJO DE VERSION ________________

# Nos quedamos solo con los genes Ensembl canonicos (prefijo ENSG),
# porque los transcritos nuevos XLOC no son comparables con datasets que
# usan anotaciones de referencia estandar. Quitamos el sufijo de version
# (.10, .5, etc.) para poder cruzar despues con anotaciones de referencia.

prostata_conteos_tbl_02 <-
  prostata_conteos_tbl_inicial %>%
  filter(grepl("^ENSG\\d+", .data[[prostata_id_genes_columna_char]])) %>%
  mutate(
    ensembl_id = str_remove(.data[[prostata_id_genes_columna_char]],
                            "\\..*$")
  ) %>%
  relocate(ensembl_id, .before = 1) %>%
  select(-all_of(prostata_id_genes_columna_char))

message("Canonical Ensembl records retained: ",
        nrow(prostata_conteos_tbl_02))
message("Unique Ensembl IDs after version stripping: ",
        length(unique(prostata_conteos_tbl_02$ensembl_id)))


# 10. COLAPSAR GENES DUPLICADOS POR MAXIMA EXPRESION MEDIA ______________________

# Cuando varios registros caen en el mismo gen tras quitar la version,
# nos quedamos con el de mayor expresion media entre muestras. Es la
# regla estandar de colapso en meta-analisis y conserva la senal mas
# informativa por gen.

prostata_conteos_tbl_curado <-
  prostata_conteos_tbl_02 %>%
  mutate(
    expresion_media = rowMeans(
      select(., -ensembl_id),
      na.rm = TRUE
    )
  ) %>%
  arrange(ensembl_id, desc(expresion_media)) %>%
  distinct(ensembl_id, .keep_all = TRUE) %>%
  select(-expresion_media)

message("Genes retained after duplicate collapse: ",
        nrow(prostata_conteos_tbl_curado))


# 11. ENSAMBLAR LA MATRIZ ESTANDARIZADA ________________________________________

# Pasamos la tabla de conteos ya curada a una matriz numerica con los
# IDs Ensembl como nombres de fila y las muestras como columnas.

prostata_genes_anotacion_tbl_final <-
  tibble(
    gene_id_native = paste0(prostata_conteos_tbl_curado$ensembl_id),
    ensembl_id     = prostata_conteos_tbl_curado$ensembl_id,
    gene_symbol    = NA_character_
  )

prostata_conteos_matriz_final <-
  prostata_conteos_tbl_curado %>%
  select(-ensembl_id) %>%
  as.matrix()

rownames(prostata_conteos_matriz_final) <-
  prostata_genes_anotacion_tbl_final$ensembl_id

# Reordenamos los metadatos para que sigan el orden de columnas de la matriz.

prostata_metadatos_tbl_final <-
  prostata_grupos_alineados_tbl %>%
  rename(sample_id = sample_id_matriz) %>%
  mutate(paciente_id = NA_character_) %>%
  select(sample_id, gsm_accession, sample_title,
         grupo_diagnostico, tejido_control_subtipo,
         gleason_score, risk_group, paciente_id) %>%
  slice(match(colnames(prostata_conteos_matriz_final), sample_id))

message("Standardized matrix assembled: ",
        nrow(prostata_conteos_matriz_final), " genes x ",
        ncol(prostata_conteos_matriz_final), " samples")


# 12. DESCRIPTOR DEL DATASET ___________________________________________________

# Descriptor autoexplicativo que recoge la procedencia biologica y las
# decisiones de preprocesado de este dataset, incluida la declaracion del
# comparador BPH y la regla de filtrado de Cufflinks.

prostata_descriptor_lista <-
  list(
    gse_id              = gse_objetivo_id_char,
    organo              = organo_objetivo_char,
    plataforma_tipo     = plataforma_tipo_char,
    gene_id_convention  = "Ensembl (version-stripped)",
    matrix_value_nature = if (prostata_valores_son_enteros_lgl)
      "Raw integer counts"
    else "Non-integer values; verify normalization",
    n_total             = ncol(prostata_conteos_matriz_final),
    n_tumor             = sum(prostata_metadatos_tbl_final$grupo_diagnostico
                              == "Tumor", na.rm = TRUE),
    n_normal            = sum(prostata_metadatos_tbl_final$grupo_diagnostico
                              == "Normal", na.rm = TRUE),
    n_sin_etiqueta      = sum(is.na(prostata_metadatos_tbl_final$grupo_diagnostico)),
    comparator_subtype  = "Benign prostatic hyperplasia (BPH)",
    paired_design       = FALSE,
    alignment_strategy  = estrategia_alineacion_char,
    feature_filtering   = "Retained canonical ENSG records; XLOC novel transcripts discarded",
    duplicate_collapse  = "Maximum mean expression per Ensembl ID",
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
message("Features retained: ",
        nrow(prostata_dataset_estandarizado_lista$conteos_matriz),
        " canonical Ensembl genes")
message("Samples retained: ",
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
message("Feature filtering: ",
        prostata_descriptor_lista$feature_filtering)
message("Duplicate collapse rule: ",
        prostata_descriptor_lista$duplicate_collapse)
message("Output: ", prostata_rds_destino_path)
