# ______________________________________________________________________________
# 02_01  Estandarizacion de GSE104131 (prostata, pareado tumor/normal)
#
# Que hace: pasa la matriz FPKM de Cufflinks y los metadatos de GSE104131 al
#           esquema estandar del proyecto (simbolo HGNC, log2) y lo guarda.
# Entra:    Prostata/GSE104131*metadatos*.csv,
#           Prostata_Counts/GSE104131*batch1*fpkms_tracking*.txt.gz
# Sale:     Datasets_Estandarizados/GSE104131_estandarizado.rds
# ______________________________________________________________________________

# __ INDEX ___
# 01. LIBRERIAS
# 02. RUTAS
# 03. LOCALIZAR ARCHIVOS
# 04. IMPORTAR METADATOS CLINICOS
# 05. IMPORTAR MATRIZ FPKM CUFFLINKS (BATCH1)
# 06. ARMONIZAR IDENTIFICADORES Y CORREGIR ERROR DE TRANSCRIPCION
# 07. ASIGNAR GRUPO DIAGNOSTICO E IDENTIFICAR PACIENTES PAREADOS
# 08. CONTROL DE CALIDAD DE LOS GENES
# 09. TRANSFORMAR A LOG2 Y ALINEAR MATRIZ
# 10. CONSTRUIR ANOTACION DE GENES
# 11. COLAPSAR SIMBOLOS DUPLICADOS
# 12. ENSAMBLAR LA MATRIZ ESTANDARIZADA
# 13. DESCRIPTOR DEL DATASET
# 14. EXPORTAR EL RDS ESTANDARIZADO
# 15. RESUMEN
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

gse_objetivo_id_char  <- "GSE104131"
organo_objetivo_char  <- "Prostate"
plataforma_tipo_char  <- "RNA-seq (Cufflinks FPKM, gene-level)"


# 03. LOCALIZAR ARCHIVOS _______________________________________________________

message("Auditing dataset ", gse_objetivo_id_char)

prostata_metadatos_archivo_path <-
  list.files(
    path       = prostata_metadatos_directorio_path,
    pattern    = paste0(gse_objetivo_id_char, ".*metadatos.*\\.csv$"),
    recursive  = TRUE,
    full.names = TRUE
  )

prostata_batch1_archivo_path <-
  list.files(
    path       = prostata_conteos_directorio_path,
    pattern    = paste0(gse_objetivo_id_char,
                        ".*batch1.*fpkms_tracking.*\\.txt\\.gz$"),
    recursive  = TRUE,
    full.names = TRUE
  )

if (length(prostata_metadatos_archivo_path) == 0 ||
    length(prostata_batch1_archivo_path) == 0) {
  message("Metadata files located: ",
          length(prostata_metadatos_archivo_path))
  message("Cufflinks batch1 files located: ",
          length(prostata_batch1_archivo_path))
  stop("Critical error: incomplete local repository for ",
       gse_objetivo_id_char)
}

message("Metadata file located: ",
        basename(prostata_metadatos_archivo_path[1]))
message("Cufflinks batch1 file located: ",
        basename(prostata_batch1_archivo_path[1]))


# 04. IMPORTAR METADATOS CLINICOS ______________________________________________

prostata_metadatos_tbl_inicial <-
  prostata_metadatos_archivo_path[1] %>%
  fread() %>%
  as_tibble()

message("Clinical metadata imported: ",
        nrow(prostata_metadatos_tbl_inicial), " samples, ",
        ncol(prostata_metadatos_tbl_inicial), " annotation columns")

message("Native characteristics_ch1 distribution:")
print(table(prostata_metadatos_tbl_inicial$characteristics_ch1))

message("Native ethnicity distribution (characteristics_ch1.1):")
print(table(prostata_metadatos_tbl_inicial$characteristics_ch1.1))


# 05. IMPORTAR MATRIZ FPKM CUFFLINKS (BATCH1) __________________________________

prostata_fpkm_tbl_inicial <-
  prostata_batch1_archivo_path[1] %>%
  fread() %>%
  as_tibble()

message("FPKM matrix imported: ",
        nrow(prostata_fpkm_tbl_inicial), " features x ",
        ncol(prostata_fpkm_tbl_inicial) - 3, " sequencing libraries")

prostata_columnas_anotacion_vector <- c("GID", "GENESYMBOL", "Locus")

prostata_libraries_columnas_vector <-
  setdiff(names(prostata_fpkm_tbl_inicial),
          prostata_columnas_anotacion_vector)


# 06. ARMONIZAR IDENTIFICADORES Y CORREGIR ERROR DE TRANSCRIPCION ______________

# La matriz nombra las muestras con codigos como '1014_C17', pero los
# metadatos guardan el mismo identificador con el prefijo 'Sample_'
# (por ejemplo 'Sample_1014_C17'). Quitamos el prefijo para poder cruzarlos.
# Ademas hay un error de transcripcion en el paciente 881: los metadatos
# ponen '881_819' / '881_826' y la matriz pone '881_B19' / '881_B26'
# (la letra B se escribio como el numero 8). Lo corregimos aqui para
# recuperar el alineamiento completo 32/32.

prostata_metadatos_tbl_armonizado <-
  prostata_metadatos_tbl_inicial %>%
  mutate(
    sample_id_matriz_corregido = title %>%
      str_remove("^Sample_") %>%
      str_replace("^881_819$", "881_B19") %>%
      str_replace("^881_826$", "881_B26")
  )

prostata_typos_corregidos_int <-
  sum(prostata_metadatos_tbl_armonizado$sample_id_matriz_corregido !=
        str_remove(prostata_metadatos_tbl_armonizado$title, "^Sample_"))

message("Transcription typos corrected (patient 881): ",
        prostata_typos_corregidos_int)

prostata_match_corregido_vector <-
  intersect(
    prostata_libraries_columnas_vector,
    prostata_metadatos_tbl_armonizado$sample_id_matriz_corregido
  )

message("Sample alignment after typo correction: ",
        length(prostata_match_corregido_vector), " / ",
        length(prostata_libraries_columnas_vector))

if (length(prostata_match_corregido_vector) !=
    length(prostata_libraries_columnas_vector)) {
  stop("Critical error: alignment incomplete despite typo correction")
}

estrategia_alineacion_char <-
  "Title field with 'Sample_' prefix stripped and patient 881 typo corrected"


# 07. ASIGNAR GRUPO DIAGNOSTICO E IDENTIFICAR PACIENTES PAREADOS ________________

# El diagnostico viene en el campo 'source_name_ch1' como 'tumor N' o
# 'normal tissue matched to tumor N', donde N es el indice del paciente
# pareado. Extraemos ese indice y lo usamos como paciente_id para enlazar
# cada tumor con su normal adyacente.

extraer_paciente_id_func <- function(texto) {
  patron <- regmatches(texto,
                       regexpr("(?<=tumor\\s)\\d+", texto, perl = TRUE))
  ifelse(length(patron) == 0, NA_character_, patron)
}

prostata_metadatos_tbl_curado <-
  prostata_metadatos_tbl_armonizado %>%
  rowwise() %>%
  mutate(
    paciente_id = extraer_paciente_id_func(source_name_ch1)
  ) %>%
  ungroup() %>%
  mutate(
    grupo_diagnostico = case_when(
      str_detect(source_name_ch1,
                 regex("^tumor\\s", ignore_case = TRUE))
      ~ "Tumor",
      str_detect(source_name_ch1,
                 regex("normal|matched", ignore_case = TRUE))
      ~ "Normal",
      TRUE
      ~ NA_character_
    ),
    tejido_control_subtipo = case_when(
      grupo_diagnostico == "Normal"
      ~ "Adjacent normal (patient-matched)",
      TRUE
      ~ NA_character_
    ),
    ethnicity_status = str_remove(characteristics_ch1.1,
                                  regex("^ethnicity:\\s*",
                                        ignore_case = TRUE))
  )

message("Standardized phenotypic distribution:")
print(table(prostata_metadatos_tbl_curado$grupo_diagnostico))

message("Paired patient structure:")
print(prostata_metadatos_tbl_curado %>%
        count(paciente_id, grupo_diagnostico) %>%
        pivot_wider(names_from = grupo_diagnostico,
                    values_from = n,
                    values_fill = 0))

message("Ethnicity distribution within diagnostic groups:")
print(table(prostata_metadatos_tbl_curado$grupo_diagnostico,
            prostata_metadatos_tbl_curado$ethnicity_status))


# 08. CONTROL DE CALIDAD DE LOS GENES __________________________________________

# Anotamos que tipo de identificador trae GSE104131: da a la vez IDs
# Ensembl canonicos (columna GID) y simbolos HGNC (GENESYMBOL).

prostata_ids_ensembl_vector <-
  prostata_fpkm_tbl_inicial[["GID"]]

prostata_ids_symbol_vector <-
  prostata_fpkm_tbl_inicial[["GENESYMBOL"]]

prostata_ids_son_ensembl_lgl <-
  mean(grepl("^ENSG\\d+$",
             head(prostata_ids_ensembl_vector, 100))) > 0.9

prostata_ids_tienen_version_lgl <-
  any(grepl("\\.\\d+$", head(prostata_ids_ensembl_vector, 20)))

message("Gene identifier convention: ",
        if (prostata_ids_son_ensembl_lgl) "Canonical Ensembl"
        else "Mixed / non-canonical",
        if (prostata_ids_tienen_version_lgl) " (with version suffix)"
        else " (no version suffix)")

message("Symbol coverage: ",
        sum(!is.na(prostata_ids_symbol_vector) &
              prostata_ids_symbol_vector != ""),
        " / ", length(prostata_ids_symbol_vector))

# Miramos el rango de valores FPKM.

prostata_valores_muestra_matriz <-
  prostata_fpkm_tbl_inicial %>%
  select(all_of(prostata_libraries_columnas_vector[
    1:min(5, length(prostata_libraries_columnas_vector))
  ])) %>%
  as.matrix()

message("FPKM value range across first 5 samples:")
message("  min = ",
        round(min(prostata_valores_muestra_matriz, na.rm = TRUE), 3))
message("  max = ",
        round(max(prostata_valores_muestra_matriz, na.rm = TRUE), 3))
message("  median = ",
        round(median(prostata_valores_muestra_matriz, na.rm = TRUE), 3))

prostata_valores_tienen_negativos_lgl <-
  any(prostata_valores_muestra_matriz < 0, na.rm = TRUE)

message("Negative values present: ", prostata_valores_tienen_negativos_lgl)


# 09. TRANSFORMAR A LOG2 Y ALINEAR MATRIZ ______________________________________

# Aplicamos log2(FPKM + 1) para dejar los datos en la escala que espera
# la rama de limma del analisis de expresion diferencial.

prostata_expresion_matriz_log2 <-
  prostata_fpkm_tbl_inicial %>%
  select(all_of(prostata_libraries_columnas_vector)) %>%
  as.matrix() %>%
  pmax(0) %>%
  {log2(. + 1)}

message("Log2(FPKM + 1) transformation applied")
message("Transformed value range:")
message("  min = ",
        round(min(prostata_expresion_matriz_log2, na.rm = TRUE), 3))
message("  max = ",
        round(max(prostata_expresion_matriz_log2, na.rm = TRUE), 3))
message("  median = ",
        round(median(prostata_expresion_matriz_log2, na.rm = TRUE), 3))


# 10. CONSTRUIR ANOTACION DE GENES _____________________________________________

# Nos quedamos solo con los genes que tienen simbolo HGNC valido, porque
# el meta-analisis trabaja sobre simbolos y los genes sin mapear no se
# pueden integrar con el resto de datasets.

prostata_expresion_tbl_inicial <-
  bind_cols(
    tibble(
      ensembl_id  = prostata_ids_ensembl_vector,
      gene_symbol = prostata_ids_symbol_vector
    ),
    as_tibble(prostata_expresion_matriz_log2)
  ) %>%
  mutate(
    gene_symbol = na_if(str_trim(gene_symbol), "")
  ) %>%
  filter(!is.na(gene_symbol))

message("Features with valid HGNC symbol: ",
        nrow(prostata_expresion_tbl_inicial), " / ",
        nrow(prostata_fpkm_tbl_inicial))


# 11. COLAPSAR SIMBOLOS DUPLICADOS _____________________________________________

# Cuando varios registros Ensembl caen en el mismo simbolo HGNC, nos
# quedamos con el de mayor expresion media entre muestras.

prostata_expresion_tbl_curado <-
  prostata_expresion_tbl_inicial %>%
  mutate(
    expresion_media = rowMeans(
      select(., -ensembl_id, -gene_symbol),
      na.rm = TRUE
    )
  ) %>%
  arrange(gene_symbol, desc(expresion_media)) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  select(-expresion_media)

message("Genes retained after duplicate collapse: ",
        nrow(prostata_expresion_tbl_curado))


# 12. ENSAMBLAR LA MATRIZ ESTANDARIZADA ________________________________________

prostata_genes_anotacion_tbl_final <-
  tibble(
    gene_id_native = prostata_expresion_tbl_curado$ensembl_id,
    ensembl_id     = prostata_expresion_tbl_curado$ensembl_id,
    gene_symbol    = prostata_expresion_tbl_curado$gene_symbol
  )

prostata_expresion_matriz_final <-
  prostata_expresion_tbl_curado %>%
  select(-ensembl_id, -gene_symbol) %>%
  as.matrix()

rownames(prostata_expresion_matriz_final) <-
  prostata_genes_anotacion_tbl_final$gene_symbol

prostata_metadatos_tbl_final <-
  prostata_metadatos_tbl_curado %>%
  select(sample_id = sample_id_matriz_corregido,
         gsm_accession = geo_accession,
         paciente_id,
         grupo_diagnostico,
         tejido_control_subtipo,
         ethnicity_status) %>%
  slice(match(colnames(prostata_expresion_matriz_final), sample_id))

message("Standardized matrix assembled: ",
        nrow(prostata_expresion_matriz_final), " genes x ",
        ncol(prostata_expresion_matriz_final), " samples")


# 13. DESCRIPTOR DEL DATASET ___________________________________________________

prostata_descriptor_lista <-
  list(
    gse_id                  = gse_objetivo_id_char,
    organo                  = organo_objetivo_char,
    plataforma_tipo         = plataforma_tipo_char,
    gene_id_convention      = "HGNC symbol (with Ensembl ID preserved)",
    matrix_value_nature     = "log2(FPKM + 1) transformed expression",
    dge_branch_recommended  = "limma on log2(FPKM + 1) with paciente_id blocking",
    n_total                 = ncol(prostata_expresion_matriz_final),
    n_tumor                 = sum(prostata_metadatos_tbl_final$grupo_diagnostico
                                  == "Tumor", na.rm = TRUE),
    n_normal                = sum(prostata_metadatos_tbl_final$grupo_diagnostico
                                  == "Normal", na.rm = TRUE),
    n_sin_etiqueta          = sum(is.na(prostata_metadatos_tbl_final$grupo_diagnostico)),
    n_patients              = length(unique(prostata_metadatos_tbl_final$paciente_id)),
    comparator_subtype      = "Adjacent normal tissue (patient-matched)",
    paired_design           = TRUE,
    ethnicity_balance       = "16 African American + 16 European American",
    alignment_strategy      = estrategia_alineacion_char,
    typos_corrected         = "Patient 881: '881_819' -> '881_B19', '881_826' -> '881_B26'",
    feature_filtering       = "Retained features with non-empty HGNC symbol",
    duplicate_collapse      = "Maximum mean expression per gene symbol",
    fecha_estandarizacion   = format(Sys.Date(), "%Y-%m-%d")
  )


# 14. EXPORTAR EL RDS ESTANDARIZADO ____________________________________________

prostata_dataset_estandarizado_lista <-
  list(
    conteos_matriz         = prostata_expresion_matriz_final,
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


# 15. RESUMEN __________________________________________________________________

message("")
message("Pipeline completed for ", gse_objetivo_id_char)
message("Organ: ",
        prostata_descriptor_lista$organo)
message("Platform: ",
        prostata_descriptor_lista$plataforma_tipo)
message("Matrix nature: ",
        prostata_descriptor_lista$matrix_value_nature)
message("DGE branch recommended: ",
        prostata_descriptor_lista$dge_branch_recommended)
message("Comparator subtype: ",
        prostata_descriptor_lista$comparator_subtype)
message("Paired design: ",
        prostata_descriptor_lista$paired_design,
        " (", prostata_descriptor_lista$n_patients, " patients)")
message("Ethnicity balance: ",
        prostata_descriptor_lista$ethnicity_balance)
message("Alignment strategy: ",
        prostata_descriptor_lista$alignment_strategy)
message("Transcription typos corrected: ",
        prostata_descriptor_lista$typos_corrected)
message("Features retained: ",
        nrow(prostata_dataset_estandarizado_lista$conteos_matriz),
        " unique gene symbols")
message("Samples retained: ",
        prostata_descriptor_lista$n_total,
        " (Tumor = ", prostata_descriptor_lista$n_tumor,
        ", Normal = ", prostata_descriptor_lista$n_normal, ")")
message("Gene identifier convention: ",
        prostata_descriptor_lista$gene_id_convention)
message("Output: ", prostata_rds_destino_path)
