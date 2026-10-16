# ______________________________________________________________________________
# 02_04  Estandarizacion de GSE214940 (prostata, estroma pareado tumor/normal)
#
# Que hace: pasa la matriz de conteos crudos de Salmon y los metadatos de
#           GSE214940 al esquema estandar del proyecto (IDs Ensembl) y lo guarda.
# Entra:    Prostata/GSE214940*metadatos*.csv,
#           Prostata_Counts/GSE214940*raw_salmon*.csv.gz
# Sale:     Datasets_Estandarizados/GSE214940_estandarizado.rds
# ______________________________________________________________________________

# __ INDEX ___
# 01. LIBRERIAS
# 02. RUTAS
# 03. LOCALIZAR ARCHIVOS
# 04. IMPORTAR METADATOS CLINICOS
# 05. IMPORTAR MATRIZ DE CONTEOS SALMON CRUDOS
# 06. CLASIFICAR COMPARTIMENTO ESTROMAL E IDENTIFICAR PACIENTES
# 07. RESTRINGIR AL SUBCONJUNTO DE PACIENTES PAREADOS
# 08. CONSERVAR VARIABLES CLINICOPATOLOGICAS
# 09. ALINEAR MUESTRAS Y REDONDEAR A ENTEROS
# 10. CONTROL DE CALIDAD DE LOS GENES
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

gse_objetivo_id_char  <- "GSE214940"
organo_objetivo_char  <- "Prostate"
plataforma_tipo_char  <- "RNA-seq (Salmon, gene-level)"


# 03. LOCALIZAR ARCHIVOS _______________________________________________________

message("Auditing dataset ", gse_objetivo_id_char)

prostata_metadatos_archivo_path <-
  list.files(
    path       = prostata_metadatos_directorio_path,
    pattern    = paste0(gse_objetivo_id_char, ".*metadatos.*\\.csv$"),
    recursive  = TRUE,
    full.names = TRUE
  )

prostata_salmon_archivo_path <-
  list.files(
    path       = prostata_conteos_directorio_path,
    pattern    = paste0(gse_objetivo_id_char,
                        ".*raw_salmon.*\\.csv\\.gz$"),
    recursive  = TRUE,
    full.names = TRUE
  )

if (length(prostata_metadatos_archivo_path) == 0 ||
    length(prostata_salmon_archivo_path) == 0) {
  message("Metadata files located: ",
          length(prostata_metadatos_archivo_path))
  message("Raw Salmon files located: ",
          length(prostata_salmon_archivo_path))
  stop("Critical error: incomplete local repository for ",
       gse_objetivo_id_char)
}

message("Metadata file located: ",
        basename(prostata_metadatos_archivo_path[1]))
message("Raw Salmon matrix located: ",
        basename(prostata_salmon_archivo_path[1]))


# 04. IMPORTAR METADATOS CLINICOS ______________________________________________

prostata_metadatos_tbl_inicial <-
  prostata_metadatos_archivo_path[1] %>%
  fread() %>%
  as_tibble()

message("Clinical metadata imported: ",
        nrow(prostata_metadatos_tbl_inicial), " samples, ",
        ncol(prostata_metadatos_tbl_inicial), " annotation columns")


# 05. IMPORTAR MATRIZ DE CONTEOS SALMON CRUDOS _________________________________

# Los conteos de Salmon a nivel de gen son estimaciones con decimales
# (no enteros estrictos). Los redondeamos al entero mas cercano antes de
# entrar en la rama de DESeq2 del analisis de expresion diferencial.

prostata_salmon_tbl_inicial <-
  prostata_salmon_archivo_path[1] %>%
  fread() %>%
  as_tibble()

names(prostata_salmon_tbl_inicial)[1] <- "ensembl_id"

message("Raw Salmon matrix imported: ",
        nrow(prostata_salmon_tbl_inicial), " features x ",
        ncol(prostata_salmon_tbl_inicial) - 1, " sequencing libraries")

prostata_libraries_columnas_vector <-
  setdiff(names(prostata_salmon_tbl_inicial), "ensembl_id")


# 06. CLASIFICAR COMPARTIMENTO ESTROMAL E IDENTIFICAR PACIENTES _________________

# Los titulos de cada libreria llevan el codigo del paciente (los digitos
# iniciales) y el compartimento estromal ('SN' = Stroma Normal, 'ST' =
# Stroma Tumor). Esos dos campos definen juntos el diseno experimental.

prostata_metadatos_tbl_curado <-
  prostata_metadatos_tbl_inicial %>%
  mutate(
    paciente_id = str_extract(title, "^\\d+"),
    compartimento_estromal = case_when(
      str_detect(title, "_SN_") ~ "Stroma Normal",
      str_detect(title, "_ST_") ~ "Stroma Tumor",
      TRUE                      ~ NA_character_
    ),
    grupo_diagnostico = case_when(
      compartimento_estromal == "Stroma Tumor"  ~ "Tumor",
      compartimento_estromal == "Stroma Normal" ~ "Normal",
      TRUE                                       ~ NA_character_
    ),
    tejido_control_subtipo = case_when(
      grupo_diagnostico == "Normal" ~ "Stroma normal (intra-patient)",
      TRUE                          ~ NA_character_
    )
  )

message("Stromal compartment distribution:")
print(table(prostata_metadatos_tbl_curado$compartimento_estromal,
            useNA = "ifany"))

prostata_n_pacientes_unicos_int <-
  length(unique(prostata_metadatos_tbl_curado$paciente_id))

message("Unique patient codes: ", prostata_n_pacientes_unicos_int)


# 07. RESTRINGIR AL SUBCONJUNTO DE PACIENTES PAREADOS __________________________

# Buscamos los pacientes que tienen los dos compartimentos (SN y ST).
# Solo esos entran en el analisis, con lo que queda un diseno pareado
# limpio. Si un paciente tiene varias muestras del mismo compartimento
# las conservamos; el modelo absorbe la correlacion dentro del paciente
# mediante el factor de bloqueo por paciente.

prostata_pareamiento_tbl <-
  prostata_metadatos_tbl_curado %>%
  count(paciente_id, compartimento_estromal) %>%
  pivot_wider(
    names_from  = compartimento_estromal,
    values_from = n,
    values_fill = 0
  )

prostata_pacientes_pareados_vector <-
  prostata_pareamiento_tbl %>%
  filter(`Stroma Normal` > 0 & `Stroma Tumor` > 0) %>%
  pull(paciente_id)

message("Patients with paired SN and ST specimens: ",
        length(prostata_pacientes_pareados_vector))

prostata_metadatos_tbl_pareado <-
  prostata_metadatos_tbl_curado %>%
  filter(paciente_id %in% prostata_pacientes_pareados_vector)

message("Libraries retained after pairing restriction: ",
        nrow(prostata_metadatos_tbl_pareado),
        " (from ", nrow(prostata_metadatos_tbl_curado), " original)")

message("Final group balance:")
print(table(prostata_metadatos_tbl_pareado$grupo_diagnostico))


# 08. CONSERVAR VARIABLES CLINICOPATOLOGICAS ___________________________________

# Extraemos Gleason, letalidad y morfologia cribiforme, quitando los
# prefijos de texto para dejar las variables clinicas limpias.

prostata_metadatos_tbl_clinico <-
  prostata_metadatos_tbl_pareado %>%
  mutate(
    gleason_score = str_remove(characteristics_ch1,
                               regex("^gleason score:\\s*",
                                     ignore_case = TRUE)),
    lethal_status = str_remove(characteristics_ch1.1,
                               regex("^lethal:\\s*",
                                     ignore_case = TRUE)),
    cribriform_morphology = str_remove(characteristics_ch1.2,
                                       regex("^cribriform:\\s*",
                                             ignore_case = TRUE)),
    used_in_original_analysis = str_remove(characteristics_ch1.3,
                                           regex("^used in analysis:\\s*",
                                                 ignore_case = TRUE))
  )


# 09. ALINEAR MUESTRAS Y REDONDEAR A ENTEROS ___________________________________

# Esperamos que las columnas de la matriz crucen directamente con el
# campo 'title' de los metadatos, porque ambos llevan la cadena paciente
# + compartimento (por ejemplo '11_SN_'). Restringimos las columnas de la
# matriz a las librerias que quedan en el subconjunto pareado.

prostata_libraries_retenidas_vector <-
  intersect(prostata_libraries_columnas_vector,
            prostata_metadatos_tbl_clinico$title)

message("Matrix columns retained after pairing restriction: ",
        length(prostata_libraries_retenidas_vector), " / ",
        length(prostata_libraries_columnas_vector))

prostata_salmon_tbl_pareado <-
  prostata_salmon_tbl_inicial %>%
  select(ensembl_id, all_of(prostata_libraries_retenidas_vector))

# Redondeamos las estimaciones de Salmon al entero mas cercano para que
# DESeq2 las acepte. Es lo que hace tximport por defecto cuando se usa
# 'countsFromAbundance'; mantiene la estructura relativa de conteos que
# necesita el modelo binomial negativo.

prostata_salmon_matriz_redondeada <-
  prostata_salmon_tbl_pareado %>%
  select(-ensembl_id) %>%
  as.matrix() %>%
  round() %>%
  apply(2, as.integer)

rownames(prostata_salmon_matriz_redondeada) <-
  prostata_salmon_tbl_pareado$ensembl_id

message("Salmon estimates rounded to nearest integer for DESeq2 compatibility")


# 10. CONTROL DE CALIDAD DE LOS GENES __________________________________________

# La matriz ya viene prefiltrada por los autores originales. Confirmamos
# que los IDs son Ensembl canonicos y que tras redondear son enteros.

prostata_ids_genes_vector <- rownames(prostata_salmon_matriz_redondeada)

prostata_ids_son_ensembl_lgl <-
  mean(grepl("^ENSG\\d+$",
             head(prostata_ids_genes_vector, 100))) > 0.9

message("Gene identifier convention: ",
        if (prostata_ids_son_ensembl_lgl) "Canonical Ensembl (no version)"
        else "Mixed / non-canonical")

prostata_valores_muestra_matriz <-
  prostata_salmon_matriz_redondeada[, 1:min(5, ncol(prostata_salmon_matriz_redondeada))]

message("Rounded count range across first 5 libraries:")
message("  min = ", min(prostata_valores_muestra_matriz, na.rm = TRUE))
message("  max = ", max(prostata_valores_muestra_matriz, na.rm = TRUE))
message("  median = ", round(median(prostata_valores_muestra_matriz, na.rm = TRUE), 2))

prostata_valores_son_enteros_lgl <-
  all(prostata_valores_muestra_matriz ==
        round(prostata_valores_muestra_matriz), na.rm = TRUE)

message("Integer-valued matrix (raw counts): ",
        prostata_valores_son_enteros_lgl)


# 11. ENSAMBLAR LA MATRIZ ESTANDARIZADA ________________________________________

prostata_genes_anotacion_tbl_final <-
  tibble(
    gene_id_native = rownames(prostata_salmon_matriz_redondeada),
    ensembl_id     = rownames(prostata_salmon_matriz_redondeada),
    gene_symbol    = NA_character_
  )

# Reordenamos los metadatos para que sigan el orden de columnas de la matriz.

prostata_metadatos_tbl_final <-
  prostata_metadatos_tbl_clinico %>%
  select(sample_id = title,
         gsm_accession = geo_accession,
         paciente_id,
         compartimento_estromal,
         grupo_diagnostico,
         tejido_control_subtipo,
         gleason_score,
         lethal_status,
         cribriform_morphology,
         used_in_original_analysis) %>%
  slice(match(colnames(prostata_salmon_matriz_redondeada), sample_id))

message("Standardized matrix assembled: ",
        nrow(prostata_salmon_matriz_redondeada), " genes x ",
        ncol(prostata_salmon_matriz_redondeada), " libraries")


# 12. DESCRIPTOR DEL DATASET ___________________________________________________

prostata_descriptor_lista <-
  list(
    gse_id                   = gse_objetivo_id_char,
    organo                   = organo_objetivo_char,
    plataforma_tipo          = plataforma_tipo_char,
    gene_id_convention       = "Canonical Ensembl (no version suffix)",
    matrix_value_nature      = "Integer counts (Salmon estimates rounded)",
    dge_branch_recommended   = "DESeq2 with ~ paciente_id + grupo_diagnostico",
    n_total                  = ncol(prostata_salmon_matriz_redondeada),
    n_tumor                  = sum(prostata_metadatos_tbl_final$grupo_diagnostico
                                   == "Tumor", na.rm = TRUE),
    n_normal                 = sum(prostata_metadatos_tbl_final$grupo_diagnostico
                                   == "Normal", na.rm = TRUE),
    n_sin_etiqueta           = sum(is.na(prostata_metadatos_tbl_final$grupo_diagnostico)),
    n_patients_paired        = length(prostata_pacientes_pareados_vector),
    biological_compartment   = "Prostate stroma (not glandular epithelium)",
    comparator_subtype       = "Intra-patient paired stroma normal",
    paired_design            = TRUE,
    n_libraries_original     = length(prostata_libraries_columnas_vector),
    n_patients_original      = prostata_n_pacientes_unicos_int,
    pairing_restriction      = "Restricted to patients with both SN and ST specimens",
    salmon_rounding_applied  = TRUE,
    pre_filtered_by_authors  = TRUE,
    n_features_pre_filtered  = nrow(prostata_salmon_matriz_redondeada),
    fecha_estandarizacion    = format(Sys.Date(), "%Y-%m-%d")
  )


# 13. EXPORTAR EL RDS ESTANDARIZADO ____________________________________________

prostata_dataset_estandarizado_lista <-
  list(
    conteos_matriz         = prostata_salmon_matriz_redondeada,
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
message("Organ: ",
        prostata_descriptor_lista$organo)
message("Platform: ",
        prostata_descriptor_lista$plataforma_tipo)
message("Biological compartment: ",
        prostata_descriptor_lista$biological_compartment)
message("Matrix nature: ",
        prostata_descriptor_lista$matrix_value_nature)
message("DGE branch recommended: ",
        prostata_descriptor_lista$dge_branch_recommended)
message("Paired design: ",
        prostata_descriptor_lista$paired_design,
        " (", prostata_descriptor_lista$n_patients_paired,
        " patients with both compartments)")
message("Original cohort: ",
        prostata_descriptor_lista$n_libraries_original,
        " libraries from ",
        prostata_descriptor_lista$n_patients_original, " patients")
message("Restricted to paired subset: ",
        prostata_descriptor_lista$n_total,
        " libraries (Tumor = ",
        prostata_descriptor_lista$n_tumor,
        ", Normal = ", prostata_descriptor_lista$n_normal,
        if (prostata_descriptor_lista$n_sin_etiqueta > 0)
          paste0(", Unlabeled = ",
                 prostata_descriptor_lista$n_sin_etiqueta)
        else "",
        ")")
message("Features retained: ",
        nrow(prostata_dataset_estandarizado_lista$conteos_matriz),
        " Ensembl genes (pre-filtered by authors)")
message("Salmon estimates rounded to integers: ",
        prostata_descriptor_lista$salmon_rounding_applied)
message("Output: ", prostata_rds_destino_path)
