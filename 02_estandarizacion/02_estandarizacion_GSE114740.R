# ______________________________________________________________________________
# 02_02  Estandarizacion de GSE114740 (prostata, pareado tumor/normal)
#
# Que hace: pasa la matriz RPKM y los metadatos de GSE114740 al esquema
#           estandar del proyecto (simbolo HGNC) y lo guarda.
# Entra:    Prostata/GSE114740*metadatos*.csv, Prostata_Counts/GSE114740*.xlsx.gz
# Sale:     Datasets_Estandarizados/GSE114740_estandarizado.rds
# ______________________________________________________________________________

# __ INDEX ___
# 01. LIBRERIAS
# 02. RUTAS
# 03. LOCALIZAR ARCHIVOS
# 04. IMPORTAR METADATOS CLINICOS
# 05. IMPORTAR MATRIZ DE EXPRESION (RPKM)
# 06. ALINEAR MUESTRAS ENTRE METADATOS Y MATRIZ
# 07. ASIGNAR GRUPO DIAGNOSTICO E IDENTIFICAR PACIENTES PAREADOS
# 08. CONTROL DE CALIDAD DE LOS GENES
# 09. CONSTRUIR ANOTACION DE GENES
# 10. ENSAMBLAR LA MATRIZ ESTANDARIZADA
# 11. DESCRIPTOR DEL DATASET
# 12. EXPORTAR EL RDS ESTANDARIZADO
# 13. RESUMEN
# ________________________________________________________________________


# 01. LIBRERIAS ________________________________________________________________

library(tidyverse)
library(data.table)
library(readxl)
library(R.utils)


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

gse_objetivo_id_char  <- "GSE114740"
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

prostata_xlsx_gz_archivo_path <-
  list.files(
    path       = prostata_conteos_directorio_path,
    pattern    = paste0(gse_objetivo_id_char, ".*\\.xlsx\\.gz$"),
    recursive  = TRUE,
    full.names = TRUE
  )

if (length(prostata_metadatos_archivo_path) == 0 ||
    length(prostata_xlsx_gz_archivo_path) == 0) {
  message("Metadata files located: ",
          length(prostata_metadatos_archivo_path))
  message("Compressed Excel files located: ",
          length(prostata_xlsx_gz_archivo_path))
  stop("Critical error: incomplete local repository for ",
       gse_objetivo_id_char)
}

message("Metadata file located: ",
        basename(prostata_metadatos_archivo_path[1]))
message("Compressed Excel located: ",
        basename(prostata_xlsx_gz_archivo_path[1]))


# 04. IMPORTAR METADATOS CLINICOS ______________________________________________

prostata_metadatos_tbl_inicial <-
  prostata_metadatos_archivo_path[1] %>%
  fread() %>%
  as_tibble()

message("Clinical metadata imported: ",
        nrow(prostata_metadatos_tbl_inicial), " samples, ",
        ncol(prostata_metadatos_tbl_inicial), " annotation columns")

# El campo 'isolate' de characteristics_ch1.2 lleva a la vez el
# identificador del paciente pareado y si la muestra es tumor o normal.

if ("characteristics_ch1.2" %in% names(prostata_metadatos_tbl_inicial)) {
  message("Native isolate distribution (characteristics_ch1.2):")
  print(head(table(prostata_metadatos_tbl_inicial$characteristics_ch1.2,
                   useNA = "ifany"), 20))
}


# 05. IMPORTAR MATRIZ DE EXPRESION (RPKM) ______________________________________

# La matriz de expresion viene como un Excel de una sola hoja comprimido
# en gzip. Lo descomprimimos a un archivo temporal antes de leerlo con
# readxl, porque R no puede leer un .xlsx directamente desde un gzip.

prostata_xlsx_temp_path <- tempfile(fileext = ".xlsx")

R.utils::gunzip(
  filename  = prostata_xlsx_gz_archivo_path[1],
  destname  = prostata_xlsx_temp_path,
  remove    = FALSE,
  overwrite = TRUE
)

prostata_excel_hojas_vector <- excel_sheets(prostata_xlsx_temp_path)

message("Excel sheets detected: ",
        paste(prostata_excel_hojas_vector, collapse = ", "))

prostata_rpkm_tbl_inicial <-
  read_excel(
    path  = prostata_xlsx_temp_path,
    sheet = 1
  )

message("Expression matrix imported: ",
        nrow(prostata_rpkm_tbl_inicial), " features x ",
        ncol(prostata_rpkm_tbl_inicial), " columns")

# Separamos las columnas de anotacion (Gene ID y GeneSymbol) de las
# columnas de muestras. Los nombres de las muestras son descriptivos y
# reflejan el diseno pareado (por ejemplo 'PC_1 prostate tumor').

prostata_columnas_anotacion_vector <- c("Gene ID", "GeneSymbol")

prostata_columnas_muestras_vector <-
  setdiff(names(prostata_rpkm_tbl_inicial),
          prostata_columnas_anotacion_vector)

message("Annotation columns: ",
        paste(prostata_columnas_anotacion_vector, collapse = ", "))
message("Biological sample columns: ",
        length(prostata_columnas_muestras_vector))


# 06. ALINEAR MUESTRAS ENTRE METADATOS Y MATRIZ ________________________________

# La matriz etiqueta las muestras con nombres descriptivos (por ejemplo
# 'PC_1 prostate tumor', 'PC_1 matched normal tissue') en vez de GSM.
# Esas etiquetas coinciden con el campo 'title' de los metadatos de GEO,
# que es la clave de cruce para este dataset.

prostata_match_titulo_vector <-
  intersect(
    prostata_columnas_muestras_vector,
    prostata_metadatos_tbl_inicial$title
  )

message("Title-based alignment: ",
        length(prostata_match_titulo_vector), " / ",
        length(prostata_columnas_muestras_vector))

if (length(prostata_match_titulo_vector) ==
    length(prostata_columnas_muestras_vector)) {
  estrategia_alineacion_char <- "GEO title field"
  prostata_metadatos_tbl_armonizado <-
    prostata_metadatos_tbl_inicial %>%
    mutate(id_alineacion = title)
} else {
  estrategia_alineacion_char <- "Partial; manual review required"
  prostata_metadatos_tbl_armonizado <-
    prostata_metadatos_tbl_inicial %>%
    mutate(id_alineacion = title)
  message("Warning: alignment incomplete; manual review required")
}

message("Alignment strategy selected: ", estrategia_alineacion_char)


# 07. ASIGNAR GRUPO DIAGNOSTICO E IDENTIFICAR PACIENTES PAREADOS ________________

# Cada title lleva el grupo diagnostico y el paciente pareado en la forma
# '<codigo_paciente> <descripcion_tejido>', donde el codigo sigue el
# patron 'PC_<digitos>'. El tejido normal es normal adyacente real del
# mismo paciente (no BPH), lo que hace de este el comparador mas limpio
# del grupo prostatico.

extraer_paciente_id_func <- function(titulo_text) {
  patron_codigo <- regmatches(
    titulo_text,
    regexpr("PC_\\d+", titulo_text)
  )
  ifelse(length(patron_codigo) == 0, NA_character_, patron_codigo)
}

prostata_metadatos_tbl_curado <-
  prostata_metadatos_tbl_armonizado %>%
  rowwise() %>%
  mutate(
    paciente_id = extraer_paciente_id_func(title)
  ) %>%
  ungroup() %>%
  mutate(
    grupo_diagnostico = case_when(
      str_detect(title,
                 regex("tumor|cancer|carcinoma|adenocarcinoma",
                       ignore_case = TRUE))
      ~ "Tumor",
      str_detect(title,
                 regex("normal|adjacent|matched",
                       ignore_case = TRUE))
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
    edad_paciente = str_remove(characteristics_ch1.1,
                               regex("^age:\\s*",
                                     ignore_case = TRUE)) %>%
      str_remove(regex("\\s*years?\\s*old\\s*$",
                       ignore_case = TRUE))
  )

message("Standardized phenotypic distribution:")
print(table(prostata_metadatos_tbl_curado$grupo_diagnostico,
            useNA = "ifany"))

message("Paired patient structure:")
print(prostata_metadatos_tbl_curado %>%
        count(paciente_id, grupo_diagnostico) %>%
        pivot_wider(names_from = grupo_diagnostico,
                    values_from = n,
                    values_fill = 0))

# Construimos la tabla de alineamiento que respeta el orden de columnas
# de la matriz de expresion.

prostata_grupos_alineados_tbl <-
  tibble(sample_id_matriz = prostata_columnas_muestras_vector) %>%
  mutate(
    grupo_diagnostico = prostata_metadatos_tbl_curado$grupo_diagnostico[
      match(sample_id_matriz, prostata_metadatos_tbl_curado$id_alineacion)
    ],
    tejido_control_subtipo = prostata_metadatos_tbl_curado$tejido_control_subtipo[
      match(sample_id_matriz, prostata_metadatos_tbl_curado$id_alineacion)
    ],
    paciente_id = prostata_metadatos_tbl_curado$paciente_id[
      match(sample_id_matriz, prostata_metadatos_tbl_curado$id_alineacion)
    ],
    edad_paciente = prostata_metadatos_tbl_curado$edad_paciente[
      match(sample_id_matriz, prostata_metadatos_tbl_curado$id_alineacion)
    ],
    gsm_accession = prostata_metadatos_tbl_curado$geo_accession[
      match(sample_id_matriz, prostata_metadatos_tbl_curado$id_alineacion)
    ]
  )

message("Final diagnostic balance after alignment:")
print(table(prostata_grupos_alineados_tbl$grupo_diagnostico,
            useNA = "ifany"))


# 08. CONTROL DE CALIDAD DE LOS GENES __________________________________________

# Anotamos que identificadores trae GSE114740: da a la vez Entrez Gene
# IDs (numericos) y simbolos HGNC, con lo que hay doble cobertura de
# anotacion frente a los datasets con un solo identificador.

prostata_ids_entrez_vector <-
  prostata_rpkm_tbl_inicial[["Gene ID"]]

prostata_ids_symbol_vector <-
  prostata_rpkm_tbl_inicial[["GeneSymbol"]]

message("Entrez ID coverage: ",
        sum(!is.na(prostata_ids_entrez_vector)), " / ",
        length(prostata_ids_entrez_vector))
message("Gene symbol coverage: ",
        sum(!is.na(prostata_ids_symbol_vector) &
              prostata_ids_symbol_vector != ""),
        " / ", length(prostata_ids_symbol_vector))

# Miramos el rango de valores para confirmar que la matriz es RPKM.

prostata_valores_muestra_matriz <-
  prostata_rpkm_tbl_inicial %>%
  select(all_of(prostata_columnas_muestras_vector[
    1:min(5, length(prostata_columnas_muestras_vector))
  ])) %>%
  as.matrix()

message("Value range across first 5 samples:")
message("  min = ",
        round(min(prostata_valores_muestra_matriz, na.rm = TRUE), 3))
message("  max = ",
        round(max(prostata_valores_muestra_matriz, na.rm = TRUE), 3))
message("  median = ",
        round(median(prostata_valores_muestra_matriz, na.rm = TRUE), 3))

prostata_valores_son_enteros_lgl <-
  all(prostata_valores_muestra_matriz ==
        round(prostata_valores_muestra_matriz), na.rm = TRUE)

prostata_valores_tienen_negativos_lgl <-
  any(prostata_valores_muestra_matriz < 0, na.rm = TRUE)

message("Integer-valued matrix (would suggest raw counts): ",
        prostata_valores_son_enteros_lgl)
message("Negative values present: ",
        prostata_valores_tienen_negativos_lgl)

# Comprobacion: marcadores clasicos de biologia prostatica.

prostata_genes_marcadores_vector <-
  c("AR", "KLK3", "FOXA1", "TMPRSS2", "NKX3-1", "ESR1", "TP53")

prostata_marcadores_presentes_vector <-
  intersect(prostata_genes_marcadores_vector,
            prostata_ids_symbol_vector)

message("Canonical prostate markers detected (by symbol): ",
        paste(prostata_marcadores_presentes_vector, collapse = ", "))


# 09. CONSTRUIR ANOTACION DE GENES _____________________________________________

# Conservamos tanto la anotacion Entrez como la de simbolo HGNC. El
# contenedor estandarizado usa el simbolo HGNC como identificador de fila
# (cuando existe) porque el meta-analisis trabaja sobre simbolos.

prostata_genes_anotacion_tbl_inicial <-
  tibble(
    gene_id_native = as.character(prostata_ids_entrez_vector),
    entrez_id      = as.character(prostata_ids_entrez_vector),
    gene_symbol    = prostata_ids_symbol_vector,
    ensembl_id     = NA_character_
  )

# Quitamos los genes sin simbolo HGNC valido, porque no se pueden cruzar
# a nivel de gen con el resto de datasets. Luego colapsamos duplicados
# por maxima expresion media, la regla estandar para integrar datos.

prostata_rpkm_tbl_02 <-
  prostata_rpkm_tbl_inicial %>%
  rename(entrez_id = `Gene ID`,
         gene_symbol = GeneSymbol) %>%
  mutate(entrez_id = as.character(entrez_id)) %>%
  filter(!is.na(gene_symbol) & gene_symbol != "")

message("Features with valid HGNC symbol: ",
        nrow(prostata_rpkm_tbl_02), " / ",
        nrow(prostata_rpkm_tbl_inicial))

prostata_rpkm_tbl_curado <-
  prostata_rpkm_tbl_02 %>%
  mutate(
    expresion_media = rowMeans(
      select(., -entrez_id, -gene_symbol),
      na.rm = TRUE
    )
  ) %>%
  arrange(gene_symbol, desc(expresion_media)) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  select(-expresion_media)

message("Genes retained after duplicate collapse: ",
        nrow(prostata_rpkm_tbl_curado))


# 10. ENSAMBLAR LA MATRIZ ESTANDARIZADA ________________________________________

# Armamos el tibble de anotacion final, la matriz numerica indexada por
# simbolo HGNC y los metadatos por paciente respetando el orden de
# columnas de la matriz.

prostata_genes_anotacion_tbl_final <-
  tibble(
    gene_id_native = prostata_rpkm_tbl_curado$entrez_id,
    entrez_id      = prostata_rpkm_tbl_curado$entrez_id,
    gene_symbol    = prostata_rpkm_tbl_curado$gene_symbol,
    ensembl_id     = NA_character_
  )

prostata_expresion_matriz_final <-
  prostata_rpkm_tbl_curado %>%
  select(-entrez_id, -gene_symbol) %>%
  as.matrix()

rownames(prostata_expresion_matriz_final) <-
  prostata_genes_anotacion_tbl_final$gene_symbol

prostata_metadatos_tbl_final <-
  prostata_grupos_alineados_tbl %>%
  rename(sample_id = sample_id_matriz) %>%
  select(sample_id, gsm_accession, paciente_id,
         grupo_diagnostico, tejido_control_subtipo,
         edad_paciente) %>%
  slice(match(colnames(prostata_expresion_matriz_final), sample_id))

message("Standardized matrix assembled: ",
        nrow(prostata_expresion_matriz_final), " genes x ",
        ncol(prostata_expresion_matriz_final), " samples")


# 11. DESCRIPTOR DEL DATASET ___________________________________________________

# Descriptor autoexplicativo que marca que este dataset va a la rama de
# limma del analisis de expresion diferencial y no a DESeq2, por ser RPKM.

prostata_descriptor_lista <-
  list(
    gse_id                  = gse_objetivo_id_char,
    organo                  = organo_objetivo_char,
    plataforma_tipo         = plataforma_tipo_char,
    gene_id_convention      = "HGNC symbol (with Entrez ID preserved)",
    matrix_value_nature     = "RPKM (pre-normalized expression units)",
    dge_branch_recommended  = "limma on log2(RPKM + 1) with paciente_id blocking",
    n_total                 = ncol(prostata_expresion_matriz_final),
    n_tumor                 = sum(prostata_metadatos_tbl_final$grupo_diagnostico
                                  == "Tumor", na.rm = TRUE),
    n_normal                = sum(prostata_metadatos_tbl_final$grupo_diagnostico
                                  == "Normal", na.rm = TRUE),
    n_sin_etiqueta          = sum(is.na(prostata_metadatos_tbl_final$grupo_diagnostico)),
    comparator_subtype      = "Adjacent normal tissue (patient-matched)",
    paired_design           = TRUE,
    n_patients              = length(unique(prostata_metadatos_tbl_final$paciente_id)),
    alignment_strategy      = estrategia_alineacion_char,
    feature_filtering       = "Retained features with non-empty HGNC symbol",
    duplicate_collapse      = "Maximum mean expression per gene symbol",
    fecha_estandarizacion   = format(Sys.Date(), "%Y-%m-%d")
  )


# 12. EXPORTAR EL RDS ESTANDARIZADO ____________________________________________

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


# 13. RESUMEN __________________________________________________________________

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
message("Alignment strategy: ",
        prostata_descriptor_lista$alignment_strategy)
message("Features retained: ",
        nrow(prostata_dataset_estandarizado_lista$conteos_matriz),
        " unique gene symbols")
message("Samples retained: ",
        prostata_descriptor_lista$n_total,
        " (Tumor = ", prostata_descriptor_lista$n_tumor,
        ", Normal = ", prostata_descriptor_lista$n_normal,
        if (prostata_descriptor_lista$n_sin_etiqueta > 0)
          paste0(", Unlabeled = ",
                 prostata_descriptor_lista$n_sin_etiqueta)
        else "",
        ")")
message("Gene identifier convention: ",
        prostata_descriptor_lista$gene_id_convention)
message("Feature filtering: ",
        prostata_descriptor_lista$feature_filtering)
message("Duplicate collapse rule: ",
        prostata_descriptor_lista$duplicate_collapse)
message("Output: ", prostata_rds_destino_path)
