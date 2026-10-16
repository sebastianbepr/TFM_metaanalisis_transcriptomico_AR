# ______________________________________________________________________________
# 02_11  Estandarizacion de TCGA_GTEX_Testis (testiculo)
#
# Que hace: pasa la cohorte armonizada TCGA-TGCT + GTEx (recomputo Toil,
#           testiculo) al esquema estandar del proyecto (Ensembl sin version,
#           escala log2(conteo + 1)) y la guarda.
# Entra:    TCGA_GTEX/grupos_testiculo.csv, gencode.v23.gene.probemap, TcgaTargetGtex_gene_expected_count.gz
# Sale:     Datasets_Estandarizados/TCGA_GTEX_Testis_estandarizado.rds
# ______________________________________________________________________________

# __ INDEX ___
# 01. LIBRERIAS
# 02. RUTAS
# 03. CARGAR GRUPOS DE MUESTRAS Y PROBEMAP
# 04. CARGAR LA MATRIZ TOIL CRUDA (SOLO TESTICULO)
# 05. ASIGNAR GRUPO DIAGNOSTICO Y ORIGEN
# 06. ANOTACION ENSEMBL A SIMBOLO
# 07. CONTROL DE CALIDAD DE LOS GENES
# 08. ENSAMBLAR LA MATRIZ ESTANDARIZADA
# 09. DESCRIPTOR DEL DATASET
# 10. EXPORTAR EL RDS ESTANDARIZADO
# 11. RESUMEN
# ________________________________________________________________________


# 01. LIBRERIAS ________________________________________________________________

library(tidyverse)
library(data.table)


# 02. RUTAS ____________________________________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/TCGA_GTEX"
setwd(ruta)
tcga_fuente_directorio_path <- ruta   # el resto del script usa este nombre sin cambios

# Carpeta del repositorio general donde se guardan los datasets estandarizados.
base_datos_directorio_path <-
  "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"

datasets_estandarizados_directorio_path <-
  file.path(base_datos_directorio_path, "Datasets_Estandarizados")

if (!dir.exists(datasets_estandarizados_directorio_path)) {
  dir.create(datasets_estandarizados_directorio_path, recursive = TRUE)
}

dataset_id_char         <- "TCGA_GTEX_Testis"
organo_objetivo_char    <- "Testis"
plataforma_tipo_char    <- "RNA-seq (Toil recompute, harmonized TCGA + GTEx)"


# 03. CARGAR GRUPOS DE MUESTRAS Y PROBEMAP _____________________________________

message("Standardizing dataset ", dataset_id_char)

# Grupos de muestras: 411 filas con sample, grupo y origen.
testiculo_grupos_archivo_path <-
  file.path(tcga_fuente_directorio_path, "grupos_testiculo.csv")

testiculo_grupos_tbl_inicial <-
  testiculo_grupos_archivo_path %>%
  fread() %>%
  as_tibble()

message("Sample groups imported: ",
        nrow(testiculo_grupos_tbl_inicial), " samples")
message("Native distribution by origin:")
print(table(testiculo_grupos_tbl_inicial$origen,
            testiculo_grupos_tbl_inicial$grupo))

# Probemap: gene_id -> simbolo HGNC (Gencode v23).
testiculo_probemap_archivo_path <-
  file.path(tcga_fuente_directorio_path, "gencode.v23.gene.probemap")

testiculo_probemap_tbl <-
  testiculo_probemap_archivo_path %>%
  fread() %>%
  as_tibble() %>%
  select(gene_id_native = id,
         gene_symbol    = gene)

message("Probemap imported: ",
        nrow(testiculo_probemap_tbl), " gene-to-symbol mappings")


# 04. CARGAR LA MATRIZ TOIL CRUDA (SOLO TESTICULO) _____________________________

# Para no cargar de memoria un equipo de 16 GB, leemos la matriz cruda
# seleccionando solo las columnas de las 411 muestras testiculares mas la
# columna del identificador de gen. Asi la matriz de trabajo se queda en unos
# 60.498 x 412, bien por debajo de 300 MB.

testiculo_matriz_cruda_archivo_path <-
  file.path(tcga_fuente_directorio_path,
            "TcgaTargetGtex_gene_expected_count.gz")

# Primero identificamos la cabecera de la columna de genes (suele ser 'sample').
testiculo_cabecera_vector <-
  testiculo_matriz_cruda_archivo_path %>%
  fread(nrows = 0) %>%
  names()

columna_id_genes_char <- testiculo_cabecera_vector[1]

message("Gene identifier column header: '", columna_id_genes_char, "'")

# Comprobamos que las 411 muestras estan en la cabecera de la matriz cruda.
samples_objetivo_vector <- testiculo_grupos_tbl_inicial$sample

samples_no_encontradas_vector <-
  setdiff(samples_objetivo_vector, testiculo_cabecera_vector)

if (length(samples_no_encontradas_vector) > 0) {
  message("Warning: ", length(samples_no_encontradas_vector),
          " sample IDs from grupos_testiculo.csv are not present in the raw matrix header")
  print(head(samples_no_encontradas_vector, 5))
}

samples_efectivas_vector <-
  intersect(samples_objetivo_vector, testiculo_cabecera_vector)

message("Sample IDs to be extracted from raw matrix: ",
        length(samples_efectivas_vector), " / ",
        length(samples_objetivo_vector))

# Carga selectiva de columnas.
message("Loading raw matrix (column-restricted load to minimize memory)...")

testiculo_matriz_cruda_tbl <-
  testiculo_matriz_cruda_archivo_path %>%
  fread(select = c(columna_id_genes_char, samples_efectivas_vector)) %>%
  as_tibble()

message("Raw matrix loaded: ",
        nrow(testiculo_matriz_cruda_tbl), " genes x ",
        ncol(testiculo_matriz_cruda_tbl) - 1, " samples")


# 05. ASIGNAR GRUPO DIAGNOSTICO Y ORIGEN _______________________________________

# El grupo diagnostico y el origen del dataset vienen directos de
# grupos_testiculo.csv. Dejamos anotado que el origen esta completamente
# confundido con el grupo.

testiculo_metadatos_tbl_curado <-
  testiculo_grupos_tbl_inicial %>%
  dplyr::rename(sample_id          = "sample",
                grupo_diagnostico  = "grupo",
                origen_dataset     = "origen") %>%
  filter(sample_id %in% samples_efectivas_vector) %>%
  mutate(
    tejido_control_subtipo = case_when(
      grupo_diagnostico == "Normal"
      ~ "Healthy testis (GTEx donor)",
      TRUE
      ~ NA_character_
    ),
    histologia_subtipo = case_when(
      grupo_diagnostico == "Tumor"
      ~ "Testicular germ cell tumor (TCGA-TGCT)",
      TRUE
      ~ NA_character_
    )
  )

message("Standardized diagnostic distribution:")
print(table(testiculo_metadatos_tbl_curado$grupo_diagnostico,
            testiculo_metadatos_tbl_curado$origen_dataset))


# 06. ANOTACION ENSEMBL A SIMBOLO ______________________________________________

# Los IDs nativos llevan sufijo de version (ENSG00000169083.15). Quitamos el
# sufijo para alinear con el espacio Ensembl comun a los otros datasets y luego
# unimos con el probemap para pegar los simbolos HGNC.

testiculo_genes_anotacion_tbl_inicial <-
  tibble(
    gene_id_native = testiculo_matriz_cruda_tbl[[columna_id_genes_char]]
  ) %>%
  mutate(
    ensembl_id = str_remove(gene_id_native, "\\.\\d+$")
  ) %>%
  left_join(
    testiculo_probemap_tbl %>%
      mutate(ensembl_id = str_remove(gene_id_native, "\\.\\d+$")) %>%
      select(ensembl_id, gene_symbol),
    by = "ensembl_id"
  )

n_genes_sin_simbolo_int <-
  sum(is.na(testiculo_genes_anotacion_tbl_inicial$gene_symbol) |
        testiculo_genes_anotacion_tbl_inicial$gene_symbol == "")

message("Genes without HGNC symbol mapping: ",
        n_genes_sin_simbolo_int, " / ",
        nrow(testiculo_genes_anotacion_tbl_inicial))


# 07. CONTROL DE CALIDAD DE LOS GENES __________________________________________

# Confirmamos que los valores estan en escala log2(conteo + 1) y que el receptor
# de androgenos (AR) esta presente.

testiculo_valores_muestra_matriz <-
  testiculo_matriz_cruda_tbl %>%
  select(2:min(6, ncol(.))) %>%
  as.matrix()

message("Value range across first 5 samples:")
message("  min = ",
        round(min(testiculo_valores_muestra_matriz, na.rm = TRUE), 2))
message("  max = ",
        round(max(testiculo_valores_muestra_matriz, na.rm = TRUE), 2))
message("  median = ",
        round(median(testiculo_valores_muestra_matriz, na.rm = TRUE), 2))

testiculo_valores_log2_lgl <-
  max(testiculo_valores_muestra_matriz, na.rm = TRUE) < 30 &&
  min(testiculo_valores_muestra_matriz, na.rm = TRUE) >= 0

message("Values consistent with log2(count + 1) scale: ",
        testiculo_valores_log2_lgl)

# Comprobacion rapida de AR.
ar_ensembl_id_char <- "ENSG00000169083"

testiculo_ar_presente_lgl <-
  any(testiculo_genes_anotacion_tbl_inicial$ensembl_id == ar_ensembl_id_char)

testiculo_ar_symbol_lgl <-
  "AR" %in% testiculo_genes_anotacion_tbl_inicial$gene_symbol

message("AR (", ar_ensembl_id_char, ") present in matrix: ",
        testiculo_ar_presente_lgl)
message("AR symbol mapped from probemap: ",
        testiculo_ar_symbol_lgl)


# 08. ENSAMBLAR LA MATRIZ ESTANDARIZADA ________________________________________

# Construimos la matriz numerica indexada por el Ensembl sin sufijo. Los simbolos
# quedan en la tabla de anotacion; los nombres de fila de la matriz usan el
# Ensembl para ser coherentes con GSE134073, GSE217260 y GSE214940, que tambien
# van indexados por Ensembl.

testiculo_expresion_matriz_final <-
  testiculo_matriz_cruda_tbl %>%
  select(-all_of(columna_id_genes_char)) %>%
  as.matrix()

rownames(testiculo_expresion_matriz_final) <-
  testiculo_genes_anotacion_tbl_inicial$ensembl_id

# Reordenamos los metadatos para que sigan el mismo orden de columnas que la matriz.
testiculo_metadatos_tbl_final <-
  testiculo_metadatos_tbl_curado %>%
  dplyr::slice(match(colnames(testiculo_expresion_matriz_final), sample_id)) %>%
  select(sample_id,
         grupo_diagnostico,
         origen_dataset,
         tejido_control_subtipo,
         histologia_subtipo)

message("Standardized matrix assembled: ",
        nrow(testiculo_expresion_matriz_final), " genes x ",
        ncol(testiculo_expresion_matriz_final), " samples")


# 09. DESCRIPTOR DEL DATASET ___________________________________________________

testiculo_descriptor_lista <-
  list(
    gse_id                       = dataset_id_char,
    organo                       = organo_objetivo_char,
    plataforma_tipo              = plataforma_tipo_char,
    gene_id_convention           = "Ensembl (Gencode v23, version-stripped)",
    matrix_value_nature          = "log2(count + 1) — requires reversal for DESeq2",
    dge_branch_recommended       = "DESeq2 with ~ grupo_diagnostico (after log2 reversal)",
    requiere_reversion_log2      = TRUE,
    n_total                      = ncol(testiculo_expresion_matriz_final),
    n_tumor                      = sum(testiculo_metadatos_tbl_final$grupo_diagnostico
                                       == "Tumor", na.rm = TRUE),
    n_normal                     = sum(testiculo_metadatos_tbl_final$grupo_diagnostico
                                       == "Normal", na.rm = TRUE),
    n_sin_etiqueta               = sum(is.na(testiculo_metadatos_tbl_final$grupo_diagnostico)),
    comparator_subtype           = "Healthy testis (GTEx donors)",
    histology                    = "Testicular germ cell tumor (TCGA-TGCT)",
    paired_design                = FALSE,
    role_in_meta_analysis        = "Primary testicular cohort",
    origen_grupo_confunded       = TRUE,
    origen_grupo_caveat          = "Tumor identity is perfectly confounded with dataset origin (TCGA = Tumor, GTEx = Normal); origen cannot be modelled as a covariate",
    harmonization_pipeline       = "UCSC Toil recompute (uniform alignment and quantification across TCGA + GTEx)",
    n_genes_without_symbol       = n_genes_sin_simbolo_int,
    feature_filtering            = "None at standardization stage; DESeq2 applies the canonical low-expression filter downstream",
    fecha_estandarizacion        = format(Sys.Date(), "%Y-%m-%d")
  )


# 10. EXPORTAR EL RDS ESTANDARIZADO ____________________________________________

testiculo_dataset_estandarizado_lista <-
  list(
    conteos_matriz         = testiculo_expresion_matriz_final,
    genes_anotacion_tbl    = testiculo_genes_anotacion_tbl_inicial,
    muestras_metadatos_tbl = testiculo_metadatos_tbl_final,
    dataset_descriptor     = testiculo_descriptor_lista
  )

testiculo_rds_destino_path <-
  file.path(
    datasets_estandarizados_directorio_path,
    paste0(dataset_id_char, "_estandarizado.rds")
  )

saveRDS(testiculo_dataset_estandarizado_lista, testiculo_rds_destino_path)

message("Standardized object exported to: ",
        basename(testiculo_rds_destino_path))


# 11. RESUMEN __________________________________________________________________

message("")
message("Standardization completed for ", dataset_id_char)
message("Organ: ", testiculo_descriptor_lista$organo)
message("Platform: ", testiculo_descriptor_lista$plataforma_tipo)
message("Role in meta-analysis: ",
        testiculo_descriptor_lista$role_in_meta_analysis)
message("Matrix nature: ",
        testiculo_descriptor_lista$matrix_value_nature)
message("DGE branch recommended: ",
        testiculo_descriptor_lista$dge_branch_recommended)
message("Harmonization pipeline: ",
        testiculo_descriptor_lista$harmonization_pipeline)
message("Comparator subtype: ",
        testiculo_descriptor_lista$comparator_subtype)
message("Histology: ",
        testiculo_descriptor_lista$histology)
message("Features retained: ",
        nrow(testiculo_dataset_estandarizado_lista$conteos_matriz),
        " Ensembl genes (version-stripped)")
message("Samples retained: ",
        testiculo_descriptor_lista$n_total,
        " (Tumor = ", testiculo_descriptor_lista$n_tumor,
        ", Normal = ", testiculo_descriptor_lista$n_normal, ")")
message("AR presence verified: ", testiculo_ar_presente_lgl,
        " (symbol mapping: ", testiculo_ar_symbol_lgl, ")")
message("LIMITATION: ", testiculo_descriptor_lista$origen_grupo_caveat)
message("Output: ", testiculo_rds_destino_path)
