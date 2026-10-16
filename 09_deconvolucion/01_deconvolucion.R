# ______________________________________________________________________________
# 09_01  Deconvolucion del microambiente tumoral (immunedeconv + entrada CIBERSORTx)
#
# Que hace: estima las fracciones de tipos celulares por dataset con immunedeconv
#           (MCP-counter, xCell, quanTIseq, EPIC) y exporta las mezclas en formato
#           listo para subir a CIBERSORTx.
# Entra:    Datasets_Remediados/*_remediado.rds
# Sale:     Deconvolucion/CIBERSORT_input/<GSE>_mixture.txt,
#           Deconvolucion/Resultados/deconvolucion_fracciones_largo.csv,
#           deconvolucion_resumen.csv
# ______________________________________________________________________________

# __ INDEX ___
# 01. Setup, librerias y parametros
# 02. Rutas y carga de las matrices remediadas
# 03. Tecnologia por dataset (que metodos aplican)
# 04. Preparacion de la matriz de mezcla (lineal, simbolos unicos)
# 05. Exportacion de las mezclas en formato CIBERSORTx
# 06. Deconvolucion con immunedeconv (metodos sin licencia)
# 07. Consolidacion tidy con la condicion tumor/normal
# 08. Exportacion de resultados
# 09. (Opcional) Ingesta de los resultados de CIBERSORTx web


# 01. SETUP ____________________________________________________________________

# Instalacion (una sola vez):
#   install.packages(c("tidyverse","remotes"))
#   remotes::install_github("omnideconv/immunedeconv")
#     -> arrastra MCPcounter, xCell, EPIC y quanTIseq.
#   (xCell puede requerir: remotes::install_github("dviraran/xCell"))

library(tidyverse)
library(immunedeconv)

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
mutate <- dplyr::mutate

parametros_deconv <- list(
  # Umbral por debajo del cual se asume que la matriz esta en escala log2 y se
  # devuelve a escala lineal (CIBERSORT y las firmas esperan valores lineales).
  umbral_log_max = 50,
  semilla        = 1234
)
set.seed(parametros_deconv$semilla)


# 02. RUTAS Y CARGA ____________________________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_remediados_chr <- file.path(ruta_base_chr, "Datasets_Remediados")
ruta_deconv_chr     <- file.path(ruta_base_chr, "Deconvolucion")
ruta_cibersort_chr  <- file.path(ruta_deconv_chr, "CIBERSORT_input")
ruta_result_chr     <- file.path(ruta_deconv_chr, "Resultados")

for (d in c(ruta_deconv_chr, ruta_cibersort_chr, ruta_result_chr)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

archivos_remediados_chr <- list.files(ruta_remediados_chr, pattern = "remediado.*\\.rds$",
                                      full.names = TRUE, ignore.case = TRUE)

cargar_remediado <- function(ruta_archivo_chr) {
  x <- readRDS(ruta_archivo_chr)
  list(
    gse_id  = x$dataset_descriptor$gse_id,
    organo  = x$dataset_descriptor$organo,
    matriz  = x$expresion_matriz,     # genes x muestras (rownames = simbolo)
    grupo   = x$grupo_binario         # 1 = Tumor, 0 = Normal
  )
}

datasets_lista <- archivos_remediados_chr %>% map(cargar_remediado)
names(datasets_lista) <- map_chr(datasets_lista, ~ .x$gse_id)

message("Datasets cargados: ", paste(names(datasets_lista), collapse = ", "))


# 03. TECNOLOGIA POR DATASET ___________________________________________________

# quanTIseq y EPIC asumen RNA-seq (TPM). MCP-counter y xCell son agnosticos a la
# plataforma y se corren en todos. Verifica esta tabla contra tus metadatos:
# la clasificacion solo afecta a que metodos se aplican, no a la mezcla exportada
# para CIBERSORT (que se genera para todos).
tecnologia_tbl <- tribble(
  ~gse_id,             ~tecnologia,
  "GSE104131",         "RNAseq",
  "GSE114740",         "RNAseq",     # <- VERIFICAR (si es microarray, cambiar)
  "GSE134073",         "RNAseq",
  "GSE214940",         "RNAseq",
  "GSE217260",         "RNAseq",
  "GSE221107",         "RNAseq",
  "GSE3218",           "microarray",
  "GSE8607",           "microarray",
  "TCGA_GTEX_Testis",  "RNAseq",
  "GSE7476",           "microarray",
  "GSE40355",          "microarray"
)

tecnologia_de <- function(gse_id_chr) {
  v <- tecnologia_tbl$tecnologia[tecnologia_tbl$gse_id == gse_id_chr]
  if (length(v) == 0) "microarray" else v
}


# 04. PREPARACION DE LA MATRIZ DE MEZCLA _______________________________________

# Devuelve una matriz genes x muestras, en escala lineal, con simbolos unicos
# (si un simbolo se repite, se conserva la fila de mayor expresion media).
preparar_mezcla <- function(matriz_num) {

  matriz_num <- as.matrix(matriz_num)

  # A escala lineal si parece log2.
  if (max(matriz_num, na.rm = TRUE) < parametros_deconv$umbral_log_max) {
    message("    matriz en escala log -> se devuelve a lineal (2^x)")
    matriz_num <- 2 ^ matriz_num
  }

  # NA -> 0 y recorte de negativos residuales.
  matriz_num[is.na(matriz_num)] <- 0
  matriz_num[matriz_num < 0]    <- 0

  # Colapso de simbolos duplicados por mayor media.
  simbolos_chr <- rownames(matriz_num)
  if (any(duplicated(simbolos_chr))) {
    medias_num <- rowMeans(matriz_num, na.rm = TRUE)
    orden_int  <- order(medias_num, decreasing = TRUE)
    matriz_num <- matriz_num[orden_int, , drop = FALSE]
    matriz_num <- matriz_num[!duplicated(rownames(matriz_num)), , drop = FALSE]
  }

  # Se eliminan genes constantes a cero (no aportan a la deconvolucion).
  matriz_num <- matriz_num[rowSums(matriz_num) > 0, , drop = FALSE]
  matriz_num
}


# 05. EXPORTACION DE LAS MEZCLAS EN FORMATO CIBERSORTx _________________________

# CIBERSORTx espera un archivo de texto tabulado: primera columna 'GeneSymbol'
# con el simbolo, seguida de una columna por muestra, en escala lineal.
exportar_cibersort <- function(matriz_num, gse_id_chr) {
  salida_tbl <- as.data.frame(matriz_num) %>%
    rownames_to_column("GeneSymbol")
  ruta_out <- file.path(ruta_cibersort_chr, paste0(gse_id_chr, "_mixture.txt"))
  write_tsv(salida_tbl, ruta_out)
  message("    mezcla CIBERSORT -> ", basename(ruta_out),
          " (", nrow(matriz_num), " genes x ", ncol(matriz_num), " muestras)")
}


# 06. DECONVOLUCION CON immunedeconv ___________________________________________

# Nota metodologica IMPORTANTE:
#   - quanTIseq y EPIC devuelven FRACCIONES verdaderas (suman ~1, comparables
#     entre muestras).
#   - MCP-counter y xCell devuelven PUNTUACIONES relativas de abundancia
#     (NO fracciones): validas para comparar un mismo tipo celular ENTRE muestras,
#     no para comparar tipos entre si dentro de una muestra.
# Se registra el tipo de salida en la columna 'escala'.

correr_metodo <- function(matriz_num, metodo_chr, es_tumor_lgl = TRUE) {
  tryCatch({
    if (metodo_chr %in% c("quantiseq", "epic")) {
      res <- immunedeconv::deconvolute(matriz_num, metodo_chr, tumor = es_tumor_lgl)
    } else {
      res <- immunedeconv::deconvolute(matriz_num, metodo_chr)
    }
    res %>%
      pivot_longer(-cell_type, names_to = "muestra", values_to = "valor") %>%
      mutate(metodo = metodo_chr,
             escala = if_else(metodo_chr %in% c("quantiseq", "epic"),
                              "fraccion", "score_relativo"))
  }, error = function(e) {
    message("    [", metodo_chr, "] FALLO: ", conditionMessage(e))
    tibble()
  })
}

deconvolucionar_dataset <- function(dataset_lista) {

  gse_id_chr <- dataset_lista$gse_id
  message("\n[", gse_id_chr, " | ", dataset_lista$organo, "] deconvolucion")

  matriz_num <- preparar_mezcla(dataset_lista$matriz)

  # Exporta la mezcla para CIBERSORTx (todos los datasets).
  exportar_cibersort(matriz_num, gse_id_chr)

  # Selecciona metodos segun tecnologia.
  tec_chr <- tecnologia_de(gse_id_chr)
  metodos_chr <- if (tec_chr == "RNAseq") {
    c("mcp_counter", "xcell", "quantiseq", "epic")
  } else {
    c("mcp_counter", "xcell")   # microarray: solo metodos agnosticos
  }
  message("    tecnologia: ", tec_chr, " -> metodos: ",
          paste(metodos_chr, collapse = ", "))

  metodos_chr %>%
    map_dfr(~ correr_metodo(matriz_num, .x)) %>%
    mutate(dataset = gse_id_chr, organo = dataset_lista$organo, .before = 1)
}

message("\n=== DECONVOLUCION (immunedeconv) ===")
fracciones_tbl <- datasets_lista %>% map_dfr(deconvolucionar_dataset)


# 07. CONSOLIDACION CON LA CONDICION TUMOR/NORMAL ______________________________

# Se construye el mapa muestra -> condicion desde grupo_binario de cada dataset.
condicion_tbl <- datasets_lista %>%
  imap_dfr(function(d, gse_id_chr) {
    tibble(dataset  = gse_id_chr,
           muestra  = names(d$grupo),
           condicion = if_else(d$grupo == 1, "Tumor", "Normal"))
  })

fracciones_tbl <- fracciones_tbl %>%
  left_join(condicion_tbl, by = c("dataset", "muestra"))


# 08. EXPORTACION ______________________________________________________________

write_csv(fracciones_tbl,
          file.path(ruta_result_chr, "deconvolucion_fracciones_largo.csv"))

# Resumen por dataset/metodo/tipo celular/condicion (media de la fraccion/score).
resumen_tbl <- fracciones_tbl %>%
  group_by(organo, dataset, metodo, escala, cell_type, condicion) %>%
  summarise(media = mean(valor, na.rm = TRUE),
            n     = n(), .groups = "drop")

write_csv(resumen_tbl,
          file.path(ruta_result_chr, "deconvolucion_resumen.csv"))

message("\nResultados exportados en: ", ruta_result_chr)
message("Mezclas para CIBERSORTx en: ", ruta_cibersort_chr)


# 09. (OPCIONAL) INGESTA DE LOS RESULTADOS DE CIBERSORTx WEB ____________________

# Tras subir cada *_mixture.txt a CIBERSORTx (Impute Cell Fractions, firma LM22;
# desactiva la normalizacion por cuantiles en los datasets RNA-seq y activala en
# los microarray), descarga los CSV de resultados a esta carpeta y ejecuta:
#
#   ingerir_cibersortx(file.path(ruta_deconv_chr, "CIBERSORTx_output"))
#
# La funcion los fusiona con la tabla tidy (metodo = "cibersort_lm22").
ingerir_cibersortx <- function(ruta_carpeta_chr) {
  archivos_chr <- list.files(ruta_carpeta_chr, pattern = "\\.csv$", full.names = TRUE)
  if (length(archivos_chr) == 0) {
    message("Sin archivos CIBERSORTx en ", ruta_carpeta_chr); return(invisible(NULL))
  }
  cols_meta_chr <- c("P-value", "Correlation", "RMSE", "P.value", "Absolute score (sig.score)")
  cx_tbl <- archivos_chr %>%
    map_dfr(function(f) {
      gse_id_chr <- str_remove(basename(f), "_.*$")
      read_csv(f, show_col_types = FALSE) %>%
        rename(muestra = 1) %>%
        select(-any_of(cols_meta_chr)) %>%
        pivot_longer(-muestra, names_to = "cell_type", values_to = "valor") %>%
        mutate(dataset = gse_id_chr, metodo = "cibersort_lm22",
               escala = "fraccion")
    })
  cx_tbl <- cx_tbl %>%
    left_join(distinct(select(datasets_lista %>% imap_dfr(~ tibble(
      dataset = .y, organo = .x$organo)), dataset, organo)), by = "dataset") %>%
    left_join(condicion_tbl, by = c("dataset", "muestra"))
  ruta_out <- file.path(ruta_result_chr, "deconvolucion_fracciones_con_cibersort.csv")
  write_csv(bind_rows(fracciones_tbl, cx_tbl), ruta_out)
  message("Fusion con CIBERSORTx -> ", ruta_out)
  invisible(ruta_out)
}

message("\n________ DECONVOLUCION LISTA ________")
message("Metodos sin licencia ejecutados. Sube las mezclas a CIBERSORTx y luego")
message("llama a ingerir_cibersortx() para anadir LM22.")
