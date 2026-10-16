# ______________________________________________________________________________
# 05_01  Construccion de los objetos DExMA por grupo y control de heterogeneidad
#
# Que hace: agrupa los datasets remediados por organo, construye un objectMA de
#           DExMA por cada grupo tumoral y evalua la heterogeneidad (I2) de cada
#           grupo antes del meta-analisis.
# Entra:    Datasets_Remediados/*_remediado.rds
# Sale:     Meta_Analisis_Objetos/objetos_ma_por_grupo.rds,
#           heterogeneidad_por_grupo.csv
# ______________________________________________________________________________

# __ INDEX ___
# 01. Setup, librerias y parametros de la etapa
# 02. Rutas y carga de los datasets remediados
# 03. Agrupacion por organo
# 04. Construccion del objectMA por grupo tumoral
# 05. Control de calidad: escala logaritmica y heterogeneidad
# 06. Definicion de la rejilla de meta-analisis
# 07. Ejecucion del meta-analisis (tres brazos, tres metodos)
# 08. Recuperacion de la direccionalidad mediante tamanos de efecto
# 09. Exportacion de resultados y objetos intermedios
# 10. Resumen final


# 01. SETUP ____________________________________________________________________

library(tidyverse)
library(DExMA)

# La preferencia de funciones se fija explicitamente porque DExMA arrastra la
# pila de Bioconductor, que enmascara varias funciones de dplyr y purrr.
select   <- dplyr::select
filter   <- dplyr::filter
count    <- dplyr::count
slice    <- dplyr::slice
rename   <- dplyr::rename
mutate   <- dplyr::mutate
arrange  <- dplyr::arrange
summarise <- dplyr::summarise
group_by <- dplyr::group_by
map      <- purrr::map
map_chr  <- purrr::map_chr
map_dbl  <- purrr::map_dbl
map_dfr  <- purrr::map_dfr
imap     <- purrr::imap
imap_dfr <- purrr::imap_dfr
pmap     <- purrr::pmap
pmap_dfr <- purrr::pmap_dfr
keep     <- purrr::keep
reduce   <- purrr::reduce

parametros_meta_lista <- list(
  # Proporcion minima de datasets en que un gen debe estar presente para el
  # brazo de "proporcion minima".
  proporcion_minima     = 2 / 3,

  # Numero de vecinos para la imputacion sampleKNN de DExMA.
  k_vecinos_imputacion  = 7,

  # Umbral de significancia sobre el p-valor combinado, tras correccion FDR.
  umbral_fdr            = 0.05,

  semilla_aleatoria     = 1234
)

set.seed(parametros_meta_lista$semilla_aleatoria)


# 02. RUTAS Y CARGA DE LOS DATASETS REMEDIADOS _________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_remediados_chr <- file.path(ruta_base_chr, "Datasets_Remediados")
ruta_salida_chr     <- file.path(ruta_base_chr, "Meta_Analisis_Resultados")
ruta_objetos_chr    <- file.path(ruta_base_chr, "Meta_Analisis_Objetos")

if (!dir.exists(ruta_objetos_chr)) {
  dir.create(ruta_objetos_chr, recursive = TRUE)
}

# Se listan los archivos con un patron laxo y se filtra la extension aparte,
# para evitar fallos de coincidencia por mayusculas o variaciones del sufijo.

archivos_remediados_chr <- list.files(ruta_remediados_chr,
                                      pattern     = "remediado",
                                      full.names  = TRUE,
                                      ignore.case = TRUE)

archivos_remediados_chr <- archivos_remediados_chr[
  grepl("\\.rds$", archivos_remediados_chr, ignore.case = TRUE)
]

stopifnot(length(archivos_remediados_chr) == 11)

message("Datasets remediados encontrados: ", length(archivos_remediados_chr))

# Cada objeto remediado se carga junto con su organo, para agruparlos despues.

cargar_remediado <- function(ruta_archivo_chr) {

  remediado_lista <- readRDS(ruta_archivo_chr)

  list(
    gse_id           = remediado_lista$dataset_descriptor$gse_id,
    organo           = remediado_lista$dataset_descriptor$organo,
    expresion_matriz = remediado_lista$expresion_matriz,
    grupo_binario    = remediado_lista$grupo_binario
  )
}

datasets_remediados_lista <- archivos_remediados_chr %>%
  map(cargar_remediado)

names(datasets_remediados_lista) <- datasets_remediados_lista %>%
  map_chr(~ .x$gse_id)


# 03. AGRUPACION POR ORGANO ____________________________________________________

organos_chr <- datasets_remediados_lista %>%
  map_chr(~ .x$organo) %>%
  unique()

message("Grupos tumorales: ", paste(organos_chr, collapse = ", "))


# 04. CONSTRUCCION DEL objectMA POR GRUPO TUMORAL ______________________________

# createObjectMA() de DExMA exige, para cada estudio, la matriz de expresion y
# un phenodata con una columna de condicion (no un vector 0/1 suelto). El
# phenodata se construye a partir del vector binario preparado en la remediacion,
# con nombres de fila iguales a las columnas de la matriz para que DExMA alinee
# las muestras por nombre.

construir_objectma_grupo <- function(organo_chr) {

  message("\n[", organo_chr, "] Construyendo objectMA")

  datasets_grupo_lista <- datasets_remediados_lista %>%
    keep(~ .x$organo == organo_chr)

  message("    Datasets en el grupo: ",
          paste(map_chr(datasets_grupo_lista, ~ .x$gse_id), collapse = ", "))

  matrices_lista <- datasets_grupo_lista %>%
    map(~ .x$expresion_matriz)

  # El vector 0/1 se convierte en un data.frame con etiquetas textuales. Los
  # nombres de fila coinciden con las columnas de la matriz correspondiente.
  phenodatas_lista <- datasets_grupo_lista %>%
    map(function(dataset_lista) {
      condicion_chr <- if_else(dataset_lista$grupo_binario == 1,
                               "Tumor", "Normal")
      data.frame(
        condicion        = condicion_chr,
        row.names        = names(dataset_lista$grupo_binario),
        stringsAsFactors = FALSE
      )
    })

  n_estudios_int <- length(matrices_lista)

  # Cada argumento de grupo se pasa como vector de longitud igual al numero de
  # estudios, para no depender del reciclado interno de DExMA.
  objeto_ma <- createObjectMA(
    listEX    = matrices_lista,
    listPheno = phenodatas_lista,
    namePheno = rep("condicion", n_estudios_int),
    expGroups = rep("Tumor",     n_estudios_int),
    refGroups = rep("Normal",    n_estudios_int)
  )

  objeto_ma
}

# Se construye un objectMA por grupo y se guardan para reutilizarlos sin
# recalcular.

objetos_ma_lista <- organos_chr %>%
  set_names() %>%
  map(construir_objectma_grupo)

saveRDS(objetos_ma_lista,
        file.path(ruta_objetos_chr, "objetos_ma_por_grupo.rds"))


# 05. CONTROL DE CALIDAD: ESCALA Y HETEROGENEIDAD ______________________________

# dataLog() verifica que los datos esten en escala logaritmica, requisito de los
# metodos de combinacion de p-valores. Como la remediacion ya dejo todo en log,
# esta llamada debe confirmarlo sin transformar nada relevante.

message("\n=== CONTROL DE CALIDAD ===")

objetos_ma_lista <- objetos_ma_lista %>%
  imap(function(objeto_ma, organo_chr) {
    message("[", organo_chr, "] Verificando escala logaritmica")
    dataLog(objeto_ma)
  })


# heterogeneityTest() cuantifica la heterogeneidad entre estudios mediante el
# estadistico I2. Se guarda para la seccion de resultados: informa si el grupo
# es homogeneo (I2 < 0.25 en la mayoria de genes) o no, aunque los metodos de
# p-valor no modelen la heterogeneidad directamente.

heterogeneidad_por_grupo <- function(objeto_ma, organo_chr) {

  message("[", organo_chr, "] Calculando heterogeneidad (I2)")

  # heterogeneityTest() devuelve un vector de cuantiles del I2 (no una lista y
  # no el I2 por gen). Se captura y se aplana a un vector nombrado.
  cuantiles_i2 <- heterogeneityTest(objeto_ma)
  cuantiles_i2 <- unlist(cuantiles_i2, use.names = TRUE)

  # Se vuelca cada cuantil como columna, prefijando con el organo.
  tibble(
    organo   = organo_chr,
    cuantil  = names(cuantiles_i2),
    valor_i2 = round(as.numeric(cuantiles_i2), 3)
  )
}

heterogeneidad_tbl <- objetos_ma_lista %>%
  imap_dfr(heterogeneidad_por_grupo)

write_csv(heterogeneidad_tbl,
          file.path(ruta_salida_chr, "heterogeneidad_por_grupo.csv"))
