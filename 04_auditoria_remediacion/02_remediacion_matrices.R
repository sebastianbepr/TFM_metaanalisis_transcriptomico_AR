# ______________________________________________________________________________
# 04_02  Remediacion de matrices (simbolo HGNC, escala, filtro por varianza)
#
# Que hace: deja los datasets estandarizados en un estado tecnicamente homogeneo:
#           reindexa a simbolo HGNC, colapsa duplicados, unifica la escala y
#           aplica un filtro unico por varianza antes del meta-analisis.
# Entra:    Datasets_Estandarizados/*.rds
# Sale:     Datasets_Remediados/<GSE>_remediado.rds, reporte_remediacion.csv,
#           espacio_genico_comun.csv
# ______________________________________________________________________________

# __ INDEX ___
# 01. Setup, librerias y parametros de la etapa
# 02. Rutas y descubrimiento de archivos
# 03. Funciones auxiliares generales
# 04. B1 - Verificacion de integridad y seleccion de clave
# 05. B2 - Reindexacion a simbolo genico y colapso de duplicados
# 06. B3 - Transformacion de escala
# 07. B4 - Filtro de expresion armonizado
# 08. Pipeline de remediacion por dataset
# 09. Ejecucion sobre los 11 datasets
# 10. Exportacion y resumen final


# 01. SETUP ____________________________________________________________________

library(tidyverse)
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Parametros de la etapa. Centralizados para que ningun umbral quede escondido
# dentro del codigo.

parametros_remediacion_lista <- list(
  # Cobertura minima aceptable al mapear filas de la matriz contra la tabla de
  # anotacion. Por debajo de este valor el dataset se marca para revision.
  cobertura_minima_mapeo = 0.80,

  # Percentil de varianza por debajo del cual se descartan los genes. Se filtra
  # por varianza y no por media porque GSE40355 contiene log-ratios de dos
  # colores centrados en cero, donde una media baja no indica baja expresion.
  percentil_varianza     = 0.10,

  # Numero maximo de valores muestreados para inferir la escala de una matriz.
  n_valores_muestreo     = 20000,

  semilla_aleatoria      = 1234
)

set.seed(parametros_remediacion_lista$semilla_aleatoria)


# 02. RUTAS Y DESCUBRIMIENTO DE ARCHIVOS _______________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_datasets_chr   <- file.path(ruta_base_chr, "Datasets_Estandarizados")
ruta_remediados_chr <- file.path(ruta_base_chr, "Datasets_Remediados")
ruta_salida_chr     <- file.path(ruta_base_chr, "Meta_Analisis_Resultados")

for (ruta_chr in c(ruta_remediados_chr, ruta_salida_chr)) {
  if (!dir.exists(ruta_chr)) {
    dir.create(ruta_chr, recursive = TRUE)
  }
}

# Se listan los archivos ignorando mayusculas y descartando cualquier directorio
# que pudiera coincidir con el patron.

listar_archivos_rds <- function(ruta_carpeta_chr) {

  if (!dir.exists(ruta_carpeta_chr)) {
    stop("La carpeta no existe: ", ruta_carpeta_chr)
  }

  archivos_chr <- list.files(ruta_carpeta_chr,
                             pattern     = "\\.rds$",
                             full.names  = TRUE,
                             ignore.case = TRUE)

  archivos_chr <- archivos_chr[!file.info(archivos_chr)$isdir]

  if (length(archivos_chr) == 0) {
    stop("No se encontraron archivos .rds en: ", ruta_carpeta_chr)
  }

  archivos_chr
}

archivos_datasets_chr <- listar_archivos_rds(ruta_datasets_chr)

message("Datasets a remediar: ", length(archivos_datasets_chr))


# 03. FUNCIONES AUXILIARES GENERALES ___________________________________________

# Extraccion segura de campos del descriptor. Devuelve NA en lugar de fallar si
# el campo no existe en algun dataset.

extraer_campo <- function(lista_origen, campo_chr) {

  valor <- lista_origen[[campo_chr]]

  if (is.null(valor) || length(valor) == 0) {
    return(NA_character_)
  }

  as.character(valor[1])
}


# Elimina el sufijo de version de identificadores Ensembl (ENSG00000141510.11
# pasa a ENSG00000141510). Es una causa frecuente de fallo de mapeo entre la
# matriz y la tabla de anotacion.

limpiar_version_ensembl <- function(identificadores_chr) {
  sub("\\.\\d+$", "", identificadores_chr)
}


# Convierte a NA las cadenas vacias y los marcadores textuales que algunos
# pipelines dejan en lugar de un valor ausente real.

normalizar_ausentes <- function(valores_chr) {

  valores_chr <- as.character(valores_chr)
  valores_chr[valores_chr %in% c("", " ", "NA", "N/A", "-", "---")] <- NA_character_

  valores_chr
}


# Muestreo de valores de la matriz para inferir su escala sin recorrerla
# completa, lo que seria costoso en las matrices de 60.000 genes.

muestrear_valores <- function(expresion_matriz) {

  n_total_int   <- length(expresion_matriz)
  n_muestra_int <- min(parametros_remediacion_lista$n_valores_muestreo,
                       n_total_int)

  indices_int <- sample.int(n_total_int, n_muestra_int)
  valores_num <- as.numeric(expresion_matriz[indices_int])

  valores_num[is.finite(valores_num)]
}


# Traduccion externa de identificadores Ensembl a simbolo HGNC. Se usa
# org.Hs.eg.db en lugar de biomaRt porque funciona sin conexion y su version
# queda registrada en sessionInfo(), lo que hace el mapeo reproducible.

mapear_ensembl_externo <- function(identificadores_chr) {

  claves_chr    <- limpiar_version_ensembl(identificadores_chr)
  resultado_chr <- rep(NA_character_, length(claves_chr))

  validas_lgl <- !is.na(claves_chr) & nzchar(claves_chr)

  if (!any(validas_lgl)) {
    return(resultado_chr)
  }

  simbolos_chr <- suppressMessages(suppressWarnings(
    AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys      = unique(claves_chr[validas_lgl]),
      column    = "SYMBOL",
      keytype   = "ENSEMBL",
      multiVals = "first"
    )
  ))

  resultado_chr[validas_lgl] <- simbolos_chr[claves_chr[validas_lgl]]

  normalizar_ausentes(resultado_chr)
}


# 04. B1 - VERIFICACION DE INTEGRIDAD Y SELECCION DE CLAVE _____________________

# Columnas de la tabla de anotacion que pueden actuar como clave. No se asume
# cual corresponde a los rownames de la matriz: la auditoria mostro que en
# varios datasets el indice es el simbolo mientras gene_id_native guarda el
# identificador Ensembl, es decir, las columnas estan cruzadas.

columnas_clave_chr <- c("gene_id_native", "gene_symbol", "ensembl_id")


# Evalua la cobertura de los rownames contra una columna concreta, con y sin
# sufijo de version Ensembl.

evaluar_clave <- function(identificadores_chr, anotacion_tbl, columna_chr) {

  if (!columna_chr %in% names(anotacion_tbl)) {
    return(tibble(columna = columna_chr, limpiar_version = FALSE, cobertura = 0))
  }

  valores_chr <- as.character(anotacion_tbl[[columna_chr]])

  tibble(
    columna         = rep(columna_chr, 2),
    limpiar_version = c(FALSE, TRUE),
    cobertura       = c(
      mean(identificadores_chr %in% valores_chr),
      mean(limpiar_version_ensembl(identificadores_chr) %in%
             limpiar_version_ensembl(valores_chr))
    )
  )
}


# Selecciona la combinacion de columna y limpieza de version que maximiza la
# cobertura del mapeo. La union posterior se hara siempre por clave y nunca por
# posicion, porque una union posicional en un dataset desalineado asignaria los
# identificadores a filas equivocadas de forma silenciosa.

verificar_integridad <- function(expresion_matriz, anotacion_tbl) {

  identificadores_chr <- rownames(expresion_matriz)

  if (is.null(identificadores_chr)) {
    stop("La matriz de expresion no tiene rownames")
  }

  alineado_lgl <- identical(identificadores_chr,
                            as.character(anotacion_tbl$gene_id_native))

  evaluacion_tbl <- columnas_clave_chr %>%
    map_dfr(~ evaluar_clave(identificadores_chr, anotacion_tbl, .x)) %>%
    arrange(desc(cobertura))

  mejor_tbl <- dplyr::slice(evaluacion_tbl, 1)

  # Proporcion de simbolos utilizables en la tabla de anotacion propia. Si es
  # cero, el dataset dependera por completo de la traduccion externa.
  simbolos_disponibles_num <- 0

  if ("gene_symbol" %in% names(anotacion_tbl)) {
    simbolos_disponibles_num <- mean(
      !is.na(normalizar_ausentes(anotacion_tbl$gene_symbol))
    )
  }

  list(
    alineado             = alineado_lgl,
    columna_clave        = mejor_tbl$columna,
    limpiar_version      = mejor_tbl$limpiar_version,
    cobertura            = mejor_tbl$cobertura,
    simbolos_disponibles = simbolos_disponibles_num
  )
}


# 05. B2 - REINDEXACION A SIMBOLO GENICO _______________________________________

# Devuelve el vector de simbolos correspondiente a cada fila de la matriz. El
# mapeo se hace por clave con match(), respetando el orden de la matriz, y se
# resuelve en tres pasos sucesivos:
#   1. Traduccion mediante la tabla de anotacion propia del dataset.
#   2. Si el rowname ya es un simbolo, se conserva tal cual.
#   3. Si el rowname es Ensembl y sigue sin simbolo, se traduce con org.Hs.eg.db.

mapear_a_simbolo <- function(expresion_matriz, anotacion_tbl, integridad_lista) {

  claves_matriz_chr <- rownames(expresion_matriz)

  claves_anotacion_chr <- as.character(
    anotacion_tbl[[integridad_lista$columna_clave]]
  )

  if (integridad_lista$limpiar_version) {
    claves_matriz_chr    <- limpiar_version_ensembl(claves_matriz_chr)
    claves_anotacion_chr <- limpiar_version_ensembl(claves_anotacion_chr)
  }

  # Paso 1: anotacion propia del dataset.
  posiciones_int <- match(claves_matriz_chr, claves_anotacion_chr)
  simbolos_chr   <- normalizar_ausentes(anotacion_tbl$gene_symbol[posiciones_int])

  # Paso 2: el rowname ya es un simbolo, no un identificador Ensembl.
  propios_lgl <- is.na(simbolos_chr) & !grepl("^ENSG", claves_matriz_chr)
  simbolos_chr[propios_lgl] <- claves_matriz_chr[propios_lgl]

  # Paso 3: traduccion externa de los Ensembl que siguen sin simbolo.
  pendientes_lgl <- is.na(simbolos_chr) & grepl("^ENSG", claves_matriz_chr)
  n_externos_int <- 0

  if (any(pendientes_lgl)) {
    message("    B2 traduciendo ", sum(pendientes_lgl),
            " identificadores Ensembl con org.Hs.eg.db")

    simbolos_chr[pendientes_lgl] <- mapear_ensembl_externo(
      claves_matriz_chr[pendientes_lgl]
    )

    n_externos_int <- sum(pendientes_lgl & !is.na(simbolos_chr))

    message("    B2 traducidos con exito: ", n_externos_int)
  }

  list(
    simbolos    = normalizar_ausentes(simbolos_chr),
    n_externos  = n_externos_int
  )
}


# Colapsa las filas que comparten simbolo conservando, para cada simbolo, la
# fila de mayor expresion media. Es la regla ya empleada en la estandarizacion
# original ("Maximum mean expression per gene symbol") y una de las recomendadas
# en la literatura de meta-analisis para resumir sondas redundantes.

colapsar_por_simbolo <- function(expresion_matriz, simbolos_chr) {

  filas_validas_lgl <- !is.na(simbolos_chr)

  expresion_matriz <- expresion_matriz[filas_validas_lgl, , drop = FALSE]
  simbolos_chr     <- simbolos_chr[filas_validas_lgl]

  if (nrow(expresion_matriz) == 0) {
    stop("Ningun gen pudo mapearse a simbolo. Revisar la tabla de anotacion.")
  }

  medias_num <- rowMeans(expresion_matriz, na.rm = TRUE)

  # Se ordena por media descendente y se conserva la primera aparicion de cada
  # simbolo, que corresponde a la fila de mayor expresion media.
  orden_int <- order(medias_num, decreasing = TRUE)

  expresion_matriz <- expresion_matriz[orden_int, , drop = FALSE]
  simbolos_chr     <- simbolos_chr[orden_int]

  conservar_lgl <- !duplicated(simbolos_chr)

  expresion_matriz <- expresion_matriz[conservar_lgl, , drop = FALSE]
  rownames(expresion_matriz) <- simbolos_chr[conservar_lgl]

  # Se reordena alfabeticamente para que todas las matrices compartan el mismo
  # criterio de ordenacion de filas.
  expresion_matriz[order(rownames(expresion_matriz)), , drop = FALSE]
}


# 06. B3 - TRANSFORMACION DE ESCALA ____________________________________________

# La decision se toma sobre las propiedades numericas de la matriz y no sobre el
# texto del descriptor, porque la auditoria mostro que ambos pueden discrepar
# (GSE114740 declaraba log2(RPKM + 1) pero contenia RPKM lineal).

determinar_transformacion <- function(expresion_matriz) {

  valores_num <- muestrear_valores(expresion_matriz)

  if (length(valores_num) == 0) {
    return("indeterminada")
  }

  maximo_num        <- max(valores_num)
  todos_enteros_lgl <- all(abs(valores_num - round(valores_num)) < 1e-8)

  # Rango compatible con datos ya logaritmados.
  if (maximo_num <= 30) {
    return("ninguna")
  }

  # Conteos crudos: enteros con rango muy amplio. Requieren normalizacion por
  # tamano de libreria antes de logaritmar.
  if (todos_enteros_lgl && maximo_num > 100) {
    return("log2_cpm")
  }

  # Unidades lineales ya normalizadas (RPKM, FPKM, TPM, conteos normalizados).
  # El tamano de libreria ya esta corregido, solo falta el logaritmo.
  "log2_mas_uno"
}


aplicar_transformacion <- function(expresion_matriz, tipo_chr) {

  if (tipo_chr == "ninguna") {
    return(expresion_matriz)
  }

  # Los valores ausentes se tratan como no detectados antes de transformar. Es
  # relevante solo en GSE221107, donde la auditoria encontro 32 NA.
  if (anyNA(expresion_matriz)) {
    message("    B3 sustituyendo ", sum(is.na(expresion_matriz)),
            " valores NA por 0")
    expresion_matriz[is.na(expresion_matriz)] <- 0
  }

  if (tipo_chr == "log2_mas_uno") {
    return(log2(expresion_matriz + 1))
  }

  if (tipo_chr == "log2_cpm") {
    # Normalizacion TMM seguida de log-CPM. Se usa prior.count = 1 para evitar
    # varianzas infladas en los genes de baja expresion.
    lista_dge <- DGEList(counts = expresion_matriz)
    lista_dge <- calcNormFactors(lista_dge, method = "TMM")

    return(cpm(lista_dge, log = TRUE, prior.count = 1))
  }

  expresion_matriz
}


# 07. B4 - FILTRO DE EXPRESION ARMONIZADO ______________________________________

# Se sustituyen los once filtros heterogeneos originales por una regla unica.
# El filtro es por varianza y no por media: en GSE40355 los valores son
# log-ratios de dos colores centrados en cero, de modo que una media baja indica
# ausencia de cambio frente a la referencia y no baja expresion. La varianza es
# interpretable de la misma forma en log-ratios, log-intensidades y log-CPM.
#
# El filtro no utiliza las etiquetas de grupo, por lo que es un filtro
# independiente y no sesga los contrastes posteriores.

filtrar_por_varianza <- function(expresion_matriz, percentil_num) {

  varianzas_num <- apply(expresion_matriz, 1, var, na.rm = TRUE)

  # Los genes sin variacion no aportan informacion y hacen fallar el ajuste
  # lineal, de modo que se descartan siempre.
  conservar_lgl <- is.finite(varianzas_num) & varianzas_num > 0

  if (percentil_num > 0) {
    umbral_num <- quantile(varianzas_num[conservar_lgl],
                           probs = percentil_num,
                           na.rm = TRUE)

    conservar_lgl <- conservar_lgl & varianzas_num >= umbral_num
  }

  expresion_matriz[conservar_lgl, , drop = FALSE]
}


# 08. PIPELINE DE REMEDIACION POR DATASET ______________________________________

remediar_dataset <- function(ruta_archivo_chr) {

  # Guarda defensiva: el pipeline debe fallar con un mensaje legible si recibe
  # una carpeta o una ruta inexistente.
  if (dir.exists(ruta_archivo_chr)) {
    stop("Se recibio una carpeta en lugar de un archivo: ", ruta_archivo_chr)
  }

  if (!file.exists(ruta_archivo_chr)) {
    stop("Archivo inexistente: ", ruta_archivo_chr)
  }

  dataset_lista <- readRDS(ruta_archivo_chr)

  descriptor_lista <- dataset_lista$dataset_descriptor
  expresion_matriz <- dataset_lista$conteos_matriz
  metadatos_tbl    <- dataset_lista$muestras_metadatos_tbl
  anotacion_tbl    <- dataset_lista$genes_anotacion_tbl

  gse_id_chr <- extraer_campo(descriptor_lista, "gse_id")

  message("\n[", gse_id_chr, "] Iniciando remediacion")

  n_genes_inicial_int <- nrow(expresion_matriz)

  # --- B1: integridad y seleccion de clave ------------------------------------

  integridad_lista <- verificar_integridad(expresion_matriz, anotacion_tbl)

  message("    B1 clave usada: ", integridad_lista$columna_clave,
          " | cobertura: ", round(100 * integridad_lista$cobertura, 1), "%")

  message("    B1 simbolos disponibles en anotacion propia: ",
          round(100 * integridad_lista$simbolos_disponibles, 1), "%")

  if (integridad_lista$limpiar_version) {
    message("    B1 se aplico limpieza de sufijos de version Ensembl")
  }

  if (integridad_lista$cobertura <
      parametros_remediacion_lista$cobertura_minima_mapeo) {
    message("    B1 AVISO: cobertura por debajo del umbral, revisar dataset")
  }

  # --- B2: reindexacion -------------------------------------------------------

  mapeo_lista <- mapear_a_simbolo(expresion_matriz, anotacion_tbl,
                                  integridad_lista)

  simbolos_chr      <- mapeo_lista$simbolos
  n_sin_simbolo_int <- sum(is.na(simbolos_chr))

  expresion_matriz <- colapsar_por_simbolo(expresion_matriz, simbolos_chr)

  n_genes_tras_colapso_int <- nrow(expresion_matriz)

  message("    B2 genes sin simbolo descartados: ", n_sin_simbolo_int)
  message("    B2 genes tras colapso por simbolo: ", n_genes_tras_colapso_int)

  # --- B3: escala -------------------------------------------------------------

  tipo_transformacion_chr <- determinar_transformacion(expresion_matriz)

  message("    B3 transformacion aplicada: ", tipo_transformacion_chr)

  expresion_matriz <- aplicar_transformacion(expresion_matriz,
                                             tipo_transformacion_chr)

  # --- B4: filtro -------------------------------------------------------------

  expresion_matriz <- filtrar_por_varianza(
    expresion_matriz,
    parametros_remediacion_lista$percentil_varianza
  )

  n_genes_final_int <- nrow(expresion_matriz)

  message("    B4 genes tras filtro de varianza: ", n_genes_final_int)

  # --- Vector de grupo para la etapa C ----------------------------------------

  # DExMA exige un vector 0/1 donde 0 es el grupo control. Se construye aqui
  # para que la etapa C no dependa de la codificacion textual original.
  grupo_binario_int <- as.integer(metadatos_tbl$grupo_diagnostico == "Tumor")
  names(grupo_binario_int) <- metadatos_tbl$sample_id

  # Se reordena el vector segun las columnas de la matriz para garantizar la
  # correspondencia muestra a muestra.
  grupo_binario_int <- grupo_binario_int[colnames(expresion_matriz)]

  # --- Objeto remediado -------------------------------------------------------

  dataset_remediado_lista <- list(
    expresion_matriz       = expresion_matriz,
    grupo_binario          = grupo_binario_int,
    muestras_metadatos_tbl = metadatos_tbl,
    dataset_descriptor     = descriptor_lista,
    remediacion_descriptor = list(
      columna_clave           = integridad_lista$columna_clave,
      limpieza_version        = integridad_lista$limpiar_version,
      cobertura_mapeo         = integridad_lista$cobertura,
      simbolos_externos       = mapeo_lista$n_externos,
      transformacion_aplicada = tipo_transformacion_chr,
      percentil_varianza      = parametros_remediacion_lista$percentil_varianza,
      n_genes_inicial         = n_genes_inicial_int,
      n_genes_final           = n_genes_final_int,
      fecha_remediacion       = as.character(Sys.Date())
    )
  )

  saveRDS(
    dataset_remediado_lista,
    file.path(ruta_remediados_chr, paste0(gse_id_chr, "_remediado.rds"))
  )

  # --- Fila del reporte -------------------------------------------------------

  tibble(
    gse_id               = gse_id_chr,
    organo               = extraer_campo(descriptor_lista, "organo"),
    convencion_original  = extraer_campo(descriptor_lista,
                                         "gene_id_convention"),
    alineado_original    = integridad_lista$alineado,
    columna_clave        = integridad_lista$columna_clave,
    limpieza_version     = integridad_lista$limpiar_version,
    cobertura_mapeo      = round(integridad_lista$cobertura, 3),
    simbolos_propios_pct = round(100 * integridad_lista$simbolos_disponibles, 1),
    simbolos_externos    = mapeo_lista$n_externos,
    n_sin_simbolo        = n_sin_simbolo_int,
    transformacion       = tipo_transformacion_chr,
    n_genes_inicial      = n_genes_inicial_int,
    n_genes_tras_colapso = n_genes_tras_colapso_int,
    n_genes_final        = n_genes_final_int,
    pct_conservado       = round(100 * n_genes_final_int / n_genes_inicial_int,
                                 1),
    valor_minimo_final   = round(min(expresion_matriz, na.rm = TRUE), 2),
    valor_maximo_final   = round(max(expresion_matriz, na.rm = TRUE), 2),
    n_muestras           = ncol(expresion_matriz),
    grupo_completo       = !anyNA(grupo_binario_int)
  )
}


# 09. EJECUCION SOBRE LOS 11 DATASETS __________________________________________

message("\n=== INICIANDO ETAPA B: REMEDIACION DE MATRICES ===")

reporte_remediacion_tbl <- archivos_datasets_chr %>%
  map_dfr(remediar_dataset) %>%
  arrange(organo, gse_id)


# 10. EXPORTACION Y RESUMEN FINAL ______________________________________________

write_csv(reporte_remediacion_tbl,
          file.path(ruta_salida_chr, "reporte_remediacion.csv"))


# Diagnostico del espacio genico comun por grupo tumoral. Es la cifra que
# determina cuanta informacion se pierde con la interseccion estricta y, por
# tanto, cuanto aporta la imputacion de genes faltantes en la etapa siguiente.

calcular_espacio_comun <- function(organo_chr) {

  archivos_grupo_chr <- reporte_remediacion_tbl %>%
    filter(organo == organo_chr) %>%
    pull(gse_id) %>%
    paste0("_remediado.rds") %>%
    file.path(ruta_remediados_chr, .)

  listas_genes <- archivos_grupo_chr %>%
    map(~ rownames(readRDS(.x)$expresion_matriz))

  genes_union_chr        <- purrr::reduce(listas_genes, union)
  genes_interseccion_chr <- purrr::reduce(listas_genes, intersect)

  tibble(
    organo             = organo_chr,
    n_datasets         = length(listas_genes),
    genes_union        = length(genes_union_chr),
    genes_interseccion = length(genes_interseccion_chr),
    pct_interseccion   = round(100 * length(genes_interseccion_chr) /
                                 length(genes_union_chr), 1)
  )
}

espacio_comun_tbl <- reporte_remediacion_tbl$organo %>%
  unique() %>%
  map_dfr(calcular_espacio_comun)

write_csv(espacio_comun_tbl,
          file.path(ruta_salida_chr, "espacio_genico_comun.csv"))


message("\n________________ RESUMEN DE LA REMEDIACION ________________\n")

message("Datasets remediados: ", nrow(reporte_remediacion_tbl))

message("\nClaves de anotacion utilizadas:")

reporte_remediacion_tbl %>%
  count(columna_clave, name = "n") %>%
  pmap(function(columna_clave, n) {
    message("  ", columna_clave, ": ", n)
  }) %>%
  invisible()

message("\nTransformaciones de escala aplicadas:")

reporte_remediacion_tbl %>%
  count(transformacion, name = "n") %>%
  pmap(function(transformacion, n) {
    message("  ", transformacion, ": ", n)
  }) %>%
  invisible()

message("\nDatasets que requirieron traduccion externa: ",
        sum(reporte_remediacion_tbl$simbolos_externos > 0))

message("Total de identificadores traducidos con org.Hs.eg.db: ",
        sum(reporte_remediacion_tbl$simbolos_externos))

message("Cobertura de mapeo minima observada: ",
        round(100 * min(reporte_remediacion_tbl$cobertura_mapeo), 1), "%")

message("\nEspacio genico por grupo tumoral:")

espacio_comun_tbl %>%
  pmap(function(organo, n_datasets, genes_union, genes_interseccion,
                pct_interseccion) {
    message("  ", organo, " (", n_datasets, " datasets): union ",
            genes_union, " | interseccion ", genes_interseccion,
            " (", pct_interseccion, "%)")
  }) %>%
  invisible()

message("\nDatasets con vector de grupo incompleto: ",
        sum(!reporte_remediacion_tbl$grupo_completo))

message("\nMatrices remediadas guardadas en: ", ruta_remediados_chr)
message("Tablas exportadas en: ", ruta_salida_chr)
message("___________________________________________________________\n")
