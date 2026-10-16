# ______________________________________________________________________________
# 06_02  Analisis de sensibilidad del grupo prostatico (6 vs 5 datasets)
#
# Que hace: repite el meta-analisis prostatico con y sin GSE134073 (comparador
#           HPB) y compara lado a lado la firma global y el comportamiento del
#           eje AR, para decidir si se mantiene o se excluye ese dataset.
# Entra:    Datasets_Remediados/*_remediado.rds
# Sale:     Meta_Analisis_Resultados/prostata_*_comparativo.csv,
#           Meta_Analisis_Objetos/prostata_sensibilidad.rds
# ______________________________________________________________________________

# __ INDEX ___
# 01. Setup, librerias y parametros
# 02. Rutas y carga de los datasets remediados prostaticos
# 03. Verificacion del comparador en todos los datasets
# 04. Construccion de los dos objectMA prostaticos
# 05. Funciones de meta-analisis reutilizadas
# 06. Ejecucion del meta-analisis en ambas versiones
# 07. Comportamiento del eje AR en ambas versiones
# 08. Tabla comparativa lado a lado
# 09. Exportacion y resumen final


# 01. SETUP ____________________________________________________________________

library(tidyverse)
library(DExMA)
library(R.utils)

select    <- dplyr::select
filter    <- dplyr::filter
count     <- dplyr::count
slice     <- dplyr::slice
arrange   <- dplyr::arrange
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
group_by  <- dplyr::group_by
map       <- purrr::map
map_chr   <- purrr::map_chr
map_dbl   <- purrr::map_dbl
map_dfr   <- purrr::map_dfr
imap      <- purrr::imap
imap_dfr  <- purrr::imap_dfr
pmap      <- purrr::pmap
pmap_dfr  <- purrr::pmap_dfr
keep      <- purrr::keep
reduce    <- purrr::reduce

parametros_lista <- list(
  proporcion_minima    = 2 / 3,
  k_vecinos_imputacion = 7,
  umbral_fdr           = 0.05,
  timeout_imputacion   = 600,
  dataset_a_excluir    = "GSE134073",
  semilla_aleatoria    = 1234
)

set.seed(parametros_lista$semilla_aleatoria)

# Genes del eje del receptor de androgenos, para rastrear su comportamiento.
genes_eje_ar_chr <- c("AR", "KLK3", "KLK2", "NKX3-1", "FKBP5",
                      "TMPRSS2", "FOLH1", "STEAP1", "SLC45A3", "HOXB13")


# 02. RUTAS Y CARGA DE LOS DATASETS REMEDIADOS PROSTATICOS ______________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_remediados_chr <- file.path(ruta_base_chr, "Datasets_Remediados")
ruta_salida_chr     <- file.path(ruta_base_chr, "Meta_Analisis_Resultados")
ruta_objetos_chr    <- file.path(ruta_base_chr, "Meta_Analisis_Objetos")

archivos_remediados_chr <- list.files(ruta_remediados_chr,
                                      pattern     = "remediado",
                                      full.names  = TRUE,
                                      ignore.case = TRUE)
archivos_remediados_chr <- archivos_remediados_chr[
  grepl("\\.rds$", archivos_remediados_chr, ignore.case = TRUE)
]

cargar_remediado <- function(ruta_archivo_chr) {
  remediado_lista <- readRDS(ruta_archivo_chr)
  list(
    gse_id             = remediado_lista$dataset_descriptor$gse_id,
    organo             = remediado_lista$dataset_descriptor$organo,
    comparador         = remediado_lista$dataset_descriptor$comparator_subtype,
    expresion_matriz   = remediado_lista$expresion_matriz,
    grupo_binario      = remediado_lista$grupo_binario
  )
}

datasets_lista <- archivos_remediados_chr %>%
  map(cargar_remediado)

names(datasets_lista) <- map_chr(datasets_lista, ~ .x$gse_id)


# 03. VERIFICACION DEL COMPARADOR EN TODOS LOS DATASETS ________________________

# Antes de decidir sobre GSE134073, se verifica que ningun otro dataset use un
# comparador distinto de tejido normal. Es el control de homogeneidad del
# contraste que la auditoria inicial no incluyo, y debe quedar documentado.

# Normaliza el comparador a un unico string, cubriendo los casos en que el campo
# no existe (NULL) o esta vacio.
normalizar_comparador <- function(valor) {
  if (is.null(valor) || length(valor) == 0) return("No especificado")
  if (is.na(valor[1]))                       return("No especificado")
  as.character(valor[1])
}

comparadores_tbl <- datasets_lista %>%
  imap_dfr(function(dataset_lista, gse_id_chr) {
    tibble(
      gse_id     = gse_id_chr,
      organo     = dataset_lista$organo,
      comparador = normalizar_comparador(dataset_lista$comparador)
    )
  }) %>%
  arrange(organo, gse_id)

write_csv(comparadores_tbl,
          file.path(ruta_salida_chr, "verificacion_comparadores.csv"))

message("--- Comparadores por dataset ---")
print(comparadores_tbl, n = Inf)


# 04. CONSTRUCCION DE LOS DOS objectMA PROSTATICOS _____________________________

# Se reutiliza el patron de construccion validado en C1: phenodata con columna
# de condicion y argumentos de grupo como vectores explicitos.

construir_objectma <- function(datasets_grupo_lista) {

  matrices_lista <- datasets_grupo_lista %>%
    map(~ .x$expresion_matriz)

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

  createObjectMA(
    listEX    = matrices_lista,
    listPheno = phenodatas_lista,
    namePheno = rep("condicion", n_estudios_int),
    expGroups = rep("Tumor",     n_estudios_int),
    refGroups = rep("Normal",    n_estudios_int)
  )
}

# Datasets prostaticos, con y sin el dataset a excluir.
datasets_prostata_lista <- datasets_lista %>%
  keep(~ .x$organo == "Prostate")

datasets_6ds_lista <- datasets_prostata_lista
datasets_5ds_lista <- datasets_prostata_lista %>%
  keep(~ .x$gse_id != parametros_lista$dataset_a_excluir)

message("\nVersion 6ds: ",
        paste(map_chr(datasets_6ds_lista, ~ .x$gse_id), collapse = ", "))
message("Version 5ds: ",
        paste(map_chr(datasets_5ds_lista, ~ .x$gse_id), collapse = ", "))

message("\nConstruyendo objectMA (6 datasets)")
objeto_6ds <- construir_objectma(datasets_6ds_lista)
objeto_6ds <- dataLog(objeto_6ds)

message("Construyendo objectMA (5 datasets)")
objeto_5ds <- construir_objectma(datasets_5ds_lista)
objeto_5ds <- dataLog(objeto_5ds)


# 05. FUNCIONES DE META-ANALISIS REUTILIZADAS __________________________________

# Prefiltro por proporcion antes de imputar, para hacer viable el KNN.

prefiltrar_por_proporcion <- function(objeto_ma, proporcion_num) {

  n_estudios_int <- length(objeto_ma)

  conteo_tbl <- objeto_ma %>%
    map(~ rownames(.x$mExpres)) %>%
    unlist() %>%
    tibble(gen = .) %>%
    count(gen, name = "n_estudios")

  umbral_int <- ceiling(proporcion_num * n_estudios_int)
  genes_chr  <- conteo_tbl %>%
    filter(n_estudios >= umbral_int) %>%
    pull(gen)

  objeto_ma %>%
    map(function(estudio_lista) {
      presentes_chr <- intersect(rownames(estudio_lista$mExpres), genes_chr)
      estudio_lista$mExpres <- estudio_lista$mExpres[presentes_chr, ,
                                                     drop = FALSE]
      estudio_lista
    })
}

# Imputacion protegida con limite de tiempo.

imputar_con_timeout <- function(objeto_ma) {

  objeto_prefiltrado <- prefiltrar_por_proporcion(
    objeto_ma, parametros_lista$proporcion_minima
  )

  tryCatch(
    withTimeout(
      {
        missGenesImput(objeto_prefiltrado,
                       k = parametros_lista$k_vecinos_imputacion)$objectMA
      },
      timeout   = parametros_lista$timeout_imputacion,
      onTimeout = "error"
    ),
    TimeoutException = function(e) { message("    imputacion abortada"); NULL },
    error            = function(e) { message("    imputacion fallida"); NULL }
  )
}

# Ejecuta los tres brazos por tres metodos para un objectMA dado.

ejecutar_meta_version <- function(objeto_ma, etiqueta_chr) {

  message("\n[", etiqueta_chr, "] Ejecutando meta-analisis")

  message("  Resolviendo imputacion")
  objeto_imputado <- imputar_con_timeout(objeto_ma)
  imputacion_viable_lgl <- !is.null(objeto_imputado)

  rejilla_tbl <- expand_grid(
    brazo  = c("comun", "proporcion", "imputacion"),
    metodo = c("Stouffer", "Fisher", "maxP")
  )

  rejilla_tbl %>%
    pmap_dfr(function(brazo, metodo) {

      if (brazo == "imputacion") {
        if (!imputacion_viable_lgl) return(NULL)
        objeto_trabajo <- objeto_imputado
        proporcion_num <- 1
      } else if (brazo == "proporcion") {
        objeto_trabajo <- objeto_ma
        proporcion_num <- parametros_lista$proporcion_minima
      } else {
        objeto_trabajo <- objeto_ma
        proporcion_num <- 1
      }

      message("    ", brazo, " | ", metodo)

      resultado_tbl <- metaAnalysisDE(objeto_trabajo,
                                      typeMethod     = metodo,
                                      proportionData = proporcion_num) %>%
        as.data.frame() %>%
        rownames_to_column("gen_simbolo") %>%
        as_tibble()

      resultado_tbl %>%
        mutate(version = etiqueta_chr, brazo = brazo, metodo = metodo,
               .before = 1)
    })
}


# 06. EJECUCION DEL META-ANALISIS EN AMBAS VERSIONES ___________________________

message("\n=== META-ANALISIS: VERSION 6 DATASETS ===")
resultados_6ds_tbl <- ejecutar_meta_version(objeto_6ds, "6ds")

message("\n=== META-ANALISIS: VERSION 5 DATASETS ===")
resultados_5ds_tbl <- ejecutar_meta_version(objeto_5ds, "5ds")

resultados_ambos_tbl <- bind_rows(resultados_6ds_tbl, resultados_5ds_tbl)

saveRDS(resultados_ambos_tbl,
        file.path(ruta_objetos_chr, "prostata_sensibilidad.rds"))


# 07. COMPORTAMIENTO DEL EJE AR EN AMBAS VERSIONES _____________________________

# Se extrae, por version, el logFC y signo por estudio de los genes del eje AR,
# usando calculateES. La comparacion muestra si excluir GSE134073 devuelve a AR
# un comportamiento concordante.

eje_ar_por_version <- function(objeto_ma, etiqueta_chr) {

  es_lista     <- calculateES(objeto_ma)
  logfc_matriz <- es_lista$logFC
  var_matriz   <- es_lista$Var
  es_matriz    <- es_lista$ES

  # logFC combinado ponderado por inverso de varianza.
  pesos_matriz <- 1 / var_matriz
  logfc_comb_num <- rowSums(logfc_matriz * pesos_matriz, na.rm = TRUE) /
    rowSums(pesos_matriz, na.rm = TRUE)
  es_comb_num <- rowSums(es_matriz * pesos_matriz, na.rm = TRUE) /
    rowSums(pesos_matriz, na.rm = TRUE)

  genes_presentes_chr <- intersect(genes_eje_ar_chr, rownames(logfc_matriz))

  tibble(
    version            = etiqueta_chr,
    gen                = genes_presentes_chr,
    logfc_combinado    = round(logfc_comb_num[genes_presentes_chr], 3),
    direccion          = if_else(es_comb_num[genes_presentes_chr] > 0,
                                 "Up", "Down"),
    n_estudios_signo_pos = map_dbl(genes_presentes_chr, function(g) {
      sum(sign(logfc_matriz[g, ]) > 0, na.rm = TRUE)
    }),
    n_estudios_signo_neg = map_dbl(genes_presentes_chr, function(g) {
      sum(sign(logfc_matriz[g, ]) < 0, na.rm = TRUE)
    })
  )
}

eje_ar_6ds_tbl <- eje_ar_por_version(objeto_6ds, "6ds")
eje_ar_5ds_tbl <- eje_ar_por_version(objeto_5ds, "5ds")

eje_ar_comparativo_tbl <- bind_rows(eje_ar_6ds_tbl, eje_ar_5ds_tbl) %>%
  arrange(gen, version)

write_csv(eje_ar_comparativo_tbl,
          file.path(ruta_salida_chr, "prostata_eje_ar_comparativo.csv"))


# 08. TABLA COMPARATIVA LADO A LADO ____________________________________________

# Conteo de genes significativos por version, brazo y metodo, para ver el efecto
# global de excluir GSE134073.

detectar_columna_fdr <- function(nombres_chr) {
  candidatos_chr <- c("FDR", "padj", "adj.P.Val", "p.adjust", "adjPval")
  interseccion_chr <- intersect(candidatos_chr, nombres_chr)
  if (length(interseccion_chr) == 0) return(NA_character_)
  interseccion_chr[1]
}

columna_fdr_chr <- detectar_columna_fdr(names(resultados_ambos_tbl))

resumen_comparativo_tbl <- resultados_ambos_tbl %>%
  mutate(fdr_valor = .data[[columna_fdr_chr]]) %>%
  group_by(version, brazo, metodo) %>%
  summarise(
    n_genes          = n(),
    n_significativos = sum(fdr_valor < parametros_lista$umbral_fdr,
                           na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = version,
    values_from = c(n_genes, n_significativos)
  )

write_csv(resumen_comparativo_tbl,
          file.path(ruta_salida_chr, "prostata_resumen_comparativo.csv"))


# 09. EXPORTACION Y RESUMEN FINAL ______________________________________________

message("\n________________ SENSIBILIDAD PROSTATA 6ds vs 5ds ________________\n")

message("--- Genes significativos por brazo y metodo ---")
print(resumen_comparativo_tbl, n = Inf)

message("\n--- Eje AR: comparacion de direccion ---")

eje_ar_comparativo_tbl %>%
  select(gen, version, logfc_combinado, direccion,
         n_estudios_signo_pos, n_estudios_signo_neg) %>%
  pmap(function(gen, version, logfc_combinado, direccion,
                n_estudios_signo_pos, n_estudios_signo_neg) {
    message("  ", gen, " [", version, "]: logFC ", logfc_combinado,
            " (", direccion, ") | +", n_estudios_signo_pos,
            " / -", n_estudios_signo_neg, " estudios")
  }) %>%
  invisible()

# Comprobacion clave: cambia la direccion de AR entre versiones?
ar_6ds <- eje_ar_comparativo_tbl %>%
  filter(gen == "AR", version == "6ds") %>% pull(direccion)
ar_5ds <- eje_ar_comparativo_tbl %>%
  filter(gen == "AR", version == "5ds") %>% pull(direccion)

message("\n--- Veredicto sobre AR ---")

if (length(ar_6ds) > 0 && length(ar_5ds) > 0) {
  if (ar_6ds != ar_5ds) {
    message("  AR CAMBIA de direccion: ", ar_6ds, " (6ds) -> ", ar_5ds,
            " (5ds)")
    message("  Confirma que GSE134073 (HPB) invertia la senal de AR.")
    message("  Recomendacion: usar la version de 5 datasets como principal.")
  } else {
    message("  AR mantiene la direccion (", ar_5ds, ") en ambas versiones.")
    message("  Excluir GSE134073 no altera la conclusion sobre AR.")
  }
} else {
  message("  AR no esta presente en alguna de las versiones; revisar.")
}

message("\nTablas exportadas en: ", ruta_salida_chr)
message("___________________________________________________________________\n")
