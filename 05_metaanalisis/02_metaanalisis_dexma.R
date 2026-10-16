# ______________________________________________________________________________
# 05_02  Meta-analisis con DExMA (tres brazos de genes faltantes x tres metodos)
#
# Que hace: ejecuta el meta-analisis de cada grupo reutilizando los objectMA ya
#           construidos, con tres estrategias de genes faltantes y tres metodos
#           de combinacion, y recupera la direccionalidad por tamanos de efecto.
# Entra:    Meta_Analisis_Objetos/objetos_ma_por_grupo.rds
# Sale:     resultados_meta_completos (.rds/.csv), direccionalidad_por_grupo.csv,
#           resumen_significativos.csv
# ______________________________________________________________________________

# __ INDEX ___
# 01. Setup, librerias y parametros de la etapa
# 02. Rutas y carga de los objectMA ya construidos
# 03. Definicion de la rejilla de meta-analisis
# 04. Prefiltro del espacio genico por proporcion
# 05. Imputacion protegida con limite de tiempo
# 06. Ejecucion de una corrida individual
# 07. Ejecucion del meta-analisis completo
# 08. Recuperacion de la direccionalidad
# 09. Exportacion de resultados
# 10. Resumen final


# 01. SETUP ____________________________________________________________________

library(tidyverse)
library(DExMA)
library(R.utils)   # para withTimeout()

# Preferencia de funciones frente al enmascaramiento de Bioconductor.
select    <- dplyr::select
filter    <- dplyr::filter
count     <- dplyr::count
slice     <- dplyr::slice
rename    <- dplyr::rename
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
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

parametros_meta_lista <- list(
  proporcion_minima    = 2 / 3,
  k_vecinos_imputacion = 7,
  umbral_fdr           = 0.05,

  # Limite de tiempo, en segundos, para la imputacion de un grupo. Si se supera,
  # el brazo de imputacion se marca como no completado para ese grupo.
  timeout_imputacion   = 600,

  semilla_aleatoria    = 1234
)

set.seed(parametros_meta_lista$semilla_aleatoria)


# 02. RUTAS Y CARGA DE LOS objectMA YA CONSTRUIDOS _____________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_salida_chr  <- file.path(ruta_base_chr, "Meta_Analisis_Resultados")
ruta_objetos_chr <- file.path(ruta_base_chr, "Meta_Analisis_Objetos")

archivo_objetos_chr <- file.path(ruta_objetos_chr, "objetos_ma_por_grupo.rds")

if (!file.exists(archivo_objetos_chr)) {
  stop("No se encontro objetos_ma_por_grupo.rds. Ejecutar antes la construccion.")
}

objetos_ma_lista <- readRDS(archivo_objetos_chr)

message("objectMA cargados: ", paste(names(objetos_ma_lista), collapse = ", "))

# Se verifica que esten en escala logaritmica. Como ya se remediaron, dataLog
# debe confirmarlo sin cambios relevantes.
objetos_ma_lista <- objetos_ma_lista %>%
  imap(function(objeto_ma, organo_chr) {
    message("[", organo_chr, "] Verificando escala logaritmica")
    dataLog(objeto_ma)
  })


# 03. DEFINICION DE LA REJILLA DE META-ANALISIS ________________________________

rejilla_meta_tbl <- expand_grid(
  brazo_faltantes = c("comun", "proporcion", "imputacion"),
  metodo          = c("Stouffer", "Fisher", "maxP")
)

message("\nCombinaciones por grupo: ", nrow(rejilla_meta_tbl))


# 04. PREFILTRO DEL ESPACIO GENICO POR PROPORCION ______________________________

# La imputacion sampleKNN sobre la union total de genes es inviable cuando el
# espacio compartido es pequeno (en prostata la interseccion es solo el 13.7%).
# Se restringe el objeto a los genes presentes en al menos la proporcion minima
# de estudios ANTES de imputar: no tiene sentido imputar genes medidos en un
# unico dataset, y esto reduce drasticamente el coste del KNN.

prefiltrar_por_proporcion <- function(objeto_ma, proporcion_num) {

  n_estudios_int <- length(objeto_ma)

  # Recuento, por gen, del numero de estudios en que aparece.
  genes_por_estudio_lista <- objeto_ma %>%
    map(~ rownames(.x$mExpres))

  conteo_genes_tbl <- tibble(gen = unlist(genes_por_estudio_lista)) %>%
    count(gen, name = "n_estudios")

  umbral_estudios_int <- ceiling(proporcion_num * n_estudios_int)

  genes_conservar_chr <- conteo_genes_tbl %>%
    filter(n_estudios >= umbral_estudios_int) %>%
    pull(gen)

  message("    Prefiltro: ", length(genes_conservar_chr),
          " genes en >= ", umbral_estudios_int, "/", n_estudios_int, " estudios")

  # Se recorta cada matriz del objeto a esos genes, respetando el orden propio.
  objeto_filtrado <- objeto_ma %>%
    map(function(estudio_lista) {
      genes_presentes_chr <- intersect(rownames(estudio_lista$mExpres),
                                       genes_conservar_chr)
      estudio_lista$mExpres <- estudio_lista$mExpres[genes_presentes_chr, ,
                                                     drop = FALSE]
      estudio_lista
    })

  objeto_filtrado
}


# 05. IMPUTACION PROTEGIDA CON LIMITE DE TIEMPO ________________________________

# Ejecuta missGenesImput sobre el objeto prefiltrado, abortando si supera el
# limite de tiempo. Devuelve el objeto imputado o NULL si no fue viable.

imputar_con_timeout <- function(objeto_ma, organo_chr) {

  objeto_prefiltrado <- prefiltrar_por_proporcion(
    objeto_ma,
    parametros_meta_lista$proporcion_minima
  )

  message("    Imputando (limite ",
          parametros_meta_lista$timeout_imputacion, "s)...")

  resultado <- tryCatch(
    withTimeout(
      {
        imputacion_lista <- missGenesImput(
          objeto_prefiltrado,
          k = parametros_meta_lista$k_vecinos_imputacion
        )
        imputacion_lista$objectMA
      },
      timeout   = parametros_meta_lista$timeout_imputacion,
      onTimeout = "error"
    ),
    TimeoutException = function(e) {
      message("    IMPUTACION ABORTADA: excedio el limite de tiempo")
      NULL
    },
    error = function(e) {
      message("    IMPUTACION FALLIDA: ", conditionMessage(e))
      NULL
    }
  )

  resultado
}


# 06. EJECUCION DE UNA CORRIDA INDIVIDUAL ______________________________________

# Devuelve la tabla de resultados de una corrida, o NULL si el brazo no pudo
# completarse (unico caso posible: imputacion inviable).

ejecutar_corrida <- function(objeto_ma, brazo_chr, metodo_chr, organo_chr) {

  objeto_trabajo <- objeto_ma
  proporcion_num <- 1

  if (brazo_chr == "proporcion") {
    proporcion_num <- parametros_meta_lista$proporcion_minima
  }

  if (brazo_chr == "imputacion") {
    objeto_trabajo <- imputar_con_timeout(objeto_ma, organo_chr)
    if (is.null(objeto_trabajo)) {
      return(NULL)
    }
    proporcion_num <- 1
  }

  resultado_tbl <- metaAnalysisDE(
    objeto_trabajo,
    typeMethod     = metodo_chr,
    proportionData = proporcion_num
  )

  resultado_tbl %>%
    as.data.frame() %>%
    rownames_to_column("gen_simbolo") %>%
    as_tibble()
}


# 07. EJECUCION DEL META-ANALISIS COMPLETO _____________________________________

# La imputacion se resuelve una sola vez por grupo y se reutiliza para los tres
# metodos, evitando repetir el KNN tres veces. Los brazos comun y proporcion no
# necesitan pre-calculo.

ejecutar_grupo <- function(objeto_ma, organo_chr) {

  message("\n[", organo_chr, "] Ejecutando meta-analisis")

  # Se pre-resuelve la imputacion una vez para el grupo.
  message("  Resolviendo brazo de imputacion para el grupo")
  objeto_imputado <- imputar_con_timeout(objeto_ma, organo_chr)
  imputacion_viable_lgl <- !is.null(objeto_imputado)

  if (!imputacion_viable_lgl) {
    message("  [", organo_chr, "] imputacion no viable; se omite ese brazo")
  }

  rejilla_meta_tbl %>%
    pmap_dfr(function(brazo_faltantes, metodo) {

      # El brazo de imputacion reutiliza el objeto ya imputado.
      if (brazo_faltantes == "imputacion") {
        if (!imputacion_viable_lgl) return(NULL)

        message("    imputacion | ", metodo)
        resultado_tbl <- metaAnalysisDE(
          objeto_imputado,
          typeMethod     = metodo,
          proportionData = 1
        ) %>%
          as.data.frame() %>%
          rownames_to_column("gen_simbolo") %>%
          as_tibble()

      } else {
        message("    ", brazo_faltantes, " | ", metodo)
        resultado_tbl <- ejecutar_corrida(objeto_ma, brazo_faltantes,
                                          metodo, organo_chr)
      }

      if (is.null(resultado_tbl)) return(NULL)

      resultado_tbl %>%
        mutate(
          organo          = organo_chr,
          brazo_faltantes = brazo_faltantes,
          metodo          = metodo,
          .before = 1
        )
    })
}

message("\n=== EJECUCION DEL META-ANALISIS ===")

resultados_meta_tbl <- objetos_ma_lista %>%
  imap_dfr(ejecutar_grupo)


# 08. RECUPERACION DE LA DIRECCIONALIDAD _______________________________________

# calculateES devuelve, por gen y estudio, el tamano de efecto (g de Hedges,
# $ES), su varianza ($Var) y el log fold-change ($logFC). El signo del efecto
# combinado ponderado da la direccion de consenso.

calcular_direccion_grupo <- function(objeto_ma, organo_chr) {

  message("[", organo_chr, "] Calculando direccionalidad")

  es_lista <- calculateES(objeto_ma)

  efectos_matriz   <- es_lista$ES
  varianzas_matriz <- es_lista$Var
  logfc_matriz     <- es_lista$logFC

  pesos_matriz <- 1 / varianzas_matriz

  efecto_combinado_num <- rowSums(efectos_matriz * pesos_matriz, na.rm = TRUE) /
    rowSums(pesos_matriz, na.rm = TRUE)

  logfc_combinado_num <- rowSums(logfc_matriz * pesos_matriz, na.rm = TRUE) /
    rowSums(pesos_matriz, na.rm = TRUE)

  signo_consenso_int <- sign(efecto_combinado_num)

  concordancia_num <- map_dbl(seq_len(nrow(efectos_matriz)), function(i) {
    signos_estudio_int <- sign(efectos_matriz[i, ])
    signos_validos_int <- signos_estudio_int[!is.na(signos_estudio_int)]
    if (length(signos_validos_int) == 0) return(NA_real_)
    mean(signos_validos_int == signo_consenso_int[i])
  })

  tibble(
    organo             = organo_chr,
    gen_simbolo        = rownames(efectos_matriz),
    efecto_combinado   = round(efecto_combinado_num, 4),
    logfc_combinado    = round(logfc_combinado_num, 4),
    direccion          = if_else(signo_consenso_int > 0, "Up", "Down"),
    concordancia_signo = round(concordancia_num, 3)
  )
}

direccionalidad_tbl <- objetos_ma_lista %>%
  imap_dfr(calcular_direccion_grupo)

write_csv(direccionalidad_tbl,
          file.path(ruta_salida_chr, "direccionalidad_por_grupo.csv"))


# 09. EXPORTACION DE RESULTADOS ________________________________________________

resultados_completos_tbl <- resultados_meta_tbl %>%
  left_join(direccionalidad_tbl, by = c("organo", "gen_simbolo"))

saveRDS(resultados_completos_tbl,
        file.path(ruta_objetos_chr, "resultados_meta_completos.rds"))

write_csv(resultados_completos_tbl,
          file.path(ruta_salida_chr, "resultados_meta_completos.csv"))


# 10. RESUMEN FINAL ____________________________________________________________

detectar_columna_fdr <- function(nombres_chr) {
  candidatos_chr   <- c("FDR", "padj", "adj.P.Val", "p.adjust", "adjPval")
  interseccion_chr <- intersect(candidatos_chr, nombres_chr)
  if (length(interseccion_chr) == 0) return(NA_character_)
  interseccion_chr[1]
}

columna_fdr_chr <- detectar_columna_fdr(names(resultados_completos_tbl))

message("\n________________ RESUMEN DEL META-ANALISIS ________________\n")

message("Columna de FDR detectada: ",
        if_else(is.na(columna_fdr_chr), "NINGUNA (revisar nombres)",
                columna_fdr_chr))

if (!is.na(columna_fdr_chr)) {

  resumen_significativos_tbl <- resultados_completos_tbl %>%
    mutate(fdr_valor = .data[[columna_fdr_chr]]) %>%
    group_by(organo, brazo_faltantes, metodo) %>%
    summarise(
      n_genes          = n(),
      n_significativos = sum(fdr_valor < parametros_meta_lista$umbral_fdr,
                             na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(pct_significativos = round(100 * n_significativos / n_genes, 1))

  write_csv(resumen_significativos_tbl,
            file.path(ruta_salida_chr, "resumen_significativos.csv"))

  resumen_significativos_tbl %>%
    pmap(function(organo, brazo_faltantes, metodo, n_genes,
                  n_significativos, pct_significativos) {
      message("  [", organo, "] ", brazo_faltantes, " | ", metodo,
              " -> ", n_significativos, "/", n_genes,
              " (", pct_significativos, "%)")
    }) %>%
    invisible()

} else {
  message("Columnas disponibles: ",
          paste(names(resultados_completos_tbl), collapse = ", "))
}

# Se registra que brazos de imputacion se completaron, para documentar en la
# metodologia los grupos donde la imputacion no fue viable.
brazos_completados_tbl <- resultados_meta_tbl %>%
  distinct(organo, brazo_faltantes) %>%
  arrange(organo, brazo_faltantes)

write_csv(brazos_completados_tbl,
          file.path(ruta_salida_chr, "brazos_completados.csv"))

message("\nBrazos completados por grupo:")

brazos_completados_tbl %>%
  group_by(organo) %>%
  summarise(brazos = paste(brazo_faltantes, collapse = ", "),
            .groups = "drop") %>%
  pmap(function(organo, brazos) {
    message("  ", organo, ": ", brazos)
  }) %>%
  invisible()

message("\nObjetos guardados en: ", ruta_objetos_chr)
message("Tablas exportadas en: ", ruta_salida_chr)
message("___________________________________________________________\n")
