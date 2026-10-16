# ______________________________________________________________________________
# 05_03  Meta-analisis Stouffer con signo, ponderaciones y robustez
#
# Que hace: combina los z firmados de cada estudio (Stouffer ponderado), aplica
#           el filtro de concordancia de signo y evalua robustez con tres
#           ponderaciones y leave-one-dataset-out.
# Entra:    Meta_Analisis_Objetos/objetos_ma_por_grupo.rds
# Sale:     meta_firmado_completo (.csv/.rds), resumen_stouffer_firmado.csv,
#           leave_one_out.csv
# ______________________________________________________________________________

# __ INDEX ___
# 01. Setup, librerias y parametros de la etapa
# 02. Rutas y carga de los objectMA
# 03. Extraccion de estadisticos por estudio (ES, Var, logFC, n)
# 04. Derivacion de z-scores firmados por estudio
# 05. Combinacion de Stouffer con signo y ponderacion
# 06. Filtro de concordancia de signo (Marot)
# 07. Meta-analisis firmado por grupo y esquema de ponderacion
# 08. Analisis leave-one-dataset-out
# 09. Exportacion de resultados
# 10. Resumen final


# 01. SETUP ____________________________________________________________________

library(tidyverse)
library(DExMA)

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

parametros_c2_lista <- list(
  umbral_fdr              = 0.05,

  # Umbrales del filtro de concordancia de signo. Se reportan ambos: el estricto
  # (todos los estudios en la misma direccion) y el mayoritario.
  concordancia_estricta   = 1.00,
  concordancia_mayoritaria = 0.75,

  semilla_aleatoria       = 1234
)

set.seed(parametros_c2_lista$semilla_aleatoria)


# 02. RUTAS Y CARGA DE LOS objectMA ____________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_salida_chr  <- file.path(ruta_base_chr, "Meta_Analisis_Resultados")
ruta_objetos_chr <- file.path(ruta_base_chr, "Meta_Analisis_Objetos")

objetos_ma_lista <- readRDS(
  file.path(ruta_objetos_chr, "objetos_ma_por_grupo.rds")
)

message("objectMA cargados: ", paste(names(objetos_ma_lista), collapse = ", "))


# 03. EXTRACCION DE ESTADISTICOS POR ESTUDIO ___________________________________

# calculateES devuelve, por gen y estudio: tamano de efecto (ES), su varianza
# (Var) y el log fold-change (logFC). El tamano muestral de cada estudio se
# cuenta directamente del vector de condicion del objectMA, garantizando que el
# n usado para ponderar es el mismo del analisis.

extraer_estadisticos_grupo <- function(objeto_ma) {

  es_lista <- calculateES(objeto_ma)

  # Tamanos muestrales por estudio, en el mismo orden que las columnas de las
  # matrices de efecto.
  tamanos_tbl <- objeto_ma %>%
    imap_dfr(function(estudio_lista, gse_id_chr) {
      condicion_int <- estudio_lista$condition
      n_tumor_int   <- sum(condicion_int == 1, na.rm = TRUE)
      n_normal_int  <- sum(condicion_int == 0, na.rm = TRUE)
      tibble(
        gse_id     = gse_id_chr,
        n_tumor    = n_tumor_int,
        n_normal   = n_normal_int,
        n_total    = n_tumor_int + n_normal_int,
        # n efectivo: media armonica de los tamanos de grupo, que gobierna la
        # varianza de un contraste de dos grupos.
        n_efectivo = (n_tumor_int * n_normal_int) / (n_tumor_int + n_normal_int)
      )
    })

  list(
    ES      = es_lista$ES,
    Var     = es_lista$Var,
    logFC   = es_lista$logFC,
    tamanos = tamanos_tbl
  )
}


# 04. DERIVACION DE Z-SCORES FIRMADOS POR ESTUDIO ______________________________

# Para cada gen y estudio se calcula un z-score firmado. El p-valor de dos colas
# se deriva del estadistico ES/sqrt(Var) (lo mismo que hace DExMA internamente),
# y el signo se toma del logFC. El resultado es un z cuya magnitud refleja la
# evidencia y cuyo signo refleja la direccion del cambio.

derivar_z_firmado <- function(estadisticos_lista) {

  efectos_matriz   <- estadisticos_lista$ES
  varianzas_matriz <- estadisticos_lista$Var
  logfc_matriz     <- estadisticos_lista$logFC

  # Estadistico de Wald por gen y estudio.
  wald_matriz <- efectos_matriz / sqrt(varianzas_matriz)

  # p-valor de dos colas correspondiente.
  p_matriz <- 2 * pnorm(-abs(wald_matriz))

  # z firmado: la magnitud viene del p de dos colas convertido a cuantil normal
  # de una cola, y el signo viene del logFC.
  z_una_cola_matriz <- qnorm(1 - p_matriz / 2)
  z_firmado_matriz  <- sign(logfc_matriz) * z_una_cola_matriz

  list(
    z_firmado = z_firmado_matriz,
    p_valor   = p_matriz,
    signo     = sign(logfc_matriz)
  )
}


# 05. COMBINACION DE STOUFFER CON SIGNO Y PONDERACION __________________________

# Combina los z firmados de todos los estudios en un unico Z por gen, con la
# formula de Stouffer ponderada:
#   Z = sum(w_i * z_i) / sqrt(sum(w_i^2))
# Los pesos dependen del esquema elegido. El signo del Z combinado es la
# direccion meta-analitica; un gen con direcciones opuestas entre estudios se
# cancela en lugar de dar un falso positivo.

combinar_stouffer_firmado <- function(z_firmado_matriz, pesos_num) {

  n_genes_int <- nrow(z_firmado_matriz)

  z_combinado_num <- map_dbl(seq_len(n_genes_int), function(i) {

    z_num <- z_firmado_matriz[i, ]

    validos_lgl <- !is.na(z_num) & is.finite(z_num)
    z_num       <- z_num[validos_lgl]
    w_num       <- pesos_num[validos_lgl]

    if (length(z_num) == 0) return(NA_real_)

    sum(w_num * z_num) / sqrt(sum(w_num^2))
  })

  # p-valor de dos colas del Z combinado.
  p_combinado_num <- 2 * pnorm(-abs(z_combinado_num))

  tibble(
    gen_simbolo     = rownames(z_firmado_matriz),
    z_combinado     = round(z_combinado_num, 4),
    p_combinado     = p_combinado_num,
    direccion       = if_else(z_combinado_num > 0, "Up", "Down")
  )
}


# Devuelve el vector de pesos segun el esquema, alineado al orden de estudios.

calcular_pesos <- function(tamanos_tbl, esquema_chr) {

  if (esquema_chr == "n_efectivo") {
    return(sqrt(tamanos_tbl$n_efectivo))
  }
  if (esquema_chr == "n_total") {
    return(sqrt(tamanos_tbl$n_total))
  }
  # sin ponderar: todos los estudios pesan igual.
  rep(1, nrow(tamanos_tbl))
}


# 06. FILTRO DE CONCORDANCIA DE SIGNO __________________________________________

# Calcula, por gen, la proporcion de estudios cuyo signo coincide con la
# direccion del Z combinado. Es el filtro de Marot: descarta genes cuya
# significancia combinada proviene de estudios que discrepan en direccion.

calcular_concordancia <- function(signo_matriz, direccion_chr) {

  direccion_int <- if_else(direccion_chr == "Up", 1, -1)

  map_dbl(seq_len(nrow(signo_matriz)), function(i) {
    signos_num <- signo_matriz[i, ]
    signos_num <- signos_num[!is.na(signos_num)]
    if (length(signos_num) == 0) return(NA_real_)
    mean(signos_num == direccion_int[i])
  })
}


# 07. META-ANALISIS FIRMADO POR GRUPO Y ESQUEMA ________________________________

# Ejecuta el Stouffer firmado para un grupo bajo los tres esquemas de
# ponderacion, anadiendo la concordancia de signo y las banderas de los dos
# umbrales de filtro.

meta_firmado_grupo <- function(objeto_ma, organo_chr) {

  message("\n[", organo_chr, "] Stouffer firmado, tres ponderaciones")

  estadisticos_lista <- extraer_estadisticos_grupo(objeto_ma)
  z_lista            <- derivar_z_firmado(estadisticos_lista)

  esquemas_chr <- c("n_efectivo", "n_total", "sin_peso")

  esquemas_chr %>%
    map_dfr(function(esquema_chr) {

      message("    Esquema: ", esquema_chr)

      pesos_num <- calcular_pesos(estadisticos_lista$tamanos, esquema_chr)

      combinado_tbl <- combinar_stouffer_firmado(z_lista$z_firmado, pesos_num)

      concordancia_num <- calcular_concordancia(z_lista$signo,
                                                combinado_tbl$direccion)

      combinado_tbl %>%
        mutate(
          organo               = organo_chr,
          esquema_ponderacion  = esquema_chr,
          fdr                  = p.adjust(p_combinado, method = "BH"),
          concordancia_signo   = round(concordancia_num, 3),
          pasa_concordancia_estricta =
            concordancia_num >= parametros_c2_lista$concordancia_estricta,
          pasa_concordancia_mayoritaria =
            concordancia_num >= parametros_c2_lista$concordancia_mayoritaria,
          .before = 1
        )
    })
}

message("\n=== STOUFFER FIRMADO POR GRUPO ===")

meta_firmado_tbl <- objetos_ma_lista %>%
  imap_dfr(meta_firmado_grupo)

write_csv(meta_firmado_tbl,
          file.path(ruta_salida_chr, "meta_firmado_completo.csv"))


# 08. ANALISIS LEAVE-ONE-DATASET-OUT ___________________________________________

# Recalcula el meta-analisis firmado de cada grupo quitando un dataset por vez,
# usando el esquema de n efectivo. Se compara el numero de genes significativos
# con el analisis completo. Una caida grande al quitar un dataset indica que ese
# dataset domina el resultado del grupo, situacion esperada para TCGA en el
# grupo testicular dado su n efectivo muy superior al de los otros estudios.

leave_one_out_grupo <- function(objeto_ma, organo_chr) {

  gse_ids_chr <- names(objeto_ma)

  # Solo tiene sentido si el grupo tiene al menos tres estudios; con dos, quitar
  # uno deja un unico estudio y no hay meta-analisis.
  if (length(gse_ids_chr) < 3) {
    message("[", organo_chr, "] leave-one-out omitido (menos de 3 estudios)")
    return(tibble())
  }

  message("\n[", organo_chr, "] Leave-one-out sobre ", length(gse_ids_chr),
          " datasets")

  # Significancia del analisis completo, como referencia.
  estadisticos_completo <- extraer_estadisticos_grupo(objeto_ma)
  z_completo            <- derivar_z_firmado(estadisticos_completo)
  pesos_completo        <- calcular_pesos(estadisticos_completo$tamanos,
                                          "n_efectivo")
  combinado_completo    <- combinar_stouffer_firmado(z_completo$z_firmado,
                                                     pesos_completo) %>%
    mutate(fdr = p.adjust(p_combinado, method = "BH"))

  n_sig_completo_int <- sum(combinado_completo$fdr <
                              parametros_c2_lista$umbral_fdr, na.rm = TRUE)

  # Recalculo quitando cada dataset.
  gse_ids_chr %>%
    map_dfr(function(gse_excluido_chr) {

      message("    Excluyendo: ", gse_excluido_chr)

      objeto_reducido <- objeto_ma[names(objeto_ma) != gse_excluido_chr]

      estadisticos_reducido <- extraer_estadisticos_grupo(objeto_reducido)
      z_reducido            <- derivar_z_firmado(estadisticos_reducido)
      pesos_reducido        <- calcular_pesos(estadisticos_reducido$tamanos,
                                              "n_efectivo")
      combinado_reducido    <- combinar_stouffer_firmado(z_reducido$z_firmado,
                                                         pesos_reducido) %>%
        mutate(fdr = p.adjust(p_combinado, method = "BH"))

      n_sig_reducido_int <- sum(combinado_reducido$fdr <
                                  parametros_c2_lista$umbral_fdr, na.rm = TRUE)

      tibble(
        organo               = organo_chr,
        dataset_excluido     = gse_excluido_chr,
        n_genes_evaluados    = nrow(combinado_reducido),
        n_sig_completo       = n_sig_completo_int,
        n_sig_sin_dataset    = n_sig_reducido_int,
        cambio_absoluto      = n_sig_reducido_int - n_sig_completo_int,
        cambio_relativo_pct  = round(100 * (n_sig_reducido_int -
                                              n_sig_completo_int) /
                                       n_sig_completo_int, 1)
      )
    })
}

message("\n=== LEAVE-ONE-DATASET-OUT ===")

leave_one_out_tbl <- objetos_ma_lista %>%
  imap_dfr(leave_one_out_grupo)

write_csv(leave_one_out_tbl,
          file.path(ruta_salida_chr, "leave_one_out.csv"))


# 09. EXPORTACION DE RESULTADOS ________________________________________________

# Se guarda un resumen por grupo y esquema con el conteo de genes significativos
# antes y despues del filtro de concordancia, que es la salida principal para
# comparar la robustez de la firma.

resumen_firmado_tbl <- meta_firmado_tbl %>%
  group_by(organo, esquema_ponderacion) %>%
  summarise(
    n_genes            = n(),
    n_sig              = sum(fdr < parametros_c2_lista$umbral_fdr, na.rm = TRUE),
    n_sig_concordancia_estricta = sum(
      fdr < parametros_c2_lista$umbral_fdr & pasa_concordancia_estricta,
      na.rm = TRUE
    ),
    n_sig_concordancia_mayoritaria = sum(
      fdr < parametros_c2_lista$umbral_fdr & pasa_concordancia_mayoritaria,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write_csv(resumen_firmado_tbl,
          file.path(ruta_salida_chr, "resumen_stouffer_firmado.csv"))

saveRDS(meta_firmado_tbl,
        file.path(ruta_objetos_chr, "meta_firmado_completo.rds"))


# 10. RESUMEN FINAL ____________________________________________________________

message("\n________________ RESUMEN DE STOUFFER FIRMADO ________________\n")

message("Genes significativos por grupo y esquema (crudo / estricto / mayoritario):")

resumen_firmado_tbl %>%
  pmap(function(organo, esquema_ponderacion, n_genes, n_sig,
                n_sig_concordancia_estricta, n_sig_concordancia_mayoritaria) {
    message("  [", organo, "] ", esquema_ponderacion, ": ",
            n_sig, " / ", n_sig_concordancia_estricta, " / ",
            n_sig_concordancia_mayoritaria, " (de ", n_genes, ")")
  }) %>%
  invisible()

message("\nEfecto del filtro de concordancia:")
message("  El descenso de 'crudo' a 'estricto' cuantifica cuantos genes")
message("  significativos provienen de estudios discordantes en direccion.")

message("\nLeave-one-dataset-out (cambio en genes significativos):")

if (nrow(leave_one_out_tbl) > 0) {
  leave_one_out_tbl %>%
    pmap(function(organo, dataset_excluido, n_genes_evaluados, n_sig_completo,
                  n_sig_sin_dataset, cambio_absoluto, cambio_relativo_pct) {
      message("  [", organo, "] sin ", dataset_excluido, ": ",
              n_sig_completo, " -> ", n_sig_sin_dataset,
              " (", cambio_relativo_pct, "%)")
    }) %>%
    invisible()
}

message("\nComparacion de esquemas de ponderacion:")
message("  Diferencias grandes en n_sig entre esquemas indican dependencia")
message("  del resultado respecto al peso de los datasets grandes (TCGA).")

message("\nTablas exportadas en: ", ruta_salida_chr)
message("___________________________________________________________\n")
