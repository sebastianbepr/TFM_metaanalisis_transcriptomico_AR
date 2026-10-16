# ______________________________________________________________________________
# 05_04  Contraste entre grupos: firma hormonal compartida (Firma 1) y divergencia (Firma 2)
#
# Que hace: identifica los genes con comportamiento comun en prostata y testiculo
#           (firma hormonal compartida) que divergen en vejiga, en tres versiones
#           prostaticas para sensibilidad, y ubica el eje del receptor de androgenos.
# Entra:    Datasets_Remediados/*_remediado.rds
# Sale:     Meta_Analisis_Objetos/contraste_grupos.rds,
#           contraste_prostata_4ds/6ds/3ds.csv, contraste_resumen.csv
# ______________________________________________________________________________

# __ INDEX ___
# 01. Setup, librerias y parametros
# 02. Rutas y carga de los datasets remediados
# 03. Definicion de las versiones de cada grupo
# 04. Motor de meta-analisis firmado (reutilizado de C2)
# 05. Calculo de la firma firmada por grupo
# 06. Interseccion prostata-testiculo (firma hormonal compartida)
# 07. Clasificacion de la divergencia en vejiga
# 08. Ubicacion del eje del receptor de androgenos
# 09. Analisis de sensibilidad entre versiones prostaticas
# 10. Exportacion y resumen final


# 01. SETUP ____________________________________________________________________

library(tidyverse)
library(DExMA)

select    <- dplyr::select
filter    <- dplyr::filter
count     <- dplyr::count
slice     <- dplyr::slice
arrange   <- dplyr::arrange
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
group_by  <- dplyr::group_by
rename    <- dplyr::rename
map       <- purrr::map
map_chr   <- purrr::map_chr
map_dbl   <- purrr::map_dbl
map_dfr   <- purrr::map_dfr
imap      <- purrr::imap
imap_dfr  <- purrr::imap_dfr
keep      <- purrr::keep
reduce    <- purrr::reduce

parametros_d_lista <- list(
  umbral_fdr            = 0.05,

  # Concordancia de signo minima para considerar consistente un gen dentro de un
  # grupo (proporcion de estudios que coinciden con la direccion combinada).
  concordancia_minima   = 0.75,

  # Umbral de logFC combinado por debajo del cual se considera que un gen no
  # cambia de forma biologicamente relevante en vejiga (para la categoria de
  # ausencia por equivalencia).
  umbral_logfc_ausencia = 0.30,

  semilla_aleatoria     = 1234
)

set.seed(parametros_d_lista$semilla_aleatoria)

genes_eje_ar_chr <- c("AR", "KLK3", "KLK2", "NKX3-1", "FKBP5",
                      "TMPRSS2", "FOLH1", "STEAP1", "SLC45A3", "HOXB13")


# 02. RUTAS Y CARGA DE LOS DATASETS REMEDIADOS _________________________________

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
    gse_id           = remediado_lista$dataset_descriptor$gse_id,
    organo           = remediado_lista$dataset_descriptor$organo,
    expresion_matriz = remediado_lista$expresion_matriz,
    grupo_binario    = remediado_lista$grupo_binario
  )
}

datasets_lista <- archivos_remediados_chr %>%
  map(cargar_remediado)

names(datasets_lista) <- map_chr(datasets_lista, ~ .x$gse_id)


# 03. DEFINICION DE LAS VERSIONES DE CADA GRUPO ________________________________

# Cada version es un conjunto de identificadores de dataset. Prostata tiene tres
# versiones por el problema del comparador HPB; testiculo y vejiga son unicas.

versiones_lista <- list(
  prostata_4ds = c("GSE104131", "GSE114740", "GSE214940", "GSE221107"),
  prostata_6ds = c("GSE104131", "GSE114740", "GSE134073", "GSE214940",
                   "GSE217260", "GSE221107"),
  prostata_3ds = c("GSE104131", "GSE114740", "GSE214940"),
  testiculo    = c("GSE3218", "GSE8607", "TCGA_GTEX_Testis"),
  vejiga       = c("GSE7476", "GSE40355")
)


# 04. MOTOR DE META-ANALISIS FIRMADO ___________________________________________

# Construye un objectMA a partir de una lista de identificadores de dataset.

construir_objectma <- function(gse_ids_chr) {

  datasets_grupo_lista <- datasets_lista[gse_ids_chr]

  matrices_lista <- datasets_grupo_lista %>%
    map(~ .x$expresion_matriz)

  phenodatas_lista <- datasets_grupo_lista %>%
    map(function(dataset_lista) {
      condicion_chr <- if_else(dataset_lista$grupo_binario == 1,
                               "Tumor", "Normal")
      data.frame(condicion = condicion_chr,
                 row.names = names(dataset_lista$grupo_binario),
                 stringsAsFactors = FALSE)
    })

  n_int <- length(matrices_lista)

  createObjectMA(
    listEX    = matrices_lista,
    listPheno = phenodatas_lista,
    namePheno = rep("condicion", n_int),
    expGroups = rep("Tumor",     n_int),
    refGroups = rep("Normal",    n_int)
  )
}


# Meta-analisis de Stouffer con signo, ponderado por raiz de n efectivo. Devuelve
# una tabla por gen con Z combinado, p, FDR, direccion y concordancia de signo.
# Es el mismo procedimiento de C2, encapsulado para reutilizarlo por grupo.

meta_firmado <- function(objeto_ma) {

  es_lista     <- calculateES(objeto_ma)
  ef_matriz    <- es_lista$ES
  var_matriz   <- es_lista$Var
  logfc_matriz <- es_lista$logFC

  # Tamanos por estudio (n efectivo) desde el vector de condicion.
  n_efectivo_num <- objeto_ma %>%
    map_dbl(function(estudio_lista) {
      cond_int <- estudio_lista$condition
      n_t <- sum(cond_int == 1, na.rm = TRUE)
      n_c <- sum(cond_int == 0, na.rm = TRUE)
      (n_t * n_c) / (n_t + n_c)
    })
  pesos_num <- sqrt(n_efectivo_num)

  # z firmado por gen y estudio.
  wald_matriz       <- ef_matriz / sqrt(var_matriz)
  p_ind_matriz      <- 2 * pnorm(-abs(wald_matriz))
  z_una_cola_matriz <- qnorm(1 - p_ind_matriz / 2)
  z_firmado_matriz  <- sign(logfc_matriz) * z_una_cola_matriz

  # Combinacion de Stouffer ponderada, gen a gen.
  n_genes_int <- nrow(z_firmado_matriz)

  z_comb_num <- map_dbl(seq_len(n_genes_int), function(i) {
    z_num <- z_firmado_matriz[i, ]
    ok_lgl <- !is.na(z_num) & is.finite(z_num)
    if (!any(ok_lgl)) return(NA_real_)
    sum(pesos_num[ok_lgl] * z_num[ok_lgl]) / sqrt(sum(pesos_num[ok_lgl]^2))
  })

  p_comb_num <- 2 * pnorm(-abs(z_comb_num))

  # logFC combinado ponderado por inverso de varianza.
  pesos_iv_matriz <- 1 / var_matriz
  logfc_comb_num  <- rowSums(logfc_matriz * pesos_iv_matriz, na.rm = TRUE) /
    rowSums(pesos_iv_matriz, na.rm = TRUE)

  # Concordancia de signo: proporcion de estudios que coinciden con el signo del
  # Z combinado.
  direccion_int <- sign(z_comb_num)
  concordancia_num <- map_dbl(seq_len(n_genes_int), function(i) {
    signos_num <- sign(logfc_matriz[i, ])
    signos_num <- signos_num[!is.na(signos_num)]
    if (length(signos_num) == 0) return(NA_real_)
    mean(signos_num == direccion_int[i])
  })

  tibble(
    gen_simbolo      = rownames(z_firmado_matriz),
    z_combinado      = round(z_comb_num, 4),
    p_combinado      = p_comb_num,
    fdr              = p.adjust(p_comb_num, method = "BH"),
    logfc_combinado  = round(logfc_comb_num, 4),
    direccion        = if_else(z_comb_num > 0, "Up", "Down"),
    concordancia     = round(concordancia_num, 3)
  )
}


# 05. CALCULO DE LA FIRMA FIRMADA POR GRUPO ____________________________________

message("=== META-ANALISIS FIRMADO POR VERSION ===")

firmas_lista <- versiones_lista %>%
  imap(function(gse_ids_chr, nombre_version_chr) {
    message("[", nombre_version_chr, "] ", length(gse_ids_chr), " datasets")
    objeto_ma <- construir_objectma(gse_ids_chr)
    objeto_ma <- dataLog(objeto_ma)
    meta_firmado(objeto_ma)
  })


# Marca de significancia y consistencia de cada gen dentro de su grupo. Un gen es
# "consistente" si es significativo y su concordancia de signo supera el umbral.

marcar_consistencia <- function(firma_tbl) {
  firma_tbl %>%
    mutate(
      significativo = fdr < parametros_d_lista$umbral_fdr,
      consistente   = significativo &
        concordancia >= parametros_d_lista$concordancia_minima
    )
}

firmas_lista <- firmas_lista %>%
  map(marcar_consistencia)


# 06. INTERSECCION PROSTATA-TESTICULO __________________________________________

# Para una version prostatica dada, se cruza su firma consistente con la de
# testiculo. La firma compartida son los genes consistentes en AMBOS grupos, con
# la misma direccion (ambos Up o ambos Down): el patron hormonal comun.

calcular_interseccion <- function(firma_prostata_tbl, firma_testiculo_tbl) {

  prostata_sig_tbl <- firma_prostata_tbl %>%
    filter(consistente) %>%
    select(gen_simbolo,
           z_prostata     = z_combinado,
           logfc_prostata = logfc_combinado,
           dir_prostata   = direccion)

  testiculo_sig_tbl <- firma_testiculo_tbl %>%
    filter(consistente) %>%
    select(gen_simbolo,
           z_testiculo     = z_combinado,
           logfc_testiculo = logfc_combinado,
           dir_testiculo   = direccion)

  # Interseccion por simbolo, exigiendo misma direccion en ambos grupos.
  prostata_sig_tbl %>%
    inner_join(testiculo_sig_tbl, by = "gen_simbolo") %>%
    filter(dir_prostata == dir_testiculo) %>%
    mutate(direccion_hormonal = dir_prostata)
}


# 07. CLASIFICACION DE LA DIVERGENCIA EN VEJIGA ________________________________

# A cada gen de la firma hormonal compartida se le anade su comportamiento en
# vejiga, y se clasifica en dos categorias:
#   (A) Inversion : significativo en vejiga con direccion opuesta a la hormonal.
#   (B) Ausencia  : no significativo en vejiga y con |logFC| por debajo del
#                   umbral de relevancia (ausencia de cambio por equivalencia).
# Los genes que no caen en ninguna (significativos y en la misma direccion que el
# eje hormonal) se marcan como "concordante", es decir, NO divergentes.

clasificar_divergencia <- function(interseccion_tbl, firma_vejiga_tbl) {

  vejiga_tbl <- firma_vejiga_tbl %>%
    select(gen_simbolo,
           z_vejiga     = z_combinado,
           logfc_vejiga = logfc_combinado,
           dir_vejiga   = direccion,
           fdr_vejiga   = fdr,
           sig_vejiga   = significativo)

  interseccion_tbl %>%
    left_join(vejiga_tbl, by = "gen_simbolo") %>%
    mutate(
      categoria_vejiga = case_when(
        # No medido en vejiga: no clasificable.
        is.na(sig_vejiga) ~ "No evaluado en vejiga",
        # (A) Inversion: significativo y con signo opuesto al hormonal.
        sig_vejiga & dir_vejiga != direccion_hormonal ~ "A_inversion",
        # (B) Ausencia por equivalencia: no significativo y cambio pequeno.
        !sig_vejiga &
          abs(logfc_vejiga) < parametros_d_lista$umbral_logfc_ausencia ~
          "B_ausencia",
        # No significativo pero con cambio apreciable: ausencia debil.
        !sig_vejiga ~ "B_ausencia_debil",
        # Significativo y misma direccion: concordante, no divergente.
        TRUE ~ "Concordante (no divergente)"
      )
    )
}


# 08. EJECUCION DEL CONTRASTE PARA CADA VERSION PROSTATICA ______________________

# Se ejecuta la interseccion y clasificacion para las tres versiones prostaticas
# contra la misma firma de testiculo y vejiga.

versiones_prostata_chr <- c("prostata_4ds", "prostata_6ds", "prostata_3ds")

contraste_lista <- versiones_prostata_chr %>%
  set_names() %>%
  map(function(version_chr) {

    message("\n[", version_chr, "] Calculando interseccion y divergencia")

    interseccion_tbl <- calcular_interseccion(
      firmas_lista[[version_chr]],
      firmas_lista[["testiculo"]]
    )

    clasificar_divergencia(interseccion_tbl, firmas_lista[["vejiga"]])
  })


# 09. UBICACION DEL EJE DEL RECEPTOR DE ANDROGENOS _____________________________

# Se localizan los genes del eje AR dentro del contraste principal (4 datasets),
# reportando si forman parte de la firma hormonal compartida y su categoria de
# divergencia. Es el puente hacia el analisis especifico de AR del paso siguiente.

eje_ar_en_contraste_tbl <- contraste_lista[["prostata_4ds"]] %>%
  filter(gen_simbolo %in% genes_eje_ar_chr) %>%
  select(gen_simbolo, direccion_hormonal, logfc_prostata, logfc_testiculo,
         logfc_vejiga, categoria_vejiga)

# Tambien se reporta el comportamiento de AR en cada grupo aunque no supere el
# filtro de interseccion, para no perderlo de vista.
eje_ar_por_grupo_tbl <- genes_eje_ar_chr %>%
  map_dfr(function(gen_chr) {
    extraer <- function(firma_tbl, etiqueta_chr) {
      fila <- firma_tbl %>% filter(gen_simbolo == gen_chr)
      if (nrow(fila) == 0) return(tibble())
      tibble(gen = gen_chr, grupo = etiqueta_chr,
             logfc = fila$logfc_combinado, dir = fila$direccion,
             fdr = round(fila$fdr, 4), consistente = fila$consistente)
    }
    bind_rows(
      extraer(firmas_lista[["prostata_4ds"]], "Prostata(4ds)"),
      extraer(firmas_lista[["testiculo"]],    "Testiculo"),
      extraer(firmas_lista[["vejiga"]],       "Vejiga")
    )
  })


# 10. EXPORTACION Y RESUMEN FINAL ______________________________________________

# Se guarda el contraste completo de cada version y los resumenes.

contraste_lista %>%
  iwalk(function(contraste_tbl, version_chr) {
    write_csv(contraste_tbl,
              file.path(ruta_salida_chr,
                        paste0("contraste_", version_chr, ".csv")))
  })

write_csv(eje_ar_en_contraste_tbl,
          file.path(ruta_salida_chr, "contraste_eje_ar.csv"))
write_csv(eje_ar_por_grupo_tbl,
          file.path(ruta_salida_chr, "eje_ar_por_grupo.csv"))

saveRDS(list(firmas = firmas_lista, contraste = contraste_lista),
        file.path(ruta_objetos_chr, "contraste_grupos.rds"))

# Resumen comparativo: tamano de la firma hormonal y reparto por categoria de
# divergencia, para las tres versiones prostaticas.

resumen_contraste_tbl <- contraste_lista %>%
  imap_dfr(function(contraste_tbl, version_chr) {
    contraste_tbl %>%
      count(categoria_vejiga, name = "n") %>%
      mutate(version = version_chr, .before = 1)
  })

write_csv(resumen_contraste_tbl,
          file.path(ruta_salida_chr, "contraste_resumen.csv"))


message("\n________________ CONTRASTE ENTRE GRUPOS ________________\n")

message("Tamano de la firma hormonal compartida (prostata ∩ testiculo):")

contraste_lista %>%
  iwalk(function(contraste_tbl, version_chr) {
    message("  ", version_chr, ": ", nrow(contraste_tbl), " genes compartidos")
  })

message("\nDivergencia en vejiga por version (categorias):")

resumen_contraste_tbl %>%
  arrange(version, categoria_vejiga) %>%
  pmap(function(version, categoria_vejiga, n) {
    message("  [", version, "] ", categoria_vejiga, ": ", n)
  }) %>%
  invisible()

message("\nEje del receptor de androgenos en el contraste principal (4ds):")

if (nrow(eje_ar_en_contraste_tbl) > 0) {
  eje_ar_en_contraste_tbl %>%
    pmap(function(gen_simbolo, direccion_hormonal, logfc_prostata,
                  logfc_testiculo, logfc_vejiga, categoria_vejiga) {
      message("  ", gen_simbolo, ": ", direccion_hormonal,
              " (P=", logfc_prostata, ", T=", logfc_testiculo,
              ", V=", logfc_vejiga, ") -> ", categoria_vejiga)
    }) %>%
    invisible()
} else {
  message("  Ningun gen del eje AR supero el filtro de interseccion.")
  message("  Ver eje_ar_por_grupo.csv para su comportamiento por grupo.")
}

message("\nTablas exportadas en: ", ruta_salida_chr)
message("_______________________________________________________\n")
