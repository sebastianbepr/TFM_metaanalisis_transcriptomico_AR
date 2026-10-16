# ______________________________________________________________________________
# 06_01  Diagnostico del peso de GSE134073 en el grupo prostatico
#
# Que hace: comprueba si la dominancia de GSE134073 en el grupo prostatico refleja
#           senal biologica concordante o un artefacto, con tres comprobaciones:
#           diseno y tamano por dataset, comportamiento del eje AR y solapamiento
#           de genes top.
# Entra:    Meta_Analisis_Objetos/objetos_ma_por_grupo.rds
# Sale:     Meta_Analisis_Resultados/diag_gse134073_*.csv
# ______________________________________________________________________________

# __ INDEX ___
# 01. Setup, librerias y parametros
# 02. Rutas y carga del objectMA prostatico
# 03. Comprobacion 1: diseno y tamano por dataset
# 04. Estadisticos por estudio para el grupo prostatico
# 05. Comprobacion 2: comportamiento de AR y su eje
# 06. Comprobacion 3: solapamiento de genes top
# 07. Exportacion y resumen final


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
map       <- purrr::map
map_chr   <- purrr::map_chr
map_dfr   <- purrr::map_dfr
imap_dfr  <- purrr::imap_dfr
pmap      <- purrr::pmap
reduce    <- purrr::reduce

# Eje del receptor de androgenos: AR y una seleccion de genes regulados por o
# asociados a la senalizacion androgenica, relevantes para el objetivo final.
genes_eje_ar_chr <- c(
  "AR",        # receptor de androgenos
  "KLK3",      # PSA, diana androgenica clasica
  "KLK2",      # diana androgenica
  "NKX3-1",    # regulado por androgenos, supresor prostatico
  "FKBP5",     # diana androgenica
  "TMPRSS2",   # diana androgenica, fusion en prostata
  "FOLH1",     # PSMA
  "STEAP1",    # diana androgenica
  "SLC45A3",   # regulado por androgenos
  "HOXB13"     # cofactor de AR en prostata
)

# Numero de genes top por dataset para el analisis de solapamiento.
n_top_genes_int <- 500


# 02. RUTAS Y CARGA DEL objectMA PROSTATICO ____________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_salida_chr  <- file.path(ruta_base_chr, "Meta_Analisis_Resultados")
ruta_objetos_chr <- file.path(ruta_base_chr, "Meta_Analisis_Objetos")

objetos_ma_lista <- readRDS(
  file.path(ruta_objetos_chr, "objetos_ma_por_grupo.rds")
)

objeto_prostata <- objetos_ma_lista[["Prostate"]]

gse_ids_chr <- names(objeto_prostata)

message("Datasets prostaticos: ", paste(gse_ids_chr, collapse = ", "))


# 03. COMPROBACION 1: DISENO Y TAMANO POR DATASET ______________________________

# Un dataset puede dominar legitimamente si es mucho mayor o si su contraste es
# mas limpio. Se cuentan tamanos de grupo directamente del vector de condicion.

diseno_tbl <- objeto_prostata %>%
  imap_dfr(function(estudio_lista, gse_id_chr) {

    condicion_int <- estudio_lista$condition
    n_tumor_int   <- sum(condicion_int == 1, na.rm = TRUE)
    n_normal_int  <- sum(condicion_int == 0, na.rm = TRUE)

    tibble(
      gse_id     = gse_id_chr,
      n_genes    = nrow(estudio_lista$mExpres),
      n_tumor    = n_tumor_int,
      n_normal   = n_normal_int,
      n_total    = n_tumor_int + n_normal_int,
      n_efectivo = round((n_tumor_int * n_normal_int) /
                           (n_tumor_int + n_normal_int), 2),
      es_dominante = gse_id_chr == "GSE134073"
    )
  }) %>%
  arrange(desc(n_efectivo))

message("\n--- Comprobacion 1: diseno por dataset ---")
print(diseno_tbl)


# 04. ESTADISTICOS POR ESTUDIO PARA EL GRUPO PROSTATICO ________________________

# calculateES devuelve las matrices de efecto, varianza y logFC por gen y
# estudio, indexadas por simbolo. Son la base de las comprobaciones 2 y 3.

es_prostata_lista <- calculateES(objeto_prostata)

efectos_matriz   <- es_prostata_lista$ES
varianzas_matriz <- es_prostata_lista$Var
logfc_matriz     <- es_prostata_lista$logFC

# calculateES ordena las columnas segun los estudios del objectMA; se les asigna
# el identificador correspondiente para poder referenciarlas por nombre.
colnames(efectos_matriz)   <- gse_ids_chr
colnames(varianzas_matriz) <- gse_ids_chr
colnames(logfc_matriz)     <- gse_ids_chr

# p-valor por gen y estudio, derivado del estadistico de Wald.
wald_matriz <- efectos_matriz / sqrt(varianzas_matriz)
p_matriz    <- 2 * pnorm(-abs(wald_matriz))


# 05. COMPROBACION 2: COMPORTAMIENTO DE AR Y SU EJE ____________________________

# Para cada gen del eje AR se extrae su logFC y p-valor en cada dataset. La
# pregunta clave: la senal de AR es concordante entre datasets, o depende
# exclusivamente de GSE134073? Si AR se comporta parecido en todos, la
# dominancia global de GSE134073 no compromete la conclusion sobre AR.

extraer_gen_eje <- function(gen_chr) {

  if (!gen_chr %in% rownames(logfc_matriz)) {
    return(tibble(gen = gen_chr, gse_id = NA_character_,
                  logfc = NA_real_, p_valor = NA_real_, presente = FALSE))
  }

  tibble(
    gen      = gen_chr,
    gse_id   = colnames(logfc_matriz),
    logfc    = round(as.numeric(logfc_matriz[gen_chr, ]), 3),
    p_valor  = signif(as.numeric(p_matriz[gen_chr, ]), 3),
    presente = TRUE
  )
}

eje_ar_tbl <- genes_eje_ar_chr %>%
  map_dfr(extraer_gen_eje) %>%
  filter(!is.na(gse_id))

# Resumen por gen: en cuantos datasets esta presente, cuantos son significativos
# y si los signos concuerdan. La columna clave es 'concordancia_signo': si es 1,
# el gen se mueve igual en todos los datasets donde se mide.
eje_ar_resumen_tbl <- eje_ar_tbl %>%
  group_by(gen) %>%
  summarise(
    n_datasets_presente = sum(!is.na(logfc)),
    n_significativos    = sum(p_valor < 0.05, na.rm = TRUE),
    signo_concordante   = length(unique(sign(logfc[!is.na(logfc)]))) == 1,
    logfc_gse134073     = logfc[gse_id == "GSE134073"][1],
    logfc_mediana_resto = round(median(logfc[gse_id != "GSE134073"],
                                       na.rm = TRUE), 3),
    .groups = "drop"
  )

message("\n--- Comprobacion 2: eje AR por dataset ---")
print(eje_ar_tbl, n = Inf)

message("\n--- Comprobacion 2: resumen del eje AR ---")
print(eje_ar_resumen_tbl)


# 06. COMPROBACION 3: SOLAPAMIENTO DE GENES TOP ________________________________

# Se comparan los genes mas significativos de GSE134073 con los del resto del
# grupo. Si GSE134073 aporta genes que nadie mas senala, es un outlier; si
# aporta los mismos genes con mas fuerza, es un dataset potente pero concordante.

# Genes top de cada dataset por p-valor.
obtener_top_genes <- function(gse_id_chr) {

  p_columna_num <- p_matriz[, gse_id_chr]
  p_columna_num <- p_columna_num[!is.na(p_columna_num)]

  names(sort(p_columna_num))[seq_len(min(n_top_genes_int,
                                         length(p_columna_num)))]
}

top_por_dataset_lista <- gse_ids_chr %>%
  set_names() %>%
  map(obtener_top_genes)

top_gse134073_chr <- top_por_dataset_lista[["GSE134073"]]

top_resto_chr <- gse_ids_chr[gse_ids_chr != "GSE134073"] %>%
  map(~ top_por_dataset_lista[[.x]]) %>%
  reduce(union)

# Fraccion de los genes top de GSE134073 que tambien aparecen en el top de algun
# otro dataset. Un valor alto indica concordancia; uno bajo, comportamiento
# outlier.
n_compartidos_int <- length(intersect(top_gse134073_chr, top_resto_chr))

solapamiento_num <- n_compartidos_int / length(top_gse134073_chr)

# Solapamiento por pares, para ver si GSE134073 se parece a algun dataset en
# particular o a ninguno.
solapamiento_pares_tbl <- gse_ids_chr[gse_ids_chr != "GSE134073"] %>%
  map_dfr(function(gse_id_chr) {
    top_otro_chr <- top_por_dataset_lista[[gse_id_chr]]
    n_comun_int  <- length(intersect(top_gse134073_chr, top_otro_chr))
    tibble(
      dataset_comparado = gse_id_chr,
      n_top_comun       = n_comun_int,
      pct_solapamiento  = round(100 * n_comun_int / length(top_gse134073_chr), 1)
    )
  }) %>%
  arrange(desc(pct_solapamiento))

message("\n--- Comprobacion 3: solapamiento de genes top ---")
message("Genes top de GSE134073 compartidos con algun otro dataset: ",
        n_compartidos_int, "/", length(top_gse134073_chr),
        " (", round(100 * solapamiento_num, 1), "%)")

print(solapamiento_pares_tbl)


# 07. EXPORTACION Y RESUMEN FINAL ______________________________________________

write_csv(diseno_tbl,
          file.path(ruta_salida_chr, "diag_gse134073_diseno.csv"))
write_csv(eje_ar_tbl,
          file.path(ruta_salida_chr, "diag_gse134073_eje_ar.csv"))
write_csv(eje_ar_resumen_tbl,
          file.path(ruta_salida_chr, "diag_gse134073_eje_ar_resumen.csv"))
write_csv(solapamiento_pares_tbl,
          file.path(ruta_salida_chr, "diag_gse134073_solapamiento.csv"))

message("\n________________ VEREDICTO DEL DIAGNOSTICO ________________\n")

# Sintesis interpretativa automatica de las tres comprobaciones.

ar_concordante_lgl <- eje_ar_resumen_tbl %>%
  filter(gen == "AR") %>%
  pull(signo_concordante)

if (length(ar_concordante_lgl) == 0) ar_concordante_lgl <- NA

message("1. Diseno: GSE134073 tiene n_efectivo de ",
        diseno_tbl$n_efectivo[diseno_tbl$gse_id == "GSE134073"],
        " (rango del grupo: ", min(diseno_tbl$n_efectivo), " a ",
        max(diseno_tbl$n_efectivo), ")")

message("2. AR concordante en signo entre datasets: ",
        if_else(isTRUE(ar_concordante_lgl), "SI",
                if_else(is.na(ar_concordante_lgl),
                        "AR no presente en el grupo", "NO")))

message("3. Solapamiento de genes top de GSE134073 con el resto: ",
        round(100 * solapamiento_num, 1), "%")

message("\nInterpretacion:")
message("  - Solapamiento alto (>50%) + AR concordante -> GSE134073 es un")
message("    dataset potente pero concordante; su dominancia no sesga AR.")
message("  - Solapamiento bajo + AR solo en GSE134073 -> outlier que requiere")
message("    ponderacion a la baja o exclusion antes del analisis de AR.")

message("\nTablas exportadas en: ", ruta_salida_chr)
message("___________________________________________________________\n")
