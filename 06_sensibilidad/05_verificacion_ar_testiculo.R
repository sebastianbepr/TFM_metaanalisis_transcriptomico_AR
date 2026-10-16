# ______________________________________________________________________________
# 06_05  Verificacion de AR por dataset individual en testiculo
#
# Que hace: comprueba si la caida de AR en el grupo testicular la reproducen los
#           tres datasets por separado (limma dataset a dataset) o depende solo
#           de TCGA/GTEx.
# Entra:    Datasets_Remediados/GSE3218_remediado.rds, GSE8607_remediado.rds,
#           TCGA_GTEX_Testis_remediado.rds
# Sale:     Meta_Analisis_Resultados/verif_ar_testiculo_*.csv
# ______________________________________________________________________________


# 01. SETUP ____________________________________________________________________

library(tidyverse)
library(limma)

select   <- dplyr::select
filter   <- dplyr::filter
arrange  <- dplyr::arrange
mutate   <- dplyr::mutate
map      <- purrr::map
map_dfr  <- purrr::map_dfr
imap_dfr <- purrr::imap_dfr
pmap     <- purrr::pmap

# Genes del eje androgenico a rastrear.
genes_eje_ar_chr <- c("AR", "KLK3", "KLK2", "NKX3-1", "FKBP5",
                      "TMPRSS2", "FOLH1", "STEAP1", "SLC45A3", "HOXB13")

# Datasets del grupo testicular.
datasets_testiculo_chr <- c("GSE3218", "GSE8607", "TCGA_GTEX_Testis")


# 02. RUTAS Y CARGA ____________________________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_remediados_chr <- file.path(ruta_base_chr, "Datasets_Remediados")
ruta_salida_chr     <- file.path(ruta_base_chr, "Meta_Analisis_Resultados")

cargar_remediado <- function(gse_id_chr) {
  ruta_chr <- file.path(ruta_remediados_chr,
                        paste0(gse_id_chr, "_remediado.rds"))
  readRDS(ruta_chr)
}

datasets_testiculo_lista <- datasets_testiculo_chr %>%
  set_names() %>%
  map(cargar_remediado)


# 03. EXPRESION DIFERENCIAL POR DATASET ________________________________________

# Se ejecuta limma en cada dataset por separado, sobre la matriz remediada (ya en
# escala log). El contraste es Tumor vs Normal, igual que el meta-analisis. El
# objetivo no es re-hacer el meta-analisis, sino ver el efecto crudo de cada
# dataset sobre los genes del eje AR.

dge_por_dataset <- function(remediado_lista, gse_id_chr) {

  message("[", gse_id_chr, "] Ejecutando limma")

  expresion_matriz  <- remediado_lista$expresion_matriz
  grupo_binario_int <- remediado_lista$grupo_binario

  # Diseno: 0 = Normal (referencia), 1 = Tumor.
  diseno_matriz <- model.matrix(~ grupo_binario_int)

  ajuste <- lmFit(expresion_matriz, diseno_matriz)
  ajuste <- eBayes(ajuste)

  # El coeficiente 2 corresponde al efecto Tumor vs Normal.
  resultado_tbl <- topTable(ajuste, coef = 2, number = Inf,
                            sort.by = "none") %>%
    as.data.frame() %>%
    rownames_to_column("gen_simbolo") %>%
    as_tibble() %>%
    mutate(gse_id = gse_id_chr)

  resultado_tbl
}

dge_testiculo_tbl <- datasets_testiculo_lista %>%
  imap_dfr(dge_por_dataset)


# 04. EXTRACCION DEL EJE AR POR DATASET ________________________________________

eje_ar_por_dataset_tbl <- dge_testiculo_tbl %>%
  filter(gen_simbolo %in% genes_eje_ar_chr) %>%
  select(gen_simbolo, gse_id,
         logFC, p_valor = P.Value, fdr = adj.P.Val) %>%
  mutate(
    logFC     = round(logFC, 3),
    fdr       = signif(fdr, 3),
    direccion = if_else(logFC > 0, "Up", "Down"),
    significativo = fdr < 0.05
  ) %>%
  arrange(gen_simbolo, gse_id)


# 05. VERIFICACION DE CONSISTENCIA ENTRE PLATAFORMAS ___________________________

# Para cada gen del eje, se comprueba si los tres datasets coinciden en la
# direccion del cambio. La pregunta clave es sobre AR: los dos microarrays
# acompanan a TCGA en Down, o TCGA arrastra el resultado en solitario?

consistencia_tbl <- eje_ar_por_dataset_tbl %>%
  group_by(gen_simbolo) %>%
  summarise(
    n_datasets       = n(),
    n_up             = sum(direccion == "Up"),
    n_down           = sum(direccion == "Down"),
    direccion_comun  = case_when(
      n_up == n_datasets   ~ "Todos Up",
      n_down == n_datasets ~ "Todos Down",
      TRUE                 ~ "Discordante"
    ),
    logfc_tcga       = logFC[gse_id == "TCGA_GTEX_Testis"][1],
    logfc_gse3218    = logFC[gse_id == "GSE3218"][1],
    logfc_gse8607    = logFC[gse_id == "GSE8607"][1],
    .groups = "drop"
  ) %>%
  arrange(gen_simbolo)


# 06. EXPORTACION Y RESUMEN FINAL ______________________________________________

write_csv(eje_ar_por_dataset_tbl,
          file.path(ruta_salida_chr, "verif_ar_testiculo_por_dataset.csv"))
write_csv(consistencia_tbl,
          file.path(ruta_salida_chr, "verif_ar_testiculo_consistencia.csv"))


message("\n________________ VERIFICACION DE AR EN TESTICULO ________________\n")

message("--- AR y dianas androgenicas por dataset ---")
print(eje_ar_por_dataset_tbl, n = Inf)

message("\n--- Consistencia de direccion entre datasets ---")
print(consistencia_tbl, n = Inf)

# Veredicto especifico sobre AR.
ar_fila <- consistencia_tbl %>% filter(gen_simbolo == "AR")

message("\n--- Veredicto sobre AR ---")

if (nrow(ar_fila) > 0) {
  message("  Direccion en TCGA:    logFC ", ar_fila$logfc_tcga)
  message("  Direccion en GSE3218: logFC ", ar_fila$logfc_gse3218)
  message("  Direccion en GSE8607: logFC ", ar_fila$logfc_gse8607)
  message("  Consenso: ", ar_fila$direccion_comun)

  if (ar_fila$direccion_comun == "Todos Down") {
    message("\n  CONFIRMADO: los tres datasets coinciden en AR Down.")
    message("  La caida de AR en testiculo NO es un artefacto de TCGA;")
    message("  esta respaldada por los dos microarrays independientes.")
  } else if (ar_fila$direccion_comun == "Discordante") {
    message("\n  ATENCION: los datasets discrepan en la direccion de AR.")
    message("  Revisar si TCGA arrastra el resultado en solitario.")
  }
}

message("\nTablas exportadas en: ", ruta_salida_chr)
message("________________________________________________________________\n")
