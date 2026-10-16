# ______________________________________________________________________________
# 07_01  Tabla maestra de la firma hormonal compartida
#
# Que hace: construye la tabla final del meta-analisis con los genes de la firma
#           compartida, su categoria de divergencia en vejiga y el marcado del
#           eje androgenico, y la exporta en CSV y XLSX.
# Entra:    Meta_Analisis_Resultados/contraste_prostata_4ds.csv
# Sale:     Meta_Analisis_Entregables/tabla_maestra_firma_compartida.csv/.xlsx
# ______________________________________________________________________________


# 01. SETUP Y RUTAS ____________________________________________________________

library(tidyverse)
library(writexl)   # para exportar a XLSX sin dependencias de Java

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_salida_chr   <- file.path(ruta_base_chr, "Meta_Analisis_Resultados")
ruta_entregables_chr <- file.path(ruta_base_chr, "Meta_Analisis_Entregables")

if (!dir.exists(ruta_entregables_chr)) {
  dir.create(ruta_entregables_chr, recursive = TRUE)
}

genes_eje_ar_chr <- c("AR", "KLK3", "KLK2", "NKX3-1", "FKBP5",
                      "TMPRSS2", "FOLH1", "STEAP1", "SLC45A3", "HOXB13")


# 02. CARGA DEL CONTRASTE PRINCIPAL ____________________________________________

contraste_tbl <- read_csv(
  file.path(ruta_salida_chr, "contraste_prostata_4ds.csv"),
  show_col_types = FALSE
)


# 03. CONSTRUCCION DE LA TABLA MAESTRA _________________________________________

# Se traduce la categoria interna a una etiqueta legible y se ordena de forma
# que los genes divergentes (el resultado de interes) aparezcan primero.

etiquetas_categoria_chr <- c(
  "A_inversion"                 = "Divergente: inversion en vejiga",
  "B_ausencia"                  = "Divergente: ausencia en vejiga",
  "B_ausencia_debil"            = "Divergente: ausencia debil en vejiga",
  "Concordante (no divergente)" = "Concordante en los tres grupos",
  "No evaluado en vejiga"       = "No evaluado en vejiga"
)

orden_categoria_chr <- names(etiquetas_categoria_chr)

tabla_maestra_tbl <- contraste_tbl %>%
  mutate(
    categoria_legible = etiquetas_categoria_chr[categoria_vejiga],
    categoria_vejiga  = factor(categoria_vejiga, levels = orden_categoria_chr),
    es_divergente     = categoria_vejiga %in%
      c("A_inversion", "B_ausencia", "B_ausencia_debil")
  ) %>%
  arrange(categoria_vejiga, desc(abs(logfc_prostata))) %>%
  select(
    gen_simbolo,
    direccion_hormonal,
    logfc_prostata,
    logfc_testiculo,
    logfc_vejiga,
    categoria       = categoria_legible,
    es_divergente
  )


# 04. MARCADO DEL EJE ANDROGENICO ______________________________________________

tabla_maestra_tbl <- tabla_maestra_tbl %>%
  mutate(eje_androgenico = gen_simbolo %in% genes_eje_ar_chr)


# 05. EXPORTACION ______________________________________________________________

write_csv(tabla_maestra_tbl,
          file.path(ruta_entregables_chr, "tabla_maestra_firma_compartida.csv"))

write_xlsx(
  list(
    "Firma compartida (287 genes)" = tabla_maestra_tbl,
    "Solo divergentes"             = filter(tabla_maestra_tbl, es_divergente),
    "Eje androgenico"              = filter(tabla_maestra_tbl, eje_androgenico)
  ),
  file.path(ruta_entregables_chr, "tabla_maestra_firma_compartida.xlsx")
)


message("________________ TABLA MAESTRA GENERADA ________________\n")

message("Total de genes en la firma compartida: ", nrow(tabla_maestra_tbl))

message("\nReparto por categoria:")
tabla_maestra_tbl %>%
  count(categoria, name = "n") %>%
  pmap(function(categoria, n) {
    message("  ", categoria, ": ", n)
  }) %>%
  invisible()

message("\nGenes divergentes (comunes en prostata/testiculo, distintos en vejiga): ",
        sum(tabla_maestra_tbl$es_divergente))

message("\nGenes del eje androgenico en la firma: ",
        sum(tabla_maestra_tbl$eje_androgenico))

message("\nEntregables guardados en: ", ruta_entregables_chr)
message("________________________________________________________\n")
