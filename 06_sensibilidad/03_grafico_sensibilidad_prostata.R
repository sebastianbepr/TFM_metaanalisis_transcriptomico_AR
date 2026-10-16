# ______________________________________________________________________________
# 06_03  Figura del analisis de sensibilidad del grupo prostatico
#
# Que hace: figura de dos paneles que compara el meta-analisis prostatico con y
#           sin GSE134073: efecto sobre la firma global y sobre el eje AR.
# Entra:    Meta_Analisis_Resultados/prostata_resumen_comparativo.csv,
#           prostata_eje_ar_comparativo.csv
# Sale:     Meta_Analisis_Graficos/grafico_sensibilidad_prostata.png/.pdf
# ______________________________________________________________________________


# 01. SETUP ____________________________________________________________________

library(tidyverse)
library(patchwork)

# Paleta sobria: rojo apagado para la version con el dataset, azul para sin el.
col_con_dataset_chr <- "#C1666B"
col_sin_dataset_chr <- "#4A7C8C"

# Tema comun a ambos paneles, minimalista y apto para documento.
tema_tfm <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(face = "bold", size = 12),
    plot.title.position = "plot",
    axis.title         = element_text(size = 10),
    legend.position    = "top",
    legend.title       = element_blank()
  )


# 02. RUTAS Y CARGA DE DATOS ___________________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_salida_chr   <- file.path(ruta_base_chr, "Meta_Analisis_Resultados")
ruta_graficos_chr <- file.path(ruta_base_chr, "Meta_Analisis_Graficos")

if (!dir.exists(ruta_graficos_chr)) {
  dir.create(ruta_graficos_chr, recursive = TRUE)
}

resumen_tbl <- read_csv(
  file.path(ruta_salida_chr, "prostata_resumen_comparativo.csv"),
  show_col_types = FALSE
)

eje_ar_tbl <- read_csv(
  file.path(ruta_salida_chr, "prostata_eje_ar_comparativo.csv"),
  show_col_types = FALSE
)


# 03. PANEL A: FIRMA GLOBAL POR BRAZO Y METODO _________________________________

# Se pasa a formato largo para tener una fila por combinacion y version, con las
# etiquetas ordenadas de forma logica (brazo x metodo).

orden_brazo_chr  <- c("comun", "proporcion", "imputacion")
orden_metodo_chr <- c("Fisher", "Stouffer", "maxP")

firma_larga_tbl <- resumen_tbl %>%
  select(brazo, metodo,
         n_significativos_6ds, n_significativos_5ds) %>%
  pivot_longer(
    cols      = starts_with("n_significativos"),
    names_to  = "version",
    values_to = "n_significativos"
  ) %>%
  mutate(
    version = if_else(version == "n_significativos_6ds",
                      "Con GSE134073 (6 datasets)",
                      "Sin GSE134073 (5 datasets)"),
    brazo   = factor(brazo, levels = orden_brazo_chr),
    metodo  = factor(metodo, levels = orden_metodo_chr),
    etiqueta_combo = paste0(brazo, "\n", metodo)
  ) %>%
  arrange(brazo, metodo) %>%
  mutate(etiqueta_combo = fct_inorder(etiqueta_combo))

panel_a <- firma_larga_tbl %>%
  ggplot(aes(x = etiqueta_combo, y = n_significativos, fill = version)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_text(
    aes(label = n_significativos),
    position = position_dodge(width = 0.8),
    vjust = -0.4, size = 2.7, color = "grey30"
  ) +
  scale_fill_manual(values = c(col_con_dataset_chr, col_sin_dataset_chr)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "A. Firma global: efecto de excluir GSE134073",
    x     = NULL,
    y     = "Genes significativos (FDR < 0.05)"
  ) +
  tema_tfm +
  theme(axis.text.x = element_text(size = 8))


# 04. PANEL B: EJE DEL RECEPTOR DE ANDROGENOS __________________________________

# Grafico de mancuernas (dumbbell): para cada gen del eje, el logFC combinado en
# cada version, unidos por una linea. AR se resalta por ser el gen de interes.

eje_ancho_tbl <- eje_ar_tbl %>%
  select(gen, version, logfc_combinado) %>%
  pivot_wider(names_from = version, values_from = logfc_combinado) %>%
  arrange(`5ds`) %>%
  mutate(
    gen        = fct_inorder(gen),
    es_ar      = gen == "AR"
  )

panel_b <- eje_ancho_tbl %>%
  ggplot(aes(y = gen)) +
  # Banda de resaltado para AR.
  geom_tile(
    data  = filter(eje_ancho_tbl, es_ar),
    aes(x = 0.8, width = Inf, height = 0.9),
    fill  = "#F0E6A8", alpha = 0.4
  ) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_segment(
    aes(x = `5ds`, xend = `6ds`, yend = gen),
    color = "grey75", linewidth = 1
  ) +
  geom_point(aes(x = `6ds`, color = "Con GSE134073 (6 datasets)"), size = 2.6) +
  geom_point(aes(x = `5ds`, color = "Sin GSE134073 (5 datasets)"), size = 2.6) +
  scale_color_manual(values = c(col_con_dataset_chr, col_sin_dataset_chr)) +
  labs(
    title = "B. Eje del receptor de androgenos",
    x     = "logFC combinado (tumor vs control)",
    y     = NULL,
    color = NULL
  ) +
  tema_tfm +
  theme(
    axis.text.y = element_text(
      face = if_else(levels(eje_ancho_tbl$gen) == "AR", "bold", "plain")
    )
  )


# 05. COMPOSICION Y EXPORTACION ________________________________________________

figura_sensibilidad <- (panel_a | panel_b) +
  plot_layout(widths = c(1.15, 1)) +
  plot_annotation(
    title = paste("Analisis de sensibilidad del grupo prostatico -",
                  "inclusion de GSE134073 (comparador HPB)"),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(
  file.path(ruta_graficos_chr, "grafico_sensibilidad_prostata.png"),
  plot   = figura_sensibilidad,
  width  = 14, height = 6.5, dpi = 300, bg = "white"
)

# Version PDF vectorial, util para insertar en el documento sin perdida.
ggsave(
  file.path(ruta_graficos_chr, "grafico_sensibilidad_prostata.pdf"),
  plot   = figura_sensibilidad,
  width  = 14, height = 6.5, bg = "white"
)

message("Grafico de sensibilidad guardado en: ", ruta_graficos_chr)

print(figura_sensibilidad)
