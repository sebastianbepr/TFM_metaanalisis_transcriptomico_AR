# ______________________________________________________________________________
# 07_02  Figuras resumen de la firma hormonal compartida
#
# Que hace: genera dos figuras del TFM, un grafico de dispersion prostata vs
#           testiculo coloreado por divergencia en vejiga y un heatmap de los
#           genes divergentes en los tres grupos.
# Entra:    Meta_Analisis_Resultados/contraste_prostata_4ds.csv
# Sale:     Meta_Analisis_Graficos/grafico_dispersion_firma.*,
#           grafico_heatmap_divergencia.*
# ______________________________________________________________________________


# 01. SETUP ____________________________________________________________________

library(tidyverse)
library(scales)

# Paleta por categoria de divergencia.
colores_categoria_chr <- c(
  "Inversion en vejiga"       = "#C1666B",
  "Ausencia en vejiga"        = "#E8A33D",
  "Ausencia debil en vejiga"  = "#F0D080",
  "Concordante (no diverge)"  = "#B4B2A9",
  "No evaluado en vejiga"     = "#D8D8D8"
)

tema_tfm <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 12),
    plot.title.position = "plot",
    legend.position  = "right",
    legend.title     = element_text(size = 10)
  )


# 02. RUTAS Y CARGA ____________________________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_salida_chr   <- file.path(ruta_base_chr, "Meta_Analisis_Resultados")
ruta_graficos_chr <- file.path(ruta_base_chr, "Meta_Analisis_Graficos")

if (!dir.exists(ruta_graficos_chr)) {
  dir.create(ruta_graficos_chr, recursive = TRUE)
}

genes_eje_ar_chr <- c("AR", "KLK3", "KLK2", "NKX3-1", "FKBP5",
                      "TMPRSS2", "FOLH1", "STEAP1", "SLC45A3", "HOXB13")

contraste_tbl <- read_csv(
  file.path(ruta_salida_chr, "contraste_prostata_4ds.csv"),
  show_col_types = FALSE
)

# Se mapea la categoria interna a etiqueta legible y compacta.
mapa_categoria_chr <- c(
  "A_inversion"                 = "Inversion en vejiga",
  "B_ausencia"                  = "Ausencia en vejiga",
  "B_ausencia_debil"            = "Ausencia debil en vejiga",
  "Concordante (no divergente)" = "Concordante (no diverge)",
  "No evaluado en vejiga"       = "No evaluado en vejiga"
)

contraste_tbl <- contraste_tbl %>%
  mutate(
    categoria_legible = factor(
      mapa_categoria_chr[categoria_vejiga],
      levels = names(colores_categoria_chr)
    ),
    es_ar = gen_simbolo %in% genes_eje_ar_chr
  )


# 03. GRAFICO 1: DISPERSION PROSTATA VS TESTICULO ______________________________

# Cada punto es un gen de la firma compartida. Los ejes son el logFC combinado en
# prostata y testiculo. Todos caen en los cuadrantes concordantes (ambos Up o
# ambos Down) por construccion. El color indica su destino en vejiga. Se etiquetan
# los genes del eje androgenico presentes.

etiquetas_ar_tbl <- contraste_tbl %>%
  filter(es_ar)

grafico_dispersion <- contraste_tbl %>%
  ggplot(aes(x = logfc_prostata, y = logfc_testiculo)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(
    aes(color = categoria_legible),
    size = 2, alpha = 0.75
  ) +
  # Resaltar y etiquetar los genes del eje androgenico.
  geom_point(
    data = etiquetas_ar_tbl,
    shape = 21, size = 3.2, stroke = 0.8, fill = NA, color = "black"
  ) +
  ggrepel::geom_text_repel(
    data = etiquetas_ar_tbl,
    aes(label = gen_simbolo),
    size = 3, fontface = "bold", box.padding = 0.5, max.overlaps = Inf
  ) +
  scale_color_manual(values = colores_categoria_chr,
                     name = "Divergencia\nen vejiga") +
  labs(
    title = "Firma hormonal compartida: prostata vs testiculo",
    x     = "logFC combinado en prostata",
    y     = "logFC combinado en testiculo"
  ) +
  tema_tfm


# 04. GRAFICO 2: HEATMAP DE DIVERGENCIA ________________________________________

# Se seleccionan los genes divergentes (inversion o ausencia) y se muestra su
# logFC en los tres grupos. El patron esperado: color intenso y concordante en
# prostata/testiculo, atenuado o invertido en vejiga.

genes_divergentes_tbl <- contraste_tbl %>%
  filter(categoria_vejiga %in% c("A_inversion", "B_ausencia")) %>%
  # Priorizar los de mayor efecto para que el heatmap sea legible.
  arrange(desc(abs(logfc_prostata))) %>%
  slice_head(n = 40)

# El orden de los genes se fija ANTES del pivot, mientras logfc_prostata sigue
# siendo una columna, y se conserva despues como niveles del factor.
orden_genes_chr <- genes_divergentes_tbl %>%
  arrange(logfc_prostata) %>%
  pull(gen_simbolo)

heatmap_largo_tbl <- genes_divergentes_tbl %>%
  select(gen_simbolo, logfc_prostata, logfc_testiculo, logfc_vejiga) %>%
  pivot_longer(
    cols      = starts_with("logfc_"),
    names_to  = "grupo",
    values_to = "logfc"
  ) %>%
  mutate(
    grupo = recode(grupo,
                   logfc_prostata  = "Prostata",
                   logfc_testiculo = "Testiculo",
                   logfc_vejiga    = "Vejiga"),
    grupo = factor(grupo, levels = c("Prostata", "Testiculo", "Vejiga")),
    gen_simbolo = factor(gen_simbolo, levels = orden_genes_chr)
  )

# Limite simetrico para la escala de color divergente.
limite_num <- max(abs(heatmap_largo_tbl$logfc), na.rm = TRUE)

grafico_heatmap <- heatmap_largo_tbl %>%
  ggplot(aes(x = grupo, y = gen_simbolo, fill = logfc)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(
    low = "#4A7C8C", mid = "white", high = "#C1666B",
    midpoint = 0, limits = c(-limite_num, limite_num),
    name = "logFC"
  ) +
  labs(
    title = "Divergencia en vejiga de los genes hormono-compartidos",
    x     = NULL,
    y     = NULL
  ) +
  tema_tfm +
  theme(
    axis.text.y = element_text(size = 7),
    panel.grid  = element_blank()
  )


# 05. EXPORTACION ______________________________________________________________

ggsave(
  file.path(ruta_graficos_chr, "grafico_dispersion_firma.png"),
  plot = grafico_dispersion, width = 9, height = 6.5, dpi = 300, bg = "white"
)
ggsave(
  file.path(ruta_graficos_chr, "grafico_dispersion_firma.pdf"),
  plot = grafico_dispersion, width = 9, height = 6.5, bg = "white"
)

ggsave(
  file.path(ruta_graficos_chr, "grafico_heatmap_divergencia.png"),
  plot = grafico_heatmap, width = 6, height = 8, dpi = 300, bg = "white"
)
ggsave(
  file.path(ruta_graficos_chr, "grafico_heatmap_divergencia.pdf"),
  plot = grafico_heatmap, width = 6, height = 8, bg = "white"
)

message("Graficos guardados en: ", ruta_graficos_chr)

print(grafico_dispersion)
print(grafico_heatmap)
