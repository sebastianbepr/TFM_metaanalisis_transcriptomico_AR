# ______________________________________________________________________________
# 06_04  Figura 8 - QC de heterogeneidad (I2) y leave-one-out
#
# Que hace: compone la figura de control de calidad del meta-analisis: panel A con
#           la distribucion de I2 por grupo tumoral y panel B con la robustez
#           leave-one-dataset-out.
# Entra:    Meta_Analisis_Resultados/heterogeneidad_por_grupo.csv, leave_one_out.csv
# Sale:     Meta_Analisis_Graficos/Figura8_QC_heterogeneidad.png/.pdf
# ______________________________________________________________________________


# 01. SETUP ____________________________________________________________________

library(tidyverse)
library(patchwork)

col_prost <- "#4A7C8C"   # prostata (teal)
col_test  <- "#8C6A9E"   # testiculo (ciruela)
col_vej   <- "#E8A33D"   # vejiga (ambar)
col_drop  <- "#C1666B"   # el dataset aporta senal (al quitarlo bajan los sig.)
col_gain  <- "#4A7C8C"   # el dataset resta (al quitarlo suben los sig.)

colores_organo <- c("Próstata" = col_prost, "Testículo" = col_test,
                    "Vejiga" = col_vej)

tema_tfm <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 12),
    plot.title.position = "plot",
    plot.subtitle    = element_text(color = "grey35", size = 9),
    axis.text        = element_text(size = 9.5),
    strip.text       = element_text(face = "bold", size = 10, hjust = 0),
    legend.position  = "none"
  )


# 02. RUTAS Y CARGA ____________________________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_salida_chr   <- file.path(ruta_base_chr, "Meta_Analisis_Resultados")
ruta_graficos_chr <- file.path(ruta_base_chr, "Meta_Analisis_Graficos")
if (!dir.exists(ruta_graficos_chr)) dir.create(ruta_graficos_chr, recursive = TRUE)

het_tbl <- read_csv(file.path(ruta_salida_chr, "heterogeneidad_por_grupo.csv"),
                    show_col_types = FALSE)
loo_tbl <- read_csv(file.path(ruta_salida_chr, "leave_one_out.csv"),
                    show_col_types = FALSE)

# Traduccion de nombres de organo.
traducir_organo <- function(x) recode(x, "Prostate" = "Próstata",
                                      "Testis" = "Testículo", "Bladder" = "Vejiga")


# 03. PANEL A: HETEROGENEIDAD I2 _______________________________________________

# Cada grupo: punto en la mediana, barra gruesa mediana->Q75, linea fina hasta Q95.
# Se anota el % de genes con heterogeneidad muy alta (I2 > 0.75).
hetA_tbl <- het_tbl %>%
  mutate(organo = factor(traducir_organo(organo),
                         levels = c("Vejiga", "Testículo", "Próstata")))

panel_a <- hetA_tbl %>%
  ggplot(aes(y = organo, color = organo)) +
  geom_vline(xintercept = c(0.50, 0.75), linetype = "dashed",
             color = "grey70", linewidth = 0.35) +
  annotate("text", x = 0.50, y = 3.65, label = "alta", size = 2.7,
           color = "grey55", hjust = -0.1) +
  annotate("text", x = 0.75, y = 3.65, label = "muy alta", size = 2.7,
           color = "grey55", hjust = -0.1) +
  geom_segment(aes(x = i2_mediana, xend = i2_q95, yend = organo),
               linewidth = 0.8, alpha = 0.55) +
  geom_segment(aes(x = i2_mediana, xend = i2_q75, yend = organo),
               linewidth = 3.2, alpha = 0.85) +
  geom_point(aes(x = i2_mediana), size = 4.2) +
  geom_text(aes(x = i2_mediana, label = sprintf("%.2f", i2_mediana)),
            vjust = -1.4, size = 3, color = "grey25") +
  geom_text(aes(x = i2_q95, label = paste0(pct_muy_alta, "% con I²>0.75")),
            hjust = -0.12, size = 3, color = "grey40") +
  scale_color_manual(values = colores_organo) +
  scale_x_continuous(limits = c(0, 1.28), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.01, 0))) +
  labs(title = "A · Heterogeneidad entre estudios (I²)",
       subtitle = "Punto = mediana · barra gruesa = mediana→Q75 · línea = →Q95",
       x = "I² por gen", y = NULL) +
  tema_tfm


# 04. PANEL B: LEAVE-ONE-DATASET-OUT ___________________________________________

# Cambio relativo (%) en genes significativos al excluir cada dataset. Valores
# muy negativos = el dataset aporta la mayor parte de la senal del grupo.
loB_tbl <- loo_tbl %>%
  mutate(
    organo  = factor(traducir_organo(organo), levels = c("Próstata", "Testículo")),
    signo   = if_else(cambio_relativo_pct < 0, "aporta", "resta"),
    dataset = fct_reorder(dataset_excluido, cambio_relativo_pct)
  )

panel_b <- loB_tbl %>%
  ggplot(aes(x = cambio_relativo_pct, y = dataset, fill = signo)) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%+.1f%%", cambio_relativo_pct),
                hjust = if_else(cambio_relativo_pct < 0, 1.1, -0.1)),
            size = 3, color = "grey25") +
  facet_grid(organo ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c("aporta" = col_drop, "resta" = col_gain)) +
  scale_x_continuous(limits = c(-90, 45), breaks = seq(-75, 25, 25),
                     labels = function(x) paste0(x, "%")) +
  labs(title = "B · Robustez leave-one-dataset-out",
       subtitle = "Cambio en genes significativos al excluir cada estudio (rojo = el estudio aporta señal)",
       x = "Cambio relativo en genes significativos", y = NULL) +
  tema_tfm +
  theme(strip.placement = "outside")


# 05. COMPOSICION Y EXPORTACION ________________________________________________

figura8 <- panel_a / panel_b +
  plot_layout(heights = c(1, 1.35)) +
  plot_annotation(
    title = "Figura 8. Control de calidad del meta-análisis integrativo",
    subtitle = paste0("La alta heterogeneidad del grupo testicular y el dominio ",
                      "de TCGA/GTEX se controlan con Stouffer firmado ponderado ",
                      "por √n_efectivo y concordancia de signo ≥0.75."),
    caption = paste0("I²: convención de Higgins (>0.50 alta, >0.75 muy alta). ",
                     "Leave-one-out sobre el meta-análisis firmado por grupo. ",
                     "Vejiga (2 datasets) no admite leave-one-out informativo."),
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(color = "grey35", size = 9.5),
                  plot.caption = element_text(color = "grey55", size = 7.5, hjust = 0))
  )

ggsave(file.path(ruta_graficos_chr, "Figura8_QC_heterogeneidad.png"),
       figura8, width = 9.5, height = 7.6, dpi = 300, bg = "white")
ggsave(file.path(ruta_graficos_chr, "Figura8_QC_heterogeneidad.pdf"),
       figura8, width = 9.5, height = 7.6, bg = "white")

message("Figura 8 (QC heterogeneidad) guardada en: ", ruta_graficos_chr)
print(figura8)
