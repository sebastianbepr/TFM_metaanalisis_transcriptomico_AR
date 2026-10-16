# ______________________________________________________________________________
# 08_02  Figura 5 - Terminos ORA significativos por firma
#
# Que hace: dibuja los terminos ORA significativos (FDR<0.05) como un
#           lollipop-dotplot faceteado por firma, con el eje X de Fold
#           Enrichment en escala log.
# Entra:    Meta_Analisis_Enriquecimiento/enriquecimiento_ORA.csv
# Sale:     Meta_Analisis_Graficos/Figura5_ORA_terminos_por_firma.png/.pdf
# ______________________________________________________________________________


# 01. SETUP ____________________________________________________________________

library(tidyverse)
library(scales)

# Paleta on-brand (misma familia que "12 Graficos resumen ...").
col_firma1 <- "#4A7C8C"   # teal
col_firma2 <- "#C1666B"   # terracota
grad_low   <- "#F0D080"   # menos significativo
grad_high  <- "#8A3B34"   # mas significativo

tema_tfm <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(face = "bold", size = 13),
    plot.title.position = "plot",
    plot.subtitle      = element_text(color = "grey35", size = 9.5),
    plot.caption       = element_text(color = "grey55", size = 7.5, hjust = 0),
    legend.position    = "right",
    strip.text         = element_text(face = "bold", hjust = 0, size = 10.5),
    axis.text.y        = element_text(size = 9)
  )


# 02. RUTAS Y CARGA ____________________________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_enriq_chr    <- file.path(ruta_base_chr, "Meta_Analisis_Enriquecimiento")
ruta_graficos_chr <- file.path(ruta_base_chr, "Meta_Analisis_Graficos")
if (!dir.exists(ruta_graficos_chr)) dir.create(ruta_graficos_chr, recursive = TRUE)

ora_tbl <- read_csv(file.path(ruta_enriq_chr, "enriquecimiento_ORA.csv"),
                    show_col_types = FALSE)


# 03. PREPARACION ______________________________________________________________

# Etiquetas legibles de firma. Se conserva el orden Firma1 -> Firma2.
etiquetas_firma_chr <- c(
  "Firma1_Down" = "Firma 1 · Down — señal compartida próstata ∩ testículo",
  "Firma2_Down" = "Firma 2 · Down — subconjunto hormono-específico"
)

ora_plot_tbl <- ora_tbl %>%
  mutate(
    firma        = factor(lista, levels = names(etiquetas_firma_chr),
                          labels = etiquetas_firma_chr[names(etiquetas_firma_chr)]),
    neglog10_fdr = -log10(p.adjust),
    # Genes del termino listos para tooltip / anotacion, si se quisieran.
    genes        = str_replace_all(geneID, "/", ", ")
  ) %>%
  # Orden de terminos por Fold Enrichment (mas enriquecido arriba). Como las
  # descripciones son unicas, un factor global sirve; el facet libre descarta
  # los niveles no usados en cada panel.
  arrange(FoldEnrichment) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Genes del eje androgenico dentro del ORA (para resaltar SRD5A2 en Firma2).
genes_eje_ar_chr <- c("AR", "KLK3", "KLK2", "NKX3-1", "FKBP5",
                      "TMPRSS2", "FOLH1", "STEAP1", "SLC45A3", "HOXB13", "SRD5A2")


# 04. FIGURA ___________________________________________________________________

fig_ora <- ora_plot_tbl %>%
  ggplot(aes(x = FoldEnrichment, y = Description)) +
  # Baston del lollipop desde 1 (sin enriquecimiento) hasta el valor.
  geom_segment(aes(x = 1, xend = FoldEnrichment,
                   y = Description, yend = Description),
               color = "grey78", linewidth = 0.5) +
  geom_point(aes(size = Count, color = neglog10_fdr)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey60", linewidth = 0.3) +
  # Alturas proporcionales por firma (12 vs 1 termino) y eje X compartido.
  facet_grid(firma ~ ., scales = "free_y", space = "free_y",
             labeller = label_wrap_gen(width = 45)) +
  scale_color_gradient(low = grad_low, high = grad_high,
                       name = expression(-log[10]~"(FDR)")) +
  scale_size_continuous(range = c(3, 9), name = "N.º de genes\n(Count)") +
  # Escala X logaritmica compartida: hace legibles a la vez el grupo muscular
  # (FE ~ 4-12) y el termino de benceno (FE ~ 85), sin distorsionar la posicion.
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 85),
                expand = expansion(mult = c(0.03, 0.10))) +
  labs(
    title = "Figura 5. Términos ORA significativos por firma",
    subtitle = paste0("13 términos con FDR < 0.05, todos en subconjuntos de ",
                      "infraexpresión (Down); ninguna firma Up dio términos ",
                      "significativos.\nLa Firma1_Down está dominada por el ",
                      "programa muscular/contráctil."),
    x = "Fold Enrichment (escala log)",
    y = NULL,
    caption = paste0("Sobre-representación (ORA, GO:BP) sobre las firmas ",
                     "consenso del meta-análisis firmado. Eje X logarítmico ",
                     "compartido. Firma2_Down: benzene-containing compound ",
                     "metabolic process (GSTM2, GSTM3, SRD5A2 — SRD5A2 = ",
                     "5α-reductasa).")
  ) +
  tema_tfm


# 05. EXPORTACION ______________________________________________________________

ggsave(file.path(ruta_graficos_chr, "Figura5_ORA_terminos_por_firma.png"),
       fig_ora, width = 10, height = 6.8, dpi = 300, bg = "white")
ggsave(file.path(ruta_graficos_chr, "Figura5_ORA_terminos_por_firma.pdf"),
       fig_ora, width = 10, height = 6.8, bg = "white")

message("Figura 5 (ORA por firma) guardada en: ", ruta_graficos_chr)
print(fig_ora)
