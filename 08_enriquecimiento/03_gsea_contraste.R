# ______________________________________________________________________________
# 08_03  Figura 6 - GSEA, contraste prostata vs testiculo
#
# Que hace: compara, para un panel curado de rutas, el NES y el FDR de la GSEA
#           en prostata y testiculo como un dotplot faceteado por bloque
#           funcional.
# Entra:    Meta_Analisis_Enriquecimiento/enriquecimiento_GSEA.csv
# Sale:     Meta_Analisis_Graficos/Figura6_GSEA_contraste_prostata_testiculo.png/.pdf
# ______________________________________________________________________________


# 01. SETUP ____________________________________________________________________

library(tidyverse)

col_neg <- "#4A7C8C"   # NES negativo (infraexpresado) — teal, como en el heatmap
col_pos <- "#C1666B"   # NES positivo (sobreexpresado) — terracota

tema_tfm <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    plot.title       = element_text(face = "bold", size = 13),
    plot.title.position = "plot",
    plot.subtitle    = element_text(color = "grey35", size = 9.5),
    plot.caption     = element_text(color = "grey55", size = 7.5, hjust = 0),
    strip.text.y     = element_text(face = "bold", angle = 0, hjust = 0, size = 9.5),
    strip.placement  = "outside",
    axis.text.x      = element_text(face = "bold", size = 11),
    axis.text.y      = element_text(size = 9),
    legend.position  = "right"
  )


# 02. RUTAS Y CARGA ____________________________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_enriq_chr    <- file.path(ruta_base_chr, "Meta_Analisis_Enriquecimiento")
ruta_graficos_chr <- file.path(ruta_base_chr, "Meta_Analisis_Graficos")
if (!dir.exists(ruta_graficos_chr)) dir.create(ruta_graficos_chr, recursive = TRUE)

gsea_tbl <- read_csv(file.path(ruta_enriq_chr, "enriquecimiento_GSEA.csv"),
                     show_col_types = FALSE)


# 03. PANEL CURADO DE RUTAS ____________________________________________________

# Se seleccionan por ID para que el contraste sea legible y reproducible. El orden
# (nivel del factor) fija la disposicion vertical; el bloque agrupa por facet.
panel_tbl <- tribble(
  ~ID,             ~bloque,                              ~etiqueta,                                        ~orden,
  "GO:0048545",    "Receptores nucleares / esteroideo",  "response to steroid hormone",                    1,
  "GO:0009755",    "Receptores nucleares / esteroideo",  "hormone-mediated signaling pathway",             2,
  "GO:0141193",    "Receptores nucleares / esteroideo",  "nuclear receptor-mediated signaling pathway",    3,
  "GO:0043401",    "Receptores nucleares / esteroideo",  "steroid hormone receptor signaling pathway",     4,
  "R-HSA-9006931", "Receptores nucleares / esteroideo",  "Signaling by Nuclear Receptors (Reactome)",      5,
  "hsa04935",      "Growth hormone / biosíntesis",       "GH synthesis, secretion & action (KEGG)",        6,
  "R-HSA-982772",  "Growth hormone / biosíntesis",       "Growth hormone receptor signaling (Reactome)",   7,
  "GO:0060396",    "Growth hormone / biosíntesis",       "growth hormone receptor signaling pathway",      8,
  "GO:0071378",    "Growth hormone / biosíntesis",       "cellular response to growth hormone stimulus",   9,
  "hsa00140",      "Growth hormone / biosíntesis",       "Steroid hormone biosynthesis (KEGG)",           10
)


# 04. UNION CON LOS RESULTADOS POR COHORTE _____________________________________

comp_tbl <- gsea_tbl %>%
  filter(ID %in% panel_tbl$ID) %>%
  select(ID, grupo, NES, p.adjust) %>%
  right_join(panel_tbl, by = "ID") %>%          # conserva rutas ausentes (NES = NA)
  filter(!is.na(grupo)) %>%                       # deja solo las combinaciones observadas
  mutate(
    cohorte      = recode(grupo, "Prostata" = "Próstata", "Testiculo" = "Testículo"),
    cohorte      = factor(cohorte, levels = c("Próstata", "Testículo")),
    bloque       = factor(bloque, levels = c("Receptores nucleares / esteroideo",
                                             "Growth hormone / biosíntesis")),
    etiqueta     = fct_reorder(etiqueta, -orden),
    neglog10_fdr = -log10(p.adjust),
    significativo = p.adjust < 0.05
  )


# 05. DOTPLOT COMPARATIVO ______________________________________________________

lim_nes <- max(abs(comp_tbl$NES), na.rm = TRUE)

fig6 <- comp_tbl %>%
  ggplot(aes(x = cohorte, y = etiqueta)) +
  geom_point(aes(color = NES, size = neglog10_fdr, shape = significativo,
                 fill = NES), stroke = 0.9) +
  scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 1),
                     labels = c("TRUE" = "FDR < 0.05", "FALSE" = "no signif."),
                     name = "Significancia") +
  scale_color_gradient2(low = col_neg, mid = "grey88", high = col_pos,
                        midpoint = 0, limits = c(-lim_nes, lim_nes),
                        name = "NES", aesthetics = c("colour", "fill")) +
  scale_size_continuous(range = c(3, 10), name = expression(-log[10]~"(FDR)")) +
  facet_grid(bloque ~ ., scales = "free_y", space = "free_y",
             switch = "y", labeller = label_wrap_gen(width = 20)) +
  labs(
    title = "Figura 6. Enriquecimiento GSEA: contraste próstata vs testículo",
    subtitle = paste0("El eje de receptores nucleares/esteroideo se enriquece ",
                      "de forma coordinada hacia la infraexpresión (NES < 0) ",
                      "solo en próstata;\nen testículo predomina, en cambio, la ",
                      "señal positiva de growth hormone."),
    x = NULL, y = NULL,
    caption = paste0("GSEA (clusterProfiler) sobre el ranking de Z combinado ",
                     "firmado por cohorte. Punto relleno = FDR < 0.05; anillo ",
                     "vacío = no significativo. Ausencia de punto = ruta no ",
                     "enriquecida en esa cohorte. NES < 0: extremo infraexpresado ",
                     "del ranking tumor–normal.")
  ) +
  tema_tfm

ggsave(file.path(ruta_graficos_chr, "Figura6_GSEA_contraste_prostata_testiculo.png"),
       fig6, width = 9.5, height = 6.6, dpi = 300, bg = "white")
ggsave(file.path(ruta_graficos_chr, "Figura6_GSEA_contraste_prostata_testiculo.pdf"),
       fig6, width = 9.5, height = 6.6, bg = "white")

message("Figura 6 guardada en: ", ruta_graficos_chr)
print(fig6)
