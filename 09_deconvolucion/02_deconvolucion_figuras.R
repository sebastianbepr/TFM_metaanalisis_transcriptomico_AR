# ______________________________________________________________________________
# 09_02  Consolidacion de la deconvolucion y figuras 9 y 10
#
# Que hace: fusiona los resultados de CIBERSORTx (LM22) con las salidas de
#           immunedeconv (script anterior) en una tabla maestra y genera la
#           Figura 9 (infiltrado inmune) y la Figura 10 (contraccion del
#           compartimento no maligno / dilucion del microambiente).
# Entra:    Deconvolucion/Resultados/deconvolucion_fracciones_largo.csv,
#           Deconvolucion/CIBERSORT_Output/*.csv
# Sale:     deconvolucion_master_todos_metodos.csv, Figura9_*, Figura10_*
# ______________________________________________________________________________

# __ INDEX ___
# 01. Setup, librerias, paleta y tema
# 02. Rutas
# 03. Carga de la tabla tidy de immunedeconv
# 04. Ingesta robusta de los CSV de CIBERSORTx
# 05. Tabla maestra y exportacion
# 06. Figura 9 - infiltrado inmune (LM22)
# 07. Figura 10 - contraccion del compartimento estromal
# 08. Exportacion de figuras


# 01. SETUP ____________________________________________________________________

library(tidyverse)
library(patchwork)

col_org <- c("Próstata" = "#4A7C8C", "Testículo" = "#8C6A9E", "Vejiga" = "#E8A33D")

tema_tfm <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.title.position = "plot",
        plot.subtitle = element_text(color = "grey35", size = 9),
        legend.position = "right")


# 02. RUTAS ____________________________________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_deconv_chr    <- file.path(ruta_base_chr, "Deconvolucion")
ruta_cx_chr        <- file.path(ruta_deconv_chr, "CIBERSORT_Output")
ruta_result_chr    <- file.path(ruta_deconv_chr, "Resultados")
ruta_graficos_chr  <- file.path(ruta_base_chr, "Meta_Analisis_Graficos")
if (!dir.exists(ruta_graficos_chr)) dir.create(ruta_graficos_chr, recursive = TRUE)

traducir_organo <- function(x) recode(x, "Prostate" = "Próstata",
                                      "Testis" = "Testículo", "Bladder" = "Vejiga")


# 03. CARGA DE immunedeconv ____________________________________________________

imm_tbl <- read_csv(file.path(ruta_result_chr, "deconvolucion_fracciones_largo.csv"),
                    show_col_types = FALSE)

# Mapa (dataset, muestra) -> organo + condicion, para etiquetar CIBERSORT.
meta_muestra_tbl <- imm_tbl %>%
  distinct(dataset, muestra, organo, condicion)


# 04. INGESTA ROBUSTA DE CIBERSORTx ____________________________________________

# Los archivos se llaman <ID>_CIBERSORTx_JobN_Results.csv. El ID puede contener
# guiones bajos (TCGA_GTEX_Testis), asi que se recorta desde "_CIBERSORTx".
cols_meta_cx <- c("P-value", "P.value", "Correlation", "RMSE",
                  "Absolute score (sig.score)")

archivos_cx <- list.files(ruta_cx_chr, pattern = "\\.csv$", full.names = TRUE)

cx_tbl <- archivos_cx %>%
  map_dfr(function(f) {
    id_chr <- str_remove(basename(f), "_CIBERSORTx.*$")
    read_csv(f, show_col_types = FALSE) %>%
      rename(muestra = 1) %>%
      select(-any_of(cols_meta_cx)) %>%
      pivot_longer(-muestra, names_to = "cell_type", values_to = "valor") %>%
      mutate(dataset = id_chr, metodo = "cibersort_lm22", escala = "fraccion")
  }) %>%
  left_join(meta_muestra_tbl, by = c("dataset", "muestra"))

message("CIBERSORT ingerido: ", n_distinct(cx_tbl$dataset), " datasets, ",
        n_distinct(cx_tbl$muestra), " muestras")


# 05. TABLA MAESTRA ____________________________________________________________

master_tbl <- bind_rows(imm_tbl, cx_tbl) %>%
  mutate(organo_es = traducir_organo(organo))

write_csv(master_tbl,
          file.path(ruta_result_chr, "deconvolucion_master_todos_metodos.csv"))

# Resumen tumor vs normal (media por organo/metodo/tipo/condicion).
resumen_tn_tbl <- master_tbl %>%
  group_by(organo_es, metodo, escala, cell_type, condicion) %>%
  summarise(media = mean(valor, na.rm = TRUE), n = n(), .groups = "drop")
write_csv(resumen_tn_tbl,
          file.path(ruta_result_chr, "deconvolucion_resumen_tumor_normal.csv"))


# 06. FIGURA 9 - INFILTRADO INMUNE (LM22) ______________________________________

linaje_chr <- c(
  "B cells naive"="Linfocitos B/plasma","B cells memory"="Linfocitos B/plasma","Plasma cells"="Linfocitos B/plasma",
  "T cells CD8"="T CD8",
  "T cells CD4 naive"="T CD4/otras","T cells CD4 memory resting"="T CD4/otras","T cells CD4 memory activated"="T CD4/otras",
  "T cells follicular helper"="T CD4/otras","T cells regulatory (Tregs)"="T CD4/otras","T cells gamma delta"="T CD4/otras",
  "NK cells resting"="NK","NK cells activated"="NK",
  "Monocytes"="Monocitos/macrófagos","Macrophages M0"="Monocitos/macrófagos","Macrophages M1"="Monocitos/macrófagos","Macrophages M2"="Monocitos/macrófagos",
  "Dendritic cells resting"="Dendríticas","Dendritic cells activated"="Dendríticas",
  "Mast cells resting"="Mastocitos","Mast cells activated"="Mastocitos",
  "Eosinophils"="Eosinófilos","Neutrophils"="Neutrófilos")

linaje_niveles <- c("Linfocitos B/plasma","T CD8","T CD4/otras","NK",
                    "Monocitos/macrófagos","Dendríticas","Mastocitos",
                    "Eosinófilos","Neutrófilos")
linaje_colores <- setNames(c("#4A7C8C","#6E8FA6","#B4C4CE","#8C6A9E","#C1666B",
                             "#E8A33D","#F0D080","#B4B2A9","#7A8B6F"), linaje_niveles)

fig9_tbl <- master_tbl %>%
  filter(metodo == "cibersort_lm22") %>%
  mutate(linaje = factor(linaje_chr[cell_type], levels = linaje_niveles),
         organo_es = factor(organo_es, levels = names(col_org)),
         condicion = factor(condicion, levels = c("Normal", "Tumor"))) %>%
  group_by(organo_es, condicion, muestra, linaje) %>%
  summarise(valor = sum(valor), .groups = "drop") %>%
  group_by(organo_es, condicion, linaje) %>%
  summarise(valor = mean(valor), .groups = "drop") %>%
  group_by(organo_es, condicion) %>%
  mutate(frac = valor / sum(valor)) %>%
  ungroup()

fig9 <- fig9_tbl %>%
  ggplot(aes(x = condicion, y = frac, fill = linaje)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.3) +
  facet_wrap(~ organo_es, nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values = linaje_colores, name = "Linaje inmune") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Figura 9. Composición del infiltrado inmune (CIBERSORT · LM22)",
       subtitle = "Fracciones relativas medias por órgano y condición; las 22 poblaciones LM22 agrupadas por linaje.",
       x = NULL, y = "Composición relativa del infiltrado inmune") +
  tema_tfm +
  theme(strip.text = element_text(face = "bold", size = 11),
        strip.placement = "outside", panel.spacing = unit(1, "lines"))


# 07. FIGURA 10 - CONTRACCION DEL COMPARTIMENTO ESTROMAL ________________________

# Panel A: score de estroma de xCell (todos los organos), tumor vs normal.
figA_tbl <- master_tbl %>%
  filter(metodo == "xcell", cell_type == "stroma score") %>%
  mutate(organo_es = factor(organo_es, levels = names(col_org)),
         condicion = factor(condicion, levels = c("Normal", "Tumor")))

panel_a <- figA_tbl %>%
  ggplot(aes(x = condicion, y = valor, fill = condicion, color = organo_es)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.7) +
  facet_wrap(~ organo_es, nrow = 1, strip.position = "bottom") +
  scale_color_manual(values = col_org, guide = "none") +
  scale_fill_manual(values = c("Normal" = "white", "Tumor" = "grey75"),
                    name = NULL) +
  labs(title = "A · Contracción del estroma (xCell)",
       x = NULL, y = "xCell · stroma score") +
  tema_tfm +
  theme(strip.text = element_text(face = "bold"), legend.position = "none")

# Panel B: compartimentos EPIC (fraccion), solo RNA-seq (prostata y testiculo).
mapa_epic <- c("uncharacterized cell" = "Contenido tumoral",
               "Endothelial cell" = "Endotelio",
               "Cancer associated fibroblast" = "Fibroblastos (CAF)")

figB_tbl <- master_tbl %>%
  filter(metodo == "epic", cell_type %in% names(mapa_epic),
         organo_es %in% c("Próstata", "Testículo")) %>%
  mutate(compartimento = factor(mapa_epic[cell_type],
                                levels = c("Contenido tumoral", "Endotelio", "Fibroblastos (CAF)")),
         condicion = factor(condicion, levels = c("Normal", "Tumor")),
         organo_es = factor(organo_es, levels = c("Próstata", "Testículo"))) %>%
  group_by(organo_es, compartimento, condicion) %>%
  summarise(media = mean(valor, na.rm = TRUE), .groups = "drop")

panel_b <- figB_tbl %>%
  ggplot(aes(x = compartimento, y = media, fill = condicion)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68,
           color = "grey40", linewidth = 0.3) +
  facet_wrap(~ organo_es, nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values = c("Normal" = "white", "Tumor" = "grey60"),
                    name = NULL) +
  labs(title = "B · Compartimentos EPIC (solo RNA-seq)",
       x = NULL, y = "EPIC · fracción celular") +
  tema_tfm +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(size = 8, angle = 20, hjust = 1))

fig10 <- panel_a + panel_b +
  plot_layout(widths = c(1, 1.15)) +
  plot_annotation(
    title = "Figura 10. Contracción del compartimento no maligno en el tumor",
    subtitle = paste0("El score de estroma (xCell) cae en tumor en los tres ",
                      "órganos y el contenido tumoral (EPIC) se expande: patrón ",
                      "de dilución del microambiente."),
    caption = paste0("Relevante para la hipótesis del AR: una caída de AR en ",
                     "bulk testicular es compatible con dilución del estroma AR+ ",
                     "por expansión de células germinales tumorales AR-negativas, ",
                     "no necesariamente con ausencia de AR estromal."),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(color = "grey35", size = 9.5),
                  plot.caption = element_text(color = "grey55", size = 7.5, hjust = 0)))


# 08. EXPORTACION DE FIGURAS ___________________________________________________

ggsave(file.path(ruta_graficos_chr, "Figura9_microambiente_inmune_CIBERSORT.png"),
       fig9, width = 9.2, height = 6.2, dpi = 300, bg = "white")
ggsave(file.path(ruta_graficos_chr, "Figura9_microambiente_inmune_CIBERSORT.pdf"),
       fig9, width = 9.2, height = 6.2, bg = "white")

ggsave(file.path(ruta_graficos_chr, "Figura10_compartimento_estromal.png"),
       fig10, width = 11, height = 5.8, dpi = 300, bg = "white")
ggsave(file.path(ruta_graficos_chr, "Figura10_compartimento_estromal.pdf"),
       fig10, width = 11, height = 5.8, bg = "white")

message("Consolidacion y figuras listas.")
message("  Tabla maestra: ", file.path(ruta_result_chr, "deconvolucion_master_todos_metodos.csv"))
message("  Figuras 9 y 10 en: ", ruta_graficos_chr)
