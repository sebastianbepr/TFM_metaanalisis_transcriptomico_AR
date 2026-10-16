# ______________________________________________________________________________
# 10_01  Asociacion de AR con el compartimento de celulas de Leydig (testiculo)
#
# Que hace: comprueba si el AR del testiculo reside en las celulas de Leydig,
#           combinando coexpresion con marcadores de linaje (Spearman por dataset),
#           la correlacion AR-estroma dentro de cada condicion y un test de
#           dilucion con un score germinal (POU5F1/KIT). Produce la Figura 11.
# Entra:    Deconvolucion/CIBERSORT_input/*_mixture.txt,
#           Deconvolucion/Resultados/deconvolucion_master_todos_metodos.csv
# Sale:     AR_coexpresion_marcadores.csv, AR_vs_estroma_por_condicion.csv,
#           AR_test_dilucion_v2.csv, Figura11_AR_Leydig_testiculo.*
# ______________________________________________________________________________

# __ INDEX ___
# 01. Setup y marcadores
# 02. Rutas y carga (expresion de testiculo + fracciones)
# 03. Analisis 1: coexpresion AR-marcador
# 04. Analisis 2: AR vs estroma por condicion
# 05. Analisis 3: test de dilucion (score germinal)
# 06. Figura 11 y exportacion


# 01. SETUP ____________________________________________________________________

library(tidyverse)
library(patchwork)

col_clase <- c("Leydig" = "#2F6B7C", "Peritubular" = "#6E9CAE",
               "Endotelio" = "#A9C7D1", "Germinal/tumor" = "#C1666B")
col_cond  <- c("Normal" = "#8C6A9E", "Tumor" = "#4A7C8C")

tema_tfm <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 10.5),
        plot.title.position = "plot", legend.position = "right")

marcadores_tbl <- tribble(
  ~gen, ~clase,
  "CYP17A1","Leydig", "STAR","Leydig", "HSD3B1","Leydig", "INSL3","Leydig", "CYP11A1","Leydig",
  "ACTA2","Peritubular", "MYH11","Peritubular", "PDGFRB","Peritubular", "DCN","Peritubular",
  "PECAM1","Endotelio", "VWF","Endotelio",
  "POU5F1","Germinal/tumor", "NANOG","Germinal/tumor", "KIT","Germinal/tumor",
  "SALL4","Germinal/tumor", "DDX4","Germinal/tumor")
leydig_chr <- marcadores_tbl$gen[marcadores_tbl$clase == "Leydig"]
germ_chr   <- c("POU5F1", "KIT")   # marcadores especificos de tumor germinal
genes_chr  <- c("AR", marcadores_tbl$gen)
datasets_chr <- c("GSE3218", "GSE8607", "TCGA_GTEX_Testis")


# 02. RUTAS Y CARGA ____________________________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_mix_chr    <- file.path(ruta_base_chr, "Deconvolucion", "CIBERSORT_input")
ruta_result_chr <- file.path(ruta_base_chr, "Deconvolucion", "Resultados")
ruta_graf_chr   <- file.path(ruta_base_chr, "Meta_Analisis_Graficos")

# Expresion (log2) de AR y marcadores por muestra.
leer_genes <- function(gse) {
  read_tsv(file.path(ruta_mix_chr, paste0(gse, "_mixture.txt")), show_col_types = FALSE) %>%
    rename(gen = 1) %>% filter(gen %in% genes_chr) %>%
    pivot_longer(-gen, names_to = "muestra", values_to = "expr") %>%
    mutate(dataset = gse, l2 = log2(expr + 1))
}
expr_long <- datasets_chr %>% map_dfr(leer_genes)
expr_wide <- expr_long %>%
  select(dataset, muestra, gen, l2) %>%
  pivot_wider(names_from = gen, values_from = l2)

# Fracciones + condicion (tabla maestra del script de consolidacion).
master_tbl <- read_csv(file.path(ruta_result_chr, "deconvolucion_master_todos_metodos.csv"),
                       show_col_types = FALSE)
frac_testis <- master_tbl %>% filter(organo == "Testis")
stroma_tbl <- frac_testis %>% filter(metodo == "xcell", cell_type == "stroma score") %>%
  transmute(dataset, muestra, stroma = valor)
tumorc_tbl <- frac_testis %>% filter(metodo == "epic", cell_type == "uncharacterized cell") %>%
  transmute(dataset, muestra, contenido_tumoral = valor)
cond_tbl <- frac_testis %>% distinct(dataset, muestra, condicion)

datos_tbl <- expr_wide %>%
  left_join(cond_tbl, by = c("dataset", "muestra")) %>%
  left_join(stroma_tbl, by = c("dataset", "muestra")) %>%
  left_join(tumorc_tbl, by = c("dataset", "muestra"))


# 03. ANALISIS 1: COEXPRESION AR-MARCADOR ______________________________________

rho_seguro <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 8 || sd(a[ok]) == 0 || sd(b[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(a[ok], b[ok], method = "spearman"))
}

coexpr_tbl <- datos_tbl %>%
  pivot_longer(all_of(marcadores_tbl$gen), names_to = "gen", values_to = "expr_m") %>%
  group_by(dataset, gen) %>%
  summarise(rho = rho_seguro(AR, expr_m), .groups = "drop") %>%
  left_join(marcadores_tbl, by = "gen")
write_csv(coexpr_tbl, file.path(ruta_result_chr, "AR_coexpresion_marcadores.csv"))
print(coexpr_tbl %>% filter(dataset == "TCGA_GTEX_Testis") %>% arrange(desc(rho)))


# 04. ANALISIS 2: AR vs ESTROMA POR CONDICION __________________________________

ar_estroma_cond_tbl <- datos_tbl %>%
  group_by(dataset, condicion) %>%
  summarise(rho = rho_seguro(AR, stroma), n = n(), .groups = "drop")
write_csv(ar_estroma_cond_tbl, file.path(ruta_result_chr, "AR_vs_estroma_por_condicion.csv"))
print(ar_estroma_cond_tbl)


# 05. ANALISIS 3: TEST DE DILUCION _____________________________________________

t_tbl <- datos_tbl %>%
  filter(dataset == "TCGA_GTEX_Testis", is.finite(contenido_tumoral)) %>%
  mutate(AR_z = as.numeric(scale(AR)),
         germ_score = rowMeans(scale(across(all_of(germ_chr)))),
         condicion = factor(condicion, levels = c("Normal", "Tumor")))

m0 <- lm(AR_z ~ condicion, data = t_tbl)
m1 <- lm(AR_z ~ condicion + germ_score, data = t_tbl)
m2 <- lm(AR_z ~ condicion + contenido_tumoral, data = t_tbl)
b <- c(m0 = coef(m0)["condicionTumor"], m1 = coef(m1)["condicionTumor"], m2 = coef(m2)["condicionTumor"])
message("beta condicion  crudo=", round(b[1],3),
        "  +germ=", round(b[2],3), " (atenua ", round(100*(1-abs(b[2])/abs(b[1])),0), "%)",
        "  +tumorEPIC=", round(b[3],3), " (atenua ", round(100*(1-abs(b[3])/abs(b[1])),0), "%)")

dilucion_tbl <- tibble(
  modelo = factor(c("Crudo", "+ germ score", "+ contenido tumoral"),
                  levels = c("Crudo", "+ germ score", "+ contenido tumoral")),
  beta = c(coef(m0)["condicionTumor"], coef(m1)["condicionTumor"], coef(m2)["condicionTumor"]),
  se   = c(summary(m0)$coefficients["condicionTumor","Std. Error"],
           summary(m1)$coefficients["condicionTumor","Std. Error"],
           summary(m2)$coefficients["condicionTumor","Std. Error"]))
write_csv(dilucion_tbl, file.path(ruta_result_chr, "AR_test_dilucion_v2.csv"))


# 06. FIGURA 11 ________________________________________________________________

# Panel A: AR vs score de Leydig (TCGA/GTEX).
tcga_tbl <- t_tbl %>%
  mutate(leydig_score = rowMeans(scale(across(all_of(leydig_chr)))))
rho_A <- rho_seguro(tcga_tbl$AR, tcga_tbl$leydig_score)

panel_a <- tcga_tbl %>%
  ggplot(aes(leydig_score, AR, color = condicion)) +
  geom_point(size = 1.5, alpha = 0.55) +
  geom_smooth(method = "lm", se = FALSE, color = "grey20", linewidth = 0.7) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
           label = sprintf("ρ = %+.2f", rho_A), fontface = "bold",
           color = "#2F6B7C", size = 4) +
  scale_color_manual(values = col_cond, name = NULL) +
  labs(title = "A · AR sigue al programa de Leydig (TCGA/GTEX)",
       x = "Score de Leydig (z)", y = "AR (log2)") +
  tema_tfm

# Panel B: fingerprint de coexpresion (TCGA/GTEX).
panel_b <- coexpr_tbl %>%
  filter(dataset == "TCGA_GTEX_Testis") %>%
  mutate(gen = fct_reorder(gen, rho)) %>%
  ggplot(aes(rho, gen, fill = clase)) +
  geom_col() + geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
  scale_fill_manual(values = col_clase, name = "Marcador") +
  labs(title = "B · Coexpresión AR–marcador (TCGA/GTEX)",
       x = "ρ de Spearman con AR", y = NULL) +
  tema_tfm + theme(axis.text.y = element_text(size = 8))

# Panel C: test de dilucion.
panel_c <- dilucion_tbl %>%
  ggplot(aes(modelo, beta, fill = modelo)) +
  geom_col(width = 0.6, color = "grey30", linewidth = 0.3) +
  geom_errorbar(aes(ymin = beta - se, ymax = beta + se), width = 0.15) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%+.2f", beta)), vjust = 1.4, size = 3) +
  scale_fill_manual(values = c("#8C6A9E", "#4A7C8C", "#B4B2A9"), guide = "none") +
  labs(title = "C · Test de dilución (germ score ≈ 49 %)",
       x = NULL, y = "β condición (Tumor) sobre AR") +
  tema_tfm + theme(axis.text.x = element_text(size = 8))

figura11 <- panel_a + panel_b + panel_c +
  plot_layout(widths = c(1, 1.05, 0.9)) +
  plot_annotation(
    title = "Figura 11. El AR reside en el compartimento de células de Leydig del testículo",
    subtitle = paste0("El AR coexpresa con el programa esteroidogénico de Leydig ",
                      "y anticorrelaciona con marcadores germinales; ~la mitad de ",
                      "su descenso en tumor se explica por la expansión germinal (POU5F1/KIT)."),
    caption = paste0("Spearman en TCGA/GTEX (n=411). Test de dilución: modelo ",
                     "lineal AR~condición ± covariable. Evidencia in silico de ",
                     "bulk: apoya, no demuestra, la localización en Leydig."),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(color = "grey35", size = 9),
                  plot.caption = element_text(color = "grey55", size = 7.5, hjust = 0)))

ggsave(file.path(ruta_graf_chr, "Figura11_AR_Leydig_testiculo.png"),
       figura11, width = 13.5, height = 5.8, dpi = 300, bg = "white")
ggsave(file.path(ruta_graf_chr, "Figura11_AR_Leydig_testiculo.pdf"),
       figura11, width = 13.5, height = 5.8, bg = "white")
message("Figura 11 (v2) y tablas exportadas.")
