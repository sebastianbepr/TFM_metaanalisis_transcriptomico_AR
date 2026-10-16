# ______________________________________________________________________________
# 07_03  Diagrama de interseccion de las firmas consenso (Firma 1 y Firma 2)
#
# Que hace: dibuja la figura de interseccion de las firmas, un panel Euler
#           anidado (Firma 1 = prostata/testiculo, Firma 2 anidada como donut) y
#           una barra apilada con la anatomia de la Firma 1. Los conteos se
#           derivan del contraste.
# Entra:    Meta_Analisis_Resultados/contraste_prostata_4ds.csv
# Sale:     Meta_Analisis_Graficos/Diagrama_interseccion_firmas_consenso.png/.pdf
# ______________________________________________________________________________


# 01. SETUP ____________________________________________________________________

library(tidyverse)
library(ggforce)     # geom_circle / geom_arc_bar (donut)
library(patchwork)   # composicion de paneles

# Paleta tomada de los scripts del TFM (12 Graficos resumen ...).
col_prost <- "#4A7C8C"   # prostata (teal)
col_test  <- "#8C6A9E"   # testiculo (ciruela)
col_lens  <- "#F2ECDD"   # firma 1 (pergamino)
col_inv   <- "#C1666B"   # inversion en vejiga
col_aus   <- "#E8A33D"   # ausencia en vejiga
col_ausd  <- "#F0D080"   # ausencia debil en vejiga
col_conc  <- "#B4B2A9"   # concordante (no diverge)
col_noev  <- "#D8D8D8"   # no evaluado en vejiga

tema_diag <- theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0))


# 02. RUTAS Y CARGA ____________________________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_salida_chr   <- file.path(ruta_base_chr, "Meta_Analisis_Resultados")
ruta_graficos_chr <- file.path(ruta_base_chr, "Meta_Analisis_Graficos")
if (!dir.exists(ruta_graficos_chr)) dir.create(ruta_graficos_chr, recursive = TRUE)

genes_eje_ar_chr <- c("AR", "KLK3", "KLK2", "NKX3-1", "FKBP5",
                      "TMPRSS2", "FOLH1", "STEAP1", "SLC45A3", "HOXB13")

contraste_tbl <- read_csv(
  file.path(ruta_salida_chr, "contraste_prostata_4ds.csv"),
  show_col_types = FALSE
)


# 03. CONTEOS DE FIRMA 1 Y FIRMA 2 _____________________________________________

# La Firma 1 es la totalidad del contraste (interseccion prostata y testiculo).
# La Firma 2 son las categorias divergentes en vejiga.
divergentes_chr <- c("A_inversion", "B_ausencia", "B_ausencia_debil")

conteos_tbl <- contraste_tbl %>%
  count(categoria_vejiga, name = "n")

n_de <- function(cat) {
  v <- conteos_tbl$n[conteos_tbl$categoria_vejiga == cat]
  if (length(v) == 0) 0L else v
}

n_firma1 <- nrow(contraste_tbl)
n_inv    <- n_de("A_inversion")
n_aus    <- n_de("B_ausencia")
n_ausd   <- n_de("B_ausencia_debil")
n_conc   <- n_de("Concordante (no divergente)")
n_noev   <- n_de("No evaluado en vejiga")
n_firma2 <- n_inv + n_aus + n_ausd
n_resto  <- n_conc + n_noev

genes_ar_firma2 <- contraste_tbl %>%
  filter(categoria_vejiga %in% divergentes_chr,
         gen_simbolo %in% genes_eje_ar_chr) %>%
  pull(gen_simbolo)

stopifnot(n_firma1 == n_firma2 + n_resto)
message("Firma 1 = ", n_firma1, "  Firma 2 = ", n_firma2,
        "  (inv ", n_inv, " / aus ", n_aus, " / aus.debil ", n_ausd, ")")


# 04. PANEL A: DIAGRAMA DE INTERSECCION ________________________________________

R  <- 1.62      # radio de los circulos de tejido
cx <- 0.72      # semiseparacion (mayor solapamiento -> lente mas ancho)

# Lente de interseccion (Firma 1) como poligono, para poder rellenarlo.
yint  <- sqrt(R^2 - cx^2)
ang_p <- atan2(yint, cx)
arco  <- function(cx0, a0, a1, n = 200) {
  a <- seq(a0, a1, length.out = n)
  tibble(x = cx0 + R * cos(a), y = R * sin(a))
}
lente_tbl <- bind_rows(
  arco(-cx, -ang_p, ang_p),                  # arco derecho de prostata
  arco( cx, pi - ang_p, pi + ang_p)          # arco izquierdo de testiculo
)

# Donut de la Firma 2 (por categoria) con geom_arc_bar.
f2_cx <- 0; f2_cy <- -0.30; f2_r0 <- 0.40; f2_r1 <- 0.70
donut_tbl <- tibble(
  cat = factor(c("Inversión", "Ausencia", "Ausencia débil"),
               levels = c("Inversión", "Ausencia", "Ausencia débil")),
  n   = c(n_inv, n_aus, n_ausd),
  col = c(col_inv, col_aus, col_ausd)
) %>%
  mutate(frac = n / sum(n),
         a1   = pi/2 - 2*pi*cumsum(lag(frac, default = 0)),
         a2   = pi/2 - 2*pi*cumsum(frac))

circulos_tbl <- tibble(
  x0 = c(-cx, cx), y0 = c(0, 0), r = R,
  tejido = c("Próstata", "Testículo"), col = c(col_prost, col_test)
)

panel_a <-
  ggplot() +
  # circulos de tejido
  geom_circle(data = circulos_tbl,
              aes(x0 = x0, y0 = y0, r = r, fill = tejido, color = tejido),
              alpha = 0.16, linewidth = 0.9) +
  scale_fill_manual(values = c("Próstata" = col_prost, "Testículo" = col_test),
                    guide = "none") +
  scale_color_manual(values = c("Próstata" = col_prost, "Testículo" = col_test),
                     guide = "none") +
  # lente = Firma 1
  geom_polygon(data = lente_tbl, aes(x = x, y = y),
               fill = col_lens, color = "#B9A97F", linewidth = 0.6) +
  # donut = Firma 2
  geom_arc_bar(data = donut_tbl,
               aes(x0 = f2_cx, y0 = f2_cy, r0 = f2_r0, r = f2_r1,
                   start = a1, end = a2, fill = cat),
               color = "white", linewidth = 0.5, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Inversión" = col_inv, "Ausencia" = col_aus,
                               "Ausencia débil" = col_ausd),
                    name = "Divergencia en vejiga") +
  # nucleo blanco del donut
  geom_circle(aes(x0 = f2_cx, y0 = f2_cy, r = f2_r0),
              fill = "white", color = "#C24E45", linewidth = 0.6) +
  # etiquetas
  annotate("text", x = -cx - 0.62, y = 1.45, label = "PRÓSTATA",
           fontface = "bold", size = 5, color = col_prost) +
  annotate("text", x = cx + 0.62, y = 1.45, label = "TESTÍCULO",
           fontface = "bold", size = 5, color = col_test) +
  annotate("text", x = 0, y = 1.15, label = "FIRMA 1",
           fontface = "bold", size = 4.5, color = "#6E5B2E") +
  annotate("text", x = 0, y = 0.78,
           label = paste0(n_firma1, " genes\ncompartidos"),
           size = 3.3, color = "#6E5B2E", lineheight = 0.9) +
  annotate("text", x = f2_cx, y = f2_cy + 0.10, label = "FIRMA 2",
           fontface = "bold", size = 3.8, color = "#B0433C") +
  annotate("text", x = f2_cx, y = f2_cy - 0.14, label = n_firma2,
           fontface = "bold", size = 5, color = "#B0433C") +
  annotate("text", x = 0, y = -1.70,
           label = paste0(n_resto, " restantes mantienen señal en vejiga\n(",
                          n_conc, " concordantes · ", n_noev, " no evaluados)"),
           size = 2.8, color = "#6b6b6b", lineheight = 0.9) +
  annotate("text", x = -1.85, y = -0.42,
           label = paste0("eje AR en Firma 2:\n",
                          paste(genes_ar_firma2, collapse = " · ")),
           fontface = "bold", size = 2.9, color = "#333333", lineheight = 0.95) +
  coord_fixed(xlim = c(-2.7, 2.7), ylim = c(-2.0, 1.7)) +
  labs(title = NULL) +
  tema_diag +
  theme(legend.position = "none")


# 05. PANEL B: ANATOMIA DE LA FIRMA 1 (BARRA APILADA) __________________________

barra_tbl <- tibble(
  cat = factor(c("Inversión en vejiga", "Ausencia en vejiga",
                 "Ausencia débil en vejiga", "Concordante (mantiene señal)",
                 "No evaluado en vejiga"),
               levels = rev(c("Inversión en vejiga", "Ausencia en vejiga",
                              "Ausencia débil en vejiga",
                              "Concordante (mantiene señal)",
                              "No evaluado en vejiga"))),
  n   = c(n_inv, n_aus, n_ausd, n_conc, n_noev)
)
paleta_barra <- c(
  "Inversión en vejiga" = col_inv, "Ausencia en vejiga" = col_aus,
  "Ausencia débil en vejiga" = col_ausd,
  "Concordante (mantiene señal)" = col_conc, "No evaluado en vejiga" = col_noev)

panel_b <-
  ggplot(barra_tbl, aes(x = 1, y = n, fill = cat)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(cat, " · ", n)),
            position = position_stack(vjust = 0.5), hjust = 0,
            x = 1.33, size = 3, color = "#333333") +
  scale_fill_manual(values = paleta_barra, guide = "none") +
  coord_cartesian(xlim = c(0.7, 2.6), clip = "off") +
  labs(title = "Anatomía de la Firma 1 (287 → filtro de divergencia)") +
  tema_diag


# 06. COMPOSICION Y EXPORTACION ________________________________________________

figura <- panel_a + panel_b +
  plot_layout(widths = c(1.7, 1)) +
  plot_annotation(
    title = "Diagrama de intersección de las firmas consenso",
    subtitle = paste0("Meta-análisis integrativo de tumores urológicos ",
                      "hormono-dependientes · señal compartida próstata ∩ ",
                      "testículo (Firma 1) y subconjunto hormono-específico ",
                      "(Firma 2)"),
    caption = paste0("Firma consenso por meta-análisis firmado (Stouffer con ",
                     "signo, FDR<0.05, concordancia ≥0.75). Firma 2 = genes de ",
                     "la Firma 1 con divergencia en vejiga urotelial."),
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )

ggsave(file.path(ruta_graficos_chr, "Diagrama_interseccion_firmas_consenso.png"),
       figura, width = 13.2, height = 7.6, dpi = 300, bg = "white")
ggsave(file.path(ruta_graficos_chr, "Diagrama_interseccion_firmas_consenso.pdf"),
       figura, width = 13.2, height = 7.6, bg = "white")

message("Figura guardada en: ", ruta_graficos_chr)
