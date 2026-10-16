# ______________________________________________________________________________
# 05_05  Diagrama de flujo del pipeline de meta-analisis
#
# Que hace: dibuja un diagrama tipo workflow (Graphviz DOT) que resume las etapas
#           del pipeline, desde los 11 datasets crudos hasta el contraste entre
#           grupos tumorales, pensado para la seccion de Metodologia.
# Entra:    (ninguno; el diagrama se define en Graphviz DOT)
# Sale:     Meta_Analisis_Graficos/diagrama_flujo_metaanalisis.png/.pdf
# ______________________________________________________________________________


# 01. SETUP ____________________________________________________________________

library(DiagrammeR)      # diagramas Graphviz
library(DiagrammeRsvg)   # exportar a SVG
library(rsvg)            # convertir SVG a PNG/PDF

# Si falta algun paquete:
#install.packages(c("DiagrammeR", "DiagrammeRsvg", "rsvg"))


# 02. RUTAS ____________________________________________________________________

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
setwd(ruta)
ruta_base_chr <- ruta

ruta_graficos_chr <- file.path(ruta_base_chr, "Meta_Analisis_Graficos")

if (!dir.exists(ruta_graficos_chr)) {
  dir.create(ruta_graficos_chr, recursive = TRUE)
}


# 03. DEFINICION DEL DIAGRAMA __________________________________________________

# El diagrama se define en lenguaje DOT (Graphviz). Cada bloque es una etapa del
# pipeline; los colores agrupan las etapas por tipo (entrada, control de calidad,
# meta-analisis, robustez, salida). Las cifras clave estan embebidas para que el
# diagrama sea auto-explicativo en el documento.

diagrama_dot_chr <- '
digraph meta_analisis {

  graph [rankdir = TB, fontname = "Helvetica", nodesep = 0.4, ranksep = 0.55,
         bgcolor = "white"]
  node  [fontname = "Helvetica", fontsize = 11, shape = box, style = "filled,rounded",
         penwidth = 0.8, margin = "0.18,0.10"]
  edge  [fontname = "Helvetica", fontsize = 9, color = "#666666", penwidth = 0.9,
         arrowsize = 0.7]

  # ---- Entrada ----
  entrada [label = "11 datasets estandarizados\\n(GEO + TCGA/GTEx)",
           fillcolor = "#E8EEF4", color = "#4A7C8C"]

  # ---- Etapa A: auditoria ----
  auditoria [label = "Etapa A — Auditoria de compatibilidad\\n8 convenciones de anotacion · 5 problemas detectados",
             fillcolor = "#E8EEF4", color = "#4A7C8C"]

  # ---- Etapa B: remediacion ----
  remediacion [label = "Etapa B — Remediacion de matrices\\nReindexacion a simbolo HGNC · escala log · filtro por varianza",
               fillcolor = "#EAF0E4", color = "#5E8C4A"]

  # ---- Agrupacion por organo ----
  grupos [label = "Agrupacion por organo",
          fillcolor = "#F1EFE8", color = "#888780", shape = box]

  prostata [label = "Prostata\\n6 datasets\\n(HPB verificado)",
            fillcolor = "#F5EDE9", color = "#C1666B"]
  testiculo [label = "Testiculo\\n3 datasets\\n(TCGA domina: -75%)",
             fillcolor = "#F5EDE9", color = "#C1666B"]
  vejiga [label = "Vejiga\\n2 datasets\\n(control no hormonal)",
          fillcolor = "#F5EDE9", color = "#C1666B"]

  # ---- Etapa C1: meta-analisis DExMA ----
  c1 [label = "Etapa C1 — Meta-analisis (DExMA)\\n3 brazos de genes faltantes x 3 metodos de combinacion\\ncomun / proporcion / imputacion  x  Stouffer / Fisher / maxP",
      fillcolor = "#EDE9F5", color = "#7A6BA8"]

  # ---- Etapa C2: robustez ----
  c2 [label = "Etapa C2 — Robustez\\nStouffer con signo · 3 ponderaciones · leave-one-out\\nConcordancia de signo (Marot)",
      fillcolor = "#EDE9F5", color = "#7A6BA8"]

  # ---- Sensibilidad ----
  sensibilidad [label = "Analisis de sensibilidad\\nGSE134073 (HPB): firma global x2, AR estable",
                fillcolor = "#FAF4E4", color = "#BA8C2A"]

  # ---- Salida ----
  contraste [label = "Contraste entre grupos\\nInterseccion Prostata n Testiculo · divergencia en Vejiga",
             fillcolor = "#E8EEF4", color = "#4A7C8C", penwidth = 1.4]

  # ---- Conexiones ----
  entrada     -> auditoria
  auditoria   -> remediacion
  remediacion -> grupos

  grupos -> prostata
  grupos -> testiculo
  grupos -> vejiga

  prostata  -> c1
  testiculo -> c1
  vejiga    -> c1

  c1 -> c2
  c2 -> sensibilidad
  sensibilidad -> contraste

  # Alinear los tres grupos en la misma fila.
  { rank = same; prostata; testiculo; vejiga }
}
'


# 04. RENDERIZADO Y EXPORTACION ________________________________________________

# Se renderiza el diagrama y se convierte a PNG de alta resolucion y a PDF
# vectorial para el documento.

diagrama <- grViz(diagrama_dot_chr)

# Exportar a SVG intermedio.
svg_chr <- export_svg(diagrama)

# PNG de alta resolucion (300 dpi aprox mediante factor de escala).
rsvg_png(
  charToRaw(svg_chr),
  file.path(ruta_graficos_chr, "diagrama_flujo_metaanalisis.png"),
  width = 1400
)

# PDF vectorial.
rsvg_pdf(
  charToRaw(svg_chr),
  file.path(ruta_graficos_chr, "diagrama_flujo_metaanalisis.pdf")
)

message("Diagrama de flujo guardado en: ", ruta_graficos_chr)

# Para ver el diagrama en el visor de RStudio:
diagrama
