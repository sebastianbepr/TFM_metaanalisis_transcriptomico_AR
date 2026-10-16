# Requisitos de R

Todo el análisis se ejecutó en R (salvo CIBERSORT, que se ejecuta en línea en CIBERSORTx).
Se recomienda R >= 4.3 y Bioconductor >= 3.18. La semilla aleatoria está fijada en `1234` en
todos los scripts que lo requieren, para reproducibilidad.

## Paquetes por etapa

**01_descarga**
- CRAN: `tidyverse`, `data.table`, `here`
- Bioconductor: `GEOquery`
- TCGA/GTEx: `UCSCXenaTools` (o descarga directa del hub Toil de UCSC Xena), `data.table`

**02_estandarizacion**
- CRAN: `tidyverse`, `data.table`, `readxl`
- Bioconductor: `GEOquery`, `limma`, `WGCNA` (colapso de sondas por gen)

**03_dge**
- Bioconductor: `DESeq2`, `limma`, `edgeR`
- CRAN: `tidyverse`, `here`

**04_auditoria_remediacion**
- CRAN: `tidyverse`
- Bioconductor: `edgeR`, `AnnotationDbi`, `org.Hs.eg.db`

**05_metaanalisis**
- Bioconductor: `DExMA`
- CRAN: `tidyverse`, `R.utils` (límite de tiempo en la imputación)

**06_sensibilidad**
- Bioconductor: `DExMA`, `limma`
- CRAN: `tidyverse`, `ggplot2`

**07_firmas**
- CRAN: `tidyverse`, `ggplot2`, `writexl`, `pheatmap` / `ggplot2` (según figura)

**08_enriquecimiento**
- Bioconductor: `clusterProfiler`, `ReactomePA`, `org.Hs.eg.db`, `enrichplot`
- CRAN: `tidyverse`, `patchwork`, `ggplot2`

**09_deconvolucion**
- `immunedeconv` (MCP-counter, xCell, quanTIseq, EPIC) y sus dependencias
- CRAN: `tidyverse`, `ggplot2`
- Externo: CIBERSORTx (plataforma web), firma LM22

**10_ar_estroma**
- CRAN: `tidyverse`, `ggplot2`

## Diagramas de flujo

Los scripts que dibujan diagramas (05_metaanalisis/05, etc.) usan `DiagrammeR` y `rsvg`
para exportar el Graphviz a PNG/PDF.

## Nota

Las versiones exactas de los paquetes deben registrarse con `sessionInfo()` al final de una
ejecución completa, y adjuntarse al TFM para dejar constancia de la versión de las anotaciones
(por ejemplo `org.Hs.eg.db`), que es de versión fija.
