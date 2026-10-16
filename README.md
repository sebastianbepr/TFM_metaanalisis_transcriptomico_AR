# Meta-análisis transcriptómico del eje androgénico en tumores urológicos

Repositorio de scripts del Trabajo Fin de Máster (Máster Universitario en Bioinformática,
Universidad Europea de Madrid).

**Autor:** Sebastián Andrés Beltrán Prieto
**Tema:** meta-análisis transcriptómico integrativo y deconvolución celular para localizar la
señal del receptor de andrógenos (AR) en adenocarcinoma de próstata y tumores germinales de
testículo, usando el carcinoma urotelial de vejiga como comparador no androgénico.

Los scripts están numerados por etapa y reproducen, en orden, toda la metodología: desde la
descarga de las cohortes públicas hasta las figuras finales del TFM.

---

## Estructura del repositorio

```
01_descarga/              Descarga de las cohortes (GEO, TCGA-TGCT + GTEx)
02_estandarizacion/       Estandarizacion de cada dataset al esquema comun (.rds)
03_dge/                   Expresion diferencial por dataset (DESeq2 / limma)
04_auditoria_remediacion/ Auditoria de los objetos y remediacion de matrices
05_metaanalisis/          Objetos DExMA, meta-analisis, Stouffer firmado y contraste de firmas
06_sensibilidad/          Analisis de sensibilidad y control de calidad (I2, leave-one-out)
07_firmas/                Tabla maestra y figuras de las firmas consenso
08_enriquecimiento/       ORA y GSEA (GO:BP, Reactome, KEGG) y sus figuras
09_deconvolucion/         Deconvolucion del microambiente (immunedeconv + CIBERSORTx)
10_ar_estroma/            Asociacion de AR con el compartimento de celulas de Leydig
```

Cada carpeta contiene scripts con un prefijo numérico que marca el orden de ejecución dentro
de la etapa (`01_...`, `02_...`).

---

## Cómo ejecutar

1. **Datos.** Los datos crudos **no** se incluyen en el repositorio (son públicos y pesados).
   Se descargan con los scripts de `01_descarga/`. El repositorio contiene solo el código.

2. **Ruta de trabajo.** Cada script empieza con una única línea editable:

   ```r
   ruta <- "C:/Users/sebas/OneDrive/Escritorio/BIOINFORMATICA/UE/TFM/DATOS/Data_TFM"
   setwd(ruta)
   ```

   Basta con cambiar esa ruta por la carpeta de trabajo de cada equipo. Los scripts de
   descarga de TCGA/GTEx usan una carpeta propia (`C:/Users/sebas/TCGA_GTEX`), también
   editable en su cabecera.

3. **Orden.** Se ejecutan las carpetas en orden (01 → 10) y, dentro de cada una, por su prefijo
   numérico. La salida de una etapa es la entrada de la siguiente (cada script indica qué lee y
   qué escribe en su cabecera, en las líneas `Entra:` / `Sale:`).

4. **CIBERSORT.** El paso de CIBERSORT (etapa 09) se ejecuta en línea sobre la plataforma
   CIBERSORTx a partir de las matrices `*_mixture.txt` que genera `09_deconvolucion/01_...R`.
   Los resultados descargados de CIBERSORTx se colocan en `Deconvolucion/CIBERSORT_Output/`
   antes de correr `09_deconvolucion/02_...R`.

---

## Flujo de datos (resumen)

```
descarga  ->  estandarizacion (*_estandarizado.rds)
          ->  DGE por dataset (*_dge_*.rds)
          ->  auditoria + remediacion (*_remediado.rds)
          ->  objetos DExMA por grupo (objetos_ma_por_grupo.rds)
          ->  meta-analisis + Stouffer firmado + robustez
          ->  contraste entre grupos  ->  Firma 1 (287 genes) y Firma 2 (113 genes)
          ->  enriquecimiento (ORA / GSEA)
          ->  deconvolucion del microambiente
          ->  asociacion AR - celulas de Leydig
```

La cohorte integrada son 11 conjuntos de datos (966 muestras: 533 tumorales, 433 controles) de
próstata, testículo y vejiga.

---

## Requisitos

Ver `requisitos_R.md` para el listado de paquetes de R (CRAN y Bioconductor) por etapa.

## Declaración de uso de IA

Ver `DECLARACION_USO_IA.md`.
