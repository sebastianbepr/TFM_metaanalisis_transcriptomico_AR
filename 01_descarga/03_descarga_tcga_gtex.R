# ______________________________________________________________________________
# 01_03  Descarga de la cohorte TCGA-TGCT + GTEx (testiculo)
#
# Que hace: descarga de UCSC Xena la matriz de conteos Toil y el fenotipo, extrae
#           las muestras de testiculo (tumor TCGA-TGCT y normal GTEx), remapea los
#           IDs Ensembl a simbolo (probemap gencode v23) y guarda la matriz lista.
# Entra:    UCSC Xena (matriz Toil de conteos, phenotype, probemap gencode v23)
# Sale:     testiculo_TCGA_GTEx_conteos.rds, grupos_testiculo.csv (carpeta TCGA_GTEX)
# ______________________________________________________________________________


# 01. DESCARGA _________________________________________________________________

install.packages(c("UCSCXenaTools", "data.table", "R.utils"))

library(UCSCXenaTools)
library(data.table)

# Ruta de trabajo: unica linea que hay que editar para reproducir el pipeline.
ruta <- "C:/Users/sebas/TCGA_GTEX"
setwd(ruta)
dir_datos <- ruta   # el resto del script usa dir_datos sin cambios

if (!dir.exists(dir_datos)) dir.create(dir_datos)

# Paso 1. Descargar la metadata de fenotipo (archivo pequeno).
cat(">>> Descargando metadata de fenotipo (TCGA TARGET GTEx)...\n")
pheno_q <- XenaGenerate(subset = XenaHostNames == "toilHub") |>
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") |>
  XenaFilter(filterDatasets = "TcgaTargetGTEX_phenotype")
XenaQuery(pheno_q) |> XenaDownload(destdir = dir_datos)

# Localizamos el archivo de fenotipo descargado.
f_pheno <- list.files(dir_datos, pattern = "phenotype", full.names = TRUE)[1]
pheno <- fread(f_pheno, sep = "\t", header = TRUE)
cat(sprintf("Metadata cargada: %d muestras, columnas: %s\n",
            nrow(pheno), paste(names(pheno), collapse=", ")))

# Paso 2. Identificar las muestras de testiculo.
# La metadata trae columnas tipo sample, _primary_site, _sample_type, _study
# (TCGA/GTEX/TARGET), detailed_category, primary disease or tissue... Los nombres
# exactos pueden variar, asi que los detectamos de forma flexible.
nm <- names(pheno)
col_sample <- nm[grepl("^sample$|sampleID|^id$", nm, ignore.case = TRUE)][1]
col_site   <- nm[grepl("primary.?site|tissue", nm, ignore.case = TRUE)][1]
col_study  <- nm[grepl("study", nm, ignore.case = TRUE)][1]
col_type   <- nm[grepl("sample.?type", nm, ignore.case = TRUE)][1]
col_cat    <- nm[grepl("category|disease|detailed", nm, ignore.case = TRUE)][1]

cat(sprintf("\nColumnas detectadas: sample=%s, site=%s, study=%s, type=%s, cat=%s\n",
            col_sample, col_site, col_study, col_type, col_cat))

# Mostramos los valores de sitio que contienen 'test' para verificar los nombres.
sitios_test <- unique(pheno[[col_site]][grepl("test", pheno[[col_site]], ignore.case = TRUE)])
cat("Valores de sitio relacionados con testículo:", paste(sitios_test, collapse=" | "), "\n")

# Tumores TCGA-TGCT: study = TCGA y sitio = testis.
es_tcga <- grepl("TCGA", pheno[[col_study]], ignore.case = TRUE)
es_testis_site <- grepl("test", pheno[[col_site]], ignore.case = TRUE)
tumores <- pheno[[col_sample]][es_tcga & es_testis_site]

# Normales GTEx Testis.
es_gtex <- grepl("GTEX", pheno[[col_study]], ignore.case = TRUE)
normales <- pheno[[col_sample]][es_gtex & es_testis_site]

cat(sprintf("\n>>> Muestras identificadas:\n"))
cat(sprintf("    Tumores  (TCGA-TGCT): %d\n", length(tumores)))
cat(sprintf("    Normales (GTEx Testis): %d\n", length(normales)))

if (length(tumores) == 0 || length(normales) == 0)
  stop("No se identificaron muestras. Revisa los nombres de columnas/valores arriba.")

# Tabla de grupos para el DGE.
grupos <- data.frame(
  sample = c(tumores, normales),
  grupo  = c(rep("Tumor", length(tumores)), rep("Normal", length(normales))),
  origen = c(rep("TCGA", length(tumores)),  rep("GTEx", length(normales))),
  stringsAsFactors = FALSE
)
fwrite(grupos, file.path(dir_datos, "grupos_testiculo.csv"))
cat(sprintf("    Tabla de grupos guardada (%d muestras)\n", nrow(grupos)))

# Paso 3. Descargar la matriz de conteos (grande) y leer solo testiculo.
cat("\n>>> Descargando matriz de conteos (~2.5 GB, una sola vez)...\n")
counts_q <- XenaGenerate(subset = XenaHostNames == "toilHub") |>
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") |>
  XenaFilter(filterDatasets = "TcgaTargetGtex_gene_expected_count")
XenaQuery(counts_q) |> XenaDownload(destdir = dir_datos)

f_counts <- list.files(dir_datos, pattern = "expected_count", full.names = TRUE)[1]

# Lectura selectiva: solo la 1a columna (genes) + las columnas de testiculo.
# Primero leemos solo la cabecera para saber que columnas existen.
cab <- fread(f_counts, nrows = 0)
todas_cols <- names(cab)
col_genes <- todas_cols[1]   # primera col = identificador de gen
cols_testiculo <- intersect(grupos$sample, todas_cols)
cat(sprintf("De %d muestras de testículo, %d están en la matriz de conteos.\n",
            nrow(grupos), length(cols_testiculo)))

# Leemos solo esas columnas para ahorrar memoria.
mat <- fread(f_counts, select = c(col_genes, cols_testiculo))
cat(sprintf("Matriz de testículo leída: %d genes × %d muestras\n",
            nrow(mat), ncol(mat)-1))

# Paso 4. Convertir IDs Ensembl a simbolo de gen (probemap Gencode v23).
cat("\n>>> Descargando mapa de genes (Ensembl -> símbolo)...\n")
probemap_url <- "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap/gencode.v23.annotation.gene.probemap"
f_map <- file.path(dir_datos, "gencode.v23.gene.probemap")
if (!file.exists(f_map)) download.file(probemap_url, f_map, mode = "wb")
probemap <- fread(f_map, header = TRUE)
# probemap: columnas tipo 'id' (Ensembl) y 'gene' (simbolo).
setnames(probemap, names(probemap)[1:2], c("id", "gene"))

# Unimos matriz y probemap.
setnames(mat, col_genes, "id")
mat <- merge(probemap[, .(id, gene)], mat, by = "id")
cat(sprintf("Mapeados a símbolo: %d genes\n", nrow(mat)))

# Paso 5. Guardar la matriz lista para el DGE.
# Estos conteos de Xena estan en log2(count+1). Para DESeq2 hay que revertir a
# conteos crudos; lo dejamos anotado para el script de DGE.
saveRDS(list(matriz = mat, grupos = grupos),
        file.path(dir_datos, "testiculo_TCGA_GTEx_conteos.rds"))
cat("\n>>> GUARDADO: testiculo_TCGA_GTEx_conteos.rds\n")
cat("    Contiene: matriz (genes×muestras, en log2(count+1)) + tabla de grupos.\n")
cat("    NOTA: para DESeq2, en el siguiente script se revierte log2 -> conteos crudos.\n")


# 02. VERIFICACION _____________________________________________________________

# Comprobaciones post-descarga de la cohorte TCGA+GTEx testiculo.

# Cargamos lo que guardo el bloque de descarga.
obj <- readRDS(file.path(dir_datos, "testiculo_TCGA_GTEx_conteos.rds"))
mat    <- obj$matriz
grupos <- obj$grupos

cat("=====================================================\n")
cat("VERIFICACIÓN 1 — Dimensiones y estructura\n")
cat("=====================================================\n")
cat(sprintf("Matriz: %d genes × %d columnas\n", nrow(mat), ncol(mat)))
cat(sprintf("Tabla de grupos: %d muestras\n", nrow(grupos)))
cat("Primeras columnas de la matriz:", paste(head(names(mat),4), collapse=", "), "...\n")

cat("\n=====================================================\n")
cat("VERIFICACIÓN 2 — Balance de grupos (tumor vs normal)\n")
cat("=====================================================\n")
print(table(grupo = grupos$grupo, origen = grupos$origen))
cat(sprintf("\nTumores (TCGA-TGCT): %d  | esperado ~150-156\n", sum(grupos$grupo=="Tumor")))
cat(sprintf("Normales (GTEx Testis): %d | esperado ~165\n", sum(grupos$grupo=="Normal")))

cat("\n=====================================================\n")
cat("VERIFICACIÓN 3 — Concordancia matriz <-> grupos\n")
cat("=====================================================\n")
muestras_mat <- setdiff(names(mat), c("id","gene"))
en_ambos <- intersect(muestras_mat, grupos$sample)
solo_mat <- setdiff(muestras_mat, grupos$sample)
solo_grp <- setdiff(grupos$sample, muestras_mat)
cat(sprintf("Muestras en matriz: %d\n", length(muestras_mat)))
cat(sprintf("Coinciden matriz<->grupos: %d\n", length(en_ambos)))
cat(sprintf("En grupos pero NO en matriz: %d %s\n", length(solo_grp),
            if(length(solo_grp)>0) paste0("(", paste(head(solo_grp,3),collapse=","), "...)") else ""))
if (length(en_ambos)==length(muestras_mat)) cat("OK: todas las columnas de la matriz tienen grupo asignado.\n")

cat("\n=====================================================\n")
cat("VERIFICACIÓN 4 — Rango de valores (¿están en log2?)\n")
cat("=====================================================\n")
# Tomamos una submuestra numerica.
num_cols <- muestras_mat[1:min(5,length(muestras_mat))]
vals <- as.matrix(mat[, ..num_cols])
cat(sprintf("Rango de valores: min=%.2f, max=%.2f, mediana=%.2f\n",
            min(vals,na.rm=TRUE), max(vals,na.rm=TRUE), median(vals,na.rm=TRUE)))
cat("Interpretación: si el max está ~15-20, están en log2(count+1) (esperado en Xena).\n")
cat("                si el max es enorme (miles+), serían conteos crudos.\n")

cat("\n=====================================================\n")
cat("VERIFICACIÓN 5 — Genes clave presentes (control de sanidad)\n")
cat("=====================================================\n")
# Comprobamos AR y genes de testiculo como control biologico.
genes_check <- c("AR","FOXA1","KLK3","NANOG","POU5F1","KIT","DDX4","CD8A","FOXP3")
col_gene <- if("gene" %in% names(mat)) "gene" else "id"
presentes <- genes_check[genes_check %in% mat[[col_gene]]]
cat("Genes clave encontrados:", paste(presentes, collapse=", "), "\n")
ausentes <- setdiff(genes_check, presentes)
if(length(ausentes)>0) cat("No encontrados (revisar mapeo):", paste(ausentes,collapse=", "), "\n")

cat("\n=====================================================\n")
cat("VERIFICACIÓN 6 — ¿Hay NAs o duplicados de genes?\n")
cat("=====================================================\n")
cat(sprintf("Genes duplicados (símbolo): %d\n", sum(duplicated(mat[[col_gene]]))))
cat(sprintf("NAs en la matriz numérica: %d\n", sum(is.na(vals))))
cat("\n>>> Si todo se ve bien, listo para el QC de batch (PCA) y el DGE.\n")


# 03. CORRECCION DE REMAPEO ____________________________________________________

# La matriz quedaba con solo 9.646 genes y AR ausente porque el merge con el
# probemap fallaba en parte: los IDs Ensembl traen sufijo de version. Rehacemos
# el mapeo quitando ese sufijo en ambos lados (matriz y probemap) antes de unir.
# Parte del archivo de conteos ya descargado, no vuelve a descargar la matriz.

# Recuperamos la tabla de grupos ya generada.
grupos <- fread(file.path(dir_datos, "grupos_testiculo.csv"))

# Releemos el archivo de conteos crudo, solo las columnas de testiculo.
f_counts <- list.files(dir_datos, pattern = "expected_count", full.names = TRUE)[1]
cat("Archivo de conteos:", basename(f_counts), "\n")
cab <- fread(f_counts, nrows = 0)
col_id_raw <- names(cab)[1]
cols_test <- intersect(grupos$sample, names(cab))
mat <- fread(f_counts, select = c(col_id_raw, cols_test))
setnames(mat, col_id_raw, "id_versionado")
cat(sprintf("Matriz cruda releída: %d genes × %d muestras\n", nrow(mat), ncol(mat)-1))

# Clave: quitar el sufijo de version del Ensembl ID.
mat[, id := sub("\\..*$", "", id_versionado)]   # ENSG00000000938.12 -> ENSG00000000938
cat("Ejemplo de limpieza:", mat$id_versionado[1], "->", mat$id[1], "\n")

# Releemos el probemap y limpiamos tambien su version.
f_map <- file.path(dir_datos, "gencode.v23.gene.probemap")
if (!file.exists(f_map)) {
  url <- "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap/gencode.v23.annotation.gene.probemap"
  download.file(url, f_map, mode = "wb")
}
probemap <- fread(f_map, header = TRUE)
setnames(probemap, names(probemap)[1:2], c("id_pm", "gene"))
probemap[, id := sub("\\..*$", "", id_pm)]   # limpiar version tambien
cat(sprintf("Probemap: %d entradas\n", nrow(probemap)))

# Merge por ID sin version.
mat_mapeada <- merge(probemap[, .(id, gene)], mat[, !"id_versionado"], by = "id")
cat(sprintf("\n>>> Tras merge correcto: %d genes (antes eran 9646)\n", nrow(mat_mapeada)))

# Verificamos los genes clave.
genes_check <- c("AR","FOXA1","KLK3","NANOG","POU5F1","KIT","DDX4","CD8A","FOXP3","NR3C4")
presentes <- genes_check[genes_check %in% mat_mapeada$gene]
cat("Genes clave AHORA encontrados:", paste(presentes, collapse=", "), "\n")
faltan <- setdiff(genes_check, presentes)
if(length(faltan)>0) cat("Aún faltan:", paste(faltan,collapse=", "), "\n")

# Colapsamos genes duplicados (varios Ensembl -> mismo simbolo): nos quedamos con
# la fila de mayor expresion media, que es la mas informativa.
num_cols <- setdiff(names(mat_mapeada), c("id","gene"))
mat_mapeada[, media_expr := rowMeans(.SD), .SDcols = num_cols]
setorder(mat_mapeada, gene, -media_expr)
mat_colapsada <- mat_mapeada[!duplicated(gene)]
mat_colapsada[, media_expr := NULL]
cat(sprintf("\nTras colapsar duplicados: %d genes únicos (antes %d con duplicados)\n",
            nrow(mat_colapsada), nrow(mat_mapeada)))
cat("Genes duplicados restantes:", sum(duplicated(mat_colapsada$gene)), "\n")

# Guardamos la matriz corregida.
saveRDS(list(matriz = mat_colapsada, grupos = grupos),
        file.path(dir_datos, "testiculo_TCGA_GTEx_conteos_CORREGIDO.rds"))
cat("\n>>> GUARDADO: testiculo_TCGA_GTEx_conteos_CORREGIDO.rds\n")
cat(sprintf("    %d genes únicos × %d muestras, en log2(count+1)\n",
            nrow(mat_colapsada), length(num_cols)))
cat("    AR presente:", "AR" %in% mat_colapsada$gene, "\n")
