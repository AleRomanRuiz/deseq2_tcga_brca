# Análisis de expresión diferencial con DESeq2 (TCGA-BRCA)

- **Estudiante:** Alejandro Román Ruiz
- **Asignatura:** Bioinformática y Análisis Genómico, Tarea 02
- **Lenguaje:** R (DESeq2, clusterProfiler, ReactomePA)
- **Repositorio:** https://github.com/AleRomanRuiz/deseq2_tcga_brca

Contraste: tejido tumoral de mama (`Primary_Solid_Tumor`) frente a tejido normal
adyacente (`Solid_Tissue_Normal`), ajustando por placa de secuenciación.
Diseño del modelo: `~ site + plate`.

## Ejecución

Desde la raíz del proyecto (las rutas del script son relativas):

```bash
conda env create -p ./.venv -f environment.yml
conda activate ./.venv
Rscript src/01_run_deseq2.R
```

## Datos de entrada

Los datos **no están en GitHub** (excluidos en `.gitignore`); se incluyen en el ZIP de la entrega.
Si se clona el repositorio, descargarlos del aula virtual y colocarlos con estos nombres exactos:

```
data/raw/counts.csv.gz     # conteos crudos: genes en filas, muestras en columnas
data/raw/metadata.csv.gz   # una fila por muestra: columnas site y plate
```

## Resultados

Se generan al ejecutar el script (tampoco están en GitHub):

- `results/deseq2_all_results.tsv`: todos los genes, ordenados por padj.
- `results/deseq2_sig_padj_0.05.tsv`: genes con padj < 0.05.
- `results/deseq2_de_padj_0.05_lfc_1.tsv`: genes DE (padj < 0.05 y |log2FC| > 1).
- `results/top10_up_genes.tsv`, `results/top10_down_genes.tsv`, `results/plot_top10_genes.png`.
- `results/plot_volcano.png`, `results/plot_ma.png`, `results/dispersion_plot.png`.
- `results/normalized_counts.tsv`: conteos normalizados por size factors.
- `results/top_gene_context.txt`: gen destacado (MKI67).
- `results/go_bp_up.tsv`, `results/go_bp_down.tsv` y sus figuras `plot_go_bp_*.png`.
- `results/reactome_up.tsv`, `results/reactome_down.tsv` y sus figuras `plot_reactome_*.png`.
- `logs/run_summary.txt`: resumen de la ejecución.
- `logs/sessionInfo.txt`: versiones de R y paquetes.

La interpretación biológica está en **`informe.pdf`**.

## Parámetros

- Filtro de genes: ≥ 10 lecturas en ≥ 20 muestras (tamaño del grupo más pequeño).
- Significación: `padj < 0.05` (Benjamini-Hochberg).
- Umbral biológico: `|log2FoldChange| > 1`.

Definidos al principio de `src/01_run_deseq2.R`.

## Resultados principales

| Métrica | Valor |
|---|---|
| Muestras | 75 |
| Genes tras el filtro | 17606 |
| Genes significativos (padj < 0.05) | 10240 |
| Genes DE (padj < 0.05 y \|log2FC\| > 1) | 4054 |
| Up en tumor / Down en tumor | 2115 / 1939 |

## Decisiones e incidencias

- **Entorno no reubicable.** El entorno conda se creó en otra ruta y, al mover la carpeta
  del proyecto, `Rscript` dejó de funcionar (las rutas absolutas quedan grabadas en la
  instalación). Se recreó desde `environment.yml`, del que se eliminaron `name:` y `prefix:`
  (rutas absolutas de la máquina original).
- **`r-ggplot2` explícito** en `environment.yml` y `library(ggplot2)` en el script: el script
  usa ggplot2 directamente y antes solo estaba instalado como dependencia de enrichplot.
- **Correcciones sobre el script de clase:** errores de sintaxis y de nombres de variables,
  y nombres de los ficheros de salida ajustados a los del enunciado.
- **Bloques añadidos:** figura top 10 (barras ordenadas por |log2FC| en cada panel),
  conteos normalizados (`counts(dds, normalized = TRUE)`) y gen destacado (MKI67, marcador
  de proliferación citado en el enunciado).
- **Top 10:** son los 10 genes DE con menor padj en cada dirección, no los de mayor |log2FC|.
- **Universo del enriquecimiento:** los genes DE (criterio del script de clase), no todos los
  genes analizados.
- **Diferencias numéricas:** mismos genes, orden y log2FC que el enunciado; los padj difieren
  ligeramente, probablemente por la versión de DESeq2 (ver `logs/sessionInfo.txt`).
- **Verificación:** el análisis se ejecutó desde cero (entorno nuevo a partir de
  `environment.yml`) y reprodujo exactamente las cifras de `logs/run_summary.txt`.
