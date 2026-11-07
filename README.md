# Análisis de expresión diferencial con DESeq2 (TCGA-BRCA)

Tarea 02 – Bioinformática y Análisis Genómico. Alejandro Román Ruiz.

Contraste: tejido tumoral de mama (`Primary_Solid_Tumor`) frente a tejido normal adyacente (`Solid_Tissue_Normal`), ajustando por placa de secuenciación.
Diseño del modelo: `~ site + plate`.

## Ejecución

Desde la raíz del proyecto (las rutas del script son relativas):

```bash
conda env create -p ./.venv -f environment.yml
conda activate ./.venv
Rscript src/01_run_deseq2.R
```

## Estructura

- `data/raw/counts.csv.gz`: conteos crudos, genes en filas y muestras en columnas.
- `data/raw/metadata.csv.gz`: una fila por muestra, columnas `site` y `plate`.
- `src/01_run_deseq2.R`: script del análisis (8 pasos).
- `results/`: tablas (`.tsv`) y figuras (`.png`) generadas.
- `logs/`: `sessionInfo.txt` (versiones de R y paquetes) y `run_summary.txt`.

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

### Interpretación biológica

Los genes sobreexpresados en tumor (2115) muestran sobre todo una **firma proliferativa**:
GO BP prioriza la segregación cromosómica y el ciclo celular (*chromosome segregation*,
padj 3.8e-17; *cell cycle process*) y Reactome lo confirma con *Cell Cycle*,
*Cell Cycle, Mitotic* y *Cell Cycle Checkpoints* (padj 9.9e-22). A nivel de gen lo apoyan
PKMYT1, NUF2 y el marcador de proliferación MKI67 (log2FC 3.39, padj 1.4e-23).
Un segundo componente up es la **remodelación de la matriz extracelular**: entre los genes
más significativos están MMP1, MMP11, COL10A1 y COL11A1, compatibles con interacción
tumor-estroma y potencial invasivo. Además aparece una **firma inmune** (*immune response*,
*adaptive immune response*, *lymphocyte activation*), compatible con infiltración linfocitaria.

Los genes con menor expresión en tumor (1939) se asocian a **metabolismo lipídico y
adipogénesis** (*lipid metabolic process*, *fat cell differentiation*; Reactome *Adipogenesis*,
*Metabolism of lipids*), a **desarrollo vascular** (*circulatory system development*,
*angiogenesis*) y a **contracción muscular**, con genes como ATP1A2, GPAM, RDH5 y MME.

En conjunto, el tumor combina proliferación activa, remodelación estromal y un componente
inmune, con una pérdida relativa de los programas metabólicos y de diferenciación del
tejido mamario normal.

**Limitación:** es RNA-seq *bulk*, así que cada muestra es una mezcla de tipos celulares.
El tejido mamario normal es rico en adipocitos, por lo que parte de la firma down (y de la
inmune up) puede reflejar diferencias de composición celular entre muestras más que cambios
de regulación dentro de las mismas células.
