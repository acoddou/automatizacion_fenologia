# Automatización de Análisis Fenológico — Landsat

[![DOI](https://zenodo.org/badge/1333526786.svg)](https://doi.org/10.5281/zenodo.23067268)

Pipeline automatizado para caracterizar la respuesta fenológica de cualquier área de interés
(bosque, matorral, cuenca, etc.) a partir de un shapefile de entrada, usando series de tiempo
NDVI de Landsat y extracción de métricas fenológicas 100% scripteada en R (sin TIMESAT).

Basado en la metodología desarrollada en
[residencia_aguas_de_ramon](https://github.com/acoddou/residencia_aguas_de_ramon),
adaptada para correr de extremo a extremo sobre cualquier área, sin intervención manual.

## Por qué Landsat

Landsat se eligió como sensor único del pipeline por:
- **Registro histórico largo (1984–presente)**, a diferencia de Sentinel-2 (2015+ con cobertura
  completa recién desde 2019), lo que permite correr el pipeline sobre cualquier área y rango de
  fechas sin quedar limitado por disponibilidad de datos.
- **Alta correlación con Sentinel-2 en métricas de productividad acumulada** (LIN R=0.95, AMP R=0.92),
  validado en el análisis de la cuenca Aguas de Ramón.
- **Resolución intermedia (30m)** que da suficiente detalle espacial sin la carga computacional de Sentinel.

### Advertencia de uso — métricas fenológicas de fecha (SOS/EOS/PET/LOS)

La correlación entre sensores para fechas específicas de eventos fenológicos es baja o negativa
(ej. SOS entre MODIS y Landsat: R=-0.69). El pipeline calcula estas métricas igual, pero el output
las marca explícitamente como **de menor confiabilidad** frente a LIN/AMP, que sí son robustas
entre sensores.

## Estructura

```
automatizacion_fenologia/
│
├── 00_input/
│   └── area.shp                      # Único input manual: el polígono del área a analizar
│
├── 01_pipeline/
│   ├── 00_inicializar_gee.R           # Autenticación rgee/Earth Engine (correr 1 vez por sesión de R)
│   ├── 01_extraccion_landsat.R       # rgee: NDVI Landsat recortado al shp, con enmascarado de nubes
│   │                                  # (QA_PIXEL) aplicado en Earth Engine antes de extraer
│   ├── 02_limpieza_series.R          # Regularización temporal (binning por periodo) y gap-filling
│   │                                  # sobre los huecos que dejó el enmascarado de nubes
│   ├── 03_fenologia_phenofit.R       # Suavizado + extracción SOS/EOS/PET/LOS/LIN/AMP (reemplaza TIMESAT)
│   ├── 04_muestreo_espacial.R        # Muestreo estratificado sobre el área (genérico, no fijo a Aguas de Ramón)
│   └── 05_outputs.R                  # Tablas, rasters y plots automáticos
│
├── 02_output/
│   └── <nombre_area>_<fecha_corrida>/
│       ├── metricas.csv
│       ├── rasters/
│       └── plots/
│
├── R/
│   └── funciones_auxiliares.R        # Funciones compartidas entre scripts
│
├── run_pipeline.R                    # Punto de entrada único
└── README.md
```

## Uso

```r
source("run_pipeline.R")

correr_analisis(
  shp_path    = "00_input/area.shp",
  fecha_inicio = "2019-01-01",
  fecha_fin    = "2024-12-31",
  nombre_area  = "mi_area"
)
```

## Requisitos

- R >= 4.2
- `rgee` (requiere cuenta de Google Earth Engine autenticada)
- `phenofit` (extracción de métricas fenológicas, reemplaza TIMESAT)
- `sf`, `terra`, `dplyr`, `ggplot2`

## Estado

🚧 En desarrollo — Paso actual: script de extracción Landsat (`01_extraccion_landsat.R`)
