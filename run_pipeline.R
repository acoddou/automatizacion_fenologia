# ==========================================================================
# Script: run_pipeline.R
# Descripción: Orquestador único del pipeline. Corre de punta a punta:
#   1. Extracción de la serie de tiempo NDVI Landsat (Earth Engine)
#   2. Regularización temporal y gap-filling
#   3. Extracción de métricas fenológicas (Savitzky-Golay + Zhang, phenofit)
#   4. Mapas y gráficos (NDVI, métricas fenológicas -- ver 06_visualizaciones.R)
#   5. (opcional) Muestreo estratificado + clasificación de cobertura
#      -- solo si se pasan conaf_path y raster_clasificacion_path
#
# Uso:
#   source("run_pipeline.R")
#   resultado <- correr_pipeline(
#     shp_path     = "00_input/limite_apoquindo.shp",
#     fecha_inicio = "2023-01-01",
#     fecha_fin    = "2024-12-31",
#     nombre_area  = "apoquindo",
#     gee_user     = "agustincoddoudiaz@gmail.com",
#     conaf_path   = "00_input/conaf_apoquindo.shp"  # opcional, agrega desglose por cobertura
#   )
# ==========================================================================

source("01_pipeline/01_extraccion_landsat.R")
source("01_pipeline/02_limpieza_series.R")
source("01_pipeline/03_fenologia_phenofit.R")
source("01_pipeline/04_muestreo_espacial.R")
source("01_pipeline/06_visualizaciones.R")

#' Corre el pipeline completo de automatización de fenología, de punta a punta
#'
#' @param shp_path Ruta al shapefile del área de interés
#' @param fecha_inicio,fecha_fin Rango de fechas a analizar ("YYYY-MM-DD")
#' @param nombre_area Nombre identificador del área (usado en carpetas y archivos de salida)
#' @param gee_user Usuario de Google Earth Engine
#' @param output_dir Carpeta de salida (default "02_output")
#' @param grid_scale Espaciamiento (m) de la grilla de puntos de muestreo (default 150)
#' @param periods_per_year Periodos por año para regularizar la serie (default 24)
#' @param sos_eos_trs Umbral de amplitud estacional para SOS/EOS (default 0.2)
#' @param conaf_path Ruta a la capa de cobertura (CONAF u otra) -- opcional, si se
#'   entrega junto con raster_clasificacion_path se agrega el Paso 4 (clasificación)
#' @param raster_clasificacion_path Ruta al raster RGB+NDVI para clasificación -- opcional
#' @param campo_cobertura Campo de la capa de cobertura a usar (default "USO")
correr_pipeline <- function(shp_path,
                             fecha_inicio,
                             fecha_fin,
                             nombre_area,
                             gee_user,
                             output_dir = "02_output",
                             grid_scale = 150,
                             periods_per_year = 24,
                             sos_eos_trs = 0.2,
                             conaf_path = NULL,
                             raster_clasificacion_path = NULL,
                             campo_cobertura = "USO",
                             fecha_evento = NULL,
                             evento_label = "Evento") {

  if (!file.exists(shp_path)) {
    stop(sprintf(
      "No encontre el shapefile del area en:\n  %s\nRevisa que el archivo (.shp, .shx, .dbf, .prj) este guardado ahi con ese nombre exacto.",
      shp_path
    ))
  }
  if (!is.null(conaf_path) && !file.exists(conaf_path)) {
    stop(sprintf(
      "No encontre la capa de cobertura (conaf_path) en:\n  %s\nSi no tenes esa capa para esta area, deja conaf_path = NULL.",
      conaf_path
    ))
  }

  incluye_clasificacion <- !is.null(conaf_path) && !is.null(raster_clasificacion_path)
  n_pasos <- if (incluye_clasificacion) 5 else 4
  t0 <- Sys.time()

  cat(strrep("=", 72), "\n")
  cat("  PIPELINE AUTOMATIZACION_FENOLOGIA\n")
  cat("  Área:", nombre_area, " | Período:", fecha_inicio, "a", fecha_fin, "\n")
  cat(strrep("=", 72), "\n\n")

  # -- Paso 1: extracción -----------------------------------------------------
  cat(sprintf(">> [Paso 1/%d] Extrayendo serie de tiempo NDVI desde Google Earth Engine (Landsat)...\n", n_pasos))
  res_extraccion <- extraer_landsat(
    shp_path = shp_path, fecha_inicio = fecha_inicio, fecha_fin = fecha_fin,
    nombre_area = nombre_area, gee_user = gee_user, output_dir = output_dir,
    grid_scale = grid_scale
  )
  cat(sprintf("   listo (%s puntos) -- %.1f min transcurridos\n\n",
              res_extraccion$n_puntos, as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  # -- Paso 2: limpieza ---------------------------------------------------------
  cat(sprintf(">> [Paso 2/%d] Regularizando la serie temporal y rellenando huecos (nubes)...\n", n_pasos))
  res_limpieza <- limpiar_series(
    input_csv = res_extraccion$csv_path, nombre_area = nombre_area,
    periods_per_year = periods_per_year
  )
  cat(sprintf("   listo -- %.1f min transcurridos\n\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  # -- Paso 3: fenología ----------------------------------------------------------
  cat(sprintf(">> [Paso 3/%d] Ajustando curvas (Savitzky-Golay + Zhang) y extrayendo métricas fenológicas...\n", n_pasos))
  cat("   (este paso es el más lento -- ajusta una curva por punto)\n")
  res_fenologia <- extraer_fenologia(
    input_csv = res_limpieza$csv_path, nombre_area = nombre_area,
    periods_per_year = periods_per_year, sos_eos_trs = sos_eos_trs
  )
  cat(sprintf("   listo (%s puntos) -- %.1f min transcurridos\n\n",
              res_fenologia$n_puntos, as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  resultado <- list(extraccion = res_extraccion, limpieza = res_limpieza, fenologia = res_fenologia)

  # -- Paso 4: mapas y gráficos --------------------------------------------------
  cat(sprintf(">> [Paso 4/%d] Generando mapas y gráficos...\n", n_pasos))
  archivos_viz <- generar_visualizaciones(
    resultado, nombre_area = nombre_area, conaf_path = conaf_path, campo_cobertura = campo_cobertura,
    fecha_fin = fecha_fin, fecha_evento = fecha_evento, evento_label = evento_label
  )
  cat(sprintf("   listo -- %.1f min transcurridos\n\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  resultado$visualizaciones <- archivos_viz

  # -- Paso 5 (opcional): muestreo + clasificación ---------------------------
  if (incluye_clasificacion) {
    cat(sprintf(">> [Paso 5/%d] Muestreo estratificado y clasificación de cobertura (RandomForest)...\n", n_pasos))
    res_clasif <- muestrear_y_clasificar(
      raster_path = raster_clasificacion_path, conaf_path = conaf_path,
      nombre_area = nombre_area, output_dir = output_dir, campo_cobertura = campo_cobertura
    )
    cat(sprintf("   listo -- %.1f min transcurridos\n\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    resultado$clasificacion <- res_clasif
  }

  t_total <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(strrep("=", 72), "\n")
  cat(sprintf("  PIPELINE COMPLETO en %.1f minutos.\n", t_total))
  cat("  Resultados en:", file.path(output_dir, paste0(nombre_area, "_", format(Sys.Date(), "%Y%m%d"))), "\n")
  cat(strrep("=", 72), "\n")

  invisible(resultado)
}
