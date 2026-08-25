# ==========================================================================
# Script: 05_extraccion_multisensor_demo.R
# Descripción: Extracción LIVIANA de NDVI para MODIS y Sentinel-2, sobre la
# misma área que Landsat. Es una DEMO para mostrar que el patrón de
# extracción se puede extender a otros sensores -- NO incluye limpieza de
# nubes rigurosa, gap-filling, ni fenología. No usar para análisis, solo
# para graficar junto a la serie Landsat ya validada.
# ==========================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, terra, sf, rgee, readr)

# Reutiliza el mismo truco que 01_extraccion_landsat.R: reduceRegions()
# llamado como método de instancia (no ee$Image$reduceRegions(image=...)),
# que es lo que evita el bug de compatibilidad con earthengine-api actual.
extraer_puntos_ee <- function(stack_ee, grilla_sf, scale) {
  grilla_ee <- rgee::sf_as_ee(grilla_sf)
  fc <- stack_ee$reduceRegions(collection = grilla_ee, reducer = ee$Reducer$mean(), scale = scale)
  info <- fc$getInfo()
  propiedades_a_df <- function(props) {
    props[sapply(props, is.null)] <- NA
    as.data.frame(props, stringsAsFactors = FALSE)
  }
  dplyr::bind_rows(lapply(info$features, function(f) propiedades_a_df(f$properties)))
}

generar_grilla <- function(area, grid_scale) {
  grilla <- sf::st_make_grid(area, cellsize = grid_scale, what = "centers")
  grilla <- sf::st_sf(geometry = grilla)
  grilla <- sf::st_filter(grilla, area)
  grilla$id_row <- seq_len(nrow(grilla))
  coords <- sf::st_coordinates(grilla)
  grilla$x <- coords[, "X"]
  grilla$y <- coords[, "Y"]
  grilla
}

# 1. MODIS (MOD13Q1) ----------------------------------------------------------
#' Extrae NDVI MODIS (MOD13Q1, ya compuesto 16 días y filtrado de nubes) --
#' extracción liviana de demo, sin limpieza adicional
extraer_modis_demo <- function(shp_path, fecha_inicio, fecha_fin, nombre_area,
                                gee_user, output_dir = "02_output", grid_scale = 250) {
  cat("Inicializando GEE (MODIS demo)...\n")
  rgee::ee_Initialize(user = gee_user, drive = FALSE)  # no usamos Drive en ningún paso

  area <- sf::st_read(shp_path, quiet = TRUE)
  area <- sf::st_sf(geometry = sf::st_geometry(area))  # solo geometria, ver nota en 01_extraccion_landsat.R
  region <- rgee::sf_as_ee(area)$geometry()

  coleccion <- ee$ImageCollection("MODIS/061/MOD13Q1")$
    filterDate(fecha_inicio, fecha_fin)$
    filterBounds(region)$
    select("NDVI")$
    map(function(img) {
      img$multiply(0.0001)$
        rename(ee$Date(img$get("system:time_start"))$format("YYYYMMdd"))$
        copyProperties(img, list("system:index", "system:time_start"))
    })

  grilla <- generar_grilla(area, grid_scale)
  stack <- coleccion$toBands()$clip(region)

  cat("Extrayendo (", nrow(grilla), "puntos, resolución", grid_scale, "m)...\n")
  valores <- extraer_puntos_ee(stack, grilla, scale = grid_scale)
  resultado <- valores %>% dplyr::relocate(id_row, x, y)

  carpeta <- file.path(output_dir, paste0(nombre_area, "_", format(Sys.Date(), "%Y%m%d")))
  dir.create(carpeta, recursive = TRUE, showWarnings = FALSE)
  ruta <- file.path(carpeta, paste0(nombre_area, "_modis_ndvi_raw.csv"))
  readr::write_csv(resultado, ruta)
  cat("✔️ MODIS demo guardado en:", ruta, "\n")
  list(csv_path = ruta, n_puntos = nrow(grilla))
}

# 2. Sentinel-2 (COPERNICUS/S2_SR_HARMONIZED) ----------------------------------
#' Extrae NDVI Sentinel-2 -- extracción liviana de demo, enmascarado de nubes
#' simple vía la banda de probabilidad de nube (umbral fijo), sin el
#' refinamiento que tiene el enmascarado QA_PIXEL de Landsat
extraer_sentinel_demo <- function(shp_path, fecha_inicio, fecha_fin, nombre_area,
                                   gee_user, output_dir = "02_output", grid_scale = 100,
                                   max_puntos = 300) {
  cat("Inicializando GEE (Sentinel-2 demo)...\n")
  rgee::ee_Initialize(user = gee_user, drive = FALSE)  # no usamos Drive en ningún paso

  area <- sf::st_read(shp_path, quiet = TRUE)
  area <- sf::st_sf(geometry = sf::st_geometry(area))  # solo geometria, ver nota en 01_extraccion_landsat.R
  region <- rgee::sf_as_ee(area)$geometry()

  enmascarar_nubes_s2 <- function(img) {
    qa <- img$select("QA60")
    nube_bit <- bitwShiftL(1L, 10)
    cirrus_bit <- bitwShiftL(1L, 11)
    mascara <- qa$bitwiseAnd(nube_bit)$eq(0)$And(qa$bitwiseAnd(cirrus_bit)$eq(0))
    img$updateMask(mascara)
  }

  coleccion <- ee$ImageCollection("COPERNICUS/S2_SR_HARMONIZED")$
    filterDate(fecha_inicio, fecha_fin)$
    filterBounds(region)$
    filter(ee$Filter$lt("CLOUDY_PIXEL_PERCENTAGE", 40))$
    map(function(img) {
      img_masc <- enmascarar_nubes_s2(img)
      ndvi <- img_masc$normalizedDifference(c("B8", "B4"))$
        rename(ee$Date(img$get("system:time_start"))$format("YYYYMMdd"))
      ndvi$copyProperties(img, list("system:index", "system:time_start"))
    })

  grilla <- generar_grilla(area, grid_scale)
  # Sentinel-2 a 10-20m sobre áreas grandes puede pesar mucho -- si la grilla
  # sale muy grande, se achica a una muestra aleatoria para que la demo sea
  # liviana (esto es a propósito, no es la extracción completa)
  if (nrow(grilla) > max_puntos) {
    set.seed(123)
    grilla <- grilla[sample(nrow(grilla), max_puntos), ]
    cat("Grilla reducida a una muestra de", max_puntos, "puntos (demo liviana)\n")
  }

  stack <- coleccion$toBands()$clip(region)

  cat("Extrayendo (", nrow(grilla), "puntos, resolución", grid_scale, "m)...\n")
  valores <- extraer_puntos_ee(stack, grilla, scale = grid_scale)
  resultado <- valores %>% dplyr::relocate(id_row, x, y)

  carpeta <- file.path(output_dir, paste0(nombre_area, "_", format(Sys.Date(), "%Y%m%d")))
  dir.create(carpeta, recursive = TRUE, showWarnings = FALSE)
  ruta <- file.path(carpeta, paste0(nombre_area, "_sentinel_ndvi_raw.csv"))
  readr::write_csv(resultado, ruta)
  cat("✔️ Sentinel-2 demo guardado en:", ruta, "\n")
  list(csv_path = ruta, n_puntos = nrow(grilla))
}

# 3. Gráfico comparativo de 3 sensores -----------------------------------------
#' Grafica NDVI promedio (crudo, sin limpiar) de Landsat, MODIS y Sentinel-2
#' juntos -- para mostrar potencial multisensor, no para validar
graficar_multisensor_demo <- function(landsat_raw_csv, modis_csv, sentinel_csv, titulo = "Comparación multisensor (demo)") {
  leer_serie <- function(ruta, sensor) {
    d <- readr::read_csv(ruta, show_col_types = FALSE)
    columnas_fecha <- setdiff(names(d), c("id_row", "x", "y"))
    d %>%
      tidyr::pivot_longer(all_of(columnas_fecha), names_to = "columna", values_to = "ndvi") %>%
      dplyr::mutate(fecha = as.Date(stringr::str_extract(columna, "\\d{8}$"), format = "%Y%m%d")) %>%
      dplyr::filter(!is.na(ndvi), !is.na(fecha)) %>%
      dplyr::group_by(fecha) %>%
      dplyr::summarise(ndvi_medio = mean(ndvi), .groups = "drop") %>%
      dplyr::mutate(sensor = sensor)
  }

  serie <- dplyr::bind_rows(
    leer_serie(landsat_raw_csv, "Landsat 8"),
    leer_serie(modis_csv, "MODIS"),
    leer_serie(sentinel_csv, "Sentinel-2")
  )

  ggplot2::ggplot(serie, ggplot2::aes(x = fecha, y = ndvi_medio, color = sensor)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1) +
    ggplot2::labs(title = titulo,
                  subtitle = "Demo -- series crudas sin limpiar, para mostrar potencial multisensor",
                  x = NULL, y = "NDVI promedio", color = "Sensor") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::scale_color_manual(values = c("Landsat 8" = "#1a7a3c", "MODIS" = "#b2182b", "Sentinel-2" = "#2166ac"))
}
