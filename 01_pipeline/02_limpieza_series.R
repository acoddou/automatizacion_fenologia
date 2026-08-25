# ==========================================================================
# Script: 02_limpieza_series.R
# Descripción: Regulariza la serie de tiempo NDVI cruda (una columna por
# fecha de paso de Landsat, espaciado irregular) a una grilla temporal fija
# de periodos por año, y rellena los huecos que dejó el enmascarado de nubes
# de 01_extraccion_landsat.R (QA_PIXEL) mediante interpolación temporal.
# Adaptado de 04_limpieza_landsat.R del repo residencia_aguas_de_ramon,
# generalizado para leer el CSV de cualquier área/rango de fechas en vez de
# un raster con años hardcodeados.
# ==========================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse,
  lubridate,
  zoo,
  readr
)

# 2. Función principal -------------------------------------------------------
#' Regulariza y rellena huecos en la serie NDVI cruda de Landsat
#'
#' @param input_csv Ruta al CSV crudo generado por extraer_landsat()
#'   (columnas: id_row, x, y, y una columna NDVI por fecha)
#' @param nombre_area Nombre identificador del área (usado en el nombre del archivo de salida)
#' @param output_dir Carpeta donde se guardará el CSV limpio (default: misma carpeta que input_csv)
#' @param periods_per_year Cantidad de periodos en que se divide cada año para regularizar
#'   la serie (default 24, ~15 días por periodo, igual que en la residencia)
limpiar_series <- function(input_csv,
                            nombre_area,
                            output_dir = dirname(input_csv),
                            periods_per_year = 24) {

  # 2.1 Cargar el CSV crudo ---------------------------------------------------
  if (!file.exists(input_csv)) {
    stop("No se encontró el CSV de entrada: ", input_csv)
  }
  datos_crudos <- readr::read_csv(input_csv, show_col_types = FALSE)
  cat("CSV crudo cargado:", nrow(datos_crudos), "puntos,", ncol(datos_crudos) - 3, "fechas.\n")

  # 2.2 Extraer fecha, año y periodo de cada columna NDVI ---------------------
  # Cada columna viene nombrada tipo "LC08_<path_row>_<fecha>_<fecha>"; se
  # extrae el primer bloque de 8 dígitos rodeado de "_" como la fecha real.
  columnas_fecha <- setdiff(names(datos_crudos), c("id_row", "x", "y"))
  info_columnas <- tibble(
    columna  = columnas_fecha,
    date_str = str_extract(columna, "(?<=_)\\d{8}(?=_)")
  ) %>%
    mutate(
      date    = ymd(date_str),
      year    = year(date),
      doy     = yday(date),
      periodo = pmin(ceiling(doy / (365.25 / periods_per_year)), periods_per_year)
    )

  target_years <- sort(unique(info_columnas$year))
  cat("Regularizando a", periods_per_year, "periodos/año, años:",
      paste(range(target_years), collapse = "-"), "\n")

  # 2.3 Agregar por periodo e interpolar los huecos ---------------------------
  datos_largo <- datos_crudos %>%
    pivot_longer(cols = all_of(columnas_fecha), names_to = "columna", values_to = "ndvi_valor") %>%
    filter(!is.na(ndvi_valor)) %>%                    # descarta píxeles enmascarados por nubes
    inner_join(info_columnas, by = "columna")

  datos_regularizados <- datos_largo %>%
    group_by(id_row, x, y, year, periodo) %>%
    # si hay más de una pasada válida en el mismo periodo, se toma el valor
    # más alto: dentro de un periodo corto, más NDVI = menos contaminación
    # residual de nubes/sombra que el enmascarado QA_PIXEL no haya sacado
    summarise(ndvi_agg = max(ndvi_valor, na.rm = TRUE), .groups = "drop") %>%
    tidyr::complete(nesting(id_row, x, y), year = target_years, periodo = 1:periods_per_year) %>%
    arrange(id_row, year, periodo) %>%
    group_by(id_row) %>%
    mutate(ndvi_interpolado = zoo::na.approx(ndvi_agg, na.rm = FALSE, rule = 2)) %>%
    # puntos sin ningún dato válido en todo el rango: no hay de dónde interpolar
    mutate(ndvi_interpolado = ifelse(is.na(ndvi_interpolado), 0, ndvi_interpolado)) %>%
    ungroup()

  # 2.4 Formatear a ancho y exportar -------------------------------------------
  cat("Pivotando a formato ancho y exportando...\n")

  datos_ancho <- datos_regularizados %>%
    mutate(periodo_nombre = sprintf("NDVI_%d_P%02d", year, periodo)) %>%
    select(id_row, x, y, periodo_nombre, ndvi_interpolado) %>%
    pivot_wider(names_from = periodo_nombre, values_from = ndvi_interpolado)

  orden_columnas <- datos_regularizados %>%
    distinct(year, periodo) %>%
    arrange(year, periodo) %>%
    mutate(periodo_nombre = sprintf("NDVI_%d_P%02d", year, periodo)) %>%
    pull(periodo_nombre)

  datos_ancho <- datos_ancho %>% select(id_row, x, y, all_of(orden_columnas))

  ruta_csv <- file.path(output_dir, paste0(nombre_area, "_landsat_ndvi_limpio.csv"))
  readr::write_csv(datos_ancho, ruta_csv)

  cat("✔️ Limpieza completa. CSV guardado en:", ruta_csv, "\n")
  return(list(csv_path = ruta_csv, n_puntos = nrow(datos_ancho), n_periodos = length(orden_columnas)))
}
