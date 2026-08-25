# ==========================================================================
# Script: 04_muestreo_espacial.R
# Descripción: Cruza el raster RGB+NDVI con la capa de cobertura de CONAF,
# hace un muestreo aleatorio estratificado (N píxeles por categoría) y
# entrena un RandomForest para clasificar cobertura a partir de las bandas
# espectrales. Generaliza el patrón de 17_muestreo_modis.R de la residencia
# (umbral de dominancia + muestreo estratificado de 60 px/categoría), usando
# la capa de CONAF en vez de una clasificación manual de Landsat.
# ==========================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse,
  terra,
  sf,
  randomForest
)

# 2. Función principal -------------------------------------------------------
#' Muestrea estratificado por cobertura CONAF y entrena un RandomForest
#'
#' @param raster_path Ruta al raster RGB+NDVI (generado por extraer_raster_clasificacion())
#' @param conaf_path Ruta al shapefile de CONAF ya recortado al área de interés
#' @param nombre_area Nombre identificador del área
#' @param output_dir Carpeta de salida
#' @param campo_cobertura Campo de CONAF a usar como categoría. Default "USO"
#'   (3 clases: Bosques / Praderas y Matorrales / Áreas Urbanas e Industriales,
#'   el mismo nivel de detalle que el esquema de la residencia). "SUBUSO" (5
#'   clases) o "USO_TIERRA" (7 clases) dan más detalle si hace falta.
#' @param n_por_categoria Píxeles a muestrear por categoría (default 60, igual
#'   que 17_muestreo_modis.R de la residencia)
muestrear_y_clasificar <- function(raster_path,
                                    conaf_path,
                                    nombre_area,
                                    output_dir = "02_output",
                                    campo_cobertura = "USO",
                                    n_por_categoria = 60,
                                    semilla = 123) {

  set.seed(semilla)
  carpeta_area <- file.path(output_dir, paste0(nombre_area, "_", format(Sys.Date(), "%Y%m%d")))
  dir.create(carpeta_area, recursive = TRUE, showWarnings = FALSE)

  # 2.1 Cargar raster espectral y capa de cobertura ----------------------------
  if (!file.exists(raster_path)) stop("No se encontró el raster en: ", raster_path)
  if (!file.exists(conaf_path)) stop("No se encontró la capa de cobertura en: ", conaf_path)
  raster_espectral <- terra::rast(raster_path)
  # Brillo (promedio R+G+B): las áreas urbanas/industriales tienen NDVI
  # similar al de la vegetación (píxeles mixtos con techos/pavimento y
  # jardines a 30m), pero son notoriamente más brillantes en las 3 bandas --
  # se agrega como predictor extra para ayudar al RF a separarlas.
  raster_espectral$BRILLO <- mean(raster_espectral[[c("R", "G", "B")]])
  cobertura_sf <- sf::st_read(conaf_path, quiet = TRUE)

  if (!campo_cobertura %in% names(cobertura_sf)) {
    stop("El campo '", campo_cobertura, "' no existe en la capa de cobertura. Campos disponibles: ",
         paste(names(cobertura_sf), collapse = ", "))
  }

  # 2.2 Rasterizar la categoría de cobertura sobre la grilla del raster --------
  cobertura_sf[[campo_cobertura]] <- as.factor(cobertura_sf[[campo_cobertura]])
  cobertura_vect <- terra::vect(sf::st_transform(cobertura_sf, terra::crs(raster_espectral)))
  r_cobertura <- terra::rasterize(cobertura_vect, raster_espectral, field = campo_cobertura)

  # 2.3 Armar tabla de píxeles con cobertura asignada --------------------------
  pixeles <- c(raster_espectral, r_cobertura)
  names(pixeles)[terra::nlyr(pixeles)] <- "cobertura"
  df_pixeles <- as.data.frame(pixeles, xy = TRUE, na.rm = TRUE)
  df_pixeles$cobertura <- as.factor(df_pixeles$cobertura)
  cat("Píxeles con cobertura asignada:", nrow(df_pixeles), "de", terra::ncell(raster_espectral), "\n")
  print(table(df_pixeles$cobertura))

  # 2.4 Muestreo aleatorio estratificado ----------------------------------------
  muestra <- df_pixeles %>%
    dplyr::group_by(cobertura) %>%
    dplyr::group_modify(~ dplyr::slice_sample(.x, n = min(n_por_categoria, nrow(.x)), replace = FALSE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cobertura = droplevels(cobertura))
  cat("\nMuestra final por categoría (objetivo:", n_por_categoria, "c/u):\n")
  print(table(muestra$cobertura))

  # 2.5 Entrenar RandomForest ----------------------------------------------------
  predictoras <- names(raster_espectral)
  modelo_rf <- randomForest::randomForest(
    x = muestra[, predictoras], y = muestra$cobertura,
    importance = TRUE, ntree = 500
  )
  cat("\n---- RandomForest: matriz de confusión (OOB) ----\n")
  print(modelo_rf$confusion)
  cat("Error OOB global:", round(100 * sum(modelo_rf$predicted != muestra$cobertura) / nrow(muestra), 1), "%\n")
  cat("\n---- Importancia de variables ----\n")
  print(randomForest::importance(modelo_rf))

  # 2.6 Clasificar el raster completo con el modelo entrenado --------------------
  r_clasificado <- terra::predict(raster_espectral, modelo_rf, na.rm = TRUE)

  ruta_modelo <- file.path(carpeta_area, paste0(nombre_area, "_rf_modelo.rds"))
  ruta_muestra <- file.path(carpeta_area, paste0(nombre_area, "_muestra_estratificada.csv"))
  ruta_raster_clasif <- file.path(carpeta_area, paste0(nombre_area, "_clasificacion.tif"))

  saveRDS(modelo_rf, ruta_modelo)
  readr::write_csv(muestra, ruta_muestra)
  terra::writeRaster(r_clasificado, ruta_raster_clasif, overwrite = TRUE)

  cat("\n✔️ Modelo, muestra y raster clasificado guardados en:", carpeta_area, "\n")
  return(list(
    modelo = modelo_rf, muestra = muestra, raster_clasificado = r_clasificado,
    ruta_raster = ruta_raster_clasif, niveles_cobertura = levels(muestra$cobertura)
  ))
}
