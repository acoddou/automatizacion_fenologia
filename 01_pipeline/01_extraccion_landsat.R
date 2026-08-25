# ==========================================================================
# Script: 01_extraccion_landsat.R
# Descripción: Inicializa el entorno (rgee) y extrae la serie de tiempo NDVI
# de Landsat para un área de interés (shp) y rango de fechas arbitrarios.
# Adaptado de 01_main.R + 05_extraccion_series_landsat.R del repo
# residencia_aguas_de_ramon, generalizado para correr sobre cualquier área.
# ==========================================================================

# 1. Instalación y carga de paquetes ---------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse,
  terra,
  sf,
  rgee,
  readr
)

# 2. Enmascarado de nubes (compartido entre extraer_landsat() y
# extraer_raster_clasificacion()) ---------------------------------------
# Usa el bitmask QA_PIXEL de Collection 2 (bit 1 = Dilated Cloud, 2 = Cirrus,
# 3 = Cloud, 4 = Cloud Shadow, 5 = Snow). Los píxeles marcados quedan como NA.
enmascarar_nubes_l8 <- function(image) {
  qa <- image$select("QA_PIXEL")
  mascara <- qa$bitwiseAnd(bitwShiftL(1L, 1))$eq(0)$    # sin nube dilatada
    And(qa$bitwiseAnd(bitwShiftL(1L, 2))$eq(0))$        # sin cirrus
    And(qa$bitwiseAnd(bitwShiftL(1L, 3))$eq(0))$        # sin nube
    And(qa$bitwiseAnd(bitwShiftL(1L, 4))$eq(0))$        # sin sombra de nube
    And(qa$bitwiseAnd(bitwShiftL(1L, 5))$eq(0))         # sin nieve
  image$updateMask(mascara)
}

# 3. Función principal -------------------------------------------------------
#' Extrae la serie de tiempo NDVI Landsat para un área de interés
#'
#' @param shp_path Ruta al shapefile del área de interés (.shp)
#' @param fecha_inicio Fecha de inicio, formato "YYYY-MM-DD"
#' @param fecha_fin Fecha de fin, formato "YYYY-MM-DD"
#' @param nombre_area Nombre identificador del área (usado en nombres de archivo)
#' @param gee_user Usuario de Google Earth Engine (tu cuenta autenticada)
#' @param output_dir Carpeta donde se guardará el CSV crudo de salida
#' @param grid_scale Espaciamiento (m) de la grilla de puntos de muestreo. Default 30 (resolución nativa Landsat)
extraer_landsat <- function(shp_path,
                            fecha_inicio,
                            fecha_fin,
                            nombre_area,
                            gee_user,
                            output_dir = "02_output",
                            grid_scale = 30) {
  
  # 3.1 Inicializar GEE ------------------------------------------------------
  cat("Inicializando sesión de Google Earth Engine...\n")
  tryCatch({
    rgee::ee_Initialize(user = gee_user, drive = FALSE)  # no usamos Drive en ningún paso, evita ese login extra
    cat("✔️ GEE inicializado correctamente.\n")
  }, error = function(e) {
    stop("Error al inicializar GEE. Ejecuta ee_Authenticate() en la consola primero.\n", e)
  })
  
  # 3.2 Leer el área de interés (genérico, no hardcodeado) -------------------
  if (!file.exists(shp_path)) {
    stop("No se encontró el shapefile en: ", shp_path)
  }
  area <- sf::st_read(shp_path, quiet = TRUE)
  # Se descartan todos los atributos del shapefile y se queda solo la
  # geometria: Earth Engine no acepta nombres de propiedad con puntos, y
  # shapefiles con muchos atributos (ej. exportados desde OpenStreetMap)
  # suelen traer columnas como "X.id" que rompen sf_as_ee() al subirlas.
  area <- sf::st_sf(geometry = sf::st_geometry(area))
  area_ee <- rgee::sf_as_ee(area)
  
  # 3.3 Preparar carpeta de salida -------------------------------------------
  carpeta_area <- file.path(output_dir, paste0(nombre_area, "_", format(Sys.Date(), "%Y%m%d")))
  dir.create(carpeta_area, recursive = TRUE, showWarnings = FALSE)
  
  # 3.4 Definir región y colección Landsat 8 (Collection 2, Tier 1, Level 2) --
  region_of_interest <- area_ee$geometry()
  landsat8 <- ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")

  calcular_ndvi_l8 <- function(image) {
    image <- enmascarar_nubes_l8(image)
    # Escala + offset correctos para reflectancia de superficie C2 L2
    # (reflectancia = DN * 0.0000275 - 0.2; el offset no se cancela en la razón NDVI)
    nir <- image$select("SR_B5")$multiply(0.0000275)$add(-0.2)
    red <- image$select("SR_B4")$multiply(0.0000275)$add(-0.2)
    ndvi <- nir$subtract(red)$divide(nir$add(red))
    # NDVI solo está matemáticamente acotado a [-1, 1] si nir+red > 0. Un
    # píxel ruidoso que pasó el enmascarado de nubes puede tener nir+red
    # cercano a 0 (o negativo), disparando el NDVI a valores absurdos
    # (se vieron casos de hasta 12). Se enmascara como inválido en vez de
    # dejar pasar el número: es mejor un hueco que interpolar en 02, que
    # un valor imposible arrastrado por toda la pipeline.
    ndvi <- ndvi$updateMask(ndvi$gte(-1)$And(ndvi$lte(1)))$
      rename(ee$Date(image$get("system:time_start"))$format("YYYYMMdd"))
    ndvi$copyProperties(image, list("system:index", "system:time_start"))
  }
  
  ndvi_l8 <- landsat8$
    filter(ee$Filter$date(fecha_inicio, fecha_fin))$
    filter(ee$Filter$intersects(".geo", region_of_interest))$
    map(calcular_ndvi_l8)
  
  reprojected_ndvi_l8 <- ndvi_l8$map(function(image) {
    image$reproject(crs = "EPSG:4326", scale = 30)
  })
  
  ndvi_stack_l8 <- reprojected_ndvi_l8$toBands()$clip(region_of_interest)
  
  # 3.5 Generar grilla de puntos de muestreo dentro del área -----------------
  # (ee_extract necesita geometrías vectoriales; no se puede "bajar el raster
  # completo" sin pasar por Drive, así que muestreamos el área con una grilla)
  cat("Generando grilla de puntos de muestreo (resolución:", grid_scale, "m)...\n")
  grilla <- sf::st_make_grid(area, cellsize = grid_scale, what = "centers")
  grilla <- sf::st_sf(geometry = grilla)              # st_as_sf() nombraría la columna "x",
                                                       # que luego se pisaría con las coordenadas
  grilla <- sf::st_filter(grilla, area)              # solo puntos dentro del área
  grilla$id_row <- seq_len(nrow(grilla))
  coords <- sf::st_coordinates(grilla)
  grilla$x <- coords[, "X"]
  grilla$y <- coords[, "Y"]
  
  # 3.6 Extracción directa a R (sin pasar por Drive) --------------------------
  # No se usa rgee::ee_extract(): su llamada interna a reduceRegions()
  # (ee$Image$reduceRegions(image = img, ...)) está rota con las versiones
  # actuales de earthengine-api, porque el primer parámetro del método ya no
  # se llama "image" sino "self". Se llama reduceRegions() directo sobre la
  # instancia de la imagen, que funciona sin depender del nombre del parámetro.
  cat("Extrayendo serie NDVI directo a R vía reduceRegions (", nrow(grilla), "puntos)...\n")
  grilla_ee <- rgee::sf_as_ee(grilla)
  extraccion_fc <- ndvi_stack_l8$reduceRegions(
    collection = grilla_ee,
    reducer    = ee$Reducer$first(),
    scale      = 30
  )
  info <- extraccion_fc$getInfo()
  propiedades_a_df <- function(props) {
    # Un punto totalmente enmascarado (nube) en una fecha vuelve como NULL
    # (largo 0) en vez de NA; hay que igualar los largos antes de armar el
    # data.frame o as.data.frame() falla al mezclar columnas de distinto largo.
    props[sapply(props, is.null)] <- NA
    as.data.frame(props, stringsAsFactors = FALSE)
  }
  valores_ndvi <- dplyr::bind_rows(lapply(info$features, function(f) propiedades_a_df(f$properties)))

  # 3.7 Ensamblar y guardar CSV crudo (mismo formato que el resto del pipeline)
  # Formato: id_row, x, y, NDVI_<fecha1>, NDVI_<fecha2>, ...
  resultado <- valores_ndvi %>% dplyr::relocate(id_row, x, y)

  ruta_csv <- file.path(carpeta_area, paste0(nombre_area, "_landsat_ndvi_raw.csv"))
  readr::write_csv(resultado, ruta_csv)

  cat("✔️ Extracción completa. CSV guardado en:", ruta_csv, "\n")
  return(list(carpeta_area = carpeta_area, csv_path = ruta_csv, n_puntos = nrow(grilla)))
}

# 4. Función auxiliar: raster RGB + NDVI para clasificación de cobertura -----
#' Descarga la mejor escena Landsat (menor CLOUD_COVER) dentro de un rango de
#' fechas corto como raster RGB + NDVI, para usar como base de una futura
#' clasificación de cobertura (ej. cruce con una capa de CONAF).
#' No pasa por Google Drive: usa la descarga directa de rgee (getDownloadURL).
#'
#' @param shp_path Ruta al shapefile del área de interés (.shp)
#' @param fecha_inicio,fecha_fin Rango de fechas donde buscar la mejor escena
#'   (recomendado: una ventana corta, ~1 mes, al inicio del período a analizar)
#' @param nombre_area Nombre identificador del área (usado en nombres de archivo)
#' @param gee_user Usuario de Google Earth Engine (tu cuenta autenticada)
#' @param output_dir Carpeta donde se guardará el TIF de salida
#' @param scale Resolución de descarga en metros (default 30, nativa Landsat)
extraer_raster_clasificacion <- function(shp_path,
                                          fecha_inicio,
                                          fecha_fin,
                                          nombre_area,
                                          gee_user,
                                          output_dir = "02_output",
                                          scale = 30) {

  # 4.1 Inicializar GEE --------------------------------------------------------
  cat("Inicializando sesión de Google Earth Engine...\n")
  rgee::ee_Initialize(user = gee_user, drive = FALSE)  # la descarga usa getDownloadURL, no Drive

  # 4.2 Leer el área de interés ------------------------------------------------
  if (!file.exists(shp_path)) {
    stop("No se encontró el shapefile en: ", shp_path)
  }
  area <- sf::st_read(shp_path, quiet = TRUE)
  area <- sf::st_sf(geometry = sf::st_geometry(area))  # solo geometria, ver nota en extraer_landsat()
  region_of_interest <- rgee::sf_as_ee(area)$geometry()

  carpeta_area <- file.path(output_dir, paste0(nombre_area, "_", format(Sys.Date(), "%Y%m%d")))
  dir.create(carpeta_area, recursive = TRUE, showWarnings = FALSE)

  # 4.3 Elegir la escena menos nubosa del rango --------------------------------
  coleccion <- ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")$
    filter(ee$Filter$date(fecha_inicio, fecha_fin))$
    filter(ee$Filter$intersects(".geo", region_of_interest))$
    sort("CLOUD_COVER")

  n_escenas <- coleccion$size()$getInfo()
  if (n_escenas == 0) {
    stop("No se encontró ninguna escena Landsat entre ", fecha_inicio, " y ", fecha_fin,
         " para esta área. Probá un rango de fechas más amplio.")
  }
  mejor_escena <- coleccion$first()
  # system:time_start es un entero grande (epoch en milisegundos); pedirlo
  # directo con getInfo() falla en Windows por overflow de C long al
  # convertir Python->R. Se formatea como string DENTRO de Earth Engine
  # (server-side) antes de traerlo, así nunca viaja el número grande.
  fecha_escena <- as.Date(ee$Date(mejor_escena$get("system:time_start"))$format("YYYY-MM-dd")$getInfo())
  nubosidad <- mejor_escena$get("CLOUD_COVER")$getInfo()
  cat(n_escenas, "escena(s) en el rango. Usando la del", format(fecha_escena),
      "(CLOUD_COVER =", nubosidad, "%).\n")

  # 4.4 Armar RGB (reflectancia real) + NDVI, con el mismo enmascarado --------
  # de nubes que extraer_landsat() -------------------------------------------
  escena_enmascarada <- enmascarar_nubes_l8(mejor_escena)
  reflectancia <- escena_enmascarada$
    select(c("SR_B2", "SR_B3", "SR_B4", "SR_B5"))$
    multiply(0.0000275)$add(-0.2)

  ndvi <- reflectancia$normalizedDifference(c("SR_B5", "SR_B4"))$rename("NDVI")
  rgb  <- reflectancia$select(c("SR_B4", "SR_B3", "SR_B2"))$rename(c("R", "G", "B"))
  # normalizedDifference() devuelve Float64 mientras que multiply()/add() dejan
  # las bandas RGB en Float32 -- toFloat() empareja el tipo de todas las bandas
  # antes de exportar/descargar (si no, GEE tira "inconsistent types").
  raster_salida <- rgb$addBands(ndvi)$clip(region_of_interest)$toFloat()

  # 4.5 Descargar el raster directo a R (sin Drive) ----------------------------
  ruta_tif <- file.path(carpeta_area, paste0(nombre_area, "_rgb_ndvi_", format(fecha_escena, "%Y%m%d"), ".tif"))
  cat("Descargando raster (RGB + NDVI) a", ruta_tif, "...\n")

  # La descarga a veces llega truncada/corrupta (falla GDAL recién al leer
  # los píxeles, no al bajar el archivo), así que se valida leyendo todos
  # los valores y se reintenta si hace falta.
  raster_descargado <- NULL
  intento <- 0
  while (is.null(raster_descargado) && intento < 3) {
    intento <- intento + 1
    raster_descargado <- tryCatch({
      r <- rgee::ee_as_rast(
        image  = raster_salida,
        region = region_of_interest,
        dsn    = ruta_tif,
        via    = "getDownloadURL",
        scale  = scale
      )
      terra::values(r)  # fuerza a leer todos los píxeles; revienta acá si el archivo está corrupto
      r
    }, error = function(e) {
      cat("  descarga intento", intento, "falló (", conditionMessage(e), "), reintentando...\n")
      NULL
    })
  }
  if (is.null(raster_descargado)) {
    stop("No se pudo descargar un raster válido después de 3 intentos.")
  }

  # La descarga puede volver etiquetada como UTM 19N (coordenadas numéricamente
  # consistentes pero con el EPSG "equivocado" para Chile) -- reproyectar acá
  # en R (no server-side en GEE, que corrompía el GeoTIFF descargado) y
  # renombrar bandas, que ee_as_rast() tampoco conserva.
  # (no se puede sobreescribir el mismo archivo del que se está leyendo:
  # se graba a un temporal y se reemplaza.)
  raster_descargado <- terra::project(raster_descargado, "EPSG:32719")
  names(raster_descargado) <- c("R", "G", "B", "NDVI")
  ruta_temporal <- tempfile(fileext = ".tif")
  terra::writeRaster(raster_descargado, ruta_temporal, overwrite = TRUE)
  file.copy(ruta_temporal, ruta_tif, overwrite = TRUE)
  file.remove(ruta_temporal)

  cat("✔️ Raster de clasificación guardado en:", ruta_tif, "\n")
  return(list(tif_path = ruta_tif, fecha_escena = fecha_escena, cloud_cover = nubosidad))
}