# ==========================================================================
# Script: correr_mi_area.R
# Descripción: Plantilla para correr el pipeline sobre cualquier área nueva.
# Editar SOLO la zona de AJUSTES de abajo y apretar "Source" (Ctrl+Shift+S).
# No hace falta tocar nada más del archivo.
# ==========================================================================

source("run_pipeline.R")

# --------------------------------------------------------------------------
# ZONA DE AJUSTES -- lo único que hay que cambiar para una área nueva
# --------------------------------------------------------------------------
shp_path     <- "00_input/plaza_sur.shp"   # polígono del área (.shp+.shx+.dbf+.prj en 00_input)
fecha_inicio <- "2022-01-01"
fecha_fin    <- "2026-07-31"
nombre_area  <- "plaza_sur"                         # nombre libre, se usa para nombrar la carpeta de salida
gee_user     <- "agustincoddoudiaz@gmail.com"
conaf_path   <- NULL                              # opcional: "00_input/conaf_mi_area.shp" -- si no hay capa
                                                   # de cobertura para esta área, dejar en NULL
# --------------------------------------------------------------------------

resultado <- correr_pipeline(
  shp_path     = shp_path,
  fecha_inicio = fecha_inicio,
  fecha_fin    = fecha_fin,
  nombre_area  = nombre_area,
  gee_user     = gee_user,
  conaf_path   = conaf_path
)
