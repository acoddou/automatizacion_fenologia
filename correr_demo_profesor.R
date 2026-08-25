# ==========================================================================
# Script: correr_demo_profesor.R
# Descripción: Corre el pipeline dos veces para la demo con el profesor guía:
# primero Plaza Sur (área chica, sale rápido -- ~1.5 min, para mostrar algo
# terminado altiro) y después Apoquindo (área más grande, queda corriendo
# de fondo en la consola mientras se sigue conversando).
# Abrir este archivo y apretar "Source" (o Ctrl+Shift+S) -- no hace falta
# tocar la consola.
# ==========================================================================

source("run_pipeline.R")

# fecha real del incendio que afectó el sector (confirmada, no la de diciembre
# 2025/2024 asumida al principio) -- se marca en el gráfico de evolución de
# métricas (06) para que la caída post-incendio no se confunda con ruido
fecha_incendio <- "2025-09-29"

# -- 1) Plaza Sur -- rápido, para mostrar un resultado completo altiro -------
resultado_plaza <- correr_pipeline(
  shp_path     = "00_input/plaza_sur.shp",
  fecha_inicio = "2022-01-01",
  fecha_fin    = "2026-07-31",
  nombre_area  = "plaza_sur_demo",
  gee_user     = "agustincoddoudiaz@gmail.com",
  conaf_path   = NULL,   # Plaza Sur no tiene recorte CONAF propio todavía
  fecha_evento = fecha_incendio, evento_label = "Incendio"
)

# -- 2) Apoquindo -- más grande, queda corriendo mientras se conversa -------
resultado_apoquindo <- correr_pipeline(
  shp_path     = "00_input/limite_apoquindo.shp",
  fecha_inicio = "2022-01-01",
  fecha_fin    = "2026-07-31",
  nombre_area  = "apoquindo_demo",
  gee_user     = "agustincoddoudiaz@gmail.com",
  conaf_path   = "00_input/conaf_apoquindo.shp",
  fecha_evento = fecha_incendio, evento_label = "Incendio"
)
