# ==========================================================================
# Script: 06_visualizaciones.R
# Descripción: Genera los mapas y gráficos "de prueba de vida" del pipeline
# a partir de los resultados de correr_pipeline() (run_pipeline.R): mapa
# NDVI de la escena más reciente, serie de tiempo de NDVI de toda el área,
# y mapas de las métricas fenológicas principales (LIN, AMP, SOS). No
# requiere una capa de cobertura -- si se entrega una (conaf_path), agrega
# además la serie de tiempo y el violín/boxplot por cobertura.
# ==========================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, terra, sf, zoo, lubridate)

#' Genera los mapas y gráficos principales a partir de un resultado de
#' correr_pipeline()
#'
#' @param resultado_pipeline La lista devuelta por correr_pipeline() (tiene
#'   $extraccion$csv_path, $limpieza$csv_path, $fenologia$csv_path)
#' @param nombre_area Nombre identificador del área (para títulos)
#' @param conaf_path Ruta a una capa de cobertura (opcional). Si se entrega,
#'   agrega serie de tiempo y violín/boxplot por cobertura
#' @param campo_cobertura Campo de la capa de cobertura a usar (default "USO")
#' @param output_dir Carpeta donde guardar los PNG (default: la misma del pipeline)
generar_visualizaciones <- function(resultado_pipeline,
                                     nombre_area,
                                     conaf_path = NULL,
                                     campo_cobertura = "USO",
                                     output_dir = NULL,
                                     fecha_fin = NULL,
                                     fecha_evento = NULL,
                                     evento_label = "Evento") {

  fecha_corrida <- format(Sys.Date(), "%Y-%m-%d")
  etiqueta_run <- sprintf("Extracción automatizada — pipeline automatizacion_fenologia, corrida %s", fecha_corrida)

  if (is.null(output_dir)) output_dir <- dirname(resultado_pipeline$fenologia$csv_path)
  carpeta_plots <- file.path(output_dir, "plots")
  dir.create(carpeta_plots, recursive = TRUE, showWarnings = FALSE)

  raw <- read_csv(resultado_pipeline$extraccion$csv_path, show_col_types = FALSE)
  limpio <- read_csv(resultado_pipeline$limpieza$csv_path, show_col_types = FALSE)
  feno <- read_csv(resultado_pipeline$fenologia$csv_path, show_col_types = FALSE)

  paleta <- hcl.colors(50, "viridis")

  # 1. Mapa NDVI de la escena mas reciente con datos completos -----------------
  cat("[1/4] Mapa NDVI...\n")
  columnas_fecha <- names(raw)[-(1:3)]
  fechas <- as.Date(str_extract(columnas_fecha, "\\d{8}$"), format = "%Y%m%d")
  cobertura_n <- sapply(raw[columnas_fecha], function(col) sum(!is.na(col)))
  bien_cubiertas <- cobertura_n >= 0.9 * nrow(raw)
  col_reciente <- columnas_fecha[bien_cubiertas][which.max(fechas[bien_cubiertas])]
  fecha_usada <- fechas[columnas_fecha == col_reciente]

  d <- raw %>% select(id_row, x, y, ndvi = all_of(col_reciente))
  puntos_v <- vect(d, geom = c("x", "y"), crs = "EPSG:32719")
  r_template <- rast(puntos_v, resolution = 150)
  r_ndvi <- rasterize(puntos_v, r_template, field = "ndvi", fun = "mean")

  png(file.path(carpeta_plots, "01_mapa_ndvi.png"), width = 1000, height = 850, res = 120, bg = "white")
  plot(r_ndvi, col = paleta, range = c(-0.2, 0.8),
       main = paste0("NDVI Landsat — ", nombre_area, " (", format(fecha_usada, "%d-%b-%Y"), ")"))
  mtext(etiqueta_run, side = 1, line = 3.5, cex = 0.7, col = "grey30", font = 3)
  dev.off()

  # 2. Serie de tiempo NDVI general (sin distincion de cobertura) --------------
  cat("[2/4] Serie de tiempo NDVI general...\n")
  columnas_periodo <- setdiff(names(limpio), c("id_row", "x", "y"))
  dias_por_periodo <- 365.25 / 24
  info_col <- tibble(columna = columnas_periodo) %>%
    mutate(year = as.integer(str_extract(columna, "(?<=NDVI_)\\d{4}")),
           periodo = as.integer(str_extract(columna, "(?<=_P)\\d{2}")),
           doy_mid = round((periodo - 0.5) * dias_por_periodo),
           fecha = as.Date(paste0(year, "-01-01")) + doy_mid - 1)

  resumen_general <- limpio %>%
    pivot_longer(all_of(columnas_periodo), names_to = "columna", values_to = "ndvi") %>%
    left_join(info_col, by = "columna") %>%
    group_by(fecha) %>%
    summarise(ndvi_medio = mean(ndvi, na.rm = TRUE), .groups = "drop") %>%
    arrange(fecha) %>%
    mutate(ndvi_suave = zoo::rollmean(ndvi_medio, k = 3, fill = NA, align = "center"))

  p_serie <- ggplot(resumen_general, aes(x = fecha, y = ndvi_medio)) +
    geom_line(color = "grey60", alpha = 0.5, linewidth = 0.4) +
    geom_line(aes(y = ndvi_suave), color = "#1a7a3c", linewidth = 1.1) +
    labs(title = paste0("Serie de tiempo de NDVI (Landsat) — ", nombre_area),
         subtitle = etiqueta_run, x = NULL, y = "NDVI promedio") +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle = element_text(size = 8, color = "grey40", face = "italic"))
  ggsave(file.path(carpeta_plots, "02_serie_tiempo_ndvi.png"), p_serie, width = 12, height = 5.5, dpi = 120, bg = "white")

  # 3. Mapas de metricas fenologicas (LIN, AMP, SOS) ----------------------------
  cat("[3/4] Mapas de métricas fenológicas...\n")
  feno_1ra_temporada <- feno %>% filter(season_n == 1, QC_valido == TRUE) %>%
    mutate(SOS_doy = yday(SOS))
  puntos_feno <- vect(feno_1ra_temporada, geom = c("x", "y"), crs = "EPSG:32719")
  r_template_feno <- rast(puntos_feno, resolution = 150)

  png(file.path(carpeta_plots, "03_mapas_fenologia.png"), width = 1500, height = 550, res = 110, bg = "white")
  par(mfrow = c(1, 3), mar = c(2, 2, 3, 4))
  for (metrica in c("LIN", "AMP", "SOS_doy")) {
    r <- rasterize(puntos_feno, r_template_feno, field = metrica, fun = "mean")
    plot(r, main = paste0(metrica, " (1ra temporada)"), col = paleta)
  }
  dev.off()

  # 4. (opcional) serie de tiempo y violin/boxplot por cobertura -----------------
  archivos_generados <- c("01_mapa_ndvi.png", "02_serie_tiempo_ndvi.png", "03_mapas_fenologia.png")

  puntos_cobertura <- NULL
  if (!is.null(conaf_path)) {
    cat("[4/5] Serie de tiempo y violín/boxplot por cobertura...\n")
    conaf <- st_read(conaf_path, quiet = TRUE)
    puntos_sf <- st_as_sf(limpio %>% distinct(id_row, x, y), coords = c("x", "y"), crs = 32719)
    puntos_cobertura <- st_join(puntos_sf, conaf[campo_cobertura]) %>%
      st_drop_geometry() %>%
      rename(cobertura = all_of(campo_cobertura)) %>%
      filter(!is.na(cobertura), cobertura != "Áreas Urbanas e Industriales")

    serie_cob <- limpio %>%
      pivot_longer(all_of(columnas_periodo), names_to = "columna", values_to = "ndvi") %>%
      left_join(info_col, by = "columna") %>%
      left_join(puntos_cobertura, by = "id_row") %>%
      filter(!is.na(cobertura)) %>%
      group_by(cobertura, fecha) %>%
      summarise(ndvi_medio = mean(ndvi, na.rm = TRUE), .groups = "drop") %>%
      arrange(cobertura, fecha) %>%
      group_by(cobertura) %>%
      mutate(ndvi_suave = zoo::rollmean(ndvi_medio, k = 3, fill = NA, align = "center")) %>%
      ungroup()

    p_cob <- ggplot(serie_cob, aes(x = fecha, y = ndvi_medio, color = cobertura)) +
      geom_line(alpha = 0.25, linewidth = 0.4) +
      geom_line(aes(y = ndvi_suave), linewidth = 1) +
      labs(title = paste0("NDVI por cobertura — ", nombre_area), subtitle = etiqueta_run,
           x = NULL, y = "NDVI promedio", color = "Cobertura") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom", plot.subtitle = element_text(size = 8, color = "grey40", face = "italic"))
    ggsave(file.path(carpeta_plots, "04_serie_tiempo_por_cobertura.png"), p_cob, width = 12, height = 5.5, dpi = 120, bg = "white")

    feno_cob <- feno %>% filter(QC_valido == TRUE) %>% inner_join(puntos_cobertura, by = "id_row")
    p_violin <- ggplot(feno_cob, aes(x = cobertura, y = LIN, fill = cobertura)) +
      geom_violin(alpha = 0.5, trim = FALSE) +
      geom_boxplot(width = 0.15, outlier.size = 0.8, alpha = 0.8) +
      labs(title = paste0("Productividad acumulada (LIN) por cobertura — ", nombre_area),
           subtitle = etiqueta_run, x = "Cobertura", y = "NDVI acumulado / periodo") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none", plot.subtitle = element_text(size = 8, color = "grey40", face = "italic"))
    ggsave(file.path(carpeta_plots, "05_violin_LIN_por_cobertura.png"), p_violin, width = 8, height = 5.5, dpi = 120, bg = "white")

    archivos_generados <- c(archivos_generados, "04_serie_tiempo_por_cobertura.png", "05_violin_LIN_por_cobertura.png")
  }

  # 5. Evolución interanual de las métricas fenológicas (SOS, LOS, AMP, LIN) ----
  # Esto es lo que distingue el análisis de una serie de NDVI cruda: no solo
  # el promedio agregado por cobertura (el violín de arriba), sino cómo se
  # mueve cada métrica año a año -- la pregunta que le importa a un lector de
  # fenología (¿se adelanta el inicio de temporada? ¿cae la amplitud?).
  cat("[5/5] Evolución interanual de métricas fenológicas...\n")

  feno_anual <- feno %>%
    filter(QC_valido == TRUE) %>%
    mutate(anio = lubridate::year(SOS), SOS_doy = lubridate::yday(SOS))

  # El detector de temporadas de phenofit a veces marca una segunda
  # "temporada" espuria de amplitud casi nula dentro del valle entre ciclos
  # (ruido residual, no un segundo crecimiento real en un clima mediterráneo
  # de una sola estación) -- se deja solo la de mayor amplitud por punto y año.
  feno_anual <- feno_anual %>%
    group_by(id_row, anio) %>%
    slice_max(AMP, n = 1, with_ties = FALSE) %>%
    ungroup()

  if (!is.null(fecha_fin)) {
    # una temporada cuyo EOS cae después del fin del rango pedido está
    # trunca (el pipeline no vio el cierre real del ciclo) -- AMP/LOS/LIN
    # de esa temporada quedan artificialmente bajos y no son comparables
    # con las demás, así que se excluyen del gráfico de tendencia. OJO: si el
    # rango pedido no llega a fin de año, el último año queda con menos
    # puntos (solo los que arrancaron temprano alcanzan a cerrar antes del
    # corte) -- eso es un sesgo de muestra real, pero NO hay que descartar el
    # año entero por eso: si ese año coincide con un evento real (ej. un
    # incendio), la caída de amplitud/LIN puede ser la señal más importante
    # del gráfico, no un artefacto. Se deja visible; el corte de muestra se
    # menciona como caveat aparte, no se oculta el dato.
    feno_anual <- feno_anual %>% filter(EOS <= as.Date(fecha_fin))
  }

  # una temporada de menos de ~90 días no es un ciclo fenológico plausible
  # para esta vegetación (las temporadas reales rondan 120-220 días) -- casi
  # siempre es un corte artificial contra el borde del rango de datos pedido
  feno_anual <- feno_anual %>% filter(LOS >= 90)

  if (!is.null(puntos_cobertura)) {
    feno_anual <- feno_anual %>% inner_join(puntos_cobertura, by = "id_row")
  } else {
    feno_anual <- feno_anual %>% mutate(cobertura = "Área completa")
  }

  feno_anual <- feno_anual %>% mutate(sos_decimal = lubridate::decimal_date(SOS))

  resumen_metricas <- feno_anual %>%
    group_by(anio, cobertura) %>%
    summarise(
      sos_decimal_media = mean(sos_decimal, na.rm = TRUE),
      across(c(SOS_doy, LOS, AMP, LIN),
             list(media = ~mean(.x, na.rm = TRUE), se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))),
             .names = "{.col}__{.fn}"),
      n = dplyr::n(), .groups = "drop"
    ) %>%
    filter(n >= 3)  # años/coberturas con muy pocos puntos no aportan una media confiable

  largo_metricas <- resumen_metricas %>%
    pivot_longer(cols = matches("__(media|se)$"),
                 names_to = c("metrica", ".value"), names_pattern = "(.*)__(media|se)$") %>%
    mutate(metrica = factor(metrica, levels = c("SOS_doy", "LOS", "AMP", "LIN"),
      labels = c("Inicio de temporada (SOS, día del año)", "Duración de temporada (LOS, días)",
                 "Amplitud estacional (AMP)", "Productividad acumulada (LIN)")))

  p_metricas <- ggplot(largo_metricas, aes(x = sos_decimal_media, y = media, color = cobertura, fill = cobertura)) +
    geom_ribbon(aes(ymin = media - se, ymax = media + se), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.8)

  if (!is.null(fecha_evento)) {
    # linea de referencia para un evento puntual conocido (ej. un incendio):
    # sin esto, una caida real coincidente con el evento es indistinguible a
    # simple vista de una caida "normal" de otro año
    evento_decimal <- lubridate::decimal_date(as.Date(fecha_evento))
    p_metricas <- p_metricas +
      geom_vline(xintercept = evento_decimal, linetype = "dashed", color = "grey30", linewidth = 0.6) +
      geom_text(data = distinct(largo_metricas, metrica),
                aes(x = evento_decimal, y = Inf, label = evento_label), inherit.aes = FALSE,
                vjust = 1.4, hjust = -0.05, size = 2.8, color = "grey30", fontface = "italic")
  }

  p_metricas <- p_metricas +
    facet_wrap(~metrica, scales = "free_y", ncol = 2) +
    labs(title = paste0("Evolución interanual de métricas fenológicas — ", nombre_area),
         subtitle = etiqueta_run, x = NULL, y = NULL, color = "Cobertura", fill = "Cobertura") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(size = 8, color = "grey40", face = "italic"),
          strip.text = element_text(face = "bold"))

  ggsave(file.path(carpeta_plots, "06_series_metricas_fenologicas.png"), p_metricas,
         width = 11, height = 8, dpi = 120, bg = "white")

  archivos_generados <- c(archivos_generados, "06_series_metricas_fenologicas.png")

  cat("\n✔️ Visualizaciones guardadas en:", carpeta_plots, "\n")
  invisible(file.path(carpeta_plots, archivos_generados))
}
