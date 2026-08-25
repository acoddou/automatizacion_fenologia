# ==========================================================================
# Script: 03_fenologia_phenofit.R
# Descripción: Reemplaza TIMESAT. Toma la serie NDVI regularizada y
# gap-filled de 02_limpieza_series.R, la suaviza con Savitzky-Golay, ajusta
# una curva doble-logística (Zhang et al. 2003) por temporada y extrae
# métricas fenológicas por punto de muestreo.
# Adaptado del flujo TIMESAT del repo residencia_aguas_de_ramon
# (08_metricas_modis.R / 10_metricas_landsat.R), generalizado para leer
# cualquier CSV regularizado en vez de un archivo TXT + parseo manual del
# índice temporal continuo de TIMESAT.
# ==========================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse,
  lubridate,
  phenofit,
  readr
)

# 2. Función principal -------------------------------------------------------
#' Ajusta curvas fenológicas (Savitzky-Golay + doble-logística de Zhang) y
#' extrae métricas de temporada por punto de muestreo
#'
#' Nombres y definiciones de métricas según Miranda et al. (2025, cap. 14 de
#' "Droughts in Chile"), que a su vez es la misma convención de TIMESAT:
#'   SOS/EOS: día juliano en que el NDVI sube/baja el `sos_eos_trs`*100 % del
#'     rango entre el mínimo previo/posterior y el máximo de temporada.
#'   LOS: EOS - SOS, en días. PET: día del valor máximo dentro de SOS-EOS.
#'   VSS/VES: valor de NDVI en SOS/EOS. MVA: valor máximo en la temporada.
#'   BVP: promedio entre VSS y VES (nivel base). AMP: MVA - BVP.
#'   LIN: integral (suma) del NDVI entre SOS y EOS.
#'   SP/FAL: promedio anual de NDVI en primavera (21-sep a 21-dic) y otoño
#'     (21-mar a 21-jun) -- no dependen del ajuste de curva ni de la temporada.
#'
#' @param input_csv Ruta al CSV regularizado generado por limpiar_series()
#'   (columnas: id_row, x, y, NDVI_<año>_P<periodo>...)
#' @param nombre_area Nombre identificador del área (usado en el nombre del archivo de salida)
#' @param output_dir Carpeta donde se guardará el CSV de métricas (default: misma carpeta que input_csv)
#' @param periods_per_year Periodos por año usados al regularizar la serie en 02_limpieza_series.R
#'   (tiene que coincidir con el valor usado ahí; default 24)
#' @param sos_eos_trs Umbral (fracción de la amplitud estacional) para definir SOS/EOS.
#'   Miranda et al. citan 20-40% para bosques mediterráneos/templados; default 0.2 (extremo
#'   más inclusivo del rango). Ajustar si la tesis usó un valor distinto.
#'
#' ADVERTENCIA METODOLÓGICA (ver README): las métricas de FECHA (SOS/EOS/LOS/
#' PET) tienen correlación baja o negativa entre sensores (ej. SOS
#' MODIS-Landsat: R=-0.69) y deben tratarse con cautela frente a LIN/AMP
#' (integral y amplitud), que sí son robustas entre sensores.
extraer_fenologia <- function(input_csv,
                               nombre_area,
                               output_dir = dirname(input_csv),
                               periods_per_year = 24,
                               sos_eos_trs = 0.2) {

  # 2.1 Cargar el CSV regularizado ---------------------------------------------
  if (!file.exists(input_csv)) {
    stop("No se encontró el CSV de entrada: ", input_csv)
  }
  datos <- readr::read_csv(input_csv, show_col_types = FALSE)
  columnas_periodo <- setdiff(names(datos), c("id_row", "x", "y"))
  cat("CSV regularizado cargado:", nrow(datos), "puntos,", length(columnas_periodo), "periodos.\n")

  # 2.2 Reconstruir la fecha real de cada columna ------------------------------
  # Mismo criterio que 02_limpieza_series.R usó para asignar el periodo
  # (ceiling(doy / (365.25/periods_per_year))): se toma el día calendario del
  # punto medio de cada periodo, para no arrastrar el problema de índice
  # continuo sin decodificar que tenía el flujo con TIMESAT.
  dias_por_periodo <- 365.25 / periods_per_year
  info_columnas <- tibble(columna = columnas_periodo) %>%
    mutate(
      year    = as.integer(str_extract(columna, "(?<=NDVI_)\\d{4}")),
      periodo = as.integer(str_extract(columna, "(?<=_P)\\d{2}")),
      doy_mid = round((periodo - 0.5) * dias_por_periodo),
      fecha   = as.Date(paste0(year, "-01-01")) + doy_mid - 1
    ) %>%
    arrange(fecha)

  columnas_ordenadas <- info_columnas$columna
  fechas <- info_columnas$fecha

  # 2.3 Ajustar curva y extraer fenología, punto por punto ---------------------
  cat("Ajustando curvas (SG + Zhang) para", nrow(datos), "puntos...\n")

  ajustar_un_punto <- function(id_row, valores) {
    tryCatch({
      INPUT <- phenofit::check_input(
        t = fechas, y = valores, w = rep(1, length(valores)),
        nptperyear = periods_per_year
      )
      brks <- phenofit::season_mov(
        INPUT,
        options = list(rFUN = "smooth_wSG", threshold_max = 0.1, threshold_min = 0.1)
      )
      fit <- phenofit::curvefits(INPUT, brks, options = list(methods = "Zhang", wFUN = "wTSM"))

      # Métricas de fecha complementarias que ofrece phenofit (umbrales
      # 20/50/60%, derivada, y las 4 fases de Zhang) -- se guardan aparte
      # como referencia/comparación, no son las columnas SOS/EOS "oficiales".
      pheno_fecha_ref <- phenofit::get_pheno(fit, method = "Zhang", IsPlot = FALSE)$date$Zhang %>%
        as.data.frame()

      fila_na <- tibble(
        SOS = as.Date(NA), EOS = as.Date(NA), LOS = NA_real_, PET = as.Date(NA),
        VSS = NA_real_, VES = NA_real_, MVA = NA_real_, BVP = NA_real_, AMP = NA_real_,
        LIN = NA_real_, Lder = NA_real_, Rder = NA_real_, Sinteg = NA_real_, QC_valido = FALSE
      )

      metricas <- purrr::map_dfr(seq_along(fit), function(i) {
        temporada    <- fit[[i]]
        curva_diaria <- temporada$model$Zhang$zs$iter2 %||% temporada$model$Zhang$zs[[1]]
        # temporada$tout NO son días desde 1970-01-01 (es un índice interno de
        # phenofit en otra escala) -- la fecha real y confiable de cada
        # temporada es brks$dt$beg/end, así que la curva diaria se fecha
        # anclándola al inicio de la temporada.
        fecha_curva <- seq(brks$dt$beg[i], by = "day", length.out = length(curva_diaria))

        sos_eos <- tryCatch(
          phenofit::PhenoTrs(curva_diaria, as.numeric(fecha_curva), trs = sos_eos_trs, IsPlot = FALSE),
          error = function(e) c(sos = NA_real_, eos = NA_real_)
        )
        if (anyNA(sos_eos)) return(fila_na)

        idx_sos <- which.min(abs(as.numeric(fecha_curva) - sos_eos["sos"]))
        idx_eos <- which.min(abs(as.numeric(fecha_curva) - sos_eos["eos"]))
        if (idx_eos <= idx_sos) return(fila_na)

        tramo_temporada <- curva_diaria[idx_sos:idx_eos]
        idx_pet <- idx_sos - 1 + which.max(tramo_temporada)  # pico DENTRO de sos:eos

        vss <- curva_diaria[idx_sos]
        ves <- curva_diaria[idx_eos]
        mva <- curva_diaria[idx_pet]
        bvp <- mean(c(vss, ves))                       # nivel base = promedio inicio/fin de temporada

        ventana <- max(3, round(periods_per_year / 24 * 3))  # ventana ~ 3 periodos para la pendiente
        lder <- (curva_diaria[min(idx_sos + ventana, length(curva_diaria))] - vss) / ventana
        rder <- (ves - curva_diaria[max(idx_eos - ventana, 1)]) / ventana
        # LIN/SIN se integran sobre curva_diaria (paso = 1 día), pero TIMESAT
        # reporta estas integrales en su eje temporal nativo -- pasos de
        # ~365.25/periods_per_year días cada uno, NUNCA convertidos a días
        # por el post-procesamiento de la residencia (a diferencia de
        # Begt/Endt/Length, que sí se multiplican por el largo del paso).
        # Se divide por dias_por_periodo para quedar en las mismas unidades
        # ("NDVI acumulado por periodo") y ser comparable con TIMESAT.
        dias_por_periodo_local <- 365.25 / periods_per_year
        sinteg <- sum(pmax(tramo_temporada - bvp, 0), na.rm = TRUE) / dias_por_periodo_local

        tibble(
          SOS = fecha_curva[idx_sos], EOS = fecha_curva[idx_eos],
          LOS = as.numeric(fecha_curva[idx_eos] - fecha_curva[idx_sos]),
          PET = fecha_curva[idx_pet],
          VSS = vss, VES = ves, MVA = mva, BVP = bvp, AMP = mva - bvp,
          LIN = sum(tramo_temporada, na.rm = TRUE) / dias_por_periodo_local,
          Lder = lder, Rder = rder, Sinteg = sinteg,
          # NDVI no puede matemáticamente superar [-1, 1]; un ajuste que se
          # dispara fuera de ese rango es un artefacto numérico del fit
          # (ej. temporadas cortas/truncadas en el borde de la serie), no un
          # valor real. Se marca en vez de descartar silenciosamente.
          QC_valido = abs(mva) <= 1 & abs(bvp) <= 1
        )
      })

      # SP/FAL: promedio anual de NDVI en primavera (21-sep a 21-dic) y otoño
      # (21-mar a 21-jun), directo de la serie observada -- no dependen de la
      # temporada ni del ajuste de curva, se asignan por año calendario.
      anios <- unique(year(fechas))
      sp_fal <- purrr::map_dfr(anios, function(a) {
        tibble(
          anio = a,
          SP  = mean(valores[fechas >= as.Date(sprintf("%d-09-21", a)) & fechas <= as.Date(sprintf("%d-12-21", a))], na.rm = TRUE),
          FAL = mean(valores[fechas >= as.Date(sprintf("%d-03-21", a)) & fechas <= as.Date(sprintf("%d-06-21", a))], na.rm = TRUE)
        )
      })

      metricas <- metricas %>% dplyr::mutate(anio = year(SOS)) %>%
        dplyr::left_join(sp_fal, by = "anio") %>% dplyr::select(-anio)

      dplyr::bind_cols(id_row = id_row, season_n = seq_len(nrow(pheno_fecha_ref)),
                        flag = pheno_fecha_ref$flag, metricas, pheno_fecha_ref %>% dplyr::select(-flag))
    }, error = function(e) {
      cat("  ! Punto", id_row, "falló:", conditionMessage(e), "\n")
      NULL
    })
  }

  resultado_lista <- purrr::pmap(
    list(datos$id_row, split(as.matrix(datos[, columnas_ordenadas]), seq_len(nrow(datos)))),
    ajustar_un_punto
  )

  fenologia <- dplyr::bind_rows(resultado_lista)

  # 2.4 Unir coordenadas y exportar --------------------------------------------
  fenologia <- datos %>%
    dplyr::select(id_row, x, y) %>%
    dplyr::right_join(fenologia, by = "id_row")

  ruta_csv <- file.path(output_dir, paste0(nombre_area, "_fenologia.csv"))
  readr::write_csv(fenologia, ruta_csv)

  cat("✔️ Fenología completa. CSV guardado en:", ruta_csv, "\n")
  cat("   Puntos procesados:", length(unique(fenologia$id_row)), "de", nrow(datos), "\n")
  return(list(csv_path = ruta_csv, n_puntos = length(unique(fenologia$id_row))))
}
