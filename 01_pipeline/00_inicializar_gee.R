# ==========================================================================
# Script: 00_inicializar_gee.R
# Descripción: Inicializa la sesión de Google Earth Engine (rgee) para esta
# máquina. Correr una vez al abrir una sesión nueva de R/RStudio, antes de
# cualquier otro script del pipeline (rgee no mantiene la sesión entre
# reinicios de R).
# ==========================================================================

rgee::ee_Initialize(user = "tu_correo@gmail.com", drive = FALSE)  # no usamos Drive en ningún paso del pipeline

# Si el token de Earth Engine expira (error "EE credential has expired"),
# correr una vez y volver a intentar:
#   rgee::ee_clean_user_credentials(user = "tu_correo@gmail.com")
#   rgee::ee_Authenticate(user = "tu_correo@gmail.com")
