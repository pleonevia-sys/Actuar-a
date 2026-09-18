alpha <- 0.05

# ==============================================================================
# 1. DETECCIÓN DE FRAUDE: PRUEBA JI-CUADRADA DE INDEPENDENCIA
# ==============================================================================

# Cargar los archivos
datos_fraude <- read.csv('retail_fraud_detection_100k.csv', 
                         encoding = "UTF-8", check.names = FALSE)
# Hacer la tabla cruzando Dispositivo contra Fraude (Valores Observados)
tabla_datos <- table(datos_fraude$device_type, datos_fraude$fraud_flag)

rownames(tabla_datos) <- c("Computadora", "Celular", "Tablet")
colnames(tabla_datos) <- c("No es Fraude", "Sí es Fraude")

# Calcular los totales
totales_filas    <- rowSums(tabla_datos) # Total por cada aparato
totales_columnas <- colSums(tabla_datos) # Total de fraudes y no fraudes
gran_total       <- sum(tabla_datos)     # Total de todas las transacciones

# Calcular los valores esperados
valores_esperados <- outer(totales_filas, totales_columnas) / gran_total

# Calcular el resultado de la prueba Ji-Cuadrada
resta_al_cuadrado <- (tabla_datos - valores_esperados)^2
resultado_final   <- sum(resta_al_cuadrado / valores_esperados)

# Calcular los grados de libertad, punto crítico y el p-value
grados_libertad <- (nrow(tabla_datos) - 1) * (ncol(tabla_datos) - 1)
punto_critico   <- qchisq(1 - alpha, df = grados_libertad)
p_valor         <- pchisq(resultado_final, df = grados_libertad, lower.tail = FALSE)

cat('=======================================================================\n',
    '                    REPORTE DE REVISIÓN DE FRAUDE                    \n',
    '=======================================================================\n\n')

cat('--- Tabla de los casos reales encontrados ---\n')
print(tabla_datos)
cat('\n-----------------------------------------------------------------------\n')

cat('--- Tabla de los casos que se esperaban tener ---\n')
print(round(valores_esperados, 2))
cat('\n-----------------------------------------------------------------------\n\n')

cat('==== Números finales de la prueba ====\n',
    'Nivel de significancia (Alpha) =               ', alpha, '\n',
    'Resultado de la prueba (Calculado) =            ', resultado_final, '\n',
    'Valor Crítico de la Tabla =                    ', punto_critico, '\n',
    'Grados de libertad =                           ', grados_libertad, '\n',
    'P-Valor                                =       ', p_valor, '\n\n',
    '=======================================================================\n')

# Decisión final explicada con palabras simples
if (resultado_final > punto_critico) {
  cat(' CONCLUSIÓN DEL ANÁLISIS:\n',
      'El tipo de aparato que usa el cliente SÍ tiene que ver con el fraude.\n',
      'Hay un patrón claro ya que el resultado superó el punto crítico.\n')
} else {
  cat(' CONCLUSIÓN DEL ANÁLISIS:\n',
      'El tipo de aparato NO influye en el fraude.\n',
      'Da igual desde qué dispositivo compren, el riesgo es el mismo.\n')
}
cat('=======================================================================\n\n')

# --- GRÁFICO DE LA DISTRIBUCIÓN NULA (PRUEBA 1) ---
# Escala fija basada en el punto crítico para ignorar el estadístico si es muy alto
limite_x <- qchisq(0.999, df = grados_libertad)
x_val <- seq(0, limite_x, length.out = 500)
y_val <- dchisq(x_val, df = grados_libertad)
plot(x_val, y_val, type = "l", lwd = 2, col = "black", main = "Distribución Chi-Cuadrada (Fraude vs Dispositivo)",
     xlab = "Valor Estadístico", ylab = "Densidad", xlim = c(0, limite_x))
# Sombrear área de rechazo
x_rechazo <- seq(punto_critico, limite_x, length.out = 100)
y_rechazo <- dchisq(x_rechazo, df = grados_libertad)
polygon(c(punto_critico, x_rechazo, limite_x), c(0, y_rechazo, 0), col = rgb(1, 0, 0, 0.3), border = NA)
# Línea del punto crítico
abline(v = punto_critico, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Distribución Nula", "Zona de Rechazo (Alpha)", "Valor Crítico"),
       col = c("black", rgb(1, 0, 0, 0.3), "red"), lty = c(1, 1, 2), lwd = 2, bty = "n")







# ==============================================================================
# 2. ABANDONO BANCARIO POR PAÍS: TABLA DE CONTINGENCIA JI-CUADRADA
# ==============================================================================

# Cargar los archivos
datos_banco <- read.csv('Churn_Modelling.csv', 
                        encoding = "UTF-8", check.names = FALSE)

# Tabla cruzada país y pertenencia al banco
tabla_churn <- table(datos_banco$Geography, datos_banco$Exited)

colnames(tabla_churn) <- c("Cliente Fiel (Se queda)", "Cliente Perdido (Churn)")

# Calcular totales
totales_filas    <- rowSums(tabla_churn)
totales_columnas <- colSums(tabla_churn)
gran_total       <- sum(tabla_churn)

# Esperado si no importará el pais
valores_esperados <- outer(totales_filas, totales_columnas) / gran_total

# Calculo del estadistico de prueba
resta_al_cuadrado <- (tabla_churn - valores_esperados)^2
resultado_final   <- sum(resta_al_cuadrado / valores_esperados)

# Sacar los grados de libertad, punto crítico y el p-valor
grados_libertad <- (nrow(tabla_churn) - 1) * (ncol(tabla_churn) - 1)
punto_critico   <- qchisq(1 - alpha, df = grados_libertad)
p_valor         <- pchisq(resultado_final, df = grados_libertad, lower.tail = FALSE)

cat('=======================================================================\n',
    '             AUDITORÍA BANCARIA: RIESGO DE ABANDONO POR PAÍS          \n',
    '=======================================================================\n\n')

cat('--- Tabla de clientes reales (Observados) ---\n')
print(tabla_churn)
cat('\n-----------------------------------------------------------------------\n')

cat('--- Tabla de lo que se esperaba (Esperados) ---\n')
print(round(valores_esperados, 2))
cat('\n-----------------------------------------------------------------------\n\n')

cat('==== Números de la prueba estadística ====\n',
    'Nivel de significancia (Alpha) =               ', alpha, '\n',
    'Resultado de la prueba (Calculado) =            ', resultado_final, '\n',
    'Valor Crítico de la Tabla =                    ', punto_critico, '\n',
    'Grados de libertad =                           ', grados_libertad, '\n',
    'P-Valor                                 =      ', p_valor, '\n\n',
    '=======================================================================\n')

# Conclusión directa con palabras sencillas
if (resultado_final > punto_critico) {
  cat(' CONCLUSIÓN DEL ANÁLISIS:\n',
      'El país de residencia SÍ influye en la decisión de abandonar el banco.\n',
      'Hay una región geográfica que está perdiendo clientes de forma alarmante.\n')
} else {
  cat(' CONCLUSIÓN DEL ANÁLISIS:\n',
      'El país de residencia NO tiene relación con que el cliente se vaya.\n',
      'La tasa de abandono se comporta igual de forma homogénea en todas partes.\n')
}
cat('=======================================================================\n\n')

# --- GRÁFICO DE LA DISTRIBUCIÓN NULA (PRUEBA 2) ---
limite_x <- qchisq(0.999, df = grados_libertad)
x_val <- seq(0, limite_x, length.out = 500)
y_val <- dchisq(x_val, df = grados_libertad)
plot(x_val, y_val, type = "l", lwd = 2, col = "black", main = "Distribución Chi-Cuadrada (Churn por País)",
     xlab = "Valor Estadístico", ylab = "Densidad", xlim = c(0, limite_x))
# Sombrear área de rechazo
x_rechazo <- seq(punto_critico, limite_x, length.out = 100)
y_rechazo <- dchisq(x_rechazo, df = grados_libertad)
polygon(c(punto_critico, x_rechazo, limite_x), c(0, y_rechazo, 0), col = rgb(1, 0, 0, 0.3), border = NA)
# Línea del punto crítico
abline(v = punto_critico, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Distribución Nula", "Zona de Rechazo (Alpha)", "Valor Crítico"),
       col = c("black", rgb(1, 0, 0, 0.3), "red"), lty = c(1, 1, 2), lwd = 2, bty = "n")




# ==============================================================================
# 3. SALUD MENTAL TRABAJOS: PRUEBA DE KRUSKAL-WALLIS (RANK)
# ==============================================================================

# Cargar el archivo
datos_salud <- read.csv('tech_mental_health_burnout.csv', 
                        encoding = "UTF-8", check.names = FALSE)

# Rankear
datos_salud$Posicion <- rank(datos_salud$stress_level)

# Separamos por tipo de trabajo
remoto     <- datos_salud$Posicion[datos_salud$work_mode == "Remote"]
hibrido    <- datos_salud$Posicion[datos_salud$work_mode == "Hybrid"]
presencial <- datos_salud$Posicion[datos_salud$work_mode == "Onsite"]

n_remoto     <- length(remoto)
n_hibrido    <- length(hibrido)
n_presencial <- length(presencial)

N_total <- n_remoto + n_hibrido + n_presencial

# Calculos:
Suma_remoto     <- sum(remoto)
Suma_hibrido    <- sum(hibrido)
Suma_presencial <- sum(presencial)

# Calculo del estadistico de prueba
parte_remoto     <- (Suma_remoto^2) / n_remoto
parte_hibrido    <- (Suma_hibrido^2) / n_hibrido
parte_presencial <- (Suma_presencial^2) / n_presencial

suma_de_partes <- parte_remoto + parte_hibrido + parte_presencial

resultado_final <- (12 / (N_total * (N_total + 1))) * suma_de_partes - (3 * (N_total + 1))

# Grados de libertad y el punto crítico (Usa distribución Ji-Cuadrada)
grados_libertad <- 2
punto_critico   <- qchisq(1 - alpha, df = grados_libertad)
p_valor         <- pchisq(resultado_final, df = grados_libertad, lower.tail = FALSE)

cat('=======================================================================\n',
    '                        SALUD MENTAL TRABAJOS                          \n',
    '=======================================================================\n\n')

cat('--- Conteo de datos por grupo ---\n',
    'Empleados en Remoto =                          ', n_remoto, '\n',
    'Empleados en Híbrido =                         ', n_hibrido, '\n',
    'Empleados en Presencial =                      ', n_presencial, '\n',
    'Total de empleados analizados =                ', N_total, '\n\n')

cat('--- Suma de posiciones obtenidas ---\n',
    'Suma Remoto =                                  ', Suma_remoto, '\n',
    'Suma Híbrido =                                 ', Suma_hibrido, '\n',
    'Suma Presencial =                              ', Suma_presencial, '\n\n')

cat('==== Números finales de la prueba hecho a mano ====\n',
    'Nivel de significancia (Alpha) =               ', alpha, '\n',
    'Resultado de la operación (Calculado) =        ', resultado_final, '\n',
    'Valor Crítico de la Tabla =                    ', punto_critico, '\n',
    'P-Valor (Probabilidad) =                       ', p_valor, '\n\n',
    '=======================================================================\n')

# Conclusión fácil de leer
if (resultado_final > punto_critico) {
  cat(' CONCLUSIÓN DEL ANÁLISIS:\n',
      'El tipo de trabajo SÍ influye de forma directa en el estrés.\n',
      'Las sumas de las posiciones son demasiado diferentes para ser una casualidad.\n')
} else {
  cat(' CONCLUSIÓN DEL ANÁLISIS:\n',
      'El tipo de trabajo NO tiene relación con el estrés.\n',
      'Las posiciones están repartidas de forma muy pareja entre todos.\n')
}
cat('=======================================================================\n\n')

# --- GRÁFICO DE LA DISTRIBUCIÓN NULA (PRUEBA 3) ---
limite_x <- qchisq(0.999, df = grados_libertad)
x_val <- seq(0, limite_x, length.out = 500)
y_val <- dchisq(x_val, df = grados_libertad)
plot(x_val, y_val, type = "l", lwd = 2, col = "black", main = "Distribución Kruskal-Wallis (Aproximación Chi-Cuadrada)",
     xlab = "Estadístico H", ylab = "Densidad", xlim = c(0, limite_x))
# Sombrear área de rechazo
x_rechazo <- seq(punto_critico, limite_x, length.out = 100)
y_rechazo <- dchisq(x_rechazo, df = grados_libertad)
polygon(c(punto_critico, x_rechazo, limite_x), c(0, y_rechazo, 0), col = rgb(1, 0, 0, 0.3), border = NA)
# Línea del punto crítico
abline(v = punto_critico, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Distribución Nula", "Zona de Rechazo (Alpha)", "Valor Crítico"),
       col = c("black", rgb(1, 0, 0, 0.3), "red"), lty = c(1, 1, 2), lwd = 2, bty = "n")




# ==============================================================================
# 4. Creditos vs Ingreso: Correlación de rangos spearman 
# ==============================================================================

# 1. Leer los datos de crédito completos
datos_credito <- read.csv('train.csv', encoding = "UTF-8", check.names = FALSE)

# Convertimos las columnas a números puros
datos_credito$Annual_Income    <- as.numeric(datos_credito$Annual_Income)
datos_credito$Outstanding_Debt <- as.numeric(datos_credito$Outstanding_Debt)

# Limpiamos las filas que no tengan datos completos
datos_limpios <- na.omit(datos_credito[, c("Annual_Income", "Outstanding_Debt")])

# 2. El truco de los rangos: Ordenar todos los datos (las 100,000 filas)
datos_limpios$Rango_Ingresos <- rank(datos_limpios$Annual_Income)
datos_limpios$Rango_Deudas   <- rank(datos_limpios$Outstanding_Debt)

# 3. Calcular la diferencia de posiciones para cada fila (d)
datos_limpios$Diferencia <- datos_limpios$Rango_Ingresos - datos_limpios$Rango_Deudas

# 4. Elevar esas diferencias al cuadrado (d^2)
datos_limpios$Diferencia_Cuadrado <- datos_limpios$Diferencia^2

# 5. Aplicar la fórmula matemática manual de Spearman para todo el bloque
Suma_d2 <- sum(datos_limpios$Diferencia_Cuadrado)
n_datos <- nrow(datos_limpios)

# Ecuación de Spearman: 1 - ( (6 * Suma_d2) / (n * (n^2 - 1)) )
coeficiente_spearman <- 1 - ((6 * Suma_d2) / (n_datos * (n_datos^2 - 1)))

# Para muestras grandes, aproximamos usando la transformación t de Student
t_calculado <- coeficiente_spearman * sqrt((n_datos - 2) / (1 - coeficiente_spearman^2))
t_critico   <- qt(1 - alpha/2, df = n_datos - 2)

cat('=======================================================================\n',
    '             AUDITORÍA FINANCIERA: CORRELACIÓN DE SPEARMAN           \n',
    '=======================================================================\n\n')

cat('--- Resumen de los datos procesados ---\n',
    'Cantidad total de clientes auditados (n) =    ', n_datos, '\n',
    'Suma de diferencias al cuadrado (Suma d^2) =   ', Suma_d2, '\n\n')

cat('==== Resultado Final del Coeficiente de Rangos ====\n',
    'Coeficiente Spearman (r_s) =                   ', coeficiente_spearman, '\n',
    'Estadístico t calculado =                     ', t_calculado, '\n',
    'Valor t crítico (Dos colas) =                 ', t_critico, '\n\n',
    '=======================================================================\n')

# Conclusión directa basada en significancia estadística formal
if (abs(t_calculado) > t_critico) {
  cat(' CONCLUSIÓN DEL ANÁLISIS:\n',
      'Existe una relación clara entre el orden de ingresos y el orden de deudas.\n',
      'El lugar que ocupa un cliente en la fila de ingresos afecta su posición en la de deudas.\n')
} else {
  cat(' CONCLUSIÓN DEL ANÁLISIS:\n',
      'No hay relación entre las posiciones.\n',
      'Estar al principio o al final de la fila de ingresos no influye en el nivel de deuda.\n')
}
cat('=======================================================================\n\n')

# --- GRÁFICO DE LA DISTRIBUCIÓN NULA (PRUEBA 4) ---
# Escala estándar fija alrededor del punto crítico de dos colas
limite_x <- t_critico * 1.5
x_val <- seq(-limite_x, limite_x, length.out = 500)
df_plot <- min(n_datos - 2, 200) 
y_val <- dt(x_val, df = df_plot)

plot(x_val, y_val, type = "l", lwd = 2, col = "black", main = "Distribución t de Student (Asociación de Spearman)",
     xlab = "Valor Estadístico t", ylab = "Densidad", xlim = c(-limite_x, limite_x))
# Sombrear dos colas de rechazo
x_izq <- seq(-limite_x, -t_critico, length.out = 100)
polygon(c(-limite_x, x_izq, -t_critico), c(0, dt(x_izq, df = df_plot), 0), col = rgb(1, 0, 0, 0.3), border = NA)
x_der <- seq(t_critico, limite_x, length.out = 100)
polygon(c(t_critico, x_der, limite_x), c(0, dt(x_der, df = df_plot), 0), col = rgb(1, 0, 0, 0.3), border = NA)
# Líneas de los puntos críticos
abline(v = c(-t_critico, t_critico), col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Distribución Nula", "Zonas de Rechazo", "Valores Críticos"),
       col = c("black", rgb(1, 0, 0, 0.3), "red"), lty = c(1, 1, 2), lwd = 2, bty = "n")









# ==============================================================================
# 5. Calificaciones antes y déspues: Prueba de Wilcoxon
# ==============================================================================

# Cargar los datos reales del archivo CSV
datos_escuela <- read.csv('student_data.csv', 
                          encoding = "UTF-8", 
                          check.names = FALSE)

# Ver los nombres de las columnas para saber cuáles usar
# En este dataset, las calificaciones están en las columnas G1, G2, G3
# Para este ejemplo, tomaremos G1 (primer periodo) y G3 (periodo final)

# Extraemos las calificaciones
Antes   <- as.numeric(datos_escuela$G1)
Despues <- as.numeric(datos_escuela$G3)

# Creamos el data frame con los datos pareados
datos_analisis <- data.frame(Antes, Despues)
datos_analisis <- na.omit(datos_analisis)  # Eliminamos filas con valores faltantes

# 2. PASO 1 DE CONOVER: Restar el 'Después' menos el 'Antes' para ver la diferencia
datos_analisis$Diferencia <- datos_analisis$Despues - datos_analisis$Antes

# 3. PASO 2 DE CONOVER: Eliminar a los alumnos que sacaron exactamente la misma nota (Diferencia = 0)
datos_filtrados <- datos_analisis[datos_analisis$Diferencia != 0, ]

# 4. PASO 3 DE CONOVER: Sacar el valor absoluto
datos_filtrados$Valor_Absoluto <- abs(datos_filtrados$Diferencia)

# 5. PASO 4 DE CONOVER: El truco del Ranking. Ordenar esas diferencias de menor a mayor
datos_filtrados$Rango <- rank(datos_filtrados$Valor_Absoluto)

# 6. PASO 5 DE CONOVER: Separar los rangos según el signo original
rangos_positivos <- datos_filtrados$Rango[datos_filtrados$Diferencia > 0]
rangos_negativos <- datos_filtrados$Rango[datos_filtrados$Diferencia < 0]

# 7. PASO 6 DE CONOVER: Sumar las posiciones
T_positiva <- sum(rangos_positivos)
T_negativa <- sum(rangos_negativos)

# El estadístico final de Conover (W)
Estadistico_W <- min(T_positiva, T_negativa)
n_final       <- nrow(datos_filtrados)

# 8. CÁLCULO DE PROBABILIDAD (Aproximación Normal según Conover)
Media_esperada      <- (n_final * (n_final + 1)) / 4
Desviacion_esperada <- sqrt((n_final * (n_final + 1) * (2 * n_final + 1)) / 24)

# Convertir nuestro resultado a un valor Z estándar y sacar punto crítico normal
Valor_Z       <- (Estadistico_W - Media_esperada) / Desviacion_esperada
punto_critico_izq <- qnorm(alpha / 2)
punto_critico_der <- -punto_critico_izq
p_valor       <- 2 * pnorm(Valor_Z)  # Dos colas

cat('=======================================================================\n',
    '       AUDITORÍA ESCOLAR: PRUEBA DE WILCOXON (RANGOS CON SIGNO)       \n',
    '=======================================================================\n\n')

cat('--- Resumen de los datos procesados ---\n',
    'Total de estudiantes analizados (n total) =    ', nrow(datos_analisis), '\n',
    'Estudiantes con cambio (diferencia ≠ 0) =      ', n_final, '\n',
    'Estudiantes que mejoraron (diferencia > 0) =   ', length(rangos_positivos), '\n',
    'Estudiantes que empeoraron (diferencia < 0) =  ', length(rangos_negativos), '\n\n')

cat('--- Resultados del conteo de Rangos ---\n',
    'Suma de rangos Positivos (Mejoraron) =         ', T_positiva, '\n',
    'Suma de rangos Negativos (Empeoraron) =         ', T_negativa, '\n',
    'Estadístico W (el menor de los dos) =           ', Estadistico_W, '\n\n')

cat('==== Números Finales de la Ecuación de Conover ====\n',
    'Nivel de significancia (Alpha) =               ', alpha, '\n',
    'Media esperada bajo H0 =                       ', Media_esperada, '\n',
    'Desviación estándar esperada =                 ', Desviacion_esperada, '\n',
    'Valor Z calculado =                            ', Valor_Z, '\n',
    'Punto Crítico Z (α/2 = 0.025) =                ±', round(punto_critico_der, 4), '\n',
    'P-Valor (Probabilidad bilateral) =             ', p_valor, '\n\n',
    '=======================================================================\n')

# Conclusión directa comparando Z contra el Punto Crítico
if (abs(Valor_Z) > punto_critico_der) {
  cat(' CONCLUSIÓN DEL ANÁLISIS:\n',
      'Se rechaza H0. Las calificaciones cambiaron de forma significativa entre exámenes.\n',
      'El valor Z (', round(Valor_Z, 4), ') cayó en la zona de rechazo (|Z| > ', round(punto_critico_der, 4), ').\n')
} else {
  cat(' CONCLUSIÓN DEL ANÁLISIS:\n',
      'No se rechaza H0. No hay cambios significativos en las calificaciones.\n',
      'Los alumnos subieron y bajaron de forma tan pareja que el rendimiento se mantiene igual.\n')
}
cat('=======================================================================\n\n')

# --- GRÁFICO DE LA DISTRIBUCIÓN NULA (PRUEBA 5) ---
# Escala estándar fija para la curva normal
limite_x <- 4
x_val <- seq(-limite_x, limite_x, length.out = 500)
y_val <- dnorm(x_val)
plot(x_val, y_val, type = "l", lwd = 2, col = "black", 
     main = "Distribución Normal Estándar Nula (Wilcoxon - Conover)",
     xlab = "Valor Estadístico Z", ylab = "Densidad", xlim = c(-limite_x, limite_x))

# Sombrear las dos colas de rechazo (alfa/2 a cada lado)
x_izq <- seq(-limite_x, punto_critico_izq, length.out = 100)
if (length(x_izq) > 1) {
  polygon(c(-limite_x, x_izq, punto_critico_izq), c(0, dnorm(x_izq), 0), 
          col = rgb(1, 0, 0, 0.3), border = NA)
}
x_der <- seq(punto_critico_der, limite_x, length.out = 100)
if (length(x_der) > 1) {
  polygon(c(punto_critico_der, x_der, limite_x), c(0, dnorm(x_der), 0), 
          col = rgb(1, 0, 0, 0.3), border = NA)
}

# Línea del estadístico Z calculado
abline(v = Valor_Z, col = "blue", lwd = 2, lty = 1)

# Líneas de los puntos críticos
abline(v = c(punto_critico_izq, punto_critico_der), col = "red", lwd = 2, lty = 2)

legend("topright", legend = c("Distribución Nula (Z)", "Zonas de Rechazo", 
                              "Puntos Críticos", "Estadístico Z calculado"),
       col = c("black", rgb(1, 0, 0, 0.3), "red", "blue"), 
       lty = c(1, 1, 2, 1), lwd = 2, bty = "n")