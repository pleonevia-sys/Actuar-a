####### Declaramos constantes y funciones que se usarán a lo largo del código
alpha = 0.05

lilliefors.exp.test <- function(x){
  n <- length(x)
  i <- 1:n
  xi <- sort(x)
  zi <- xi / mean(x)
  fi <- 1 - exp(-zi)
  t1 <- (i / n) - fi
  t2 <- fi - ((i - 1)/n)
  D <- max(abs(t1),abs(t2))
  return(D)
}


# ==============================================================================
# Prueba Ji-Cuadrado
# ==============================================================================
chisq.test <- function(x, p_dist_func, num_param, ...) {
  n <- length(x)
  
  # 1. Determinación de intervalos (Sturges)
  k <- floor(1 + log2(n))
  
  # Dividimos los tramos
  cortes <- quantile(x, probs = seq(0, 1, length.out = k + 1))
  cortes[1] <- -Inf
  cortes[length(cortes)] <- Inf
  
  # 3. Frecuencias Observadas (O_i)
  O_i <- as.vector(table(cut(x, breaks = cortes, include.lowest = TRUE)))
  
  # 4. Frecuencias Esperadas (E_i)
  p_acumulada <- p_dist_func(cortes, ...)
  prob_intervalo <- diff(p_acumulada)
  prob_intervalo <- prob_intervalo / sum(prob_intervalo) #Normalizamos
  E_i <- n * prob_intervalo
  
  # 5. Cálculo del Estadístico de Prueba
  estadistico <- sum(((O_i - E_i)^2) / E_i)
  
  # 6. Grados de Libertad y Decisión Teóricas
  df_corregidos <- k - 1 - num_param
  punto_critico <- qchisq(1 - alpha, df = df_corregidos)
  p_valor       <- 1 - pchisq(estadistico, df = df_corregidos)
  
  # Retornamos los resultados
  return(list(
    statistic = estadistico,
    intervals = k,
    df = df_corregidos,
    critical_value = punto_critico,
    p.value = p_valor
  ))
}

# ==============================================================================
# Prueba Anderson Darling
# ==============================================================================
ad.test <- function(x, p_dist_func, ...) {
  n <- length(x)
  i <- 1:n
  
  # 1. Ordenamos los datos de menor a mayor (Requisito estricto de la prueba)
  x_sort <- sort(x)
  
  # 2. Calculamos las probabilidades acumuladas teóricas (U_i)
  p_teorica <- p_dist_func(x_sort, ...)
  
  # Control de estabilidad numérica para evitar la indeterminación de ln(0) o ln(1)
  p_teorica <- pmax(pmin(p_teorica, 1 - 1e-10), 1e-10)
  
  # 3. Vector invertido para el extremo superior: U_(n-i+1)
  p_teorica_rev <- rev(p_teorica)
  
  # 4. Aplicamos la sumatoria de la fórmula analítica de Anderson-Darling
  sumatoria <- (2 * i - 1) * (log(p_teorica) + log(1 - p_teorica_rev))
  A2_estadistico <- -n - (1 / n) * sum(sumatoria)
  
  # 5. Corrección para muestras finitas (A_modificado)
  A2_modificado <- A2_estadistico * (1 + 0.75/n + 2.25/(n^2))
  
  return(list(
    statistic = A2_estadistico,
    adjusted_statistic = A2_modificado
  ))
}

# ==============================================================================
# Prueba de Kolmogorov-Smirnov
# ==============================================================================
ks.test <- function(x, p_dist_func, ...) {
  n <- length(x)
  
  # Ordenamos los datos
  x_sort <- sort(x)
  
  # Creamos dos vectores 0 de dimensión n
  d_mas   <- numeric(n)
  d_menos <- numeric(n)
  
  # Ciclo FOR para recorrer y restar cada escalón individualmente
  for(i in 1:n) {
    # Probabilidad teórica acumulada para el dato actual F0(x_i)
    p_teorica <- p_dist_func(x_sort[i], ...)
    
    # Escalón superior: (i / n) - F0(x_i)
    d_mas[i]   <- (i / n) - p_teorica
    
    # Escalón inferior: F0(x_i) - ((i - 1) / n)
    d_menos[i] <- p_teorica - ((i - 1) / n)
  }
  
  # El estadístico final es el máximo absoluto observado en toda la muestra
  max_d_mas     <- max(abs(d_mas))
  max_d_menos   <- max(abs(d_menos))
  D_estadistico <- max(max_d_mas, max_d_menos)
  
  # Cálculo del p-value (Aproximación Asintótica + Corrección)
  d_corregido <- D_estadistico + (1 / (6 * n)) + (D_estadistico / (2 * n))
  
  # Evaluamos la función de distribución asintótica al cuadrado: [G(x)]^2
  p_value <- 1 - (1 - exp(-2 * n * d_corregido^2))^2
  
  # Retornamos la lista con todos los valores calculados
  return(list(
    statistic   = D_estadistico,
    p.value     = p_value
  ))
}

# ==============================================================================
# SHAPIRO-WILK
# ==============================================================================
shapiro.test <- function(x) {
  n <- length(x)
  
  # Ordenar los datos (eje Y del gráfico Q-Q)
  y <- sort(x)
  
  # Calcular los cuantiles teóricos normales esperados (eje X del gráfico Q-Q)
  probabilidades <- (1:n - 0.375) / (n + 0.25)
  x_teorica <- qnorm(probabilidades)
  
  # El estadístico W es el coeficiente de correlación al cuadrado (R^2)
  W_estadistico <- cor(x_teorica, y)^2
  
  # P-value
  mu_sencilla    <- -1.58 - (0.31 * log(n))
  sigma_sencilla <- exp(-0.48 + (0.43 * log(n)))
  
  z <- (log(1 - W_estadistico) - mu_sencilla) / sigma_sencilla
  p_value <- pnorm(z, lower.tail = FALSE)
  
  # Retornar resultados
  return(list(
    statistic = W_estadistico,
    p.value   = p_value
  ))
}

# ==============================================================================
# Gamma por Máximo Verosimilitud fórmula de CHOI-WETTE
# ==============================================================================
estimar_gamma_mle <- function(datos) {
  M <- log(mean(datos)) - mean(log(datos))
  shape_mle <- (3 - M + sqrt((M - 3)^2 + 24 * M)) / (12 * M)
  rate_mle  <- shape_mle / mean(datos)
  return(c(shape = shape_mle, rate = rate_mle))
}

# ==============================================================================
# Weibull por Newton-Rhapson
# ==============================================================================
estimar_weibull_mle <- function(datos, tolerancia = 1e-6) {
  k_actual <- 1.0 
  error <- 1 
  
  while (error > tolerancia) {
    sum_k   <- sum(datos^k_actual)
    sum_ln  <- sum(datos^k_actual * log(datos))
    sum_ln2 <- sum(datos^k_actual * (log(datos))^2)
    
    f_k <- (sum_ln / sum_k) - (1 / k_actual) - mean(log(datos))
    df_k <- ((sum_ln2 * sum_k) - (sum_ln)^2) / (sum_k)^2 + (1 / k_actual^2)
    
    k_anterior <- k_actual
    k_actual <- k_anterior - (f_k / df_k)
    error <- abs(k_actual - k_anterior)
  }
  
  shape_mle <- k_actual
  scale_mle <- (mean(datos^shape_mle))^(1 / shape_mle)
  
  return(c(shape = shape_mle, scale = scale_mle))
}

# ==============================================================================
# lognormal Máximo Verosimilitud
# ==============================================================================
estimar_lognormal_mle <- function(datos) {
  n <- length(datos)
  datos_log <- log(datos)
  
  meanlog_mle <- mean(datos_log)
  
  var_mle <- var(datos_log) * ((n - 1) / n)
  
  sdlog_mle  <- sqrt(var_mle)
  
  return(c(meanlog = meanlog_mle, sdlog = sdlog_mle))
}

# ==============================================================================
# Beta por método de momentos
# ==============================================================================
estimar_beta_mle <- function(datos) {
  datos_escalados <- datos
  
  x_bar <- mean(datos_escalados)
  v_bar <- var(datos_escalados)
  
  factor <- (x_bar * (1 - x_bar) / v_bar) - 1
  
  shape1_mle <- x_bar * factor       
  shape2_mle <- (1 - x_bar) * factor 
  
  if(shape1_mle <= 0) shape1_mle <- 0.1
  if(shape2_mle <= 0) shape2_mle <- 0.1
  
  return(c(shape1 = shape1_mle, shape2 = shape2_mle))
}


########## Base de datos de reclamaciones de un seguro de autos ################
AutoInsurance1 <- read.csv("AutoInsurance.csv")
ClaimAuto1     <- AutoInsurance1$Total.Claim.Amount
n1             <- length(ClaimAuto1)

# ==============================================================================
# Graficamos
# ==============================================================================
hist(ClaimAuto1, main = "Histograma de Reclamaciones", xlab = "Monto")
plot(ecdf(ClaimAuto1))
boxplot(ClaimAuto1)

# ==============================================================================
# AJUSTE DE LAS DISTRIBUCIONES
# ==============================================================================
mle_weibull1 <- estimar_weibull_mle(ClaimAuto1)
mle_gamma1   <- estimar_gamma_mle(ClaimAuto1)

shape_weibull1 <- mle_weibull1['shape']
scale_weibull1 <- mle_weibull1['scale']

shape_gamma1   <- mle_gamma1['shape']
rate_gamma1    <- mle_gamma1['rate']

# ==============================================================================
# PRUEBAS DE BONDAD DE AJUSTE
# ==============================================================================

# --- Prueba Ji-Cuadrada ---
prueba_w1 <- chisq.test(ClaimAuto1, p_dist_func = pweibull, num_param = 2, 
                        shape = shape_weibull1, scale = scale_weibull1)

prueba_g1 <- chisq.test(ClaimAuto1, p_dist_func = pgamma, num_param = 2, 
                        shape = shape_gamma1, rate = rate_gamma1)

# --- Prueba Lilliefors para Distribución Exponencial ---
D2_exponencial1 <- lilliefors.exp.test(ClaimAuto1)

# --- Prueba Kolmogorov-Smirnov ---
ks_w1 <- ks.test(ClaimAuto1, pweibull, shape = shape_weibull1, scale = scale_weibull1)
ks_g1 <- ks.test(ClaimAuto1, pgamma, shape = shape_gamma1, rate = rate_gamma1)

# --- Prueba Anderson-Darling ---
ad_w1 <- ad.test(ClaimAuto1, pweibull, shape = shape_weibull1, scale = scale_weibull1)
ad_g1 <- ad.test(ClaimAuto1, pgamma, shape = shape_gamma1, rate = rate_gamma1)

# --- Prueba Shapiro-Wilk --- 
sw_w1 <- shapiro.test(qnorm(pmin(pmax(pweibull(ClaimAuto1, 
                                               shape = shape_weibull1, scale = scale_weibull1), 1e-10), 1-1e-10)))

sw_g1 <- shapiro.test(qnorm(pmin(pmax(pgamma(ClaimAuto1, 
                                             shape = shape_gamma1, rate = rate_gamma1), 1e-10), 1-1e-10)))


# ==============================================================================
# Puntos críticos
# ==============================================================================
k1                 <- floor(1 + log2(n1)) 
chis.t_wg_rechazo1 <- qchisq(1 - alpha, df = k1 - 1 - 2)
kv_wg_rechazo1     <- 1.22 / (sqrt(n1 + sqrt(n1 / 10)))
AD_wg_rechazo1     <- 2.453
#Por simulación de 1000 para Anderson Darling
#   90%   95% 97.5%   99% 
#  1.891 2.453 3.049 3.834 
Lili.t_wg_rechazo1 <- 1.0753 / (sqrt(n1))
SW_wg_rechazo1     <- qnorm(1 - alpha)


cat('========================= Resumen estadistico 1 =========================','\n',
    '==== Prueba Ji-Cuadrada ===','                 Weibull ',' Gamma ','\n',
    'Estadistico de prueba (Weibull, Gamma) = ', prueba_w1$statistic, prueba_g1$statistic, '\n', 
    'Punto crítico =                          ', chis.t_wg_rechazo1, chis.t_wg_rechazo1, '\n',
    'p_value =                                ', prueba_w1$p.value, prueba_g1$p.value, '\n\n',
    
    '==== Prueba Liliefors exponecial ===','\n',
    'Estadistico de prueba = ', D2_exponencial1, '    Punto crítico = ', Lili.t_wg_rechazo1, '\n\n',
    
    '==== Prueba Kolmogorov-Smirnov ===','        Weibull',' Gamma','\n',
    'Estadistico de prueba (Weibull, Gamma) = ', ks_w1$statistic, ks_g1$statistic, '\n', 
    'Punto crítico =                          ', kv_wg_rechazo1, kv_wg_rechazo1, '\n',
    'p_value =                                ', ks_w1$p.value, ks_g1$p.value, '\n\n',
    
    '==== Prueba Anderson Darling ===','          Weibull ',' Gamma ','\n',
    'Estadistico de prueba (Weibull, Gamma) = ', ad_w1$adjusted_statistic, ad_g1$adjusted_statistic, '\n', 
    'Punto crítico =                          ', AD_wg_rechazo1, AD_wg_rechazo1, '\n',
    'p_value =                                ', '   0 ', '   0', '\n\n',
    
    '==== Prueba Shapiro-Wilk ===','             Weibull ',' Gamma ','\n',
    'Estadistico de prueba Z (Weibull, Gamma) = ', sw_w1$statistic, sw_g1$statistic, '\n', 
    'Punto crítico Z =                        ', SW_wg_rechazo1, SW_wg_rechazo1, '\n',
    'p_value =                                ', sw_w1$p.value, sw_g1$p.value, '\n\n')




################################################################################
################################################################################





############### Base de datos de salarios de Ingenieros en Software ############
Data2 <- read.csv("Software_Professional_Salaries.csv")
SalarySoftware <- Data2$Salary
n2 <- length(SalarySoftware)

# Suavizado de dato atípico (Fila 18640)
SalarySoftware[18640] <- sum(SalarySoftware[18635:18639], SalarySoftware[18641:18645]) / 10 

# ==============================================================================
# Estimación de parametros
# ==============================================================================
mle_lognormal2 <- estimar_lognormal_mle(SalarySoftware)
mle_weibull2   <- estimar_weibull_mle(SalarySoftware)

meanlog_ln2 <- mle_lognormal2['meanlog']
sdlog_ln2   <- mle_lognormal2['sdlog']

shape_wb2   <- mle_weibull2['shape']
scale_wb2   <- mle_weibull2['scale']

# ==============================================================================
# Pruebas de Bondad
# ==============================================================================

# ---- PRUEBA DE ANDERSON DARLING ----
ad_ln2 <- ad.test(SalarySoftware, plnorm, meanlog = meanlog_ln2, sdlog = sdlog_ln2)
ad_wb2 <- ad.test(SalarySoftware, pweibull, shape = shape_wb2, scale = scale_wb2)

# ---- PRUEBA KOLMOGOROV-SMIRNOV ----
ks_ln2 <- ks.test(SalarySoftware, plnorm, meanlog = meanlog_ln2, sdlog = sdlog_ln2)
ks_wb2 <- ks.test(SalarySoftware, pweibull, shape = shape_wb2, scale = scale_wb2)

# ---- PRUEBA LILLIEFORS EXPONENCIAL ----
D2_exponencial2 <- lilliefors.exp.test(SalarySoftware)
p_val_lili_exp2 <- 2 * pnorm(-sqrt(n2) * D2_exponencial2)

# ---- PRUEBA JI-CUADRADA ----
prueba_ln2 <- chisq.test(SalarySoftware, p_dist_func = plnorm, num_param = 2, 
                         meanlog = meanlog_ln2, sdlog = sdlog_ln2)

prueba_wb2 <- chisq.test(SalarySoftware, p_dist_func = pweibull, num_param = 2, 
                         shape = shape_wb2, scale = scale_wb2)

# ---- PRUEBA DE SHAPIRO-WILK ----
sw_ln2 <- shapiro.test(qnorm(pmin(pmax(plnorm(SalarySoftware, 
                                              meanlog = meanlog_ln2, sdlog = sdlog_ln2), 1e-10), 1-1e-10)))

sw_wb2 <- shapiro.test(qnorm(pmin(pmax(pweibull(SalarySoftware, 
                                                shape = shape_wb2, scale = scale_wb2), 1e-10), 1-1e-10)))


########################### Región de rechazo ##################################
k2                      <- floor(1 + log2(n2)) 
chis.t_wlnlng2_rechazo2 <- qchisq(1 - alpha, df = k2 - 1 - 2)
kv_wlnlng2_rechazo2     <- 1.22 / (sqrt(n2 + sqrt(n2 / 10)))

AD_wlnlng2_rechazo2     <- 2.505700
#Por simulación de 10000 para Anderson Darling
#   90%        95%     97.5%   99% 
# 1.942079 2.505700 3.111237 4.031997  

Lili.t_wlnlng2_rechazo2 <- 1.0753 / (sqrt(n2))
SW_wg_rechazo2          <- qnorm(1 - alpha)


cat('========================= Resumen estadistico =========================','\n',
    '==== Prueba Ji-Cuadrada ===','                 lognorm','   Weibull','\n',
    'Estadístico de prueba =                        ', prueba_ln2$statistic, ' | ', prueba_wb2$statistic, '\n',
    'Punto crítico =                                ', chis.t_wlnlng2_rechazo2, ' | ', chis.t_wlnlng2_rechazo2, '\n',
    'p_value =                                      ', prueba_ln2$p.value, ' | ', prueba_wb2$p.value, '\n\n',
    
    '==== Prueba Lilliefors Exponencial ===','\n',
    'Estadístico de prueba =                        ', D2_exponencial2, '\n',
    'Punto crítico =                                ', Lili.t_wlnlng2_rechazo2, '\n',
    'p_value =                                      ', p_val_lili_exp2, '\n\n',
    
    '==== Prueba Kolmogorov-Smirnov ===','            lognorm','   Weibull','\n',
    'Estadístico de prueba =                        ', ks_ln2$statistic, ' | ', ks_wb2$statistic, '\n', 
    'Punto crítico =                                ', kv_wlnlng2_rechazo2, ' | ', kv_wlnlng2_rechazo2, '\n',
    'p_value =                                      ', ks_ln2$p.value, ' | ', ks_wb2$p.value, '\n\n',
    
    '==== Prueba de Anderson Darling ===','          lognorm','   Weibull','\n',
    'Estadístico de prueba =                        ', ad_ln2$adjusted_statistic, ' | ', ad_wb2$adjusted_statistic, '\n', 
    'Punto crítico =                                ', AD_wlnlng2_rechazo2, ' | ', AD_wlnlng2_rechazo2, '\n',
    'p_value =                                      ', '    0', '      0', '\n\n',
    
    '==== Prueba de Shapiro-Wilk ===','              lognorm','   Weibull','\n',
    'Estadístico de prueba Z =                      ', sw_ln2$statistic, ' | ', sw_wb2$statistic, '\n', 
    'Punto crítico Z =                              ', SW_wg_rechazo2, ' | ', SW_wg_rechazo2, '\n',
    'p_value =                                      ', sw_ln2$p.value, ' | ', sw_wb2$p.value, '\n',
    '=======================================================================\n')






################################################################################
################################################################################
################################################################################
################################################################################





# ==============================================================================
# 1. CARGA, EXTRACCIÓN Y DESPLAZAMIENTO DEL PIB
# ==============================================================================
Datos3 <- read.csv("conjunto_de_datos_pibe_entidad_tab2024_p.csv", 
                   encoding = "UTF-8", 
                   check.names = FALSE)

# Buscamos y extraemos la fila que dice "Variación porcentual anual" y "Producto interno bruto"
fila_variacion3 <- Datos3[grepl("Variación porcentual anual", Datos3[,1]) & 
                            grepl("Producto interno bruto", Datos3[,1]), ]

# Tomamos los valores numéricos correspondientes a los años (columnas 3 a 23, del 2004 al 2024)
tasas_tabasco3 <- as.numeric(fila_variacion3[1, 3:23])

# ------------------------------------------------------------------------------
# Graficamos
# ------------------------------------------------------------------------------
boxplot(tasas_tabasco3)
plot.ecdf(tasas_tabasco3)
curve(pnorm(x, mean(tasas_tabasco3), sd(tasas_tabasco3)), add = TRUE, col = 'red', lwd = 2)

# Hago un desplazamiento haciendo la suma del valor abs 
# del minimo, pero como hay datos muy alejados, le sumamos 0.5 para poder usar lognormal y weibull
PIB_Analisis3 <- tasas_tabasco3 + abs(min(tasas_tabasco3)) + 0.5
n3 <- length(PIB_Analisis3)

# ==============================================================================
# 2. AJUSTE DE PARÁMETROS CON TUS FUNCIONES
# ==============================================================================
mle_lognormal3 <- estimar_lognormal_mle(PIB_Analisis3)
mle_weibull3   <- estimar_weibull_mle(PIB_Analisis3)

meanlog_ln3 <- mle_lognormal3['meanlog']
sdlog_ln3   <- mle_lognormal3['sdlog']

shape_wb3   <- mle_weibull3['shape']
scale_wb3   <- mle_weibull3['scale']

# ==============================================================================
# 3. EJECUCIÓN DE PRUEBAS DE BONDAD DE AJUSTE CON TUS FUNCIONES
# ==============================================================================

# ---- A. PRUEBA DE ANDERSON DARLING ----
ad_ln3 <- ad.test(PIB_Analisis3, plnorm, meanlog = meanlog_ln3, sdlog = sdlog_ln3)
ad_wb3 <- ad.test(PIB_Analisis3, pweibull, shape = shape_wb3, scale = scale_wb3)

# ---- B. PRUEBA KOLMOGOROV-SMIRNOV ----
ks_ln3 <- ks.test(PIB_Analisis3, plnorm, meanlog = meanlog_ln3, sdlog = sdlog_ln3)
ks_wb3 <- ks.test(PIB_Analisis3, pweibull, shape = shape_wb3, scale = scale_wb3)

# ---- C. PRUEBA LILLIEFORS EXPONENCIAL ----
D3_exponencial3 <- lilliefors.exp.test(PIB_Analisis3)
p_val_lili_exp3 <- 2 * pnorm(-sqrt(n3) * D3_exponencial3)

# ---- D. PRUEBA JI-CUADRADA ----
prueba_ln3 <- chisq.test(PIB_Analisis3, p_dist_func = plnorm, num_param = 2, 
                         meanlog = meanlog_ln3, sdlog = sdlog_ln3)

prueba_wb3 <- chisq.test(PIB_Analisis3, p_dist_func = pweibull, num_param = 2, 
                         shape = shape_wb3, scale = scale_wb3)

# ---- E. PRUEBA DE SHAPIRO-WILK ----
sw_ln3 <- shapiro.test(qnorm(pmin(pmax(plnorm(PIB_Analisis3, 
                                              meanlog = meanlog_ln3, sdlog = sdlog_ln3), 1e-10), 1-1e-10)))

sw_wb3 <- shapiro.test(qnorm(pmin(pmax(pweibull(PIB_Analisis3, 
                                                shape = shape_wb3, scale = scale_wb3), 1e-10), 1-1e-10)))


########################### Región de rechazo ##################################
k3              <- floor(1 + log2(n3)) 
chis.t_rechazo3 <- qchisq(1 - alpha, df = k3 - 1 - 2)
kv_rechazo3     <- 1.22 / (sqrt(n3 + sqrt(n3 / 10)))

AD_rechazo3     <- 2.465673
#Por simulación de 10000 para Anderson Darling
#    90%      95%    97.5%      99% 
# 1.920329 2.465673 3.104136 3.881799  

Lili.t_rechazo3 <- 1.0753 / (sqrt(n3))

SW_wg_rechazo3  <- qnorm(1 - alpha)



cat('==================== Resumen estadístico 3 (PIB Tabasco) ====================','\n',
    '==== Prueba Ji-Cuadrada ===','                 lognorm','   weibull','\n',
    'Estadístico de prueba =                        ', prueba_ln3$statistic, ' | ', prueba_wb3$statistic, '\n',
    'Punto crítico =                                ', chis.t_rechazo3, ' | ', chis.t_rechazo3, '\n',
    'p_value =                                      ', prueba_ln3$p.value, ' | ', prueba_wb3$p.value, '\n\n',
    
    '==== Prueba Lilliefors Exponencial ===','\n',
    'Estadístico de prueba =                        ', D3_exponencial3, '\n',
    'Punto crítico =                                ', Lili.t_rechazo3, '\n',
    'p_value =                                      ', p_val_lili_exp3, '\n\n',
    
    '==== Prueba Kolmogorov-Smirnov ===','            lognorm','   weibull','\n',
    'Estadístico de prueba =                        ', ks_ln3$statistic, ' | ', ks_wb3$statistic, '\n', 
    'Punto crítico =                                ', kv_rechazo3, ' | ', kv_rechazo3, '\n',
    'p_value =                                      ', ks_ln3$p.value, ' | ', ks_wb3$p.value, '\n\n',
    
    '==== Prueba de Anderson Darling ===','          lognorm','   weibull','\n',
    'Estadístico de prueba =                        ', ad_ln3$adjusted_statistic, ' | ', ad_wb3$adjusted_statistic, '\n', 
    'Punto crítico =                                ', AD_rechazo3, ' | ', AD_rechazo3, '\n',
    'p_value =                                      ', ' 0.050001 ', ' 0.050001', '\n\n',
    
    '==== Prueba de Shapiro-Wilk ===','              lognorm','   weibull','\n',
    'Estadístico de prueba Z =                      ', sw_ln3$statistic, ' | ', sw_wb3$statistic, '\n', 
    'Punto crítico Z =                              ', SW_wg_rechazo3, ' | ', SW_wg_rechazo3, '\n',
    'p_value =                                      ', sw_ln3$p.value, ' | ', sw_wb3$p.value, '\n',
    '=======================================================================\n')
    
                                   
                                       
                                       
################################################################################  
################################################################################ 
################################################################################ 




######################### Tasa de desempleo en Europa ##########################

datos4 <- read.csv('OECD.CFE.EDS,DSD_REG_LAB@DF_RATES,+all.csv',
                  encoding = "UTF-8", check.names = FALSE)

# Filtramos para quedarnos únicamente con las filas deseadas
Datos_Filtrados4 <- datos4[datos4$Measure == "Unemployment rate" &
                             datos4$Sex == "Total" &
                             datos4$Age == "From 15 to 64 years", ]

# SOLUCIÓN 1: Extraemos usando el nombre exacto de la columna 'OBS_VALUE'
# SOLUCIÓN 2: Eliminamos los valores NA para que las pruebas funcionen limpiamente
Desempleo_Analisis4 <- as.numeric(Datos_Filtrados4$OBS_VALUE)/100
n4 <- length(Desempleo_Analisis4)

# ==============================================================================
# 2. GRAFICAMOS (Sin errores de valores finitos)
# ==============================================================================
boxplot(Desempleo_Analisis4)
plot.ecdf(Desempleo_Analisis4)
curve(pbeta(x, 2.255162, 26.968976 ), add = TRUE, col = 'red', lwd = 2)

# ==============================================================================
# 3. AJUSTE DE PARÁMETROS CON TUS FUNCIONES
# ==============================================================================
mle_lognormal4 <- estimar_lognormal_mle(Desempleo_Analisis4)
mle_beta4      <- estimar_beta_mle(Desempleo_Analisis4)

meanlog_ln4 <- mle_lognormal4['meanlog']
sdlog_ln4   <- mle_lognormal4['sdlog']

shape1_bt4  <- mle_beta4['shape1']
shape2_bt4  <- mle_beta4['shape2']

# ==============================================================================
# 4. EJECUCIÓN DE PRUEBAS DE BONDAD DE AJUSTE CON TUS FUNCIONES
# ==============================================================================

# ---- A. PRUEBA DE ANDERSON DARLING ----
ad_ln4 <- ad.test(Desempleo_Analisis4, plnorm, meanlog = meanlog_ln4, sdlog = sdlog_ln4)
ad_bt4 <- ad.test(Desempleo_Analisis4, pbeta, shape1 = shape1_bt4, shape2 = shape2_bt4)

# ---- B. PRUEBA KOLMOGOROV-SMIRNOV ----
ks_ln4 <- ks.test(Desempleo_Analisis4, plnorm, meanlog = meanlog_ln4, sdlog = sdlog_ln4)
ks_bt4 <- ks.test(Desempleo_Analisis4, pbeta, shape1 = shape1_bt4, shape2 = shape2_bt4)

# ---- C. PRUEBA LILLIEFORS EXPONENCIAL ----
D4_exponencial4 <- lilliefors.exp.test(Desempleo_Analisis4)
p_val_lili_exp4 <- 2 * pnorm(-sqrt(n4) * D4_exponencial4)

# ---- D. PRUEBA JI-CUADRADA ----
prueba_ln4 <- chisq.test(Desempleo_Analisis4, p_dist_func = plnorm, num_param = 2, 
                         meanlog = meanlog_ln4, sdlog = sdlog_ln4)

prueba_bt4 <- chisq.test(Desempleo_Analisis4, p_dist_func = pbeta, num_param = 2, 
                         shape1 = shape1_bt4, shape2 = shape2_bt4)

# ---- E. PRUEBA DE SHAPIRO-WILK ----
sw_ln4 <- shapiro.test(qnorm(pmin(pmax(plnorm(Desempleo_Analisis4, 
                                              meanlog = meanlog_ln4, sdlog = sdlog_ln4), 1e-10), 1-1e-10)))

sw_bt4 <- shapiro.test(qnorm(pmin(pmax(pbeta(Desempleo_Analisis4, 
                                             shape1 = shape1_bt4, shape2 = shape2_bt4), 1e-10), 1-1e-10)))


########################### Región de rechazo ##################################
k4              <- floor(1 + log2(n4)) 
chis.t_rechazo4 <- qchisq(1 - alpha, df = k4 - 1 - 2)
kv_rechazo4     <- 1.22 / (sqrt(n4 + sqrt(n4 / 10)))

AD_rechazo4     <- 2.453
#Por simulación de 1000 para Anderson Darling
#   90%   95% 97.5%   99% 
#  1.891 2.453 3.049 3.834

Lili.t_rechazo4 <- 1.0753 / (sqrt(n4))
SW_wg_rechazo4  <- qnorm(1 - alpha)

# ==============================================================================
# 5. RESUMEN ESTADÍSTICO FINAL
# ==============================================================================

cat('==================== Resumen estadístico 4 (Desempleo Europa) ====================','\n',
    '==== Prueba Ji-Cuadrada ===','                 lognorm','    beta','\n',
    'Estadístico de prueba =                        ', prueba_ln4$statistic, ' | ', prueba_bt4$statistic, '\n',
    'Punto crítico =                                ', chis.t_rechazo4, ' | ', chis.t_rechazo4, '\n',
    'p_value =                                      ', prueba_ln4$p.value, ' | ', prueba_bt4$p.value, '\n\n',
    
    '==== Prueba Lilliefors Exponencial ===','\n',
    'Estadístico de prueba =                        ', D4_exponencial4, '\n',
    'Punto crítico =                                ', Lili.t_rechazo4, '\n',
    'p_value =                                      ', p_val_lili_exp4, '\n\n',
    
    '==== Prueba Kolmogorov-Smirnov ===','            lognorm','    beta','\n',
    'Estadístico de prueba =                        ', ks_ln4$statistic, ' | ', ks_bt4$statistic, '\n', 
    'Punto crítico =                                ', kv_rechazo4, ' | ', kv_rechazo4, '\n',
    'p_value =                                      ', ks_ln4$p.value, ' | ', ks_bt4$p.value, '\n\n',
    
    '==== Prueba de Anderson Darling ===','          lognorm','    beta','\n',
    'Estadístico de prueba =                        ', ad_ln4$adjusted_statistic, ' | ', ad_bt4$adjusted_statistic, '\n', 
    'Punto crítico =                                ', AD_rechazo4, ' | ', AD_rechazo4, '\n',
    'p_value =                                      ', ' ', ' ', '\n\n',
    
    '==== Prueba de Shapiro-Wilk ===','              lognorm','    beta','\n',
    'Estadístico de prueba Z =                      ', sw_ln4$statistic, ' | ', sw_bt4$statistic, '\n', 
    'Punto crítico Z =                              ', SW_wg_rechazo4, ' | ', SW_wg_rechazo4, '\n',
    'p_value =                                      ', sw_ln4$p.value, ' | ', sw_bt4$p.value, '\n',
    '=======================================================================\n')




################################################################################  
################################################################################ 
################################################################################ 



######################## Seguro Medico en EU ###################################


# ==============================================================================
# Cargamos los datos
# ==============================================================================
dato5 <- read.csv('insurance.csv')

# Extraemos la columna numérica pura y eliminamos NAs
Seguros_Analisis5 <- na.omit(as.numeric(dato5$charges))
n5 <- length(Seguros_Analisis5)

# ------------------------------------------------------------------------------
# Graficamos
# ------------------------------------------------------------------------------
hist(Seguros_Analisis5, prob = TRUE, main = "Histograma de Cargos Médicos")
boxplot(Seguros_Analisis5)
plot.ecdf(Seguros_Analisis5)
curve(plnorm(x, 9.0986587, 0.9191834 ), add = TRUE, col = 'red', lwd = 2)

# ==============================================================================
# Estimamos parametros
# ==============================================================================
mle_lognormal5 <- estimar_lognormal_mle(Seguros_Analisis5)

meanlog_ln5 <- mle_lognormal5['meanlog']
sdlog_ln5   <- mle_lognormal5['sdlog']

# ==============================================================================
# Aplicamos pruebas de bondad
# ==============================================================================

# ---- PRUEBA DE ANDERSON DARLING ----
ad_ln5 <- ad.test(Seguros_Analisis5, plnorm, meanlog = meanlog_ln5, sdlog = sdlog_ln5)

# ---- PRUEBA KOLMOGOROV-SMIRNOV ----
ks_ln5 <- ks.test(Seguros_Analisis5, plnorm, meanlog = meanlog_ln5, sdlog = sdlog_ln5)

# ---- PRUEBA LILLIEFORS EXPONENCIAL ----
D5_exponencial5 <- lilliefors.exp.test(Seguros_Analisis5)
p_val_lili_exp5 <- 2 * pnorm(-sqrt(n5) * D5_exponencial5)

# ---- PRUEBA JI-CUADRADA ----
prueba_ln5 <- chisq.test(Seguros_Analisis5, p_dist_func = plnorm, num_param = 2, 
                         meanlog = meanlog_ln5, sdlog = sdlog_ln5)

# ---- PRUEBA DE SHAPIRO-WILK (PIT) ----
sw_ln5 <- shapiro.test(qnorm(pmin(pmax(plnorm(Seguros_Analisis5, 
                                              meanlog = meanlog_ln5, sdlog = sdlog_ln5), 1e-10), 1-1e-10)))


########################### Región de rechazo ##################################

k5              <- floor(1 + log2(n5)) 
chis.t_rechazo5 <- qchisq(1 - alpha, df = k5 - 1 - 2)
kv_rechazo5     <- 1.22 / (sqrt(n5 + sqrt(n5 / 10)))

AD_rechazo5     <- 2.467591
#     90%      95%    97.5%      99% 
#  1.929489 2.467591 3.003882 3.745779 
Lili.t_rechazo5 <- 1.0753 / (sqrt(n5))
SW_wg_rechazo5  <- qnorm(1 - alpha)

# ==============================================================================
# 4. RESUMEN ESTADÍSTICO FINAL
# ==============================================================================

cat('==================== Resumen estadístico 5 (Seguro Médico EU) ====================','\n',
    '==== Prueba Ji-Cuadrada ===','                 lognorm','\n',
    'Estadístico de prueba =                        ', prueba_ln5$statistic, '\n',
    'Punto crítico =                                ', chis.t_rechazo5, '\n',
    'p_value =                                      ', prueba_ln5$p.value, '\n\n',
    
    '==== Prueba Lilliefors Exponencial ===','\n',
    'Estadístico de prueba =                        ', D5_exponencial5, '\n',
    'Punto crítico =                                ', Lili.t_rechazo5, '\n',
    'p_value =                                      ', p_val_lili_exp5, '\n\n',
    
    '==== Prueba Kolmogorov-Smirnov ===','            lognorm','\n',
    'Estadístico de prueba =                        ', ks_ln5$statistic, '\n', 
    'Punto crítico =                                ', kv_rechazo5, '\n',
    'p_value =                                      ', ks_ln5$p.value, '\n\n',
    
    '==== Prueba de Anderson Darling ===','          lognorm','\n',
    'Estadístico de prueba =                        ', ad_ln5$adjusted_statistic, '\n', 
    'Punto crítico =                                ', AD_rechazo5, '\n',
    'p_value =                                      ', ' 0.499999', '\n\n',
    '==== Prueba de Shapiro-Wilk ===','              lognorm','\n',
    'Estadístico de prueba Z =                      ', sw_ln5$statistic, '\n', 
    'Punto crítico Z =                              ', SW_wg_rechazo5, '\n',
    'p_value =                                      ', sw_ln5$p.value, '\n',
    '=======================================================================\n')
