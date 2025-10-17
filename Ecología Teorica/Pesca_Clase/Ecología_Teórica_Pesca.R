#Ecología Teoríca - Caza y Pesca
#Maldonado Martínez Pablo Adahir
#Aplicación de un Modelo Logístico para calcular el Rendimiento Máximo Pesquero
#05/09/2025
#424062646@iztacala.unam.mx
#Notas: 

####################################
library(ggplot2)#Para graficar
library(tidyr)#Para ordenar
library(dplyr)#Para manejo de bases de datos
####################################

#-----------Parámetros
# Parámetros del modelo logístico
r <- 0.3     # tasa intrínseca de crecimiento
K <- 1000    # capacidad de carga
N0 <- 100    # población inicial
time <- seq(0, 50, by = 1)#define 51 unidades de tiempo (de 0 a 50) 
print(time)

#-----------Función logistic_growth()
# 1. Definimos una función llamada logistic_growth que recibe tres parámetros:
# N = tamaño poblacional en el tiempo actual
# r = tasa intrínseca de crecimiento
# K = capacidad de carga del ambiente

# 2. Se calcula la tasa de cambio de la población (dN) usando la ecuación logística:
# dN/dt = r * N * (1 - N/K)
# Este valor representa cuánto cambiará la población en la siguiente unidad de tiempo

# 3. Se devuelve el valor calculado de dN

# Función logística de crecimiento poblacional
logistic_growth <- function(N, r, K) {
dN <- r * N * (1 - N / K)
return(dN)
}
  
#------------Simulación sin extracción
# Vector para almacenar la población
population <- numeric(length(time))
population[1] <- N0 #Cuando haces simulaciones dinámicas (como el modelo logístico), necesitas partir de un valor inicial conocido. Ese valor es el tamaño poblacional al inicio del tiempo simulado, en este caso, N0.
  
# Simulación de crecimiento sin pesca
# Bucle for que simula el crecimiento poblacional a lo largo del tiempo, aplicando el modelo de crecimiento logístico definido antes.
  
#1. Inicia un bucle que va desde t = 2 hasta el final del vector time (en este caso, t = 51, ya que time va de 0 a 50).
#for (t in 2:length(time)) {
  
#2.  population[t] <- population[t - 1] + logistic_growth(population[t - 1], r, K)
#calcula el tamaño poblacional en el tiempo t, utilizando el valor anterior (t - 1) y aplicando la fórmula del modelo logístico.
  
#3.evita que la población supere la capacidad de carga K
# if (population[t] > K) population[t] <- K
  
  
for (t in 2:length(time)) {
    population[t] <- population[t - 1] + logistic_growth(population[t - 1], r, K)
    if (population[t] > K) population[t] <- K
  }
print(population)
  
#------------Cálculo del Rendimiento Máximo Sostenible (RMS)
# Fórmula RMS = (r * K) / 4
RMS <- (r * K) / 4
cat("Rendimiento Máximo Sostenible (RMS):", RMS)

#----------- Simulación con extracción constante
# Función de simulación con pesca (extracción constante)
simular_con_pesca <- function(N0, r, K, F, time) {
  N <- numeric(length(time))
  N[1] <- N0
  for (t in 2:length(time)) {
    crecimiento <- logistic_growth(N[t - 1], r, K)
    N[t] <- N[t - 1] + crecimiento - F
    if (N[t] < 0) N[t] <- 0
    if (N[t] > K) N[t] <- K
  }
  return(N)
}

# Simulamos con diferentes tasas de extracción
F_bajo <- 30
F_optimo <- RMS
F_alto <- 100

pop_bajo <- simular_con_pesca(N0, r, K, F_bajo, time)
pop_optimo <- simular_con_pesca(N0, r, K, F_optimo, time)
pop_alto <- simular_con_pesca(N0, r, K, F_alto, time)
print(pop_bajo)

#------------Visualización con ggplot2
# Construimos un data frame para graficar
df <- data.frame(
  tiempo = time,
  Sin_extraccion = population,
  Extraccion_baja = pop_bajo,
  Extraccion_RMS = pop_optimo,
  Extraccion_alta = pop_alto
)

# Transponemos la base de datos
df_long <- pivot_longer(df, cols = -tiempo, names_to = "Escenario", values_to = "Poblacion")

# Graficamos
ggplot(df_long, aes(x = tiempo, y = Poblacion, color = Escenario)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = K / 2, linetype = "dotted", color = "black") +
  annotate("text", x = 40, y = K / 2 + 30, label = "K/2 (máximo crecimiento)", size = 4) +
  labs(title = "Simulación del crecimiento poblacional con distintas tasas de pesca",
       subtitle = "Modelo logístico y extracción constante",
       x = "Tiempo",
       y = "Tamaño poblacional",
       color = "Escenario") +
  scale_color_manual(values = c(
    "Sin_extraccion" = "purple",
    "Extraccion_RMS" = "deepskyblue3",
    "Extraccion_baja" = "forestgreen",
    "Extraccion_alta" = "tomato"
  )) +
  theme_minimal()


#Librerías para R-Marckdown
install.packages(c("knitr", "rmarkdown", "xfun", "htmltools", "yaml", "tinytex"), dependencies = TRUE)
install.packages("devtools")
devtools::install_github("rstudio/rmarkdown")

install.packages("remotes")
remotes::install_cran("rmarkdown")

library(rmarkdown)

# Paso 1: Asegura que xfun se instala primero
install.packages("xfun", dependencies = TRUE)

# Paso 2: Instala knitr
install.packages("knitr", dependencies = TRUE)

# Paso 3: Instala tinytex
install.packages("tinytex", dependencies = TRUE)

# Paso 4: Instala rmarkdown
install.packages("rmarkdown", dependencies = TRUE)
# Paso 5: 
install.packages("remotes")
remotes::install_github("yihui/tinytex")
tinytex::install_tinytex()

# Actualiza xfun a la versión más reciente
install.packages("xfun")

# Luego instala rmarkdown y sus dependencias
install.packages("rmarkdown")

# Especificar que quieres las versiones binarias (no source)
install.packages("rmarkdown", type = "binary")
