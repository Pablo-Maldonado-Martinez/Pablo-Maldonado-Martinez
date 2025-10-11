#Ecología Teórica
#Especies Invasoras: Taller de Modelación

#El presente script se modífico con fines académicos; para su ejemplificación del uso de modelos, análisis y gráficos del lenguaje de programación R, en el tema ecológico de las especies invasoras y sus efectos. 
#Árticulo original: Local coexistence of native and invasive ant species is associated with micro-spatial shifts in foraging activity
#Autores: Gippet, George & Bertelsmeier (2021)
#DOI: https://doi.org/10.1007/s10530-021-02678-2

#Equipo  a cargo:
#Farias Ruiz Sara Fernanda
#Garcia Morales Eva Mariani 
#González López Luz Elena
#Guzmán Alvarez Olin
#Maldonado Martinez Pablo Adahir

#------------LIBRERÍAS

# Gráficos
library(ggplot2) #Sistema elegante para crear gráficos
library(vioplot) #Gráficos de violín (boxplot + densidad)
library(gridExtra) #Combinar múltiples gráficos en uno
library(scales) #Formatear ejes y leyendas de gráficos
library(randomcoloR) #Genera paleta de colores aleatorias
library(viridis) #Paleta de colores específica

#Modelos Estadísticos
library(glmmTMB) #Modelos mixtos generalizados (con varianza)
library(performance) #Diagnosticar y comparar modelos
library(DHARMa) #Diagnosticar modelos mixtos (residuos)
library(buildmer) #Selección automática de modelos
library(ggeffects) #Visualizar efectos de modelos


#Análisis estadísticos
library(car) #Análisis de regresión y ANOVA
library(emmeans) #Medias marginales y comparaciones múltiples
library(multcomp) #Comparaciones múltiples post-hoc

#Análisis Multivariados
library(vegan) #Ecología numérica y ordenaciones
library(ade4) #Análisis de datos ecológicos
library(factoextra) #Visualizar análisis multivariados


#-------------PARTE 1

# PALETA DE COLORES DE LIVIE <3
livie_Cardinals <- c("firebrick", "purple", "slateblue", "violetred")
# North: "firebrick"
# East: "purple"
# South: "slateblue"
# West: "violetred"

livie_Invasion <- c("purple", "firebrick")
# Non-invaded: "purple"
# Invaded: "firebrick"


#-------------PARTE 2
# Cargar Base de datos

# Datos de los cebos
dataTI <- read.table("Gippet2021_BioInv_dataBaits.txt", h=T, sep="\t") 
head(dataTI) #Primeras Filas
dim(dataTI) #Dimensiones

# Datos de los edificios (a nivel de cebo >> muchos valores duplicados, eso es normal)
dataBuildings <- read.table("Gippet2021_BioInv_dataBuildings.txt", h=T, sep="\t")
head(dataBuildings) #Primeras Filas
dim(dataBuildings) #Dimensiones

#-------------PARTE 3
# Manipulación de datos y preparación para análisis

# Fusión de datos de cebos y edificios
dataTI2 <- cbind(dataTI, dataBuildings[,7:12][match(dataTI$ID, dataBuildings$ID),])
head(dataTI2)  #Primeras Filas
dim(dataTI2)   #Dimensiones: 3840  29(23 + 6 )

# Crea una nueva columna que indique si el cebo estaba cerca (posición A) o lejos (posición B) del edificio
dataTI2$positionAB <- rep(c("A","B"), 1920) #1920*2=3840, Número de Filas(obs)

for (i in 1:dim(dataTI2)[1]){
  if (grepl("A",dataTI2$position[i])){ 
    dataTI2$positionAB[i] <- "A"
  } else {dataTI2$positionAB[i] <- "B"}
}                                        

# Niveles de factores de reordenamiento
dataTI2$time <- as.factor(dataTI2$time)
dataTI2$time = factor(dataTI2$time, levels(dataTI2$time)[c(2,3,1)]) 
dataTI2$Bface <- as.factor(dataTI2$Bface)
dataTI2$Bface = factor(dataTI2$Bface, levels(dataTI2$Bface)[c(2,1,3,4)]) 
dataTI2$Date <- as.factor(dataTI2$Date)
dataTI2$building <- as.factor(dataTI2$building)
dataTI2$building = factor(dataTI2$building, levels(dataTI2$building)[c(1,9,10,11,12,13,14,15,16,2,3,4,5,6,7,8)]) 

# Diseño de muestreo (número de cebos por eventos de muestreo y día en que se realizó cada evento de muestreo, por edificio)
table(dataTI2$Date, dataTI2$time, dataTI2$building)

# Número total de hormigas obreras nativas reclutadas en cada cebo
dataTI2$Ab_allNatives <- rowSums(dataTI2[ ,10:23])
head(dataTI2) #Primeras filas
dim(dataTI2) #Dimensiones

# Cambiando la columna ID con una nueva columna IDbait (sin información de atún y miel)
dataTI2$IDbait <- as.factor(gsub('.{4}$', '', dataTI2$ID))
dataTI2 <- dataTI2[,2:32]
head(dataTI2) #Primeras filas
dim(dataTI2) #Dimensiones

# Fusionar líneas de atún y miel por cebo para tener un conjunto de datos a nivel de cebo
dataTI3_part1 <- dataTI2[seq(2,3840, by=2),c(31,27,28,29,1,2,3,4,5,7,23,24,25,26)]
head(dataTI3_part1) #Primeras filas
dim(dataTI3_part1) #Dimensiones

dataTI3_part2 <- aggregate(dataTI2[,c(8:22,30)], by=list(dataTI2$IDbait), sum) 
head(dataTI3_part2) #Primeras filas
dim(dataTI3_part2) #Dimensiones
colnames(dataTI3_part2)[1] <- "IDbait"

dataTI3 <- cbind(dataTI3_part1, dataTI3_part2[match(dataTI3_part1$IDbait, dataTI3_part2$IDbait),])

# Comprobando si hay desajustes
for (i in 1:dim(dataTI3)[1]){
  if (dataTI3[i,1] != dataTI3[i,15]){
    print("MISMATCH PROBLEM")} else{}
} #Comprueba la consistencia de los datos

dataTI3 <- dataTI3[,-15] #Esto sirve de ejemplo para eliminar la columna 15

dataTI3$presence_Natives <- 0
for (i in 1:1920){if(dataTI3$Ab_allNatives[i] > 0){dataTI3$presence_Natives[i] <- 1}} #Asignar 0 o 1 a los datos de "Ab_allNatives"
head(dataTI3)

# Niveles de factores de reordenamiento
dataTI3$positionAB <- as.factor(dataTI3$positionAB)
dataTI3$zone <- as.factor(dataTI3$zone)
dataTI3$building <- as.factor(dataTI3$building)
dataTI3$position <- as.factor(dataTI3$position)
head(dataTI3) #Primeras filas
dim(dataTI3)  #Dimensiones


#-------------PARTE 4

# Variaciones microclimáticas inducidas por las condiciones de sombreado
dataTI3_Temperature <- dataTI3
dataTI3_Temperature$IDface <- paste0(dataTI3_Temperature$building, dataTI3_Temperature$time, dataTI3_Temperature$Bface)
dataTI3_Temperature <-  dataTI3_Temperature[!duplicated(dataTI3_Temperature$IDface),]
dim(dataTI3_Temperature)
head(dataTI3_Temperature)

Tmean_test <- glmmTMB(Tground_mean   ~ 1 +
                        Bface + 
                        time + 
                        time:Bface +
                        (1|building) + (1|Date), 
                      data=dataTI3_Temperature, 
                      family=gaussian)
ef_Tmean_test <- ggemmeans(Tmean_test, c("time", "Bface"), type = "fixed")

# FIGURA 1.D
plot(ef_Tmean_test, colors = livie_Cardinals) + 
  theme_gray() +
  labs(
    title="Temperatura del Suelo respecto al Estado de Tiempo",
    x="Temperatura de Suelo(C°)", 
    y="Estado del Tiempo", 
    caption="Figura 1.D",
  )

# FIGURA 1.C
vioplot(Tground_mean ~ zone, data=dataTI3_Temperature, col=livie_Invasion, 
        main="Temperatura del Suelo respecto a la Zona",
        xlab ="Zona", 
        ylab ="Temperatura de Suelo(C°)", 
        sub ="Figura 1.C"
        )

#-------------PARTE 5
# Efecto de T. magnum en las comunidades de hormigas nativas

# ¿Tapinoma influye en la abundancia, riqueza y diversidad de especies de hormigas nativas?
dataTI3presence <- dataTI3
for (i in 1:1920) {
  for (j in 15:29) {
    if(dataTI3presence[i,j] > 0){
      dataTI3presence[i,j] <- 1}
  }
} 
head(dataTI3presence)
dim(dataTI3presence)

# Agregado a nivel de edificio
dataTI3presence_v2 <- aggregate(dataTI3presence[,15:29], 
                                by = list(dataTI3presence$zone, 
                                          dataTI3presence$building, 
                                          dataTI3presence$time,
                                          dataTI3presence$Date,
                                          dataTI3presence$Age_building), 
                                sum)
colnames(dataTI3presence_v2)[1:5] <- c("zone", "building", "time", "Date", "age_building")
head(dataTI3presence_v2)
dim(dataTI3presence_v2)
dataTI3presence_v2 <- dataTI3presence_v2[order(dataTI3presence_v2$zone, dataTI3presence_v2$building, dataTI3presence_v2$time),]

dataTI3presence_v2$nbBaits <- 40
dataTI3presence_v2$allNative <-  rowSums(dataTI3presence_v2[,7:20])
dataTI3presence_v2$allANTS <-  rowSums(dataTI3presence_v2[,6:20])

# Riqueza de especies por evento de muestreo
dataTI3presence_v2$richness <- NA
for (i in 1:dim(dataTI3presence_v2)[1]){
  dataTI3presence_v2$richness[i] = 14 - length(which(dataTI3presence_v2[i,7:20]==0))
}

# Diversidad de Shannon mediante evento de muestreo
dataTI3presence_v2$diversity <- diversity(dataTI3presence_v2[,7:20], index = "shannon")
hist(dataTI3presence_v2$diversity)
head(dataTI3presence_v2)
dim(dataTI3presence_v2)


#dataTI3presence_v2$building
#colorBuilding <- c(rep(distinctColorPalette(8), each=3),rep(distinctColorPalette(8), each=3))
#Formas
shapeBuilding <- c(rep(c(0:7), each=3), rep(c(0:7), each=3))

# FIGURA 2.A
p_richness <- ggplot(dataTI3presence_v2, aes(y=richness, x=zone, fill=building)) +
  geom_boxplot(width=0.5, alpha=0.9, fill=livie_Invasion, show.legend = FALSE) +
  geom_point(position=position_jitter(width = 0.25, height = 0, seed=3), shape=shapeBuilding, size=2.5, stroke=1.5, aes(color=zone), show.legend = FALSE) +
  geom_text(position=position_jitter(width = 0.25, height = 0,seed=3), aes(label=building), size=2.5) +
  scale_color_manual(values = livie_Invasion) +
  #ylim(c(0,10)) +
  scale_y_continuous( breaks=c(0,1,2,3,4,5,6,7,8,9), labels=c(0,1,2,3,4,5,6,7,8,9), limits=c(0.75,9.25)) +
  labs( 
    title = "Riqueza por Zona",
    x = "Zona",
    y = "Riqueza",
    caption = "Figura 2.A"
  ) +
  theme_classic()
p_richness


# FIGURA 2.C
p_diversity <- ggplot(dataTI3presence_v2, aes(y=diversity, x=zone, fill=building)) + 
  geom_boxplot(width=0.5, alpha=0.9, fill=livie_Invasion, show.legend = FALSE) + 
  geom_point(position=position_jitter(width = 0.25, height = 0, seed=3), shape=shapeBuilding, size=2.5, stroke=1.5, aes(color=zone), show.legend = FALSE) +
  geom_text(position=position_jitter(width = 0.25, height = 0,seed=3), aes(label=building), size=2.5) +
  scale_color_manual(values = livie_Invasion) +
  #ylim(c(0,10)) +
  scale_y_continuous( breaks=c(0,0.25,0.5,0.75,1,1.25,1.5,1.75), labels=c(0,0.25,0.5,0.75,1,1.25,1.5,1.75), limits=c(0,1.85)) +
  labs( 
  title = "Diversidad por Zona",
  x = "Zona",
  y = "Índice de Diversidad de Shannon",
  caption = "Figura 2.C"
  ) +
  theme_classic()
p_diversity

# FIGURA 2.B
p_abundance <- ggplot(dataTI3presence_v2, aes(y=allNative, x=zone, fill=building)) + 
  geom_boxplot(width=0.5, alpha=0.9, fill=livie_Invasion, show.legend = FALSE) + 
  geom_point(position=position_jitter(width = 0.25, height = 0, seed=3), shape=shapeBuilding, size=2.5, stroke=1.5, aes(color=zone), show.legend = FALSE) +
  geom_text(position=position_jitter(width = 0.25, height = 0,seed=3), aes(label=building), size=2.5) +
  scale_color_manual(values = livie_Invasion) +
  #ylim(c(0,42)) +
  scale_y_continuous( breaks=c(0,8,16,24,32,40), labels= round(c(0,8,16,24,32,40)/40*100)) +
  labs( 
    title = "Abundancia por Zona",
    x = "Zona",
    y = "Abundancia Relativa",
    caption = "Figura 2.B"
  ) +
  theme_classic()
p_abundance

# FIGURA COMPLETA (A,B Y C)
figure2 <- function(){
  grid.arrange(p_richness, p_abundance, p_diversity, nrow = 1)
}
figure2()


# Trazando el reclutamiento de especies nativas en cebos

# Nuevo conjunto de datos para agregar información
head(dataTI3)
dataTI3_2 <- dataTI3
for (i in 15:29) {
  for (j in 1: dim(dataTI3_2)[1]){
    if (dataTI3_2[j,i] > 0) {dataTI3_2[j,i] <- 1}
  }
}
head(dataTI3_2)
dataTI3_2$Ab_allNatives[dataTI3_2$Ab_allNatives>0] <- 1
dataTI3_2$Ab_baits_allspecies <- as.numeric(as.character(dataTI3_2$Ab_allNatives)) + dataTI3_2$tap_mag
dataTI3_2$Ab_baits_allspecies <- as.numeric(as.character(dataTI3_2$Ab_baits_allspecies))
range(dataTI3_2$Ab_baits_allspecies)
dataTI3_2$Ab_baits_allspecies[dataTI3_2$Ab_baits_allspecies>0] <- 1

# Crea un conjunto de datos para trazar el diagrama de barras.
ALLSpecies_occBaits <- aggregate(as.numeric(as.character(dataTI3_2$Ab_baits_allspecies)), by=list(dataTI3_2$building), mean)
colnames(ALLSpecies_occBaits) <- c("building", "meanOccBaits")
ALLSpecies_occBaits$zone <- c(rep("Invaded", 8), rep("Non-Invaded", 8))
w = wilcox.test(ALLSpecies_occBaits$meanOccBaits ~ ALLSpecies_occBaits$zone)
w2 <- c(round(as.numeric(w$statistic, digits=0)), round(w$p.value, digits=4))
mean_species_occBaits <- aggregate(ALLSpecies_occBaits$meanOccBaits, by=list(ALLSpecies_occBaits$zone), mean)
colnames(mean_species_occBaits) <- c("zone", "meanOccBaits")
mean_species_occBaits$median <- aggregate(ALLSpecies_occBaits$meanOccBaits, by=list(ALLSpecies_occBaits$zone), median)[,2]
mean_species_occBaits$se <- aggregate(ALLSpecies_occBaits$meanOccBaits, by=list(ALLSpecies_occBaits$zone), sd)[,2]/sqrt(8)
mean_species_occBaits$species <- "All species"
mean_species_occBaits$statsWilcox <- w2
Stats_occBaits <- mean_species_occBaits

NativeSpecies_occBaits <- aggregate(as.numeric(as.character(dataTI3_2$Ab_allNatives)), by=list(dataTI3_2$building), mean)
colnames(NativeSpecies_occBaits) <- c("building", "meanOccBaits")
NativeSpecies_occBaits$zone <- c(rep("Invaded", 8), rep("Non-Invaded", 8))
w = wilcox.test(NativeSpecies_occBaits$meanOccBaits ~ NativeSpecies_occBaits$zone)
w2 <- c(round(as.numeric(w$statistic, digits=0)), round(w$p.value, digits=4))
mean_species_occBaits <- aggregate(NativeSpecies_occBaits$meanOccBaits, by=list(NativeSpecies_occBaits$zone), mean)
colnames(mean_species_occBaits) <- c("zone", "meanOccBaits")
mean_species_occBaits$median <- aggregate(NativeSpecies_occBaits$meanOccBaits, by=list(NativeSpecies_occBaits$zone), median)[,2]
mean_species_occBaits$se <- aggregate(NativeSpecies_occBaits$meanOccBaits, by=list(NativeSpecies_occBaits$zone), sd)[,2]/sqrt(8)
mean_species_occBaits$species <- "Native only"
mean_species_occBaits$statsWilcox <- w2
Stats_occBaits <- rbind(Stats_occBaits, mean_species_occBaits)


for (i in 15:29){
  species_occBaits <- aggregate(dataTI3_2[,i], by=list(dataTI3_2$building), mean)
  colnames(species_occBaits) <- c("building", "meanOccBaits")
  species_occBaits$zone <- c(rep("Invaded", 8), rep("Non-Invaded", 8))
  w = wilcox.test(species_occBaits$meanOccBaits ~ species_occBaits$zone)
  w2 <- c(round(as.numeric(w$statistic, digits=0)), round(w$p.value, digits=4))
  
  mean_species_occBaits <- aggregate(species_occBaits$meanOccBaits, by=list(species_occBaits$zone), mean)
  colnames(mean_species_occBaits) <- c("zone", "meanOccBaits")
  mean_species_occBaits$median <- aggregate(species_occBaits$meanOccBaits, by=list(species_occBaits$zone), median)[,2]
  mean_species_occBaits$se <- aggregate(species_occBaits$meanOccBaits, by=list(species_occBaits$zone), sd)[,2]/sqrt(8)
  mean_species_occBaits$species <- colnames(dataTI3_2)[i]
  mean_species_occBaits$statsWilcox <- w2
  
  Stats_occBaits <- rbind(Stats_occBaits, mean_species_occBaits)
}

order_NI <- Stats_occBaits[Stats_occBaits$zone=="Non-Invaded",]
order_species_NI <- as.character(order_NI$species[order(order_NI$meanOccBaits, decreasing=T) ])

# Ploteo
Stats_occBaits$zone <- ordered(Stats_occBaits$zone, levels = c("Non-Invaded","Invaded"))
Stats_occBaits$species <- ordered(Stats_occBaits$specie, levels = order_species_NI)

# FIGURA 3
figure3 <- ggplot(data=Stats_occBaits[-c(1:4),], aes(x=species, y=meanOccBaits, fill=zone)) +
  geom_bar(stat = "identity", position="dodge",lwd=0.65, col="black", width=0.75) +
  geom_errorbar( aes(x=species, ymin=(meanOccBaits-se), ymax=(meanOccBaits+se)), 
                 position=position_dodge(0.75),
                 width=0.25, colour="black", alpha=0.9, linewidth =0.5) +
  geom_point(aes(x=species, y=median, fill=zone),
             position=position_dodge(0.75),
             color= "magenta", size=1.25) +
  labs(title = "Proporción de Cebos Ocupados por Especie",
       y="Proportion de Cebos Ocupados",
       x="Especies",
       caption = "Figura 3") +
  scale_fill_manual(values=livie_Invasion) + 
  scale_y_continuous(trans = 'pseudo_log') +
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"),
        axis.text.x = element_text(size=12, angle=45, hjust=1, vjust=1.15),
        axis.text.y = element_text(size=14))
figure3


#-------------PARTE 6
# Efecto de T. magnum y las condiciones de sombreado en la composición de la comunidad de hormigas nativas

# Composición de la comunidad

# Diferencias laterales del edificio
dataTI3recruit_vBface <- aggregate(dataTI3[,15:29],
                                    by = list(dataTI3presence$zone,
                                              dataTI3presence$building,
                                              dataTI3presence$time,
                                              dataTI3presence$Bface,
                                              dataTI3presence$Date),
                                    median)
colnames(dataTI3recruit_vBface)[1:5] <- c("zone", "building", "time", "Bface", "Date")
head(dataTI3recruit_vBface)
dim(dataTI3recruit_vBface)
colnames(dataTI3recruit_vBface)[6:20] <- paste0(colnames(dataTI3recruit_vBface)[6:20], "_recruit")

dataTI3presence_vBface <- aggregate(dataTI3presence[,15:29], 
                                by = list(dataTI3presence$zone, 
                                          dataTI3presence$building, 
                                          dataTI3presence$time,
                                          dataTI3presence$Bface,
                                          dataTI3presence$Date), 
                                sum)
colnames(dataTI3presence_vBface)[1:5] <- c("zone", "building", "time", "Bface", "Date")
head(dataTI3presence_vBface)
dim(dataTI3presence_vBface)  # 192 20
dataTI3presence_vBface <- dataTI3presence_vBface[order(dataTI3presence_vBface$zone, dataTI3presence_vBface$building, dataTI3presence_vBface$time, dataTI3presence_vBface$Bface),] # 

dataTI3presence_vBface$nbBaits <- 10
dataTI3presence_vBface$allNative <-  rowSums(dataTI3presence_vBface[,7:20])
dataTI3presence_vBface$allANTS <-  rowSums(dataTI3presence_vBface[,6:20])


# Riqueza de especies por evento de muestreo
dataTI3presence_vBface$richness <- NA
for (i in 1:dim(dataTI3presence_vBface)[1]){
  dataTI3presence_vBface$richness[i] = 14 - length(which(dataTI3presence_vBface[i,7:20]==0))
}

rowSums(dataTI3presence_vBface[,6:20])
test1 <- dataTI3presence_vBface[,7:20]+1

nmds1 <- metaMDS(test1, k=2, trymax=20) 

df_discri <- as.data.frame(nmds1$points)
colnames(df_discri) <- c("Axis1", "Axis2")
df_discri$Bface <- dataTI3presence_vBface$Bface 
df_discri$time <- dataTI3presence_vBface$time 
df_discri$zone <- dataTI3presence_vBface$zone
df_discri$building <- dataTI3presence_vBface$building

discri_sp  <- as.data.frame(nmds1$species)
colnames(discri_sp) <- c("Axis1", "Axis2")

species_names <- c("T. erraticum", 
                   "L. niger", 
                   "L. emarginatus", 
                   "M. specioides", 
                   "M. sabuleti", 
                   "M. schencki", 
                   "My. graminicola", 
                   "Tetramorium sp.", 
                   "S. fugax", 
                   "F. cunicularia", 
                   "F. rufibarbis", 
                   "F. fusca", 
                   "P. pygmeae", 
                   "Temnothorax sp.")

stat_Bface <- adonis2(test1 ~ dataTI3presence_vBface$Bface, permutations = 999, distance = "bray") 

# FIGURA 4.C
p_discri_Bface <-ggplot(df_discri, aes(x=Axis1, y=Axis2, color=Bface)) + 
                              geom_point(size=3) + 
                              stat_ellipse(aes(x=Axis1, y=Axis2, color=Bface), type = "norm", level=0.95, lwd=1.25) +
                              labs(title = "Comunidades por Lado",
                                   color = "Lado",
                                   caption = "Figura 4.C")+
                              theme_light() + 
                              scale_color_manual(values=livie_Cardinals) +
                              xlim(-0.35, 0.32) +
                              ylim(-0.32, 0.32)
p_discri_Bface

stat_time <- adonis2(test1 ~ dataTI3presence_vBface$time, permutations = 999, distance = "bray") 

# FIGURA 4.D
p_discri_time <-ggplot(df_discri, aes(x=Axis1, y=Axis2, color=time)) + 
                              geom_point(size=3) + 
                              stat_ellipse(aes(x=Axis1, y=Axis2, color=time),type = "norm", level=0.95, lwd=1.25) +
  labs(
    title = "Comunidades respecto al Estado de Tiempo",
    color = "Estado de Tiempo",
    caption = "Figura 4.D"
  )+
                              theme_light() +  
                              scale_color_manual(values=c("magenta4", "darkorchid", "deeppink4")) +
# Morning: "magenta4"
# Noon: "darkorchid"
# Afternoon: "deeppink4"
                              xlim(-0.35, 0.32) +
                              ylim(-0.32, 0.32)
p_discri_time

stat_zone <- adonis2(test1 ~ dataTI3presence_vBface$zone, permutations = 999, distance = "bray") 

# FIGURA 4.B
p_discri_zone <-ggplot(df_discri, aes(x=Axis1, y=Axis2, color=zone)) + 
                            geom_point(size=3) + 
                            stat_ellipse(aes(x=Axis1, y=Axis2, color=zone),type = "norm", level=0.95, lwd=1.25) +
  labs(
    title = "Comunidades por Zona",
    color = "Zona",
    caption = "Figura 4.B"
  )+
                            theme_light() + 
                            scale_color_manual(values=livie_Invasion) +
                            xlim(-0.35, 0.32) +
                            ylim(-0.32, 0.32)
p_discri_zone

# FIGURA 4.A
mult1 <- 1.2
p_discri_species <- ggplot(df_discri, aes(x=Axis1, y=Axis2, color=Bface)) + 
  labs(
    title = "Composición de la Comunidad de Hormigas Nativas",
    caption = "Figura 4.A",
    x = "Axis 1",
    y = "Axis 2"
  ) +
                              theme_light() + 
                              geom_segment(data = discri_sp, 
                                           aes(x = 0, y = 0, xend = (Axis1*mult1), yend = (Axis2*mult1)), 
                                           arrow = arrow(length = unit(1/2, "picas")),
                                           color = "black",
                                           lwd=1) +
                              annotate("text", x = (discri_sp$Axis1*(mult1+0.05)), y = (discri_sp$Axis2*(mult1+0.05)),
                                       label = species_names) +
                              xlim(-0.35, 0.32) +
                              ylim(-0.32, 0.32)
p_discri_species

# FIGURA 4 UNIDA (A, B, C Y D)
marge = 1
grid.arrange(p_discri_species, p_discri_zone, p_discri_Bface, p_discri_time, nrow = 2)

#-------------PARTE 7
# Efecto de T. magnum y las condiciones de sombreado en la actividad de alimentación de las hormigas nativas

# Lasius niger
binomial_lasnig <- glmmTMB(cbind(las_nig, nbBaits) ~ 1 +
                             zone + 
                             Bface +
                             time +
                             zone:Bface +
                             #zone:time +
                             #time:Bface +
                             (1|building) +
                             (1|Date), # remove this random effect to get non-NA result at the r2 computation below
                           data=dataTI3presence_vBface, 
                           family=binomial)

ef_binomial_lasnig_side <- ggemmeans(binomial_lasnig, c("Bface", "zone"), type = "fixed")

# FIGURA 5
plot(ef_binomial_lasnig_side,dodge=0.25) +
  labs (
    title = "Efecto de la Invasión de T. magnum",
    x = "Lado",
    y = "Proporción de Cebos Ocupados", 
    caption = "Figura 5",
    color = "Zona"
  )


