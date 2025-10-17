# Práctica Epidemias
# Ecología teórica I
# 2/10/2025
# Artículo:  Long-term benefits of nonpharmaceutical interventions for
# endemicinfections are shaped by respiratory pathogen dynamics
# Autores: Rachel E. Baker, Chadi M. Saad-Roy, Sang Woo Park, Jeremy Farrar, C. Jessica E. Metcalf , and Bryan T. Grenfell 
####Código
require("deSolve")
require("scales")
require('zoo')
library(tsiR)
library("doBy")

### Parte 1: Simulación con parámetros fijos (Gripa; R0 bajo)
sirsmod_cs = function(t, y, parameters) {
  S = y[1]
  I = y[2]
  R = y[3]
  mu = parameters[["mu"]]
  gamma = parameters[["gamma"]]
  omega = parameters[["omega"]]### que tanto tarda un individuo en perder su inmunidad
  N = parameters[["N"]]
  controlstart = parameters[["controlstart"]]## que día inicia el uso de NPI
  controlend = parameters[["controlend"]]### que dái acaba el uso de NPI
  controlstart2 = parameters[["controlstart2"]]### que día inicia el uso de NPI longterm
  controlend2 = parameters[["controlend2"]]### que día acaba
  betachange = parameters[["betachange"]]### que tanto cambia b si se usan NPI
  
  beta=R0.list[t]*gamma ##Calcula el valor de b para cada tiempo
  if(t >= controlstart & t < controlend){
    beta = betachange*beta ## si se están utilizando NPI disminuye b 
  }
  if(t >= controlstart2 & t < controlend2){
    beta = betachange*beta ## si se están utilizando NPI disminuye b 
  }###Modelo SIRS con población constante
  dS = omega*R + mu * (N - S) - beta * S * I/N ##Calcular el cambio de S para cada tiempo
  dI = beta * S * I/N - (mu + gamma) * I ### Calcular el cambio de I para cada tiempo
  dR = gamma*I - omega*R - mu * R ### Calcular el cambio de R para cada tiempo
  res = c(dS, dI, dR) ###Agrupa los cambios en un vector
  list(res)
}
###Calcular y simular el valor de beta para cada semana dependiendo de la humedad ambiental

week <- seq(1,52,1)### has una secuencia de semanas del 1 al 52
qsim <-(0.006/2)*sin(2*pi*week/52 + 10.5) + 0.006 # Obtenemos la humedad ambiental de las semanas con una función sinoidad
qsim[qsim < 0] <- 0##Si en alguna semana el valor de la humedad es menoro a 0 corrigelo a 0 
plot(qsim) ### vamos a graficarlo
R0min = 1.2 ### la virulencia mínina
R0max = 3 ### la virulencia máxima
R0 = exp(-180*qsim + log(R0max - R0min)) + R0min ### calcula la virulencia para cada humedad dentro de un año
mean(R0) #obten la media

times = seq(1, 52 * 100,by = 1)### genera una secuencia con todas las semanas que queremos simular la epidemia
R0.list <- rep(R0,length=length(times))###repite los valores de la virulencia todos los años que vamos a simular

####Simular la epidemia
out = as.data.frame(ode(y = c(S = 0.19, I = 0.01, R = 0.8), times = times, 
                        func = sirsmod_cs, parms = list(mu = 1/(50 * 52), N = 1, R0.list = R0.list, gamma = 1,
                                                        omega = 1/40,
                                                        controlstart = (52*43 + 11), 
                                                        controlend =(52*43 + 11 + 52), betachange = 0.8, 
                                                        controlstart2 = (52*43 + 11 + 104) ,  controlend2 = (52*400))))


###Parámetros necesarios para la gráfica
subsetI <- out$I[times > (39*52)]
subsetS <- out$S[times > (39*52)]

timeuse <- seq(2016,2100,1/52)[1:length(subsetI)]
controlWeekStart = 11

### Lo graficamos
{
par(mar=c(3,3,1,3))
cols <- c(dichromat_pal("DarkRedtoBlue.12")(12)[1:4],dichromat_pal("DarkRedtoBlue.12")(12)[9:12])
cols2 <- c(dichromat_pal("BluetoOrange.10")(10))

timeuse <- seq(2016,2100,1/52)[1:length(subsetI)]
plot(timeuse,subsetI,type="n",col="#F64740",lwd =2,xlab="",ylab=
       "",bty = "n",xaxs="i",xlim=c(2016,2035),  yaxt ="n", xaxt= "n", ylim=c(0,0.3))
polygon(c(timeuse[52*4+controlWeekStart],timeuse[52*4+controlWeekStart],timeuse[52*4+controlWeekStart + 52],
          timeuse[52*4+controlWeekStart + 52]), c(-1,1,1,-1), border = NA,col="grey88")
polygon(c(timeuse[52*4+controlWeekStart + 104],timeuse[52*4+controlWeekStart + 104],max(timeuse),
          max(timeuse)), c(-1,1,1,-1), border = NA,col="grey88")
axis(side = 2, col.axis= cols[7], at = c(0.000,0.04, 0.08), label = c("0.00","0.04","0.08"))
lines(timeuse,subsetI,col=cols[7],lwd =2,xlab="",ylab=
        "")
subsetI2 <- out$I[times > (38.5*52)]
rm <- rollmean(subsetI2, k = 52 )[1:length(timeuse)]
lines(timeuse, rm, col="grey30", lwd = 2, lty = 2)
mtext(side= 2, "I/N", line = 2,col = cols[7], at = 0.001)
mtext(side= 2, "mean(I/N)", line = 2,col="grey30", at =0.06)
axis(side = 1, at = seq(2015,2035,5), labels =  seq(2015,2035,5))
abline(v = seq(2016,2100,1),lty=2,col="gray")

text(2030, 0.3,expression('"Influenza like": Low R'[0]*', SIRS'))


par(new = TRUE)
plot(timeuse, subsetS, col=cols[2],type="l",xlab="",ylab="", 
     axes = F,xaxs="i", ylim=c(0,1), lwd = 2, lty = 1,xlim=c(2016,2035))
axis(side = 4, col.axis=cols[2], at = c(0.4, 0.6), label = c(0.4, 0.6))
title(xlab="Year", line = 2)
mtext(side = 4, line = 2, 'S/N',col=cols[2])

par(new = TRUE)
subsetS2 <- out$S[times > (38.5*52)]
beta2 = c(R0[26:52] , rep(R0, length = length(subsetS2)))[1:length(subsetS2)]
beta2[245:length(beta2)] <- beta2[245:length(beta2)]*0.8

reff = beta2*subsetS2

rmreff <- rollmean(reff, k = 52 )[1:length(timeuse)]
plot(timeuse, rmreff, ylab = " ", xlab = "", yaxt = "n", xaxt = "n", ylim=c(-2.5,1.5), type="l", lwd =2, col =cols2[9],xlim=c(2016,2035),bty = "n",
     xaxs="i", lty = 2)
axis(side = 2, col.axis=cols2[9], at = c(0.5,1,1.5), label = c(0.5,1,1.5))

mtext(side = 2, line = 2, expression('mean(R'[e]*')'),col=cols2[9], at= 1, lty  = 1)
}

####Parte 2: Estimación de los parámetros através de series de tiempo
###Parte 2.1: Estiimación de los parámetros

load("TX.RData")
data
tsir <- runtsir(data=data, IP = 1, xreg ="cumcases",regtype= "spline",userYhat = NULL,
                alpha=0.97, family='poisson',link='log',method='negbin') ###

# Grafica el número de casos observados y estimados los el modelo TSIR
{
par(mar=c(3,3,1,1))
plot(tsir$res$time, tsir$res$cases, pch = 16, col="grey64", xlab="", ylab = "",ylim=c(0,1200))
polygon(c(tsir$res$time, rev(tsir$res$time))  ,c(tsir$res$mean + 1.96*tsir$res$sd, rev(tsir$res$mean - 1.96*tsir$res$sd)), border = NA, col=rgb(0,0,1,0.2))
title(xlab = "Year", line = 2)
title(ylab = "Cases", line = 2)
lines(tsir$res$time, tsir$res$mean, col="navy", lwd = 2)
}
#Gráfica de la beta estimada con los intervalos de confianza par cada semana
{
  par(mar=c(3,3,1,1))
  plot(tsir$beta*mean(data$pop), type="n", xlab="", ylab = "", bty="n")
  abline(h = mean(tsir$beta*mean(data$pop)),col="grey1", lty = 3)
  polygon(c(seq(1,52,1), rev(seq(1,52,1)))  ,c(tsir$contact$betahigh*mean(data$pop), rev(tsir$contact$betalow)*mean(data$pop)), border = NA, col="grey64")
  title(xlab="Week", line = 2)
  title(ylab="Beta (t)", line = 2)
  lines(tsir$beta*mean(data$pop), col="black", lwd = 2)
}

#Cálculo de R0
mean(tsir$beta*mean(data$pop))

####### Parte 2.2 Simulación 
source('predtsirMultiControlVax.R', encoding = 'UTF-8')###cargamos la función predtsirMultControlVax
# now we are going to predict forward, we need a dataset with the seasonal transmission rate and average births/pop 
# we can also use as seasonal birth rate here 
data2 <- data.frame(seas_beta_order <- seq(1,52,1), sea_beta = tsir$beta, births = 368190/52, pop = 29145505)

controlWeekStart = 11 # week in year that NPIs go into place

times <- seq(1,100,1/52) 
controlStart = 43*52 + controlWeekStart # week in simulation controls go in place #NB for endemic infections important to remove burn in period/transient dynamics
controlEnd =  controlStart + 52 # first controls in place for a year i.e. COVID-19 controls
p2controlStart = controlEnd + 52 # add controls back in after a year
p2controlEnd = length(times)
#p2controlEnd = p2controlStart + 8

pred <- predtsirMultiControlVax(times = times, births = rep(data2$births, length = length(times)), beta = data2$sea_beta, alpha = 0.97, 
                                S0 =floor(0.8*data2$pop[1]), I0 = floor(0.2*data2$pop[1]), nsim = 10, stochastic = F, 
                                controlStart = controlStart, controlEnd =controlEnd, betachange = 0.8,
                                p2controlStart =  p2controlStart, p2controlEnd = p2controlEnd , p2betachange = 0.8)
# I wrote this function, but it just does the forward simulation and allows you to put in a control period

subsetI <- pred$I$mean[times > 40]
subsetS <- pred$S$mean[times > 40]
timeuse <- seq(2016,2100,1/52)[1:length(subsetI)]

{
  par(mar=c(3,3,1,3))
  cols <- c(dichromat_pal("DarkRedtoBlue.12")(12)[1:4],dichromat_pal("DarkRedtoBlue.12")(12)[9:12])
  cols2 <- c(dichromat_pal("BluetoOrange.10")(10))
  
  popuse <- mean(data$pop)
  timeuse <- seq(2016,2100,1/52)[1:length(subsetI)]
  plot(timeuse,subsetI/popuse,type="n",col="#F64740",lwd =2,ylim=c(0,0.005),xlab="",ylab=
         "",bty = "n",xaxs="i",xlim=c(2016,2035), main="Texas", yaxt ="n", xaxt= "n")
  polygon(c(timeuse[52*4+controlWeekStart],timeuse[52*4+controlWeekStart],timeuse[52*4+controlWeekStart + 52],
            timeuse[52*4+controlWeekStart + 52]), c(-1,1,1,-1), border = NA,col="grey88")
  polygon(c(timeuse[52*4+controlWeekStart + 104],timeuse[52*4+controlWeekStart + 104],max(timeuse),
            max(timeuse)), c(-1,1,1,-1), border = NA,col="grey88")
  axis(side = 2, col.axis= cols[7], at = c(0.000,0.001,0.002,0.003), label = c("0.000","0.001","0.002","0.003"))
  lines(timeuse,subsetI/popuse,col=cols[7],lwd =2,ylim=c(0,0.01),xlab="",ylab=
          "")
  subsetI2 <- pred$I$mean[times > 39.5]/data2$pop[1]
  rm <- rollmean(subsetI2, k = 52 )[1:length(timeuse)]
  lines(timeuse, rm, col="grey30", lwd = 2, lty = 2)
  mtext(side= 2, "I/N", line = 2,col = cols[7], at = 0.001)
  mtext(side= 2, "mean(I/N)", line = 2,col="grey30", at =0.002)
  axis(side = 1, at = seq(2015,2035,5), labels =  seq(2015,2035,5))
  abline(v = seq(2016,2100,1),lty=2,col="gray")
  
  text(2030, 0.005,expression('"RSV like": High R'[0]*', SIR'))
  
  par(new = TRUE)
  plot(timeuse, subsetS/popuse, col=cols[2],type="l",xlab="",ylab="", 
       axes = F,xaxs="i", ylim=c(0.1,0.3), lwd = 2, lty = 1,xlim=c(2016,2035))
  axis(side = 4, col.axis=cols[2], at = c(0.15,0.20,0.25), label = c(0.15,"0.20",0.25))
  title(xlab="Year", line = 2)
  mtext(side = 4, line = 2, 'S/N',col=cols[2])
  
  par(new = TRUE)
  subsetS2 <- pred$S$mean[times > 39.5]
  beta2 = c(data2$sea_beta[26:52] , rep(data2$sea_beta, length = length(subsetS2)))[1:length(subsetS2)]
  beta2[245:length(beta2)] <- beta2[245:length(beta2)]*0.8
  reff = beta2*subsetS2
  rmreff <- rollmean(reff, k = 52 )[1:length(timeuse)]
  plot(timeuse, rmreff, ylab = " ", xlab = "", yaxt = "n", xaxt = "n", ylim=c(-2.5,1.5), type="l", lwd =2, col =cols2[9],xlim=c(2016,2035),bty = "n",
       xaxs="i", lty = 2)
  axis(side = 2, col.axis=cols2[9], at = c(0.5,1,1.5), label = c(0.5,1,1.5))
  
  mtext(side = 2, line = 2, expression('mean(R'[e]*')'),col=cols2[9], at= 1, lty  = 1)
  
}
