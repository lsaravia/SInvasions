
## 
orig_data <- read.table("sl02Pomac.out", header=T)
attach(orig_data)

orig_data <- read.table("sl01Pomac.out", header=T)
attach(orig_data)

orig_data <- read.table("sl03Pomac.out", header=T)
attach(orig_data)

orig_data <- read.table("sl04Pomac.out", header=T)
attach(orig_data)

## Eliminar datos con 0 en densidades 
orig_data0 <- subset(orig_data, AdGrad1>0); nrow(orig_data0)
attach(orig_data0)
layout(matrix(c(1,2,3,4,5,6,7,8), 2, 4, byrow = TRUE))
hist(g)
hist(dj0)
hist(dj1)
hist(dj2)
hist(alfa)
hist(b)
hist(da)
hist(densIni)



hist(JuvGrad1)
log_JuvGrad1 <- log(JuvGrad1)
hist(log_JuvGrad1)

model <- lm(log_JuvGrad1 ~ g+dj0+dj1+alfa+b+da+densIni)
summary(model)

library(car)
vif(model)

y_hats <- fitted(model)
residuals <- resid(model)
plot(y_hats, residuals)

model <- lm(AdOcupado14 ~ g+dj0+dj1+alfa+b+da+densIni)

sq_AdOcupado14 <- sqrt(AdOcupado14)
hist(sq_AdOcupado14)

model <- lm(sq_AdOcupado14 ~ g+dj0+dj1+alfa+b+da+densIni)
summary(model)

## Escalar variables dependientes e independiente al intervalo [0,1]
##
g_esc <- (g - min(g))/max(g)
dj0_esc <- (dj0 - min(dj0))/max(dj0)
dj1_esc <- (dj1 - min(dj1))/max(dj1)
alfa_esc <- (alfa - min(alfa))/max(alfa)
b_esc <- (b - min(b))/max(b)
da_esc <- (da - min(da))/max(da)
densIni_esc <- (densIni - min(densIni))/max(densIni)
sq_AdOcupado14_esc <- (sq_AdOcupado14 - min(sq_AdOcupado14))/max(sq_AdOcupado14)
model <- lm(sq_AdOcupado14_esc ~ g_esc+dj0_esc+dj1_esc+alfa_esc+b_esc+da_esc+densIni_esc)
summary(model)

log_JuvGrad1_esc <- (log_JuvGrad1 - min(log_JuvGrad1))/max(log_JuvGrad1)
model <- lm(log_JuvGrad1_esc ~ g_esc+dj0_esc+dj1_esc+alfa_esc+b_esc+da_esc+densIni_esc)
summary(model)


## Subset para seleccionar datos de los criterios

crit1 <- subset(orig_data, AdGrad1>205 & AdGrad1<620); nrow(crit1)
crit2 <- subset(orig_data, AdGrad2>1100 & AdGrad2<2270); nrow(crit2)
crit3 <- subset(orig_data,JuvGrad1>20 & JuvGrad1<450); nrow(crit3)
crit4 <- subset(orig_data,JuvGrad2>600 & JuvGrad2<2280); nrow(crit4)
crit5 <- subset(orig_data,AdOcupado14>283 & AdOcupado14<301); nrow(crit5)
crit6 <- subset(orig_data,AdOcupado23>349 & AdOcupado23<366); nrow(crit6)
critSin35 <- subset(orig_data, (AdGrad1>205 & AdGrad1<620) & (AdGrad2>1100 & AdGrad2<2270) & (JuvGrad2>600 & JuvGrad2<2280) & (AdOcupado23>349 & AdOcupado23<366)); nrow(critSin35)
critSin5 <- subset(orig_data, (AdGrad1>205 & AdGrad1<620) & (AdGrad2>1100 & AdGrad2<2270) & (JuvGrad1>20 & JuvGrad1<450) & (JuvGrad2>600 & JuvGrad2<2280) & (AdOcupado23>349 & AdOcupado23<366)); nrow(critSin5)
critSin3 <- subset(orig_data, (AdGrad1>205 & AdGrad1<620) & (AdGrad2>1100 & AdGrad2<2270) & (JuvGrad2>600 & JuvGrad2<2280) & (AdOcupado14>283 & AdOcupado14<301) & (AdOcupado23>349 & AdOcupado23<366)); nrow(critSin3)
crit56 <- subset(orig_data, (AdOcupado14>283 & AdOcupado14<301) & (AdOcupado23>349 & AdOcupado23<366)); nrow(crit56)


## Guardar tabla de datos con un criterio
write.table(critSin35, file="critSin35.txt", quote=FALSE, sep="\t")
write.table(crit3, file="crit3.txt", quote=FALSE, sep="\t")
write.table(crit5, file="crit5.txt", quote=FALSE, sep="\t")

write.table(crit6, file="crit6.txt")

## Grafico de boxplot de subconjunto de parámetros segun criterios
##
## Primera simulacion con termino mortalidad juveniles linear
ciInva01 <- read.table("ciInva01.txt", header=T)
attach(ciInva01)
svg("figPomac01.svg", 9,7)
layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
x11(,9,7)
svg("figPomac02.svg", 9,7)
layout(matrix(c(1,2,3,4,5,6,7,8), 2, 4, byrow = TRUE))
boxplot(g ~crit,main='g',col='bisque')
boxplot(dj0 ~crit,main='d0',col='bisque')
boxplot(dj1 ~crit,main='d1',col='bisque')
boxplot(alfa ~crit,main='alfa',col='bisque')
boxplot(b ~crit,main='B',col='bisque')
boxplot(da ~crit,main='D',col='bisque')
boxplot(densIni ~crit,main='densIni',col='bisque')
dev.off()

## Interaccion entre parametros para los criterios C5 y CS35
##
library(lattice)
ciInva01a <- subset(ciInva01, crit!="C3"); nrow(ciInva01a)
svg("fig03Inter.svg", 9,7)
trellis.par.set(superpose.symbol=list(pch=c(2, 3))) 
sup.sym <- list(col=trellis.par.get()[["superpose.symbol"]]$col[2:3],pch = c(3,2))
splom(ciInva01a[2:7], groups=crit, data = ciInva01, xlab = "", pscales=0, 
varnames = c("g", "d0","d1","alfa", "B", "D"),
key = list(
points = sup.sym,
text = list(c("C5", "CS35")), space="right"))
dev.off()

## Interaccion entre parametros para todos los criterios
##
trellis.par.set(superpose.symbol=list(pch=c(1,2,3))) 
splom(ciInva01[2:7], groups=crit, data = ciInva01, xlab = "", pscales=0, 
varnames = c("g", "d0","d1","alfa", "B", "D"),
auto.key = list(space = "right", points = TRUE))


##
##
## Segunda simulacion con termino mortalidad juveniles Allee
ciInva02 <- read.table("ciInva02.txt", header=T)
attach(ciInva02)
svg("figPomac04.svg", 9,7)
x11(,9,7)
layout(matrix(c(1,2,3,4,5,6,7,8), 2, 4, byrow = TRUE))
boxplot(g ~crit,main='g',col='bisque')
boxplot(dj0 ~crit,main='d0',col='bisque')
boxplot(dj1 ~crit,main='d1',col='bisque')
boxplot(dj2 ~crit,main='d2',col='bisque')
boxplot(alfa ~crit,main='alfa',col='bisque')
boxplot(b ~crit,main='B',col='bisque')
boxplot(da ~crit,main='D',col='bisque')
boxplot(densIni ~crit,main='densIni',col='bisque')
dev.off()

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
boxplot(g ~crit,main='g',col='bisque')
boxplot(dj0 ~crit,main='dj0',col='bisque')
boxplot(dj1 ~crit,main='dj1',col='bisque')
boxplot(dj2 ~crit,main='dj2',col='bisque')


## Tercera simulacion con termino mortalidad juveniles Allee + dispersion larga distancia potencial inversa
ciInva03 <- read.table("ciInva03.txt", header=T)
attach(ciInva03)
svg("fig05.svg", 9,7)
x11(,9,7)
layout(matrix(c(1,2,3,4,5,6,7,8), 2, 4, byrow = TRUE))
boxplot(g ~crit,main='g',col='bisque')
boxplot(dj0 ~crit,main='d0',col='bisque')
boxplot(dj1 ~crit,main='d1',col='bisque')
boxplot(dj2 ~crit,main='d2',col='bisque')
boxplot(alfa ~crit,main='alfa',col='bisque')
boxplot(alfa1 ~crit,main='alfa1',col='bisque')
boxplot(b ~crit,main='B',col='bisque')
boxplot(da ~crit,main='D',col='bisque')


## Interaccion entre parametros para los criterios C5 y CS35
##
library(lattice)
ciInva03a <- subset(ciInva03, crit!="C3"); nrow(ciInva03a)
svg("fig06Inter.svg", 9,7)
trellis.par.set(superpose.symbol=list(pch=c(2,3))) 
sup.sym <- list(col=trellis.par.get()[["superpose.symbol"]]$col[2:3],pch = c(3,2))
splom(ciInva03a[2:9], groups=crit, data = ciInva01, xlab = "", pscales=0, 
varnames = c("g", "d0","d1","d2","alfa","alfa1", "B", "D"),
key = list(
points = sup.sym,
text = list(c("C5", "CS35")), space="right"))
dev.off()

## Cuarta simulacion con termino mortalidad juveniles Allee + dispersion larga distancia uniforme
##
ciInva04 <- read.table("ciInva04.txt", header=T)
attach(ciInva04)
svg("fig06.svg", 9,7)
x11(,9,7)
layout(matrix(c(1,2,3,4,5,6,7,8), 2, 4, byrow = TRUE))
boxplot(g ~crit,main='g',col='bisque')
boxplot(dj0 ~crit,main='d0',col='bisque')
boxplot(dj1 ~crit,main='d1',col='bisque')
boxplot(dj2 ~crit,main='d2',col='bisque')
boxplot(alfa ~crit,main='alfa',col='bisque')
boxplot(alfa1 ~crit,main='beta',col='bisque')
boxplot(b ~crit,main='B',col='bisque')
boxplot(da ~crit,main='D',col='bisque')
dev.off()

## Comparacion de modelos 
##
library(lattice)
compMod <- read.table("compMod.txt", header=T)
attach(compMod)
trellis.par.set(superpose.symbol=list(pch=c(1,2,3))) 
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
dotplot(JuvGrad1 ~ factor(modelo), data = compMod, groups = crit, 
auto.key = list(space = "right", points = TRUE))

dotplot(AdOcupado14 ~ factor(modelo), data = compMod, groups = crit, 
auto.key = list(space = "right", points = TRUE))

trellis.par.set(superpose.symbol=list(pch=c(1,2,3))) 
dotplot(PorcSimul ~ factor(modelo), data = compMod, groups = crit, 
auto.key = list(space = "right", points = TRUE))

trellis.par.set(superpose.symbol=list(pch=c(1,2,3))) 
dotplot(AdOcupado23 ~ factor(modelo), data = compMod, groups = crit, 
auto.key = list(space = "right", points = TRUE))
