# NUSTATYMAI -------------------------------------------------------------------

rm(list = ls())
Sys.setlocale(locale = "Lithuanian")

library("rotations")
library(CircStats)
library(circular)
require("caTools")

# VON MISES DENSITY FUNCTION ---------------------------------------------------

pdf("vonmises density.pdf", paper = "special")
par(mar=c(0,0,0,0))
r <- seq(-pi, pi, length = 500)
plot(r, dvmises(r, kappa = 4, Haar = FALSE), type = "l", ylab=' ', 
     ylim=c(0,0.8), axes=FALSE, xlab=' ')
axis(1, pos=0, at=c(-pi, -pi/2, 0, pi/2, pi), labels=c("-180", "-90", "0",  "90", "180"))
axis(2, pos=0, labels=TRUE, at=c(0.2,0.4,0.6,0.8))
abline(v=0, h=0)
par(new=TRUE)
plot(r, dvmises(r, kappa = 2, Haar = FALSE), type = "l", ylab = " ",
     ylim=c(0,0.8), axes=FALSE, col="red", xlab=' ')
par(new=TRUE)
plot(r, dvmises(r, kappa = 1, Haar = FALSE), type = "l", ylab = " ",
     ylim=c(0,0.8), axes=FALSE, col="blue", xlab=' ')
par(new=TRUE)
plot(r, dvmises(r, kappa = 0.5, Haar = FALSE), type = "l", ylab = " ",
     ylim=c(0,0.8), axes=FALSE, col="green", xlab=' ')
dev.off()

curve(dvmises(x, kappa = 4, Haar = FALSE), -pi,pi)
curve(dvmises(x, kappa = 0.5, Haar = FALSE), -pi,pi, add=TRUE)

# SKLAIDOS GRAFIKO PAVYZDYS --------------------------------------------------

pdf("roses.pdf", paper = "special")
par(mar=c(0,0,0,0))
windc <- circular(wind, type="angles",units="radians",template="geographics")
plot(windc, cex=1.5, bin=720, stack=TRUE, sep=0.035, shrink=1.3)
ticks.circular(circular(seq(0,2*pi,pi/8)), zero=pi/2, rotation="clock", tcl=0.075)
rose.diag(windc, bins=16, col="darkgrey", cex=1.5, prop=1.3, add=TRUE)
dev.off()

# DUOMENYS APIE VËJO KRYPTÁ -------------------------------------------------
setwd("C:/Users/User/Desktop/bakalauras")
duom <- read.csv("Preila2010_meteo.csv", header = T)
head(duom)
glod <- duom[,c(1,2,6)]
v.k <- duom[,6]
v.k <- as.numeric(as.vector(v.k), na.rm=T)
a <- which(is.na(v.k)==TRUE)
v.kdeg <- v.k[-a]
length(v.kdeg)
v.k <- v.k[-a]
v.k <- v.k*pi/180 #pasiverciam i radianus
# circular laipsniais 
v.kcl <- circular(v.kdeg, type = "angles", units="degrees", template = "geographics", zero=pi/2)
# circular radianais
v.kc <- circular(v.k, type = "angles", units="radians", template = "geographics", zero=0)

#patikrinam ar turime kryptinius duomenis
is.circular(v.kc) 
is.circular(v.kcl) 

# GRAFINIS VAIZDAVIMAS ------------------------------------------------------

# Kryptine statistika 
pdf("KSgrafikai.pdf", paper = "special")
par(mfrow=c(1,1))
par(mar=c(0,0,0,0)+0.1)
plot(v.kc, cex=0.15, bin=720, stack = TRUE, sep=0.035, shrink = 1.3, col="black")
rose.diag(v.kc, bins=16, col="darkgrey", cex=1.5, prop=1.3, add=TRUE)
lines(density.circular(v.kcl, bw=40), lwd=2)
dev.off()

# Klasikinë statistika
h<-10
sl.vid <- runmean(v.kcl, h, alg = "C")
sl.vidr <- sl.vid*pi/180
b <- range(v.k)/8
b[2]
br <- seq(0, 6.28144, b[2])

labi1 <- c("0",
expression(paste(pi,"/4")),
expression(paste(pi,"/2")),
expression(paste("3",pi,"/4")),
expression(paste(pi)),
expression(paste("5",pi, "/4")),
expression(paste("3",pi, "/2")),
expression(paste("7",pi,"/4")),
expression(paste("2",pi)))

pdf("TSgrafikaislvid.pdf",width=8,height=4,paper='special')
par(mar=c(0,0,0,0)+2)
par(mfrow=c(1,2))
plot(sl.vidr, cex=0.1,
     type = "l",frame.plot=FALSE,  pch=16, yaxt="n",xlab="Stebinio numeris", ylab="Vëjo kryptis radianais", xlim = c(0,20000))
axis(side=2, at=c(0,1,2,3,4,5,6),  labels=c("0", expression(paste(pi,"/3")), expression(paste("2",pi, "/3")), 
                          expression(pi), expression(paste("4",pi, "/3")), expression(paste("5",pi,"/3")), expression(paste("2",pi))))
hist(v.k, freq = FALSE, ylim = c(0,0.4),
     main=" ",ylab="Tankis", xlab="Vëjo kryptis radianais", breaks=br, xaxt="n")
axis(side=1, at=br, cex.axis=0.9,
     labels=labi1)
lines(density(sl.vidr), col="red")
dev.off()


# Apraðomoji statistika -----------------------------------------------------

# SUSIKURIAM REIKALINGAS FUNKCIJAS:

# asimetrija
hatsC <- function(bbar2,V) {bbar2/(V**(3/2))}
# ekscesas
hatkC <- function(abar2,V, Rbar) {(abar2-Rbar**4)/(V**2)}
# kryptinë dispersija
delhatf <- function(R, R2) { (1-R2)/(2*R^2)}


# REZULTATAI -----------------------------------------------------------

theta <- mean(v.kcl)
Rbar <- rho.circular(v.kcl)
phi <- median.circular(v.kcl)
V <- 1-rho.circular(v.kcl)
v <- sd.circular(v.kcl)
w <- range.circular(v.kcl)
D0 <- meandeviation(v.kcl)
delhat <- delhatf(rho.circular(v.kcl),trigonometric.moment(v.kcl, p=2)$rho)
theta1 <- trigonometric.moment(v.kcl, p=1)$mu
a1 <- trigonometric.moment(v.kcl, p=1)$cos
b1 <- trigonometric.moment(v.kcl, p=1)$sin
R1 <- trigonometric.moment(v.kcl, p=1)$rho
m_1 <- trigonometric.moment(v.kcl, p=1)$cos+(-1)*trigonometric.moment(v.kcl, p=1)$sin
theta2 <- trigonometric.moment(v.kcl, p=2)$mu
a2 <- trigonometric.moment(v.kcl, p=2)$cos
b2 <- trigonometric.moment(v.kcl, p=2)$sin
R2 <- trigonometric.moment(v.kcl, p=2)$rho
m_2 <- trigonometric.moment(v.kcl, p=2)$cos+(-1)*trigonometric.moment(v.kcl, p=2)$sin
abar1 <- trigonometric.moment(v.kcl, p=1, control.circular = list(units="degrees"), center = TRUE)$cos
bbar1 <- trigonometric.moment(v.kcl, p=1, control.circular = list(units="degrees"), center = TRUE)$sin
m1 <- trigonometric.moment(v.kcl, p=1, control.circular = list(units="degrees"), center = TRUE)$cos+(-1)*trigonometric.moment(v.kcl, p=1, control.circular = list(units="degrees"), center = TRUE)$sin
abar2 <- trigonometric.moment(v.kcl, p=2, control.circular = list(units="degrees"), center = TRUE)$cos
bbar2 <- trigonometric.moment(v.kcl, p=2, control.circular = list(units="degrees"), center = TRUE)$sin
m2 <- trigonometric.moment(v.kcl, p=2, control.circular = list(units="degrees"), center = TRUE)$cos+(-1)*trigonometric.moment(v.kcl, p=2, control.circular = list(units="degrees"), center = TRUE)$sin
hats <- hatsC(trigonometric.moment(v.kcl, p=2, control.circular = list(units="degrees"), center = TRUE)$sin, 1-rho.circular(v.kcl))
hatk <- hatkC(trigonometric.moment(v.kcl, p=2, control.circular = list(units="degrees"), center = TRUE)$cos, 1-rho.circular(v.kcl), rho.circular(v.kcl))

theta # vidutinio kampo kryptis
Rbar # vidutinio kampo krypties vektoriaus ilgis
phi # medianos kryptis
V # imties kryptinis nuokrypis (klas. st. atitinka dispersija)
v # imties kryptinis standartinis nuokrypis
delhat # imties kryptine dispersija
D0 # kyrptinio vidurkio skirtumas
w # kryptinis lankas
theta1 # 1trig. momentas apie 0
a1 # 1trig. momentas apie 0
b1 # 1trig. momentas apie 0
R1 # 1trig. momentas apie 0
m_1 # 1trig. momentas apie 0
theta2 # 2trig. momentas apie 0
a2 # 2trig. momentas apie 0
b2 # 2trig. momentas apie 0
R2 # 2trig. momentas apie 0
m_2 # 2trig. momentas apie 0
abar1  # 1trig. momentas apie vid.
bbar1 # 1trig. momentas apie vid.
m1 # 1trig. momentas apie vid.
abar2 # 2trig. momentas apie vid.
bbar2 # 2trig. momentas apie vid.
m2 # 2trig. momentas apie vid.
hats # asimetrijos koeficientas 
hatk # eksceso koeficientas

# PASISKIRSTYMAS ------------------------------------------------------------
kuiper.test(v.kcl, alpha = 0.05)
watson.test(v.kcl, dist = "vonmises", alpha = 0.05)
rao.spacing.test(v.kcl, alpha=0.05) 

# Rayleigh testas
R <- 0.12
l <- 2*n*R^2
S <- (1-1/(2*n))*2*n*R^2+(n*R^4)/2
S
l

# Tankio grafikas
plot(v.kcl, shrink=1.2, stack=TRUE, pch=16, bins=720, cex=0.2)
lines(density.circular(v.kcl, bw=40), lwd=2)
rose.diag(v.kcl, bins=24, cex=1.5, prop=2.3, col="grey", add=TRUE)

# PAPILDOMI GRAFIKAI ---------------------------------------------------------

# grafikai su rodyklëmis
pdf("rodykles.pdf", paper = "special")
par(mfrow=c(1,1))
par(mar=c(2,0,2,0)+2)
plot(v.kcl, shrink=1.2, stack=TRUE, pch=16, bins=720, cex=0.2)
lines(density.circular(v.kcl, bw=40), lwd=2)
rose.diag(v.kcl, bins=24, cex=1.5, prop=2.3, col="grey", add=TRUE)
# vidutinio kampo kryptis 117,71 raudona
arrows.circular(mean(v.kcl), y=1, lwd=2, col="red")
# vidurkis 173,78 melyna
arrows.circular(circular(mean(v.kdeg), type = "angles", units="degrees", template = "geographics", zero=pi/2), y=1, lwd=2, col="blue")
# mediana Kryp 112.5 zalia
arrows.circular(median.circular(v.kcl), y=1, lwd=2, col="green")
# mediana tiesineje statistikoje 158,6 juoda
arrows.circular(circular(median(v.kdeg), type = "angles", units="degrees", template = "geographics", zero=pi/2), y=1, lwd=2, col="black")
dev.off()


# SLENKANCIOS CHARAKTERISTIKOS -----------------------------------------------

# funkcija
h <- 480 # 10 dienø
n <- length(v.kcl)

# slenkancios charakteristikos funkcija
karpymas <- function(i, h, x){
  w <- x[i:(i+h-1)]
}

l <- lapply(1:(n-h+1), karpymas, h=480, x=v.kdeg)

# vidutinë kampo kryptis
a<-head(sort(vid.sl, decreasing = TRUE))
b <- vector()
for(i in 1:6){
  b[i] <- which(vid.sl==a[i])
}
a1<-head(sort(vid.sl, decreasing = FALSE))
b1 <- vector()
for(i in 1:6){
  b1[i] <- which(vid.sl==a1[i])
}
vid.sl <- sapply(l, function(x)  mean(circular(x, type = "angles",
          units="degrees", template = "geographics", zero=pi/2, modulo="2pi")))
length(vid.sl)
at.men <- c(0,31*48, 31*48+28*48, 31*48+31*48+28*48, 30*48+31*48+31*48+28*48,
            31*48+30*48+31*48+31*48+28*48,30*48+31*48+30*48+31*48+31*48+28*48,
            31*48+30*48+31*48+30*48+31*48+31*48+28*48,
            31*48+31*48+30*48+31*48+30*48+31*48+31*48+28*48,
            30*48+31*48+31*48+30*48+31*48+30*48+31*48+31*48+28*48,
            31*48+30*48+31*48+31*48+30*48+31*48+30*48+31*48+31*48+28*48,
            30*48+31*48+30*48+31*48+31*48+30*48+31*48+30*48+31*48+31*48+28*48,
            31*48+30*48+31*48+30*48+31*48+31*48+30*48+31*48+30*48+31*48+31*48+28*48)
labels.men <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                "11", "12", "")

pdf("vidsl.pdf", paper = "special")
plot(vid.sl, type = "l",xaxt="n", frame.plot=FALSE, xlim=c(0, 17000),
     xlab = "Mënesis",
     ylab="Vidutinio kampo kryptis laipsniais", ylim=c(-0, 360), yaxt="n")
axis(side=1,at=at.men, labels = labels.men)
axis(side=2,at=seq(0, 360, 90))
abline(h=180, col="red", lwd=2)
abline(h=c(360, 0), col="blue", lwd=2)
abline(h=270, col="violet", lwd=2)
abline(h=90, col="green", lwd=2)
dev.off()

# Vidutinio kampo krypties vektoriaus ilgis
rho.sl <- sapply(l, function(x)  rho.circular(circular(x, type = "angles",
                  units="degrees", template = "geographics", zero=pi/2)))

rho.max <- which(rho.sl==max(rho.sl))
rho.min <- which(rho.sl==min(rho.sl))

pdf("rhosl.pdf",paper='special')
plot(rho.sl,xaxt="n", type = "l", frame.plot=FALSE, xlim=c(0, 17000),
     xlab = "Mënesis", ylab="Vidutinio kampo krypties vektoriaus ilgis")
axis(side=1,at=at.men, labels = labels.men)
points(x=rho.max, y=max(rho.sl), col = "red", pch=19)
points(x=rho.min, y=min(rho.sl), col = "blue", pch=19)
dev.off()

# medianos kryptis
med.sl <- sapply(l, function(x)  
  median.circular(circular(x, type = "angles", units="degrees", 
             template = "geographics", zero=pi/2, modulo="2pi")))
pdf("medsl.pdf",paper='special')
plot(med.sl,xaxt="n", type = "l", frame.plot=FALSE, xlim=c(0, 17000),
     ylim=c(0,360), xlab = "Mënesis",ylab="Medianos kryptis laipsniais", 
     yaxt="n")
axis(at=at.men,labels = labels.men, side=1)
axis(at=seq(0,360,60), side=2)
abline(h=180, col="red", lwd=2)
abline(h=c(360, 0), col="blue", lwd=2)
abline(h=270, col="violet", lwd=2)
abline(h=90, col="green", lwd=2)
dev.off()

# kryptinis nuokrypis 
V.sl <- 1-rho.sl
V.max <- which(V.sl==max(V.sl))
V.min <- which(V.sl==min(V.sl))

pdf("Vsl.pdf", paper='special')
plot(V.sl, type = "l",xaxt="n", frame.plot=FALSE, xlim=c(0, 17000),
     xlab = "Mënesis", ylab = "Kryptinis nuokrypis")
axis(at=at.men,labels = labels.men, side=1)
points(x=V.max, y=max(V.sl), col = "red", pch=19)
points(x=V.min, y=min(V.sl), col = "blue", pch=19)
dev.off()

# vidrukiu skirtumas 
vidsk.sl <- sapply(l, function(x) meandeviation(circular(x,
            type = "angles", units="degrees", template = "geographics",
            zero=pi/2)))
vidsk.max <- which(vidsk.sl==max(vidsk.sl))
vidsk.min <- which(vidsk.sl==min(vidsk.sl))

pdf("vidsksl.pdf",paper='special')
plot(vidsk.sl, xaxt="n", xlim=c(0, 17000), type = "l", frame.plot=FALSE,
     xlab="Mënesis", ylab="Kryptinio vidurkio skirtumas", 
     ylim=c(0, 1.6), yaxt="n")
axis(side=1,at=at.men,labels = labels.men)
axis(side=2,at=seq(0, 1.6, 0.2))
points(x=vidsk.max, y=max(vidsk.sl), col = "red", pch=19)
points(x=vidsk.min, y=min(vidsk.sl), col = "blue", pch=19)
dev.off()

# kryptinis lankas 
krl.sl <- sapply(l, function(x) range.circular(circular(x, type = "angles",
                    units="degrees", template = "geographics", zero=pi/2)))
krl.max <- which(krl.sl==max(krl.sl))
krl.min <- which(krl.sl==min(krl.sl))

pdf("krlsl.pdf" ,paper='special')
plot(krl.sl,xaxt="n", xlim=c(0, 17000), type = "l", frame.plot=FALSE, 
     xlab="Mënesis", ylab="Kryptinis lankas laipsniais", ylim=c(0,360), yaxt="n")
axis(side=1,at=at.men,labels = labels.men)
axis(side=2,at=seq(0, 360, 60))
points(x=krl.max, y=rep(max(krl.sl),length(krl.max)) , col = "red", pch=19)
points(x=krl.min, y=rep(min(krl.sl),length(krl.min)), col = "blue", pch=19)
dev.off()

# SAUSIS IR LIEPA ----------------------------------------------------------
glod$Data <- as.Date(glod$Data, format = "%m/%d/%Y") 
df <- data.frame(data = glod$Data,
                 metai = as.numeric(format(glod$Data , format = "%Y")),
                 men = as.numeric(format(glod$Data, format = "%m")),
                 diena = as.numeric(format(glod$Data, format = "%d")))
glod <- cbind(glod, metai=0)
glod$metai <-df[match(glod$Data, df$data),2]
glod <- cbind(glod, men=0)
glod$men <-df[match(glod$Data, df$data),3]
glod <- cbind(glod, diena=0)
glod$diena <-df[match(glod$Data, df$data),4]
names(glod)[3] <- "vejok"
glod$vejok <- as.numeric(as.vector(glod$vejok), na.rm=TRUE)

vk.men <- glod[,c(3, 5)]
# Imam sausio ir liepos men
vk.sausis <- subset(vk.men, men==1, vejok)
vk.sausis <- vk.sausis[,1]
vk.liepa <- subset(vk.men, men==7, vejok)
vk.liepa <- vk.liepa[,1]
vk.liepa <- vk.liepa[-which(is.na(vk.liepa)==TRUE)]
vk.liepaC <- circular(vk.liepa, type = "angles", units="degrees", template = "geographics", zero=pi/2, modulo = "2pi")
vk.sausisC <- circular(vk.sausis, type = "angles", units="degrees", template = "geographics", zero=pi/2, modulo = "2pi")

# Grafikai

# SAUSIS
pdf("KSsausis.pdf",width = 10, height = 8 ,paper = "special")
par(mfrow=c(1,1))
par(mar=c(0,0,0,0)+0.1)
plot(vk.sausisC, cex=0.7, axes=FALSE, bin=720, stack = TRUE, sep=0.035, shrink = 1.3, col="black")
rose.diag(vk.sausisC, bins=16, col="darkgrey", cex=1.5, prop=1.3, add=TRUE)
lines(density.circular(vk.sausisC, bw=40), lwd=2)
dev.off()

# LIEPA
pdf("KSliepa.pdf", width = 10, height = 8 ,paper = "special")
par(mfrow=c(1,1))
par(mar=c(0,0,0,0)+0.1)
plot(vk.liepaC,cex=0.7,axes=FALSE, bin=720, stack = TRUE, sep=0.035, shrink = 1.3, col="black")
rose.diag(vk.liepaC, bins=16, col="darkgrey", cex=1.5, prop=1.3, add=TRUE)
lines(density.circular(vk.liepaC, bw=40), lwd=2)
dev.off()

# Sausios ir liepos kartu
pdf("KSliepasausis.pdf", width = 8, height = 12 ,paper = "special")
par(mfrow=c(2,1))
par(mar=c(0,0,0,0)+1)
plot(vk.sausisC,main="SAUSIS", axes=FALSE,cex=0.7, bin=720, stack = TRUE, sep=0.035, shrink = 1.3, col="black")
rose.diag(vk.sausisC, bins=16, col="darkgrey", cex=1.5, prop=1.3, add=TRUE)
lines(density.circular(vk.sausisC, bw=40), lwd=2)
plot(vk.liepaC, main="LIEPA",cex=0.7,axes=FALSE, bin=720, stack = TRUE, sep=0.035, shrink = 1.3, col="black")
rose.diag(vk.liepaC, bins=16, col="darkgrey", cex=1.5, prop=1.3, add=TRUE)
lines(density.circular(vk.liepaC, bw=40), lwd=2)
dev.off()

pdf("ls.pdf", width = 8, height = 6 ,paper = "special")
par(mfrow=c(1,2))
par(mar=c(0,0,0,2))
plot(vk.sausisC, axes=FALSE,cex=0.7, bin=720, stack = TRUE, sep=0.035, shrink = 1.3, col="black")
rose.diag(vk.sausisC, bins=16, col="darkgrey", cex=1.5, prop=1.3, add=TRUE)
lines(density.circular(vk.sausisC, bw=40), lwd=2)
plot(vk.liepaC,cex=0.7,axes=FALSE, bin=720, stack = TRUE, sep=0.035, shrink = 1.3, col="black")
rose.diag(vk.liepaC, bins=16, col="darkgrey", cex=1.5, prop=1.3, add=TRUE)
lines(density.circular(vk.liepaC, bw=40), lwd=2)
dev.off()


# Skaièiuojam charakteristikas:

# priskiriam tà, kuriam mënesiui norim skaièiuoti
x <- vk.liepaC
x <- vk.sausisC

theta <- mean(x)
Rbar <- rho.circular(x)
phi <- median.circular(x)
V <- 1-rho.circular(x)
v <- sd.circular(x)
w <- range.circular(x)
D0 <- meandeviation(x)
delhat <- delhatf(rho.circular(x),trigonometric.moment(x, p=2)$rho)
theta1 <- trigonometric.moment(x, p=1)$mu
a1 <- trigonometric.moment(x, p=1)$cos
b1 <- trigonometric.moment(x, p=1)$sin
R1 <- trigonometric.moment(x, p=1)$rho
m_1 <- trigonometric.moment(x, p=1)$cos+(-1)*trigonometric.moment(x, p=1)$sin
theta2 <- trigonometric.moment(x, p=2)$mu
a2 <- trigonometric.moment(x, p=2)$cos
b2 <- trigonometric.moment(x, p=2)$sin
R2 <- trigonometric.moment(x, p=2)$rho
m_2 <- trigonometric.moment(x, p=2)$cos+(-1)*trigonometric.moment(x, p=2)$sin
abar1 <- trigonometric.moment(x, p=1, control.circular = list(units="degrees"), center = TRUE)$cos
bbar1 <- trigonometric.moment(x, p=1, control.circular = list(units="degrees"), center = TRUE)$sin
m1 <- trigonometric.moment(x, p=1, control.circular = list(units="degrees"), center = TRUE)$cos+(-1)*trigonometric.moment(x, p=1, control.circular = list(units="degrees"), center = TRUE)$sin
abar2 <- trigonometric.moment(x, p=2, control.circular = list(units="degrees"), center = TRUE)$cos
bbar2 <- trigonometric.moment(x, p=2, control.circular = list(units="degrees"), center = TRUE)$sin
m2 <- trigonometric.moment(x, p=2, control.circular = list(units="degrees"), center = TRUE)$cos+(-1)*trigonometric.moment(x, p=2, control.circular = list(units="degrees"), center = TRUE)$sin
hats <- hatsC(trigonometric.moment(x, p=2, control.circular = list(units="degrees"), center = TRUE)$sin, V)
hatk <- hatkC(trigonometric.moment(x, p=2, control.circular = list(units="degrees"), center = TRUE)$cos, 1-rho.circular(x), rho.circular(x))

theta # vidutinio kampo kryptis
Rbar # vidutinio kampo krypties vektoriaus ilgis
phi # medianos kryptis
V # imties kryptinis nuokrypis (klas. st. atitinka dispersija)
v # imties kryptinis standartinis nuokrypis
delhat # imties kryptine dispersija
D0 # kyrptinio vidurkio skirtumas
w # kryptinis lankas
theta1 # 1trig. momentas apie 0
a1 # 1trig. momentas apie 0
b1 # 1trig. momentas apie 0
R1 # 1trig. momentas apie 0
m_1 # 1trig. momentas apie 0
theta2 # 2trig. momentas apie 0
a2 # 2trig. momentas apie 0
b2 # 2trig. momentas apie 0
R2 # 2trig. momentas apie 0
m_2 # 2trig. momentas apie 0
abar1  # 1trig. momentas apie vid.
bbar1 # 1trig. momentas apie vid.
m1 # 1trig. momentas apie vid.
abar2 # 2trig. momentas apie vid.
bbar2 # 2trig. momentas apie vid.
m2 # 2trig. momentas apie vid.
hats # asimetrijos koeficientas 
hatk # eksceso koeficientas








