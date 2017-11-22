############################################################################################
#Replication code for Figure 2 (main model) in Manatschal, Anita and Julian Bernauer (2016):  
#"Consenting to Exclude? Empirical Patterns of Democracy and Immigrant Integration Policy", 
#West European Politics 39(2): 183-204.
#Published on GitHub on November 22, 2017

setwd("...")

library(foreign)
library(R2WinBUGS)
library(R2jags)
library(rjags)

mipex_av <- read.dta("mipex_av.dta")
mipex_inst <- read.dta("mipex_institutions.dta")
attach(mipex_av)
attach(mipex_inst)
ls(mipex_av)
ls(mipex_inst)

N <- length(cab2)
dir <- dir-mean(dir)
J <- length(mipex10)

mipex.data <- list(J=J, N=N, idmip=idmip, par=par, elec=elec, cab2=cab2, fed=fed, dec=dec, bic=bic, const=const, jud=jud, mipex10=mipex10, dir=dir, 
                   loggdp=loggdp, hist_settler=hist_settler,  right_gov=right_gov, pop_elec=pop_elec, immigrants=immigrants, exeleg=exeleg)

mipex.parameters <- c("sigma.elec", "sigma.par", "sigma.dec", "sigma.dir", "sigma.exeleg", "alpha.exeleg", "alpha.dir", "alpha.elec", "alpha.par", 
                      "alpha.dec", "gamma.elec1", "gamma.dec2", "gamma.par1", "beta.fed2", "beta.bic2", "beta.const2", "beta.jud2", "beta.cab1", 
                      "gamma.dir3", "gamma.exeleg1","cx1", "cx2", "cx3", "b0", "b.exepar", "b.feduni", "b.dir", "b.loggdp", "b.rightg", "b.popelec", 
                      "b.immi", "b.settl", "sigma.mipex", "tau.mipex", "tau.fed", "tau.bic", "tau.const", "tau.jud", "tau.cab")

mipex.model1 <- "model{

#outcome model

for(k in 1:J){
mipex10[k] ~ dnorm(mu.mipex[k],tau.mipex)
mu.mipex[k] <- b0 + b.exepar*cx1[k] + b.feduni*cx2[k] + b.dir*cx3[k] + b.settl*hist_settler[k] +  b.loggdp*loggdp[k] + b.immi*immigrants[k] + 
b.rightg*right_gov[k] + b.popelec*pop_elec[k] 
}

b0 ~ dnorm(0,.001)
b.exepar ~ dnorm(0,.001)
b.feduni ~ dnorm(0,.001)
b.dir ~ dnorm(0,.001)
b.immi ~ dnorm(0,.001)
b.settl ~ dnorm(0,.001)
b.rightg ~ dnorm(0,.001)
b.popelec ~ dnorm(0,.001)
b.loggdp ~ dnorm(0,.001)
sigma.mipex ~ dunif(0,100) 
tau.mipex <- pow(sigma.mipex,-2)  


#measurement institutions
for(i in 1:N){

elec[i] ~ dnorm(mu.elec[i],tau.elec)
mu.elec[i] <- alpha.elec + gamma.elec1*x1[i]

par[i] ~ dnorm(mu.par[i],tau.par)
mu.par[i] <- alpha.par + gamma.par1*x1[i]

exeleg[i] ~ dnorm(mu.exeleg[i],tau.exeleg)
mu.exeleg[i] <- alpha.exeleg + gamma.exeleg1*x1[i]

dec[i] ~ dnorm(mu.dec[i],tau.dec)
mu.dec[i] <- alpha.dec + gamma.dec2*x2[i]

dir[i] ~ dnorm(mu.dir[i],tau.dir)
mu.dir[i] <- alpha.dir + gamma.dir3*x3[i]

cab2[i] ~ dcat(p.cab[i,1:4])
mu.cab[i] <- beta.cab1*x1[i]
logit(Q.cab[i,1]) <- tau.cab[1] - mu.cab[i]
p.cab[i,1] <- Q.cab[i,1]
logit(Q.cab[i,2]) <- tau.cab[2] - mu.cab[i]
p.cab[i,2] <- Q.cab[i,2] - Q.cab[i,1] 
logit(Q.cab[i,3]) <- tau.cab[3] - mu.cab[i]
p.cab[i,3] <- Q.cab[i,3] - Q.cab[i,2] 
p.cab[i,4] <- 1 - Q.cab[i,3]

fed[i] ~ dcat(p.fed[i,1:3])
mu.fed[i] <- beta.fed2*x2[i]
logit(Q.fed[i,1]) <- tau.fed[1] - mu.fed[i]
p.fed[i,1] <- Q.fed[i,1]
logit(Q.fed[i,2]) <- tau.fed[2] - mu.fed[i]
p.fed[i,2] <- Q.fed[i,2] - Q.fed[i,1] 
p.fed[i,3] <- 1 - Q.fed[i,2]

bic[i] ~ dcat(p.bic[i,1:4])
mu.bic[i] <- beta.bic2*x2[i]
logit(Q.bic[i,1]) <- tau.bic[1]-mu.bic[i]
p.bic[i,1] <- Q.bic[i,1]
logit(Q.bic[i,2]) <- tau.bic[2] - mu.bic[i]
p.bic[i,2] <- Q.bic[i,2] - Q.bic[i,1]
logit(Q.bic[i,3]) <- tau.bic[3] - mu.bic[i]
p.bic[i,3] <- Q.bic[i,3] - Q.bic[i,2]
p.bic[i,4] <- 1 - Q.bic[i,3]

const[i] ~ dcat(p.const[i,1:3])
mu.const[i] <- beta.const2*x2[i]
logit(Q.const[i,1]) <- tau.const[1] - mu.const[i]
p.const[i,1] <- Q.const[i,1]
logit(Q.const[i,2]) <- tau.const[2] - mu.const[i]
p.const[i,2] <- Q.const[i,2] - Q.const[i,1] 
p.const[i,3] <- 1 - Q.const[i,2]

jud[i] ~ dcat(p.jud[i,1:3])
mu.jud[i] <- beta.jud2*x2[i]
logit(Q.jud[i,1]) <- tau.jud[1] - mu.jud[i]
p.jud[i,1] <- Q.jud[i,1]
logit(Q.jud[i,2]) <- tau.jud[2] - mu.jud[i]
p.jud[i,2] <- Q.jud[i,2] - Q.jud[i,1] 
p.jud[i,3] <- 1 - Q.jud[i,2]

}

tau.elec <- pow(sigma.elec, -2)
sigma.elec ~ dunif(0, 50)
tau.par <- pow(sigma.par, -2)
sigma.par ~ dunif(0, 50)
tau.exeleg <- pow(sigma.exeleg, -2)
sigma.exeleg ~ dunif(0, 50)
tau.dec <- pow(sigma.dec, -2)
sigma.dec ~ dunif(0, 50)
tau.dir <- pow(sigma.dir, -2)
sigma.dir ~ dunif(0, 50)

alpha.elec ~ dnorm(0, .001)
alpha.par ~ dnorm(0, .001)
alpha.exeleg ~ dnorm(0, .001)
alpha.dec ~ dnorm(0, .001)
alpha.dir ~ dnorm(0, .001)

gamma.elec1 ~ dnorm(0, .001)  
gamma.par1 ~ dnorm(0, .001) I(0,)
gamma.exeleg1 ~ dnorm(0, .001) 
gamma.dec2 ~ dnorm(0, .001) 
gamma.dir3 ~ dnorm(0, .001) I(0,)

tau.cab[1] ~ dnorm(0,.01)
for (j in 1:2){
delta.cab[j] ~ dexp(2)
tau.cab[j+1] <- tau.cab[j] + delta.cab[j]
}

tau.fed[1] ~ dnorm(0,.01)
for (j in 1:1){
delta.fed[j] ~ dexp(2)
tau.fed[j+1] <- tau.fed[j] + delta.fed[j]
}

tau.bic[1] ~ dnorm(0,.01)
for (j in 1:2){
delta.bic[j] ~ dexp(2)
tau.bic[j+1] <- tau.bic[j] + delta.bic[j]
}

tau.jud[1] ~ dnorm(0,.01)
for (j in 1:1){
delta.jud[j] ~ dexp(2)
tau.jud[j+1] <- tau.jud[j] + delta.jud[j]
}

tau.const[1] ~ dnorm(0,.01)
for (j in 1:1){
delta.const[j] ~ dexp(2)
tau.const[j+1] <- tau.const[j] + delta.const[j]
}

beta.fed2 ~ dnorm(0, .001)  I(0,)
beta.bic2 ~ dnorm(0, .001) 
beta.const2 ~ dnorm(0, .001) 
beta.jud2 ~ dnorm(0, .001) 
beta.cab1 ~ dnorm(0, .001) I(0,)

for(i in 1:N){
x1[i] ~ dnorm(cx1[idmip[i]],1) 
x2[i] ~ dnorm(cx2[idmip[i]],1)  
x3[i] ~ dnorm(cx3[idmip[i]],1)
}

for(k in 1:5){
cx1[k] ~ dnorm(0,1) 
cx2[k] ~ dnorm(0,1)  
cx3[k] ~ dnorm(0,1) 
}

#restrictions Switzerland
cx1[6] ~ dnorm(0,1) I(0,)
cx2[6] ~ dnorm(0,1) I(0,)
cx3[6] ~ dnorm(0,1) I(0,)

for(k in 7:30){
cx1[k] ~ dnorm(0,1) 
cx2[k] ~ dnorm(0,1)  
cx3[k] ~ dnorm(0,1) 
}


}"

write(mipex.model1, file="mipex.model1.jags")


jags.mipex <- jags.model(file="mipex.model1.jags", data = mipex.data, n.chains = 3, n.adapt = 100)

sampleshelp <- coda.samples(jags.mipex, mipex.parameters, n.iter=100, thin=1)

samplesburn <- coda.samples(jags.mipex, mipex.parameters, n.iter=4800, thin=48)

samples <- coda.samples(jags.mipex, mipex.parameters, n.iter=5000, thin=50)


plot(sampleshelp, ask=TRUE)

plot(samples, ask=TRUE)


summary(samples[,c("gamma.elec1", "gamma.dec2", "gamma.par1", "beta.fed2", "beta.bic2", "beta.const2", "beta.jud2", "beta.cab1", "gamma.dir3", 
                   "gamma.exeleg1", "b0", "b.exepar", "b.feduni", "b.dir", "b.loggdp", "b.rightg", "b.popelec", "b.immi", "b.settl", "tau.mipex")])

dic.mod <- dic.samples(jags.mipex, 1000, "pD")
dic.mod


##################
#Creating Figure 2

kette <- as.matrix(samples)

#outcome
b0 <- kette[,"b0"] 
mb0 <- mean(b0)
sb0 <- sd(b0)

bexepar <- kette[,"b.exepar"] 
mbexepar <- mean(bexepar)
sbexepar <- sd(bexepar)

bfeduni <- kette[,"b.feduni"] 
mbfeduni <- mean(bfeduni)
sbfeduni <- sd(bfeduni)

bdir <- kette[,"b.dir"] 
mbdir <- mean(bdir)
sbdir <- sd(bdir)

bloggdp <- kette[,"b.loggdp"] 
mbloggdp <- mean(bloggdp)
sbloggdp <- sd(bloggdp)

brightg <- kette[,"b.rightg"] 
mbrightg <- mean(brightg)
sbrightg <- sd(brightg)

bpopelec <- kette[,"b.popelec"] 
mbpopelec <- mean(bpopelec)
sbpopelec <- sd(bpopelec)

bimmi <- kette[,"b.immi"] 
mbimmi <- mean(bimmi)
sbimmi <- sd(bimmi)

bsettl <- kette[,"b.settl"] 
mbsettl <- mean(bsettl)
sbsettl <- sd(bsettl)

taumipex <- kette[,"tau.mipex"] 
mtaumipex <- mean(taumipex)
staumipex <- sd(taumipex)

#patterns
gamma.elec1 <- kette[,"gamma.elec1"]
mgel1 <- mean(gamma.elec1)
sgel1 <- sd(gamma.elec1)
gamma.par1 <- kette[,"gamma.par1"]
mgp1 <- mean(gamma.par1)
sgp1 <- sd(gamma.par1)
beta.cab1 <- kette[,"beta.cab1"]
mbc1 <- mean(beta.cab1)
sbc1 <- sd(beta.cab1)
gamma.exeleg1 <- kette[,"gamma.exeleg1"]
mgex1 <- mean(gamma.exeleg1)
sgex1 <- sd(gamma.exeleg1)

beta.fed2 <- kette[,"beta.fed2"]
mbf2 <- mean(beta.fed2)
sbf2 <- sd(beta.fed2)
gamma.dec2 <- kette[,"gamma.dec2"]
mgd2 <- mean(gamma.dec2)
sgd2 <- sd(gamma.dec2)
beta.bic2 <- kette[,"beta.bic2"]
mbb2 <- mean(beta.bic2)
sbb2 <- sd(beta.bic2)
beta.jud2 <- kette[,"beta.jud2"]
mbj2 <- mean(beta.jud2)
sbj2 <- sd(beta.jud2)
beta.const2 <- kette[,"beta.const2"]
mbc2 <- mean(beta.const2)
sbc2 <- sd(beta.const2)

gamma.dir3 <- kette[,"gamma.dir3"]
mgd3 <- mean(gamma.dir3)
sgd3 <- sd(gamma.dir3)

mbsettl10 <- mbsettl/10
sbsettl10 <- sbsettl/10

mcoeffs_inds <- c(mbexepar, mbfeduni, mbdir,  mbloggdp, mbrightg, mbpopelec, mbimmi, mbsettl10,  mtaumipex, mgel1,mgp1,mbc1,mgex1,mbf2,mgd2,mbb2,mbj2,mbc2,mgd3)
scoeffs_inds <- c(sbexepar, sbfeduni, sbdir,  sbloggdp, sbrightg, sbpopelec, sbimmi, sbsettl10,  staumipex, sgel1,sgp1,sbc1,sgex1,sbf2,sgd2,sbb2,sbj2,sbc2,sgd3)

vlabels <- c("Prop. power disp." ,"Veto power disp." , "Direct power disp." ,  "GDP (log)", "Cabinet share 'right'", "Seat share 'populist'", 
             "Immigrants (share)", "Legacy: 'settler'/10", "Variance", "Disproportionality", "Number of parties", "Cabinet type", "Power of parliament", 
             "Federalism", "Decentralization","Bicameralism","Judicial review","Const. rigidity", "Direct democracy")

var.names <- c(vlabels)
m.v <- mcoeffs_inds
sd.v <- scoeffs_inds

y.axis <- length(var.names):1 
layout(matrix(c(2,1),1,2),  
    widths = c(1.5, 5)) 

par(mar=c(2,6,.5,1), lheight = .8) 
plot(m.v, y.axis, type = "p", axes = F, xlab = "", ylab = "", pch = 19, xlim = c(-6,10), cex=1, ylim = c(min(y.axis), max(y.axis)), main = "")
axis(1,at = seq(-6,10, by = 2), label = seq(-6,10, by = 2), cex.axis=.9)
axis(2, at = y.axis, label = var.names, las = 1, tick = T, font=1, cex.axis=.9)
abline(h = y.axis, lty = 2, lwd = .5, col = "grey")
segments(m.v-qnorm(.975)*sd.v, y.axis, m.v+qnorm(.975)*sd.v, y.axis, lwd =  1.5)
segments(m.v-qnorm(.9)*sd.v, y.axis -.2, m.v-qnorm(.9)*sd.v, y.axis +.2, lwd = 1.5) 
segments(m.v+qnorm(.9)*sd.v, y.axis -.2, m.v+qnorm(.9)*sd.v, y.axis +.2, lwd = 1.5)
abline(v=0, lty = 2)

par(mar=c(2,0,.5,0)) 
plot(seq(0,1,length=length(var.names)), y.axis, type = "n", axes = F,  xlab = "", ylab = "")

left.side <- .6 
segments(left.side,19,left.side,11) 
segments(left.side,19,left.side+.1,19) 
segments(left.side,11,left.side+.1,11)
text(.5, 15, "Integration policy", srt = 90, cex=.8)
segments(left.side,10,left.side,1) 
segments(left.side,10,left.side+.1,10) 
segments(left.side,1,left.side+.1,1)
text(.5, 5, "Patterns", srt = 90, cex=.8)


#Creating Figure A1 (Appendix)
#Scores of proportional power dispersion estimated 
cx11 <- kette[,"cx1[1]"]
cx12 <- kette[,"cx1[2]"]
cx13 <- kette[,"cx1[3]"]
cx14 <- kette[,"cx1[4]"]
cx15 <- kette[,"cx1[5]"]
cx16 <- kette[,"cx1[6]"]
cx17 <- kette[,"cx1[7]"]
cx18 <- kette[,"cx1[8]"]
cx19 <- kette[,"cx1[9]"]
cx110 <- kette[,"cx1[10]"]
cx111 <- kette[,"cx1[11]"]
cx112 <- kette[,"cx1[12]"]
cx113 <- kette[,"cx1[13]"]
cx114 <- kette[,"cx1[14]"]
cx115 <- kette[,"cx1[15]"]
cx116 <- kette[,"cx1[16]"]
cx117 <- kette[,"cx1[17]"]
cx118 <- kette[,"cx1[18]"]
cx119 <- kette[,"cx1[19]"]
cx120 <- kette[,"cx1[20]"]
cx121 <- kette[,"cx1[21]"]
cx122 <- kette[,"cx1[22]"]
cx123 <- kette[,"cx1[23]"]
cx124 <- kette[,"cx1[24]"]
cx125 <- kette[,"cx1[25]"]
cx126 <- kette[,"cx1[26]"]
cx127 <- kette[,"cx1[27]"]
cx128 <- kette[,"cx1[28]"]
cx129 <- kette[,"cx1[29]"]
cx130 <- kette[,"cx1[30]"]

mcx11 <- mean(cx11)
mcx12 <- mean(cx12)
mcx13 <- mean(cx13)
mcx14 <- mean(cx14)
mcx15 <- mean(cx15)
mcx16 <- mean(cx16)
mcx17 <- mean(cx17)
mcx18 <- mean(cx18)
mcx19 <- mean(cx19)
mcx110 <- mean(cx110)
mcx111 <- mean(cx111)
mcx112 <- mean(cx112)
mcx113 <- mean(cx113)
mcx114 <- mean(cx114)
mcx115 <- mean(cx115)
mcx116 <- mean(cx116)
mcx117 <- mean(cx117)
mcx118 <- mean(cx118)
mcx119 <- mean(cx119)
mcx120 <- mean(cx120)
mcx121 <- mean(cx121)
mcx122 <- mean(cx122)
mcx123 <- mean(cx123)
mcx124 <- mean(cx124)
mcx125 <- mean(cx125)
mcx126 <- mean(cx126)
mcx127 <- mean(cx127)
mcx128 <- mean(cx128)
mcx129 <- mean(cx129)
mcx130 <- mean(cx130)

scx11 <- sd(cx11)
scx12 <- sd(cx12)
scx13 <- sd(cx13)
scx14 <- sd(cx14)
scx15 <- sd(cx15)
scx16 <- sd(cx16)
scx17 <- sd(cx17)
scx18 <- sd(cx18)
scx19 <- sd(cx19)
scx110 <- sd(cx110)
scx111 <- sd(cx111)
scx112 <- sd(cx112)
scx113 <- sd(cx113)
scx114 <- sd(cx114)
scx115 <- sd(cx115)
scx116 <- sd(cx116)
scx117 <- sd(cx117)
scx118 <- sd(cx118)
scx119 <- sd(cx119)
scx120 <- sd(cx120)
scx121 <- sd(cx121)
scx122 <- sd(cx122)
scx123 <- sd(cx123)
scx124 <- sd(cx124)
scx125 <- sd(cx125)
scx126 <- sd(cx126)
scx127 <- sd(cx127)
scx128 <- sd(cx128)
scx129 <- sd(cx129)
scx130 <- sd(cx130)

mcx1scores <- c(mcx11, mcx12, mcx13, mcx14, mcx15, mcx16, mcx17, mcx18, mcx19, mcx110, mcx111, mcx112, mcx113, mcx114, mcx115, 
                mcx116, mcx117, mcx118, mcx119, mcx120, mcx121, mcx122, mcx123, mcx124, mcx125, mcx126, mcx127, mcx128, mcx129, mcx130)

scx1scores <- c(scx11, scx12, scx13, scx14, scx15, scx16, scx17, scx18, scx19, scx110, scx111, scx112, scx113, scx114, scx115, 
                scx116, scx117, scx118, scx119, scx120, scx121, scx122, scx123, scx124, scx125, scx126, scx127, scx128, scx129, scx130)


clabels <- c("Australia", "Austria", "Belgium", "Bulgaria", "Canada", "Switzerland", "Czech Republic", "Germany", "Denmark", 
             "Spain", "Estonia", "Finland", "France", "United Kingdom", "Greece", "Hungary" , "Ireland", "Italy", "Lithuania", 
             "Luxembourg",  "Latvia", "Netherlands",  "Norway", "Poland", "Portugal", "Romania", "Slovakia", "Slovenia",  "Sweden", "United States")

var.names <- c(clabels)

m.v <- mcx1scores
sd.v <- scx1scores

pic <- data.frame(var.names,m.v,sd.v)
pic.sort <- pic[order(m.v) , ]
pic.sort

y.axis <- length(var.names):1 
layout(matrix(c(2,1),1,2),  
    widths = c(1.5, 5)) 

par(mar=c(2,6,.5,1), lheight = .8) 
plot(pic.sort$m.v, y.axis, type = "p", axes = F, xlab = "", ylab = "", pch = 19, xlim = c(-4,6), cex=1, ylim = c(min(y.axis), max(y.axis)), main = "")
axis(1,at = seq(-4,6, by = 1), label = seq(-4,6, by = 1), cex.axis=.9)
axis(2, at = y.axis, label = pic.sort$var.names, las = 1, tick = T, font=1, cex.axis=.8)
abline(h = y.axis, lty = 2, lwd = .5, col = "grey")
segments(pic.sort$m.v-qnorm(.975)*pic.sort$sd.v, y.axis, pic.sort$m.v+qnorm(.975)*pic.sort$sd.v, y.axis, lwd =  1.5)

par(mar=c(2,0,.5,0)) 
plot(seq(0,1,length=length(var.names)), y.axis, type = "n", axes = F,  xlab = "", ylab = "")

left.side <- .9 
segments(left.side,30,left.side,1) 
segments(left.side,30,left.side+.1,30) 
segments(left.side,1,left.side+.1,1)
text(.5, 15, "Posterior means proportional power dispersion", srt = 90, cex=.8)

