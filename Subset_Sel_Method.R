install.packages("readr")
install.packages("leaps")
install.packages("caret")
library(leaps)
library(readr)
library(caret)
dat = read_csv("Work10_final.csv", col_types = cols(Abb = col_skip(), polymer = col_skip()))
View(dat)

preproc1 <- preProcess(dat, method=c("center", "scale"))
dat_norm <- predict(preproc1, dat)
dat = dat_norm

MSE_store = matrix(as.vector(0), nrow=9, ncol=8, byrow = TRUE)
rownames(MSE_store) = c("Sol", "Tm", "ST", "MV", "density", "Ecoh", "Tg", "HC", "EMW")
colnames(MSE_store) = c(1:8)
View(MSE_store)

attach(dat)

# For Solubility (Sol)
set.seed(123)
train = sample(102, 75)

regfit.full=regsubsets(Sol~.,dat)
summary(regfit.full)
reg.summary=summary(regfit.full)

lm.fit2 = lm(Sol~density, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Sol-pred)[-train]^2)
MSE_store [1,1] = MSE 
lm.fit2 = lm(Sol~MV+Ecoh, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Sol-pred)[-train]^2)
MSE_store [1,2] = MSE 
lm.fit2 = lm(Sol~MV+density+Ecoh, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Sol-pred)[-train]^2)
MSE_store [1,3] = MSE 
lm.fit2 = lm(Sol~MV+Ecoh+Tg+HC, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Sol-pred)[-train]^2)
MSE_store [1,4] = MSE 
lm.fit2 = lm(Sol~MV+density+Ecoh+Tg+HC, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Sol-pred)[-train]^2)
MSE_store [1,5] = MSE 
lm.fit2 = lm(Sol~ST+MV+density+Ecoh+Tg+HC, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Sol-pred)[-train]^2)
MSE_store [1,6] = MSE 
lm.fit2 = lm(Sol~ST+MV+density+Ecoh+Tg+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Sol-pred)[-train]^2)
MSE_store [1,7] = MSE 
lm.fit2 = lm(Sol~Tm+ST+MV+density+Ecoh+Tg+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Sol-pred)[-train]^2)
MSE_store [1,8] = MSE

plot(MSE_store[1,], type = 'l', lwd=2, xlab = "Number of predictors", ylab = "Test MSE")
View(MSE_store)

# For Melting Temp (Tm)

regfit.full=regsubsets(Tm~.,dat)
summary(regfit.full)
reg.summary=summary(regfit.full)

################### MSE
train = sample(102, 75)
lm.fit2 = lm(Tm~Tg, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Tm-pred)[-train]^2)
MSE_store [2,1] = MSE 
lm.fit2 = lm(Tm~Tg+HC, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Tm-pred)[-train]^2)
MSE_store [2,2] = MSE 
lm.fit2 = lm(Tm~MV+Tg+HC, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Tm-pred)[-train]^2)
MSE_store [2,3] = MSE 
lm.fit2 = lm(Tm~MV+Tg+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Tm-pred)[-train]^2)
MSE_store [2,4] = MSE 
lm.fit2 = lm(Tm~MV+density+Tg+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Tm-pred)[-train]^2)
MSE_store [2,5] = MSE 
lm.fit2 = lm(Tm~Sol+MV+density+Tg+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Tm-pred)[-train]^2)
MSE_store [2,6] = MSE 
lm.fit2 = lm(Tm~Sol+MV+density+Ecoh+Tg+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Tm-pred)[-train]^2)
MSE_store [2,7] = MSE 
lm.fit2 = lm(Tm~Sol+ST+MV+density+Ecoh+Tg+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Tm-pred)[-train]^2)
MSE_store [2,8] = MSE 
plot(MSE_store[2,], type = 'l', lwd=2, xlab = "Number of predictors", ylab = "Test MSE")

View(MSE_store)
###########################

# For Surface Tension (ST)

regfit.full=regsubsets(ST~.,dat)
summary(regfit.full)
reg.summary=summary(regfit.full)

################### MSE
train = sample(102, 75)
lm.fit2 = lm(ST~Sol, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((ST-pred)[-train]^2)
MSE_store [3,1] = MSE 
lm.fit2 = lm(ST~Sol+density, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((ST-pred)[-train]^2)
MSE_store [3,2] = MSE 
lm.fit2 = lm(ST~Sol+density+Tg, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((ST-pred)[-train]^2)
MSE_store [3,3] = MSE 
lm.fit2 = lm(ST~Sol+density++Ecoh+Tg, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((ST-pred)[-train]^2)
MSE_store [3,4] = MSE 
lm.fit2 = lm(ST~Sol+density++Ecoh+Tg+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((ST-pred)[-train]^2)
MSE_store [3,5] = MSE 
lm.fit2 = lm(ST~Sol+density++Ecoh+Tg+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((ST-pred)[-train]^2)
MSE_store [3,6] = MSE 
lm.fit2 = lm(ST~Sol+Tm+density++Ecoh+Tg+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((ST-pred)[-train]^2)
MSE_store [3,7] = MSE 
lm.fit2 = lm(ST~Sol+Tm+MV+density++Ecoh+Tg+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((ST-pred)[-train]^2)
MSE_store [3,8] = MSE 

plot(MSE_store[3,], type = 'l', lwd=2, xlab = "Number of predictors", ylab = "Test MSE")
View(MSE_store)

# for Molar Volume 
regfit.full=regsubsets(MV~.,dat)
summary(regfit.full)
reg.summary=summary(regfit.full)

################### MSE
train = sample(102, 75)
lm.fit2 = lm(MV~HC, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((MV-pred)[-train]^2)
MSE_store [4,1] = MSE 
lm.fit2 = lm(MV~Sol+Ecoh, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((MV-pred)[-train]^2)
MSE_store [4,2] = MSE 
lm.fit2 = lm(MV~Sol+Ecoh+HC, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((MV-pred)[-train]^2)
MSE_store [4,3] = MSE 
lm.fit2 = lm(MV~Sol+Ecoh+Tg+HC, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((MV-pred)[-train]^2)
MSE_store [4,4] = MSE 
lm.fit2 = lm(MV~Sol+Tm+Ecoh+Tg+HC, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((MV-pred)[-train]^2)
MSE_store [4,5] = MSE 
lm.fit2 = lm(MV~Sol+Tm+density+Ecoh+Tg+HC, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((MV-pred)[-train]^2)
MSE_store [4,6] = MSE 
lm.fit2 = lm(MV~Sol+Tm+density+Ecoh+Tg+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((MV-pred)[-train]^2)
MSE_store [4,7] = MSE 
lm.fit2 = lm(MV~Sol+Tm+ST+density+Ecoh+Tg+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((MV-pred)[-train]^2)
MSE_store [4,8] = MSE 

plot(MSE_store[4,], type = 'l', lwd=2, xlab = "Number of predictors", ylab = "Test MSE")
View(MSE_store)

# For Density (density)

regfit.full=regsubsets(density~.,dat)
summary(regfit.full)
reg.summary=summary(regfit.full)

################### MSE
train = sample(102, 75)
lm.fit2 = lm(density~Sol, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((density-pred)[-train]^2)
MSE_store [5,1] = MSE 
lm.fit2 = lm(density~Sol+ST, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((density-pred)[-train]^2)
MSE_store [5,2] = MSE 
lm.fit2 = lm(density~Sol+ST+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((density-pred)[-train]^2)
MSE_store [5,3] = MSE 
lm.fit2 = lm(density~Sol+Tm+Tg+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((density-pred)[-train]^2)
MSE_store [5,4] = MSE 
lm.fit2 = lm(density~Sol+Tm+ST+Tg+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((density-pred)[-train]^2)
MSE_store [5,5] = MSE 
lm.fit2 = lm(density~Sol+Tm+ST+Ecoh+Tg+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((density-pred)[-train]^2)
MSE_store [5,6] = MSE 
lm.fit2 = lm(density~Sol+Tm+ST+MV+Ecoh+Tg+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((density-pred)[-train]^2)
MSE_store [5,7] = MSE 
lm.fit2 = lm(density~Sol+Tm+ST+MV+Ecoh+Tg+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((density-pred)[-train]^2)
MSE_store [5,8] = MSE 

plot(MSE_store[5,], type = 'l', lwd=2, xlab = "Number of predictors", ylab = "Test MSE")
View(MSE_store)


# For Energy of Cohesion (Ecoh)

regfit.full=regsubsets(Ecoh~.,dat)
summary(regfit.full)
reg.summary=summary(regfit.full)

################### MSE
train = sample(102, 75)
lm.fit2 = lm(Ecoh~MV, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Ecoh-pred)[-train]^2)
MSE_store [6,1] = MSE 
lm.fit2 = lm(Ecoh~Sol+MV, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Ecoh-pred)[-train]^2)
MSE_store [6,2] = MSE 
lm.fit2 = lm(Ecoh~Sol+MV+HC, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Ecoh-pred)[-train]^2)
MSE_store [6,3] = MSE 
lm.fit2 = lm(Ecoh~Sol+MV+Tg+HC, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Ecoh-pred)[-train]^2)
MSE_store [6,4] = MSE 
lm.fit2 = lm(Ecoh~Sol+MV+Tg+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Ecoh-pred)[-train]^2)
MSE_store [6,5] = MSE 
lm.fit2 = lm(Ecoh~Sol+MV+density+Tg+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Ecoh-pred)[-train]^2)
MSE_store [6,6] = MSE 
lm.fit2 = lm(Ecoh~Sol+ST+MV+density+Tg+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Ecoh-pred)[-train]^2)
MSE_store [6,7] = MSE 
lm.fit2 = lm(Ecoh~Sol+Tm+ST+MV+density+Tg+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Ecoh-pred)[-train]^2)
MSE_store [6,8] = MSE 

plot(MSE_store[6,], type = 'l', lwd=2, xlab = "Number of predictors", ylab = "Test MSE")
View(MSE_store)

# for Glass Transition Temperature (Tg)
regfit.full=regsubsets(Tg~.,dat)
summary(regfit.full)
reg.summary=summary(regfit.full)

train=sample(102,75)
lm.fit2 = lm(Tg~Tm, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Tg-pred)[-train]^2)
MSE_store [7,1] = MSE 
lm.fit2 = lm(Tg~Sol+Tm, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Tg-pred)[-train]^2)
MSE_store [7,2] = MSE 
lm.fit2 = lm(Tg~Sol+MV+HC, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Tg-pred)[-train]^2)
MSE_store [7,3] = MSE 
lm.fit2 = lm(Tg~Sol+MV+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Tg-pred)[-train]^2)
MSE_store [7,4] = MSE 
lm.fit2 = lm(Tg~Sol+Tm+MV+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Tg-pred)[-train]^2)
MSE_store [7,5] = MSE 
lm.fit2 = lm(Tg~Sol+Tm+MV+density+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Tg-pred)[-train]^2)
MSE_store [7,6] = MSE 
lm.fit2 = lm(Tg~Sol+Tm+MV+density+Ecoh+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Tg-pred)[-train]^2)
MSE_store [7,7] = MSE 
lm.fit2 = lm(Tg~Sol+Tm+ST+MV+density+Ecoh+HC+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((Tg-pred)[-train]^2)
MSE_store [7,8] = MSE

plot(MSE_store[7,], type = 'l', lwd=2, xlab = "Number of predictors", ylab = "Test MSE")
View(MSE_store)


# for Molar Heat Capacity (HC)

regfit.full=regsubsets(HC~.,dat)
summary(regfit.full)
reg.summary=summary(regfit.full)

################### MSE
train = sample(102, 75)
lm.fit2 = lm(HC~MV, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((HC-pred)[-train]^2)
MSE_store [8,1] = MSE 
lm.fit2 = lm(HC~MV+Tg, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((HC-pred)[-train]^2)
MSE_store [8,2] = MSE 
lm.fit2 = lm(HC~Tm+MV+Tg, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((HC-pred)[-train]^2)
MSE_store [8,3] = MSE 
lm.fit2 = lm(HC~Sol+Tm+MV+Tg, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((HC-pred)[-train]^2)
MSE_store [8,4] = MSE 
lm.fit2 = lm(HC~Sol+Tm+MV+Ecoh+Tg, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((HC-pred)[-train]^2)
MSE_store [8,5] = MSE 
lm.fit2 = lm(HC~Sol+Tm+MV+Ecoh+Tg+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((HC-pred)[-train]^2)
MSE_store [8,6] = MSE 
lm.fit2 = lm(HC~Sol+Tm+ST+MV+Ecoh+Tg+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((HC-pred)[-train]^2)
MSE_store [8,7] = MSE 
lm.fit2 = lm(HC~Sol+Tm+ST+MV+density+Ecoh+Tg+EMW, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((HC-pred)[-train]^2)
MSE_store [8,8] = MSE 

plot(MSE_store[8,], type = 'l', lwd=2, xlab = "Number of predictors", ylab = "Test MSE")
View(MSE_store)


# for Entanglement Molecular Weight (EMW)

regfit.full=regsubsets(EMW~.,dat)
summary(regfit.full)
reg.summary=summary(regfit.full)

################### MSE

train = sample(102, 75)
lm.fit2 = lm(EMW~Tg, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((EMW-pred)[-train]^2)
MSE_store [9,1] = MSE 
lm.fit2 = lm(EMW~Ecoh+Tg, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((EMW-pred)[-train]^2)
MSE_store [9,2] = MSE 
lm.fit2 = lm(EMW~MV+density+Tg, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((EMW-pred)[-train]^2)
MSE_store [9,3] = MSE 
lm.fit2 = lm(EMW~MV+density+Tg+HC, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((EMW-pred)[-train]^2)
MSE_store [9,4] = MSE 
lm.fit2 = lm(EMW~Sol+density+Ecoh+Tg+HC, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((EMW-pred)[-train]^2)
MSE_store [9,5] = MSE 
lm.fit2 = lm(EMW~Sol+Tm+density+Ecoh+Tg+HC, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((EMW-pred)[-train]^2)
MSE_store [9,6] = MSE 
lm.fit2 = lm(EMW~Sol+Tm+ST+density+Ecoh+Tg+HC, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((EMW-pred)[-train]^2)
MSE_store [9,7] = MSE 
lm.fit2 = lm(EMW~Sol+Tm+ST+MV+density+Ecoh+Tg+HC, data = dat, subset = train)
pred=predict(lm.fit2, dat)
MSE = mean((EMW-pred)[-train]^2)
MSE_store [9,8] = MSE 

plot(MSE_store[9,], type = 'l', lwd=2, xlab = "Number of predictors", ylab = "Test MSE")
View(MSE_store)


