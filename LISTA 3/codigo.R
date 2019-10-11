# Lista 3 multivariada

setwd("/home/gui/Área de Trabalho/Lista 3 M")

library(psych)
library(expm) # Para cálculo de Matriz^1/2
library(CCA)
library(candisc)

# Exercício 1

# item a

r <- matrix(c(1,-0.488,0.15,-0.488,1,-0.13,0.15,-0.13,1),3,3) # matriz de correlação
L <- matrix(c(0.75,-0.65,0.2),3,1) # matrix de cargas fatoriais

L.Lt <- round(L%*%t(L),3)

Psi <- r - L.Lt # especificidade

# item b

# Comunalidades
L[1]^2
L[2]^2
L[3]^2

# Exercício 2

# item a

auto_v <- eigen(r)

round(t(auto_v$values),3) # autovalores
round(auto_v$vectors,3) # autovetores

# item b

L.hat <- as.matrix(sqrt(round(t(auto_v$values),3)[1])*round(auto_v$vectors,3)[,1])

LLt.hat <- L.hat %*% t(L.hat)

Psi.hat <- r - LLt.hat # especificidade

# item c

auto_v$values[1]/dim(r)[1]

# Exercício 3

F1 <- c(0.602,0.467,0.926,1,0.874,0.894)
F2 <- c(0.2,0.154,0.143,0,0.476,0.327)
F1.rot <- c(0.484,0.375,0.603,0.519,0.861,0.744)
F2.rot <- c(0.411,0.319,0.717,0.855,0.499,0.594)

C <- numeric() # vetor de comunalidades

for (i in 1:6){
  C[i] <- F1[i]^2+F2[i]^2
}

C.rot <- numeric() # vetor de comunalidades rotacioanda 

for (i in 1:6){
  C.rot[i] <- F1.rot[i]^2+F2.rot[i]^2
}

# item c

# Exercício 4

dados <- read.table("T1-9.dat")
names(dados) <- c("Pais","100m","200m","400m","800m","1500m","3000m","Maratona")

# item a

screeplot(princomp(dados[-1], cor=F),
          type="lines",
          pch=16,
          col="turquoise4",
          main=" Scree Plot")
abline(h=1,lty=3)

Cov <- round(cov(dados[2:8]),3) # Matriz de covariâncias

Fat <- principal(Cov, nfactors=2, rotate="none",covar = T)
Fat

# item b

screeplot(princomp(dados[-1], cor=T),
          type="lines",
          pch=16,
          col="turquoise4",
          main=" Scree Plot")
abline(h=1,lty=3)

Cor <- round(cor(dados[2:8]),3) # Matriz de correlações

Fat2 <- principal(Cor, nfactors=1, rotate="none")
Fat2

# Exercício 5

dados.m <- dados

# conversão dos dados para velocidades em (m/s)
dados.m$`100m` <- round(100/dados.m$`100m`,2)
dados.m$`200m` <- round(200/dados.m$`200m`,2)
dados.m$`400m` <- round(400/dados.m$`400m`,2)
dados.m$`800m` <- round(800/(dados.m$`800m`/60),2)
dados.m$`1500m` <- round((1500/(dados.m$`1500m`/60)),2)
dados.m$`3000m` <- round((3000/(dados.m$`3000m`/60)),2)
dados.m$Maratona <- round((42.195/(dados.m$Maratona/60)),2)

Cov <- round(cov(dados.m[2:8]),3) # matrix de covariâncias 

screeplot(princomp(dados.m[-1], cor=F),
          type="lines",
          pch=16,
          col="turquoise4",
          main=" Scree Plot")
abline(h=1,lty=3)

# Exercício 6

Cov <- matrix(c(8,2,3,1,
                2,5,-1,3,
                3,-1,6,-2,
                1,3,-2,7),4,4)

r11 <- Cov[1:2,1:2]
r12 <- Cov[1:2,-(1:2)]
r21 <- Cov[-(1:2),1:2]
r22 <- Cov[3:4,3:4]

E1 <- sqrtm(solve(r11))%*%r12%*%solve(r22)%*%r21%*%sqrtm(solve(r11))
E2 <- sqrtm(solve(r22))%*%r21%*%solve(r11)%*%r12%*%sqrtm(solve(r22))

VVM <- eigen(E1)$vectors
VVN <- eigen(E2)$vectors

sqrt(eigen(E1)$values)

# item b

# Coeficientes de U
A <- -sqrtm(solve(r11))%*%VVM

# Coeficientes de V
B <- -sqrtm(solve(r22))%*%VVN

# item c

M <- solve(r11)%*%r12%*%solve(r22)%*%r21
eigen(M)$values

M.sqrt <- sqrtm(solve(r11))%*%r12%*%solve(r22)%*%r21%*%sqrtm(solve(r11))
eigen(M.sqrt)$values

# Exercício 7

# item a

Cor <- matrix(c(1,0.6328,0.2412,0.0586,
                0.6328,1,-0.0553,0.0655,
                0.2412,-0.0553,1,0.4248,
                0.0586,0.0655,0.4248,1),4,4)

r11 <- Cor[1:2,1:2]
r12 <- Cor[1:2,-(1:2)]
r21 <- Cor[-(1:2),1:2]
r22 <- Cor[3:4,3:4]

E1 <- sqrtm(solve(r11))%*%r12%*%solve(r22)%*%r21%*%sqrtm(solve(r11))
E2 <- sqrtm(solve(r22))%*%r21%*%solve(r11)%*%r12%*%sqrtm(solve(r22))

VVM <- eigen(E1)$vectors
VVN <- eigen(E2)$vectors

sqrt(eigen(E1)$values)

# item b

p <- dim(r12)[1]
q <- dim(r12)[2]
s_11 <- det(r11)
s_22 <- det(r22)
S <- det(Cor)
lb <- 140*log((s_11*s_22)/S) # estatística do teste RV

1-pchisq(lb,p*q) # p-value

# teste se a primeira correlação canônica é igual a zero

r1 <- sqrt(eigen(E1)$values)[1] # primeira correlação
gl <- 2 # graus de liberdade
n <- 140 # tamanho da amostra

t.test <- r1/sqrt((1-r1^2)/(n-gl)) # estatística do teste 

2*(1-pnorm(t.test)) # p-value bilateral

# item c

# Coeficientes de U
A <- -sqrtm(solve(r11))%*%VVM

# Coeficientes de V
B <- -sqrtm(solve(r22))%*%VVN

# item d

corU1 <- A%*%r11
corU2 <- A%*%r12
corV1 <- B%*%r21
corV2 <- B%*%r22

# Exercício 8

dados <- read.table("T7-7.dat")
names(dados) <- c("BL","EM","SF","BS","AFL","LFF","FFF","ZST")

corcan <- cc(scale(dados[,1:4]),scale(dados[,5:8]))

round(t(corcan$cor),3) # correlação canonica

corcan2 <- cancor(scale(dados[,1:4]),scale(dados[,5:8]))
