# Lista 5 - Multivariada - MAE0330

library(ggplot2)
library(ca)
library(lavaan)
library(expm)
library(tidyr)

# exercício 1

# leitura dos dados
M1 <- t(matrix(c(8968,3289,11933,228,2,
                   2881,464,3187,109,63,
                   3033,5468,7919,367,264,
                   2415,1646,3328,380,350,
                   1042,277,1061,198,57,
                   522,2496,2250,382,421,
                   258,1483,77,1016,655,
                   952,2112,2200,305,13,
                   479,613,186,764,149,
                   1541,280,230,1369,228,
                   169,1355,1214,290,5,
                   2300,256,2357,125,27,
                   328,40,305,65,7,
                   4049,94,2515,1660,20,
                   16524,45,6791,8986,824,
                   2031,41,2024,40,0,
                   1834,41,1818,38,1),5,17))

dados1 <- as.data.frame(M1)
names(dados1)<- c("I", "TR", "AL", "TE", "O")
row.names(dados1) <- c("Ob","Gin","C.Cir","Hemato","Hemod",
                       "CMG","CTIA","P","UTI1","UTI2","UCI",
                       "CC","CMO","PA","PP","DQ","PPsq")

# teste qui-quadrado
chisq.test(dados1)

# Analise de correspondecnia
A.corresp <- ca(dados1)
summary(A.corresp)

# para gráfico

n <- sum(M1)
P <- M1/n

vr <- P%*%rep(1,ncol(M1))
vc <-  t(P)%*%rep(1,nrow(M1))
Dr <- diag(c(vr))
Dc <- diag(c(vc))

S <-  sqrtm(solve(Dr))%*%(P-vr%*%t(vc))%*%sqrtm(solve(Dc))
dvs_S <- svd(S)
M_X1 <-  sqrtm(solve(Dr))%*%dvs_S$u[,1:2]%*%diag(dvs_S$d[1:2])
M_Y2 <- sqrtm(solve(Dc))%*%dvs_S$v[,1:2]%*%diag(dvs_S$d[1:2])

dados_grafico <- data.frame(rbind(M_X1,M_Y2))
dados_grafico$tipo <- factor(c(rep(1,17),rep(2,5)))
row.names(dados_grafico) <- c("Ob","Gin","C.Cir","Hemato","Hemod",
                       "CMG","CTIA","P","UTI1","UTI2","UCI",
                       "CC","CMO","PA","PP","DQ","PPsq","I",
                       "TR", "AL", "TE", "O")

dados_grafico %>%
  ggplot(aes(x = X1, y = X2, shape=factor(tipo), colour=factor(tipo))) + 
  geom_point() +
  geom_hline(yintercept = 0, colour = "gray70") + 
  geom_vline(xintercept = 0, colour = "gray70") +
  geom_text(aes(y=X2+0.08,label= rownames(dados_grafico))) +
  guides(colour=FALSE, shape=FALSE) +
  labs(x="Dim 1 (59.1%)",
       y="Dim 2 (36.8%)",
       title="Clinicas e atendimentos realizados")

# Exercício 2

# Entrada dos dados

c1 <- matrix(c("7-5","Padrao","2011"),417,3,byrow = TRUE)
c2 <- matrix(c("4","Padrao","2011"),36,3,byrow = TRUE)
c3 <- matrix(c("7-5","Teste","2011"),164,3,byrow = TRUE)
c4 <- matrix(c("4","Teste","2011"),176,3,byrow = TRUE)
c5 <- matrix(c("<4","Teste","2011"),90,3,byrow = TRUE)

c6 <- matrix(c("7-5","Padrao","2012"),357,3,byrow = TRUE)
c7 <- matrix(c("4","Padrao","2012"),27,3,byrow = TRUE)
c8 <- matrix(c("7-5","Teste","2012"),169,3,byrow = TRUE)
c9 <- matrix(c("4","Teste","2012"),161,3,byrow = TRUE)
c10 <- matrix(c("<4","Teste","2012"),54,3,byrow = TRUE)

c11 <- matrix(c("7-5","Padrao","2013"),800,3,byrow = TRUE)
c12 <- matrix(c("4","Padrao","2013"),240,3,byrow = TRUE)
c13 <- matrix(c("<4","Padrao","2013"),103,3,byrow = TRUE)
c14 <- matrix(c("7-5","Teste","2013"),412,3,byrow = TRUE)
c15 <- matrix(c("4","Teste","2013"),458,3,byrow = TRUE)
c16 <- matrix(c("<4","Teste","2013"),274,3,byrow = TRUE)

c17 <- matrix(c("7-5","Padrao","2014"),273,3,byrow = TRUE)
c18 <- matrix(c("4","Padrao","2014"),176,3,byrow = TRUE)
c19 <- matrix(c("<4","Padrao","2014"),39,3,byrow = TRUE)
c20 <- matrix(c("7-5","Teste","2014"),185,3,byrow = TRUE)
c21 <- matrix(c("4","Teste","2014"),220,3,byrow = TRUE)
c22 <- matrix(c("<4","Teste","2014"),83,3,byrow = TRUE)

c23 <- matrix(c("7-5","Padrao","2015"),1521,3,byrow = TRUE)
c24 <- matrix(c("4","Padrao","2015"),1794,3,byrow = TRUE)
c25 <- matrix(c("<4","Padrao","2015"),585,3,byrow = TRUE)
c26 <- matrix(c("7-5","Teste","2015"),1420,3,byrow = TRUE)
c27 <- matrix(c("4","Teste","2015"),1681,3,byrow = TRUE)
c28 <- matrix(c("<4","Teste","2015"),635,3,byrow = TRUE)

dados2 <- rbind(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,
                c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,
                c21,c22,c23,c24,c25,c26,c27,c28)

dados2 <- data.frame(dados2)
names(dados2) <- c("Tam_Bulbilhos","Tratamento","Ano")

table(dados2)

# Analise de correspondecnia
A.corresp.M <- mjca(dados2,lambda="indicator")

summary(A.corresp.M)

g <- plot(A.corresp.M, labels=c(2,2), map="symmetric")

dados_grafico1 <- data.frame(g$cols)
row.names(dados_grafico1) <- row.names(g$cols)

dados_grafico1 %>%
  ggplot(aes(x = Dim1, y = Dim2)) + 
  geom_point() +
  geom_hline(yintercept = 0, colour = "gray70") + 
  geom_vline(xintercept = 0, colour = "gray70") +
  scale_x_continuous(limits = c(-2,1.5)) +
  geom_text(aes(y=Dim2+0.08,label= rownames(dados_grafico1))) +
  guides(colour=FALSE, shape=FALSE) +
  labs(x="Dim 1 (17.8%)",
       y="Dim 2 (15.3%)",
       title="Bulbilhos de alho por tamanho tratamento e ano de plantio")


# Exercício 3

library(lavaan)
M.cor <- data.frame(matrix(c(1.00,0.37,0.42,0.53,0.38,0.81,0.35,0.42,0.40,0.24,
                  0.37,1.00,0.33,0.14,0.10,0.34,0.65,0.32,0.14,0.15,
                  0.42,0.33,1.00,0.38,0.20,0.49,0.20,0.75,0.39,0.17,
                  0.53,0.14,0.38,1.00,0.24,0.58,-0.04,0.46,0.73,0.15,
                  0.38,0.10,0.20,0.24,1.00,0.32,0.11,0.26,0.19,0.43,
                  0.81,0.34,0.49,0.58,0.32,1.00,0.34,0.46,0.55,0.24,
                  0.35,0.65,0.20,-0.04,0.11,0.34,1.00,0.18,0.06,0.15,
                  0.42,0.32,0.75,0.46,0.26,0.46,0.18,1.00,0.54,0.20,
                  0.40,0.14,0.39,0.73,0.19,0.55,0.06,0.54,1.00,0.16,
                  0.24,0.15,0.17,0.15,0.43,0.24,0.15,0.20,0.16,1.00),10,10))

names(M.cor) <- c("Vt1","St1","Rt1","Nt1","Wt1","Vt2","St2","Rt2","Nt2","Wt2")
rownames(M.cor) <- c("Vt1","St1","Rt1","Nt1","Wt1","Vt2","St2","Rt2","Nt2","Wt2")

lower <- '
                  1.00,
                  0.37,1.00,
                  0.42,0.33,1.00,
                  0.53,0.14,0.38,1.00,
                  0.38,0.10,0.20,0.24,1.00,
                  0.81,0.34,0.49,0.58,0.32,1.00,
                  0.35,0.65,0.20,-0.04,0.11,0.34,1.00,
                  0.42,0.32,0.75,0.46,0.26,0.46,0.18,1.00,
                  0.40,0.14,0.39,0.73,0.19,0.55,0.06,0.54,1.00,
                  0.24,0.15,0.17,0.15,0.43,0.24,0.15,0.20,0.16,1.00'
# convert to a full symmetric covariance matrix with names
wheaton.cov <- getCov(lower, names=c("Vt1","St1","Rt1","Nt1","Wt1",
                                     "Vt2","St2","Rt2","Nt2","Wt2"))
model <- '
 # measurement model
 Ano1 =~ Vt1 + St1 + Rt1 + Nt1 + Wt1
 Ano2 =~ Vt2 + St2 + Rt2 + Nt2 + Wt2
 # regressions
 Ano2 ~ Ano1
 # residual correlations
 Vt1 ~~ Vt2
 St1 ~~ St2
 Rt1 ~~ Rt2
 Nt1 ~~ Nt2
 Wt1 ~~ Wt2
'

fit <- sem(model, sample.cov=wheaton.cov, sample.nobs=100)
summary(fit, standardized=TRUE)

library(semPlot)
semPaths(fit)