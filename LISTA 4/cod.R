library(ggplot2)
library(MASS)
library(ISLR)
library(caret)
library(tree)
library(RcmdrMisc)
library(smacof)
library(expm)

# Exercício 1

# item a

data <- data.frame(x=c(-1, -0.5, 0, 0.5, 1, 1.5),
                   y=c(0,0,1,1,0,0),
                   l=c("f1","f2","f1","f2","f1","f2"))

ggplot(data, aes(x = x,y=y,group=l)) +
  geom_line(aes(linetype=l))+
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(-1,1.5)) +
  scale_linetype_manual(name = "funções", values = c("solid","dashed")) +
  labs(x = "x", 
       y = "y",
       title = expression(paste("Gráfico das funções ",f[1]," e ",f[2])))

# Exercício 2

# item a

mu1 <- c(0,0); mu2 <- c(0,-1); mu3 <- c(1,0)
Sigma <- matrix(c(5,-2,-2,1),2,2)

m.g <- c(1/3,-1/3)
B <- (mu1-m.g)%*%t((mu1-m.g)) + (mu2-m.g)%*%t((mu2-m.g)) + (mu3-m.g)%*%t((mu3-m.g))

eigen(solve(Sigma)%*%B)

WBW <- eigen(sqrtm(solve(Sigma))%*%B%*%sqrtm(solve(Sigma)))
f1 <- WBW$vectors[,1]
f2 <- WBW$vectors[,2]
f1t <- sqrtm(solve(Sigma))%*%f1
f2t <- sqrtm(solve(Sigma))%*%f2

# Exercício 3

data("iris")
attach(iris)

# item a

ggplot(iris,aes(x=Sepal.Width, y=Petal.Width)) + 
  geom_point(aes(shape=Species, color=Species)) +
  scale_color_manual(values=c('black','blue','red')) +
  labs(x="largura da sépala",
       y="largura da pétala",
       title="Diagrama de disperção Disperção",
       colour="Espécies",
       shape="Espécies")

# item b

iris.qda <- qda(Species~Sepal.Width+Petal.Width, data = iris,
                prior=c(1,1,1)/3,CV=FALSE)
qda.pred <- predict(iris.qda)
table(iris$Species,qda.pred$class)

dados <- data.frame(Sepal.Width,Petal.Width,Species)

# Matrizes de covariâncias
cov.set <- cov(subset(x = dados,Species=="setosa")[-3])
cov.ver <- cov(subset(x = dados,Species=="versicolor")[-3])
cov.vir <- cov(subset(x = dados,Species=="virginica")[-3])

# Médias
aux1 <- tapply(dados[,1],Species,mean)
aux2 <- tapply(dados[,2],Species,mean)

m.set <- c(aux1[1],aux2[1])
m.ver <- c(aux1[2],aux2[2])
m.vir <- c(aux1[3],aux2[3])

# Predição
Scores <- function(x,mu,s,p=1/3){
  as.numeric(-1/2*log(det(s))-1/2*t(x-mu)%*%solve(s)%*%(x-mu) + log(p))
}

x0 <- c(3.5,1.75)
list("setosa"=Scores(x0,m.set,cov.set),"versicolor"=Scores(x0,m.ver,cov.ver),
     "virginica"=Scores(x0,m.vir,cov.vir))

# item c

iris.lda <- lda(Species~Sepal.Width+Petal.Width, data = iris,prior=c(1,1,1)/3)
iris.lda

lda.pred <- predict(iris.lda)
table(iris$Species,lda.pred$class)

iris1 <- data.frame(iris$Sepal.Width,iris$Petal.Width)
lda.data <- cbind(iris1, lda.pred)
ggplot(lda.data, aes(x.LD1, x.LD2)) + 
  geom_point(aes(color = Species)) +
  scale_color_manual(values=c('black','blue','red')) +
  labs(x="LD1",
       y="LD2",
       title="Scores Discriminantes",
       colour="Espécies")

x0 <- data.frame(3.5,1.75)
names(x0) <- c("Sepal.Width","Petal.Width")
lda.predx0 <- predict(iris.lda,newdata=x0)
lda.predx0 # predição

# Gráfico de separação
mu.k <- iris.lda$means
mu <- colMeans(mu.k)
dscores <- scale(iris1[,1:2], center=mu, scale=F) 
partimat(x=dscores[,2:1], grouping=iris$Species, method="lda")

detach(iris)

# Exercício 4

set.seed(5678)

sigma <- matrix(c(20,5,5,10),2,2)
mu1 <- c(5,5)
mu2 <- c(0,0)
n <- 100
sim_bnorm1 <- mvrnorm(n, mu1, sigma)
sim_bnorm2 <- mvrnorm(n, mu2, sigma)

y <- c(rep("1",100),rep("2",100))

df <- data.frame(rbind(sim_bnorm1,sim_bnorm2),y)

# item a

sample <- createDataPartition(y=df$y, p=0.7, list=F)
train <- df[sample, ]
test <- df[-sample, ]

# item b

disc.linear1 <- lda(y~X1+X2, data = train, prior=c(1/2,1/2), CV=FALSE)
disc.linear1

lda.pred <- predict(disc.linear1,newdata=test)
table(test$y,lda.pred$class)

# gráfico de separação
mu.k <- disc.linear1$means
mu <- colMeans(mu.k)
dscores <- scale(train[,1:2], center=mu, scale=F) 
partimat(x=dscores[,2:1], grouping=train$y, method="lda")

# item c

disc.linear2 <- lda(y~X1+X2, data = train, prior=c(1,8)/9, CV=FALSE)

lda.pred <- predict(disc.linear2,newdata=test)
table(test$y,lda.pred$class)

# item d

set.seed(5678)

n <- 20
p1 <- 0.2

y <- data.frame(y=rbinom(n,1,p1))
y[which(y[,1]==0),] <- 2

sim_bnorm.d <- matrix(NA,n,2)

for (i in 1:n){
  if (y[i,]=="1"){
    sim_bnorm.d[i,] <- mvrnorm(1, mu1, sigma)
  }
  else{
    sim_bnorm.d[i,] <- mvrnorm(1, mu2, sigma)
  }
}

df1 <- data.frame(sim_bnorm.d,y)

lda.pred.d <- predict(disc.linear1,newdata=df1)
table(df1$y,lda.pred.d$class)

lda.pred.d2 <- predict(disc.linear2,newdata=df1)
table(df1$y,lda.pred.d2$class)

# Exercício 5

primate <- read.csv("primate.scapulae.txt", sep="")
attach(primate)
set.seed(5678)

sample <- createDataPartition(y=primate$classdigit, p=0.7,list=FALSE)
train <- primate[sample, ]
test <- primate[-sample, ]

model1 <- lda(class~AD.BD+AD.CD+EA.CD+Dx.CD+SH.ACR+EAD+beta, data = train, 
              prior=c(1,1,1,1,1)/5)
model1

# classificações
plot(model1,col=as.integer(train$classdigit))

# matriz de confusão
lda.pred <- predict(model1,newdata=test)
table(test$class,lda.pred$class)

# Modelo quadrático

model2 <- qda(class~AD.BD+AD.CD+EA.CD+Dx.CD+SH.ACR+EAD+beta, data = train, 
              prior=c(1,1,1,1,1)/5)
model2

# matriz de confusão
qda.pred <- predict(model2,newdata=test)
table(test$class,qda.pred$class)

# Exercício 6

data("Carseats")
attach(Carseats)

# item a

set.seed(5678)
indice_treino <- createDataPartition(y=Carseats$Sales, p=0.7, list=FALSE)
treino <- Carseats[indice_treino,]
teste <- Carseats[-indice_treino,]

# item b

tree.reg <- tree(Sales ~ ., treino)
summary(tree.reg)

# gráfico da árvore
plot(tree.reg)
text(tree.reg, pretty = 0)

# valores de preditos e real (teste)
y.hat.teste <- predict(tree.reg, teste)
y.teste <- teste$Sales

# valores de preditos e real (treino)
y.hat.treino <- predict(tree.reg, treino)
y.treino <- treino$Sales

# soma de quadrado  dados teste
sum((y.hat.teste - y.teste)^2)

# soma de quadrado  dados teste
sum((y.hat.treino - y.treino)^2)

# item d

mod1 <- lm(Sales ~ ., treino)
step <- stepwise(mod1) # stepwise 

mod.final <- lm(formula = Sales ~ CompPrice + Income + Advertising + Price + 
                  ShelveLoc + Age, data = treino)
summary(mod.final)

# soma de quadrados arvore
sum((y.hat.teste - y.teste)^2)

# soma de quadrados regressão
pred.reg <- predict(mod.final,teste)
sum((pred.reg - y.teste)^2)

# exercício 7

M <- matrix(NA,51,51)

# Matriz de similaridades
for (i in 1:51) { 
  for (j in 1:51) {
    if (i==j) M[i,j]=9
    else if (abs(i-j) >=1 && abs(i-j) <=3) M[i,j]=8
    else if (abs(i-j) >=4 && abs(i-j) <=6) M[i,j]=7
    else if (abs(i-j) >=7 && abs(i-j) <=9) M[i,j]=6
    else if (abs(i-j) >=10 && abs(i-j) <=12) M[i,j]=5
    else if (abs(i-j) >=13 && abs(i-j) <=15) M[i,j]=4
    else if (abs(i-j) >=16 && abs(i-j) <=18) M[i,j]=3
    else if (abs(i-j) >=19 && abs(i-j) <=21) M[i,j]=2
    else if (abs(i-j) >=22 && abs(i-j) <=24) M[i,j]=1
    else M[i,j]=0
  }
}

Delta <- matrix(NA,51,51)

# Matriz de dissimilaridades
for (i in 1:51) {
  for (j in 1:51) {
    Delta[i,j] = sqrt(M[i,i]+M[j,j]-2*M[i,j])
  }
}

EM.m <- cmdscale(Delta,k=2, eig=TRUE) 
summary(EM.m)

df <- data.frame(x=EM.m$points[,1],y=EM.m$points[,2])

ggplot(df,aes(x=x,y=y)) +
  geom_point() + 
  labs(x='Coordenada 1',
       y='Coordenada 2',
       title='Escalonamento multidimensional com dim=2')

# Exercício 8

# criando uma matriz simétrica
D.aux <- c(0    ,
          2.202 , 0     , 
          1.004 , 2.025 , 0     ,
          1.108 , 1.943 , 0.233 , 0     ,
          1.122 , 1.870 , 0.719 , 0.541 , 0     ,
          0.914 , 2.070 , 0.719 , 0.679 , 0.539 , 0     ,
          0.914 , 2.186 , 0.452 , 0.681 , 1.102 , 0.916 , 0     ,
          2.056 , 2.055 , 1.986 , 1.990 , 1.963 , 2.056 , 2.027 , 0     ,
          1.608 , 1.722 , 1.358 , 1.168 , 0.681 , 1.005 , 1.719 , 1.991 , 0)
D <- matrix(0,9,9)
D[upper.tri(D, diag=TRUE)] <- D.aux
D[lower.tri(D,diag=FALSE)] <- D[upper.tri(D,diag=FALSE)]

# item a

# Método Shepard-Kruskal
EM.nm_3 <- isoMDS(D, k=3, trace=FALSE)
EM.nm_4 <- isoMDS(D, k=4, trace=FALSE)
EM.nm_5 <- isoMDS(D, k=5, maxit = 60, trace=FALSE)
iso_stress <- data.frame(stress=c(EM.nm_3$stress,EM.nm_4$stress,EM.nm_5$stress),
                         dim=c(3,4,5))

ggplot(iso_stress , aes(x=dim,y=stress)) +
  geom_point() + 
  labs(x="Dimensão",
       y="Stress",
       title="Dimensão x Stress")

# item b

EM.nm_2 <- isoMDS(D, k=2)
df1 <- data.frame(EM.nm_2$points)

ggplot(df1, aes(x=X1,y=X2)) + 
  geom_point() + 
  scale_x_continuous(limits = c(-2,2)) +
  geom_text(aes(label=as.character(seq(1:9)), hjust=0.5, vjust=-0.4)) +
  labs(x="Coordenada 1",
       y="Coordenada 2",
       title="Coordenadas pelo método Não-Métrico")

# item c

EM.C <- cmdscale(D, k=2, eig=TRUE)

df2 <- data.frame(EM.C$points)

ggplot(df2, aes(x=X1,y=X2)) + 
  geom_point() + 
  scale_x_continuous(limits = c(-2,2)) +
  geom_text(aes(label=as.character(seq(1:9)), hjust=0.5, vjust=-0.4)) +
  labs(x="Coordenada 1",
       y="Coordenada 2",
       title="Coordenadas pelo método Clássico")
