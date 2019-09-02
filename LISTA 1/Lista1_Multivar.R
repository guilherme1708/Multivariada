setwd('/home/gui/Área de Trabalho/LISTA 1 MAE 330')
library(car)
library(ICS)
library(heplots)
library(multcomp)
library(HoRM)
library(ggplot2)


#1)
tanques <- c(6.2, 432, 12.7, 431, 7,   522, 8.3,  600,
  4.8, 405, 11.3, 426, 4.4, 513, 7.1,  513,
  3,   324, 9.3,  438, 3.8, 507, 4.7,  539,
  5.6, 310, 9.5,  312, 5.0, 410, 10,   456,
  7.4, 326, 11.7, 326, 5.5, 350, 8.5,  504,
  4.8, 375, 15.3, 447, 3.2, 547, 12.4, 548)

tanquey1 <- vector()
tanquey2 <- vector()

for(i in 1:(length(tanques)/2)){
  tanquey1 <- c(tanquey1,tanques[2*i-1])
  tanquey2 <- c(tanquey2,tanques[2*i])
}

reagentes <- rep(c(1,2,3,4),6)

tanq <- data.frame(tanquey1,tanquey2,reagentes)

summary(tanquey1)
summary(tanquey2)

med_tanq1 <- tapply(tanquey1,reagentes,mean)
med_tanq2 <- tapply(tanquey2,reagentes,mean)

tapply(tanquey1,reagentes,sd)
tapply(tanquey2,reagentes,sd)

Y <- cbind(tanquey1,tanquey2)

fit.lm <- lm(Y~factor(reagentes), data=tanq)
fit.manova <- Manova(fit.lm)
summary(fit.manova, multivariate=TRUE)

mnv <- manova(Y~reagentes,data=tanq)

# comparação de matrizes de covariância
boxM(tanq[1:2],reagentes)

mvnorm.skew.test(mnv$res) # teste de normalidade multivariada, baseada na Assimetria

anova1 <- aov(tanquey1~factor(reagentes),data=tanq)
anova2 <- aov(tanquey2~factor(reagentes),data=tanq)

boxplot(tanquey2~reagentes)
boxplot(tanquey1~reagentes)

# diferença dos tratamentos de medida de clorofila

vi <- qt(1-0.05/12,20)*sqrt(3*diag(SSCP.fn(mnv)$SSCPE)[1]/20) # variacao do intervalo
dif12 <- med_tanq1[1] - med_tanq1[2]
dif13 <- med_tanq1[1] - med_tanq1[3]
dif14 <- med_tanq1[1] - med_tanq1[4]
dif23 <- med_tanq1[2] - med_tanq1[3]
dif24 <- med_tanq1[2] - med_tanq1[4]
dif34 <- med_tanq1[3] - med_tanq1[4]
c(dif12-vi,dif12+vi,use.names=F) # T1-T2
c(dif13-vi,dif13+vi,use.names=F) # T1-T3
c(dif14-vi,dif14+vi,use.names=F) # T1-T4
c(dif23-vi,dif23+vi,use.names=F) # T2-T3
c(dif24-vi,dif24+vi,use.names=F) # T2-T4

# diferença dos tratamentos de oxigênio dissolvido na água

vi <- qt(1-0.05/12,20)*sqrt(3*diag(SSCP.fn(mnv)$SSCPE)[2]/20) # variacao do intervalo
dif12 <- med_tanq2[1] - med_tanq2[2]
dif13 <- med_tanq2[1] - med_tanq2[3]
dif14 <- med_tanq2[1] - med_tanq2[4]
dif23 <- med_tanq2[2] - med_tanq2[3]
dif24 <- med_tanq2[2] - med_tanq2[4]
dif34 <- med_tanq2[3] - med_tanq2[4]
c(dif12-vi,dif12+vi,use.names=F) # T1-T2
c(dif13-vi,dif13+vi,use.names=F) # T1-T3
c(dif14-vi,dif14+vi,use.names=F) # T1-T4
c(dif23-vi,dif23+vi,use.names=F) # T2-T3
c(dif24-vi,dif24+vi,use.names=F) # T2-T4

#2
#a
dados2 <- t(matrix( c( 195.3,  153.1,  51.4,
                       194.3,  167.7,  53.7,
                       189.7,  139.5,  55.5,
                       180.4,  121.1,  44.4,
                       203.0,  156.8,  49.8,
                       195.9,  166.0,  45.8,
                       202.7,  166.1,  60.4,
                       197.6,  161.8,  54.1,
                       193.5,  164.5,  57.8,
                       187.0,  165.1,  58.6,
                       201.5,  166.8,  65.0,
                       200.0,  173.8,  67.2),3,12))
fat1 <- factor(c(1,1,2,2,1,1,2,2,1,1,2,2))
fat2 <- factor(c(5,5,5,5,6,6,6,6,8,8,8,8))
dat2 <- data.frame(trat=factor(rep(1:6,each=2)),dados2,fat1,fat2)
attach(dat2)

summary(X1)
summary(X2)
summary(X3)
tapply(X1,fat1,mean)
tapply(X2,fat1,mean)
tapply(X3,fat1,mean)
tapply(X1,fat2,mean)
tapply(X2,fat2,mean)
tapply(X3,fat2,mean)

tapply(X1,fat1,sd)
tapply(X2,fat1,sd)
tapply(X3,fat1,sd)
tapply(X1,fat2,sd)
tapply(X2,fat2,sd)
tapply(X3,fat2,sd)

interaction.plot(fat1, fat2,dados2[,1], main="Efeito de Interação dos Fatores A e B em Y")
interaction.plot(fat1, fat2,dados2[,2], main="Efeito de Interação dos Fatores A e B em Y")
interaction.plot(fat1, fat2,dados2[,3], main="Efeito de Interação dos Fatores A e B em Y")

Y <- cbind(X1,X2,X3)

fit.lm <- lm(Y~fat1+fat2+fat1*fat2, data=dat2)
fit.manova <- Manova(fit.lm)
summary(fit.manova, multivariate=TRUE)

mnv2 <- manova(cbind(X1,X2,X3)~fat1+fat2+fat1*fat2,data=dat2)
manv2 <- summary(mnv2)
par(mfrow=c(1,2))
boxplot(X1~fat1)
boxplot(X1~fat2)
boxplot(X2~fat1)
boxplot(X2~fat2)
boxplot(X3~fat1)
boxplot(X3~fat2)

#b

boxM(Y, as.factor(trat))

mvnorm.skew.test(mnv2$res) # teste de normalidade multivariada, baseada na Assimetria



par(mfrow=c(2,2))
plot(mod1$fit, mod1$res, ylab="Resíduos", xlab="Preditos", main="Resíduos x Preditos")
mnv2res <- mod1$res # extraindo residuos
s21 <- anov1$"Mean Sq"[4]
respad <- (mnv2res/sqrt(s21)) # residuos padronizados
hist(mod1$res, main="Histograma", xlab="Resíduos") #4
qqnorm(mnv2res) #5 normalidade
qqline(mnv2res)
plot(mnv2res,main="Indepêndencia",xlab="Índice", ylab="Resíduos Pad") #6 independencia
shapiro.test(mnv2res)
bartlett.test(X1~trat)

par(mfrow=c(2,2))
plot(mod2$fit, mod2$res, ylab="Resíduos", xlab="Preditos", main="Resíduos x Preditos")
mod2res <- mod1$res # extraindo residuos
s21 <- anov1$"Mean Sq"[4]
respad <- (mod2res/sqrt(s21)) # residuos padronizados
hist(mod1$res, main="Histograma", xlab="Resíduos") #4
qqnorm(mod2res) #5 normalidade
qqline(mod2res)
plot(mod2res,main="Indepêndencia",xlab="Índice", ylab="Resíduos Pad") #6 independencia
shapiro.test(mod2res)
bartlett.test(X2~trat)


par(mfrow=c(2,2))
plot(mod3$fit, mod3$res, ylab="Resíduos", xlab="Preditos", main="Resíduos x Preditos")
mod3res <- mod3$res # extraindo residuos
s21 <- anov1$"Mean Sq"[4]
respad <- (mod3res/sqrt(s21)) # residuos padronizados
hist(mod3$res, main="Histograma", xlab="Resíduos") #4
qqnorm(mod3res) #5 normalidade
qqline(mod3res)
plot(mod3res,main="Indepêndencia",xlab="Índice", ylab="Resíduos Pad") #6 independencia
shapiro.test(mod3res)
bartlett.test(X3~trat)

#c
mod0 <- aov(X1 ~ trat)
mod1 <- aov(X1~fat1+fat2+fat1*fat2,data=dat2)
mod2 <- aov(X2~fat1+fat2+fat1*fat2,data=dat2)
mod3 <- aov(X3~fat1+fat2+fat1*fat2,data=dat2)
anov1 <- anova(mod1)
anov2 <- anova(mod2)
anov3 <- anova(mod3)

#3

dados <- read.table("T1-5.dat", quote="\"", comment.char="")
names(dados) <- c("vento","radiacao","CO","NO","NO2","O3","HC")
attach(dados)
summary(dados)
pairs(dados)

#a
ml1 <- lm(NO2~vento + radiacao, data = dados)
summary(ml1)
z <- data.frame(vento = 10, radiacao=80)
p <- predict(ml1, z,interval =  "prediction")
knitr::kable(caption = "Intervalo de Predicao", p)

#b
mlm1 <- lm(cbind(NO2, O3) ~ vento + radiacao, data = dados)
summary(mlm1)
z <- data.frame(vento = 10, radiacao=80)
p1 <- predict(mlm1, z,interval =  "confidence")
p1
predictionEllipse(mod = mlm1, newdata = z)

# Funções

Cov.Mtest=function(x,ina,a=0.05) {
  ## x is the data set
  ## ina is a numeric vector indicating the groups of the data set
  ## a is the significance level, set to 0.05 by default
  x=as.matrix(x)
  p=ncol(x) ## dimension of the data set
  n=nrow(x) ## total sample size
  k=max(ina) ## number of groups
  nu=rep(0,k)
  ## the sample size of each group will be stored here later
  pame=rep(0,k)
  ## the determinant of each covariance will be stored here
  ## the next "for" function calculates the covariance matrixof each group
  nu=as.vector(table(ina))
  mat=mat1=array(dim=c(p,p,k))
  for (i in 1:k) {
    mat[,,i]=cov(x[ina==i,])
    pame[i]=det(mat[,,i]) ## the detemirnant of each covarian
    mat1[,,i]=(nu[i]-1)*cov(x[ina==i,])
  }
  ## the next 2 lines calculate the pooled covariance matrix
  Sp=apply(mat1,1:2,sum)
  Sp=Sp/(n-k)
  pamela=det(Sp)
  ## determinant of the pooled covariance matrix
  test1=sum((nu-1)*log(pamela/pame))
  gama1=(2*(p^2)+3*p-1)/(6*(p+1)*(k-1))
  gama2=(sum(1/(nu-1))-1/(n-k))
  gama=1-gama1*gama2
  test=gama*test1
  ## this is the M (test statistic)
  df=0.5*p*(p+1)*(k-1)
  pvalue=1-pchisq(test,df) ## p-value of the test statistic
  crit=qchisq(1-a,df)
  ## critical value of the chi-square distribution
  list(M.test=test,degrees=df,critical=crit,p.value=pvalue)
}

predictionEllipse <- function(mod, newdata, level = 0.95, ggplot = TRUE){
  # labels
  lev_lbl <- paste0(level * 100, "%")
  resps <- colnames(mod$coefficients)
  title <- paste(lev_lbl, "confidence ellipse for", resps[1], "and", resps[2])
  
  # prediction
  p <- predict(mod, newdata)
  
  # center of ellipse
  cent <- c(p[1,1],p[1,2])
  
  # shape of ellipse
  Z <- model.matrix(mod)
  Y <- mod$model[[1]]
  n <- nrow(Y)
  m <- ncol(Y)
  r <- ncol(Z) - 1
  S <- crossprod(resid(mod))/(n-r-1)
  
  # radius of circle generating the ellipse
  tt <- terms(mod)
  Terms <- delete.response(tt)
  mf <- model.frame(Terms, newdata, na.action = na.pass, 
                    xlev = mod$xlevels)
  z0 <- model.matrix(Terms, mf, contrasts.arg = mod$contrasts)
  rad <- sqrt((m*(n-r-1)/(n-r-m))*qf(level,m,n-r-m)*z0%*%solve(t(Z)%*%Z) %*% t(z0))
  
  # generate ellipse using ellipse function in car package
  ell_points <- car::ellipse(center = c(cent), shape = S, radius = c(rad), draw = FALSE)
  
  # ggplot2 plot
  if(ggplot){
    require(ggplot2, quietly = TRUE)
    ell_points_df <- as.data.frame(ell_points)
    ggplot(ell_points_df, aes(x, y)) +
      geom_path() +
      geom_point(aes(x = NO2, y = O3 ), data = data.frame(p)) +
      labs(x = resps[1], y = resps[2], 
           title = title)
  } else {
    # base R plot
    plot(ell_points, type = "l", xlab = resps[1], ylab = resps[2], main = title)
    points(x = cent[1], y = cent[2])
  }
}

pottery <- read.csv("pottery.csv")
attach(pottery)
regiao <- as.character(c(rep(1,21),rep(2,14),rep(3,10)))

pottery <- cbind(pottery,regiao)

boxplot(Al2O3~regiao)
boxplot(Fe2O3~regiao)
boxplot(MgO~regiao)
boxplot(CaO~regiao)
boxplot(Na2O~regiao)
boxplot(K2O~regiao)
boxplot(TiO2~regiao)
boxplot(MnO~regiao)
boxplot(BaO~regiao)

reg1 <- pottery[1:21,]
reg2 <- pottery[22:35,]
reg3 <- pottery[36:45,]

boxplot(reg1$Al2O3,reg1$Fe2O3,reg1$MgO,reg1$CaO,reg1$Na2O,reg1$K2O,reg1$TiO2,reg1$MnO,reg1$BaO,ylim=c(0,20))
boxplot(reg2$Al2O3,reg2$Fe2O3,reg2$MgO,reg2$CaO,reg2$Na2O,reg2$K2O,reg2$TiO2,reg2$MnO,reg2$BaO,ylim=c(0,20))
boxplot(reg3$Al2O3,reg3$Fe2O3,reg3$MgO,reg3$CaO,reg3$Na2O,reg3$K2O,reg3$TiO2,reg3$MnO,reg3$BaO,ylim=c(0,20))

ggplot(pottery, aes(x=Al2O3,y=Fe2O3,colour=regiao))+geom_point()+
  scale_colour_discrete(labels=c("1"="1","2"="2","3"="3"), limits=c("1","2","3"))+ labs(colour="Região")+
  ggtitle(paste("Simple Scatter Plot"))

ppi=300

png("Matriz.png", width=7*ppi, height=6*ppi, res=ppi)
pairs(pottery[,c(2:10)], pch=16, col=regiao, cex=0.6)
dev.off()
