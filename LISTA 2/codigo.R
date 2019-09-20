setwd("/home/gui/Área de Trabalho/LISTA 2 MAE 330")

library(ggplot2)
library(factoextra)
library(cluster)

# Exercício 2

# item b

x1 <- c(2,5,1,8,7)
x2 <- c(0,2,4,4,4)
dados <- data.frame(x1,x2)
rownames(dados) <- c("A","B","C","D","E")

d_pad <- dist(scale(dados))
clus <- hclust(d_pad, method="ward.D")

fviz_dend(clus, k = 2,                 # Cut in four groups
          cex = 0.5,                 # label size
          k_colors = c("#2E9FDF", "#FC4E07"),
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_gray()     # Change theme
) + labs(title="Dendograma",y="Distâncias",x="Agrupamentos") 

# Exercício 3

m.cor <- matrix(c(1,0.63,0.51,0.12,0.16,
                  0.63,1,0.57,0.32,0.21,
                  0.51,0.57,1,0.18,0.15,
                  0.12,0.32,0.18,1,0.68,
                  0.16,0.21,0.15,0.68,1),5,5)

d_pad <- dist(scale(m.cor))

clus <- hclust(d_pad, method="complete")
clus$labels <- setNames(c("JP_Morgan","Citibank","Wells_F.","Royal_D.","Exxon_M."), 1:5)

fviz_dend(clus, k = 2,                 # Cut in four groups
          cex = 0.5,                 # label size
          k_colors = c("#2E9FDF", "#FC4E07"),
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_gray()     # Change theme
) + labs(title="Dendograma",y="Distâncias",x="Agrupamentos") 

# Exercício 4

dados <- read.csv("primate.scapulae.txt",sep=" ")

d_pad <- dist(scale(dados[2:8]))

clus_c <- hclust(d_pad, method="complete")
fviz_dend(clus_c, k = 5,                 # Cut in four groups
          cex = 0.5,                 # label size
          k_colors = c("#0000FF", "#00FF00", "#990000", "#FF6600", "#FF33FF"),
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_gray()     # Change theme
) + labs(title="Dendograma: Vizinho mais longe",y="Distâncias",x="Agrupamentos") 

clus_s <- hclust(d_pad, method="single")
fviz_dend(clus_s, k = 5,                 # Cut in four groups
          cex = 0.5,                 # label size
          k_colors = c("#0000FF", "#00FF00", "#990000", "#FF6600", "#FF33FF"),
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_gray()     # Change theme
) + labs(title="Dendograma: Vizinho mais próximo",y="Distâncias",x="Agrupamentos") 

clus_a <- hclust(d_pad, method="average")
fviz_dend(clus_a, k = 5,                 # Cut in four groups
          cex = 0.5,                 # label size
          k_colors = c("#0000FF", "#00FF00", "#990000", "#FF6600", "#FF33FF"),
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_gray()     # Change theme
) + labs(title="Dendograma: Média das Distâncias",y="Distâncias",x="Agrupamentos") 

clus_d <- diana(d_pad, diss=TRUE)
fviz_dend(clus_d, k = 5,                 # Cut in four groups
          cex = 0.5,                 # label size
          k_colors = c("#0000FF", "#00FF00", "#990000", "#FF6600", "#FF33FF"),
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_gray()     # Change theme
) + labs(title="Dendograma: DIANA",y="Distâncias",x="Agrupamentos") 

previsoes_c <- cutree(clus_c,5)
previsoes_s <- cutree(clus_s,5)
previsoes_a <- cutree(clus_a,5)
previsoes_d <- cutree(clus_d,5)
taxa_c <- 0
taxa_s <- 0
taxa_a <- 0
taxa_d <- 0
for(i in 1:105){
  if(previsoes_c[i]==dados$classdigit[i]){
    taxa_c <- taxa_c+1
  }
  if(previsoes_s[i]==dados$classdigit[i]){
    taxa_s <- taxa_s+1
  }
  if(previsoes_a[i]==dados$classdigit[i]){
    taxa_a <- taxa_a+1
  }
  if(previsoes_d[i]==dados$classdigit[i]){
    taxa_d <- taxa_d+1
  }
}
clusplot(dados[,2:8],previsoes_c,main="Vizinho mais longe")
round(taxa_c/length(dados[,1]),3)*100
clusplot(dados[,2:8],previsoes_s,main="Vizinho mais próximo")
round(taxa_s/length(dados[,1]),3)*100
clusplot(dados[,2:8],previsoes_a,main="Média das Distâncias")
round(taxa_a/length(dados[,1]),3)*100
clusplot(dados[,2:8],previsoes_d,main="DIANA")
round(taxa_d/length(dados[,1]),3)*100
clusplot(dados[,2:8],dados$classdigit,main="Classificação correta")

# Exercicio 5

S <- matrix(c(5,2,2,2),2,2)

# item a

PCA <- princomp(S, cor=F)
summary(PCA, loadings=TRUE)

# item b
cor <- c(S[2]/(sqrt(S[1])*sqrt(S[4])))

C <- matrix(c(1,cor,cor,1),2,2)

PCA_c <- princomp(C)

summary(PCA_c, loadings=TRUE)


# Exercício 6

S <- matrix(c(7476.45,303.62,303.62,26.19),2,2)

# item a

PCA <- princomp(S, cor=F)
summary(PCA, loadings=TRUE)

# item b

autovalores <- eigen(S)$values
autovalores[1]/sum(autovalores)

corr1 <- sqrt(eigen(S)$values[1])*eigen(S)$vectors[1,1]/sqrt(S[1])
corr2 <- sqrt(eigen(S)$values[2])*eigen(S)$vectors[1,2]/sqrt(S[4])
c(corr1,corr2)

# item c

cor <- c(S[2]/(sqrt(S[1])*sqrt(S[4])))

R <- matrix(c(1,cor,cor,1),2,2)

PCA_R <- princomp(R, cor=F)
summary(PCA_R, loadings=TRUE)

autovalores_R <- eigen(R)$values
autovalores_R[1]/sum(autovalores_R)

corr1_R <- sqrt(eigen(R)$values[1])*eigen(R)$vectors[1,1]
corr2_R <- sqrt(eigen(R)$values[2])*eigen(R)$vectors[1,2]
c(corr1_R,corr2_R)


# Exercício 7

# item a

dat <- read.table("T8-6.DAT")
names(dat) <- c("Pais","100m","200m","400m","800m","1500m","5000m","10000m","Maratona")

cor(dat[2:9]) # matriz de covariancias
knitr::kable(caption="matrix de cor", cov(dat[2:9]))

auto_v <- eigen(cor(dat[2:9]))

round(auto_v$values,3)
knitr::kable(caption="Autovalores", round(t(auto_v$values),3),col.names = c( "$\\lambda_1$","$\\lambda_2$","$\\lambda_3$","$\\lambda_4$","$\\lambda_5$","$\\lambda_6$","$\\lambda_7$","$\\lambda_8$"))

round(auto_v$vectors,3)
knitr::kable(caption="Autovetores", round(auto_v$vectors,3),col.names = c("$e_1$","$e_2$","$e_3$","$e_4$","$e_5$","$e_6$","$e_7$","$e_8$"))

# item b

dat_padr <- scale(dat[2:9])
dat_padr <- data.frame(dat_padr)

auto_v_p <- eigen(cor(dat_padr))

duas_pc <- auto_v_p$vectors[,1:2]

# Correlações entre variáveis componentes principais
corr_1comp <- sqrt(eigen(cor(dat_padr))$values[1])*eigen(cor(dat_padr))$vectors[,1]
corr_2comp <- sqrt(eigen(cor(dat_padr))$values[2])*eigen(cor(dat_padr))$vectors[,2]

acum <- vector()
for(i in 1:8){
  acum[i] <- sum(auto_v_p$values[1:i])/sum(auto_v_p$values)
}

acum

# Exercício 8

dado <- read.csv("oleo.csv")
rownames(dado) <- dado$Nation

cor(dado[-1])

PCA <- princomp(dado[-1], cor=T)
summary(PCA, loadings=TRUE)

screeplot(PCA,
         type="lines",
         pch=16,
         col="turquoise4",
         main=" Scree Plot",
         ylim=c(0,2))
abline(h=1,lty=3)

# BIPLOT

biplot(PCA, choices=1:2,main="Biplot")
       
