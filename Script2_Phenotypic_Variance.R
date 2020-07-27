setwd("C:/Users/User/Desktop/2020_Corona_Fun/Coding/Rm1_ToadRDA/Code")

library(dplyr)
library(readxl)
library(gdata)
library(nlme)
library(ggplot2)
library(vegan)
library(factoextra)
library(FD)
library(tibble)
library("ggpubr")
#install.packages('Rmisc', dependencies = TRUE)
library("Rmisc")



######
######Ammending fdisp function to report per sample distances
######View(fdisp) to get the original function
######

panda <- function (d, a, tol = 1e-07) 
{
  if (!inherits(d, "dist")) 
    stop("'d' must be a 'dist' object.")
  n <- attr(d, "Size")
  if (is.null(attr(d, "Labels"))) 
    stop("'d' must have labels.", "\n")
  else sn.d <- attr(d, "Labels")
  if (missing(a)) {
    ab.names <- list("Community1", sn.d)
    a <- matrix(1, 1, n, dimnames = ab.names)
  }
  com <- nrow(a)
  if (!is.matrix(a)) 
    stop("'a' must be a matrix.")
  if (ncol(a) != n) 
    stop("Number of columns in 'a' must be equal to the number of objects in 'd'.")
  if (is.null(colnames(a))) 
    stop("'a' must have column names", "\n")
  else sn.a <- colnames(a)
  if (any(sn.d != sn.a)) 
    stop("Species labels in 'd' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
         "\n")
  a[which(is.na(a))] <- 0
  abun.sum <- apply(a, 1, sum)
  if (any(abun.sum == 0)) 
    stop("At least one community has zero-sum abundances (no species).", 
         "\n")
  abun.sum2 <- apply(a, 2, sum)
  if (any(abun.sum2 == 0)) 
    stop("At least one species does not occur in any community (zero total abundance across all communities).", 
         "\n")
  if (any(is.na(d))) 
    stop("NA's in the distance matrix.", "\n")
  A <- matrix(0, ncol = n, nrow = n)
  A[row(A) > col(A)] <- -0.5 * d^2
  A <- A + t(A)
  G <- bicenter.wt(A)
  e <- eigen(G, symmetric = TRUE)
  vectors <- e$vectors
  eig <- e$values
  w0 <- eig[n]/eig[1]
  if (w0 > -tol) 
    r <- sum(eig > (eig[1] * tol))
  else r <- length(eig)
  vectors <- vectors[, 1:r, drop = FALSE] %*% diag(sqrt(abs(eig <- eig[1:r])), 
                                                   r)
  dimnames(vectors) <- list(colnames(a), NULL)
  pos <- eig > 0
  avg.dist.cent <- rep(NA, nrow(a))
  names(avg.dist.cent) <- row.names(a)
  for (i in 1:com) {
    pres <- which(a[i, ] > 0)
    nb.sp <- nrow((unique(vec <- vectors[pres, , drop = F])))
    if (nb.sp >= 2) {
      w <- a[i, pres]
      centroid <- apply(vec, 2, weighted.mean, w = w)
      dist.pos <- sweep(vec[, pos, drop = F], 2, centroid[pos])
      dist.pos <- rowSums(dist.pos^2)
      if (any(!pos)) {
        dist.neg <- sweep(vec[, !pos, drop = F], 2, 
                          centroid[!pos])
        dist.neg <- rowSums(dist.neg^2)
      }
      else dist.neg <- 0
      zij <- sqrt(abs(dist.pos - dist.neg))
      avg.dist.cent[i] <- weighted.mean(zij, w)
    }
    else avg.dist.cent[i] <- 0
  }
  return(list(FDis = avg.dist.cent, zij = zij, eig = eig, vectors = vectors))
}


######
###### Wild caught toads (with Sex)
###### Functional Dispersion ANOVA and plot
######


###### Importing the data
FDis_MASTER_FG_F <- read.csv("FDis_MASTER_FG_F.csv",stringsAsFactors=TRUE,sep=",",row.names=1)
fdisp(gowdis(FDis_MASTER_FG_F))$`FDis`

FDis_MASTER_HI_F <- read.csv("FDis_MASTER_HI_F.csv",stringsAsFactors=TRUE,sep=",",row.names=1)
fdisp(gowdis(FDis_MASTER_HI_F))$`FDis`

FDis_MASTER_QLD_F <- read.csv("FDis_MASTER_QLD_F.csv",stringsAsFactors=TRUE,sep=",",row.names=1)
fdisp(gowdis(FDis_MASTER_QLD_F))$`FDis`

FDis_MASTER_NT_F <- read.csv("FDis_MASTER_NT_F.csv",stringsAsFactors=TRUE,sep=",",row.names=1)
fdisp(gowdis(FDis_MASTER_NT_F))$`FDis`

FDis_MASTER_WA_F <- read.csv("FDis_MASTER_WA_F.csv",stringsAsFactors=TRUE,sep=",",row.names=1)
fdisp(gowdis(FDis_MASTER_WA_F))$`FDis`


FDis_MASTER_FG_M <- read.csv("FDis_MASTER_FG_M.csv",stringsAsFactors=TRUE,sep=",",row.names=1)
fdisp(gowdis(FDis_MASTER_FG_M))$`FDis`

FDis_MASTER_HI_M <- read.csv("FDis_MASTER_HI_M.csv",stringsAsFactors=TRUE,sep=",",row.names=1)
fdisp(gowdis(FDis_MASTER_HI_M))$`FDis`

FDis_MASTER_QLD_M <- read.csv("FDis_MASTER_QLD_M.csv",stringsAsFactors=TRUE,sep=",",row.names=1)
fdisp(gowdis(FDis_MASTER_QLD_M))$`FDis`

FDis_MASTER_NT_M <- read.csv("FDis_MASTER_NT_M.csv",stringsAsFactors=TRUE,sep=",",row.names=1)
fdisp(gowdis(FDis_MASTER_NT_M))$`FDis`

FDis_MASTER_WA_M <- read.csv("FDis_MASTER_WA_M.csv",stringsAsFactors=TRUE,sep=",",row.names=1)
fdisp(gowdis(FDis_MASTER_WA_M))$`FDis`


###### Getting PhenoDisp values
ind12.FG_F <- as.data.frame(panda(gowdis(FDis_MASTER_FG_F[,2:11]))$`zij`)
names(ind12.FG_F)[1]<-paste("distance")
ind12.HI_F <- as.data.frame(panda(gowdis(FDis_MASTER_HI_F[,2:11]))$`zij`)
names(ind12.HI_F)[1]<-paste("distance")
ind12.QLD_F <- as.data.frame(panda(gowdis(FDis_MASTER_QLD_F[,2:11]))$`zij`)
names(ind12.QLD_F)[1]<-paste("distance")
ind12.NT_F <- as.data.frame(panda(gowdis(FDis_MASTER_NT_F[,2:11]))$`zij`)
names(ind12.NT_F)[1]<-paste("distance")
ind12.WA_F <- as.data.frame(panda(gowdis(FDis_MASTER_WA_F[,2:11]))$`zij`)
names(ind12.WA_F)[1]<-paste("distance")


ind12.FG_M <- as.data.frame(panda(gowdis(FDis_MASTER_FG_M[,2:11]))$`zij`)
names(ind12.FG_M)[1]<-paste("distance")
ind12.HI_M <- as.data.frame(panda(gowdis(FDis_MASTER_HI_M[,2:11]))$`zij`)
names(ind12.HI_M)[1]<-paste("distance")
ind12.QLD_M <- as.data.frame(panda(gowdis(FDis_MASTER_QLD_M[,2:11]))$`zij`)
names(ind12.QLD_M)[1]<-paste("distance")
ind12.NT_M <- as.data.frame(panda(gowdis(FDis_MASTER_NT_M[,2:11]))$`zij`)
names(ind12.NT_M)[1]<-paste("distance")
ind12.WA_M <- as.data.frame(panda(gowdis(FDis_MASTER_WA_M[,2:11]))$`zij`)
names(ind12.WA_M)[1]<-paste("distance")


#dist.ALL_F <- rbind(ind12.FG, ind12.HI, ind12.QLD, ind12.NT, ind12.WA)
#dist.ALL_M <- rbind(ind12.FG, ind12.HI, ind12.QLD, ind12.NT, ind12.WA)
#groups.ALL <- cbind(c(rep("1-FG",83),rep("2-HI",499),rep("3-QLD",270), rep("4-NT",428),rep("5-WA",452))) #F
#groups.ALL <- cbind(c(rep("1-FG",214),rep("2-HI",509),rep("3-QLD",235), rep("4-NT",497),rep("5-WA",391))) #M
#groups.ALL <- cbind(c(rep("1-FG-F",83),rep("1-FG-M",214),rep("2-HI-F",499),rep("2-HI-M",509),rep("3-QLD-F",270),rep("3-QLD-M",235),rep("4-NT-F",428),rep("4-NT-M",497),rep("5-WA-F",452),rep("5-WA-M",391))) #F+M
#dataframe.ALL <- data.frame(dist.ALL, groups.ALL)

dist.ALL <- rbind(ind12.FG_F, ind12.FG_M, ind12.HI_F, ind12.HI_M, ind12.QLD_F, ind12.QLD_M, ind12.NT_F, ind12.NT_M, ind12.WA_F, ind12.WA_M)
str(dist.ALL)
groups.ALL.POP <- cbind(c(rep("1-FG",83),rep("1-FG",214),rep("2-HI",499),rep("2-HI",509),rep("3-QLD",270),rep("3-QLD",235),rep("4-NT",428),rep("4-NT",497),rep("5-WA",452),rep("5-WA",391))) 
groups.ALL.SEX <- cbind(c(rep("Female",83),rep("Male",214),rep("Female",499),rep("Male",509),rep("Female",270),rep("Male",235),rep("Female",428),rep("Male",497),rep("Female",452),rep("Male",391))) #F+M

dataframe.ALL <- data.frame(dist.ALL, groups.ALL.POP, groups.ALL.SEX)


fdisp.summary <- summarySE(dataframe.ALL, measurevar="distance", groupvars=c("groups.ALL.POP","groups.ALL.SEX"))
fdisp.summary

res.aov2 <- aov(distance ~ groups.ALL.POP * groups.ALL.SEX, data = dataframe.ALL)
summary(res.aov2)

res.aov2 <- aov(distance ~ groups.ALL, data = dataframe.ALL)
summary(res.aov2)

TukeyHSD(res.aov2)


######
###### Wild caught toads (No Sex)
###### Functional Dispersion ANOVA and plot
######

ind12.FG_F <- as.data.frame(panda(gowdis(FDis_MASTER_FG_F[,2:11]))$`zij`)
names(ind12.FG_F)[1]<-paste("distance")
ind12.HI_F <- as.data.frame(panda(gowdis(FDis_MASTER_HI_F[,2:11]))$`zij`)
names(ind12.HI_F)[1]<-paste("distance")
ind12.QLD_F <- as.data.frame(panda(gowdis(FDis_MASTER_QLD_F[,2:11]))$`zij`)
names(ind12.QLD_F)[1]<-paste("distance")
ind12.NT_F <- as.data.frame(panda(gowdis(FDis_MASTER_NT_F[,2:11]))$`zij`)
names(ind12.NT_F)[1]<-paste("distance")
ind12.WA_F <- as.data.frame(panda(gowdis(FDis_MASTER_WA_F[,2:11]))$`zij`)
names(ind12.WA_F)[1]<-paste("distance")


ind12.FG_M <- as.data.frame(panda(gowdis(FDis_MASTER_FG_M[,2:11]))$`zij`)
names(ind12.FG_M)[1]<-paste("distance")
ind12.HI_M <- as.data.frame(panda(gowdis(FDis_MASTER_HI_M[,2:11]))$`zij`)
names(ind12.HI_M)[1]<-paste("distance")
ind12.QLD_M <- as.data.frame(panda(gowdis(FDis_MASTER_QLD_M[,2:11]))$`zij`)
names(ind12.QLD_M)[1]<-paste("distance")
ind12.NT_M <- as.data.frame(panda(gowdis(FDis_MASTER_NT_M[,2:11]))$`zij`)
names(ind12.NT_M)[1]<-paste("distance")
ind12.WA_M <- as.data.frame(panda(gowdis(FDis_MASTER_WA_M[,2:11]))$`zij`)
names(ind12.WA_M)[1]<-paste("distance")


dist.ALL.NoSex <- rbind(ind12.FG_F, ind12.FG_M, ind12.HI_F, ind12.HI_M, ind12.QLD_F, ind12.QLD_M, ind12.NT_F, ind12.NT_M, ind12.WA_F, ind12.WA_M)
str(dist.ALL.NoSex)
groups.ALL.POP.NoSex <- cbind(c(rep("1-FG",83),rep("1-FG",214),rep("2-HI",499),rep("2-HI",509),rep("3-QLD",270),rep("3-QLD",235),rep("4-NT",428),rep("4-NT",497),rep("5-WA",452),rep("5-WA",391))) 
groups.ALL.SEX.NoSex <- cbind(c(rep("Adult",83),rep("Adult",214),rep("Adult",499),rep("Adult",509),rep("Adult",270),rep("Adult",235),rep("Adult",428),rep("Adult",497),rep("Adult",452),rep("Adult",391))) #F+M

dataframe.ALL.NoSex <- data.frame(dist.ALL.NoSex, groups.ALL.POP.NoSex, groups.ALL.SEX.NoSex)

head(dataframe.ALL.NoSex)


fdisp.summary <- summarySE(dataframe.ALL.NoSex, measurevar="distance", groupvars=c("groups.ALL.POP.NoSex"))
fdisp.summary

res.aov2 <- aov(distance ~ groups.ALL.POP.NoSex, data = dataframe.ALL.NoSex)
summary(res.aov2)

TukeyHSD(res.aov2)



######
###### Common-garden-raised progeny
###### Functional Dispersion ANOVA and plot
######

FDis_QLD_all <- read.csv("FDis_QLD_All.csv",stringsAsFactors=TRUE,sep=",",row.names=1)
fdisp(gowdis(FDis_QLD_all[1:5]))$`FDis` #Morphology

FDis_NT_all <- read.csv("FDis_NT_All.csv",stringsAsFactors=TRUE,sep=",",row.names=1)
fdisp(gowdis(FDis_NT_all[1:5]))$`FDis` #Morphology

FDis_WA_all <- read.csv("FDis_WA_All.csv",stringsAsFactors=TRUE,sep=",",row.names=1)
panda(gowdis(FDis_WA_all[1:5]))$`FDis` #Morphology

ind12.QLD_morph <- as.data.frame(panda(gowdis(FDis_QLD_all[1:5]))$`zij`)
names(ind12.QLD_morph)[1]<-paste("distance")
ind12.NT_morph <- as.data.frame(panda(gowdis(FDis_NT_all[1:5]))$`zij`)
names(ind12.NT_morph)[1]<-paste("distance")
ind12.WA_morph <- as.data.frame(panda(gowdis(FDis_WA_all[1:5]))$`zij`)
names(ind12.WA_morph)[1]<-paste("distance")

dist.ALL.CG <- rbind(ind12.QLD_morph, ind12.NT_morph, ind12.WA_morph)
str(dist.ALL.CG)

groups.ALL.POP.CG <- cbind(c(rep("3-QLD",34), rep("4-NT",19),rep("5-WA",23)))
groups.ALL.TREATMENT.CG <- cbind(c(rep("Juvenile",76)))

dataframe.ALL.CG <- data.frame(dist.ALL.CG, groups.ALL.POP.CG, groups.ALL.TREATMENT.CG)
dataframe.ALL.CG


fdisp.summary <- summarySE(dataframe.ALL, measurevar="distance", groupvars=c("groups.ALL.POP","groups.ALL.TREATMENT"))
fdisp.summary

res.aov2 <- aov(distance ~ groups.ALL.POP.CG, data = dataframe.ALL.CG)
summary(res.aov2)

######
###### Figure 1 plot
###### Combine data sets for the plot
######

library(dplyr)
library(plyr)

str(dataframe.ALL)
str(dataframe.ALL.CG)
str(dataframe.ALL.NoSex)

colnames(dataframe.ALL.CG) <- c("distance","groups.ALL.POP",
                          "groups.ALL.SEX")

colnames(dataframe.ALL.NoSex) <- c("distance","groups.ALL.POP",
                                "groups.ALL.SEX")

str(dataframe.ALL.CG)
str(dataframe.ALL.NoSex)

bind2 <- rbind(dataframe.ALL, dataframe.ALL.CG, dataframe.ALL.NoSex)

str(bind2)

fdisp.summary <- summarySE(bind2, measurevar="distance", groupvars=c("groups.ALL.POP","groups.ALL.SEX"))
fdisp.summary

###### plotting figure 1
Figure1<- ggplot(fdisp.summary, aes(x=groups.ALL.POP, y=distance, colour=groups.ALL.SEX, group=groups.ALL.SEX)) + 
  geom_errorbar(aes(ymin=distance-se, ymax=distance+se), colour="black", width=.1) +
  geom_line(size=1.0) +
  geom_point(size=3, shape=21, fill="white") +
  xlab("Population") +
  ylab("Phenotypic Dispersion") +
  scale_color_manual(values=c("#000000", "#A9A9A9" ,"#FF0000","663300")) +
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(0.95,0.78))  +
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text(colour="black", size=12, 
                                   face="bold")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Figure1

png(filename = paste("Figure1_Variance_Final.png", sep = ""), width = 800, height = 600, units = "px", res = NA)
print(Figure1)
dev.off()


######
###### ADONIS2 Wild caught cane toads
###### Analysis of variance using distance matrices
######

CaneToad_Pheno_SC <- read.csv("Cane_Toad_Master_Sheet_SC_adult_NoNSW_transformed.csv",stringsAsFactors=TRUE,sep=",",row.names=1)
str(CaneToad_Pheno_SC)

dune <- CaneToad_Pheno_SC[,c(9:18)]
str(dune)

dune.env <- CaneToad_Pheno_SC[,c(1:8)]
str(dune.env)

adonis2(dune ~ State, data = CaneToad_Pheno_SC)

#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#adonis2(formula = dune ~ State, data = CaneToad_Pheno_SC)
#Df SumOfSqs      R2      F Pr(>F)    
#State       4   0.6087 0.06012 57.138  0.001 ***
#  Residual 3573   9.5153 0.93988                  
#Total    3577  10.1240 1.00000                  
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



######
###### ADONIS2 common-garden raised progeny
###### Analysis of variance using distance matrices
######


CaneToad_Pheno_SC <- read.csv("CaneToad_Pheno_SC_slim_transformed.csv",stringsAsFactors=TRUE,sep=",",row.names=1)
str(CaneToad_Pheno_SC)

dune <- CaneToad_Pheno_SC[,c(3:9)]
str(dune)

dune.env <- CaneToad_Pheno_SC[,c(1:2)]
str(dune.env)


adonis2(dune ~ Population.Origin, data = CaneToad_Pheno_SC)

#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#adonis2(formula = dune ~ Population.Origin, data = CaneToad_Pheno_SC)
#Df SumOfSqs      R2      F Pr(>F)  
#Population.Origin  2   0.2552 0.07314 2.8802  0.015 *
#  Residual          73   3.2345 0.92686                
#Total             75   3.4897 1.00000                
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1




