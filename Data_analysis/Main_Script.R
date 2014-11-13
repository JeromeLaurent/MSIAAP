

#######################################
########       Projet AD       ########
########       Sujet n°4       ########
########    encoding : UTF8    ########
#######################################

# Packages nécessaires

require(ggplot2)
require(plyr)
require(dplyr)
require(pls) # PLS
require(glmnet) # LASSO & elastic net
require(caret) # Feature selection. http://cran.r-project.org/web/packages/caret/vignettes/caret.pdf
require(corrplot)


# Définition du répertoire de travail

setwd("~/R/MSIAAP/Data_analysis")

### Lecture des données ###

# Lecture des fichiers, regroupés en un fichier commun dataset

dataset <- do.call(cbind, lapply(list.files("~/R/MSIAAP/Data_analysis/Data", 
                   full.names = TRUE), read.delim, header = TRUE))

dataset <- dataset[,!(names(dataset) %in% "Ligand")] # enlève les colonnes "Ligands" qui étaient répétées

# NB : grep() = pattern matching and replacement

dim(dataset) # 279 variables : 278 descripteurs moléculaires + 1 variable explicative

# Sélection des variables non numériques
dataset[!sapply(dataset,is.numeric)]  # Aucune n'est retenue, les variables ont toutes été bien reconnues


str(dataset)
na <- is.na(dataset[[ncol(dataset)]])
as.numeric(na)


### Création des sets d'entraînement et de test ###

### Sélection des variables ###

# p-value : mauvais
# AIC & BIC : à éviter
# LASSO & Elastic Net : OUi 

# Centrage-réduction de toutes les variables explicatives
dat.scale<- scale(dataset[1:(ncol(dataset)-1)],center=TRUE,scale=TRUE) 

colnames(dat.scale) <- (1:(ncol(dataset)-1))
# Nom des colonnes remplacés par des chiffres afin de simplifier la représentation graphique

corMat <- cor(dat.scale)
# Calcul des corrélations

corrplot(corMat, order = "hclust", tl.cex = 0.1)
#visualize the matrix, clustering features by correlation index.
# Il y a des groupes de variables fortement corrélées.

highCor <- findCorrelation(corMat, 0.70)
#Apply correlation filter at 0.70,
#then we remove all the variable correlated with more 0.7.
datFiltered.scale <- dat.scale[,-highCor]
dim(datFiltered.scale)
corMat.Filt <- cor(datFiltered.scale)
corrplot(corMat.Filt, order = "hclust", tl.cex = 0.1)



### Création du modèle ###


# Test PLS
res.plsr <- plsr(Vss ~ ., data = dataset, validation = "CV", segments = 2)

res.plsr$validation$PRESS
#plot(1:(ncol(dataset)-1), res.plsr$validation$PRESS)

#plot(1:100, res.plsr$validation$PRESS[,1:100])
min(res.plsr$validation$PRESS)
# Beaucoup trop de variables, le PRESS augmente extremement vite.

# Test LASSO

xfactors <- model.matrix(Vss ~ ., data = dataset)#[,-1]

glmmod <- glmnet(xfactors, y = dataset$Vss,alpha=1,family='multinomial')

#plot variable coefficients vs. shrinkage parameter lambda.
#plot(glmmod,xvar="lambda")
#grid()



### Validation du modèle ###
