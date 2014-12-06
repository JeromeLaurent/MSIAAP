

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
require(glmnet) # http://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
require(caret) # Feature selection. http://cran.r-project.org/web/packages/caret/vignettes/caret.pdf
# http://topepo.github.io/caret/featureselection.html
require(corrplot)
require(elasticnet)
require(neuralnet)


# Définition du répertoire de travail

setwd("~/R/MSIAAP/Data_analysis")

### Lecture des données ###

# Lecture des fichiers, regroupés en un fichier commun dat

dat <- do.call(cbind, lapply(list.files("~/R/MSIAAP/Data_analysis/Data", 
                   full.names = TRUE), read.delim, header = TRUE))

dat <- dat[,!(names(dat) %in% "Ligand")] # enlève les colonnes "Ligands" qui étaient répétées

# NB : grep() = pattern matching and replacement

dim(dat) # 279 variables : 278 descripteurs moléculaires + 1 variable explicative

# Sélection des variables non numériques
dat[!sapply(dat,is.numeric)]  # Aucune n'est retenue, les variables ont toutes été bien reconnues


str(dat)
na <- is.na(dat[[ncol(dat)]])
as.numeric(na)


dataset <- dat
dataset$Vss <- log(dataset$Vss) 
names(dataset)[names(dataset)=="Vss"] <- "logvss" # renomme Vss en logvss

### Création des sets d'entraînement et de test ###

### Sélection des variables ###

# p-value : mauvais
# AIC & BIC : à éviter
# LASSO & Elastic Net : OUi 

# Centrage-réduction de toutes les variables explicatives
dat.scale<- scale(dataset[1:(ncol(dat)-1)],center=TRUE,scale=TRUE) 

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
res.plsr <- plsr(logvss ~ ., data = dataset, validation = "CV", segments = 2)

res.plsr$validation$PRESS
#plot(1:(ncol(dataset)-1), res.plsr$validation$PRESS)

#plot(1:100, res.plsr$validation$PRESS[,1:100])
min(res.plsr$validation$PRESS)
# Beaucoup trop de variables, le PRESS augmente extremement vite.


# fit model
res.pls <- plsr(logvss ~ ., data=dataset, validation="CV")
# summarize the fit
summary(res.pls)
# make predictions
predictions <- predict(res.pls, dataset, ncomp=6)
# summarize accuracy
rmse <- mean((dataset$logvss - predictions)^2)
print(rmse)
min(res.pls$validation$PRESS)


# Test LASSO

xfactors <- model.matrix(logvss ~ ., data = dataset)#[,-1]

glmmod <- glmnet(xfactors, y = dataset$logvss,alpha=1,family='multinomial')

#plot variable coefficients vs. shrinkage parameter lambda.
#plot(glmmod,xvar="lambda")
#grid()

# Test Elastic Net


x <- as.matrix(dataset[,1:278])
y <- dataset[,279]

# Optimisation des paramètres


## Optimisation du paramètre pour la correction de l'erreur, C

lambda <- 10^((-5):1) # Séquences de puissances de 10 à tester pour lambda
alpha <- seq(from = 0, to = 1, by = 0.1)

enrmoyen = NULL
for (i in lambda)
{
        for (k in alpha)
        {
                        pred.CV <- matrix(NA, nrow = 633, ncol = 30)
                        rmse <- rep(NA, 30)
                        
                        for (j in 1:30)
                        {
                                ID <- sample(1:nrow(dataset), round((nrow(dataset))/2))
                                
                                
                                fit <- glmnet(x[-ID,], y[-ID], family="gaussian", alpha=alpha, lambda= lambda)
                                pred.glm <- predict(fit, x[ID, ])
                                
                                pred.CV[ID, j] <- pred.glm # on a construit le modèle sans ID, donc on prédit ID
                                
                                fit <- glmnet(x[ID,], y[ID], family="gaussian", alpha=alpha, lambda= lambda)
                                pred.glm <- predict(fit, x[-ID, ])
                                
                                pred.CV[- ID, j] <- pred.glm
                                
                                
                                rmse[j] <- mean((y - pred.CV[, j])^2) #pourcentage de prédiction pour chaque modèle
                                
                        }
                        rmsemoy <- c(rmsemoy, mean(rmse)) # Evolution du pourcentage de prédiction selon C
        }                
}



fit <- glmnet(x[-ID,], y[-ID], family="gaussian", alpha=0.5, lambda= 0.001)




# Validation croisée pour obtenir le meilleur lambda
cvfit <- cv.glmnet(x, y)

fit <- glmnet(x, y, family="gaussian", alpha=0.5, lambda= 0.001)

print(fit)

coef(fit,s=0.1)

fit <- glmnet(x, y)

cvfit <- cv.glmnet(x, y)

plot(cvfit)
# summarize the fit
summary(fit)
# make predictions
predictions <- predict(fit, x, type="link")
# summarize accuracy
rmse <- mean((y - predictions)^2)
print(rmse)

### Validation du modèle ###



#### Utilisation du package "caret" ###

set.seed(157) # Permet d'obtenir des résultats reproductibles en fixant l'aléatoire

# stratiﬁed random split of the data

inTrain <- createDataPartition(y = dataset$logvss, p = .75,list = FALSE)
training <- dataset[ inTrain,]
testing <- dataset[-inTrain,]

plsFit <- train(logvss ~ ., data = training,  method = "pls",
                        preProc = c("center", "scale")) ## Center and scale the predictors for the training set and all future samples.

plsPred <- predict(plsFit, newdata = testing)


ridgeFit1 <- train(logvss ~ ., data = training, method = 'ridge', 
                   preProc = c("center", "scale"))
plot(ridgeFit1)

neuralFit <- train(logvss ~ ., data = training, method = 'neuralnet', 
                   preProc = c("center", "scale"))
