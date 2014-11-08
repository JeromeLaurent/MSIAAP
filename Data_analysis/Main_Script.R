

#######################################
########       Projet AD       ########
########       Sujet n°4       ########
########    encoding : UTF8    ########
#######################################


# Définition du répertoire de travail

setwd("~/R/MSIAAP/Data_analysis")

# Lecture des fichiers, regroupés en un fichier commun

dataset <- do.call(cbind, lapply(list.files("~/R/MSIAAP/Data_analysis/Data", 
                full.names = TRUE), read.delim, header = TRUE))

head(dataset['Vss'])

str(dataset)
na <- is.na(dataset[[ncol(dataset)]])
class(na)
as.numeric(na)
