# TD 3 : Données d'expression
# 19/10/2015
# Dalmasso

# 1) Importation des données
# a)
setwd("/home/etudiant/Bureau")
getwd()
wang05 <- read.table("wang05.txt", h=T)
annot <- read.table("annot.txt", h=T)
# Données déjà normalisées et prétraitées

# b)
dim(wang05) # 22277 lignes (les transcrits) et 286 colonnes (les individus)
summary(wang05)
summary(annot)
str(wang05)
head(wang05)
head(annot)

# c)
# Histogramme des ages
hist(annot$age)

plot(annot$e.rfs~annot$age)

# Nombre de rechute
sum(annot$e.rfs)
hist(unlist(wang05))
which(row.names(wang05)=="213627_at")

# 2) Etude préliminaire
# a)
summary(annot[,"age"]) # Age médian 52 ans, plus agés etc...
hist(annot[,"age"])
hist(annot[,"age"], nclass = 50) # On augmente notre nombre de classe


# b)
table(annot[,"e.rfs"]) # Rechute : variable relapse. 179 avec rechute, 107 sans rechute. On a 2 groupes différents

# c)
hist(wang05[,1]) # Distribution du niveau d'expression des gènes pour le patient 1
hist(wang05[,123], breaks = 100) # Niveau pour le patient 123, avec plus de plots (sur nos 20000 genes)

hist(as.numeric(wang05[3,])) # Données normales, bien pour la normalité par la suite
# Pour savoir à quelle ligne se situe la 13627... on utilise la fonction wich
boxplot(as.numeric(wang05["213627_at",])~ as.factor(annot$e.rfs)) # Boxplot de 

# Test de Student
# On se pose la question si il y a une différence entre les individus avec ou sans rechute pour le transcrit particulier
# On fait un test de welsch (pas d'homoscédasticité des variances)
fcttemp<- function(vect){ # Fonction pour test de student par ligne à l'ensemble de la matrice
  vect <- as.numeric(vect)
  res <-t.test(vect~as.factor(annot$e.rfs))$p.value
  return(res)
}
resp <- apply(wang05,1,fcttemp)
# Restemp on stocke dans un vecteur l'ensemble des valeurs
hist(resp) # Histogramme des pvalues
# Si on veut les 10 pvalues les plus significatives
names(sort(resp))[1:10]
# Comment savoir combien de pvalues plusn significatives prendre ? On fait un test multiple (avec multtest)
library(multtest)
install.packages("multtest")
install.packages("biocLite")
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("multtest")
# On ajuste les pvalues a cause du test multiple (5% * le nbombre de transcrits ca fait beaucoup de faux positifs)
# Fonction multtest qui transforme des pvalues brutes en pvalues ajustées. Méthode BH Benjamini et Hochberg
adjp <- mt.rawp2adjp(resp, proc="BH")
adjp$adjp[1:10,] 
sum(adjp$adjp[,2]<= 0.05) # Au niveau 5%, BH donne 117 genes différentiellement exprimés entre les patients avec ou sans rechute (117 pvalues inférieures à 5%)

# 3) Analyse différentielle
