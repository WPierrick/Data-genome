# Analyse données d'expression


# TD 1 20/09/2015

# Exercice 1

x<-c(-4:4) # Matrice de -4 à 4 ou x<= -4 : 4 ou fonction seq :
x<- seq (-4,4, 1) # intervalle des valeurs puis "pas"
which(x>0) # Positions des valeurs supérieures à 0. Renvoie les indices des valeurs vérifiant la condition
(x>=0)*x # Multiplication par FALSE (0) ou TRUE (1) de x. Produit de x par l'indication de x positif ou nul (met à 0 les valeurs ne vérifiant pas la condition)
sum(x>0) # Somme comptant le nombre de valeurs vérifiant la condition (TRUE)
mean(x>2) # Moyenne
x[which(x>0)] # Crochets : indices dans la matrice (positions)
sqrt(x[which(x>0)]) # Racine carrée des positifs
any(x>4) # Regarde si une valeur est dans le data set. Effectue un test logique qui prend la valeur si une des valeur vérifie la condition


# Exercice 2

# 1)
#a)
runif(1,0,1) # On génère des valeurs selon une loi uniforme continue entre 0 et 1 : nb de valeurs, bornes inf, borne sup
#b)
sample(0:1, 1) # Tirage au hasard entre l'intervalle, 1 = nombre de tirage. Tire seulement des variables DISCRETES
# ou
rbinom(1,1,0.5) # Bernouille
# ou
runif(1,0,1) > 0.5 # Loi uniforme sur 1 lancer
#c)
sum (sample(c(-1,1), 1000, replace=TRUE)) # Somme des tirages avec perte d'un euro ou gain sur 1000 lancers
# ou
sum( rbinom(1000, 1, 0.5)*2-1) # Binomiale

# 2)
x<- sample(1:6, 1000) # x le résultat du lancer sur 1000 lancers
y <- 3*(x>4) +1*(x==4)-2.5*(x<4) # y le gain 

# ou on calcule l'espérance E(Y)=-2.5/2 + 1/6 + (3x2)/6 = -0.08 = -80 € sur 1000 lancers
x <- (-2.5,-2.5, -2.5, 1, 3, 3)
sum(sample(x,1000, replace=TRUE))

# Exercice 3

# 1)
extra <- c(0.7, -1.6, -0.2, -1.2, -0.1, 3.4, 3.7, 0.8, 0, 2, 1.9, 0.8,1.1,0.1,-0.1,4.4, 5.5,1.4,4.6,3.4)
group <- c(rep(1,10),rep(2,10)) # Rep = répéter un chiffre n fois autre possibilité repn (1,2, each = 10)
matrix <- matrix(c(extra, group), ncol = 2) # On met les deux vecteurs dans une matrice
?stem
stem(extra)

# 2)
plot(sort(extra),1:length(extra)/length(extra), type = "s") # Fonction de répartition empirique :  représenter les probablitlités d'avoir des des données inférieures à la moyenne
# On peut aussi utiliser abline pour faire une courbe
plot(ecdf(extra))# fonction de répartition empirique

# 3)
summary (extra [group=1]) # Pour avoir tous les infos (médiane, variance etc,) pour chaque groupe
summary (extra [group=2])
var (extra[group == 1])
var (extra[group == 2])

# Exercice 4

# 1)
# Xsp ~ N(gaussienne)(esperance 1000, écart type 100²)
# Xsn ~ N(400, 150²)
pnorm (700,1000,100) # Pour le spot d'un gène exprimé

# 2)
# Xgp = 1/4 Somme I=1 à 4 Xsp ~ N (1000, 50²)
pnorm(700,1000,50) # Pareil que en 1) mais avec un écart type de 50 (variance 100² / 4)
# P(Xgp <=700)


# 3) 
# On cherche t tel que : 
# P(Xgp<= t) = P (Xgn<=t)
# P((Xgp-1000))/50 <= (t-1000)/50) = P((Xgn-400)/75)>= (t-400)/75)
# On centre C((t-1000)/50) = 1-C((y-400)/75)
# (t-1000)/50 = - (t-400)/75
# 125 t = 95000 / 125
# t = 760
# On vient de déterminer une valeur t telle que :
# La probabliité d'être plus grand que t sous la condition non exprimée est
# La même que la probilité d'être plus petit que t pour les gènes exprimés
# On peut avoir une règle de décision :
# Faux négatifs, on dit qu'il y a pas d'expression alors qu'en réalité ils sont exprimés

# 4)
# Donc probabilité d'être P(Xgp <= 760) :
pnorm(760,1000,50)

# 5)
# Même réponse qu'au dessus, on la même proba d'avoir un faux positif qu'un faux négatif (intersection des deux gaussiennes)

# 6)
# On fait une loi normale
xp<-rnorm(1000,1000,100)
xn<-rnorm(1000,400,150)
x<-c(xp,xn)

hist(x)
hist(x, nclass= 50)
# Breaks on définit les points qui définissent les classes ou identique nclasses, on fait des tailles de classes équiprobables

# Pour superposer la densité : pour chacune des valeurs observées
# On a un mélange de deux gaussiennes, on fait la somme des densités pour chacune des gaussiennes
# ex, densité de 600 : f(600) = 0.5*1/(Sqr(2pi*sigma²))e^(-(x-400/2pi)²)
y<-0.5*dnorm(x,400,150)+0.5*dnorm(x,1000,100) # Densité : commande dnorm
hist(x, nclass= 50, probability = TRUE, xlim=c(-100,1300),ylim = c(0,0.0025))
par(new=TRUE) # Avoir 2 graphes sur le même plot
plot(x[order(x)],y[order(x)],type="l",col="red",lwd=3)


# TD 2 05/10/2015

# Bioconductor : page avec pleins d'outils pour l'analyse de données génomiques
# Pour installer un package, installer la fonction d'installation (ex Bioclite) puis ligne de code sous R

source("http://bioconductor.org/biocLite.R")
source("http://bioconductor.org/marray.R")
biocLite("marray")
library(marray)
help(swirl) # On a les données de l'expérience dans help(swirl)
datadir<-system.file("swirldata",package="marray")
# En regardant le fichier que datadir a généré. 

# 1)
swirltargets<- readTargets(file.path(datadir, "SwirlSample.txt"))# On récupère dans R les données dans le fichier généré
swirltargets # Descriptif de l'expérience, qu'est ce qui est hybridé, comment, de quelle manière

# 2)
?read.maimages # Liste Red/Green des microarrays
RG<-read.maimages((maInfo=swirltargets)$Names, source="spot", path = datadir) # On récupères les intensités
# Les données sur les intensités
show(RG) # Donne des infos sur les données
summary(RG$R) # Statistiques usuelles

# 3)
# Identifier les transcrits