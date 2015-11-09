# Analyse données d'expression

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

# a)
swirltargets<- readTargets(file.path(datadir, "SwirlSample.txt"))# On récupère dans R les données dans le fichier généré
swirltargets # Descriptif de l'expérience, qu'est ce qui est hybridé, comment, de quelle manière

# b)
?read.maimages # Liste Red/Green des microarrays
RG<-read.maimages((maInfo=swirltargets)$Names, source="spot", path = datadir) # On récupères les intensités
# Les données sur les intensités
show(RG) # Donne des infos sur les données. Dans RG on a toutes les infos sur les intensités et leur localisaton sur la puce
summary(RG$R) # Statistiques usuelles

# c)
# Identifier les transcrits, lecture de l'information sur les sondes
RG$genes<-readGAL(galfile="fish.gal",path=datadir)
RG$genes[1:5,]
names(RG$genes)
RG$printer<-getLayout(RG$genes)
RG$printer # Montre comment est organisé la grille de la lame

# Identifier les contrôles
spottypes<-readSpotTypes("SpotTypes.csv", path=".") # On peut identifier les spots qui identifient les contrôles
RG$genes$Status<-controlStatus(spottypes, RG$genes)

# 2) Contrôle qualité
for (i in 1:4) {
  X11()
  par(mfrow=c(2,2))
  imageplot(RG$Rb[,i], layout = RG$printer, low="white", high = "red")
  imageplot(RG$Gb[,i], layout = RG$printer, low="white", high = "green")
  imageplot(RG$G[,i], layout = RG$printer, low="white", high = "green")
  imageplot(RG$R[,i], layout = RG$printer, low="white", high = "red")
}

# On regarde les logs des deux. Base 2 pour binaire (mais peu important)
imageplot(log(RG$R[,1]/RG$G[,1],base = 2), layout = RG$printer, low="green", high = "red")
imageplot(log(RG$R[,2]/RG$G[,2],base = 2), layout = RG$printer, low="green", high = "red")
imageplot(log(RG$R[,3]/RG$G[,3],base = 2), layout = RG$printer, low="green", high = "red")
imageplot(log(RG$R[,4]/RG$G[,4],base = 2), layout = RG$printer, low="green", high = "red")
# En rouge on demande les petites valeurs, en vert les grandes

?imageplot
?par
# On a les taches au milieu : bruit de fond. On peut faire le bruit de fond pour chaque lame

# Boxplot
# Faux boxplot(RG$Rb[,1],RG$genes[,1], main="Variabilité intrabloc du bruit de fond rouge dela puce 1", xlab="Bloc", ylab="Intensité")
boxplot(RG$Gb[,1],RG$genes[,1], main="Variabilité intrabloc du bruit de fond rouge dela puce 1", xlab="Bloc", ylab="Intensité")
boxplot(RG$Gb[,1]~RG$genes[,1], main="Variabilité intrabloc du bruit de fond rouge dela puce 1", xlab="Bloc", ylab="Intensité")
boxplot(RG$Gb[,1]~RG$genes[,1], main="Variabilité intrabloc du bruit de fond rouge dela puce 1", xlab="Bloc", ylab="Intensité")
boxplot(RG$Gb[,1]~RG$genes[,1], main="Variabilité intrabloc du bruit de fond rouge dela puce 1", xlab="Bloc", ylab="Intensité")
boxplot(RG$Gb[,1]~RG$genes[,1], main="Variabilité intrabloc du bruit de fond rouge dela puce 1", xlab="Bloc", ylab="Intensité")
boxplot(RG$Gb[,1], RG$Gb[,2], RG$Gb[,3], RG$Gb[,4], main="Variabilité interplaque du bruit de fond vert", xlab="Bloc", ylab="Intensité", outline = FALSE)
boxplot(RG$R[,1], RG$R[,2], RG$R[,3], RG$R[,4], main="Variabilité interplaque de l'expression en rouge", xlab="Bloc", ylab="Intensité", outline = FALSE)
boxplot(RG$G[,1], RG$G[,2], RG$G[,3], RG$G[,4], main="Variabilité interplaque de l'expression en vert", xlab="Bloc", ylab="Intensité", outline = FALSE)

# 3) Normalisation
# Principe de la normalisation de Loess sur des lames bicolores : corriger un biais qui soit non proportionnel et non constant sur l'intensité mesuré
# Se pase sur des MAPlot : des log ration (log Rouge - Log Vert). Si pas de biais systématique on a un nuage de point en ligne horizontale
# Souvent biais on a un biais entre Rouge et vert qui sera différent selon les intensités faibles ou fortes
# Sur le MAPlot, on a le (log R + log V) / 2. On cherche à estimer notre biais (ex sur l'incorporation des 2 fluos)
# La normalisation de Loess estime le biais et le soustrait aux intensités (logs ratios).
#Peut être faite globale ou locale (pour blocs sur la puce)
# On normalise nos données (avec ou sans bruit de fond etc...)
# Régression par bloc printiploess
# Régression globale : loess

#plot(normalizeBetweenArrays(RG$G) # Régression globale : loess
# normalizeWithinArrays() # Régression locale par bloc : printiploess

# Sans normalisation
plotMA(RG, array = 1) # On regarde le biais lié à l'effet fluo par lame
plotMA(RG, array = 2)
plotMA(RG, array = 3)
plotMA(RG, array = 4)

# Loess locale
MA<-normalizeWithinArrays(RG) 
plotMA(MA)

MA<-normalizeWithinArrays(RG, method = "none") # Par défaut ou soustrait le bruit de fond. Il vaut l'enlever car ça augmente la variablilité du signal
plotMA(MA)

# Normalisation en soustrayant le bruit de fond
boxplot(MA$M~col(MA$M),names= swirltargets$Names)

MA.2<-normalizeWithinArrays(RG, method = "none") # Par défaut ou soustrait le bruit de fond. Il vaut l'enlever car ça augmente la variablilité du signal
plotMA(MA.2)
boxplot(MA.2$M~col(MA.2$M),names= swirltargets$Names)

# Loess globale
MABG<-normalizeWithiArrays(RG, layout = RG$printer, methods = none=
?plotMA
for (i in 1:4) {
  X11()
  par(mfrow=c(2,2))
  plotMA(RG$R[,i])
}
# Si on fait une normalisaiton interlame, et qu'on rajoute des lame,s ça pose problème.
# Différentes facons d'appliquer la normalisation quantile des données (sur l'intensité moyenne, log ratio, rouge, vert etc...)
# Normalisation par quantile sur intensité moyenne

# Between Array
MA<-normalizeBetweenArrays(MAbg)
MA <- normalizeBetweenArrays(RG, method = "quantile")

boxplot(MA$A~col(MA$A, names=swirltargets$Names)
boxplot(MA$M~col(MA$M, names=swirltargets$Names) # Différence des logs
plotMA(MA, array = 4)        

# On a une variabilité plus importante sur la 1ère lame par rapport aux autres, donc on préfère une normalisation interlame.
# Quand on fait une correction interlame, le biais fluo :
# On a fait un dye swap, on a inversé rouge et vert.  Les effets fluos vont se compenser deux à deux
# Dans la normalisation quantile, quand on fait des moyennes, les biais s'annulent entre les lames


# Une fois les données normalisées, on peut faire des tests d'hypothèse etc...