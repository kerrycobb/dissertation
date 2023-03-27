library(ggplot2)
library(plotly)
library(adegenet)
library(dplyr)
library(tidyr)
library(ggmap)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(htmlwidgets)

# Load sample data and set row names to match alignment names
sampleData <- read.csv("~/Desktop/toad-project/toad-data.csv")
sampleData <- sampleData[!duplicated(sampleData["rename"]),]
rownames(sampleData) <- sampleData[,"rename"]

# Load alignment data
infile <- "clust-90-indel-16-samples-52-include-americanus-2.str"
nsamples <- 70
nloci <- 7009

# infile <- "clust-90-indel-16-samples-57-include-americanus-1.str"
# nsamples <- 76 
# nloci <- 5854 

input <- read.structure(infile, n.ind=nsamples, n.loc=nloci, ask=FALSE, 
                        onerowperind=FALSE, row.marknames=0, NA.char="-9", col.pop=0, col.lab=1)
print("Missing data:")
sum(is.na(input$tab))

# Subset samples from sample data that is present in alignment
sampleNames <- unlist(dimnames(input@tab)[1])
sampleData <- sampleData[rownames(sampleData) %in% sampleNames,]

# PCA
X <- tab(input, freq=TRUE, NA.method="mean")
pca <- dudi.pca(X, scale=FALSE, scannf=FALSE, nf=3)
barplot(pca$eig[1:50], main="PCA eigenvalues", col=heat.colors(50))
ggsave("americanus-eigenvalues.svg", device="svg")

ggplot(pca$li, aes(x=Axis1, y=Axis2)) +
  geom_point(size=1)
ggsave("americanus-pca.svg", device="svg")
p <- plot_ly(x=pca$li$Axis1, y=pca$li$Axis2, z=pca$li$Axis3, type="scatter3d", mode="markers")
p
saveWidget(p, "americanus-pca.html")

# DAPC
grp <- find.clusters(X, max.n.clust=10)
dapc1 <- dapc(X, grp$grp)
scatter(dapc1) 

sampleData["group"] <- NA
for (i in 1: length(dapc1$grp)){
  sampleData[names(dapc1$grp[i]), "group"] <- as.numeric(dapc1$grp[i]) 
}

# Maps
xlim1 <- c(-96.5, -72)
ylim1 <- c(31.5, 47)
scale <- 110
linewidth <- 1
# countries <- ne_countries(scale=scale, returnclass="sf")
# coast <- ne_coastline(scale=scale, returnclass="sf")
states <- ne_states(country="united states of america", returnclass="sf")


map1 <- ggplot() +
  geom_sf(data=states) +
  # geom_sf(data=coast) +
  geom_point(data=sampleData, mapping=aes(x=longitude, y=latitude, color=as.factor(group))) +
  coord_sf(xlim=xlim1, ylim=ylim1, expand=FALSE)
map1
p <- ggplotly(map1)
p
saveWidget(p, "americanus-dapc-map.html")

xlim2 <- c(-87, -85)
ylim2 <- c(32.5, 34)
ggplot() +
  geom_sf(data=states) +
  geom_point(data=sampleData, mapping=aes(x=longitude, y=latitude, color=as.factor(group))) +
  coord_sf(xlim=xlim2, ylim=ylim2, expand=FALSE)











