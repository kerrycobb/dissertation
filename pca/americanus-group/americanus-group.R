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

infile <- "clust-90-indel-16-samples-166-include-americanus-group-1.str"
nsamples <- 223
nloci <- 2843

input <- read.structure(infile, n.ind=nsamples, n.loc=nloci, ask=FALSE, 
                        onerowperind=FALSE, row.marknames=0, NA.char="-9", 
                        col.pop=0, col.lab=1)

# https://tomjenkins.netlify.app/tutorials/r-popgen-getting-started/
propTyped(input, by="loc")
propTyped(input, by="ind")
which(propTyped(input, by="ind") < 0.75)


print("Missing data:")
sum(is.na(input$tab))

# Subset samples from sample data that is present in alignment
sampleNames <- unlist(dimnames(input@tab)[1])
sampleData <- sampleData[rownames(sampleData) %in% sampleNames,]

# PCA
X <- tab(input, freq=TRUE, NA.method="mean")
pca <- dudi.pca(X, scale=FALSE, scannf=FALSE, nf=3)

barplot(pca$eig[1:50], main="PCA eigenvalues", col=heat.colors(50))
ggsave("americanus-group-eigenvalues.svg", device="svg")

# ggplot PCA plot
split <- strsplit(row.names(pca$li), '_')
df <- data.frame(pca$li, do.call(rbind, split))
ggplot(df, aes(x=Axis1, y=Axis2)) +
  geom_point(aes(color=X2), size=1)
ggsave("americanus-group-pca.svg", device="svg")

# Plotly PCA Plot
p <- plot_ly(x=df$Axis1, y=df$Axis2, z=df$Axis3, type="scatter3d", 
             mode="markers", color=df$X2, text=rownames(df))
p
saveWidget(p, "americanus-group-pca.html")

# dapc
grp <- find.clusters(X, max.n.clust=10)
dapc1 <- dapc(X, grp$grp)
dapc1$grp
scatter(dapc1)

sampleData["group"] <- NA
for (i in 1: length(dapc1$grp)){
  sampleData[names(dapc1$grp[i]), "group"] <- as.numeric(dapc1$grp[i]) 
}

write.csv(sampleData, "americanus-group-dapc-groups.csv", row.names=FALSE)

# Maps
xlim1 <- c(-120, -72)
ylim1 <- c(30, 47)
scale <- 110
linewidth <- 1
# countries <- ne_countries(scale=scale, returnclass="sf")
# coast <- ne_coastline(scale=scale, returnclass="sf")
states <- ne_states(country="united states of america", returnclass="sf")


map1 <- ggplot() +
  geom_sf(data=states) +
  # geom_sf(data=coast) +
  geom_point(data=sampleData, mapping=aes(x=longitude, y=latitude, text=rownames(sampleData), color=as.factor(group))) +
  coord_sf(xlim=xlim1, ylim=ylim1, expand=FALSE)
map1
p <- ggplotly(map1)
p
saveWidget(p, "americanus-group-dapc-map.html")

xlim2 <- c(-87, -85)
ylim2 <- c(32.5, 34)
ggplot() +
  geom_sf(data=states) +
  geom_point(data=sampleData, mapping=aes(x=longitude, y=latitude, color=as.factor(group))) +
  coord_sf(xlim=xlim2, ylim=ylim2, expand=FALSE)
