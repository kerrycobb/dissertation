library(pophelper)
library(gridExtra)
# library(RColorBrewer)
# library(ggplot2)

setwd("~/Desktop/toad-project/structure/")
name <- "clust-90-indel-16-samples-198-include-americanus-group-1"
# name <- "clust-90-indel-16-samples-179-include-americanus-group-2"
dir <- paste0("~/Desktop/toad-project/structure/", name)

colors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#fdbf6f","#e31a1c")
colors <- c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF","#0089B2","#FDF060","#FFA6B2","#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368","#A4A4A4","#610B5E")

files <- list.files(dir, full.names=TRUE)
qlist <- readQ(files=files, filetype="structure", indlabfromfile=T)
t <- tabulateQ(qlist=qlist)
t$iter <- as.numeric(sapply(strsplit(t$file, "\\-|\\."), getElement, 14)) # Make column with iteration number
t <- t[order(t$k, t$iter), ] # sort dataframe on k and iteration
t
s <- summariseQ(t)
s

# Evanno
evanno <- evannoMethodStructure(data=s, exportplot=F, returnplot=T, returndata=F,
                                basesize=12, linesize=0.7)
grid.arrange(evanno)

# Align and merge
aligned <- alignK(qlist)
merged <- mergeQ(aligned)

write.csv(merged[[2]], paste0("qmat-", name, "-K-3.csv"))

# Show any samples with at least % fowleri/woodhousii to exclude from hybrid zone analysis
k3Df <- merged[[2]]
print(k3Df[which(k3Df$Cluster2 > 0.1),])



# stripPanelnames <- paste0("K=",sapply(slist1, function(x) attr(x,"k")))

# groupLabels <- data.frame(Species=sapply(strsplit(rownames(merged[[1]]), "_"), 
                                       # getElement, 2))

# # Plot All Runs
# plotQ(
#   aligned,
#   imgoutput="join",
#   returnplot=F,
#   exportplot=T,
#   exportpath=getwd(),
#   outputfilename=paste0("plot-", name, "-all-runs"),
#   clustercol=colors,
#   splab=paste0("K=", t$k, "\nRun ", t$iter),
#   basesize=11
# )

# Plot merged, all K
plotQ(
  merged,
  imgoutput="join",
  outputfilename=paste0("plot-", name, "-merged-all-k"),
  exportpath=getwd(),
  imgtype="pdf",
  exportplot=T,
  returnplot=F,
  clustercol=colors,
  splab=paste0("K=", names(merged)),
  sortind="all",
  sharedindlab=F,
  #   exportplot=F, 
  #   basesize=11,
  #   clustercol=colors,
  #   
  # #  showindlab=T,
  # #  sharedindlab=T
)



# Plot merged K=
k <- 5
# legendLabels <- c("terrestris", "fowleri/woodhousii", "americanus")

plotQ(
  merged[k-1], 
  exportplot=T, 
  outputfilename=paste0("plot-", name, "-K-", k),
  imgtype="pdf",
  exportpath=getwd(),
  clustercol=colors, 
  returnplot=F, 
  showsp=F, # Doesn't work
  splab="",
  # splab=stripPanelnames[2],
  # legendlab=legendLabels,
  showlegend=T,
  # legendtextsize=4,
  sortind="all",
  showticks=T,
  # ticksize=0.1,
  showindlab=T,
  useindlab=T,
  indlabsize=1,
  indlabspacer=0.25,
  barbordercolour="black",
  width=10,
  height=3,
   
)
