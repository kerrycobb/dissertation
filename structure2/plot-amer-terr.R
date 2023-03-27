library(pophelper)
library(gridExtra)
# library(RColorBrewer)
# library(ggplot2)

setwd("~/Desktop/toad-project/structure2/")
name <- "clust-90-missing-0.9-exclude-1"

dir <- paste0("~/Desktop/toad-project/structure2/", name)

colors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#fdbf6f","#e31a1c")
colors <- c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF","#0089B2","#FDF060","#FFA6B2","#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368","#A4A4A4","#610B5E")

files <- list.files(dir, full.names=TRUE)
qlist <- readQ(files=files, filetype="structure", indlabfromfile=T)
t <- tabulateQ(qlist=qlist)
t$iter <- as.numeric(sapply(strsplit(t$file, "\\-|\\."), getElement, 11)) # Make column with iteration number
t <- t[order(t$k, t$iter), ] # sort dataframe on k and iteration
t
s <- summariseQ(t)
s

# Evanno
evanno <- evannoMethodStructure(data=s, exportplot=F, returnplot=T, returndata=F, basesize=12, linesize=0.7)
grid.arrange(evanno)

# Align and merge
aligned <- alignK(qlist)
merged <- mergeQ(aligned)

write.csv(merged[[1]], paste0("qmat-", name, "-K-2.csv"))

# Plot All Runs
plotQ(
  aligned,
  imgoutput="join",
  returnplot=F,
  exportplot=T,
  imgtype="pdf",
  exportpath=getwd(),
  outputfilename=paste0("plot-", name, "-all-runs"),
  clustercol=colors,
  splab=paste0("K=", t$k, "\nRun ", t$iter),
  basesize=11
)

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
    sortind="Cluster1",
    sharedindlab=F,
  #   exportplot=F, 
  #   basesize=11,
  #   clustercol=colors,
  #   
  # #  showindlab=T,
  # #  sharedindlab=T
)


# Plot merged, K=2
legendLabels <- c("terrestris", "americanus")
plotQ(
  merged[1], 
  exportplot=T, 
  outputfilename=paste0("plot-", name, "-K-2"),
  imgtype="pdf",
  exportpath=getwd(),
  clustercol=colors, 
  returnplot=F, 
  showsp=F, # Doesn't work
  splab="",
  # splab=stripPanelnames[2],
  legendlab=legendLabels,
  showlegend=T,
  # legendtextsize=4,
  sortind="Cluster1",
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

