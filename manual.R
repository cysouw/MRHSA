#' ---
#' title: "Historical reconstruction of MRHSA data"
#' author: "Michael Cysouw"
#' date: "`r Sys.Date()`"
#' ---

# produce output
# rmarkdown::render("manual.R")

#' ### necessary libraries

# require(qlcMatrix)
# require(qlcData)
# require(qlcVisualize)
# require(showtext)
# require(apcluster)

#' ### read data

# special function adapted to the details of the data
source("code/readData.R")

old <- read_mrhsa("sources/aeltere-generation-ipa.tsv")

# help function for visualisation
# - draw.cluster based on "limage"
# - draw.line to add separation lines into the plots
# - plot.cluster to send drawings to PDF
source("code/visualizeClusters.R")

#' ### visualise reconstructions

sounds <- sapply(strsplit(colnames(old$data), "_"), function(x) x[2])
draw.cluster("t", data = old$data, clusters = sounds)

# save all images as PDF
# sapply(sort(unique(sounds)), plot.cluster, data = old$data, clusters = sounds)

#' ### Similarity between alignments

# do not count shared gaps as similarity
# because then completely different alignments with many gaps get similar

tmp <- t(old$data)
tmp[tmp == "-"] <- NA
sim <- qlcMatrix::sim.obs(tmp, method = "res")
rm(tmp)

#' ### Clustering of alignments

# looking for clusters of alignments
clusters1 <- cutree(hclust(as.dist(-sim)),h = -0.01)
max(clusters1)

# This is an interesting alternative option
p <- apcluster::apcluster(sim)
clusters2 <- apcluster::labels(p, type = "enum")
max(clusters2)
rm(p)

# compare clusterings
compare <- table(clusters1, clusters2)
heatmap( -compare)

# relation of clustering to proposed reconstruction in the data

compare <- table(sounds, clusters1)
heatmap( -compare^.2)

compare <- table(sounds, clusters2)
heatmap( -compare^.2)

#' ### inspection of clusters

clusters <- clusters1

# most frequent sounds per cluster
stats(clusters, old$data)

# add one image to the output
draw.cluster(9, data = old$data, clusters)

# save all images as PDF
# sapply(1:max(clusters), plot.cluster, data = old$data, clusters = clusters)

#' ### compare clustering p/b/f with induced 2/3/12

# combinations of clusters
# plot.cluster(c(2,3,12), data = old$data, clusters)

# compare
# plot.cluster(c("b", "p", "f") , old$data, sounds)

#' ### phones per village

uniquePhones <- as.character(sort(unique(old$raw$ipa)))
phoneFreq <- apply(old$data, 1, function(x) { table(x)[uniquePhones] } )
phoneFreq[is.na(phoneFreq)] <- 0
rownames(phoneFreq) <- uniquePhones


