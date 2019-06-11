library(deldir)
library(sna)

make_graphdist <- function(data) {
	
	graph <- deldir::deldir(data$loc[,"long"], data$loc[,"lat"])[[1]]

	edges <- rbind(graph[,5:6], graph[,6,5])
	edges <- as.matrix(edges)
	edges <- rbind(edges, edges[,2:1])
	edges <- cbind(edges, dist = rep(1,times=dim(edges)[1]))

	attr(edges,"n") <- nrow(data$loc)
	attr(edges,"vnames") <- rownames(data$data)

	graphdist <- sna::geodist(edges)[[2]]

	return(graphdist)
	
}

sample_neighbor <- function(x) {
	
	possible <- which(x != 0)
	neighbor <- sample(possible, 1)
	return(neighbor)
	
}

select_pairs <- function(graphdist, distance, size) {
	
	graphdist[graphdist != distance] <- 0
	selection <- sample(dim(graphdist)[1], size)
	neighbor <- apply(graphdist[selection,], 1, sample_neighbor)
	return(cbind(selection, neighbor))
	
}

source("code/readData.R")

loc <- read_loc("sources/mrhsa/mrhsa-gid-wkt.tsv")
old <- read_mrhsa("sources/mrhsa/aeltere-generation-ipa.tsv", loc)

d <- make_graphdist(old)
s <- select_pairs(d,2,30)

phoible <- read.csv("sources/phoible/cldf/parameters.csv")

# remove diacritica to match phones to phoible

old$data <- gsub("͡","",old$data)
old$data <- gsub("̠","",old$data)

avail <- table(old$data[,5])
sel <- phoible$Name %in% names(avail)
tmp <- phoible[sel,-c(1:4)]
rownames(tmp) <- as.character(phoible$Name)[sel]
t(tmp)