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
#	selection <- sample(which(old$loc[,1]<7.5), size)
	neighbor <- apply(graphdist[selection,], 1, sample_neighbor)
	return(cbind(selection, neighbor))
	
}

# MASLOVA

p <- function(x, feature, sounds, graph, dist = 1, samplesize = 50) {

	feature[!(feature %in% sounds)] <- NA

	pairs <- select_pairs(graph, dist, samplesize)	
	available <- !is.na(feature[pairs[,1]]) & !is.na(feature[pairs[,2]])
	pairs <- pairs[available,]
		
	fA <- sum(feature[pairs[,1]] == sounds[1], na.rm = TRUE)
	fB <- sum(feature[pairs[,1]] == sounds[2], na.rm = TRUE)
	pA <- fA / (fA + fB)
	
	fAB <- feature[pairs[,1]] == feature[pairs[,2]]	
	pD <- sum(!fAB, na.rm = TRUE) / sum(!is.na(fAB))

	return(c(pA = pA, pD = pD))
}

pred <- function(x, nsamples = 50, ...) {
	
	tmp <- sapply(1:nsamples, p, ...)
	coef <- lm(tmp[2,] ~ tmp[1,])$coefficients
	cor <- cor.test(tmp[2,], tmp[1,])$estimate
	
	a <- coef[2]/2
	b <- coef[1]/2
	
	pAB <- (1 + a - sqrt( (1-a)^2 -4*b) )/2
	pBA <- (1 - a - sqrt( (1-a)^2 -4*b) )/2

	result <- c(a,b,cor,pAB,pBA)
	names(result) <- c("alpha", "beta", "cor", "pAB", "pBA")	
	return(result)
}

attempts_at_dist <- function(distance, nattempts = 20, ...) {
	
	all <- sapply(1:nattempts, pred, dist = distance, ...)
	med <- apply(all, 1, median)
	return(med)

} 

source("code/readData.R")

loc <- read_loc("sources/mrhsa/mrhsa-gid-wkt.tsv")
old <- read_mrhsa("sources/mrhsa/aeltere-generation-ipa.tsv", loc)

d <- make_graphdist(old)

estimates <- sapply(1:5
			, attempts_at_dist
			, feature = old$data[,"gibt_b_(52.5)"]
			, sounds = c("p","-")
			, graph = d
			)

estimates
plot(c(0,estimates["pAB",]),c(0,estimates["pBA",]))
lm(estimates["pBA",] ~ 0 + estimates["pAB",])$coefficients

plot(old$loc, type = "n")
text(old$loc, labels = old$data[,"gibt_b_(52.5)"])

compare("Kupfer_p_(98.3)")
compare("Pfennig_p_(30.4)")
compare("pfeifen_p_(171.2)")

compare("gestorben_b_(140.3)")
compare("aufräumen_p_(27.31)")

compare("fünfzig_f_(2.3)")
[1] "ab_b_(91.5)"         "Arbeit_b_(51.3)"     "Gabel_b_(33.2)"      "geben_b_(28.3)"     
[5] "gestorben_b_(140.3)" "gibt_b_(52.5)"       "Körbe_b_(106.3)"     "taub_b_(128.2)"     
[9] "Taube_b_(166.5)"   

# PHOIBLE

phoible <- read.csv("sources/phoible/cldf/parameters.csv")

# remove diacritica to match phones to phoible

old$data <- gsub("͡","",old$data)
old$data <- gsub("̠","",old$data)

avail <- table(old$data[,5])
sel <- phoible$Name %in% names(avail)
tmp <- phoible[sel,-c(1:4)]
rownames(tmp) <- as.character(phoible$Name)[sel]
t(tmp)