require(deldir)
require(sna)
require(corHMM)
require(qlcMatrix)
require(qlcVisualize)
require(markovchain)

make_graphdist <- function(data) {
	
	graph <- deldir::deldir(data$loc[,"long"], data$loc[,"lat"])[[1]]

	edges <- rbind(graph[,5:6], graph[,6,5])
	edges <- as.matrix(edges)
	edges <- rbind(edges, edges[,2:1])
	edges <- cbind(edges, dist = rep(1,times=dim(edges)[1]))

	attr(edges,"n") <- nrow(data$loc)
	attr(edges,"vnames") <- rownames(data$data)

	graphdist <- sna::geodist(edges)[[2]]
	rownames(graphdist) <- colnames(graphdist) <- rownames(data$data)

	return(graphdist)
	
}

sample_neighbor <- function(x) {
	
	possible <- which(x != 0)
	neighbor <- sample(possible, 1)
	return(neighbor)
	
}

select_pairs <- function(graphdist, distance, size) {
	
	graphdist[graphdist != distance] <- 0
	selection <- sample(nrow(graphdist), size)
#	selection <- sample(which(old$loc[,1]<7.5), size)
	neighbor <- apply(graphdist[selection,], 1, sample_neighbor)
	return(cbind(selection, neighbor))
	
}

# Maslova estimation

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

maslova_pair_estimates <- function(feature, sounds, full = FALSE) {
	
	estimates <- sapply(1:5
			, attempts_at_dist
			, feature = feature
			, sounds = sounds
			, graph = d[names(feature),names(feature)]
			)
	
	if (full) {
		return(estimates)
	} else {
		tAB <- lm(estimates["pAB",] ~ 0 + c(1:5))$coefficients
		tBA <- lm(estimates["pBA",] ~ 0 + c(1:5))$coefficients

		return(c(tAB,tBA))
	}
}

maslova_estimates <- function(feature) {
	
	sounds <- names(table(feature))
	s <- length(sounds)
	result <- matrix(NA, s, s, dimnames = list(sounds, sounds))
	
	for ( soundA in 1:(s-1) ) {
		for ( soundB in (soundA+1):s ) {
			selectSounds <- c(sounds[soundA], sounds[soundB])
			estims <- maslova_pair_estimates(feature, selectSounds)
			result[soundA, soundB] <- estims[1]
			result[soundB, soundA] <- estims[2]
		}
	}

	return(result)
}

markov_estimates <- function(feature, root = NA) {
	
	# village tree based on general similarity, subsetted to selected villages
	tree <- ape::nj(as.dist(1-qlcMatrix::sim.obs(old$data[names(feature),])))
	# branch lengths fixed to 1 to allow comparability with graph distance
	tree$edge.length <- rep(1, times = length(tree$edge.length))
	# optionally reroot
	if (!is.na(root)) {
		tree <- root(tree, root)
	}
	
	# specail dataformat as expected by rayDISC
	feature[feature == "-"] <- "_"
	data <- cbind(names(feature), feature)
	data[is.na(data)] <- "?"

	model <- corHMM::rayDISC(tree, data, model = "ARD", node.state = "marginal")

	result <- model$solution
	rownames(result)[rownames(result) == "_"] <- "-"
	colnames(result)[colnames(result) == "_"] <- "-"
	return(result)
}

clean_feature <- function(feature, cutoff = 10, long = NA) {

	feature <- old$data[,feature]
	names(feature) <- rownames(old$data)
	
	if (is.numeric(long)) {
		feature <- feature[old$loc[,1] < long]
	}
	
	values <- table(feature)
	rare_values <- which(values < cutoff)
	feature[feature %in% names(values)[rare_values]] <- NA
	return(feature)
}

compare <- function(sound) {
		
	if (is.numeric(sound)) {
		sound <- colnames(old$data)[sound]
	}
	
	if (sound %in% colnames(new$data)) {
		
		new_sound <- new$data[,sound]
		old_sound <- old$data[names(new_sound),sound]
	
		freq <- table(old_sound,new_sound)
		perc <- round(prop.table(freq,1),3)*100
	
		return(list(sound = sound, frequency = freq, percentage = perc))
	
	} else {
		warning("Sound is not available in data for younger generation")
	}
}

stable <- function(model) {
	
	diag(model) <- -rowSums(model, na.rm = TRUE)
	model[is.na(model)] <- 0
	Q <- new("ctmc", generator = model)
	markovchain::steadyStates(Q)
	
}

# ==========================

source("code/readData.R")

loc <- read_loc("sources/mrhsa/mrhsa-gid-wkt.tsv")
old <- read_mrhsa("sources/mrhsa/aeltere-generation-ipa.tsv", loc)
new <- read_mrhsa("sources/mrhsa/juengere-generation-ipa.tsv", loc)

# graph-based geographical distance between villages
d <- make_graphdist(old)

sound <- 215
qlcVisualize::lmap(old$loc, old$data[,sound], levels = .6)
compare(sound)

feature <- clean_feature(sound, 20)
table(feature)/sum(table(feature))

(maslova_estimates(feature) -> q)
(markov_estimates(feature) -> r)

stable(q)
stable(t(q))


# PHOIBLE

phoible <- read.csv("sources/phoible/cldf/parameters.csv")

# remove diacritica to match phones to phoible

old$data <- gsub("͡","",old$data)
old$data <- gsub("̠","",old$data)

avail <- table(old$data[,sound])
sel <- phoible$Name %in% names(avail)
tmp <- phoible[sel,-c(1:4)]
rownames(tmp) <- as.character(phoible$Name)[sel]
diff_test <- function(feat){max(table(feat)/nrow(tmp)) != 1}
interesting <- which(apply(tmp, 2, diff_test))

