# visualisation of clusters
# using function "limage" with different defaults

draw.cluster <- function(cluster, data, clusters, window = TRUE) {

	if(window) graphics.off()
	columns <- data[ , clusters %in% cluster]				
	qlcVisualize::limage(columns
						, col = rainbow(8)
						, order = "R2E"
						, method = "res"
						, show.remaining = TRUE
						, cex.axis = .5
						, cex.remaining = .5
						, cex.legend = 1
						, font = "Charis SIL"
						)
}

# simple function to draw lines in limage-plot

draw.line <- function(h, lwd = 1) {
	segments(rep(0,length(h))
			, rep(h, length(h))
			, rep(183, length(h))
			, rep(h, length(h))
			, lwd = lwd
			)
}

# drawing to PDF

plot.cluster <- function(cluster, data, clusters) {
	
	# PDF device
	
	name <- paste0(cluster, collapse = "+")
	name <- gsub("/","",name)
	
	filename <- paste0("images/", name, ".pdf")
	dev.new(file = filename 
			, width  =  nrow(data)/10 + 3
			, height =  sum(clusters %in% cluster)/10 + 3
			, type = "pdf"
			)

	# plotting the image
		
	sysfonts::font_add(family = "Charis SIL", regular = "CharisSIL-R.ttf")
	showtext::showtext_begin()

	draw.cluster(cluster, data = data, clusters, window = FALSE)
	
	showtext::showtext_end()
	dev.off()
}


# function to list most frequent sounds per cluster

stats <- function(clusters, alignments) {
	sapply(unique(clusters), function(clus) {
		c(sort(table(alignments[,clusters == clus]), decreasing = TRUE)[1:5]
		, cols = sum(clusters == clus)
		)
	}, simplify = FALSE)
}

