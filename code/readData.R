# read data and make some basic transformations
# this is not very optimized, so it might be a bit slow

read_loc <- function(file) {
	
	loc <- read.delim(file)
	
	l <- gsub("POINT\\(","",loc$X.wkt)
	l <- gsub("\\)","",l)
	l <- strsplit(l, split = " ")
	l <- sapply(l, as.numeric)
	l <- t(l)
	
	colnames(l) <- c("long", "lat")
	rownames(l) <- loc$X.gid

	return(l)
}

read_mrhsa <- function(file, locations) {
	
	raw <- read.delim(file, quote="", comment.char="#", na.strings = "Phone")
	colnames(raw) <- c("word", "alignment", "ipa", "place", "gid", "ID", "text")
	
	# extract indivisual correspondences:
	# combinatino of a sound in a word in a sentence
	
	sentence <- sub("^.+\\(", "(", raw$text)
	sound <- sub("^.+ ", "", raw$alignment)
	align <- factor(paste(raw$word, sound, sentence, sep = "_"))

	# there are some different villages with the same name...
	
	places <- as.factor(paste(raw$place, raw$gid))

	# make 'wide' datafile

	data <- matrix(
				  nrow = nlevels(places)
				, ncol = nlevels(align)
				, dimnames = list(levels(places), levels(align))
				)

	p <- as.numeric(places)
	a <- as.numeric(align)
	for (i in seq_along(raw$ipa)) {
		data[p[i],a[i]] <- as.character(raw$ipa[i])
		}
	
	# make unique and match to locations
	
	n <- gsub("^.+? ", "", rownames(data))
	selection <- which(n %in% rownames(locations))
	
	data <- data[selection,]
	data <- unique(data)
	
	n <- gsub("^.+? ", "", rownames(data))
	loc <- locations[n,]
	
	# return both the raw data and the wide datafile
	
	return(list(data = data, loc = loc, raw = raw))
	
}

