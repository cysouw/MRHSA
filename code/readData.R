# read data and make some basic transformations
# this is not very optimized, so it might be a bit slow

read_mrhsa <- function(file) {
	
	raw <- read.delim(file, quote="", comment.char="#", na.strings = "Phone")
	colnames(raw) <- c("word", "alignment", "ipa", "place", "gid", "ID", "text")
	
	# remove duplicated data, which only differ in 'gid'
	
	raw <- raw[!duplicated(raw[,-6]),]
	
	# extract indivisual correspondences:
	# combinatino of a sound in a word in a sentence
	
	sentence <- sub("^.+\\(", "(", raw$text)
	sound <- sub("^.+ ", "", raw$alignment)
	align <- factor(paste(raw$word, sound, sentence, sep = "_"))

	# make 'wide' datafile

	data <- matrix(
				  nrow = nlevels(raw$place)
				, ncol = nlevels(align)
				, dimnames = list(levels(raw$place), levels(align))
				)

	p <- as.numeric(raw$place)
	a <- as.numeric(align)
	for (i in seq_along(raw$ipa)) {
		data[p[i],a[i]] <- as.character(raw$ipa[i])
		}
	
	# return both the raw data and the wide datafile
	
	return(list(data = data, raw = raw))
	
}
