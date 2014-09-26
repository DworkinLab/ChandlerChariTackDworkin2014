#############################################
# Change these parameters as desired
#############################################

bin.size <- 50
overlap.size <- 5

#############################################
# Load the data
#############################################


column.names <- c("chrom", "pos", "sam.S", "sam.O", "ore.S", "ore.O", "bc.long.S", "bc.long.O", "bc.short1.S", "bc.short1.O", "bc.short2.S", "bc.short2.O", "bc.short3.S", "bc.short3.O", "bc.short4.S", "bc.short4.O", "bs.long.1.S", "bs.long.1.O", "bs.long.2.S", "bs.long.2.O", "bs.short.1.S", "bs.short.1.O", "bs.short.2.S", "bs.short.2.O")
column.types <- c("character", rep("numeric", 23))
snp.data <- read.csv("../../sequencing_scripts/snp_data/all_allele_freq_plot_data.csv", header=FALSE, col.names=column.names, colClasses=column.types)


#############################################
# Calculate sliding window allele frequencies
#############################################
calc.sliding.window.freqs <- function(cur.chrom, bin.size=100, overlap.size=10) {
	#Trim the data down to the sections we want
	cur.data <- snp.data[snp.data$chrom==cur.chrom,]
	
	#Figure out now many bins we'll need to create
	#It should be the number of SNPs in the current chromosome divided by the interval size [rounded up]
	#Interval size is the size of the bins (windows) minus the part that's
	positions.x <- rep(NA, ceiling(nrow(cur.data) / (bin.size - overlap.size - 1)))
	freq.data <- data.frame(chrom=cur.chrom, bin.start=positions.x, bin.end=positions.x, sam=positions.x, ore=positions.x, bc.long=positions.x, bc.short1=positions.x, bc.short2=positions.x, bc.short3=positions.x, bc.short4=positions.x, bc.short.avg=positions.x)
	
	#Loop through all the bins
	bin.start.index <- 1 #Initialize the starting index of the bin
	for (b in 1:(length(positions.x))) {

		#Calculate the index of the last SNP to include in this bin
		#Then save the coordinates of the start and end positions of this bin
		bin.end.index <- bin.start.index + bin.size - 1
		#Make sure this is a complete bin -- if there are not enough SNPs to make a complete bin here, we'll skip it
		if (bin.end.index <= nrow(cur.data)) {
			freq.data$bin.start[b] <- cur.data$pos[bin.start.index]
			freq.data$bin.end[b] <- cur.data$pos[bin.end.index]
			
			#Pull out the SNPs that fall within this bin
			bin.data <- cur.data[bin.start.index:bin.end.index,]
			
			#Get the stats on this bin
			freq.data$sam[b] <- sum(bin.data$sam.O) / (sum(bin.data$sam.O) + sum(bin.data$sam.S))
			freq.data$ore[b] <- sum(bin.data$ore.O) / (sum(bin.data$ore.O) + sum(bin.data$ore.S))
			freq.data$bc.short1[b] <- sum(bin.data$bc.short1.O) / (sum(bin.data$bc.short1.O) + sum(bin.data$bc.short1.S))
			freq.data$bc.short2[b] <- sum(bin.data$bc.short2.O) / (sum(bin.data$bc.short2.O) + sum(bin.data$bc.short2.S))
			freq.data$bc.short3[b] <- sum(bin.data$bc.short3.O) / (sum(bin.data$bc.short3.O) + sum(bin.data$bc.short3.S))
			freq.data$bc.short4[b] <- sum(bin.data$bc.short4.O) / (sum(bin.data$bc.short4.O) + sum(bin.data$bc.short4.S))
			freq.data$bc.long[b] <- sum(bin.data$bc.long.O) / (sum(bin.data$bc.long.O) + sum(bin.data$bc.long.S))
			short.S <- na.omit(c(bin.data$bc.short1.S, bin.data$bc.short2.S, bin.data$bc.short3.S, bin.data$bc.short4.S))
			short.O <- na.omit(c(bin.data$bc.short1.O, bin.data$bc.short2.O, bin.data$bc.short3.O, bin.data$bc.short4.O))
			freq.data$bc.short.avg[b] <- sum(short.O) / (sum(short.O) + sum(short.S))
		}
		bin.start.index <- bin.end.index - overlap.size
	}
	
	#Remove any empty rows from the data frame, which may have resulted if there weren't enough SNPs in the last bin to
	# make a full bin
	freq.data <- freq.data[is.finite(freq.data$bin.end),]
	
	return(freq.data)
}


############################################################################
############################################################################
############################################################################

bc.results <- NULL
chroms <- c("X", "2L", "2R", "3L", "3R")
for (c in chroms) {
	print(paste("Doing chromosome ", c, "...", sep=""))
	bc.results <- rbind(bc.results, calc.sliding.window.freqs(c, bin.size=bin.size, overlap.size=overlap.size))
}

bin.data.file <- paste("freq_data_bins_", as.character(bin.size), "_", as.character(overlap.size), ".csv", sep="")
write.csv(bc.results, bin.data.file, row.names=FALSE)

bc.results$bin.size <- bc.results$bin.end - bc.results$bin.start + 1
bc.results$pos <- round((bc.results$bin.start + bc.results$bin.end) / 2)



############################################################################
############################################################################
############################################################################

# Now do some plotting

mean(bc.results$bin.size)


##############################################################################
# Accessory functions for plotting
##############################################################################

get.chrom.size <- function(chrom) {
	result <- max(bc.results$pos[bc.results$chrom==chrom])
	return(result)
}


##############################################################################
# Important function: does the plotting!
##############################################################################

#Plots a specific backcross line
plot.line <- function(cur.bc, col.name, plot.color, adj.y, side.text, caption.x, sidebar.cex) {
	cur.color <- plot.color
	freq.data <- cur.bc[,c("pos", col.name)]
	colnames(freq.data) <- c("pos", "freq")
	freq.data$freq <- adj.y + freq.data$freq
	lines(freq ~ pos, data=na.omit(freq.data[,c("pos", "freq")]), cex=0.25, col="grey70", lty=3)
	points(freq ~ pos, data=freq.data, cex=0.5, col=cur.color)
	lines(x=c(-3e5, 3e10), y=c(adj.y, adj.y), lty=2, col="black")
	
	#par(srt=90, cex=sidebar.cex)
	#mtext(x=caption.x, y=0.5+adj.y, labels=side.text, col=cur.color)
	
	par(srt=90, cex=sidebar.cex)
	mtext(side.text, side=2, at=0.5+adj.y, col=cur.color, cex=sidebar.cex)
}

#Plots all the backcross line data for a particular chromosome region
plot.region <- function(chrom, xlim=NULL, gene.cex=0.5, sidebar.cex=0.8) {
	
	y.space <- 1.2
	
	#Long = 0 to 1
	#Short1 = 1.1 to 2.1
	#Short2 = 2.2 to 3.2
	#Short3 = 3.3 to 4.3
	#Short4 = 4.4 to 5.4
	
	if (is.null(xlim)) {
		xlim <- c(0, get.chrom.size(chrom))
	}
	
	cur.bc <- bc.results[bc.results$chrom==chrom,]
			
	#Set up the plot window
	par(mar=c(5,4,1,2)+0.1, cex=0.8)
	xlab <- paste("Position along chromosome", chrom)
	plot(x=NULL, y=NULL, xlim=xlim, ylim=c(0, y.space*4+1), xlab=xlab, ylab="Frequency of the ORE allele", yaxt="n")
	
	caption.x <- xlim[1] - (xlim[2]-xlim[1])/50

	plot.line(cur.bc, "bc.long", "blue", 0, "Long\nBC", caption.x, sidebar.cex)
	plot.line(cur.bc, "bc.short1", "red", y.space*1, "Short\nBC 1", caption.x, sidebar.cex)
	plot.line(cur.bc, "bc.short2", "red", y.space*2, "Short\nBC 2", caption.x, sidebar.cex)
	plot.line(cur.bc, "bc.short3", "red", y.space*3, "Short\nBC 3", caption.x, sidebar.cex)
	plot.line(cur.bc, "bc.short4", "red", y.space*4, "Short\nBC 4", caption.x, sidebar.cex)

}

##############################################################################
# Make the plots
##############################################################################

plot.chrom.pdf <- function(filename, width, height, chrom, xlim=NULL) {
	pdf(filename, width=width, height=height)
	plot.region(chrom, xlim=xlim)
	dev.off()
}


pdf.width <- 5
pdf.height <- 5

#Make a PDF for each chromosome
chroms <- c("X", "2L", "2R", "3L", "3R")
for (c in chroms) {
	cur.file <- paste("chrom_", c, "_",  as.character(bin.size), "_", as.character(overlap.size),".pdf", sep="")
	plot.chrom.pdf(cur.file, width=pdf.width, height=pdf.height, chrom=c)
}