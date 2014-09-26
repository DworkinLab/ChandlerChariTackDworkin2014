############################################################################
############################################################################
############################################################################

num.bins <- 250
bin.width <- (120e6 / 5) / num.bins
blank.col <- rep(NA, num.bins)
bc.results <- data.frame(		chrom=rep("2L", num.bins),
								bin=seq(from=1, to=num.bins),
								pos=seq(from=1, by=bin.width, length.out=num.bins),
								bc.short1=blank.col,
								bc.short2=blank.col,
								bc.short3=blank.col,
								bc.short4=blank.col,
								bc.long=blank.col
					    )


simulate.short.results <- function(num.bins, peak.locations, peak.widths) {
	#By default, short backcross lines should have all SAM alleles except for at QTLs
	#Create the baseline results
	short.results <- rep(0, num.bins)
	
	#Now add peaks
	for (c in 1:length(peak.locations)) {
		cur.peak <- peak.locations[c]
		short.results[cur.peak] <- 1
		for (i in 1:peak.widths[c]) {
			cur.freq <- 1 - (i / peak.widths[c])
			sample.size <- 16
			short.results[cur.peak + i] <- rbinom(1, size=sample.size, prob=cur.freq) / sample.size
			short.results[cur.peak - i] <- rbinom(1, size=sample.size, prob=cur.freq) / sample.size
		}
	}
	
	return(short.results)
}

peak.locations <- c(0.23, 0.66) * num.bins
peak.widths <- c(8, 7)

width.mean <- 5
width.sd <- 1.5

bc.results$bc.short1 <- simulate.short.results(num.bins, peak.locations, peak.widths=rnorm(2, mean=width.mean, sd=width.sd))
bc.results$bc.short2 <- simulate.short.results(num.bins, peak.locations, peak.widths=rnorm(2, mean=width.mean, sd=width.sd))
bc.results$bc.short3 <- simulate.short.results(num.bins, peak.locations, peak.widths=rnorm(2, mean=width.mean, sd=width.sd))
bc.results$bc.short4 <- simulate.short.results(num.bins, peak.locations, peak.widths=rnorm(2, mean=width.mean, sd=width.sd))
bc.results$bc.long <- 1 - simulate.short.results(num.bins, peak.locations, peak.widths=rnorm(2, mean=width.mean, sd=width.sd))


############################################################################
############################################################################
############################################################################

# Now do some plotting

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
	lines(freq ~ pos, data=freq.data, cex=0.5, col=cur.color)
	lines(x=c(-3e5, 3e10), y=c(adj.y, adj.y), lty=2, col="black")
	
	#par(srt=90, cex=sidebar.cex)
	#text(x=caption.x, y=0.5+adj.y, labels=side.text, col=cur.color)
	
	par(srt=90, cex=sidebar.cex)
	mtext(side.text, side=2, at=0.5+adj.y, col=cur.color, cex=sidebar.cex)
}

#Plots all the backcross line data for a particular chromosome region
plot.region <- function(xlim=NULL, gene.cex=0.5, sidebar.cex=0.8) {
	
	y.space <- 1.2
	
	#Long = 0 to 1
	#Short1 = 1.1 to 2.1
	#Short2 = 2.2 to 3.2
	#Short3 = 3.3 to 4.3
	#Short4 = 4.4 to 5.4
	
	if (is.null(xlim)) {
		xlim <- c(0, max(bc.results$pos))
	}
	
	#Set up the plot window
	par(mar=c(5,4,1,2)+0.1, cex=0.8)
	xlab <- paste("Position along chromosome")
	plot(x=NULL, y=NULL, xlim=xlim, ylim=c(0, y.space*4+1), xlab=xlab, ylab="Frequency of the ORE allele", yaxt="n")
	caption.x <- xlim[1] - (xlim[2]-xlim[1])/50

	plot.line(bc.results, "bc.long", "blue", 0, "Long\nBC", caption.x, sidebar.cex)
	plot.line(bc.results, "bc.short1", "red", y.space*1, "Short\nBC 1", caption.x, sidebar.cex)
	plot.line(bc.results, "bc.short2", "red", y.space*2, "Short\nBC 2", caption.x, sidebar.cex)
	plot.line(bc.results, "bc.short3", "red", y.space*3, "Short\nBC 3", caption.x, sidebar.cex)
	plot.line(bc.results, "bc.short4", "red", y.space*4, "Short\nBC 4", caption.x, sidebar.cex)

}

##############################################################################
# Make the plots
##############################################################################

plot.chrom.pdf <- function(filename, width, height, xlim=NULL) {
	pdf(filename, width=width, height=height)
	plot.region(chrom, xlim=xlim)
	dev.off()
}




pdf.width <- 5
pdf.height <- 5


plot.chrom.pdf(file="Mock_results.pdf", pdf.width, pdf.height)