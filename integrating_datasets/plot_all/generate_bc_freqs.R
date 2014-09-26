column.names <- c("chrom", "pos", "sam.S", "sam.O", "ore.S", "ore.O", "bc.long.S", "bc.long.O", "bc.short1.S", "bc.short1.O", "bc.short2.S", "bc.short2.O", "bc.short3.S", "bc.short3.O", "bc.short4.S", "bc.short4.O", "bs.long.1.S", "bs.long.1.O", "bs.long.2.S", "bs.long.2.O", "bs.short.1.S", "bs.short.1.O", "bs.short.2.S", "bs.short.2.O")
column.types <- c("character", rep("numeric", 23))
snp.data <- read.csv("../../sequencing_scripts/snp_data/all_allele_freq_plot_data.csv", header=FALSE, col.names=column.names, colClasses=column.types)



#############################################
# Calculate sliding window allele frequencies
#############################################
calc.sliding.window.freqs <- function(cur.chrom, window.size=10000, step.size=2000, p.thresh=0.001, xlim=NULL) {
	#Trim the data down to the sections we want
	cur.data <- snp.data[snp.data$chrom==cur.chrom,]
	if (!is.null(xlim)) {
		cur.data <- cur.data[cur.data$pos >= xlim[1] & cur.data$pos <= xlim[2],]
	}
		
	#Set up our bins and data table
	if (is.null(xlim)) {
		positions.x <- seq(from=1, to=max(cur.data$pos), by=step.size)
	} else {
		positions.x <- seq(from=xlim[1], to=xlim[2], by=step.size)
	}
	empties <- rep(NA, length(positions.x))
	freq.data <- data.frame(chrom=cur.chrom, pos=positions.x, sam=empties, ore=empties, bc.long=empties, bc.short1=empties, bc.short2=empties, bc.short3=empties, bc.short4=empties, bc.short.avg=empties)
	
	if (is.null(xlim)) {
		xlim=c(1, max(cur.data$pos))
	}
	
	#Loop through all the bins
	for (b in 1:(length(positions.x)-1)) {

		#Calculate the boundaries for this bin
		pos.min = positions.x[b] - (window.size / 2)
		pos.max = positions.x[b] + (window.size + 2)
		
		#Pull out the SNPs that fall within this bin
		bin.data <- cur.data[(cur.data$pos >= pos.min) & (cur.data$pos < pos.max),]
		
		#Make sure there are SNPs within this bin
		if (nrow(bin.data) > 0) {
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
	}
	
	return(freq.data)
}


############################################################################
############################################################################
############################################################################

all.data <- NULL
chroms <- c("X", "2L", "2R", "3L", "3R")
for (c in chroms) {
	print(paste("Doing chromosome ", c, "...", sep=""))
	all.data <- rbind(all.data, calc.sliding.window.freqs(c))
}

write.csv(all.data, "freq_data_10kbwindows_2kbsteps.csv", row.names=FALSE)