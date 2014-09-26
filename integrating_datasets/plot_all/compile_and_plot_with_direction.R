##############################################################################
# Read in flybase (list of genes with fbgn/names/symbols, and chromosomal locations)
##############################################################################

flybase.cols <- c("fbgn", "symbol", "chrom", "start", "stop")
flybase.colClasses <- c("character", "character", "character", "numeric", "numeric")
flybase.data <- read.csv("../../expression_scripts/flybase/flybase_map_table.csv", header=FALSE, colClasses=flybase.colClasses)
#flybase.data <- read.csv("flybase_map_table.csv", header=FALSE, colClasses=flybase.colClasses)
colnames(flybase.data) <- flybase.cols
newstart <- ifelse(flybase.data$start < flybase.data$stop, flybase.data$start, flybase.data$stop)
newstop <- ifelse(flybase.data$stop > flybase.data$start, flybase.data$stop, flybase.data$start)
flybase.data$start <- as.numeric(newstart)
flybase.data$stop <- as.numeric(newstop)
rm(newstart)
rm(newstop)
#Want the flybase data to be sorted by chromosome and then starting coordinates
flybase.data <- flybase.data[order(flybase.data$start),]
flybase.data <- flybase.data[order(flybase.data$chrom),]
flybase.data$y <- rep(seq(from=0.1, to=0.9, by=0.1), length.out=nrow(flybase.data)) #These values will be used for plotting



##############################################################################
# Read in the DGRC, Illumina, DGE datasets
# Here we're interested in background-dependence, so we want to use the p-values associated with the interactions
##############################################################################

dgrc.results <- read.csv('../../expression_scripts/microarray_data/DGRC_summary_results.csv', header=TRUE)
#dgrc.results <- read.csv('DGRC_summary_results.csv', header=TRUE)
dgrc.results <- dgrc.results[,2:ncol(dgrc.results)]
dgrc.results <- dgrc.results[grep("FBgn", dgrc.results$fbgn),]
dgrc.results$fbgn <- as.character(dgrc.results$fbgn)

illumina.results <- read.csv('../../expression_scripts/microarray_data/illumina_summary_results.csv', header=TRUE)
#illumina.results <- read.csv('illumina_summary_results.csv', header=TRUE)
colnames(illumina.results)[1] <- "gene"
illumina.results <- illumina.results[grep("FBgn", illumina.results$fbgn),]
illumina.results$fbgn <- as.character(illumina.results$fbgn)

dge.results <- read.csv('../../expression_scripts/dge_r_analysis/dge_allele_specific_results.csv', header=TRUE)
#dge.results <- read.csv('dge_allele_specific_results.csv', header=TRUE)
dge.results$fbgn <- as.character(dge.results$fbgn)
dge.results$gene <- as.character(dge.results$gene)
dge.results <- dge.results[is.finite(dge.results$p.2),] #Remove those that we couldn't test for an interaction

dge.genotype.results <- read.csv('../../expression_scripts/dge_r_analysis/dge_genotype_results.csv')
#dge.genotype.results <- read.csv('dge_genotype_results.csv')

##############################################################################
# Read in SRC modifier screen results
##############################################################################

src.cols <- c("stock", "p.deletion", "p.background", "p.interaction", "chrom", "start", "stop", "q.deletion", "q.interaction")
src.colclasses <- c("character", "numeric", "numeric", "numeric", "character", "numeric", "numeric", "numeric", "numeric")
src.results <- read.csv("../../src_scripts/raw_deletion_results.csv", header=TRUE, colClasses=src.colclasses)
#src.results <- read.csv("raw_deletion_results.csv", header=TRUE, colClasses=src.colclasses)
colnames(src.results) <- src.cols
#Sort and add y-values for plotting
src.results <- src.results[order(src.results$start),]
src.results <- src.results[order(src.results$chrom),]
src.results$y <- rep(seq(from=0.1, to=0.9, by=0.1), length.out=nrow(src.results))

################################################################################
### Read in binding site data
################################################################################

binding.cols <- c("fbgn", "llr.sam", "llr.ore", "llr.diff", "gc", "size", "empirical.p", "empirical.p.sam", "empirical.p.ore", "empirical.q.sam", "empirical.q.ore", "empirical.q")
binding.colclasses <- c("character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric")
binding.results <- read.csv("../../sequencing_scripts/binding_sites_custom_meme/SAM_ORE_binding_results_bygene_custom_meme.csv", header=TRUE, colClasses=binding.colclasses)
#binding.results <- read.csv("SAM_ORE_binding_results_bygene_custom_meme.csv", header=TRUE, colClasses=binding.colclasses)
colnames(binding.results) <- binding.cols
binding.results$min.q <- ifelse(binding.results$empirical.q.ore < binding.results$empirical.q.sam, binding.results$empirical.q.ore, binding.results$empirical.q.sam)

###############################################################################
## Read in backcross mapping data
###############################################################################

#Across all short backcross replicates

bc.results <- read.csv("freq_data_10kbwindows_2kbsteps.csv")


##############################################################################
# Gets overlapping deletions for a particular gene
##############################################################################
get.overlapping.deletions <- function(cur.chrom, cur.start, cur.stop) {
	
	cur.src.results <- src.results[src.results$chrom==cur.chrom,]
	#A deletion overlaps the desired region if ONE of the following is met:
	#  -Deletion contains entire region
	#  -Entire region contains deletion
	#  -Region starts between deletion's start and end
	#  -Region ends between deletion's start and end
	keeper.conditions <-     ((cur.src.results$start <= cur.start) & (cur.src.results$stop >= cur.stop)) |
	                         ((cur.start <= cur.src.results$start) & (cur.stop >= cur.src.results$stop)) |
	                         ((cur.start >= cur.src.results$start) & (cur.start <= cur.src.results$stop)) |
	                         ((cur.stop >= cur.src.results$start)  & (cur.stop <= cur.src.results$stop))
	cur.src.results <- cur.src.results[keeper.conditions,]
	return(cur.src.results)
}


##############################################################################
# Accessory functions for plotting
##############################################################################

get.chrom.size <- function(chrom) {
	result <- max(bc.results$pos[bc.results$chrom==chrom])
	return(result)
}


##############################################################################
# Custom plotting function
##############################################################################

#Draws a bar on the graph to represent a gene, deletion, etc.
draw.bar <- function(x, y, col, border) {
	new.lwd <- 0.025
	if (length(x) >= 2) {
		for (i in 1:(length(x)-1)) {
			rect(xleft=x[i], ybottom=y[i] - new.lwd, xright=x[i+1], ytop=y[i+1] + new.lwd, col=col, lwd=0.5, border=border)
		}
	}
}


##############################################################################
# Important function: does the plotting!
##############################################################################

#Plots all the results for a chromosome region
#1. Backcross data
#2. Genes that are significant for DGE 
#3. Genes that are significant for Illumina
#4. Genes that are significant for DGRC
#5. Significant modifier deletions 
#6. Binding sites
plot.region <- function(chrom, xlim=NULL, thresh=0.01, thresh.illumina=NULL, thresh.dgrc=NULL, thresh.dge.genotype=NULL, thresh.dge.interaction=NULL, thresh.src=NULL, thresh.binding=NULL, plot.gene.names=FALSE, gene.cex=0.5, sidebar.cex=0.8) {
	
	if (is.null(xlim)) {
		xlim <- c(0, get.chrom.size(chrom))
	}
	
	#Only draw those genes that meet a certain threshold for being significant in multiple datasets
	# AND for which expression effects are consistent across datasets
	meets.thresh <- sorted.flybase.data$sig.tally >= min.sig.tally
	all.expression.consistent <- is.na(sorted.flybase.data$inconsistent.expression.effects) | (sorted.flybase.data$inconsistent.expression.effects == FALSE)
	sig.fbgn.list <- sorted.flybase.data$fbgn[meets.thresh & all.expression.consistent]
	
	#Narrow down our focus to features in the requested region
	cur.flybase <- flybase.data[(flybase.data$chrom==chrom) & (flybase.data$fbgn %in% sig.fbgn.list),]
	cur.src <- get.overlapping.deletions(chrom, xlim[1], xlim[2])
	cur.dge <- dge.results[(dge.results$fbgn %in% cur.flybase$fbgn),]
	cur.dge.genotype <- dge.genotype.results[(dge.genotype.results$fbgn %in% cur.flybase$fbgn),]
	cur.illumina <- illumina.results[(illumina.results$fbgn %in% cur.flybase$fbgn),]
	cur.dgrc <- dgrc.results[(dgrc.results$fbgn %in% cur.flybase$fbgn),]
	cur.bc <- bc.results[bc.results$chrom==chrom,]
	cur.binding <- binding.results[(binding.results$fbgn %in% cur.flybase$fbgn),]
	
	#Some plotting options
	gene.y <- 0.06 #How far above the line to plot gene names
	
	#Where to plot each of the different sources of data, and what colors to use...
	#Colors consist of a triplet -- genotype-effect-only color, interaction-effect-only color, and color for genes with both genotype & interaction effects
	alpha.channels <- c(0, 0.5, 1.0)
	
	src.y <- 1
	cur.color <- c(255, 0, 0)/255
	src.color <- c(NA, NA, NA)
	for (i in 1:3) {
		src.color[i] <- rgb(cur.color[1], cur.color[2], cur.color[3], alpha.channels[i])
	}
	#src.color <- c(	rgb(0.9, 0.1, 0.1, 0.0), rgb(0.9, 0.1, 0.1, 0.33), rgb(0.9, 0.1, 0.1, 1.0) ) #c("#FFFFFF", "#888888", "#000000") #2
	src.borders <- c(rgb(cur.color[1], cur.color[2], cur.color[3]), NA, NA)
	
	illumina.y <- 4
	cur.color <- c(0, 100, 0)/255
	illumina.color <- c(NA, NA, NA)
	for (i in 1:3) {
		illumina.color[i] <- rgb(cur.color[1], cur.color[2], cur.color[3], alpha.channels[i])
	}
	#illumina.color <- c(	rgb(0.4, 0.4, 0.1, 0.0), rgb(0.4, 0.4, 0.1, 0.33), rgb(0.4, 0.4, 0.1, 1.0) ) #c("white", "gray50", "black") #3
	illumina.borders <- c(rgb(cur.color[1], cur.color[2], cur.color[3]), NA, NA)
	
	dgrc.y <- 3
	cur.color <- c(139, 28, 98)/255
	dgrc.color <- c(NA, NA, NA)
	for (i in 1:3) {
		dgrc.color[i] <- rgb(cur.color[1], cur.color[2], cur.color[3], alpha.channels[i])
	}
	#dgrc.color <- c(	rgb(0.4, 0.1, 0.4, 0.0), rgb(0.4, 0.1, 0.4, 0.33), rgb(0.4, 0.1, 0.4, 1.0) )  #c("white", "gray50", "black") #4
	dgrc.borders <- c(rgb(cur.color[1], cur.color[2], cur.color[3]), NA, NA)
	
	dge.y <- 2
	cur.color <- c(205, 102, 0)/255
	dge.color <- c(NA, NA, NA)
	for (i in 1:3) {
		dge.color[i] <- rgb(cur.color[1], cur.color[2], cur.color[3], alpha.channels[i])
	}
	#dge.color <- c(	rgb(0.1, 0.4, 0.4, 0.0), rgb(0.1, 0.4, 0.4, 0.33), rgb(0.1, 0.4, 0.4, 1.0) ) #c("white", "gray50", "black") #6
	dge.borders <- c(rgb(cur.color[1], cur.color[2], cur.color[3]), NA, NA)
	
	binding.y <- 5
	cur.color <- c(16, 78, 139)/255
	binding.color <- c(NA, NA, NA)
	for (i in 1:3) {
		binding.color[i] <- rgb(cur.color[1], cur.color[2], cur.color[3], alpha.channels[i])
	}
	#binding.color <- c(	rgb(0.1, 0.1, 0.9, 0.0), rgb(0.1, 0.1, 0.9, 0.33), rgb(0.1, 0.1, 0.9, 1.0) ) #c("white", "gray50", "black") #"orange"
	binding.borders <- c(rgb(cur.color[1], cur.color[2], cur.color[3]), NA, NA)
	
	
	#Set up the plot window
	par(mar=c(5,4,1,2)+0.1, cex=0.8)
	xlab <- paste("Position along chromosome", chrom)
	plot(x=NULL, y=NULL, xlim=xlim, ylim=c(0, 6), xlab=xlab, ylab="", yaxt="n")
	caption.x <- xlim[1] - (xlim[2]-xlim[1])/50
	
	#Plot the backcross results
	lines(bc.short.avg ~ pos, data=na.omit(cur.bc[,c("pos", "bc.short.avg")]), cex=0.25, col="grey70", lty=3)
	points(bc.short.avg ~ pos, data=cur.bc, cex=0.5)
	par(srt=90, cex=sidebar.cex)
	#text(x=caption.x, y=0.5, labels="Backcross")
	mtext(text="Backcross", side=2, at=0.5, line=1)
	
	#Now plot any significant deletions
	#In the zone 1:2
	thresh.src <- ifelse(is.null(thresh.src), thresh, thresh.src)
	cur.src.sig.i <- cur.src[(cur.src$q.interaction < thresh.src) & (cur.src$q.deletion >= thresh.src),]
	cur.src.sig.g <- cur.src[(cur.src$q.interaction >= thresh.src) & (cur.src$q.deletion < thresh.src),]
	cur.src.sig.b <- cur.src[(cur.src$q.interaction < thresh.src) & (cur.src$q.deletion < thresh.src),]
	for (i in 1:nrow(cur.src.sig.i)) {
		draw.bar(x=c(cur.src.sig.i$start[i], cur.src.sig.i$stop[i]), y=rep(src.y+cur.src.sig.i$y[i],2), col=src.color[2], border=src.borders[2])
	}
	for (i in 1:nrow(cur.src.sig.g)) {
		draw.bar(x=c(cur.src.sig.g$start[i], cur.src.sig.g$stop[i]), y=rep(src.y+cur.src.sig.g$y[i],2), col=src.color[1], border=src.borders[1])
	}
	for (i in 1:nrow(cur.src.sig.b)) {
		draw.bar(x=c(cur.src.sig.b$start[i], cur.src.sig.b$stop[i]), y=rep(src.y+cur.src.sig.b$y[i],2), col=src.color[3], border=src.borders[3])
	}
	mtext(text="Modifier\nDeletions", side=2, at=0.63+src.y, col=src.color, line=0.5)
	
	
	#Now plot any significant Illumina results
	thresh.illumina <- ifelse(is.null(thresh.illumina), thresh, thresh.illumina)
	if (which.illumina=="q") {
		cur.illumina.sig.genes.i <- cur.illumina$fbgn[(cur.illumina$genotype.wt.q >= thresh.illumina) & (cur.illumina$interaction.wt.sam.q < thresh.illumina)]
		cur.illumina.sig.genes.g <- cur.illumina$fbgn[(cur.illumina$genotype.wt.q < thresh.illumina) & (cur.illumina$interaction.wt.sam.q >= thresh.illumina)]
		cur.illumina.sig.genes.b <- cur.illumina$fbgn[(cur.illumina$genotype.wt.q < thresh.illumina) & (cur.illumina$interaction.wt.sam.q < thresh.illumina)]
	}
	if (which.illumina=="p") {
		cur.illumina.sig.genes.i <- cur.illumina$fbgn[(cur.illumina$genotype.wt.p >= thresh.illumina) & (cur.illumina$interaction.wt.sam.p < thresh.illumina)]
		cur.illumina.sig.genes.g <- cur.illumina$fbgn[(cur.illumina$genotype.wt.p < thresh.illumina) & (cur.illumina$interaction.wt.sam.p >= thresh.illumina)]
		cur.illumina.sig.genes.b <- cur.illumina$fbgn[(cur.illumina$genotype.wt.p < thresh.illumina) & (cur.illumina$interaction.wt.sam.p < thresh.illumina)]
	}
	#Do the interaction-effect-only genes
	cur.illumina.sig.info <- cur.flybase[cur.flybase$fbgn %in% cur.illumina.sig.genes.i,]
	for (i in 1:nrow(cur.illumina.sig.info)) {
		draw.bar(x=c(cur.illumina.sig.info$start[i], cur.illumina.sig.info$stop[i]), y=rep(illumina.y+cur.illumina.sig.info$y[i],2), col=illumina.color[2], border=illumina.borders[2])
		if (plot.gene.names & (length(cur.illumina.sig.info$symbol[i]) > 0)) {
			par(srt=0, cex=gene.cex)
			text(x=cur.illumina.sig.info$start[i], y=illumina.y+cur.illumina.sig.info$y[i]+gene.y, labels=cur.illumina.sig.info$symbol[i])
			par(srt=90, cex=sidebar.cex)
		}
	}
	#Do the genotype-effect-only genes
	cur.illumina.sig.info <- cur.flybase[cur.flybase$fbgn %in% cur.illumina.sig.genes.g,]
	for (i in 1:nrow(cur.illumina.sig.info)) {
		draw.bar(x=c(cur.illumina.sig.info$start[i], cur.illumina.sig.info$stop[i]), y=rep(illumina.y+cur.illumina.sig.info$y[i],2), col=illumina.color[1], border=illumina.borders[1])
		if (plot.gene.names & (length(cur.illumina.sig.info$symbol[i]) > 0)) {
			par(srt=0, cex=gene.cex)
			text(x=cur.illumina.sig.info$start[i], y=illumina.y+cur.illumina.sig.info$y[i]+gene.y, labels=cur.illumina.sig.info$symbol[i])
			par(srt=90, cex=sidebar.cex)
		}
	}
	#Do the genotype & interaction genes
	cur.illumina.sig.info <- cur.flybase[cur.flybase$fbgn %in% cur.illumina.sig.genes.b,]
	for (i in 1:nrow(cur.illumina.sig.info)) {
		draw.bar(x=c(cur.illumina.sig.info$start[i], cur.illumina.sig.info$stop[i]), y=rep(illumina.y+cur.illumina.sig.info$y[i],2), col=illumina.color[3], border=illumina.borders[3])
		if (plot.gene.names & (length(cur.illumina.sig.info$symbol[i]) > 0)) {
			par(srt=0, cex=gene.cex)
			text(x=cur.illumina.sig.info$start[i], y=illumina.y+cur.illumina.sig.info$y[i]+gene.y, labels=cur.illumina.sig.info$symbol[i])
			par(srt=90, cex=sidebar.cex)
		}
	}
	mtext(text="Illumina", side=2, at=0.5+illumina.y, col=illumina.color, line=1)
	
	#Same for DGRC
	thresh.dgrc <- ifelse(is.null(thresh.dgrc), thresh, thresh.dgrc)
	cur.dgrc.sig.genes.g <- cur.dgrc$fbgn[(cur.dgrc$interaction.wt.sam.q >= thresh.dgrc) & (cur.dgrc$genotype.wt.q < thresh.dgrc)]
	cur.dgrc.sig.genes.i <- cur.dgrc$fbgn[(cur.dgrc$interaction.wt.sam.q < thresh.dgrc) & (cur.dgrc$genotype.wt.q >= thresh.dgrc)]
	cur.dgrc.sig.genes.b <- cur.dgrc$fbgn[(cur.dgrc$interaction.wt.sam.q < thresh.dgrc) & (cur.dgrc$genotype.wt.q < thresh.dgrc)]
	#Do the genotype-effect-only genes
	cur.dgrc.sig.info <- cur.flybase[cur.flybase$fbgn %in% cur.dgrc.sig.genes.g,]
	for (i in 1:nrow(cur.dgrc.sig.info)) {
		draw.bar(x=c(cur.dgrc.sig.info$start[i], cur.dgrc.sig.info$stop[i]), y=rep(dgrc.y+cur.dgrc.sig.info$y[i],2), col=dgrc.color[1], border = dgrc.borders[1])
		if (plot.gene.names & (length(cur.dgrc.sig.info$symbol[i]) > 0)) {
			par(srt=0, cex=gene.cex)
			text(x=cur.dgrc.sig.info$start[i], y=dgrc.y+cur.dgrc.sig.info$y[i]+gene.y, labels=cur.dgrc.sig.info$symbol[i])
			par(srt=90, cex=sidebar.cex)
		}
	}
	#Do the interaction-effect-only genes
	cur.dgrc.sig.info <- cur.flybase[cur.flybase$fbgn %in% cur.dgrc.sig.genes.i,]
	for (i in 1:nrow(cur.dgrc.sig.info)) {
		draw.bar(x=c(cur.dgrc.sig.info$start[i], cur.dgrc.sig.info$stop[i]), y=rep(dgrc.y+cur.dgrc.sig.info$y[i],2), col=dgrc.color[2], border = dgrc.borders[2])
		if (plot.gene.names & (length(cur.dgrc.sig.info$symbol[i]) > 0)) {
			par(srt=0, cex=gene.cex)
			text(x=cur.dgrc.sig.info$start[i], y=dgrc.y+cur.dgrc.sig.info$y[i]+gene.y, labels=cur.dgrc.sig.info$symbol[i])
			par(srt=90, cex=sidebar.cex)
		}
	}
	#Do the genotype&interaction genes
	cur.dgrc.sig.info <- cur.flybase[cur.flybase$fbgn %in% cur.dgrc.sig.genes.b,]
	for (i in 1:nrow(cur.dgrc.sig.info)) {
		draw.bar(x=c(cur.dgrc.sig.info$start[i], cur.dgrc.sig.info$stop[i]), y=rep(dgrc.y+cur.dgrc.sig.info$y[i],2), col=dgrc.color[3], border=dgrc.borders[3])
		if (plot.gene.names & (length(cur.dgrc.sig.info$symbol[i]) > 0)) {
			par(srt=0, cex=gene.cex)
			text(x=cur.dgrc.sig.info$start[i], y=dgrc.y+cur.dgrc.sig.info$y[i]+gene.y, labels=cur.dgrc.sig.info$symbol[i])
			par(srt=90, cex=sidebar.cex)
		}
	}
	mtext(text="DGRC", side=2, at=0.5+dgrc.y, col=dgrc.color, line=1)
	
	
	#Same for DGE
	thresh.dge.genotype <- ifelse(is.null(thresh.dge.genotype), thresh, thresh.dge.genotype)
	thresh.dge.interaction <- ifelse(is.null(thresh.dge.interaction), thresh, thresh.dge.interaction)
	#Get the genes that have a significant genotype or interaction effect in the DGE dataset
	cur.dge.sig.genes.g.temp <- cur.dge.genotype$fbgn[(cur.dge.genotype$q < thresh.dge.genotype)]
	cur.dge.sig.genes.i.temp <- cur.dge$fbgn[(cur.dge$q.2 < thresh.dge.interaction)]
	#Now get the genes that have a significant effect for both genotype and interaction
	cur.dge.sig.genes.b <- cur.dge.sig.genes.i.temp[cur.dge.sig.genes.i.temp %in% cur.dge.sig.genes.g.temp]
	#Now get the genes that have a significant effect for only one or the other
	cur.dge.sig.genes.g <- cur.dge.sig.genes.g.temp[!(cur.dge.sig.genes.g.temp %in% cur.dge.sig.genes.b)]
	cur.dge.sig.genes.i <- cur.dge.sig.genes.i.temp[!(cur.dge.sig.genes.i.temp %in% cur.dge.sig.genes.b)]
	#Do the genotype-effect-only genes
	cur.dge.sig.info <- cur.flybase[cur.flybase$fbgn %in% cur.dge.sig.genes.g,]
	for (i in 1:nrow(cur.dge.sig.info)) {
		draw.bar(x=c(cur.dge.sig.info$start[i], cur.dge.sig.info$stop[i]), y=rep(dge.y+cur.dge.sig.info$y[i],2), col=dge.color[1], border=dge.borders[1])
		if (plot.gene.names & (length(cur.dge.sig.info$symbol[i]) > 0)) {
			par(srt=0, cex=gene.cex)
			text(x=cur.dge.sig.info$start[i], y=dge.y+cur.dge.sig.info$y[i]+gene.y, labels=cur.dge.sig.info$symbol[i])
			par(srt=90, cex=sidebar.cex)
		}
	}
	#Do the interaction-effect-only genes
	cur.dge.sig.info <- cur.flybase[cur.flybase$fbgn %in% cur.dge.sig.genes.i,]
	for (i in 1:nrow(cur.dge.sig.info)) {
		draw.bar(x=c(cur.dge.sig.info$start[i], cur.dge.sig.info$stop[i]), y=rep(dge.y+cur.dge.sig.info$y[i],2), col=dge.color[2], border=dge.borders[2])
		if (plot.gene.names & (length(cur.dge.sig.info$symbol[i]) > 0)) {
			par(srt=0, cex=gene.cex)
			text(x=cur.dge.sig.info$start[i], y=dge.y+cur.dge.sig.info$y[i]+gene.y, labels=cur.dge.sig.info$symbol[i])
			par(srt=90, cex=sidebar.cex)
		}
	}
	#Do the genotype&interaction genes
	cur.dge.sig.info <- cur.flybase[cur.flybase$fbgn %in% cur.dge.sig.genes.b,]
	for (i in 1:nrow(cur.dge.sig.info)) {
		draw.bar(x=c(cur.dge.sig.info$start[i], cur.dge.sig.info$stop[i]), y=rep(dge.y+cur.dge.sig.info$y[i],2), col=dge.color[3], border=dge.borders[3])
		if (plot.gene.names & (length(cur.dge.sig.info$symbol[i]) > 0)) {
			par(srt=0, cex=gene.cex)
			text(x=cur.dge.sig.info$start[i], y=dge.y+cur.dge.sig.info$y[i]+gene.y, labels=cur.dge.sig.info$symbol[i])
			par(srt=90, cex=sidebar.cex)
		}
	}
	mtext(text="DGE", side=2, at=0.5+dge.y, col=dge.color, line=1)

	#Now plot the binding site results
	thresh.binding <- ifelse(is.null(thresh.binding), thresh, thresh.binding)
	thresh.binding <- ifelse(is.null(thresh.binding), thresh, thresh.binding)
	#Get the genotype-effect genes	
		#Count a gene as significant for the overall main effect if it meets the threshold for both SAM and ORE,
		# or if it meets the threshold for one and is close in the other (e.g., we want to keep the genes that have
		# p-values of 0.049 and 0.059, whereas these would be discarded if we required them both to meet the threshold)
		sig.sam <- cur.binding$empirical.p.sam < thresh.binding
		sig.ore <- cur.binding$empirical.p.ore < thresh.binding
		similar.p <- (cur.binding$empirical.p.sam - cur.binding$empirical.p.ore) < 0.01
		cur.binding.sig.genes.g.temp <- cur.binding$fbgn[(sig.sam & sig.ore) | (sig.sam & similar.p) | (sig.ore & similar.p)]
	
	#Get the interaction-effect genes
	cur.binding.sig.genes.i.temp <- cur.binding$fbgn[cur.binding$empirical.p < thresh.binding]
	
	#Now get the genes that are in both...
	cur.binding.sig.genes.b <- cur.binding$fbgn[(cur.binding$fbgn %in% cur.binding.sig.genes.g.temp) & (cur.binding$fbgn %in% cur.binding.sig.genes.i.temp)]
	
	#Now get the sets of genes that ONLY have a genotype effect or ONLY have an interaction effect
	cur.binding.sig.genes.g <- cur.binding.sig.genes.g.temp[!(cur.binding.sig.genes.g.temp %in% cur.binding.sig.genes.b)]
	cur.binding.sig.genes.i <- cur.binding.sig.genes.i.temp[!(cur.binding.sig.genes.i.temp %in% cur.binding.sig.genes.b)]

	#Do the genotye-effect-only genes	
	cur.binding.sig.info <- cur.flybase[cur.flybase$fbgn %in% cur.binding.sig.genes.g,]
	for (i in 1:nrow(cur.binding.sig.info)) {
		draw.bar(x=c(cur.binding.sig.info$start[i], cur.binding.sig.info$stop[i]), y=rep(binding.y+cur.binding.sig.info$y[i],2), col=binding.color[1], border=binding.borders[1])
		if (plot.gene.names & (length(cur.binding.sig.info$symbol[i]) > 0)) {
			par(srt=0, cex=gene.cex)
			text(x=cur.binding.sig.info$start[i], y=binding.y+cur.binding.sig.info$y[i]+gene.y, labels=cur.binding.sig.info$symbol[i])
			par(srt=90, cex=sidebar.cex)
		}
	}
	#Do the interaction-effect-only genes	
	cur.binding.sig.info <- cur.flybase[cur.flybase$fbgn %in% cur.binding.sig.genes.i,]
	for (i in 1:nrow(cur.binding.sig.info)) {
		draw.bar(x=c(cur.binding.sig.info$start[i], cur.binding.sig.info$stop[i]), y=rep(binding.y+cur.binding.sig.info$y[i],2), col=binding.color[2], border=binding.borders[2])
		if (plot.gene.names & (length(cur.binding.sig.info$symbol[i]) > 0)) {
			par(srt=0, cex=gene.cex)
			text(x=cur.binding.sig.info$start[i], y=binding.y+cur.binding.sig.info$y[i]+gene.y, labels=cur.binding.sig.info$symbol[i])
			par(srt=90, cex=sidebar.cex)
		}
	}
	#Do the genotype-and-interaction genes	
	cur.binding.sig.info <- cur.flybase[cur.flybase$fbgn %in% cur.binding.sig.genes.b,]
	for (i in 1:nrow(cur.binding.sig.info)) {
		draw.bar(x=c(cur.binding.sig.info$start[i], cur.binding.sig.info$stop[i]), y=rep(binding.y+cur.binding.sig.info$y[i],2), col=binding.color[3], border=binding.borders[3])
		if (plot.gene.names & (length(cur.binding.sig.info$symbol[i]) > 0)) {
			par(srt=0, cex=gene.cex)
			text(x=cur.binding.sig.info$start[i], y=binding.y+cur.binding.sig.info$y[i]+gene.y, labels=cur.binding.sig.info$symbol[i])
			par(srt=90, cex=sidebar.cex)
		}
	}
	mtext(text="Binding\nPredictions", side=2, at=0.63+binding.y, col=binding.color, line=0.5)
	
}

##############################################################################
# Analysis
##############################################################################

thresh.illumina <- 0.05
which.illumina <- "p" #set to "p" or "q" to switch between p- and q-values for Illumina analysis significance cutoff
thresh.dgrc <- 0.05
thresh.src <- 0.05
thresh.dge.interaction <- 0.05
thresh.dge.genotype <- 0.05
thresh.binding <- 0.05
pdf.width <- 12
pdf.height <- 6
min.sig.tally <- 4

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#
#   Assemble all the results into a combined dataset
#
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################


#Returns information about the datasets in which a gene showed a significant signal
get.gene.results <- function(i) {
	fbgn <- flybase.data$fbgn[i]
	cur.chrom <- flybase.data$chrom[i]
	cur.start <- flybase.data$start[i]
	cur.stop <- flybase.data$stop[i]
	
	#Information to return
	sig.tally <- 0
	sig.expression.tally <- 0
	sig.genotype.effect.tally <- 0
	sig.interaction.effect.tally <- 0
	dge.genotype.q <- NA
	dge.interaction.q <- NA
	illumina.genotype.pq <- NA
	illumina.interaction.pq <- NA
	dgrc.genotype.q <- NA
	dgrc.interaction.q <- NA
	min.deletion.q <- NA
	min.deletion.interaction.q <- NA
	binding.p <- NA
	binding.interaction.p <- NA
	max.bcshort.freq <- NA
	inconsistent.expression.effects <- NA
	
	genotype.effect.signs <- NULL
	
	#DGE (genotype)
	cur.dge.genotype.results <- dge.genotype.results[dge.genotype.results$fbgn==fbgn,]
	cur.dge.genotype.results <- cur.dge.genotype.results[complete.cases(cur.dge.genotype.results),]
	if (nrow(cur.dge.genotype.results) > 0) {
		dge.genotype.q <- min(na.omit(cur.dge.genotype.results$q))
		min.q.index <- which.min(cur.dge.genotype.results$q)
		if (dge.genotype.q < thresh.dge.genotype) {
			sig.tally <- sig.tally + 1
			sig.expression.tally <- sig.expression.tally + 1
			sig.genotype.effect.tally <- sig.genotype.effect.tally + 1
			cur.dge.genotype.effect <- (cur.dge.genotype.results$fc5.count + cur.dge.genotype.results$fc6.count - cur.dge.genotype.results$fc3.count - cur.dge.genotype.results$fc4.count)[min.q.index]
			genotype.effect.signs <- c(genotype.effect.signs, sign(cur.dge.genotype.effect))
		}
	}
	
	#DGE (interaction)
	cur.dge.interaction.results <- dge.results[dge.results$fbgn==fbgn,]
	if (nrow(cur.dge.interaction.results) > 0) {
		dge.interaction.q <- min(na.omit(cur.dge.interaction.results$q.2))
		if (dge.interaction.q < thresh.dge.interaction) {
			sig.tally <- sig.tally + 1
			sig.expression.tally <- sig.expression.tally + 1
			sig.interaction.effect.tally <- sig.interaction.effect.tally + 1
		}
	}
	
	#Illumina (genotype & interaction)
	cur.illumina.results <- illumina.results[illumina.results$fbgn==fbgn,]
	cur.illumina.results <- cur.illumina.results[complete.cases(cur.illumina.results),]
	if (nrow(cur.illumina.results) > 0) {
		
		#What we do next depends on whether we are set to use p-values or q-values for the Illumina dataset
		if (which.illumina=="q") {
			illumina.genotype.pq <- min(na.omit(cur.illumina.results$genotype.wt.q))
			if (illumina.genotype.pq < thresh.illumina) {
				sig.tally <- sig.tally + 1
				sig.expression.tally <- sig.expression.tally + 1
				sig.genotype.effect.tally <- sig.genotype.effect.tally + 1
				cur.illumina.effect <- cur.illumina.results$genotype.wt.effect[which.min(cur.illumina.results$genotype.wt.q)]
				genotype.effect.signs <- c(genotype.effect.signs, sign(cur.illumina.effect))
			}
			illumina.interaction.pq <- min(na.omit(cur.illumina.results$interaction.wt.sam.q))
			if (illumina.interaction.pq < thresh.illumina) {
				sig.tally <- sig.tally + 1
				sig.expression.tally <- sig.expression.tally + 1
				sig.interaction.effect.tally <- sig.interaction.effect.tally + 1
			}
		}
		if (which.illumina=="p") {
			illumina.genotype.pq <- min(na.omit(cur.illumina.results$genotype.wt.p))
			if (illumina.genotype.pq < thresh.illumina) {
				sig.tally <- sig.tally + 1
				sig.expression.tally <- sig.expression.tally + 1
				sig.genotype.effect.tally <- sig.genotype.effect.tally + 1
				cur.illumina.effect <- cur.illumina.results$genotype.wt.effect[which.min(cur.illumina.results$genotype.wt.p)]
				genotype.effect.signs <- c(genotype.effect.signs, sign(cur.illumina.effect))
			}
			illumina.interaction.pq <- min(na.omit(cur.illumina.results$interaction.wt.sam.p))
			if (illumina.interaction.pq < thresh.illumina) {
				sig.tally <- sig.tally + 1
				sig.expression.tally <- sig.expression.tally + 1
				sig.interaction.effect.tally <- sig.interaction.effect.tally + 1
			}
		}

	}
	
	#DGRC (genotype & interaction)
	cur.dgrc.results <- dgrc.results[dgrc.results$fbgn==fbgn,]
	cur.dgrc.results <- cur.dgrc.results[complete.cases(cur.dgrc.results),]
	if (nrow(cur.dgrc.results) > 0) {
		dgrc.genotype.q <- min(na.omit(cur.dgrc.results$genotype.wt.q))
		if (dgrc.genotype.q < thresh.dgrc) {
			sig.tally <- sig.tally + 1
			sig.expression.tally <- sig.expression.tally + 1
			sig.genotype.effect.tally <- sig.genotype.effect.tally + 1
			cur.dgrc.effect <- cur.dgrc.results$genotype.wt.effect[which.min(cur.dgrc.results$genotype.wt.q)]
			genotype.effect.signs <- c(genotype.effect.signs, sign(cur.dgrc.effect))
		}
		dgrc.interaction.q <- min(na.omit(cur.dgrc.results$interaction.wt.sam.q))
		if (dgrc.interaction.q < thresh.dgrc) {
			sig.tally <- sig.tally + 1
			sig.expression.tally <- sig.expression.tally + 1
			sig.interaction.effect.tally <- sig.interaction.effect.tally + 1
		}
	}
	
	#Now use the genotype effect signs to flag genes that are inconsistent
	genotype.effect.signs <- na.omit(genotype.effect.signs)
	if (length(genotype.effect.signs > 1)) {
		inconsistent.expression.effects <- !(all(genotype.effect.signs == genotype.effect.signs[1]))
	}
	
	#Modifiers
	cur.src.results <- src.results[src.results$chrom==cur.chrom,]
	keeper.conditions <- ((cur.src.results$start <= cur.start) & (cur.src.results$stop >= cur.stop)) |
	                     ((cur.start <= cur.src.results$start) & (cur.stop >= cur.src.results$stop)) |
	                     ((cur.start >= cur.src.results$start) & (cur.start <= cur.src.results$stop)) |
	                     ((cur.stop >= cur.src.results$start)  & (cur.stop <= cur.src.results$stop))
	cur.src.results <- cur.src.results[keeper.conditions,]
	if (nrow(cur.src.results) > 0) {
		min.deletion.q <- min(na.omit(cur.src.results$q.deletion))
		min.deletion.interaction.q <- min(na.omit(cur.src.results$q.interaction))
		if (min.deletion.q < thresh.src) {
			sig.tally <- sig.tally + 1
			sig.genotype.effect.tally <- sig.genotype.effect.tally + 1
		}
		if (min.deletion.interaction.q < thresh.src) {
			sig.tally <- sig.tally + 1
			sig.interaction.effect.tally <- sig.interaction.effect.tally + 1	
		} 
	}
	
	#Binding sites -- only for putative *polymorphic* targets
	#Here, we will use p-values rather than q-values, since none of the q-values are below ~0.9
	cur.binding.results <- binding.results[binding.results$fbgn==fbgn,]
	if (nrow(cur.binding.results) > 0) {
		binding.p <- min(na.omit(c(cur.binding.results$empirical.p.sam, cur.binding.results$empirical.p.ore)))
		binding.interaction.p <- min(na.omit(cur.binding.results$empirical.p))
		if (binding.p < thresh.binding) sig.tally <- sig.tally + 1
		if (binding.interaction.p < thresh.binding) sig.tally <- sig.tally + 1
	}
	
	#Mapping	
	cur.bc.results <- bc.results[bc.results$chrom==cur.chrom,]
	cur.bc.results$start <- cur.bc.results$pos
	cur.bc.results$stop <- cur.bc.results$pos + 1999
	keeper.conditions <-     ((cur.bc.results$start <= cur.start) & (cur.bc.results$stop >= cur.stop)) |
	                         ((cur.start <= cur.bc.results$start) & (cur.stop >= cur.bc.results$stop)) |
	                         ((cur.start >= cur.bc.results$start) & (cur.start <= cur.bc.results$stop)) |
	                         ((cur.stop >= cur.bc.results$start)  & (cur.stop <= cur.bc.results$stop))
	cur.bc.results <- cur.bc.results[keeper.conditions,]
	if (nrow(cur.bc.results) > 0) {
		max.bcshort.freq <- max(na.omit(cur.bc.results$bc.short.avg))
		if (!is.finite(max.bcshort.freq)) max.bcshort.freq <- NA
		#Count it as "significant" if it has an average frequency of at least 0.7
		if (any(na.omit(cur.bc.results$bc.short.avg) > 0.7)) sig.tally <- sig.tally + 1
	}

	return(c(sig.tally, sig.expression.tally, sig.genotype.effect.tally, sig.interaction.effect.tally, dge.genotype.q, dge.interaction.q, illumina.genotype.pq, illumina.interaction.pq, dgrc.genotype.q, dgrc.interaction.q, min.deletion.q, min.deletion.interaction.q, binding.p, binding.interaction.p, max.bcshort.freq, inconsistent.expression.effects))
}

flybase.data$sig.tally <- rep(NA, nrow(flybase.data))
flybase.data$sig.expression.tally <- rep(NA, nrow(flybase.data))
flybase.data$sig.genotype.effect.tally <- rep(NA, nrow(flybase.data))
flybase.data$sig.interaction.effect.tally <- rep(NA, nrow(flybase.data))
flybase.data$dge.genotype.q <- rep(NA, nrow(flybase.data))
flybase.data$dge.interaction.q <- rep(NA, nrow(flybase.data))
flybase.data$illumina.genotype.pq <- rep(NA, nrow(flybase.data))
flybase.data$illumina.interaction.pq <- rep(NA, nrow(flybase.data))
flybase.data$dgrc.genotype.q <- rep(NA, nrow(flybase.data))
flybase.data$dgrc.interaction.q <- rep(NA, nrow(flybase.data))
flybase.data$min.deletion.q <- rep(NA, nrow(flybase.data))
flybase.data$min.deletion.interaction.q <- rep(NA, nrow(flybase.data))
flybase.data$binding.p <- rep(NA, nrow(flybase.data))
flybase.data$binding.interaction.p <- rep(NA, nrow(flybase.data))
flybase.data$max.bcshort.freq <- rep(NA, nrow(flybase.data))
flybase.data$inconsistent.expression.effects <- rep(NA, nrow(flybase.data))
for (i in 1:nrow(flybase.data)) {
	flybase.data[i, c("sig.tally", "sig.expression.tally", "sig.genotype.effect.tally", "sig.interaction.effect.tally", "dge.genotype.q", "dge.interaction.q", "illumina.genotype.pq", "illumina.interaction.pq", "dgrc.genotype.q", "dgrc.interaction.q", "min.deletion.q", "min.deletion.interaction.q", "binding.p", "binding.interaction.p", "max.bcshort.freq", "inconsistent.expression.effects")] <- get.gene.results(i)
	if (i %% 100 == 0) print(i)
}

sorted.flybase.data <- flybase.data[order(flybase.data$sig.tally, decreasing=TRUE),]

#We can save the results so we don't have to re-compute them every time we run the script if we want
if (which.illumina == "p") write.csv(sorted.flybase.data, file="integrated_results_p.csv", row.names=FALSE)
if (which.illumina == "q") write.csv(sorted.flybase.data, file="integrated_results_q.csv", row.names=FALSE)

#Now we can restrict plotting to only those that have a significant data source tally of at least 4

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#
#   Do the plotting
#
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################


plot.chrom.pdf <- function(filename, width, height, thresh, chrom, xlim=NULL, plot.gene.names=FALSE, thresh.dge.genotype=NULL, thresh.dge.interaction=NULL, thresh.dgrc=NULL, thresh.illumina=NULL, thresh.src=NULL, thresh.binding=NULL) {
	pdf(filename, width=width, height=height)
	plot.region(chrom, thresh=thresh, xlim=xlim, plot.gene.names=plot.gene.names, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction, thresh.dgrc=thresh.dgrc, thresh.illumina=thresh.illumina, thresh.src=thresh.src, thresh.binding=thresh.binding)
	dev.off()
}



	
	#Make a PDF for each chromosome
	chroms <- c("X", "2L", "2R", "3L", "3R")
	for (c in chroms) {
		cur.file <- paste("output_withbinding_interaction_plus_genotype/chrom_", c, "_", ".pdf", sep="")
		plot.chrom.pdf(cur.file, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding, chrom=c)
	}
	
	#Zoom in on the following regions...
	cur.file <- paste("output_withbinding_interaction_plus_genotype/2L_zoom_1_", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="2L", xlim=c(4.5e6, 1.1e7), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	cur.file <- paste("output_withbinding_interaction_plus_genotype/2L_zoom_1a_", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="2L", xlim=c(6e6, 6.5e6), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	cur.file <- paste("output_withbinding_interaction_plus_genotype/2L_zoom_1b_", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="2L", xlim=c(8.6e6, 9.5e6), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	cur.file <- paste("output_withbinding_interaction_plus_genotype/2L_zoom_1c_", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="2L", xlim=c(9.7e6, 1.05e7), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	
	cur.file <- paste("output_withbinding_interaction_plus_genotype/2R_zoom_1_", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="2R", xlim=c(6e6, 9e6), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	
	cur.file <- paste("output_withbinding_interaction_plus_genotype/3L_zoom_1_", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="3L", xlim=c(0, 1.05e6), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	cur.file <- paste("output_withbinding_interaction_plus_genotype/3L_zoom_2_", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="3L", xlim=c(1.5e6, 2.2e6), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	cur.file <- paste("output_withbinding_interaction_plus_genotype/3L_zoom_3_", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="3L", xlim=c(4.0e6, 6.0e6), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	cur.file <- paste("output_withbinding_interaction_plus_genotype/3L_zoom_4_", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="3L", xlim=c(1.7e7, 1.95e7), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	
	cur.file <- paste("output_withbinding_interaction_plus_genotype/3R_zoom_1_", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="3R", xlim=c(5e6, 1.1e7), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	cur.file <- paste("output_withbinding_interaction_plus_genotype/3R_zoom_1a_", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="3R", xlim=c(5.2e6, 6e6), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	cur.file <- paste("output_withbinding_interaction_plus_genotype/3R_zoom_1b_", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="3R", xlim=c(9.6e6, 1.04e7), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	cur.file <- paste("output_withbinding_interaction_plus_genotype/3R_zoom_2_", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="3R", xlim=c(1.15e7, 1.43e7), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	cur.file <- paste("output_withbinding_interaction_plus_genotype/3R_zoom_3_", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="3R", xlim=c(2.4e7, 2.5e7), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	
	cur.file <- paste("output_withbinding_interaction_plus_genotype/X_zoom_1_", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="X", xlim=c(0, 1e6), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	cur.file <- paste("output_withbinding_interaction_plus_genotype/X_zoom_2_", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="X", xlim=c(5.5e6, 6.5e6), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)

	#########
	########
	#######

	cur.file <- paste("output_custom_zooms/2L_zoom1", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="2L", xlim=c(6150000, 7355000), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	
	cur.file <- paste("output_custom_zooms/2L_zoom2", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="2L", xlim=c(9880000, 10050000), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	
	cur.file <- paste("output_custom_zooms/2L_zoom3", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="2L", xlim=c(18685000, 18715000), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	
	cur.file <- paste("output_custom_zooms/2R_zoom1", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="2R", xlim=c(4018000, 7395000), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	
	cur.file <- paste("output_custom_zooms/2R_zoom2", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="2R", xlim=c(8771700, 10029000), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)

	cur.file <- paste("output_custom_zooms/2R_zoom3", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="2R", xlim=c(14140000, 14256000), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	
	cur.file <- paste("output_custom_zooms/2R_zoom4", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="2R", xlim=c(20690000, 20735000), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	
	cur.file <- paste("output_custom_zooms/3L_zoom1", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="3L", xlim=c(1486000, 1515000), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	
	cur.file <- paste("output_custom_zooms/3L_zoom2", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="3L", xlim=c(8099000, 8125000), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	
	cur.file <- paste("output_custom_zooms/3L_zoom3", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="3L", xlim=c(16570000, 16578000), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	
	cur.file <- paste("output_custom_zooms/3R_zoom1", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="3R", xlim=c(1665000, 1693000), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	
	cur.file <- paste("output_custom_zooms/3R_zoom2", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="3R", xlim=c(5383000, 7300000), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	
	cur.file <- paste("output_custom_zooms/3R_zoom3", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="3R", xlim=c(17646000, 17892000), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
	
	cur.file <- paste("output_custom_zooms/3R_zoom4", ".pdf", sep="")
	plot.chrom.pdf(cur.file, chrom="3R", xlim=c(22810000, 25610000), plot.gene.names=TRUE, width=pdf.width, height=pdf.height, thresh.illumina=thresh.illumina, thresh.dgrc=thresh.dgrc, thresh.src=thresh.src, thresh.dge.genotype=thresh.dge.genotype, thresh.dge.interaction=thresh.dge.interaction,thresh.binding=thresh.binding)
