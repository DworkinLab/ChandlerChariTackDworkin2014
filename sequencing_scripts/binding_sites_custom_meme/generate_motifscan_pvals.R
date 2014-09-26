#Read and prepare the SAM and ORE binding site results
col.names <- c("discard", "gene", "llr", "motifweight")
col.classes <- c("character", "character", "numeric", "numeric")
binding.results.sam <- read.table("output_sam_custom_meme/sam_binding_seqs.fasta.fen", col.names=col.names, colClasses=col.classes)[,2:4]
binding.results.ore <- read.table("output_ore_custom_meme/ore_binding_seqs.fasta.fen", col.names=col.names, colClasses=col.classes)[,2:4]

#Make sure the genes line up...
if (all(binding.results.sam$gene==binding.results.ore$gene)) {
	binding.results <- data.frame(gene=binding.results.sam$gene, llr.sam=binding.results.sam$llr, llr.ore=binding.results.ore$llr)
	binding.results$llr.diff <- abs(binding.results$llr.sam - binding.results$llr.ore)
} else {
	print("ERROR: SAM and ORE rows do not match up")
}


#Read in the GC-content stuff
gc.col.names <- c("gene", "gc", "size")
gc.data <- read.table("sam_gc.txt", col.names=gc.col.names) #Just use SAM's GC-content; differences from ORE are negligible for this analysis

if (all(as.character(gc.data$gene) == as.character(binding.results$gene))) {
	binding.results$gc <- gc.data$gc
	binding.results$size <- gc.data$size
	#binding.results <- data.frame(gene=binding.results$gene, llr=binding.results$llr, gc=gc.data$gc)
} else {
	print("ERROR: Cannot merge GC content with binding results")
}


#Calculate empirical p-values
num.to.compare <- 1000

binding.results$empirical.p <- rep(NA, nrow(binding.results))
binding.results$empirical.p.sam <- rep(NA, nrow(binding.results))
binding.results$empirical.p.ore <- rep(NA, nrow(binding.results))
resorted.results <- binding.results[order(binding.results$gc),]
for (cur.row in 1:nrow(binding.results)) {
	if (cur.row %% 500 == 0) {
		print(paste("On row: ", cur.row))
	}
	
	cur.gene <- as.character(binding.results$gene[cur.row])
	obs.llr.diff <- binding.results$llr.diff[cur.row]
	obs.llr.sam <- binding.results$llr.sam[cur.row]
	obs.llr.ore <- binding.results$llr.ore[cur.row]
	
	#Which row in resorted.results corresponds to the current one in binding.results?
	sorted.index <- which(as.character(resorted.results$gene)==cur.gene)

	#Which rows represent genes with the most similar gene content to the current one?
	similar.min <- sorted.index - round(num.to.compare/2)
	similar.max <- sorted.index + round(num.to.compare/2) - 1
	if (similar.min < 1) {
		similar.min = 1
		similar.max = 1000
	}
	if (similar.max > nrow(resorted.results)) {
		similar.max <- nrow(resorted.results)
		similar.min <- similar.max - 999
	}
	similar.results <- resorted.results[similar.min:similar.max,]

	empirical.p <- sum(similar.results$llr.diff >= obs.llr.diff) / num.to.compare
	empirical.p.sam <- sum(similar.results$llr.sam >= obs.llr.sam) / num.to.compare
	empirical.p.ore <- sum(similar.results$llr.ore >= obs.llr.ore) / num.to.compare

	##Now rank them by LLR so that we can empirically estimate the p-value
	#similar.results <- similar.results[order(similar.results$llr.diff),]
	#found.index <- which(as.character(similar.results$gene) == cur.gene)
	#empirical.p <- 1 - (found.index / nrow(similar.results))
	
	binding.results$empirical.p[cur.row] <- empirical.p
	binding.results$empirical.p.sam[cur.row] <- empirical.p.sam
	binding.results$empirical.p.ore[cur.row] <- empirical.p.ore
}

#Add q-values
library(qvalue)

binding.results$empirical.q.sam <- rep(NA, nrow(binding.results))
full.rows <- which((!is.na(binding.results$empirical.p.sam)) & (is.finite(binding.results$empirical.p.sam)))
binding.results$empirical.q.sam[full.rows] <- qvalue(binding.results$empirical.p.sam[full.rows])$qvalues

binding.results$empirical.q.ore <- rep(NA, nrow(binding.results))
full.rows <- which((!is.na(binding.results$empirical.p.ore)) & (is.finite(binding.results$empirical.p.ore)))
binding.results$empirical.q.ore[full.rows] <- qvalue(binding.results$empirical.p.ore[full.rows])$qvalues

binding.results$empirical.q <- rep(NA, nrow(binding.results))
full.rows <- which((!is.na(binding.results$empirical.p)) & (is.finite(binding.results$empirical.p)))
binding.results$empirical.q[full.rows] <- qvalue(binding.results$empirical.p[full.rows])$qvalues

write.csv(binding.results, file="SAM_ORE_binding_results_bygene_custom_meme.csv", row.names=FALSE, quote=FALSE)

#NOTE: In the output file,
#	empirical.p is the empirical p-value that the gene is predicted to be differentially regulated by SAM & ORE
#	 (i.e., difference between SAM & ORE binding scores is at the tail of the distribution)
#	empirical.p.sam & empirical.p.ore are the empricial p-values that the gene is predicted to be a scalloped
#	 binding target in each

#Double-binding site dataset definitely has bigger differences, notably, many genes where one is >= 3 and the other is not...

binding.results <- read.csv(file="SAM_ORE_binding_results_bygene_custom_meme.csv")
sig.thresh <- 0.01
#How many predicted targets in SAM?
sprintf("Number of predicted binding targets in SAM: %d", sum(binding.results$empirical.p.sam <= sig.thresh))
#How many predicted targets in ORE?
sprintf("Number of predicted binding targets in ORE: %d", sum(binding.results$empirical.p.ore <= sig.thresh))
#How many predicted to be targets in both?
sprintf("Number of predicted binding targets in both: %d", sum((binding.results$empirical.p.sam <= sig.thresh) & (binding.results$empirical.p.ore <= sig.thresh)))
#How many predicted to be targets in either?
sprintf("Number of predicted binding targets in either: %d", sum((binding.results$empirical.p.sam <= sig.thresh) | (binding.results$empirical.p.ore <= sig.thresh)))
sprintf("Number of predicted binding targets in SAM only: %d", sum((binding.results$empirical.p.sam <= sig.thresh) & !(binding.results$empirical.p.ore <= sig.thresh)))
sprintf("Number of predicted binding targets in ORE only: %d", sum(!(binding.results$empirical.p.sam <= sig.thresh) & (binding.results$empirical.p.ore <= sig.thresh)))


#How many are predicted to be differentially regulated in SAM & ORE -- total?
sprintf("Number of predicted differential targets: %d", sum(binding.results$empirical.p <= sig.thresh))
#Now, of those that are predicted to be significant targets in either SAM or ORE, how many are predicted to be differential targets?
sig.targets <- binding.results[(binding.results$empirical.p.sam <= sig.thresh) | (binding.results$empirical.p.ore <= sig.thresh),]
sprintf("Number of predicted differential targets, considering only those predicted to be a target: %d", sum(sig.targets$empirical.p <= sig.thresh))
