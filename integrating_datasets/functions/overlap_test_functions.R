#This function does the randomization test, shuffling over ALL genes in each dataset, rather than
# limiting to the analysis to the genes that are already common to both datasets
overlap.test.all <- function(data.1.orig, data.2.orig, p.col.1, p.col.2, thresh.1, thresh.2, n.reps=1000) {
	data.1 <- na.omit(data.1.orig[,c("fbgn", p.col.1)])
	data.2 <- na.omit(data.2.orig[,c("fbgn", p.col.2)])
	data.1.sig <- unique(data.1$fbgn[data.1[,2] < thresh.1])
	data.2.sig <- unique(data.2$fbgn[data.2[,2] < thresh.2])
	obs.overlap <- sum(data.1.sig %in% data.2.sig)
	sim.overlap <- rep(NA, n.reps)
	sim.overlap[1] <- obs.overlap
	for (i in 2:n.reps) {
		data.1$fbgn <- sample(data.1$fbgn, nrow(data.1), replace=FALSE)
		data.2$fbgn <- sample(data.2$fbgn, nrow(data.2), replace=FALSE)
		data.1.sig <- unique(data.1$fbgn[data.1[,2] < thresh.1])
		data.2.sig <- unique(data.2$fbgn[data.2[,2] < thresh.2])
		sim.overlap[i] <- sum(data.1.sig %in% data.2.sig)
	}
	p <- sum(sim.overlap >= obs.overlap) / length(sim.overlap)
	return(p)
}

#This function does the randomization test, but first trims each dataset down to the genes that are present in both
overlap.test.trimmed <- function(data.1.orig, data.2.orig, p.col.1, p.col.2, thresh.1, thresh.2, n.reps=1000) {
	data.1 <- na.omit(data.1.orig[,c("fbgn", p.col.1)])
	data.2 <- na.omit(data.2.orig[,c("fbgn", p.col.2)])
	data.1 <- data.1[data.1$fbgn %in% data.2$fbgn,]
	data.2 <- data.2[data.2$fbgn %in% data.1$fbgn,]
	data.1.sig <- unique(data.1$fbgn[data.1[,2] < thresh.1])
	n.data.1.sig <- length(data.1.sig)
	data.2.sig <- unique(data.2$fbgn[data.2[,2] < thresh.2])
	n.data.2.sig <- length(data.2.sig)
	obs.overlap <- sum(data.1.sig %in% data.2.sig)
	sim.overlap <- rep(NA, n.reps)
	sim.overlap[1] <- obs.overlap
	for (i in 2:n.reps) {
		data.1$fbgn <- sample(data.1$fbgn, nrow(data.1), replace=FALSE)
		data.2$fbgn <- sample(data.2$fbgn, nrow(data.2), replace=FALSE)
		data.1.sig <- unique(data.1$fbgn[data.1[,2] <= thresh.1])
		data.2.sig <- unique(data.2$fbgn[data.2[,2] <= thresh.2])
		sim.overlap[i] <- sum(data.1.sig %in% data.2.sig)
	}
	p <- sum(sim.overlap >= obs.overlap) / length(sim.overlap)
	
	print("Number of significant genes in dataset 1")
	print(n.data.1.sig)
	print("Number of significant genes in dataset 2")
	print(n.data.2.sig) 
	print("Number of genes common to both:")
	n.common.genes <- length(unique(data.1$fbgn))
	print(n.common.genes)
	print("Observed overlap")
	print(obs.overlap)
	
	return(p)
}