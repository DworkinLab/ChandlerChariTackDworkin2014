#####################################################################################
#####################################################################################
# Functions

#Takes a raw microarray dataset (each column representing an array)
#Does a log base 2 transformation, and then performs a median normalization
median_normalize <- function(dataset, cols=3:18) {
	new.data <- dataset

	#First do a log-base-2 transformation
	new.data[,cols] <- log2(new.data[,cols])
	
	#Now for each column, subtract the median from each cell
	new.data[,cols] <- apply(X=new.data[,cols], MARGIN=2, FUN= function(x) x - median(x))
	
	return(new.data)
}

#Does a log base 2 transformation, and then normalizes by subtracting the mean
#    value of a set of control genes
#Identities of control genes are supplied in the file pointed to by control.file
control_normalize <- function(dataset, cols=3:18, control.file) {
	new.data <- dataset
	
	#First, do a log-base-2 transformation
	new.data[,cols] <- log2(new.data[,cols])
	
	#Extract the control gene data
	control.genes <- read.csv(file=control.file, colClasses=rep("character", 6))
	control.data <- new.data[as.character(new.data$SYMBOL) %in% control.genes$Symbol,]
	
	#Normalize by subtracting the mean of the control genes from each value
	#Sure, using an apply-like function might be a bit faster but we are not working with many
	#  columns here, and it would be complicated to implement because the control genes
	#  are in a separate data frame
	for (i in cols) {
		new.data[,i] <- new.data[,i] - mean(control.data[,i])
	}
	
	return(new.data)
}

#Takes an Illumina microarray data frame in the wide format and converts
#    it to long format; also generates columns describing the background, genotype,
#    etc. of each array in the dataset
convert.to.long <- function(dataset, label.file) {
	
	translation.table <- read.csv(label.file, colClasses=rep("character", 5))
	#Convert to long format
	#ill.data <- reshape(dataset, idvar="ProbeID", direction="long", varying=3:18,
	#    times=translation.table$Array, v.names="log2expr", timevar="Array")
	ill.data <- reshape(dataset, idvar="ProbeID", direction="long", varying=3:18,
	    times=names(dataset)[3:18], v.names="log2expr", timevar="Array")


	#Add the new columns and fill them in
	ill.data$Background <- rep("", nrow(ill.data))
	ill.data$Genotype <- rep("", nrow(ill.data))
	ill.data$biol_rep <- rep("", nrow(ill.data))
	ill.data$label_rep <- rep("", nrow(ill.data))
	for (i in 1:nrow(translation.table)) {
		which.rows <- (as.character(ill.data$Array) == as.character(translation.table$Array[i]))
		ill.data[which.rows, c("Background", "Genotype", "biol_rep", "label_rep")] <-
		    translation.table[i,c("Background", "Genotype", "biol_rep", "label_rep")]
	}
	
    return(ill.data)
}

#####################################################################################
#####################################################################################

#Read and normalize the data in wide format
ill.data.wide <- read.csv("../microarray_data/sd_array_illumina_R.csv")
ill.data.c <- control_normalize(ill.data.wide, cols=3:18, control.file="../microarray_data/sd_array_illumina_control_list_R.csv")
ill.data.m <- median_normalize(ill.data.wide, cols=3:18)

#Convert the data to long format
ill.data.c <- convert.to.long(ill.data.c, "../microarray_data/sd_array_illumina_labels_R.csv")
ill.data.m <- convert.to.long(ill.data.m, "../microarray_data/sd_array_illumina_labels_R.csv")

#Use the table of array definitions to generate Background, Genotype, biol_rep, and label_rep columns:

#Run gene-specific models:
#Can run probe-specific models instead by replacing INDICES=SYMBOL with INDICES=ProbeID
gene.model <- function(data) {
	return(lm(log2expr ~ Background + Genotype + Background:Genotype, data=data))
}
#The following two lines do a per-gene analysis, but some of the probes suck
#model.results.c <- by(ill.data.c, INDICES=ill.data.c$SYMBOL, FUN=gene.model)
#model.results.m <- by(ill.data.m, INDICES=ill.data.m$SYMBOL, FUN=gene.model)
#So, let's do a probe-level analysis instead
model.results.c <- by(ill.data.c, INDICES=ill.data.c$ProbeID, FUN=gene.model)
model.results.m <- by(ill.data.m, INDICES=ill.data.m$ProbeID, FUN=gene.model)


#Generate a table of q-values from p-values:
library(qvalue)
#Extract the p-values into a matrix
pvals.m <- sapply(model.results.m, FUN= function(x) summary(x)$coefficients[,4], simplify=TRUE)
pvals.m <- t(pvals.m)
pvals.c <- sapply(model.results.c, FUN= function(x) summary(x)$coefficients[,4], simplify=TRUE)
pvals.c <- t(pvals.c)
#Use the qvalue package to obtain q-values from p-values
qvals.m <- apply(X=pvals.m, MARGIN=2, FUN= function(x) qvalue(x)$qvalues)
qvals.c <- apply(X=pvals.c, MARGIN=2, FUN= function(x) qvalue(x)$qvalues)



#Generate a summary table with all the output data
#We want the gene name, FBGN, genotype/background/interaction effect sizes/t-values/p-values
model.results <- model.results.m #Replace with model.results.c if you want to use the control-normalized dataset
extract.data <- function(x) {
	model.summary <- summary(x)
	
	intercept <- model.summary$coefficients[1,1]
	
	genotype.wt.effect <- model.summary$coefficients[3,1]
	genotype.wt.t <- model.summary$coefficients[3, 3]
	genotype.wt.p <- model.summary$coefficients[3, 4]
	
	background.sam.effect <- model.summary$coefficients[2,1]
	background.sam.t <- model.summary$coefficients[2, 3]
	background.sam.p <- model.summary$coefficients[2, 4]
	
	interaction.wt.sam.effect <- model.summary$coefficients[4,1]
	interaction.wt.sam.t <- model.summary$coefficients[4, 3]
	interaction.wt.sam.p <- model.summary$coefficients[4, 4]
	
	result <- cbind(	intercept=intercept,
						genotype.wt.effect=genotype.wt.effect,
						genotype.wt.t=genotype.wt.t,
						genotype.wt.p=genotype.wt.p,
						background.sam.effect=background.sam.effect,
						background.sam.t=background.sam.t,
						background.sam.p=background.sam.p,
						interaction.wt.sam.effect=interaction.wt.sam.effect,
						interaction.wt.sam.t=interaction.wt.sam.t,
						interaction.wt.sam.p=interaction.wt.sam.p
			       )
	
	return(result)
}
summary.table <- sapply(model.results, FUN= extract.data, simplify=TRUE)
summary.table <- data.frame(t(summary.table))
colnames(summary.table) <- c("intercept", "genotype.wt.effect", "genotype.wt.t", "genotype.wt.p", "background.sam.effect", "background.sam.t", "background.sam.p", "interaction.wt.sam.effect", "interaction.wt.sam.t", "interaction.wt.sam.p")
#Generate a "dictionary" that allows us to look up the FBgn's that correspond to gene symbols
fbgn.key <- read.table("../microarray_data/illumina_fbgn.txt", header=TRUE)
fbgn.key[,1] <- as.character(fbgn.key[,1])
fbgn.key[,2] <- as.character(fbgn.key[,2])
#Generate a "dictionary" that allows us to look up the gene symbols that correspond to Probe ID's
probe.gene.key <- ill.data.m[,c("SYMBOL", "ProbeID")]
probe.gene.key <- probe.gene.key[!duplicated(probe.gene.key$ProbeID),]
get.fbgn <- function(x) {
  match <- which(probe.gene.key$ProbeID == x)
  cur.gene <- probe.gene.key$SYMBOL[match]
  match <- which(fbgn.key$symbol == cur.gene)
  return(fbgn.key$fbgn[match])
}
fbgn <- sapply(X=rownames(summary.table), FUN=get.fbgn)
summary.table$fbgn <- fbgn

#Add q-values
summary.table$genotype.wt.q <- rep(NA, nrow(summary.table))
full.rows <- which(!is.na(summary.table$genotype.wt.p))
summary.table$genotype.wt.q[full.rows] <- qvalue(na.exclude(summary.table$genotype.wt.p))$qvalues

summary.table$background.sam.q <- rep(NA, nrow(summary.table))
full.rows <- which(!is.na(summary.table$background.sam.p))
summary.table$background.sam.q[full.rows] <- qvalue(na.exclude(summary.table$background.sam.p))$qvalues

summary.table$interaction.wt.sam.q <- rep(NA, nrow(summary.table))
full.rows <- which(!is.na(summary.table$interaction.wt.sam.p))
summary.table$interaction.wt.sam.q[full.rows] <- qvalue(na.exclude(summary.table$interaction.wt.sam.p))$qvalues

write.csv(summary.table, file="../microarray_data/illumina_summary_results.csv")



#####################################################################################
#####################################################################################
# Some other outputs that might be interesting but which we will not save to include
# in the output file

#Which genes show background dependence in expression?
background.m <- names(which(qvals.m[,2] < 0.1))
background.c <- names(which(qvals.c[,2] < 0.1))

#Which genes are affected by sd[E3]?
geno.m <- names(which(qvals.m[,3] < 0.05))
geno.c <- names(which(qvals.c[,3] < 0.05))

#Which genes show a genotype-by-background interaction?
interaction.m <- names(which(qvals.m[,4] < 0.05))
interaction.c <- names(which(qvals.c[,4] < 0.05))

#For genotype: which ones are common to both normalization methods? Which ones are unique to just one
#  normalization method?
print("Common to both normalization methods:")
geno.m[geno.m %in% geno.c]
print("Unique to median normalization")
geno.m[!(geno.m %in% geno.c)]
print("Unique to control normalization")
geno.c[!(geno.c %in% geno.m)]
