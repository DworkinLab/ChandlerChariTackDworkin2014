column.types <- c("character", "character", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric")
column.names <- c("gene", "fbgn", "gene_transcript_start", "fc3.count.ore", "fc3.count.sam", "fc4.count.ore", "fc4.count.sam", "fc5.count.ore", "fc5.count.sam", "fc6.count.ore", "fc6.count.sam")
all.data <- read.csv("../dge_data/dge_data_allele_specific.csv", header=FALSE, col.names=column.names, colClasses=column.types)


#What we need to do is look at replicates individually and compare models
#	Null: 0.5 prob for each allele for both mutant and wild-type (nothing going on)
#	A: Uneven prob for each allele, but same across both mutant and wild-type (allele-specific expression, analogous to a background effect)
#	B: Uneven prob for each allele, and unequal between mutant and wild-type (allele-by-genotype interaction)

#Get the negative log likelihood for the null model for each gene
nll.null <- -dbinom(all.data$fc3.count.ore, all.data$fc3.count.ore+all.data$fc3.count.sam, prob=0.5, log=TRUE)
nll.null <- nll.null - dbinom(all.data$fc4.count.ore, all.data$fc4.count.ore+all.data$fc4.count.sam, prob=0.5, log=TRUE)
nll.null <- nll.null - dbinom(all.data$fc5.count.ore, all.data$fc5.count.ore+all.data$fc5.count.sam, prob=0.5, log=TRUE)
nll.null <- nll.null - dbinom(all.data$fc6.count.ore, all.data$fc6.count.ore+all.data$fc6.count.sam, prob=0.5, log=TRUE)
all.data$nll.null <- nll.null

#Get the negative log likelihood for model A for each gene
#p.A = maximum likelihood estimate for the allelic imbalance
p.A <- (all.data$fc3.count.ore + all.data$fc4.count.ore + all.data$fc5.count.ore + all.data$fc6.count.ore) / 
       (all.data$fc3.count.ore + all.data$fc4.count.ore + all.data$fc5.count.ore + all.data$fc6.count.ore + all.data$fc3.count.sam + all.data$fc4.count.sam + all.data$fc5.count.sam + all.data$fc6.count.sam)
nll.A <- -dbinom(all.data$fc3.count.ore, all.data$fc3.count.ore+all.data$fc3.count.sam, prob=p.A, log=TRUE)
nll.A <- nll.A - dbinom(all.data$fc4.count.ore, all.data$fc4.count.ore+all.data$fc4.count.sam, prob=p.A, log=TRUE)
nll.A <- nll.A - dbinom(all.data$fc5.count.ore, all.data$fc5.count.ore+all.data$fc5.count.sam, prob=p.A, log=TRUE)
nll.A <- nll.A - dbinom(all.data$fc6.count.ore, all.data$fc6.count.ore+all.data$fc6.count.sam, prob=p.A, log=TRUE)
all.data$nll.A <- nll.A

#Get the negative log likelihood for model B for each gene
#p.B.wt and p.B.mut = maximum likelihood estimates for the allelic imbalance for wild-type and mutant genotypes, respectively
p.B.wt <- (all.data$fc5.count.ore + all.data$fc6.count.ore) / (all.data$fc5.count.ore + all.data$fc6.count.ore + all.data$fc5.count.sam + all.data$fc6.count.sam)
p.B.mut <- (all.data$fc3.count.ore + all.data$fc4.count.ore) / (all.data$fc3.count.ore + all.data$fc4.count.ore + all.data$fc3.count.sam + all.data$fc4.count.sam)
nll.B <- -dbinom(all.data$fc3.count.ore, all.data$fc3.count.ore+all.data$fc3.count.sam, prob=p.B.mut, log=TRUE)
nll.B <- nll.B - dbinom(all.data$fc4.count.ore, all.data$fc4.count.ore+all.data$fc4.count.sam, prob=p.B.mut, log=TRUE)
nll.B <- nll.B - dbinom(all.data$fc5.count.ore, all.data$fc5.count.ore+all.data$fc5.count.sam, prob=p.B.wt, log=TRUE)
nll.B <- nll.B - dbinom(all.data$fc6.count.ore, all.data$fc6.count.ore+all.data$fc6.count.sam, prob=p.B.wt, log=TRUE)
all.data$nll.B <- nll.B

p.1 <- pchisq(q=2*(all.data$nll.null - all.data$nll.A), df=1, lower.tail=FALSE)
all.data$p.1 <- p.1
p.2 <- pchisq(q=2*(all.data$nll.A - all.data$nll.B), df=1, lower.tail=FALSE)
all.data$p.2 <- p.2
p.3 <- pchisq(q=2*(all.data$nll.null - all.data$nll.B), df=2, lower.tail=FALSE)
all.data$p.3 <- p.3

interactions <- all.data[(all.data$p.2 < 0.001) & is.finite(all.data$p.2),]
interactions <- interactions[order(interactions$p.2),]
allele.specific <- all.data[(all.data$p.1 < 0.001) & is.finite(all.data$p.1),]
allele.specific <- allele.specific[order(allele.specific$p.1),]

library(qvalue)

#Interesting -- there are a lot more genes that have genotype-specific allelic imbalance than genes that have allelic imbalance on average
#How many genes show evidence of allelic imbalance in each genotype?
#FC5&6 = WT
#FC3&4 = mut
get.allele.specific.rows <- function(s1.ore, s1.sam, s2.ore, s2.sam, thresh=0.001) {
	#Null model
	nll.null <- -dbinom(s1.ore, s1.ore+s1.sam, prob=0.5, log=TRUE)
	nll.null <- nll.null - dbinom(s2.ore, s2.ore+s2.sam, prob=0.5, log=TRUE)
	
	#Model 1: allele-specific
	p.A <- (s1.ore + s2.ore) / (s1.ore + s2.ore + s1.sam + s2.sam)
	nll.A <- -dbinom(s1.ore, s1.ore+s1.sam, prob=p.A, log=TRUE)
	nll.A <- nll.A - dbinom(s2.ore, s2.ore+s2.sam, prob=p.A, log=TRUE)

	p.1 <- pchisq(q=2*(nll.null - nll.A), df=1, lower.tail=FALSE)
	q.1 <- qvalue(na.exclude(p.1))$qvalues
	return(which(q.1 <= thresh))
}

thresh <- 0.001
#Row indices for transcripts that show allelic imbalance, on average
as.avg <- which(all.data$p.1 < thresh)
print(length(as.avg))

#Row indices that show genotype-dependent allelic imbalance (model 2 better than model 1)
as.gd.1 <- which(all.data$p.2 < thresh)
print(length(as.gd.1))

#Row indices that show genotype-dependent allelic imbalance (model 2 better than null model)
as.gd.2 <- which(all.data$p.3 < thresh)
print(length(as.gd.2))

#Wild-types
thresh <- 0.001
as.wt <- get.allele.specific.rows(all.data$fc5.count.ore, all.data$fc5.count.sam, all.data$fc6.count.ore, all.data$fc6.count.sam, thresh=thresh)
print("Number of tags with allelic imbalance in WT:")
print(length(as.wt))


#Mutants
as.mut <- get.allele.specific.rows(all.data$fc3.count.ore, all.data$fc3.count.sam, all.data$fc4.count.ore, all.data$fc4.count.sam, thresh=thresh)
print("Number of tags with allelic imbalance in sdE3 mutants:")
print(length(as.mut))

print("Number they have in common:")
print(length(intersect(as.wt, as.mut)))

#Add q-values

full.rows <- which(!is.na(all.data$p.1))
all.data$q.1 <- rep(NA, nrow(all.data))
all.data$q.1[full.rows] <- qvalue(na.exclude(all.data$p.1))$qvalues

full.rows <- which(!is.na(all.data$p.2))
all.data$q.2 <- rep(NA, nrow(all.data))
all.data$q.2[full.rows] <- qvalue(na.exclude(all.data$p.2))$qvalues

full.rows <- which(!is.na(all.data$p.3))
all.data$q.3 <- rep(NA, nrow(all.data))
all.data$q.3[full.rows] <- qvalue(na.exclude(all.data$p.3))$qvalues


write.csv(all.data, file="dge_allele_specific_results.csv", quote=FALSE, row.names=FALSE)