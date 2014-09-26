#This script reads in the results from all the dge processing scripts in the ../dge directory
# (these results are in the form: 
#   gene symbol, fbgn, gene_transcript_startpos, fc3count, fc4count, fc5count, fc6count
#  just for reference, fc3 and fc4 are the mutant samples, and fc5 and fc6 are the wild-type samples)
#Then analyzes them to look for genes that show evidence of differences between wild-type and mutant genotypes

column.types <- c("character", "character", "character", "numeric", "numeric", "numeric", "numeric")
column.names <- c("gene", "fbgn", "gene_transcript_start", "fc3.count", "fc4.count", "fc5.count", "fc6.count")
all.data <- read.csv("../dge_data/dge_data_all.csv", header=FALSE, col.names=column.names, colClasses=column.types)



#####################################################
#####################################################

# First, analyze the data using edgeR

library(edgeR)
dge.counts <- all.data[,c("fc3.count", "fc4.count", "fc5.count", "fc6.count")]
library.sizes <- apply(all.data[,c("fc3.count", "fc4.count", "fc5.count", "fc6.count")], MARGIN=2, FUN=sum)
dge.data <- DGEList(counts=dge.counts, group=c(1,1,2,2), lib.size=library.sizes)

# Here, allow overdispersion for each tag
prior.n <- 10 #Use a fixed smoothing parameter as recommended by the authors
d <- estimateCommonDisp(dge.data)
d2 <- estimateTagwiseDisp(d, prior.n)

de.tagwise <- exactTest(d2, common.disp=FALSE)

library(qvalue)
p.vals <- de.tagwise$table$p.value
p.vals <- ifelse(p.vals >= 1, 1, p.vals)
q.vals <- qvalue(p.vals)$qvalues

all.data$p <- p.vals
all.data$q <- q.vals


# Now with DEGseq
library(DEGseq)
DEGseq.results <- DEGexp(geneExpMatrix1 = all.data, geneCol1 = 3, expCol1 = c(6,7), groupLabel1="nonmutant",
                         geneExpMatrix2 = all.data, geneCol2 = 3, expCol2 = c(4,5), groupLabel2="mutant",
                         outputDir="DEGseq_output")


# Now compare edgeR to DEGseq
degseq.results <- read.table(file="DEGseq_output/output_score.txt", header=TRUE)
#Need to sort them both first according to transcript
all.data$gene_transcript_start <- as.character(all.data$gene_transcript_start)
all.data <- all.data[order(all.data$gene_transcript_start),]
degseq.results$GeneNames <- as.character(degseq.results$GeneNames)
degseq.results <- degseq.results[order(degseq.results$GeneNames),]
all(as.character(degseq.results$GeneNames) == all.data$gene_transcript_start) #This should be TRUE


axis.lim <- c(0, 1)
png.width <- 480*1.75
png.height <- 480*1.75
png(file="edgeR_DEGseq_comparison.png", width=png.width, height=png.height)
par(cex.axis=2, cex.lab=2, mar=c(6,6,1,1), mgp=c(4,1,0))
plot(all.data$p, degseq.results$p.value, xlim=c(0,1), ylim=axis.lim, xlab="edgeR p-value", ylab="DEGseq p-value", cex=1, pch=16, col="#00000033") 
dev.off()
png(file="edgeR_DEGseq_comparison_log.png", width=png.width, height=png.height)
plot(log(all.data$p), log(degseq.results$p.value), xlab="edgeR log-p", ylab="DEGseq log-p") 
dev.off()

#Most of the points are on or below the diagonal, suggesting that edgeR is more conservative
#plot(all.data$q, degseq.results$q.value.Storey.et.al..2003) #Same
#plot(all.data$q, degseq.results$q.value.Benjamini.et.al..1995) #Same

#lines(x=c(0,1), y=c(0,1), lty=3)
#p.thresh <- 0.001
#lines(x=c(0, 1), y=c(p.thresh, p.thresh), lty=2)
#lines(x=c(p.thresh, p.thresh), y=c(0, 1), lty=2)

thresh <- 0.001
# number significant by edgeR and degseq
length(which(all.data$p < thresh))
length(which(degseq.results$p.value < thresh))
# by q-value
length(which(all.data$q < thresh))
length(which(degseq.results$q.value.Benjamini.et.al..1995 < thresh))


##edgeR is CLEARLY a more conservative approach -- let's go with that one...



#Sort by increasing q-value
all.data <- all.data[order(all.data$q),]

save.image(file="mutant_nonmutant_image.RData")
load(file="mutant_nonmutant_image.RData")

#Write the output results
write.csv(all.data, "dge_genotype_results.csv", row.names=FALSE, quote=FALSE)