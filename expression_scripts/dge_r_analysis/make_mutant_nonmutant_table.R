#This script reads in the results from the edgeR analysis of mutant vs. non-mutant expression profiles
#Then outputs a CSV file containing information on the genes that are significant at some threshold,
# including log fold-change in expression levels as estimated by edgeR

thresh <- 0.001

column.types <- c("character", "character", "character", "numeric", "numeric", "numeric", "numeric")
column.names <- c("gene", "fbgn", "gene_transcript_start", "fc3.count", "fc4.count", "fc5.count", "fc6.count")


load(file="mutant_nonmutant_image.RData")


#Publication stuff
pub.data <- read.csv("../dge_data/dge_data_all.csv", header=FALSE, col.names=column.names, colClasses=column.types)
pub.data$p <- p.vals
pub.data$q <- q.vals
pub.data$logConc <- de.tagwise$table$logConc
pub.data$logFC <- de.tagwise$table$logFC

pub.data <- pub.data[order(pub.data$q),]

pub.data <- pub.data[pub.data$q < thresh,]

#Write the output results
write.csv(pub.data, "publication_genotype_results.csv", row.names=FALSE, quote=FALSE)

