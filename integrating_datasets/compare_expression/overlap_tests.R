##############################################################################
# Read in flybase (list of genes with fbgn/names/symbols, and chromosomal locations)##############################################################################

flybase.cols <- c("fbgn", "symbol", "chrom", "start", "stop")
flybase.colClasses <- c("character", "character", "character", "numeric", "numeric")
flybase.data <- read.csv("../../expression_scripts/flybase/flybase_map_table.csv", header=FALSE, colClasses=flybase.colClasses)
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
##############################################################################

#DGRC
dgrc.results <- read.csv('../../expression_scripts/microarray_data/DGRC_summary_results.csv', header=TRUE)
dgrc.results <- dgrc.results[,2:ncol(dgrc.results)]
dgrc.results <- dgrc.results[grep("FBgn", dgrc.results$fbgn),]
dgrc.results$fbgn <- as.character(dgrc.results$fbgn)

#Illumina
illumina.results <- read.csv('../../expression_scripts/microarray_data/illumina_summary_results.csv', header=TRUE)
colnames(illumina.results)[1] <- "gene"
illumina.results <- illumina.results[grep("FBgn", illumina.results$fbgn),]
illumina.results$fbgn <- as.character(illumina.results$fbgn)

#All dge results -- to test for genotype effects
dge.results.all <- read.csv('../../expression_scripts/dge_r_analysis/dge_genotype_results.csv', header=TRUE)
dge.results.all$fbgn <- as.character(dge.results.all$fbgn)
dge.results.all$gene <- as.character(dge.results.all$gene)

#Allele-specific DGE results
dge.results.as <- read.csv('../../expression_scripts/dge_r_analysis/dge_allele_specific_results.csv', header=TRUE)
dge.results.as$fbgn <- as.character(dge.results.as$fbgn)
dge.results.as$gene <- as.character(dge.results.as$gene)


##############################################################################
##############################################################################
# Now do some comparisons
##############################################################################
##############################################################################

# IV. Compare overlap between datasets in the "significant" genes

source('../functions/overlap_test_functions.R')

#print("************DGRC and Illumina -- Genotype effect overlap -- all genes************")
#overlap.test.all(dgrc.results, illumina.results, p.col.1="genotype.wt.p", p.col.2="genotype.wt.p", thresh.1=0.001, thresh.2=0.01, n.reps=10000)

#print("************DGE and Illumina -- Genotype effect overlap -- all genes************")
#overlap.test.all(dge.results.all, illumina.results, p.col.1="q", p.col.2="genotype.wt.p", thresh.1=0.01, thresh.2=0.01, n.reps=10000)

#print("************DGE and DGRC -- Genotype effect overlap -- all genes************")
#overlap.test.all(dge.results.all, dgrc.results, p.col.1="q", p.col.2="genotype.wt.p", thresh.1=0.01, thresh.2=0.001, n.reps=10000)

#print("************DGRC and Illumina -- Background effect overlap -- all genes************")
#overlap.test.all(dgrc.results, illumina.results, p.col.1="background.sam.p", p.col.2="background.sam.p", thresh.1=0.001, thresh.2=0.01, n.reps=10000)

#print("************DGRC and Illumina -- Interaction effect overlap -- all genes************")
#overlap.test.all(dgrc.results, illumina.results, p.col.1="interaction.wt.sam.p", p.col.2="interaction.wt.sam.p", thresh.1=0.001, thresh.2=0.01, n.reps=10000)




print("************DGRC and Illumina -- Genotype effect overlap -- common genes************")
overlap.test.trimmed(dgrc.results, illumina.results, p.col.1="genotype.wt.q", p.col.2="genotype.wt.p", thresh.1=0.05, thresh.2=0.05, n.reps=10000)

print("************DGE and Illumina -- Genotype effect overlap -- common genes************")
overlap.test.trimmed(dge.results.all, illumina.results, p.col.1="q", p.col.2="genotype.wt.p", thresh.1=0.05, thresh.2=0.05, n.reps=10000)

print("************DGE and DGRC -- Genotype effect overlap -- common genes************")
overlap.test.trimmed(dge.results.all, dgrc.results, p.col.1="q", p.col.2="genotype.wt.q", thresh.1=0.05, thresh.2=0.05, n.reps=10000)

print("************DGE and DGRC -- Interaction effect overlap -- common genes**************")
overlap.test.trimmed(dge.results.as, dgrc.results, p.col.1="q.2", p.col.2="interaction.wt.sam.q", thresh.1=0.05, thresh.2=0.05, n.reps=10000)

print("************DGE and Illumina -- Intearction effect overlap -- common genes**************")
overlap.test.trimmed(dge.results.as, illumina.results, p.col.1="q.2", p.col.2="interaction.wt.sam.p", thresh.1=0.05, thresh.2=0.05, n.reps=10000)


#print("************DGRC and Illumina -- Background effect overlap -- common genes************")
#overlap.test.trimmed(dgrc.results, illumina.results, p.col.1="background.sam.p", p.col.2="background.sam.p", thresh.1=0.001, thresh.2=0.01, n.reps=10000)

#print("************DGRC and Illumina -- Interaction effect overlap -- common genes************")
#overlap.test.trimmed(dgrc.results, illumina.results, p.col.1="interaction.wt.sam.p", p.col.2="interaction.wt.sam.p", thresh.1=0.001, thresh.2=0.01, n.reps=10000)

