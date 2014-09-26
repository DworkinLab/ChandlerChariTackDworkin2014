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
## Read in the datasets
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

##All dge results -- to test for genotype effects
#dge.results.all <- read.csv('../../expression_scripts/dge_r_analysis/dge_genotype_results.csv', header=TRUE)
#dge.results.all$fbgn <- as.character(dge.results.all$fbgn)
#dge.results.all$gene <- as.character(dge.results.all$gene)

#Allele-specific DGE results
dge.results.as <- read.csv('../../expression_scripts/dge_r_analysis/dge_allele_specific_results.csv', header=TRUE)
dge.results.as$fbgn <- as.character(dge.results.as$fbgn)
dge.results.as$gene <- as.character(dge.results.as$gene)

##Binding site results
#binding.results <- read.csv('../../sequencing_scripts/binding_sites/SAM_ORE_binding_results_bygene_double.csv', header=TRUE)
#colnames(binding.results)[1] <- "fbgn"
#binding.results$fbgn <- as.character(binding.results$fbgn)


##############################################################################
##############################################################################
# Now do some comparisons
##############################################################################
##############################################################################


source('../functions/overlap_test_functions.R')

#Evidence of allelic imbalance overall -- compared to evidence of background differences from expression arrays
#                                      -- compared to genes showing evidence of 



#Compare genes showing evidence of differences between WT and mutant in allelic imbalance, to genes showing
# evidence of g-by-b interaction in arrays

#overlap.test.all(dge.results.as, illumina.results, p.col.1="p.2", p.col.2="interaction.wt.sam.p", thresh.1=0.01, thresh.2=0.05, n.reps=1000)

overlap.test.trimmed(dge.results.as, illumina.results, p.col.1="p.2", p.col.2="interaction.wt.sam.p", thresh.1=0.01, thresh.2=0.05, n.reps=10000)

#overlap.test.all(dge.results.as, dgrc.results, p.col.1="p.2", p.col.2="interaction.wt.sam.p", thresh.1=0.01, thresh.2=0.001, n.reps=1000)

overlap.test.trimmed(dge.results.as, dgrc.results, p.col.1="p.2", p.col.2="interaction.wt.sam.p", thresh.1=0.01, thresh.2=0.001, n.reps=10000)