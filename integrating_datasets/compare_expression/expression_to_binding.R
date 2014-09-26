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
# Read in the datasets
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

#Binding site results
#binding.results <- read.csv('../../sequencing_scripts/binding_sites/SAM_ORE_binding_results_bygene_double.csv', header=TRUE)
#binding.results <- read.csv('../../sequencing_scripts/binding_sites_jaspar/SAM_ORE_binding_results_bygene_jaspar.csv', header=TRUE)
binding.results <- read.csv('../../sequencing_scripts/binding_sites_custom_meme/SAM_ORE_binding_results_bygene_custom_meme.csv', header=TRUE)
colnames(binding.results)[1] <- "fbgn"
binding.results$fbgn <- as.character(binding.results$fbgn)


##############################################################################
##############################################################################
# Now do some comparisons
##############################################################################
##############################################################################


source('../functions/overlap_test_functions.R')

#Evidence of binding = high LLR score in binding site prediction
#                    = regulation by sd evidence from any of the expression arrays (wt-mutant comparisons)
#Do this comparison elsewhere! (in comparing expression folder)

#First, to get overall predicted binding targets, we need to find genes that are predicted to bind scalloped in either SAM or ORE
binding.results$min.empirical.p <- ifelse(binding.results$empirical.p.ore < binding.results$empirical.p.sam, binding.results$empirical.p.ore, binding.results$empirical.p.sam)

#Compare genes predicted to bind scalloped, to genes showing evidence of being mis-regulated by scalloped in the presence of sd[E3]
#Illumina
print("Genotype effects from Illumina, compared to predicted scalloped binding targets")
overlap.test.trimmed(illumina.results, binding.results, p.col.1="genotype.wt.p", p.col.2="min.empirical.p", thresh.1=0.05, thresh.2=0.05, n.reps=10000)

#DGRC
print("Genotype effects from DGRC, compared to predicted scalloped binding targets")
overlap.test.trimmed(dgrc.results, binding.results, p.col.1="genotype.wt.q", p.col.2="min.empirical.p", thresh.1=0.05, thresh.2=0.05, n.reps=10000)

#DGE
print("Genotype effects from DGE, compared to predicted scalloped binding targets")
overlap.test.trimmed(dge.results.all, binding.results, p.col.1="q", p.col.2="min.empirical.p", thresh.1=0.05, thresh.2=0.05, n.reps=10000)

#Now compare genotype-by-background effects, to genes showing evidence of sd-binding differences between SAM and ORE
#Illumina
print("Interaction effects from Illumina, compared to sd-binding differences between SAM & ORE")
overlap.test.trimmed(illumina.results, binding.results, p.col.1="interaction.wt.sam.p", p.col.2="empirical.p", thresh.1=0.05, thresh.2=0.05, n.reps=10000)

#DGRC
print("Interaction effects from DGRC, compared to sd-binding differences between SAM & ORE")
overlap.test.trimmed(dgrc.results, binding.results, p.col.1="interaction.wt.sam.q", p.col.2="empirical.p", thresh.1=0.05, thresh.2=0.05, n.reps=10000)


#Lastly, compare genes showing evidence of differences between WT and mutant in allelic imbalance, to genes showing
# evidence of sd-binding differences between SAM and ORE
print("Genotype-dependent allelic imbalance compared to sd-binding differences between SAM & ORE")
#overlap.test.all(dge.results.as, binding.results, p.col.1="p.2", p.col.2="empirical.p", thresh.1=0.01, thresh.2=0.01, n.reps=1000)
overlap.test.trimmed(dge.results.as, binding.results, p.col.1="q.2", p.col.2="empirical.p", thresh.1=0.05, thresh.2=0.05, n.reps=10000)


#Compare genes showing evidence of allelic imbalance, to genes showing
# evidence of sd-binding differences between SAM and ORE
#Note: this comparison is probably not really useful; we don't necessarily *expect* any significant overlap here
#      because lots of genes might show allelic imbalance between SAM and ORE having nothing to do with scalloped at all
# print("Allelic imbalance compared to sd-binding differences between SAM & ORE")
# overlap.test.all(dge.results.as, binding.results, p.col.1="p.1", p.col.2="empirical.p", thresh.1=0.01, thresh.2=0.05, n.reps=1000)
# overlap.test.trimmed(dge.results.as, binding.results, p.col.1="p.1", p.col.2="empirical.p", thresh.1=0.01, thresh.2=0.05, n.reps=1000)

