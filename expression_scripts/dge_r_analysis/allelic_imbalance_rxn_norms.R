ai.data <- read.csv("dge_allele_specific_results.csv")

#Trim down to those tags having at least X reads of each
min.count <- 8

keepers <- (ai.data$fc3.count.ore >= min.count)
keepers <- keepers & (ai.data$fc3.count.sam >= min.count)
keepers <- keepers & (ai.data$fc4.count.ore >= min.count)
keepers <- keepers & (ai.data$fc4.count.sam >= min.count)
keepers <- keepers & (ai.data$fc5.count.ore >= min.count)
keepers <- keepers & (ai.data$fc5.count.sam >= min.count)
keepers <- keepers & (ai.data$fc6.count.ore >= min.count)
keepers <- keepers & (ai.data$fc6.count.sam >= min.count)

ai.data <- ai.data[keepers,]

ai.data$ratio.mut <- ((ai.data$fc3.count.sam/ai.data$fc3.count.ore) + (ai.data$fc4.count.sam/ai.data$fc4.count.ore)) / 2
ai.data$ratio.mut <- log(ai.data$ratio.mut)
#ai.data$ratio.mut <- ifelse(ai.data$ratio.mut < 1, -1/ai.data$ratio.mut, ai.data$ratio.mut)
ai.data$ratio.wt <- ((ai.data$fc5.count.sam/ai.data$fc5.count.ore) + (ai.data$fc6.count.sam/ai.data$fc6.count.ore)) / 2
ai.data$ratio.wt <- log(ai.data$ratio.wt)
#ai.data$ratio.wt <- ifelse(ai.data$ratio.wt < 1, -1/ai.data$ratio.wt, ai.data$ratio.wt)

pdf(file="PublicationAIReactionNorms.pdf", width=4, height=4)

par(mar=c(3,4,1,1))
x.locations <- c(0.05, 0.95)
plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(min(c(ai.data$ratio.mut, ai.data$ratio.wt)) , max(c(ai.data$ratio.mut, ai.data$ratio.wt))), xlab="", ylab="log(SAM:ORE ratio)", xaxt="n")
axis(side=1, at=x.locations, labels=c("WT", "Mutant"))

for (i in 1:nrow(ai.data)) {
	lines(x=x.locations, y=c(ai.data$ratio.wt[i], ai.data$ratio.mut[i]), col="#00000044", lwd=2)
}

dev.off()