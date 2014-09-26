col.types <- c("character", "character", "character", "character", "character", "character", "character", "character")
all.data <- read.csv(file="all_semiqt.csv", colClasses = col.types)
valid.numbers <- c("1", "2", "3", "4", "5", "6", "7", "8", "9")
all.data$Numeric_Phenotype <- ifelse(all.data$Numeric_Phenotype %in% valid.numbers, as.numeric(all.data$Numeric_Phenotype), NA)

#Convert all the characters to lowercase so that case differences don't screw up our attempts to match things
text.cols <- c(2, 5, 8)
for (cur.col in text.cols) {
	all.data[,cur.col] <- tolower(all.data[,cur.col])
}

#Numeric_Phenotype ~ Deletion + Background + Deletion:Background

#For each stock:
#	Pull out the data and the appropriate control
#	Run the model

#Pull out all the control data
control.rows <- which(all.data$Deletion=="control")
control.data <- all.data[control.rows,]

#Pull out all the experimental data
exp.rows <- which(all.data$Deletion!="control")
exp.data <- all.data[exp.rows,]


stock.names <- unique(exp.data$Stock_Number)
omitted.count <- 0
#Set up the empty data frame
p.deletion <- rep(NA, length(stock.names))
p.background <- rep(NA, length(stock.names))
p.interaction <- rep(NA, length(stock.names))
deletion.effect <- rep(NA, length(stock.names))
background.effect <- rep(NA, length(stock.names))
interaction.effect <- rep(NA, length(stock.names))
results.table <- data.frame(stock=stock.names, p.deletion=p.deletion, p.background=p.background, p.interaction=p.interaction)
omitted.stocks <- c()
completed.stocks <- c()
for (cur.stock in stock.names) {

	completed.stocks <- c(cur.stock, completed.stocks)

	#Pull out the current experimental data
	cur.exp.data <- exp.data[exp.data$Stock_Number==cur.stock,]
	
	#Pull out the relevant control data
	#We want it to be from the same set (i.e., BSC, DrosDel, or Exelixis) and have the same blocks as the exp data
	wanted.set <- unique(cur.exp.data$Set)
	wanted.blocks <- unique(cur.exp.data$Block)
	wanted.rows <- (control.data$Set == wanted.set) & (control.data$Block %in% wanted.blocks)
	cur.control.data <- control.data[wanted.rows,]
	
	cur.data <- rbind(cur.exp.data, cur.control.data)
	
	#Check if the current stock was measured in both backgrounds...
	measured.in.both <- (length(unique(cur.exp.data$Background)) == 2) 

	
	#Check if we have control data
	has.controls <- (nrow(cur.control.data) > 0) & (length(unique(cur.control.data$Background))==2)
	if (measured.in.both & !has.controls) {
		#If we can't find the proper block-matched controls, let's try to use some control data from other blocks,
		# matching the sample sizes of both SAM and ORE in it
		
		sam.control.rows <- which(control.data$Set==wanted.set & control.data$Background=="SAM")
		num.sam <- sum(cur.exp.data$Background=="SAM")
		sam.control.rows <- sample(sam.control.rows, size=num.sam, replace=FALSE)
		
		ore.control.rows <- which(control.data$Set==wanted.set & control.data$Background=="ORE")
		num.ore <- sum(cur.exp.data$Background=="ORE")
		ore.control.rows <- sample(ore.control.rows, size=num.ore, replace=FALSE)
		
		cur.control.data <- control.data[c(sam.control.rows, ore.control.rows),]
		cur.data <- rbind(cur.exp.data, cur.control.data)
		has.controls <- TRUE
	}
	
	#Make sure ALL of our conditions are met
	if (all(c(measured.in.both, has.controls))) {
		print("**********************************************")
		print(cur.stock)
		cur.model <- lm(Numeric_Phenotype ~ Deletion + Background + Deletion:Background, data=cur.data)
		#print(summary(cur.model))
		
		cur.p.deletion <- summary(cur.model)$coefficients[2,4]
		cur.p.background <- summary(cur.model)$coefficients[3,4]
		cur.p.interaction <- summary(cur.model)$coefficients[4,4]
		cur.deletion.effect <- summary(cur.model)$coefficients[2,1]
		cur.background.effect <- summary(cur.model)$coefficients[3,1]
		cur.interaction.effect <- summary(cur.model)$coefficients[4,1]		
		
		cur.row <- which(results.table$stock==cur.stock)
		results.table$p.deletion[cur.row] <- cur.p.deletion
		results.table$p.background[cur.row] <- cur.p.background
		results.table$p.interaction[cur.row] <- cur.p.interaction
		results.table$deletion.effect[cur.row] <- cur.deletion.effect
		results.table$background.effect[cur.row] <- cur.background.effect
		results.table$interaction.effect[cur.row] <- cur.interaction.effect
	} else {
		omitted.count <- omitted.count + 1
		omitted.stocks <- c(omitted.stocks, cur.stock)
		print(paste("Omitted ", cur.stock))
	}
	
}

print(omitted.count)
print(omitted.stocks)

results.table <- na.omit(results.table)
results.table$stock <- as.character(results.table$stock)

#Now use the dictionary to fill in the start and stop coordinates of each deletion
deletion.dict <- read.delim(file="deletion_stock_dictionary_short_deletion_names.txt", header=FALSE, colClasses=c("character", "character", "character", "numeric", "numeric"))
colnames(deletion.dict) <- c("stock", "df", "chrom", "start", "stop")
deletion.dict <- na.omit(deletion.dict)

#Trim each down to the entries that are common to both, then sort them so they are in the same order
deletion.dict <- deletion.dict[deletion.dict$stock %in% results.table$stock,]
results.table <- results.table[results.table$stock %in% deletion.dict$stock,]
deletion.dict <- deletion.dict[order(deletion.dict$stock),]
results.table <- results.table[order(results.table$stock),]
results.table$chrom <- deletion.dict$chrom
results.table$start <- deletion.dict$start
results.table$stop <- deletion.dict$stop

#Add q-values
library(qvalue)

results.table$q.deletion <- rep(NA, nrow(results.table))
full.rows <- which((!is.na(results.table$p.deletion)) & (is.finite(results.table$p.deletion)))
results.table$q.deletion[full.rows] <- qvalue(results.table$p.deletion[full.rows])$qvalues

#results.table$q.background <- rep(NA, nrow(results.table))
#full.rows <- which((!is.na(results.table$p.background)) & (is.finite(results.table$p.background)))
#results.table$q.background[full.rows] <- qvalue(results.table$p.background[full.rows])$qvalues

results.table$q.interaction <- rep(NA, nrow(results.table))
full.rows <- which((!is.na(results.table$p.interaction)) & (is.finite(results.table$p.interaction)))
results.table$q.interaction[full.rows] <- qvalue(results.table$p.interaction[full.rows])$qvalues

results.table <- results.table[,c("stock", "p.deletion", "p.background", "p.interaction", "chrom", "start", "stop", "q.deletion", "q.interaction", "deletion.effect", "background.effect", "interaction.effect")]

write.csv(results.table, file="raw_deletion_results.csv", row.names=FALSE)
