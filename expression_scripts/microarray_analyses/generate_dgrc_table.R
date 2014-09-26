##########################################################################
library(lme4)


#Performs a mixed model normalization
mixed.model.normalize <- function(dataset) {
  new.data <- na.omit(dataset)
  mm <- lmer(	log2intensity ~ Dye + (1|Array) + (1|(Array:Dye)),
  				data=new.data)
  new.data$resid <- resid(mm)
  return(new.data)
}


#Open and normalize the dataset
raw.dgrc <- read.csv("../microarray_data/DGRC_data_jan19_2006.csv")
#Log-2 transform the data
raw.dgrc$log2intensity <- log2(raw.dgrc$Intensity)
norm.dgrc <- mixed.model.normalize(raw.dgrc)
#Don't need the old dataset anymore, let's free it from memory
rm(raw.dgrc)
norm.dgrc <- norm.dgrc[,c("Spot", "resid", "genotype", "background", "Dye", "Array")]

#Run spot-specific models
spot.model <- function(data) {
  data$Array <- data$Array[1:nrow(data),drop=TRUE]
  #Some of the spot-specific models fail because missing data/NA cause some of the variables to be confounded with one another
  #So we will just skip over those
  success <- try( result <- lmer(resid ~ genotype + background + genotype:background + Dye + (1|Array),
              data=data), silent=TRUE) 
  if (class(success) != "try-error") {
    return(result)
  } else {
    return(-1)
  }
}
spot.ids <- sort(unique(norm.dgrc$Spot))
model.results <- vector("list", length=length(spot.ids))
for (i in 1:length(spot.ids)) {
	#Print a status update about where we are
	if (i %% 200 == 0) print(i)
	cur.spot <- spot.ids[i]
	cur.data <- norm.dgrc[norm.dgrc$Spot==cur.spot,]
	model.results[[i]] <- spot.model(cur.data)
}


#save.image("DGRC_image_temp.Rdata")
#load("DGRC_image_temp.Rdata")

#Generate a summary table with the results:
#Gene; genotype, background, and interaction effect sizes and t-values, and blank spots for p-values; and FBGN
extract.results <- function(x, cur.spot) {
	result <- cbind(	spot=cur.spot, fbgn="NA", genotype.wt.effect = -1, genotype.wt.t = -1, genotype.wt.p = -1,
						background.sam.effect = -1, background.sam.t = -1, background.sam.p = -1,
						interaction.wt.sam.effect = -1, interaction.wt.sam.t = -1, interaction.wt.sam.p = -1)
	if (class(x) != "mer") return(result)
	summary.results <- summary(x)
	
	genotype.wt.effect <- summary.results@coefs[2,1]
	genotype.wt.t <- summary.results@coefs[2,3]
	
	background.sam.effect <- summary.results@coefs[3,1]
	background.sam.t <- summary.results@coefs[3,3]
	
	interaction.wt.sam.effect <- summary.results@coefs[5,1]
	interaction.wt.sam.t <- summary.results@coefs[5,3]
	
	result <- cbind(	spot=cur.spot, fbgn="NA",
						genotype.wt.effect = genotype.wt.effect, genotype.wt.t = genotype.wt.t, genotype.wt.p = 1,
						background.sam.effect = background.sam.effect, background.sam.t = background.sam.t, background.sam.p = 1,
						interaction.wt.sam.effect = interaction.wt.sam.effect, interaction.wt.sam.t = interaction.wt.sam.t, interaction.wt.sam.p = 1)
	return(result)
}
blank.column <- rep(NA, length(model.results))
dgrc.summary.table <- data.frame(	spot=blank.column,
								fbgn=blank.column,
								genotype.wt.effect=blank.column,
								genotype.wt.t=blank.column,
								genotype.wt.p=blank.column,
								background.sam.effect=blank.column,
								background.sam.t=blank.column,
								background.sam.p=blank.column,
								interaction.wt.sam.effect=blank.column,
								interaction.wt.sam.t=blank.column,
								interaction.wt.sam.p=blank.column
							)
for (i in 1:nrow(dgrc.summary.table)) {
	#Print a status update about where we are
	if (i %% 200 == 0) print(i)
	cur.spot <- spot.ids[i]
	dgrc.summary.table[i,] <- extract.results(model.results[[i]], cur.spot)
}

#save.image("DGRC_image_temp_2.Rdata")
#load("DGRC_image_temp_2.Rdata")

#Estimating p-values:
#lmer refuses to do it (because of statistical issues)
#but we can come up with an estimate, flawed as it may be, from the t-value
# number of observations: nrow(model.results[[x]]@frame)

#df = (number of observations) - (number of parameters estimated)
#	We've estimated: intercept, dye effect, genotype effect, background effect, and a bunch of array effects
#		and # array effects = (# observations / 2) - 1
#       and don't forget we estimated the residual variance

#So degrees of freedom = number of observations - ([# obs/2] - 1) - 5 = (# observations / 2) - 4
#Note, this should be yield conservative estimate of p-values, because we are using a lower bound
# estimate for the degrees of freedom

#Convert all the numeric columns to numeric, otherwise R gives an error message
for (i in 3:11) dgrc.summary.table[,i] <- as.numeric(dgrc.summary.table[,i])

#Now estimate the p-values for ourselves
genotype.wt.p <- rep(1, nrow(dgrc.summary.table))
background.sam.p <- rep(1, nrow(dgrc.summary.table))
interaction.wt.sam.p <- rep(1, nrow(dgrc.summary.table))
counter <- 1
for (n in 1:nrow(dgrc.summary.table)) {
	if (counter %% 200 == 0) {
		print(counter)
	}
	if (class(model.results[[n]]) == "mer") {
		df = (nrow(model.results[[n]]@frame) / 2) - 4
		genotype.wt.p[n] <- 2*pt(q=abs(dgrc.summary.table$genotype.wt.t[n]), df=df, lower.tail=FALSE)
		background.sam.p[n] <- 2*pt(q=abs(dgrc.summary.table$background.sam.t[n]), df=df, lower.tail=FALSE)
		interaction.wt.sam.p[n] <- 2*pt(q=abs(dgrc.summary.table$interaction.wt.sam.t[n]), df=df, lower.tail=FALSE)
	}
	counter = counter+1
}
dgrc.summary.table$genotype.wt.p <- genotype.wt.p
dgrc.summary.table$background.sam.p <- background.sam.p
dgrc.summary.table$interaction.wt.sam.p <- interaction.wt.sam.p

#Add the intercepts to the table
intercept <- rep(NA, nrow(dgrc.summary.table))
for (n in 1:nrow(dgrc.summary.table)) {
	if (class(model.results[[n]]) == "mer") {
		intercept[n] <- summary(model.results[[n]])@coefs[1,1]
	}
	if (n %% 200 == 0) print(n)
}
dgrc.summary.table$intercept <- intercept

#Find the genes that correspond to the spots
spot.data <- read.csv("../microarray_data/DGRC_spot_gene_names.csv")
spot.data <- spot.data[,c("GENE_SYMBOL", "spot")]

#Now fill in the FBgn ID's that correspond to the spots
fbgn.key <- read.csv("../microarray_data/DGRC_spot_gene_names.csv")
fbgn.key[,1] <- as.character(fbgn.key[,1])
fbgn.key[,2] <- as.character(fbgn.key[,2])
fbgn.key[,3] <- as.character(fbgn.key[,3])
fbgn.key[,4] <- as.character(fbgn.key[,4])
fbgn.lookup <- function(x) {
	match <- which(fbgn.key$spot == x)
	return(fbgn.key$FLYBASE_PRIMARY[match])
}
gene.lookup <- function(x) {
	match <- which(fbgn.key$spot == x)
	return(fbgn.key$GENE_SYMBOL[match])
}
for (i in 1:nrow(dgrc.summary.table)) {
	if (i %% 200 == 0) print(i)
	dgrc.summary.table$fbgn[i] <- fbgn.lookup(dgrc.summary.table$spot[i])
	dgrc.summary.table$gene[i] <- gene.lookup(dgrc.summary.table$spot[i])
}


#Add q-values
library(qvalue)
#Make sure to take care of NAs in the results
full.rows <- which(!is.na(dgrc.summary.table$genotype.wt.p))
dgrc.summary.table$genotype.wt.q <- rep(NA, nrow(dgrc.summary.table))
dgrc.summary.table$genotype.wt.q[full.rows] <- qvalue(na.exclude(dgrc.summary.table$genotype.wt.p))$qvalues

full.rows <- which(!is.na(dgrc.summary.table$background.sam.p))
dgrc.summary.table$background.sam.q <- rep(NA, nrow(dgrc.summary.table))
dgrc.summary.table$background.sam.q[full.rows] <- qvalue(na.exclude(dgrc.summary.table$background.sam.p))$qvalues

full.rows <- which(!is.na(dgrc.summary.table$interaction.wt.sam.p))
dgrc.summary.table$interaction.wt.sam.q <- rep(NA, nrow(dgrc.summary.table))
dgrc.summary.table$interaction.wt.sam.q[full.rows] <- qvalue(na.exclude(dgrc.summary.table$interaction.wt.sam.p))$qvalues


#Save the results
#save.image("DGRC_image.Rdata")
#load("DGRC_image.Rdata")
write.csv(dgrc.summary.table, file="../microarray_data/DGRC_summary_results.csv")