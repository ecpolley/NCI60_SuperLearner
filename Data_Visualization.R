# Code to run visualizations on approved and investigational agents
# created 2016-03-07
# Eric Polley

## load libraries
library(networkD3)
library(gplots)
library(ggplot2)
library(d3heatmap)

corDist <- function(X) {
	tmp <- (1 - cor(t(X), use = "pairwise.complete.obs"))
	tmp[is.na(tmp)] <- 1
	return(tmp)
}

## running updated variant list with updated GI50 data for the approved and investigational drugs
# ALT: load(file.path("Data", "VariantTableByGene.RData")) # contains VariantTable
VariantTable <- read.csv(file.path("Data", "VariantTableByGene.csv"), row.names = 1)

## Load drug/agent GI50 data
DrugDat <- read.csv(file.path("Data", "AOD_IOA_GI50.csv"))

# find overlap in cell lines
setdiff(rownames(VariantTable), unique(DrugDat$CELL))
setdiff(unique(DrugDat$CELL), rownames(VariantTable))

# fix names in exome data
rownames(VariantTable) <- sub("A549_ATCC", "A549/ATCC", rownames(VariantTable))
rownames(VariantTable) <- sub("COLO-205", "COLO 205", rownames(VariantTable))
rownames(VariantTable) <- sub("DU145", "DU-145", rownames(VariantTable))
rownames(VariantTable) <- sub("HCC2998", "HCC-2998", rownames(VariantTable))
rownames(VariantTable) <- sub("HL-60", "HL-60(TB)", rownames(VariantTable))
rownames(VariantTable) <- sub("Hs_578T", "HS 578T", rownames(VariantTable))
rownames(VariantTable) <- sub("IGR-OV1", "IGROV1", rownames(VariantTable))
rownames(VariantTable) <- sub("K562", "K-562", rownames(VariantTable))
rownames(VariantTable) <- sub("LOX_IMVI", "LOX IMVI", rownames(VariantTable))
rownames(VariantTable) <- sub("MDA-MB-231", "MDA-MB-231/ATCC", rownames(VariantTable))
rownames(VariantTable) <- sub("NCI-ADR-RES", "NCI/ADR-RES", rownames(VariantTable))
rownames(VariantTable) <- sub("RXF-393", "RXF 393", rownames(VariantTable))

# restrict drug data to overlap with exome data
DrugDat <- DrugDat[which(DrugDat$CELL %in% rownames(VariantTable)), ]

GI50Wide <- reshape(DrugDat[, c("NSC", "CELL", "NLOGGI50")], direction = "wide", timevar = "NSC", idvar = "CELL")
colnames(GI50Wide) <- sub("NLOGGI50.", "NSC", colnames(GI50Wide))

# use cell line name as rowname
rownames(GI50Wide) <- GI50Wide[, 1]
# remove cell line name, and line up cell lines with Variant Table
GI50Wide <- GI50Wide[rownames(VariantTable), -1]
all.equal(rownames(VariantTable), rownames(GI50Wide))

## filter VariantTable
CountVar <- colSums(VariantTable)
VariantTable_sub <- as.data.frame(VariantTable[, CountVar > 4])

## clean variable names in variant table
colnames(VariantTable_sub) <- make.names(colnames(VariantTable_sub)) # some genes have "-" in the name, but this creates problems for regression models/formulas, replace with "."

## Table with drug names and NSC codes for reference
DrugNames <- unique(DrugDat[, c("NSC", "DRUGNAME")])
dim(DrugNames) # 181 agents

## filter drug with low variablility?
DrugRange <- apply(GI50Wide, 2, function(x) diff(range(x, na.rm = TRUE)))
# GI50Wide[, which(DrugRange <= 0.5)]

GI50Wide_sub <- GI50Wide[, which(DrugRange > 0.5)]

# How many missing?
NumCellLinesMissing <- apply(GI50Wide_sub, 2, function(x) sum(is.na(x))) 
table(NumCellLinesMissing)  # should we remove some with more than a quarter missing?

## add drug names?
A <- data.frame(NSC = colnames(GI50Wide_sub), stringsAsFactors = FALSE)
B <- data.frame(NSC = paste0("NSC", DrugNames$NSC), DrugName = as.character(DrugNames$DRUGNAME), stringsAsFactors = FALSE)
C <- merge(A, B, sort = FALSE)
all.equal(A$NSC, C$NSC)  # must be true
GI50Wide_sub_names <- GI50Wide_sub
colnames(GI50Wide_sub_names) <- paste0(C$NSC, "_", C$DrugName)

CorDat <- cor(as.matrix(GI50Wide_sub_names), use = "pairwise.complete.obs", method = "pearson")

AllCors <- CorDat[upper.tri(CorDat)]
summary(AllCors)
hist(AllCors)


dim(CorDat)
# the NAs (sd = 0)
data.frame(A = rep(rownames(CorDat), times = nrow(CorDat))[which(is.na(CorDat))], B = rep(rownames(CorDat), each = nrow(CorDat))[which(is.na(CorDat))]) # mostly from missing data and then no variation

# search for cut-off
# test a few values for CorCut
CorCut.seq <- seq(from = 0.60, to = 0.90, by = 0.005)
CorCut.out <- data.frame(CorCut = CorCut.seq, r.squared = NA)

for(jj in seq(nrow(CorCut.out))) {
	tNetworkData <- data.frame(Source = NULL, Target = NULL)
	
	for(ii in seq(nrow(CorDat))) {
		tmpLink <- colnames(CorDat)[which(CorDat[ii, ] > CorCut.seq[jj])]
		tmpNetworkData <- data.frame(Source = rownames(CorDat)[ii], Target = tmpLink)
		tNetworkData <- rbind(tNetworkData, tmpNetworkData)
	}
	dFrame <- as.data.frame(table(table(tNetworkData$Source)))
	dFrame$Var2 <- as.numeric(as.character(dFrame$Var1))
	dFrame$Prob <- dFrame$Freq/sum(dFrame$Freq)
	# plot(y = log10(dFrame$Prob), x = log10(dFrame$Var2))
	fit <- lm(log(Prob)~0+log(Var2), data = dFrame)
	CorCut.out[jj, 2] <- summary(fit)$r.squared
}
plot(x = CorCut.out$CorCut, y = CorCut.out$r.squared, type = "b")


CorCut <- 0.80
NetworkData <- data.frame(Source = NULL, Target = NULL, Value = NULL)
for(ii in seq(nrow(CorDat))) {
	tmpLink <- colnames(CorDat)[which(CorDat[ii, ] > CorCut)]
	tmpCor <- CorDat[ii, which(CorDat[ii, ] > CorCut)]
	tmpNetworkData <- data.frame(Source = rownames(CorDat)[ii], Target = tmpLink, Value = tmpCor)
	NetworkData <- rbind(NetworkData, tmpNetworkData)
}

NodesData <- data.frame(Name = unique(NetworkData$Source), NodeNumber = seq.int(from = 0, to = length(unique(NetworkData$Source)) - 1), Group = 1)

NetworkData$Source <- as.numeric(factor(NetworkData$Source, levels = NodesData$Name, labels = NodesData$NodeNumber)) - 1
NetworkData$Target <- as.numeric(factor(NetworkData$Target, levels = NodesData$Name, labels = NodesData$NodeNumber)) - 1
rownames(NetworkData) <- NULL

forceNetwork(Links = NetworkData, Nodes = NodesData, Source = "Source", Target = "Target", Value = "Value", NodeID = "Name", Group = "Group", height = 1600, width = 1600, zoom = FALSE)


CorCut <- 0.70
NetworkData <- data.frame(Source = NULL, Target = NULL, Value = NULL)
for(ii in seq(nrow(CorDat))) {
	tmpLink <- colnames(CorDat)[which(CorDat[ii, ] > CorCut)]
	tmpCor <- CorDat[ii, which(CorDat[ii, ] > CorCut)]
	tmpNetworkData <- data.frame(Source = rownames(CorDat)[ii], Target = tmpLink, Value = tmpCor)
	NetworkData <- rbind(NetworkData, tmpNetworkData)
}

NodesData <- data.frame(Name = unique(NetworkData$Source), NodeNumber = seq.int(from = 0, to = length(unique(NetworkData$Source)) - 1), Group = 1)

NetworkData$Source <- as.numeric(factor(NetworkData$Source, levels = NodesData$Name, labels = NodesData$NodeNumber)) - 1
NetworkData$Target <- as.numeric(factor(NetworkData$Target, levels = NodesData$Name, labels = NodesData$NodeNumber)) - 1
rownames(NetworkData) <- NULL

forceNetwork(Links = NetworkData, Nodes = NodesData, Source = "Source", Target = "Target", Value = "Value", NodeID = "Name", Group = "Group", height = 1000, width = 1000, zoom = TRUE, fontSize = 12, bounded = TRUE)


## d3heatmap
d3heatmap(t(GI50Wide_sub_names), dendrogram = "both", hclustfun = function(x) hclust(x, method = "average"), distfun = function(x) as.dist(corDist(x)), colors = "Blues", width = 1500, height = 2500, yaxis_font_size = "8px", theme = "dark", show_grid = TRUE)