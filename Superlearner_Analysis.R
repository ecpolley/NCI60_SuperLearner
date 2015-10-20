# Code to run SuperLearner on approved and investigational agents
# created 2012-01-11
# Eric Polley

## load libraries
library(glmnet)
library(rpart)
library(SuperLearner)
library(gbm)
library(gplots)
library(mboost)

## running updated variant list with updated GI50 data for the approved and investigational drugs
# ALT: load(file.path("Data", "VariantTableByGene.RData")) # contains VariantTable
VariantTable <- read.csv(file.path("Data", "VariantTableByGene.csv"))

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

######## Now fit the regression models
# SL.glmnet1se <- function(...) SL.glmnet(..., useMin = FALSE)
SL.rpartPrune2 <- function(..., cp = 0.001, minsplit = 10, maxdepth = 3, xval = 20, minbucket = 4) SL.rpartPrune(..., cp = cp, minsplit = minsplit, xval = xval, maxdepth = maxdepth, minbucket = minbucket)

SL.gbmOOB <- function(Y, X, newX, family, obsWeights, gbm.trees = 1000, interaction.depth = 2, ...) 
{
    require("gbm")
    gbm.model <- as.formula(paste("Y~", paste(colnames(X), collapse = "+")))
    if (family$family == "gaussian") {
        fit.gbm <- gbm(formula = gbm.model, data = X, distribution = "gaussian", 
            n.trees = gbm.trees, interaction.depth = interaction.depth, 
            cv.folds = 0, keep.data = TRUE, weights = obsWeights, 
            verbose = FALSE)
    }
    if (family$family == "binomial") {
        fit.gbm <- gbm(formula = gbm.model, data = X, distribution = "bernoulli", 
            n.trees = gbm.trees, interaction.depth = interaction.depth, 
            cv.folds = 0, keep.data = TRUE, verbose = FALSE, 
            weights = obsWeights)
    }
    best.iter <- gbm.perf(fit.gbm, method = "OOB", plot.it = FALSE)
    pred <- predict(fit.gbm, newdata = newX, best.iter, type = "response")
    fit <- list(object = fit.gbm, n.trees = best.iter)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("SL.gbm")
    return(out)
}
SL.glmboost <- function(Y, X, newX, family, obsWeights, mstop = 1000, centerBoost = FALSE, ...) 
{
    require("mboost")
    if (family$family == "gaussian") {
        fit.gbm <- glmboost(y = Y, x = as.matrix(X), family = GaussReg(), control = boost_control(mstop = mstop, nu = 0.1, risk = "inbag"), center = centerBoost)
    }
    if (family$family == "binomial") {
        stop("not yet implemented")
    }
    pred <- predict(fit.gbm, newdata = as.matrix(newX), type = "response")
    fit <- list(object = fit.gbm)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("SL.glmboost")
    return(out)
}
predict.SL.glmboost <- function(object, newdata, family, ...) {
  require("mboost")
  if (family$family == "binomial") stop("not yet implemented")
  pred <- predict(object$object, newdata = as.matrix(newdata), type = "response")
  return(pred)
}
create.SL.glmnet <- function(alpha = c(0.05, 0.25, 0.50, 0.75)) {
  for(mm in seq(length(alpha))){
    eval(parse(text = paste('SL.glmnet.', alpha[mm], '<- function(..., alpha = ', alpha[mm], ') SL.glmnet(..., alpha = alpha)', sep = '')), envir = .GlobalEnv)
  }
  invisible(TRUE)
}
create.SL.glmnet()
screen.corRank5 <- function(..., rank = 5) screen.corRank(..., rank = rank)
screen.corRank10 <- function(..., rank = 10) screen.corRank(..., rank = rank)
screen.corRank20 <- function(..., rank = 20) screen.corRank(..., rank = rank)


SL.library <- list(c("SL.gbmOOB", "All", "screen.corRank10", "screen.corRank20"), 
  c("SL.glmboost", "All", "screen.corRank10", "screen.corRank20"), 
  c("SL.ipredbagg", "All", "screen.corRank10", "screen.corRank20"), 
  c("SL.mean", "All"), 
  c("SL.glmnet", "All", "screen.corRank20"), 
  c("SL.glmnet.0.75", "All", "screen.corRank20"),
  c("SL.glmnet.0.5", "All", "screen.corRank20"), 
  c("SL.glmnet.0.25", "All", "screen.corRank20"), 
  c("SL.glmnet.0.05", "All", "screen.corRank20"), 
  c("SL.rpartPrune2", "All", "screen.corRank10", "screen.corRank20"), 
  c("SL.randomForest", "All", "screen.corRank10", "screen.corRank20"), 
  c("SL.glm", "screen.corRank5", "screen.corRank10", "screen.corRank20"),
  c("SL.nnet", "screen.corRank5", "screen.corRank10", "screen.corRank20"),
  c("SL.svm", "screen.corRank5", "screen.corRank10", "screen.corRank20"))

dir.create("RegOutput") # folder to save output
for(ii in seq_along(GI50Wide)) {
  Y <- GI50Wide[, ii]
  X <- VariantTable_sub[!is.na(Y), ]
  Y <- Y[!is.na(Y)]
  N <- length(Y) # need to save N so we can do LOO CV
  outSL <- SuperLearner(Y = Y, X = X, newX = VariantTable_sub, SL.library = SL.library, verbose = TRUE, cvControl = list(V = N))
  save(outSL, file = paste("./RegOutput/outSL", ii, ".RData", sep = ""))
	print(Sys.time())
  cat(ii, "\n")
  rm(outSL)
  gc()
}


# ## Files are saved in folder ./RegOutput, read and and put the list together
# coefSL <- matrix(NA, nrow = 35, ncol = ncol(GI50Wide))
# rownames(coefSL) <- c("SL.gbmOOB_All", "SL.gbmOOB_screen.corRank10", "SL.gbmOOB_screen.corRank20", 
# "SL.glmboost_All", "SL.glmboost_screen.corRank10", "SL.glmboost_screen.corRank20", 
# "SL.ipredbagg_All", "SL.ipredbagg_screen.corRank10", "SL.ipredbagg_screen.corRank20", 
# "SL.mean_All", "SL.glmnet_All", "SL.glmnet_screen.corRank20", 
# "SL.glmnet.0.75_All", "SL.glmnet.0.75_screen.corRank20", "SL.glmnet.0.5_All", 
# "SL.glmnet.0.5_screen.corRank20", "SL.glmnet.0.25_All", "SL.glmnet.0.25_screen.corRank20", 
# "SL.glmnet.0.05_All", "SL.glmnet.0.05_screen.corRank20", "SL.rpartPrune2_All", 
# "SL.rpartPrune2_screen.corRank10", "SL.rpartPrune2_screen.corRank20", 
# "SL.randomForest_All", "SL.randomForest_screen.corRank10", "SL.randomForest_screen.corRank20", 
# "SL.glm_screen.corRank5", "SL.glm_screen.corRank10", "SL.glm_screen.corRank20", 
# "SL.nnet_screen.corRank5", "SL.nnet_screen.corRank10", "SL.nnet_screen.corRank20", 
# "SL.svm_screen.corRank5", "SL.svm_screen.corRank10", "SL.svm_screen.corRank20")
# colnames(coefSL) <- names(GI50Wide)
# Zlist <- vector("list", ncol(GI50Wide))
# names(Zlist) <- names(GI50Wide)
# 
# for(nn in seq(ncol(GI50Wide))) {
#   load(paste("./RegOutput/outSL", nn, ".RData", sep = ""))
#   coefSL[, nn] <- coef(outSL)
#   Zlist[[nn]] <- outSL$Z
#   rm(outSL)
#   cat(nn, "\n")
# }
# outSLlist <- lapply(seq(ncol(GI50Wide)), function(n) get(paste("outSL", n, sep = "")))
# names(outSLlist) <- names(GI50Wide)
# save(outSLlist, file = "outSLlist5.RData")
