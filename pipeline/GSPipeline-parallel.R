#!/usr/bin/Rscript
# Pipeline for analysis of GS models

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) 
{
  messages = unlist (list (...))
  cat ("\n>>>>", messages, "\n")
}
#-------------------------------------------------------------


options (width=300)
cmmArgs = commandArgs(trailingOnly = TRUE)
cmmArgs = c ("TRAIN47K-100.csv", "TARGET47K-100.csv", "YieldBLUE-tbl.csv")

USAGE="GSPipeline.R <Training genotype filename> <Target genotype filename> <Phenotype filename>"
if (length (cmmArgs) != 3) {
	message (paste0("\n\n", USAGE, "\n\n"))
	quit ()
}

source ("lglib01.R")

library (rrBLUP)
library (BGLR)
library (brnn)
library (glmnet)
library (e1071) 
library (randomForest)

library (BWGS)
library (parallel)

# Read arguments
genoTrainFile  = cmmArgs [1]
genoTargetFile = cmmArgs [2]
phenoFile      = cmmArgs [3] 

# Load data
genoTrain   = read.table (genoTrainFile, check.names=F)
genoTarget  = read.table (genoTargetFile, check.names=F)
phenoTbl    = read.table (phenoFile, check.names=F, header=T)
phenoVector = phenoTbl [,2]
names (phenoVector) = phenoTbl [,1]
hd (phenoVector)

# Set algorithm parameters
NFOLDS  = 3
NTIMES  = 3
NCORES  = 8
TRAIT   = colnames (phenoTbl)[2]

METHODS = c ("gblup", "EGBLUP", "BA", "BL", "RF", "RKHS")
#METHODS = c ("EN", "RR", "gblup")

#################### Prediction ability comparison of methods ##################
########################## a.Execute prediction methods ######################
msg ("Cross Validation...")

cvFun <- function (method) {
	cvOutput <- bwgs.cv (genoTrain, phenoVector, predict.method= method, 
			  geno.impute.method="mni", nFolds=NFOLDS, nTimes=NTIMES )
	return (cvOutput)
}

msg ("Run cross validation...")
listCVs         = mclapply (METHODS, cvFun, mc.cores=NCORES)
names (listCVs) = METHODS

msg ("Get info from cross validation...")
outFilename             = "out-cross-validation"
predictedAbilitiesList  = mclapply (listCVs, function (cvObject ){return (cvObject$cv)})
predictedTimesList      = mclapply (listCVs, function (cvObject ){return (cvObject$summary[[3]])})
predictedAbilitiesMeans = mclapply (listCVs, function (cvObject ){return (cvObject$summary[[1]])})
bestMethodName          = METHODS [which.max (predictedAbilitiesMeans)]

predictedAbilitiesList
predictedAbilitiesMeans

msg ("Outputs for predicted ability by methods...")
compareM    = as.data.frame (predictedAbilitiesList)
write.csv (compareM, paste0 (outFilename, "-predictiveAbility.csv"), quote=F, row.names=F)
pdf (paste0 (outFilename, "-predictiveAbility.pdf"), width=7, height=7)
mainText = paste0 ("Predictive ability by methods.\nTrait: ", TRAIT)
boxplot(compareM, xlab="Prediction method", ylab="predictive ability", main=mainText)
dev.off()

msg ("Outputs for times for each method...")
compareM    = as.data.frame (predictedTimesList)
write.csv (compareM, paste0 (outFilename, "-times.csv"), quote=F, row.names=F)
pdf (paste0 (outFilename, "-times.pdf"), width=7, height=7)
mainText = paste0 ("Predictive ability time by methods.\nTrait: ", TRAIT)
barplot (unlist (predictedTimesList), main=mainText, ylab="Time (mins)", xlab="Methods")
dev.off()

############### Compare predicted value with diffetent methods ################
msg ("Comparing predicted value...")
predictFun <- function (method) {
	testPredicted  = bwgs.predict (predict.method=method,
						geno_train=genoTrain,pheno_train=phenoVector, geno_target=genoTarget,
						MAXNA=0.2,MAF=0.05,geno.reduct.method="NULL",reduct.size="NULL",
						r2="NULL",pval="NULL",MAP="NULL",geno.impute.method="MNI")
	return (testPredicted [,1])
}

listPREDS = mclapply (METHODS, predictFun, mc.cores=NCORES)
## put histograms on the diagonal
     panel.hist <- function(x, ...) {
         usr = par("usr"); on.exit(par(usr))
         par(usr = c(usr[1:2], 0, 1.5) )
         h      = hist(x, plot = FALSE)
         breaks = h$breaks; nB <- length(breaks)
         y      = h$counts; y <- y/max(y)
         rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
     }
## put (absolute) correlations on the upper panels, with size proportional to the correlations.
     panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
         usr = par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r   = abs(cor(x, y))
         txt = format(c(r, 0.123456789), digits = digits)[1]
         txt = paste0(prefix, txt)
         if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex.cor * r)
     }
msg ("Writing ouputs for comparisons of predicted value...")
outFilename = "out-predicted-value-comparison-methods"
ComparePRED = as.data.frame (listPREDS, col.names=METHODS)
write.csv (ComparePRED, paste0 (outFilename, ".csv"), quote=F, row.names=F)

pdf (paste0 (outFilename, ".pdf"), width=7, height=7)
mainText = paste0 ("Predicted value by methods.\nTrait ", TRAIT)
pairs(ComparePRED,lower.panel = panel.smooth, 
	  upper.panel = panel.cor, diag.panel=panel.hist, main=mainText)
dev.off()

##################### c. Compare predicted value varying training size ###############
msg ("Predicted value by size")
nSamples  = length (phenoVector)
N         = 5   # number of chunks
chunkSize = nSamples / N
sizes     = seq (chunkSize, nSamples, chunkSize)
sizes     = sizes [1:N]

cvFunSize <- function (sampleSize, method) {
	cvOutput <-bwgs.cv (sample.pop.size=sampleSize, 
					    predict.method=method,
					    geno=genoTrain, pheno=phenoVector, pop.reduct.method="RANDOM", 
					    geno.impute.method="mni", nFolds=NFOLDS, nTimes=NTIMES ) 
	return (cvOutput$cv)
}

listCVsSize = mclapply (sizes[1:(N-1)], cvFunSize, bestMethodName, mc.cores=NCORES)
listCVsSize = append (listCVsSize, list(listCVs[[bestMethodName]]$cv))
names (listCVsSize) = paste0("N=", sizes)

msg ("Writing outputs for predicted value by size...")
outFilename = "out-predicted-value-by-size-bestMethod"
CompareSize = as.data.frame (listCVsSize, col.names=paste0("N=",as.character(sizes)),check.names=F)
write.csv (CompareSize, paste0 (outFilename, ".csv"), quote=F, row.names=F)

pdf (paste0 (outFilename, ".pdf"), width=7, height=7)
mainText = paste0 ("Effect of training population size\n Trait: ", TRAIT, ", Method: ", bestMethodName )
boxplot(CompareSize,xlab="Training POP size",ylab="Predictive avility", main=mainText)
dev.off()

#################### Util functions ##################################################

