#!/usr/bin/Rscript
# Pipeline for analysis of GS models

options (width=300)
cmmArgs = commandArgs(trailingOnly = TRUE)
cmmArgs = c ("TRAIN47K-100.csv", "TARGET47K-100.csv", "YieldBLUE-tbl.csv")

USAGE="GSPipeline.R <Training genotype filename> <Target genotype filename> <Phenotype filename>"
if (length (cmmArgs) != 3) {
	message (paste0("\n\n", USAGE, "\n\n"))
	quit ()
}

source ("lglib01.R")

library(rrBLUP)
library(BGLR)
library(brnn)
library(glmnet)
library(e1071) 
library(randomForest)

library(BWGS)

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
NFOLDS  = 2
NTIMES  = 1
TRAIT   = colnames (phenoTbl)[2]

METHODS = c ("gblup", "EGBLUP", "BA", "BL", "RF", "RKHS")

########################## a.Execute prediction methods ######################

for (method in METHODS)

YieldGBLUP  <-bwgs.cv (genoTrain, phenoVector, predict.method= "gblup", 
			  		   geno.impute.method="mni", nFolds=NFOLDS, nTimes=NTIMES )
YieldEGBLUP <-bwgs.cv (genoTrain, phenoVector, predict.method= "EGBLUP",
					   geno.impute.method="mni", nFolds=NFOLDS, nTimes=NTIMES )

YieldBA     <-bwgs.cv (genoTrain, phenoVector, predict.method= "BA", 
					   geno.impute.method="mni", nFolds=NFOLDS, nTimes=NTIMES )
YieldBL     <-bwgs.cv (genoTrain, phenoVector, predict.method= "BL", 
					   geno.impute.method="mni", nFolds=NFOLDS, nTimes=NTIMES )

YieldRF     <-bwgs.cv (genoTrain, phenoVector, predict.method= "RF", 
					   geno.impute.method="mni", nFolds=NFOLDS, nTimes=NTIMES )
YieldRKHS   <-bwgs.cv (genoTrain, phenoVector, predict.method= "RKHS", 
					   geno.impute.method="mni", nFolds=NFOLDS, nTimes=NTIMES )

#################### Prediction ability comparison of methods ##################
outFilename = "out-predictive-ability-comparison-methods"
compareM    = cbind(YieldGBLUP$cv, YieldEGBLUP$cv, 
		  	        YieldBA$cv, YieldBL$cv, 
			        YieldRF$cv, YieldRKHS$cv)

colnames(compareM) = c("GBLUP","EGBLUP","BayesA","BayesL","RF","RKHS")

write.csv (compareM, paste0 (outFilename, ".csv"), quote=F, row.names=F)
pdf (paste0 (outFilename, ".pdf"), width=7, height=7)
mainText = paste0 ("Predictive ability of 6 methods for trait ", TRAIT)
boxplot(compareM, xlab="Prediction method", ylab="predictive ability", main=mainText)
dev.off()

############### Compare predicted value with diffetent methods ################
outFilename = "out-predicted-value-comparison-methods"
testGBLUP  = bwgs.predict (predict.method="GBLUP",
						geno_train=genoTrain,pheno_train=phenoVector, geno_target=genoTarget,
						MAXNA=0.2,MAF=0.05,geno.reduct.method="NULL",reduct.size="NULL",
						r2="NULL",pval="NULL",MAP="NULL",geno.impute.method="MNI")

testEGBLUP = bwgs.predict (predict.method="EGBLUP",
						geno_train=genoTrain,pheno_train=phenoVector,geno_target=genoTarget,
						MAXNA=0.2,MAF=0.05,geno.reduct.method="NULL",reduct.size="NULL",
						r2="NULL",pval="NULL",MAP="NULL",geno.impute.method="MNI")

testBayesA = bwgs.predict (predict.method="BA",
						geno_train=genoTrain,pheno_train=phenoVector,geno_target=genoTarget,
						MAXNA=0.2,MAF=0.05,geno.reduct.method="NULL",reduct.size="NULL",
						r2="NULL",pval="NULL",MAP="NULL",geno.impute.method="MNI")

testBayesL = bwgs.predict (predict.method="BL",
						geno_train=genoTrain,pheno_train=phenoVector,geno_target=genoTarget,
						MAXNA=0.2,MAF=0.05,geno.reduct.method="NULL",reduct.size="NULL",
						r2="NULL",pval="NULL",MAP="NULL",geno.impute.method="MNI")

testRF     = bwgs.predict (predict.method="RF",
						geno_train=genoTrain,pheno_train=phenoVector,geno_target=genoTarget,
						MAXNA=0.2,MAF=0.05,geno.reduct.method="NULL",reduct.size="NULL",
						r2="NULL",pval="NULL",MAP="NULL",geno.impute.method="MNI")

testRKHS   = bwgs.predict (predict.method="RKHS",
						geno_train=genoTrain,pheno_train=phenoVector,geno_target=genoTarget,
						MAXNA=0.2,MAF=0.05,geno.reduct.method="NULL",reduct.size="NULL",
						r2="NULL",pval="NULL",MAP="NULL",geno.impute.method="MNI")

# Join comparisons
ComparePRED = cbind(testGBLUP[,1], testEGBLUP[,1], 
					testBayesA[,1], testBayesL[,1], 
					testRF[,1], testRKHS[,1])


colnames(ComparePRED) = c("GBLUP", "EGBLUP", "BayesA", "BayesL", "RF", "RKHS")
write.csv (ComparePRED, paste0 (outFilename, ".csv"), quote=F, row.names=F)

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

pdf (paste0 (outFilename, ".pdf"), width=7, height=7)
mainText = paste0 ("Predicted value of 6 methods for trait ", TRAIT)
pairs(ComparePRED,lower.panel = panel.smooth, 
	  upper.panel = panel.cor, diag.panel=panel.hist, main=mainText)
dev.off()

##################### c. Compare predicted value varying training size ###############
outFilename = "out-predicted-value-by-size-bestMethod"
YieldGBLUP100 <-bwgs.cv (sample.pop.size=100, 
						 predict.method="gblup",
						 geno=genoTrain, pheno=phenoVector, pop.reduct.method="RANDOM", 
						 geno.impute.method="mni", nFolds=NFOLDS, nTimes=NTIMES ) 

YieldGBLUP300 <- bwgs.cv (sample.pop.size=300, 
						 predict.method="gblup",
						 geno=genoTrain, pheno=phenoVector, pop.reduct.method="RANDOM", 
						 geno.impute.method="mni", nFolds=NFOLDS, nTimes=NTIMES ) 

YieldGBLUP500 <- bwgs.cv (sample.pop.size=500, 
						 predict.method="gblup",
						 geno=genoTrain, pheno=phenoVector, pop.reduct.method="RANDOM", 
						 geno.impute.method="mni", nFolds=NFOLDS, nTimes=NTIMES ) 

CompareSize=cbind(YieldGBLUP100$cv, YieldGBLUP300$cv, YieldGBLUP500$cv, YieldGBLUP$cv)
colnames(CompareSize)=c("N=100","N=300","N=500","N=700")
write.csv (CompareSize, paste0 (outFilename, ".csv"), quote=F, row.names=F)

pdf (paste0 (outFilename, ".pdf"), width=7, height=7)
mainText = paste0 ("Predicte value by effect of training population size for trait ", TRAIT)

boxplot(CompareSize,xlab="Training POP size",ylab="Predictive avility",
		main=mainText)
dev.off()


