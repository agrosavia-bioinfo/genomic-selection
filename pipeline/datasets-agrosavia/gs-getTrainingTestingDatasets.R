#!/usr/bin/Rscript
# Get training and testing genotype datasets
# Input: genotype, penotype
# Ouput: training genotype, training phenotype, and testing genotype
# Check for common samples and add uncommon samples to testing genotype

source ("lglib01.R")

options (width=300)
cmdArgs = commandArgs(trailingOnly = TRUE)
cmdArgs = c ("geno.csv", "pheno.csv", "80")

genotypeFile  = cmdArgs [1]
phenotypeFile = cmdArgs [2]

#-------------------------------------------------------------
# Add label to filename and new extension (optional)
#-------------------------------------------------------------
addLabel <- function (filename, label, newExt=NULL)  {
	nameext = strsplit (filename, split="[.]")
	name    = nameext [[1]][1] 
	if (is.null (newExt))
		ext     = nameext [[1]][2] 
	else
		ext     = newExt
	newName = paste0 (nameext [[1]][1], "-", label, ".", ext )
	return (newName)
}
#-------------------------------------------------------------
# Filter by common sample names
#-------------------------------------------------------------
#filterByCommonNames <- function (genotypeFile, phenotypeFile) 
#{
	geno  = read.csv (file=genotypeFile, header=T, check.names=F, row.names=1)
	pheno = read.csv (file=phenotypeFile, header=T, check.names=F, row.names=1)

	genoSamples       = colnames (geno)[-1]
	phenoSamples      = rownames (pheno)
	commonSamples     = intersect (genoSamples, phenoSamples) 
	uncommonSamples   = setdiff (genoSamples, commonSamples)
	nTesting          = round (0.2*length (commonSamples))
	randomSamples     = sample (commonSamples, size=nTesting)
	trainingSamples   = setdiff (commonSamples, randomSamples)
	testingSamples    = c (uncommonSamples, randomSamples) 

	genoTraining   = t (geno  [,trainingSamples])
	genoTesting    = t (geno  [,testingSamples])
	phenoTraining  = cbind (Samples=trainingSamples, Gota=pheno [trainingSamples,])

	genoTrainingFile  = addLabel (genotypeFile, "TRAINING")
	genoTestingFile   = addLabel (genotypeFile, "TESTING")
	phenoTrainingFile = addLabel (phenotypeFile, "TRAINING")
	write.csv (file=genoTrainingFile, genoTraining, quote=F, row.names=T)
	write.csv (file=genoTestingFile, genoTesting, quote=F, row.names=T)
	write.csv (file=phenoTrainingFile, phenoTraining, quote=F, row.names=F)

	dim (pheno);dim(geno);dim(genoTraining);dim(genoTesting);dim(phenoTraining)
	message("")
	length(phenoSamples);length (genoSamples);length(commonSamples);length(uncommonSamples);
	nTesting;length(randomSamples);length(trainingSamples);length(testingSamples)

#	return (list (genotypeFile=genoCommonFile, phenotypeFile=phenoCommonFile))
#}
#-------------------------------------------------------------

