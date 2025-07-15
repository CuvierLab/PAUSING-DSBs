######################################################################################################################################################
######################################################################################################################################################
# Cuvier's Lab
# David Depierre - Thomas De Freitas
######################################################################################################################################################
######################################################################################################################################################


# Load LIBRARY
library(gplots)
'%ni%' = Negate('%in%')
require(GenomicRanges)
require(BiocGenerics)
require(parallel)


workdir = "" 

profdir = paste0(workdir, "PROFILE_MATRIX/")
outdir = paste0(workdir, "QUANTIF/")

# QUANTIFICATION OF READS ON A GIVEN WINDOW
# 500 = TSS (position 0)
# 1 = 10 bp

# QUANTIFICATION ON GENE BODY
readSumWindow = function(prof, start=400, end=800){
	Nameprof = deparse(substitute(prof))
	print(Nameprof)
	prof = rowSums(prof[,c(start:end)])
	filenameRDS = paste0(outdir, "Q_", Nameprof ,  "_readsCounts_GB_SCALED.RDS")
	saveRDS(prof, file=filenameRDS)
}

# LOAD DATA 
chipseq_data_profmat = readRDS(paste0(profdir, "chipseq_data_profmat.RDS"))

# RUN AND SAVE QUANTIFICATION
readSumWindow(chipseq_data_profmat)

