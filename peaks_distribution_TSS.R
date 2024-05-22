#Load libraries#
library(ChIPseeker)
library (ensembldb)
#Set Working directory to where the files are: these are the file of ChIP-seq peaks retrieved from MACS#
setwd("C:/OsM2 manuscript as of 25.3.24/ChIP/peaks distribution")
#Read the IP over mock peak files#
IP_B5 <- readPeakFile("OsM2_B5_peaks.xls")
IP_B6 <- readPeakFile("OsM2_B6_peaks.xls")
#Read the mock over IP peak files# they are in a different directory#
setwd("c:/OsM2 manuscript as of 25.3.24/ChIP/peaks distribution/Control")
Mock_B5 <- readPeakFile("OsM2_B5_peaks.xls")
Mock_B6 <- readPeakFile("OsM2_B6_peaks.xls")

#Make Txdb from the GTF file# From this Txdb you will get the coordinates of the defined promotor sequences#
dbf <- ensDbFromGtf("C:/OsM2 manuscript as of 25.3.24/ChIP/peaks distribution/Oryza_sativa.IRGSP-1.0.53.gtf.gz")
txdb <- EnsDb(dbf)
ens_txdb <-txdb
promoters <- getPromoters(TxDb= ens_txdb, upstream=3000, downstream=3000,
                         by = "transcript")

#Confirm that names of seqlevels are matched in the peak files and in the promotors#
#if they are not matched, then rename either one of them to match them#
promoters@seqnames
IP_B5@seqnames
promoters <- renameSeqlevels(promoters, c ("Chr1", "Chr11", "Chr12", "Chr8", "Chr6",
                                           "Chr3", "Chr10", "Chr4", "Chr2", "Chr7",
                                           "Chr9", "Chr5", "Mt", "Pt"))
#Create your tag matrices# 
#ignore the warning due to presence of ChrSy and ChrUn are sedumolecules and has no Loc_IDs#
tagMatrix1 <- getTagMatrix(IP_B5, windows=promoters)
tagMatrix2 <- getTagMatrix(IP_B6, windows=promoters)
tagMatrix3<- getTagMatrix(Mock_B5, windows=promoters)
tagMatrix4 <- getTagMatrix(Mock_B6, windows=promoters)

tagMatrix = list (tagMatrix1, tagMatrix2, tagMatrix3, tagMatrix4)

names(tagMatrix) <- c ("IP B5", "IP B6", "Mock B5", "Mock B6")

plotAvgProf(tagMatrix , xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
