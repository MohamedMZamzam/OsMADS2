#Set Working directory#
setwd("C:/transgenic_plants_ analyses 25.3.24 original transferred/osm2_ko/Received from Raghavram/ChIPSeq/macs3_keep_dup_AUTO/macs3_keep_dup_AUTO")
#Load libraries#
library(GenomicFeatures)
library(GenomicRanges)
library (tidyverse)
library (officer)
library(seqinr)
library (universalmotif)
library(Biostrings)
library(BSgenome.Osativa.MSU.MSU7)
library (ChIPpeakAnno)
library(plyranges)
# load and filter the dataset of ChIP-Seq peaks (enriched bound sites) for Replicate 1 (B5)#
# B5 is the dataset for the peaks enriched in OsMADS2 IP replicate 1 over mock IP replicate 1#
B5 <- read.table ("OsM2_B5_peaks.xls")
B5 <- B5 [, c (1:3, 7:9)]
#FILTER WITH QVALUE OF 0.05 and foldenrichment of 2#
B5<- dplyr::filter (B5, B5$V7 >= 1.30103, B5$V9 >= 1.30103, B5$V8 >=2)
B5 <- B5 [1:3]
Names <- c ("Chr", "start", "end")
B5 <- `colnames<-`(B5, Names)

# B5_Mock6 is the dataset for the peaks enriched in OsMADS2 IP replicate 1 over mock IP replicate 2#

B5_Mock6 <- read.table ("OsM2_B5with6_peaks.xls")
B5_Mock6 <- B5_Mock6 [, c (1:3, 7:9)]
B5_Mock6 <- filter (B5_Mock6, B5_Mock6$V7 >= 1.30103, B5_Mock6$V9 >= 1.30103, B5_Mock6$V8 >=2)
B5_Mock6 <- B5_Mock6 [1:3]
B5_Mock6 <- `colnames<-`(B5_Mock6, Names)

#Prepare the Genomic ranges for B5 and for B5_Mock6#
B5_GR<- makeGRangesFromDataFrame(B5,
                                 keep.extra.columns=FALSE,
                                 ignore.strand=FALSE,
                                 seqinfo=NULL,
                                 seqnames.field=c("seqname",
                                                  "chromosome", "chrom",
                                                  "chr", "chromosome_name",
                                                  "seqid"),
                                 start.field="start",
                                 end.field=c("end", "stop"),
                                 strand.field="strand",
                                 starts.in.df.are.0based=FALSE)

B5_Mock6_GR<- makeGRangesFromDataFrame(B5_Mock6,
                                       keep.extra.columns=FALSE,
                                       ignore.strand=FALSE,
                                       seqinfo=NULL,
                                       seqnames.field=c("seqname",
                                                        "chromosome", "chrom",
                                                        "chr", "chromosome_name",
                                                        "seqid"),
                                       start.field="start",
                                       end.field=c("end", "stop"),
                                       strand.field="strand",
                                       starts.in.df.are.0based=FALSE)

#Find the overlapping peaks for the two datasets of replicate 1 (B5 and B5_Mock)#
Overlap_B5 <- findOverlapsOfPeaks(B5_GR, B5_Mock6_GR)
#retreive the GRanges constituting all merged overlapping peaks for replicate 1#
Overlap_B5_mergedpeaks <- Overlap_B5$mergedPeaks
View (as.data.frame (Overlap_B5_mergedpeaks))
# load and filter the dataset of ChIP-Seq peaks (enriched bound sites) for Replicate 2 (designated as B6)#
B6 <- read.table ("OsM2_B6_peaks.xls")
B6 <- B6 [, c (1:3, 7:9)]
B6 <- filter (B6, B6$V7 >= 1.30103, B6$V9 >= 1.30103, B6$V8 >=2)
B6 <- B6 [1:3]
B6 <- `colnames<-`(B6, Names)
#Prepare the GRanges for B6#
B6_GR <- makeGRangesFromDataFrame(B6,
                                  keep.extra.columns=FALSE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field=c("seqname",
                                                   "chromosome", "chrom",
                                                   "chr", "chromosome_name",
                                                   "seqid"),
                                  start.field="start",
                                  end.field=c("end", "stop"),
                                  strand.field="strand",
                                  starts.in.df.are.0based=FALSE)


#Find the overlapping peaks between the merged peaks of replicate 1 and the peaks of replicate 2#

Overlap_B5_B6 <- findOverlapsOfPeaks(Overlap_B5_mergedpeaks_GR, B6_GR)
#retreive the GRanges constituting all merged overlapping peaks between R1 AND R2#
Overlap_B5_B6_Mergedpeaks <- Overlap_B5_B6$mergedPeaks
View (as.data.frame (Overlap_B5_B6_Mergedpeaks))

#Find the confident peaks in B5 dataset by filtering by overlap B5 GR against the merged overlapped peaks from R1 and R2#

B5_confident_peaks <- filter_by_overlaps(B5_GR, as_granges(Overlap_B5_B6_Mergedpeaks))
B5_confident_peaks <- as.data.frame(B5_confident_peaks)
# merge the other columns from the B5 dataset to this list by left joining#
B5 <- read.table ("OsM2_B5_peaks.xls")
B5 <- B5 [, c (1:3, 5, 7:9)]
#FILTER WITH QVALUE OF 0.05 and foldenrichment of 2#
B5<- filter (B5, B5$V7 >= 1.30103, B5$V9 >= 1.30103, B5$V8 >=2)
B5 <- B5 [1:4]
Names <- c ("Chr", "start", "end", "abs_summit")
B5 <- `colnames<-`(B5, Names)
B5$start <- as.integer(B5$start)

B5_confident_peaks <- left_join(B5_confident_peaks, B5, by = "start")

#Find the confident peaks in B5_Mock6 dataset by filtering by overlap B5_Mock6 GR against the merged overlapped peaks from R1 and R2#
B5_Mock6_confident_peaks <- filter_by_overlaps(B5_Mock6_GR, as_granges(Overlap_B5_B6_Mergedpeaks))

B5_Mock6_confident_peaks <- as.data.frame(B5_Mock6_confident_peaks)  
# merge the other columns from the B5_mock6 dataset to this list by left joining#
B5_Mock6 <- read.table ("OsM2_B5with6_peaks.xls")
B5_Mock6 <- B5_Mock6 [, c (1:3, 5, 7:9)]
#FILTER WITH QVALUE OF 0.05 and foldenrichment of 2#
B5_Mock6<- filter (B5_Mock6, B5_Mock6$V7 >= 1.30103, B5_Mock6$V9 >= 1.30103, B5_Mock6$V8 >=2)
B5_Mock6 <- B5_Mock6 [1:4]
Names <- c ("Chr", "start", "end", "abs_summit")
B5_Mock6 <- `colnames<-`(B5_Mock6, Names)
B5_Mock6$start <- as.integer(B5_Mock6$start)

B5_Mock6_confident_peaks <- left_join(B5_Mock6_confident_peaks, B5_Mock6, by = "start")

#Find the confident peaks in B6 dataset by filtering by overlap B6 GR against the merged overlapped peaks from R1 and R2#

B6_confident_peaks <- filter_by_overlaps(B6_GR, as_granges(Overlap_B5_B6_Mergedpeaks))
B6_confident_peaks <- as.data.frame(B6_confident_peaks)
# merge the other columns from the B6 dataset to this list by left joining#
B6 <- read.table ("OsM2_B6_peaks.xls")
B6 <- B6 [, c (1:3, 5, 7:9)]
#FILTER WITH QVALUE OF 0.05 and foldenrichment of 2#
B6<- filter (B6, B6$V7 >= 1.30103, B6$V9 >= 1.30103, B6$V8 >=2)
B6 <- B6 [1:4]
Names <- c ("Chr", "start", "end", "abs_summit")
B6 <- `colnames<-`(B6, Names)
B6$start <- as.integer(B6$start)

B6_confident_peaks <- left_join(B6_confident_peaks, B6, by = "start")


#Motif analysis#
## 1. EXTRACTING 250 bp flanking the SUMMITS for MEME analysis##
B5_confident_peaks$startt <- as.numeric (B5_confident_peaks$abs_summit) -250
B5_confident_peaks$endd <- as.numeric (B5_confident_peaks$abs_summit) +249
B5_confident_peaks_250_summ <- B5_confident_peaks [, c (1, 9:10)]

names <- 2:1600
rownames (B5_confident_peaks_250_summ) <- names

B5_flankingsummit<- makeGRangesFromDataFrame(B5_confident_peaks_250_summ,
                                keep.extra.columns=FALSE,
                                ignore.strand=FALSE,
                                seqinfo=NULL,
                                seqnames.field=c("seqnames", "seqname",
                                                 "chromosome", "chrom",
                                                 "chr", "chromosome_name",
                                                 "seqid"),
                                start.field="startt",
                                end.field=c("endd", "stop"),
                                strand.field="strand",
                                starts.in.df.are.0based=FALSE)

seqsez <- getSeq (Osativa, B5_flankingsummit)

writeXStringSet (seqsez, "B5_confident_peaks")

## 2. EXTRACTING 250 bp flanking the peak midpoit for MEME analysis##
B5_confident_peaks$startt <- ((B5_confident_peaks$start + B5_confident_peaks$end.x)/2) - 250
B5_confident_peaks$endd <- B5_confident_peaks$startt + 499
B5_confident_peaks_250_mid <- B5_confident_peaks [, c (1, 9:10)]

names <- 2:1600
rownames (B5_confident_peaks_250_mid) <- names

B5_flanking_peakmidpoint<- makeGRangesFromDataFrame(B5_confident_peaks_250_mid,
                                keep.extra.columns=FALSE,
                                ignore.strand=FALSE,
                                seqinfo=NULL,
                                seqnames.field=c("seqnames", "seqname",
                                                 "chromosome", "chrom",
                                                 "chr", "chromosome_name",
                                                 "seqid"),
                                start.field="startt",
                                end.field=c("endd", "stop"),
                                strand.field="strand",
                                starts.in.df.are.0based=FALSE)

seqsez <- getSeq (Osativa, B5_flanking_peakmidpoint)

writeXStringSet (seqsez, "B5_confident_peaks_mid")

## 3. EXTRACTING 250 bp flanking the peak summit for MEME analysis##
B5_Mock6_confident_peaks$startt <- as.numeric (B5_Mock6_confident_peaks$abs_summit) -250
B5_Mock6_confident_peaks$endd <- as.numeric (B5_Mock6_confident_peaks$abs_summit) +249
B5_Mock6_confident_peaks_250_summ <- B5_Mock6_confident_peaks [, c (1, 9:10)]

names <- 2:1635
rownames (B5_Mock6_confident_peaks_250_summ) <- names

B5_Mock6_flankingsummit<- makeGRangesFromDataFrame(B5_Mock6_confident_peaks_250_summ,
                                keep.extra.columns=FALSE,
                                ignore.strand=FALSE,
                                seqinfo=NULL,
                                seqnames.field=c("seqnames", "seqname",
                                                 "chromosome", "chrom",
                                                 "chr", "chromosome_name",
                                                 "seqid"),
                                start.field="startt",
                                end.field=c("endd", "stop"),
                                strand.field="strand",
                                starts.in.df.are.0based=FALSE)

seqsez <- getSeq (Osativa, B5_Mock6_flankingsummit)

writeXStringSet (seqsez, "B5_mock6_confident_peaks")

## 4. EXTRACTING 250 bp flanking the peak midpoint for MEME analysis##
B5_Mock6_confident_peaks$startt <- ((B5_Mock6_confident_peaks$start + B5_Mock6_confident_peaks$end.x)/2) - 250
B5_Mock6_confident_peaks$endd <- B5_Mock6_confident_peaks$startt + 499
B5_Mock6_confident_peaks_250_mid <- B5_Mock6_confident_peaks [, c (1, 9:10)]

names <- 2:1635
rownames (B5_Mock6_confident_peaks_250_mid) <- names

B5_mock5_flanking_peakmidpoint<- makeGRangesFromDataFrame(B5_Mock6_confident_peaks_250_mid,
                                keep.extra.columns=FALSE,
                                ignore.strand=FALSE,
                                seqinfo=NULL,
                                seqnames.field=c("seqnames", "seqname",
                                                 "chromosome", "chrom",
                                                 "chr", "chromosome_name",
                                                 "seqid"),
                                start.field="startt",
                                end.field=c("endd", "stop"),
                                strand.field="strand",
                                starts.in.df.are.0based=FALSE)

seqsez <- getSeq (Osativa, B5_mock5_flanking_peakmidpoint)

writeXStringSet (seqsez, "B5_Mock6_confident_peaks_mid")



## 5. EXTRACTING 250 bp flanking the peak summit for MEME analysis##
B6_confident_peaks$startt <- as.numeric (B6_confident_peaks$abs_summit) -250
B6_confident_peaks$endd <- as.numeric (B6_confident_peaks$abs_summit) +249
B6_confident_peaks_250_summ <- B6_confident_peaks [, c (1, 9:10)]

names <- 2:1677
rownames (B6_confident_peaks_250_summ) <- names

B6_flankingsummit<- makeGRangesFromDataFrame(B6_confident_peaks_250_summ,
                                keep.extra.columns=FALSE,
                                ignore.strand=FALSE,
                                seqinfo=NULL,
                                seqnames.field=c("seqnames", "seqname",
                                                 "chromosome", "chrom",
                                                 "chr", "chromosome_name",
                                                 "seqid"),
                                start.field="startt",
                                end.field=c("endd", "stop"),
                                strand.field="strand",
                                starts.in.df.are.0based=FALSE)

seqsez <- getSeq (Osativa, B6_flankingsummit)

writeXStringSet (seqsez, "B6_confident_peaks")

## 5. EXTRACTING 250 bp flanking the peak midpoint for MEME analysis##
B6_confident_peaks$startt <- ((B6_confident_peaks$start + B6_confident_peaks$end.x)/2) - 250
B6_confident_peaks$endd <- B6_confident_peaks$startt + 499
B6_confident_peaks_250_mid <- B6_confident_peaks [, c (1, 9:10)]

names <- 2:1677
rownames (B6_confident_peaks_250_mid) <- names

B6_flankingpeakmidpoint<- makeGRangesFromDataFrame(B6_confident_peaks_250_mid,
                                keep.extra.columns=FALSE,
                                ignore.strand=FALSE,
                                seqinfo=NULL,
                                seqnames.field=c("seqnames", "seqname",
                                                 "chromosome", "chrom",
                                                 "chr", "chromosome_name",
                                                 "seqid"),
                                start.field="startt",
                                end.field=c("endd", "stop"),
                                strand.field="strand",
                                starts.in.df.are.0based=FALSE)

seqsez <- getSeq (Osativa, B6_flankingpeakmidpoint)

writeXStringSet (seqsez, "B6_confident_peaks_mid")


#plotting the distribution of the enriched motifs retrieved by MEME around the peak summit#
B5_confident_peaks$startt <- as.numeric (B5_confident_peaks$abs_summit) -3000
B5_confident_peaks$endd <- as.numeric (B5_confident_peaks$abs_summit) +2999
B5_confident_peaks_3000_summ <- B5_confident_peaks [, c (1, 9:10)]

names <- 2:1600
rownames (B5_confident_peaks_3000_summ) <- names

B5_flankingsummit_3000<- makeGRangesFromDataFrame(B5_confident_peaks_3000_summ,
                                keep.extra.columns=FALSE,
                                ignore.strand=FALSE,
                                seqinfo=NULL,
                                seqnames.field=c("seqnames", "seqname",
                                                 "chromosome", "chrom",
                                                 "chr", "chromosome_name",
                                                 "seqid"),
                                start.field="startt",
                                end.field=c("endd", "stop"),
                                strand.field="strand",
                                starts.in.df.are.0based=FALSE)

B5_flankingsummit_3000_DNAset <- getSeq (Osativa, B5_flankingsummit_3000)


#plotting the distribution of the enriched motifs retrieved by MEME around the peak midpoit#
B5_confident_peaks$startt <- ((B5_confident_peaks$start + B5_confident_peaks$end.x)/2) - 3000
B5_confident_peaks$endd <- B5_confident_peaks$startt + 5999
B5_confident_peaks_3000_mid <- B5_confident_peaks [, c (1, 9:10)]

names <- 2:1600
rownames (B5_confident_peaks_3000_mid) <- names

B5_peakmidpoint_3000<- makeGRangesFromDataFrame(B5_confident_peaks_3000_mid,
                                keep.extra.columns=FALSE,
                                ignore.strand=FALSE,
                                seqinfo=NULL,
                                seqnames.field=c("seqnames", "seqname",
                                                 "chromosome", "chrom",
                                                 "chr", "chromosome_name",
                                                 "seqid"),
                                start.field="startt",
                                end.field=c("endd", "stop"),
                                strand.field="strand",
                                starts.in.df.are.0based=FALSE)

B5_peakmidpoint_3000_DNAset <- getSeq (Osativa, B5_peakmidpoint_3000)

#Function for plotting the motif distribution#

DistributionOfMotifInFasta <- function(Motif, Fasta)
{
  patt <- vmatchPattern(Motif,  Fasta, fixed = "subject")
  patt <- as.data.frame(patt)
  
  Toplot <<- matrix(0,15,2)
  
  a = 0
  b = 400
  for(i in 1:15)
  {
    
    dis_3500_4000 <- filter(patt, patt$start > a, patt$start <= b)
    Toplot[i,1] <<- a-3000
    
    Toplot[i,2] <<- nrow(dis_3500_4000)
    
    a = a+ 400
    b = b + 400
  }
  plot(Toplot[,1],
       Toplot[,2],
       xaxt = "n",
       type = "b",
       col = "red",
       xlab = "Peak and 3kb Flank from Center",
       ylab = "Count of the Motif ",main = paste0("Motif = " ,Motif))
  axis(1, at = seq(-3000, 2999, by = 400), las = 2)
}


DistributionOfMotifInFasta(Motif = "GCGGCGGCGGCGRC", Fasta = B5_flankingsummit_3000_DNAset)
DistributionOfMotifInFasta(Motif = "GCGGCGGCGGCGRC", Fasta = B5_peakmidpoint_3000_DNAset)

DistributionOfMotifInFasta(Motif = "TATACAAAGTTTGYA", Fasta = B5_flankingsummit_3000_DNAset)
DistributionOfMotifInFasta(Motif = "TATACAAAGTTTGYA", Fasta = B5_peakmidpoint_3000_DNAset)


