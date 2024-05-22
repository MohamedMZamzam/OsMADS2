# OsMADS2
Custom R script to reproduce results in manuscript Zamzam et al. (2023). [Unequal genetic redundancy among the rice transcription factors OsMADS2 and OsMADS4 reveals distinct roles in floret organ development. bioRxiv 2023.08.05.552136](https://doi.org/10.1101/2023.08.05.552136)
## ChIP-Seq data analysis
1- The provided custom R script in the file [peaks_distribution_TSS.R](https://github.com/MohamedMZamzam/OsMADS2/blob/main/peaks_distribution_TSS.R) demostartes the steps of plotting distribution of enriched ChIP-Seq peaks with respect to the TSS and the flanks (-3000 to +3000).

2- The provided custom R script in the file [Defining_confident_peaks_and_Motif_analysis](https://github.com/MohamedMZamzam/OsMADS2/blob/main/Defining_confident_peaks_and_Motif_analysis.R) deomnstrate the steps of handling the OsMADS2 enriched peaks retrevied  by MACS2 from the two biological replicate for defining the confident peak list and for extracting the DNAsets of 500 bp flanking the confident peak medpoints that were used as input to the MEME analysis. The file also provides the code for visualizing the distribution of the enriched motifs in peak region and the flanking 3000 bp.



## Downloads
1- The GTF file <b>Oryza_sativa.IRGSP-1.0.53.gtf.gz</b> for creating the TxDb can be dowloaded from [here](https://plants.ensembl.org/info/data/ftp/index.html)



## Contributors
### Mohamed Zamzam
### Raghavaram Peespati 
