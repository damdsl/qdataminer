#QBiC Tuebingen - Damien Dos Santos Luis
#15/08/11

#Script for ChIPseq analysis

library("ChIPseeker", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.2")
library("BSgenome.Mmusculus.UCSC.mm10", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.2")
require("org.Mm.eg.db")
require('GenomicFeatures')
library("TxDb.Mmusculus.UCSC.mm10.knownGene", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.2")

setwd('/home/damiendsl/data_analysis_toolbox/datasets/GSE71250/macs_SRR2124925/')
files=list.files()
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene  #Reference for annotation
peaks= readPeakFile(files[[8]])
covplot(peaks, weightCol = "V5")

#Profile of ChIP binding to TSS regions

promoter = getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tagMatrix = getTagMatrix(peaks, windows = promoter)
peakHeatmap(files[[8]], TxDb = txdb, upstream = 3000, downstream = 3000, color = "red")

#Average profile of ChIP peaks binding to TSS region

plotAvgProf(tagMatrix, xlim = c(-3000, 3000), xlab = "Genomic Region (5'->3')", ylab = "Read Count Frequency")
#plotAvgProf(tagMatrix, xlim = c(-3000, 3000), conf = 0.95, resample = 500)

#Peaks annotation

peakAnno = annotatePeak(files[[8]], tssRegion = c(-3000, 3000), TxDb = txdb, annoDb="org.Mm.eg.db")
write.table(peakAnno, file = "peakAnno.txt")

source("/home/damiendsl/Rscripts/convert_tojson.R")
require('BSgenome')
require("BSgenome.Mmusculus.UCSC.mm10")
df_peakAnno=read.table('peakAnno.txt')

#load table containing Macs statistics
SRR2124925_xls=read.table(files[[10]], skip=23, header=TRUE)
#creation of metadata fields
SRR2124925_xls$sample <- rep("SRR2124925",nrow(SRR2124925_xls)) 
SRR2124925_xls$type_experiment <- rep("ChIP-seq",nrow(SRR2124925_xls)) 
SRR2124925_xls$TF_study <- rep("Hnf1b",nrow(SRR2124925_xls)) 
SRR2124925_xls$replicate <- rep(2,nrow(SRR2124925_xls)) 

df = merge(df_peakAnno, SRR2124925_xls)
write.csv(df, 'df.csv')
write.csv(df_peakAnno, 'peakAnno.csv')

json_peakAnno=toJSONarray(df_peakAnno)
write(json_peakAnno, "peakAnnotation.json")

#Finding DNA sequences corresponding to window
genome = BSgenome.Mmusculus.UCSC.mm10
seqPeaks=data.frame(chr=df_peakAnno$geneChr, start=df_peakAnno$start, end=df_peakAnno$end, strand=df_peakAnno$strand)
Peak_DNA_sequence=as.data.frame(getSeq(genome, seqPeaks$chr, start=seqPeaks$start, end=seqPeaks$end, strand=seqPeaks$strand))
write.csv(Peak_DNA_sequence, 'Peak_DNA_sequence.csv')
seqGenes=data.frame(chr=df_peakAnno$geneChr, start=df_peakAnno$geneStart, end=df_peakAnno$geneEnd)
Gene_sequence=getSeq(genome, seqGenes$chr, start=seqGenes$start, end=seqGenes$end)

pdf(file = 'peekSeeker_report.pdf', width=8.50, height=11)
par(mfrow=(c(1,4)))
#ChIP peaks coverage plot
#This value can change depending on the species
#V5 has been controlled for mice and humans - V4 for Drosophila
covplot(peaks, weightCol = "V5")
#Average profile of ChIP peaks binding to TSS region
plotAvgProf(tagMatrix, xlim = c(-3000, 3000), xlab = "Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAnnoPie(peakAnno)
plotDistToTSS(peakAnno, title = "Distribution of transcrition factor-binding loci \n relative to TSS")
dev.off()
