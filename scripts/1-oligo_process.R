library(oligo) # oligo 
library(pd.hugene.2.0.st) # array info library
library(hugene20sttranscriptcluster.db) # annotation package
# library(limma) # diffrential expression analysis

# read data

setwd("/home/panxq/projects/BME1063_HW2/data")
celFiles <- list.celfiles()

theData <- data.frame(Key=rep(c("control", "A", "B"), each=3))
rownames(theData) <- basename(celFiles)
lvls <- c("exprs", "_ALL_")
vMtData <- data.frame(channel=factor('exprs', levels=lvls), 
                      labelDescription="Sample type")

pd <- new("AnnotatedDataFrame", data=theData, varMetadata=vMtData)

affyRaw <- read.celfiles(celFiles, pkgname="pd.hugene.2.0.st", phenoData=pd)

# modify sample name
sns <- sampleNames(affyRaw)
sns <- gsub('\\.CEL$', '', sns)
sampleNames(affyRaw) <- sns
# rm(sns, celFiles)
# It would be much more convenient if we can extract metadata from 
# file names. However, the files are named in a stupid way.

# QC
# MA plot
svg('../FinalReport/images/MAplot.svg')
MAplot(affyRaw)
# box plot
svg('../FinalReport/images/boxplot.svg')
oligo::boxplot(affyRaw, 'all')

# Normalization
# the normalization is done by RMA

eset <- rma(affyRaw)

# MA plot after normalization
pdf('../FinalReport/images/MAplot_after_rma.pdf')
MAplot(eset)
# box plot after Normalization
svg('../FinalReport/images/boxplot_after_rma.svg')
oligo::boxplot(eset)

# save data

write.csv(exprs(eset), file='../FinalReport/result/data.csv')
write.table(exprs(eset), file="../FinalReport/result/data.txt", sep='\t')

# add gene annotation

my_frame <- data.frame(exprs(eset))
Annot <- data.frame(
    ACCNUM=sapply(contents(hugene20sttranscriptclusterACCNUM), paste, collapse=", "), 
    SYMBOL=sapply(contents(hugene20sttranscriptclusterSYMBOL), paste, collapse=", "), 
    DESC=sapply(contents(hugene20sttranscriptclusterGENENAME), paste, collapse=", ")
)

all <- merge(Annot, my_frame, by.x=0, by.y=0, all=T)
write.csv(all,file="../FinalReport/result/data.ann.csv")
write.table(all, file="../FinalReport/result/data.ann.txt", sep='\t')

# differential expression analysis

# design <- model.matrix(~factor(eset[["Key"]]))
# fit <- lmFit(eset, design)
# ebayes <- eBayes(fit)
# lod <- -log10(ebayes[["p.value"]][,2])
# mtstat<- ebayes[["t"]][,2]