.libPaths(c("/home/shenlm/R/x86_64-pc-linux-gnu-library/4.1","/opt/R/4.1.2/lib/R/library"))
#input
library(getopt)
library(ggbiplot)
spec=matrix(c('work_dir','w',1,"character",
              'peak_dir','p',1,"character",
              'species','s',1,"character"),
              byrow=TRUE,ncol=4)
opt=getopt(spec)
setwd(opt$work_dir)
inputDir1 <- opt$peak_dir
species <- opt$species
scan="FALSE"
#library
library(clusterProfiler)
library(ChIPseeker)
library(ggplot2)
library(tidyverse)
if (species=="human") {
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  require(org.Hs.eg.db)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  org <- org.Hs.eg.db
  ori <- "hsa"
}else if (species=="mouse") {
  require(TxDb.Mmusculus.UCSC.mm10.knownGene)
  require(org.Mm.eg.db)
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  org <- org.Mm.eg.db
  ori <- "mmu"
}else if (species=="melanogaster") {
  require(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
  require(org.Dm.eg.db)
  txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
  org <- org.Dm.eg.db
  ori <- "Dm"
}

#readin
paths1 <- list.files(inputDir1,full.names = T)
pathsToInputFiles1 <- paths1[grep(pattern = "*_tochipseeker.txt", paths1)]
pathsToInputFiles <- c(pathsToInputFiles1)
names <- sub("_tochipseeker.txt","",basename(pathsToInputFiles1))
chip <- lapply(pathsToInputFiles1, function(x)
  readPeakFile(x))
names(chip) <- names


#tss enrichment
promoter <- getPromoters(TxDb = txdb,upstream = 3000,downstream = 3000)
tagMatrix <- lapply(chip, function(x)
  getTagMatrix(x,windows = promoter))
png("1tss_enrichment_3000.png")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()

promoter <- getPromoters(TxDb = txdb,upstream = 300,downstream = 300)
tagMatrix <- lapply(chip, function(x)
  getTagMatrix(x,windows = promoter))
png("1tss_enrichment_300.png")
plotAvgProf(tagMatrix, xlim=c(-300, 300),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()

#annotation_bar_plot
peakAnnoList <- lapply(chip, annotatePeak,
                   tssRegion = c(-3000,3000),TxDb = txdb)
p <- plotAnnoBar(peakAnnoList)
ggsave(plot = p,"2feature_distribution.pdf",width = 16,height = 16)
for (ii in 1:length(peakAnnoList)) {
  p <- upsetplot(peakAnnoList[[ii]])
  ggsave(plot = p,paste0("2_",names[ii],"_upset.pdf"),width = 20,height = 10)
}

#KEGG and go
genes <- lapply(chip, function(x)
  seq2gene(x,c(-1000,3000),3000,TxDb = txdb))
gene2 <- list()
for (i in 1:length(genes)) {
  tmp <- genes[[i]]
  gene2[[i]] <- tmp[which(str_detect(tmp,"ENS")==FALSE)]
  names(gene2)[i] <- names[i]
}
#KEGG
cc <- compareCluster(geneClusters = gene2,fun = "enrichKEGG",organism=ori)
pdf("3KEGG.pdf",height = 20,width = length(chip)*5)
dotplot(cc,showCategory=20)
dev.off()

#go
go <- compareCluster(geneClusters = gene2,fun = "enrichGO",OrgDb=org)
pdf("4go.pdf",height = 20,width = length(chip)*5)
dotplot(go,showCategory=20)
dev.off()

