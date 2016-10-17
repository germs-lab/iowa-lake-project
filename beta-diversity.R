#beta diversity
#rm(list=ls())

#setwd("~/")

##########   Load library   ############
library(phyloseq)
library(ape)
library(DESeq2)
library("ggplot2")
library("plyr")
library("vegan")

##########   read data   ############
#http://joey711.github.io/phyloseq-demo/import-biom-sd-example.html
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')
library(phyloseq)
biom_file <- "~/otu_table_mc2_w_tax.biom"
biom_file <- "~/Box Sync/2016/5May/16s/otu_table_mc2_w_tax.biom"
biom_otu_tax <- import_biom(biom_file)

# read map file
#meta  <- read.csv("~/meta_w_new_group.csv",header=T, row.names=1)
meta  <- read.csv("~/Box Sync/2016/6June/16s/meta_w_new_group.csv",header=T, row.names=1)
# Check that the rownames match the sample names
all(sample_names(biom_otu_tax) %in% rownames(meta) )
# Convert to "sample_data" class
sampledata = sample_data(meta)

#tree file
library(ape)
#tree = read.tree("~/rep_set.tre")
tree = read.tree("~/Box Sync/2016/5May/16s/rep_set.tre")

# merge
physeq2 = merge_phyloseq(sampledata, biom_otu_tax, tree)

## remove outlier
physeq3 = prune_samples(sample_names(physeq2) != "2015000001", physeq2)
physeq4 = prune_samples(sample_names(physeq3) != "2014195002", physeq3)
#physeq4 = prune_samples(sample_names(physeq4) != "2014216006", physeq4)
#physeq4 = prune_samples(sample_names(physeq4) != "2014209005", physeq4)

##########  Normalizatoin   ############
library(DESeq2)
countData <- as.data.frame(as.matrix(otu_table(physeq4)))
colData <- as.data.frame(sample_data(physeq4))

## Create a DESeqDataSet object
deseq2_deseqds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ new_group)

## Differential expression analysis.
#deseq2_de <- DESeq(deseq2_deseqds)
deseq2_de <- DESeq(deseq2_deseqds, fitType = "local")

count <- counts(deseq2_de, normalized=TRUE)
otu_table(physeq4) <- otu_table(count, taxa_are_rows=TRUE)

#saveRDS(physeq4, "~/Box Sync/Github/germs-lab/iowa_lake_project/water-project-objext.rds")
physeq4 <- readRDS("~/Box Sync/Github/germs-lab/iowa_lake_project/water-project-objext.rds")


#plot
qiime.ord <- ordinate(physeq4, "PCoA","unifrac",weighted = TRUE)
p1 = plot_ordination(physeq4, qiime.ord, color="new_group")
p1 + stat_ellipse()
cor = qiime.ord$vectors
met = sample_data(physeq4)[,c(4,5,6,7,8,9,10,11,12,13,14)]
datEF = envfit(cor, met)
datEF.df = as.data.frame(datEF$vectors$arrows*sqrt(datEF$vectors$r))
datEF.df$species <- rownames(datEF.df)
p1+ geom_segment(data = datEF.df, aes(x=0, xend=Axis.1, y=0, yend=Axis.2), arrow = arrow(length=unit(0.2, "cm")), colour="grey")+geom_text(data= datEF.df, aes(x= Axis.1, y=Axis.2, label=species), size=5,colour="black") + stat_ellipse()




pdf("pcoa.pdf", width=6, height=6)
print(p1)
dev.off()

#NMDS
GP.ord <- ordinate(physeq4, "NMDS", "bray")
p1 = plot_ordination(physeq4, GP.ord, color="new_group")
p1 + stat_ellipse()

pdf("NMDS.pdf", width=6, height=6)
print(p1)
dev.off()

library("ggplot2")
library("plyr")
library("vegan")
#d = phyloseq::distance(physeq4, method="unifrac")
d = phyloseq::distance(physeq4, method="bray")
df <- data.frame(sample_data(physeq4))
ado <- adonis(d ~ new_group, data=df)
capture.output(ado, file="adonis.txt")

b_result <- betadisper(d, df$new_group, "median")
ano_re <- anova(b_result)
capture.output(ano_re, file="anova.txt")
tuk <- TukeyHSD(b_result)
capture.output(tuk, file="tukeyHSD.txt")

# Plot dispersions , this is why I choose local
#png("qc-dispersions.png", 1000, 1000, pointsize=20)
pdf("qc-dispersions.pdf", width=6, height=6)
#des <- estimateDispersions(deseq2_de, fitType="parametric")
des <- estimateDispersions(deseq2_de, fitType="local")
#plotDispEsts(des)
plotDispEsts(deseq2_de, main="Dispersion plot")
dev.off()

##PCA plot from DESeq2
rld <- rlog(deseq2_de, blind=FALSE)
head(assay(rld))
hist(assay(rld))
pcadata <- plotPCA(rld,intgroup=c("new_group", "new_group"))
pdf("pcoa_deseq2.pdf", width=6, height=6)
pcadata
dev.off()