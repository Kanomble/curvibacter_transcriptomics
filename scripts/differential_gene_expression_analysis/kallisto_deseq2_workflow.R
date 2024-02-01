library(tximport)
library(rhdf5)
library(DESeq2)
library(DEGreport)
library(reshape2)
library(ggplot2)
library(viridis)
library(tidyr)

dir <- "data/kallisto/kallistoCounts/"
files <- list.dirs(dir)
samples <- read.csv('data/kallisto/samples.csv',header = TRUE)
target_ids <- "data/kallisto/kallistoCounts/G1_S11_L001/abundance.tsv"
cds_ids <- "data/kallisto/cds_ids.txt"

id_df <- read.table(target_ids,header=TRUE)
cds_ids <- read.table(cds_ids)
#TXNAME GENEID
tx2gene <- data.frame(TXNAME=id_df$target_id,GENEID=cds_ids$V1)

files <- file.path(dir, samples$ID, "abundance.h5")
names(files) <- samples$ID

all(file.exists(files))
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene)
#replace ID to treatment for deseq2 analysis

sampleTable <- data.frame(condition = factor(samples$treatment), sample = factor(samples$sample))

#colnames(txi.kallisto.tsv$counts) <- samples$sample
rownames(sampleTable) <- colnames(txi.kallisto.tsv$counts)

dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, sampleTable, design=~condition)

dds$sample <- factor(c("G1","G1","G1","G1","G2","G2","G2","G2","G3","G3","G3","G3",
                       "Hydra1","Hydra1","Hydra1","Hydra1",
                       "Hydra2","Hydra2","Hydra2","Hydra2",
                       "Hydra3","Hydra3","Hydra3","Hydra3"))
dds$run <- samples$sample
ddsColl <- collapseReplicates(dds, dds$sample, dds$run)
#ddsColl$condition <- factor(c("control","control","control","treatment","treatment","treatment"))


matchFirstLevel <- dds$sample == levels(dds$sample)[1]
stopifnot(all(rowSums(counts(dds[,matchFirstLevel])) == counts(ddsColl[,1])))

#pre filtering
keep <- rowSums(counts(ddsColl)) >= 10
ddsColl <- ddsColl[keep,]

#dds
dds <- DESeq(ddsColl)
res <- results(dds, contrast=c("condition","control","treated"))
#how many adjusted p-values less than 0.1:
sum(res$padj < 0.1, na.rm=TRUE)
summary(res)
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), 
          file="condition_treated_results.csv")
plotMA(res, ylim=c(-2,2))

resNorm <- lfcShrink(dds, coef=2, type="normal")

xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")

# check expression strength of specific genes
# which(rownames(res) == 'WP_087496557.1')
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
plotCounts(dds, gene=which(rownames(res) == 'WP_087494968.1'), intgroup="condition")

plotDispEsts(dds)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

plotPCA(rld, intgroup="condition")
#head(assay(vsd), 3)
#ntd <- normTransform(dds)
degPlot(dds = dds, res = res, n = 6, xs = "condition")

png(file="results/kallisto/volcano_plot_curvibacter_dataset.png",width=800, height=550)
resOrdered <- res[order(res$pvalue),]
resOrdered[["id"]] <- row.names(resOrdered)
show <- as.data.frame(resOrdered[1:10, c("log2FoldChange", "padj","id")])#"id"
degVolcano(resOrdered[,c("log2FoldChange", "padj")],plot_text = show)
dev.off()

deseq2VST <- vst(dds)
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)
sigGenes <- rownames(resOrdered[resOrdered$padj <= .05 & abs(resOrdered$log2FoldChange) > 4,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

deseq2VST_wide <- deseq2VST
deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))
deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))

heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + 
  geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_text())
heatmap

counts <- counts(dds, normalized = TRUE)
design <- as.data.frame(colData(dds))
degObj(counts, design, "results/kallisto/degObj_curvibacter.rda")
library(shiny)
shiny::runGitHub("lpantano/shiny", subdir="expression")


targets <- c('WP_087494995.1','WP_087494990.1','WP_087494996.1','WP_087494997.1')
#lipidA
targets <- c('WP_087497042.1',
             'WP_087497043.1',
             'WP_087497454.1',
             'WP_087496301.1',
             'WP_087497041.1',
             'WP_087497065.1',
             'WP_087494236.1',
             'WP_087496466.1',
             'WP_087496467.1')
#coreOligos
targets <- c('WP_087497225.1',
              'WP_232459952.1',
              'WP_232459814.1',
              'WP_087496959.1',
              'WP_087496958.1',
              'WP_087495184.1',
              'WP_087496960.1',
              'WP_087496961.1',
              'WP_087496464.1')

#epsCluster
targets <- c('WP_087496564.1',
              'WP_087496563.1',
              'WP_087496562.1',
              'WP_087496561.1',
              'WP_087496560.1',
              'WP_087496559.1',
              'WP_087496558.1',
              'WP_087496557.1',
              'WP_087496556.1',
              'WP_087496555.1')

#lacZ Operon
targets <- c(
  "WP_087497086.1","WP_087495748.1",
             "WP_087495749.1",
             "WP_087495750.1",
             "WP_087495751.1",
             "WP_087495752.1")

#glycosyltransferases cluster 3
targets <- c('WP_087494995.1','WP_087494990.1','WP_087494996.1','WP_087494997.1','WP_087494994.1','WP_087494993.1','WP_087494988.1','WP_087494987.1')
#additional genes in cluster 3
targets <- c('WP_087494994.1','WP_087494993.1','WP_087494988.1','WP_087494987.1')


#transcription factors
targets <- c("WP_087496040.1","WP_232459820.1","WP_087495013.1","WP_087497330.1","WP_232459930.1", "WP_232460005.1", "WP_087495910.1","WP_087496208.1","WP_087496231.1")


#glycosyltransferases cluster 5
targets <- c("WP_087496214.1","WP_087496217.1","WP_157673191.1","WP_087496223.1")

tcounts <- t(log2((counts(dds[targets, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(targets)+1):ncol(.))

ggplot(tcounts, aes(condition, expression, fill=condition)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales="free_y") + 
  labs(x="free living vs. host associated", 
       y="Expression (log normalized counts)", 
       fill="Control: Free Curvibacter\nTreated: Curvibacter on Hydra", 
       title="Cluster 3 additional genes")
