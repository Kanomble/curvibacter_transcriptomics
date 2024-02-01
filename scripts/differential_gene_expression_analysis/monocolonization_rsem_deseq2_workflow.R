library(DESeq2)
library(tximport)

experiments <- read.csv("data/CurvibacterMonocolonization/samples.csv")
files <- file.path("data/CurvibacterMonocolonization/rsem_results/",
                   paste0(experiments$id,".genes.results"))

names(files) <- experiments$id
all(file.exists(files))
sampleTable <- data.frame(condition = factor(experiments$treatment))


id_df <- read.table(files[1], header=TRUE)
tx2gene <- data.frame(TXNAME=id_df$gene_id,GENEID=id_df$gene_id)

txi.rsem <- tximport(files, type="rsem",tx2gene = tx2gene, txIn = FALSE, txOut = FALSE)
txi.rsem$counts <- txi.rsem$counts + 1
txi.rsem$abundance <- txi.rsem$abundance + 1
txi.rsem$length <- txi.rsem$length + 1
rownames(sampleTable) <- colnames(txi.rsem$counts)


# recommended tximport vignette cmd
threshold <- 10
dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)
keep <- rowSums(counts(dds)) >= threshold
dds <- dds[keep,]


rld <- rlog(dds, blind=TRUE)
cat(paste("[*] PRODUCING PCA",file.path(outp,paste0("curvibacter_mono_pca.png")),"\n"))
png(file=file.path(outp,paste0("curvibacter_mono_pca.png")), width=800, height=550)
pca <- plotPCA(rld, intgroup="condition")
dev.off()

dds <- DESeq(dds)

res <- results(dds, contrast=c("condition","liquid","solid"))
resOrdered <- res[order(res$pvalue),]

## writing output files ordered by pvalue
write.csv(as.data.frame(resOrdered), 
          file="results/curvibacter/curvibacter_mono_liquid_vs_solid.csv")

res <- results(dds, contrast=c("condition","liquid","hydra"))
resOrdered <- res[order(res$pvalue),]

write.csv(as.data.frame(resOrdered), 
          file="results/curvibacter/curvibacter_mono_liquid_vs_hydra.csv")

res <- results(dds, contrast=c("condition","solid","hydra"))
resOrdered <- res[order(res$pvalue),]

write.csv(as.data.frame(resOrdered), 
          file="results/curvibacter/curvibacter_mono_solid_vs_hydra.csv")

