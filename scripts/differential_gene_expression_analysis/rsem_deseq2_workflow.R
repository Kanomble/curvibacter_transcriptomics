library(DESeq2)
library(tximport)
library(DEGreport)
library(writexl)

output_dir = file.path("results/curvibacter/rsem/")
experiments <- read.csv("data/Curvibacter_rsem/samples.csv")
files <- file.path("data/Curvibacter_rsem/rsem_results/",
                   paste0(experiments$ID,".genes.results"))
                   
names(files) <- experiments$id
all(file.exists(files))
sampleTable <- data.frame(condition = factor(experiments$treatment), batch = factor(experiments$batch), sample = experiments$sample)

id_df <- read.table(files[1], header=TRUE)
tx2gene <- data.frame(TXNAME=id_df$gene_id,GENEID=id_df$gene_id)

txi.rsem <- tximport(files, type="rsem",tx2gene = tx2gene, txIn = FALSE, txOut = FALSE)
txi.rsem$counts <- txi.rsem$counts + 1
txi.rsem$abundance <- txi.rsem$abundance + 1
txi.rsem$length <- txi.rsem$length + 1
rownames(sampleTable) <- colnames(txi.rsem$counts)


# recommended tximport vignette cmd
threshold <- 5
dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)

dds$run <- experiments$replicate
ddsColl <- collapseReplicates(dds, dds$sample, dds$run)

keep <- rowSums(counts(ddsColl)) >= threshold
dds <- ddsColl[keep,]


rld <- rlog(dds, blind=TRUE)
cat(paste("[*] PRODUCING PCA",file.path(outp,paste0("curvibacter_rsem_pca.png")),"\n"))
png(file=file.path(outp,paste0("curvibacter_rsem_pca.png")), width=800, height=550)
pca <- plotPCA(rld, intgroup="condition")
dev.off()

dds <- DESeq(dds)

counts <- counts(dds, normalized = TRUE)
design <- as.data.frame(colData(dds))
degCheckFactors(counts[, 1:6])

degQC(counts, design[["condition"]], pvalue = res[["pvalue"]])
resCov <- degCovariates(log2(counts(dds)+0.5),
                        colData(dds))
degPlotWide(dds, rownames(dds)[1:5], group="condition")
degObj(counts, design, "degObj.rda")
library(shiny)
shiny::runGitHub("lpantano/shiny", subdir="expression")

for (exp in unique(experiments$treatment)){
  for(nexp in unique(experiments$treatment)){
    if(nexp != exp){
      res <- results(dds, contrast=c("condition", exp, nexp))
      resOrdered <- res[order(res$pvalue),]
      write.csv(as.data.frame(resOrdered), 
                file=file.path(output_dir,paste0(exp,"_vs_",nexp,".csv")))
      resOrdered$ID <- rownames(resOrdered)
      write_xlsx(as.data.frame(resOrdered), file.path(output_dir,paste0(exp,"_vs_",nexp,".xlsx")))
      
      png(file=file.path(output_dir,paste0(paste0(exp,"_vs_",nexp),
                                           "_volcano_plot.png")),width=800, height=550)
      print(degVolcano(resOrdered[,c("log2FoldChange", "padj")]))
      dev.off()
    }
  }
}
