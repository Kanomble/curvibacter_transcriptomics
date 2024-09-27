# use clusterProfiler for gene enrichment analysis
library(clusterProfiler)
output_dir <- "results/gene_enrichment_analysis/"

# read WP to AEP translation table and deseq2 dataframe
translation_table <- read.csv("data/curvibacter_annotation_files/translation_table_corrected.csv", sep="\t")
deseq_df <- read.csv("results/rsem/liquid_mono_culture_orgint_vs_metatranscriptome.csv")

# read data from the BLASTKoala tool
ko_mapping_df <- read.csv("data/curvibacter_annotation_files/kegg_curvibacter/curvibacter_wps_to_kegg_2.txt", sep = "\t")
ko_universe <- ko_mapping_df$ko

# rename deseq2 identifier to correct AEP identifier (without the gene: prefix)
deseq_df$X <- sapply(strsplit(as.character(deseq_df$X), "gene:"), function(x) x[2])
colnames(deseq_df)[colnames(deseq_df) == "X"] <- "old_locus_tag"
# filter all significant genes
deseq_df <- deseq_df[deseq_df$padj <= 0.05, ]

# merge translation table and deseq2 dataframe to get a WP to AEP identifier mapping in the deseq2 dataframe
merged_df <- merge(deseq_df, translation_table, by = "old_locus_tag")
# get KO numbers for all WP identifier
merged_df <- merge(merged_df, ko_mapping_df, by="protein_id")

# symbiotic genes
symbiotic_genes <- read.csv("./results/gene_enrichment_analysis/symbiotic_genes/merged_translation_table.csv")
symbiotic_merged <- merge(symbiotic_genes, ko_mapping_df, by="protein_id")

transcriptome_symbiotic_merged_df <- merge(symbiotic_genes,merged_df,by="protein_id")
up_sym <- transcriptome_symbiotic_merged_df[transcriptome_symbiotic_merged_df$log2FoldChange <= -1, ]
down_sym <- transcriptome_symbiotic_merged_df[transcriptome_symbiotic_merged_df$log2FoldChange >= 1, ]

kegg_enrich_result <- enrichKEGG(gene = symbiotic_merged$ko, universe = ko_universe, organism = 'ko')
write.csv(kegg_enrich_result@result,"results/gene_enrichment_analysis/kegg_enrich_symbiotic.csv")

png(file=file.path(output_dir,paste0("curvibacter_kegg_symbiotic.png")), width=800, height=550)
dotplot <- dotplot(kegg_enrich_result, showCategory=20)
print(dotplot)
dev.off()

kegg_enrich_result <- enrichKEGG(gene = up_sym$ko, universe = ko_universe, organism = 'ko')
write.csv(kegg_enrich_result@result,"results/gene_enrichment_analysis/kegg_enrich_up_symbiotic.csv")

png(file=file.path(output_dir,paste0("curvibacter_kegg_up_symbiotic.png")), width=800, height=550)
dotplot <- dotplot(kegg_enrich_result, showCategory=20)
print(dotplot)
dev.off()

kegg_enrich_result <- enrichKEGG(gene = down_sym$ko, universe = ko_universe, organism = 'ko')
write.csv(kegg_enrich_result@result,"results/gene_enrichment_analysis/kegg_enrich_down_symbiotic.csv")

png(file=file.path(output_dir,paste0("curvibacter_kegg_down_symbiotic.png")), width=800, height=550)
dotplot <- dotplot(kegg_enrich_result, showCategory=20)
print(dotplot)
dev.off()


# genes upregulated on host
upregulated <- merged_df[merged_df$log2FoldChange <= -1, ]

# writing upregulated into file
combined <- cbind(upregulated[upregulated$ko != "", ]$protein_id, upregulated[upregulated$ko != "", ]$ko)
write.csv(combined, "results/gene_enrichment_analysis/upregulated_kegg_ids.txt")

kegg_enrich_result <- enrichKEGG(gene = upregulated$ko, universe = ko_universe, organism = 'ko')
write.csv(kegg_enrich_result@result,"results/gene_enrichment_analysis/kegg_enrich_upregulated.csv")

png(file=file.path(output_dir,paste0("curvibacter_kegg_upregulated_on_host.png")), width=800, height=550)
dotplot <- dotplot(kegg_enrich_result, showCategory=20)
print(dotplot)
dev.off()

go_enrich_result <- enrichGO(gene=upregulated$old_locus_tag, OrgDb="org.CAEP13.eg.db", ont="ALL", keyType = "GID")

png(file=file.path(output_dir,paste0("curvibacter_go_upregulated_on_host.png")), width=800, height=550)
dotplot <- dotplot(go_enrich_result, showCategory=20)
print(dotplot)
dev.off()

# genes downregulated on host
downregulated <- merged_df[merged_df$log2FoldChange >= 1, ]
kegg_enrich_result <- enrichKEGG(gene = downregulated$ko, universe = ko_universe, organism = 'ko')

write.csv(kegg_enrich_result@result,"results/gene_enrichment_analysis/kegg_enrich_downregulated.csv")


png(file=file.path(output_dir,paste0("curvibacter_kegg_downregulated_on_host.png")), width=800, height=550)
dotplot <- dotplot(kegg_enrich_result, showCategory=20)
print(dotplot)
dev.off()


upregulated <- merged_df[merged_df$log2FoldChange <= -1, ]
highly_upregulated <- merged_df[merged_df$log2FoldChange <= -3, ]


go_enrich_result_highly_up <- enrichGO(gene=highly_upregulated$old_locus_tag, OrgDb="org.CAEP13.eg.db", ont="ALL", keyType = "GID")
go_enrich_result_up <- enrichGO(gene=upregulated$old_locus_tag, OrgDb="org.CAEP13.eg.db", ont="ALL", keyType = "GID")

result_list <- list(highly_up=go_enrich_result_highly_up, normal_upregulated=go_enrich_result_up)
go_enrich_result <- merge_result(result_list)

png(file=file.path(output_dir,paste0("highly_upregulated_go.png")), width=800, height=550)
dotplot <- dotplot(go_enrich_result, showCategory=20)
print(dotplot)
dev.off()


go_enrich_result <- enrichGO(gene=downregulated$old_locus_tag, OrgDb="org.CAEP13.eg.db", ont="ALL", keyType = "GID")

png(file=file.path(output_dir,paste0("curvibacter_go_downregulated_on_host.png")), width=800, height=550)
dotplot <- dotplot(go_enrich_result, showCategory=20)
print(dotplot)
dev.off()