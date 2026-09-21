# Load packages
library(mygene)
library(readxl)
library(DESeq2)
library(WGCNA)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(EnhancedVolcano)
library(dplyr)
library(ggplot2)

# Load data & do a preprocessig
# count data
expression = read.csv("data/gene_count.csv",sep="\t")
rownames(expression) = expression$gene_id

# load sequecing depth and aligment info
reads_mapped = read_xlsx("data/Mapping_total_reads_fernando.xlsx")
reads_mapped = as.data.frame(reads_mapped)
rownames(reads_mapped) = reads_mapped$sample

# preparing dataset
cts <- expression[,c(2:27)] # only counts variables from each sample
coldata <- reads_mapped[,c(1,2)] # sample name and total reads

coldata$condition[grep("HDM[0-9]_HIP",reads_mapped$sample)] = "HDM_HIP"
coldata$condition[grep("HDM[0-9]_PFC",reads_mapped$sample)] = "HDM_PFC"
coldata$condition[grep("CTR[0-9]_HIP",reads_mapped$sample)] = "CTR_HIP"
coldata$condition[grep("CTR[0-9]_PFC",reads_mapped$sample)] = "CTR_PFC"
coldata$condition = as.factor(coldata$condition)

coldata$type = factor("paired-end")

coldata$experiment = NA
coldata$tissue = NA

coldata$experiment[grep("CTR",coldata$condition)] = "CTR"
coldata$experiment[grep("HDM",coldata$condition)] = "HDM"
coldata$tissue[grep("HIP",coldata$condition)] = "HIP"
coldata$tissue[grep("PFC",coldata$condition)] = "PFC"

coldata$experiment = as.factor(coldata$experiment)
coldata$tissue = as.factor(coldata$tissue)

# DESeq2 - differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

# removing low count genes
smallestGroupSize = 3
keep = rowSums(counts(dds) >= 10) >= smallestGroupSize
dds = dds[keep,]
# how many genes had low count and were filtered-out 
cat("initial gene count:",nrow(cts)) # 56748
cat("\nafter filtering:",sum(keep)) # 24055

dds <- DESeq(dds) # default settings
res <- results(dds)

# comparing  HDM_HIP vs. CTR_HIP
res_hip <- results(dds, contrast=c("condition","HDM_HIP","CTR_HIP")) 
resOrdered_hip <- res_hip[order(res_hip$pvalue),]
# comparing  HDM_PFC vs. CTR_PFC
res_pfc <- results(dds, contrast=c("condition","HDM_PFC","CTR_PFC"))
resOrdered_pfc <- res_pfc[order(res_pfc$pvalue),]

table_deg = as.data.frame(res_hip[which(abs(res_hip$stat)>1.5),])
table_deg_sig = as.data.frame(res_hip[which(res_hip$padj<0.05),])

# get gene symbol
all_transcripts = as.data.frame(res_hip)
gene_symbol_all = queryMany(rownames(all_transcripts), fields = "symbol", species = "mouse", size = 1)
all_transcripts$gene_symbol = gene_symbol_all$symbol

increased = all_transcripts$stat > 1.5
decreased = all_transcripts$stat < -1.5
change = vector(mode = "character",length = length(all_transcripts$stats))

change[increased] = "increased"
change[decreased] = "decreased"
change[is.na(change)] = "-"

all_transcripts$effect = change

# order all transcripts table
all_transcritos <- bind_rows(
  all_transcritos |> slice(which(rownames(all_transcritos) %in% rownames(table_deg))),
  all_transcritos |> slice(-which(rownames(all_transcritos) %in% rownames(table_deg))) |> arrange(factor(effect, levels = c("increased", "decreased", "-")))
)

gene_symbol_deg_sig = queryMany(rownames(table_deg_sig), fields = "symbol", species = "mouse", size = 1)
rownames(table_deg_sig) = gene_symbol_deg_sig$symbol

# Volcano plot
toptable = as.data.frame(res_hip)

gene_symbol = mapIds(org.Mm.eg.db,
                    keys=rownames(res_hip),
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
gene_symbol[which(is.na(gene_symbol))] = names(which(is.na(gene_symbol)))
# duplicated
gene_symbol[which(duplicated(gene_symbol))] = paste(gene_symbol[which(duplicated(gene_symbol))],"_1",sep="")

rownames(toptable) = gene_symbol

keyvals <- ifelse(
    abs(toptable$stat) < 1.5, 'darkgrey',
    ifelse(
    toptable$log2FoldChange < 0, '#2c7fb8',
      ifelse(toptable$log2FoldChange > 0, '#de2d26',
        'black')))

names(keyvals)[keyvals == '#2c7fb8'] <- 'reduction'
names(keyvals)[keyvals == '#de2d26'] <- 'increased'
names(keyvals)[keyvals == 'darkgrey'] <- '|t-statistic| < 1.5'

FC_cutoff = NA
x_limit = c(-5,5)
p_value_threshold = mean(res_hip$pvalue[which(abs(round(res_hip$stat,2)) == 1.5)])

plot = EnhancedVolcano(toptable,
    title = 'Volcano plot',
    subtitle = "CTRL vs. HDM",
    lab = rownames(toptable),
    x = 'log2FoldChange',
    y = 'pvalue',
    colCustom = keyvals,
    pointSize = 2.0,
    colAlpha = 0.4,
    labSize = NA,
    pCutoff = p_value_threshold,
    FCcutoff = FC_cutoff,
    drawConnectors = TRUE,
    boxedLabels = TRUE,
    ylim = c(0,11),
    max.overlaps = Inf,
    xlim = x_limit)

options(repr.plot.width = 5.5, repr.plot.height = 6, repr.plot.res = 150)
#svg("volcano_plot.svg",width = 5.5, height = 6) # save 
plot
#dev.off()

# co-expression WGCNA
deg_hip = as.data.frame(res_hip)

count_table = t(counts(dds))
HIP_samples = coldata[coldata$tissue == "HIP","sample"]

WGCNA_data_hip = count_table[HIP_samples,rownames(deg_hip[abs(deg_hip$stat)>1.5,])]

WGCNA_data_hip = apply(WGCNA_data_hip,2,as.numeric)
rownames(WGCNA_data_hip) = HIP_samples

# recommended setup in WGCNA documentation
options(stringsAsFactors = FALSE);
enableWGCNAThreads(nThreads = 12)

# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(11,20, by=1))

# Call the network topology analysis function for each set in turn
powerTables = list(data = pickSoftThreshold(WGCNA_data_hip, powerVector=powers,verbose = F)[[2]]) # 
collectGarbage()

hip_network = blockwiseModules(WGCNA_data_hip, power = 9, networkType = "signed",
                        TOMType = "signed", minModuleSize = 372,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        verbose = 1)

# Convert labels to colors for plotting
mergedColors = labels2colors(hip_network$colors)
names(mergedColors) = names(hip_network$colors)
clusters_hip_network = mergedColors

adjacency = adjacency(WGCNA_data_hip, power = 9, type = "signed") #Calculating the adjacency matrix
#help(adjacency )
TOM = TOMsimilarity(adjacency) #Calculating the topological overlap matrix
dissTOM = 1-TOM ##Calculating the dissimilarity

# Define numbers of genes and samples
nGenes = ncol(WGCNA_data_hip)
nSamples = nrow(WGCNA_data_hip)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(WGCNA_data_hip, mergedColors)$eigengenes
MEs_hip = orderMEs(MEs0)

module_mship_hip = signedKME(datExpr = WGCNA_data_hip, datME = MEs0)

Connectivity_hip = softConnectivity(WGCNA_data_hip,power=9,type = "signed")
names(Connectivity_hip) = rownames(WGCNA_data_hip)

# Plot the dendrogram and the module colors underneath
options(repr.plot.width = 5.5, repr.plot.height = 2.5, repr.plot.res = 150)
#svg("WGCNA_modules_detected_smaller_modules.svg", width = 5.5, height = 2.5)
plotDendroAndColors(hip_network$dendrograms[[1]], 
                    mergedColors[hip_network$blockGenes[[1]]],
                    "Module colors",
                    main = "Co-expression modules",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

table(clusters_hip_network) # check modules size

df_clean <- as.data.frame(count_table[1:13,]) %>% 
  select(where(~n_distinct(.) > 1))

PCA_filtered_genes = prcomp(WGCNA_data_hip,scale. = TRUE, center = TRUE)

PCA_dataset = cbind(PCA_filtered_genes$x[,c(1,2)],Condition = c(rep("CTR",7),rep("HDM",6)))
PCA_dataset = as.data.frame(PCA_dataset) %>% mutate(across(c(PC1,PC2),as.numeric))

options(repr.plot.width = 5.5, repr.plot.height = 2.5, repr.plot.res = 150)
#svg("PCA_all_gene_selected.svg", width = 5.5, height = 2.5)
ggplot(PCA_dataset, aes(y=PC1,x=PC2, color = Condition)) + 
geom_point(size=2.5) + theme_bw() + xlab("PC2 (23.10%)") + ylab("PC1 (30.67%)") +
geom_vline(xintercept = 0, linetype = "dotted",linewidth = 0.5) +
scale_color_manual(values = c("#67bed9","#ff7e79"))
#dev.off()

# plot boxplot
temp = MEs0[,-which(colnames(MEs0) == "MEgrey")]; temp$id = rownames(temp) 
melted_MEs = reshape2::melt(temp,id.vars = "id")

melted_MEs$Condition = NA
melted_MEs$Condition[grep("CTR",melted_MEs$id) ] = "CTR"
melted_MEs$Condition[grep("HDM",melted_MEs$id) ] = "HDM"

options(repr.plot.width = 5.5, repr.plot.height = 3.5, repr.plot.res = 150)

# blue = #67bed9
# red = #ff7e79

#svg("Eigengene_modules.svg", width = 5.5, height = 3.5)
ggplot(melted_MEs, aes(y = value, shape = Condition, x = Condition, fill = Condition)) + 
geom_violin(trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 0.8) +
scale_fill_manual(values = c("#67bed9","#ff7e79")) +
#geom_boxplot(width = 0.2, color = "black", fill = "white", alpha = 0.8, outlier.shape = NA) +
       geom_jitter(width = 0.1, alpha = 0.8, size=2) + scale_shape_manual(values = c(16, 15))+
       facet_grid(~variable) + theme_bw() + ylab("Module eigengene (ME)") +
       geom_hline(yintercept = 0, linetype = "dotted",linewidth = 0.5)
#dev.off()

options(repr.plot.width = 3.8, repr.plot.height = 3.5, repr.plot.res = 150)

# blue = #67bed9
# red = #ff7e79

#svg("Eigengene_blue_module.svg", width = 3.8, height = 3.5)
ggplot(melted_MEs[melted_MEs$variable == "MEblue",], aes(y = value, shape = Condition, x = Condition, fill = Condition)) + 
geom_violin(trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 0.8) +
scale_fill_manual(values = c("#67bed9","#ff7e79")) +
#geom_boxplot(width = 0.2, color = "black", fill = "white", alpha = 0.8, outlier.shape = NA) +
       geom_jitter(width = 0.1, alpha = 0.8, size=2) + scale_shape_manual(values = c(16, 15))+
       theme_bw() + ylab("Module eigengene (ME)") +
       geom_hline(yintercept = 0, linetype = "dotted",linewidth = 0.5) + ggtitle("Blue module")
#dev.off()

#svg("Eigengene_turquoise_module.svg", width = 3.8, height = 3.5)
ggplot(melted_MEs[melted_MEs$variable == "MEturquoise",], aes(y = value, shape = Condition, x = Condition, fill = Condition)) + 
geom_violin(trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 0.8) +
scale_fill_manual(values = c("#67bed9","#ff7e79")) +
#geom_boxplot(width = 0.2, color = "black", fill = "white", alpha = 0.8, outlier.shape = NA) +
       geom_jitter(width = 0.1, alpha = 0.8, size=2) + scale_shape_manual(values = c(16, 15))+
       theme_bw() + ylab("Module eigengene (ME)") +
       geom_hline(yintercept = 0, linetype = "dotted",linewidth = 0.5) + ggtitle("Turquoise module")
#dev.off()

 
### SEQ DEPTH

sequencing_depth = coldata[1:13,c(1,2)]
colnames(sequencing_depth) = c("id","seq_depth")
melted_MEs = merge(melted_MEs,sequencing_depth)
melted_MEs$Condition = factor(melted_MEs$Condition, levels = c("CTR","HDM"))
melted_MEs$seq_depth = scale(melted_MEs$seq_depth)

CTR_blue = melted_MEs$value[melted_MEs$variable == "MEblue" & melted_MEs$Condition == "CTR"]
HDM_blue = melted_MEs$value[melted_MEs$variable == "MEblue" & melted_MEs$Condition == "HDM"]
t.test(CTR_blue,HDM_blue)

MEblue_model = lm("value ~ Condition + seq_depth",melted_MEs[melted_MEs$variable == "MEblue",])
summary(MEblue_model)

CTR_turq = melted_MEs$value[melted_MEs$variable == "MEturquoise" & melted_MEs$Condition == "CTR"]
HDM_turq = melted_MEs$value[melted_MEs$variable == "MEturquoise" & melted_MEs$Condition == "HDM"]
t.test(CTR_turq,HDM_turq)

MEturquoise_model = lm("value ~ Condition + seq_depth",melted_MEs[melted_MEs$variable == "MEturquoise",])
summary(MEturquoise_model)

convert_to_entrez_id <- function(ensembl_genes){
    require(AnnotationDbi)
    entrez = mapIds(org.Mm.eg.db,
                    keys=ensembl_genes,
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
    return(entrez)}

blue = names(mergedColors[mergedColors == "blue"])
blue = rownames(module_mship_hip[blue,])[module_mship_hip[blue,"kMEblue"]>0.8]

turquoise = names(mergedColors[mergedColors == "turquoise"])
turquoise = rownames(module_mship_hip[turquoise,])[module_mship_hip[turquoise,"kMEturquoise"]>0.8]

blue_GO <- enrichGO(blue, keyType="ENSEMBL", OrgDb=org.Mm.eg.db,
                ont           = "all",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 1)
blue_GO_s = summary(blue_GO)

turquoise_GO <- enrichGO(turquoise, keyType="ENSEMBL", OrgDb=org.Mm.eg.db,
                ont           = "all",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 1)
turquoise_GO_s = summary(turquoise_GO)

blue_kegg = enrichKEGG(na.omit(convert_to_entrez_id(blue)),
                            organism = "mmu",
                            keyType = "kegg",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2)
blue_kegg_s = summary(blue_kegg)

turquoise_kegg = enrichKEGG(na.omit(convert_to_entrez_id(turquoise)),
                            organism = "mmu",
                            keyType = "kegg",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2)
turquoise_kegg_s = summary(turquoise_kegg)

blue_kegg@result$Description <- gsub(" - Mus musculus (house mouse)","",blue_kegg@result$Description,fixed = TRUE)

indx_to_remove = which(blue_kegg@result$category == "Human Diseases")
indx_to_remove_name = rownames(blue_kegg@result[indx_to_remove,])
blue_kegg@result <- blue_kegg@result[-indx_to_remove, ]
blue_kegg@geneSets = blue_kegg@geneSets[!(names(blue_kegg@geneSets) %in% indx_to_remove_name)]

turquoise_kegg@result$Description <- gsub(" - Mus musculus (house mouse)","",turquoise_kegg@result$Description,fixed = TRUE)

indx_to_remove = which(turquoise_kegg@result$category == "Human Diseases")
indx_to_remove_name = rownames(turquoise_kegg@result[indx_to_remove,])
turquoise_kegg@result <- turquoise_kegg@result[-indx_to_remove, ]
turquoise_kegg@geneSets = turquoise_kegg@geneSets[!(names(turquoise_kegg@geneSets) %in% indx_to_remove_name)]

options(repr.plot.width = 7, repr.plot.height = 5, repr.plot.res = 150)

# jpeg("KEGG_enrichment.jpeg",units ="in",width =7,height = 5, res=300)
#svg("KEGG_enrichment_blue_module.svg", width = 7, height = 5)
dotplot(blue_kegg,showCategory=15,decreasing=TRUE) + ggtitle("Blue module enrichment KEGG")
ggsave("KEGG_enrichment_blue_module.svg", width = 7, height = 5) 
#dev.off()
#svg("KEGG_enrichment_turq_module.svg", width = 7, height = 5)
dotplot(turquoise_kegg,showCategory=15,decreasing=TRUE) + ggtitle("Turquoise module enrichment KEGG")
ggsave("KEGG_enrichment_turq_module.svg", width = 7, height = 5) 
#dev.off()

options(repr.plot.width = 10, repr.plot.height = 7, repr.plot.res = 150)
p_facet <- barplot(
  blue_GO,
  x = "Count",
  showCategory = 10,
  split = "ONTOLOGY"
) +
  aes(fill = ONTOLOGY) +
  scale_fill_brewer(palette = "Set2") +
  enrichplot::autofacet(by = "row", scales = "free")

# jpeg("GO_enrichment.jpeg",units ="in",width =10,height = 7, res=300)
p_facet +
  geom_text(aes(label = Count), nudge_x = 2) +
  scale_y_discrete() + ggtitle("Blue module enrichment GO")
# dev.off()

p_facet <- barplot(
  turquoise_GO,
  x = "Count",
  showCategory = 10,
  split = "ONTOLOGY"
) +
  aes(fill = ONTOLOGY) +
  scale_fill_brewer(palette = "Set2") +
  enrichplot::autofacet(by = "row", scales = "free")

# jpeg("GO_enrichment.jpeg",units ="in",width =10,height = 7, res=300)
p_facet +
  geom_text(aes(label = Count), nudge_x = 2) +
  scale_y_discrete() + ggtitle("Turquoise module enrichment GO")
# dev.off()

MAPK = mapIds(org.Mm.eg.db,
                    keys=unlist(strsplit(blue_kegg@result$geneID[2],"/",)),
                    column="ENSEMBL",
                    keytype="ENTREZID",
                    multiVals="first")

MAPK_symbol = mapIds(org.Mm.eg.db,
                    keys=unlist(strsplit(blue_kegg@result$geneID[2],"/",)),
                    column="SYMBOL",
                    keytype="ENTREZID",
                    multiVals="first")

CALC = mapIds(org.Mm.eg.db,
                    keys=unlist(strsplit(blue_kegg@result$geneID[12],"/",)),
                    column="ENSEMBL",
                    keytype="ENTREZID",
                    multiVals="first")
CALC_symbol = mapIds(org.Mm.eg.db,
                    keys=unlist(strsplit(blue_kegg@result$geneID[12],"/",)),
                    column="SYMBOL",
                    keytype="ENTREZID",
                    multiVals="first")

count_MAPK = count_table[1:13,MAPK]
colnames(count_MAPK) = MAPK_symbol

count_CALC = count_table[1:13,CALC]
colnames(count_CALC) = CALC_symbol

MAPK_pca = prcomp(count_MAPK,scale. = TRUE, center = TRUE)
CALC_pca = prcomp(count_CALC,scale. = TRUE, center = TRUE)

MAPK_df = data.frame(values = -MAPK_pca$x[,1], Condition = c(rep("CTR",7),rep("HDM",6)))
CALC_df = data.frame(values = CALC_pca$x[,1], Condition = c(rep("CTR",7),rep("HDM",6)))

MEs0$Condition = c(rep("CTR",7),rep("HDM",6))

blue_plot = ggplot(MEs0,aes(y= MEblue, x= Condition, fill = condition))+geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.05, alpha = 0.8)
mapk_plot = ggplot(MAPK_df,aes(y= values, x= Condition, fill = condition))+geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.05, alpha = 0.8)
calc_plot = ggplot(CALC_df,aes(y= values, x= Condition, fill = condition))+geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.05, alpha = 0.8)

melted_enriched_genes = rbind(MAPK_df,CALC_df)
melted_enriched_genes$gene = c(rep("MAPK enriched genes",13),rep("CALC enriched genes",13))

mapk_calc_plot = ggplot(melted_enriched_genes,
                        aes(y= values, shape = Condition, x= Condition, fill = Condition)) + 
geom_violin(trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 0.8) + ylab("Eigengene")+#,
scale_fill_manual(values = c("#67bed9","#ff7e79")) + scale_shape_manual(values = c(16, 15))+
geom_jitter(width = 0.1, alpha = 0.8, size =2) + facet_grid(~gene) + theme_bw() +
geom_hline(yintercept = 0, linetype = "dotted",linewidth = 0.5)

options(repr.plot.width = 5.5, repr.plot.height = 3.5, repr.plot.res = 150)
#svg("Eigengene_enriched_subset.svg", width = 5.5, height = 3.5)
mapk_calc_plot #+ ggtitle("Blue module gene subsets")
#dev.off()

sequencing_depth = coldata[1:13,c(1,2)]
colnames(sequencing_depth) = c("id","seq_depth")
melted_enriched_genes$id = rownames(melted_enriched_genes)
melted_enriched_genes$id = gsub("HIP1","HIP",melted_enriched_genes$id)
melted_enriched_genes = merge(melted_enriched_genes,sequencing_depth, all=FALSE)
melted_enriched_genes$seq_depth = scale(melted_enriched_genes$seq_depth)

MAPK_model = lm("values ~ Condition + seq_depth",
                  melted_enriched_genes[melted_enriched_genes$gene == "MAPK enriched genes",])
summary(MAPK_model)

CALC_model = lm("values ~ Condition + seq_depth",
                  melted_enriched_genes[melted_enriched_genes$gene == "CALC enriched genes",])
summary(CALC_model)




wgcna_dend = hip_network$dendrograms[[1]]
WGCNA_data_hip_scaled = scale(WGCNA_data_hip)

order_clusters <- function(cluster, membership, order){
    ordered_cluster = factor(cluster, levels=order, ordered=T)
    ordered_cluster = sort(ordered_cluster)
    final = list()
    for(c in order){
        members = names(ordered_cluster[ordered_cluster==c])
        membership_value = membership[members,paste("kME",c,sep="")]
        names(membership_value) = members
        membership_value = sort(membership_value,decreasing = TRUE)
        
        final[[c]] = names(membership_value) 
    }
    return(final)
}

order = c("blue","turquoise","grey")
ordered_clusters = order_clusters(mergedColors,module_mship_hip,order)

options(repr.plot.width = 5, repr.plot.height = 8, repr.plot.res = 150)
#svg("Heatmap_samples_zscores.svg", width = 5, height = 8)
Heatmap(t(WGCNA_data_hip_scaled), name = "z-score", 
        column_order = sort(colnames(t(WGCNA_data_hip_scaled))),
        column_split = c(rep("CONTROL", 7), rep("HDM", 6)),
        row_order = c(ordered_clusters[[1]],ordered_clusters[[2]],ordered_clusters[[3]]),
        #row_split = factor(gsub("[0-9]","",names(unlist(ordered_clusters))),levels=c("blue","turquoise","grey")),
        show_row_names = FALSE,
        show_column_names = FALSE,
        column_gap = unit(3, "mm"),
        #cluster_row_slices = FALSE,
        #cluster_rows = FALSE,
        border = TRUE)
#dev.off()

