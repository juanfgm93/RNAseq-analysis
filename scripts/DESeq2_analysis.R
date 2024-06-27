if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("goseq")
BiocManager::install("AnnotationDbi")
BiocManager::install("pcaExplorer")
BiocManager::install("ComplexHeatmap")
BiocManager::install("GOSemSim")
BiocManager::install("clusterProfiler")
BiocManager::install("edgeR")
install.packages("venneuler")
install.packages("gtools")
install.packages("tibble")
install.packages("pheatmap")
install.packages("gplots")
install.packages("stringr")
install.packages("corpcor")
install.packages("DataCombine")
install.packages("dplyr")
install.packages("grid.Extra")
install.packages("checkmate")


library(limma)
library(edgeR)
library(goseq)
library(GO.db)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(gplots)
library(gtools)
#library(venneuler)
library(checkmate)
library('DESeq2')
library(tibble)
library("pheatmap")
library(ggplot2)
library('RColorBrewer')
library(stringr) 
library(genefilter)
library("GOSemSim")
library("corpcor")
library(clusterProfiler)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)
library(DataCombine)
library(gridExtra)

# Get the gene counts from htseq
setwd("C:/Users/juan_/OneDrive/Escritorio/Lab/R/RNAseq results/Raw counts")

rawcount1 <- read.table("PJC1-Luc.txt", header=FALSE )
colnames(rawcount1)<- c("id","PJC1_Luc")
rawcount2 <- read.table("PJC4-Luc.txt", header=FALSE )
colnames(rawcount2)<- c("id","PJC4_Luc")
rawcount3 <- read.table("PJC7-Luc.txt", header=FALSE )
colnames(rawcount3)<- c("id","PJC7_Luc")

rawcount4 <- read.table("PJC2-DL2.txt", header=FALSE )
colnames(rawcount4)<- c("id","PJC2_DL2")
rawcount5 <- read.table("PJC5-DL2.txt", header=FALSE )
colnames(rawcount5)<- c("id","PJC5_DL2")
rawcount6 <- read.table("PJC8-DL2.txt", header=FALSE )
colnames(rawcount6)<- c("id","PJC8_DL2")

rawcount7 <- read.table("PJC3-DL2-TUTs.txt", header=FALSE )
colnames(rawcount7)<- c("id","PJC3_DL2_TUTs")
rawcount8 <- read.table("PJC6-DL2-TUTs.txt", header=FALSE )
colnames(rawcount8)<- c("id","PJC6_DL2_TUTs")
rawcount9 <- read.table("PJC9-DL2-TUTs.txt", header=FALSE)
colnames(rawcount9)<- c("id","PJC9_DL2_TUTs")
# Merge the samples in one table per condition VS control

rawcountPCA_DL2<- merge(rawcount4, merge(rawcount5, merge(rawcount6, merge(rawcount1, merge(rawcount2,rawcount3,  by = "id", all = TRUE), by = "id", all = TRUE), by = "id", all = TRUE), by = "id", all = TRUE), by = "id", all = TRUE)

rawcountPCA_DL2_TUTs<-merge(rawcount7,
                            merge(rawcount8, merge(rawcount9, merge(rawcount1, merge(rawcount2,rawcount3,  by = "id", all = TRUE), by = "id", all = TRUE), by = "id", all = TRUE), by = "id", all = TRUE), by = "id", all = TRUE)  #DF is ok, ignore warnings and remove column names


#More than 10 counts genes only
rawcountPCA_DL2_10<-rawcountPCA_DL2[rawcountPCA_DL2[2:length(colnames(rawcountPCA_DL2))]>10,] 

rawcountPCA_DL2_TUTs_10<-rawcountPCA_DL2_TUTs[rawcountPCA_DL2_TUTs[2:length(colnames(rawcountPCA_DL2_TUTs))]>10,] 


rawcountPCA_DL2_10<-na.omit(rawcountPCA_DL2_10)
rawcountPCA_DL2_TUTs_10<-na.omit(rawcountPCA_DL2_TUTs_10)


# Remove alignment stats rows

rawcountPCA_DL2<-rawcountPCA_DL2[!grepl("__no_feature", rawcountPCA_DL2$id),]
rawcountPCA_DL2<-rawcountPCA_DL2[!grepl("__ambiguous", rawcountPCA_DL2$id),]
rawcountPCA_DL2<-rawcountPCA_DL2[!grepl("__too_low_aQual", rawcountPCA_DL2$id),]
rawcountPCA_DL2<-rawcountPCA_DL2[!grepl("__not_aligned", rawcountPCA_DL2$id),]
rawcountPCA_DL2<-rawcountPCA_DL2[!grepl("__alignment_not_unique", rawcountPCA_DL2$id),]

rawcountPCA_DL2_TUTs<-rawcountPCA_DL2_TUTs[!grepl("__no_feature", rawcountPCA_DL2_TUTs$id),]
rawcountPCA_DL2_TUTs<-rawcountPCA_DL2_TUTs[!grepl("__ambiguous", rawcountPCA_DL2_TUTs$id),]
rawcountPCA_DL2_TUTs<-rawcountPCA_DL2_TUTs[!grepl("__too_low_aQual", rawcountPCA_DL2_TUTs$id),]
rawcountPCA_DL2_TUTs<-rawcountPCA_DL2_TUTs[!grepl("__not_aligned", rawcountPCA_DL2_TUTs$id),]
rawcountPCA_DL2_TUTs<-rawcountPCA_DL2_TUTs[!grepl("__alignment_not_unique", rawcountPCA_DL2_TUTs$id),]

rawcountPCA_DL2_10<-rawcountPCA_DL2_10[!grepl("__no_feature", rawcountPCA_DL2_10$id),]
rawcountPCA_DL2_10<-rawcountPCA_DL2_10[!grepl("__ambiguous", rawcountPCA_DL2_10$id),]
rawcountPCA_DL2_10<-rawcountPCA_DL2_10[!grepl("__too_low_aQual", rawcountPCA_DL2_10$id),]
rawcountPCA_DL2_10<-rawcountPCA_DL2_10[!grepl("__not_aligned", rawcountPCA_DL2_10$id),]
rawcountPCA_DL2_10<-rawcountPCA_DL2_10[!grepl("__alignment_not_unique", rawcountPCA_DL2_10$id),]

rawcountPCA_DL2_TUTs_10<-rawcountPCA_DL2_TUTs_10[!grepl("__no_feature", rawcountPCA_DL2_TUTs_10$id),]
rawcountPCA_DL2_TUTs_10<-rawcountPCA_DL2_TUTs_10[!grepl("__ambiguous", rawcountPCA_DL2_TUTs_10$id),]
rawcountPCA_DL2_TUTs_10<-rawcountPCA_DL2_TUTs_10[!grepl("__too_low_aQual", rawcountPCA_DL2_TUTs_10$id),]
rawcountPCA_DL2_TUTs_10<-rawcountPCA_DL2_TUTs_10[!grepl("__not_aligned", rawcountPCA_DL2_TUTs_10$id),]
rawcountPCA_DL2_TUTs_10<-rawcountPCA_DL2_TUTs_10[!grepl("__alignment_not_unique", rawcountPCA_DL2_TUTs_10$id),]

#Merge all samples (to make the PCA and Heatmap with all samples)

#rawcountPCA_all<-merge(rawcountPCA_DL2,rawcountPCA_DL2_TUTs, by = "id", all = TRUE)   ->[As? lo hace Marcelo]
rawcountPCA_all<-merge(rawcountPCA_DL2_10,rawcountPCA_DL2_TUTs_10, by = "id", all = TRUE)

rawcountPCA_all[,13] <- NULL
rawcountPCA_all[,12] <- NULL
rawcountPCA_all[,11] <- NULL

colnames(rawcountPCA_all)<- c("id",           "PJC2-DL2",       "PJC5-DL2",       "PJC8-DL2",       "PJC1-Luc",     "PJC4-Luc",     "PJC7-Luc",     "PJC3-DL2-TUTs",  "PJC6-DL2-TUTs", "PJC9-DL2-TUTs")

rawcountPCA_all<-rawcountPCA_all[,c(1:4,8,9,10,5:7)]

#function by parts
x <- rawcountPCA_DL2_TUTs_10
cond <-"DL2_TUTs"
x<-rawcountPCA_DL2_10 
cond<-"DL2"

#x<-rawcountPCA_all 
# Tables of raw gene counts
write.table(x, paste(cond,".tsv", sep =""), quote=FALSE, sep='\t', row.names=FALSE)
write.table(x, paste(cond,".csv", sep =""), quote=FALSE, sep=',', row.names=FALSE)
table<-read.csv(paste(cond,".tsv", sep =""),sep="\t",row.names="id")
colnames(table)<-colnames(x)[2:length(colnames(x))]
rawcountPCAcol <- as.matrix(table)
rawcountPCAcol <- na.omit(rawcountPCAcol)
### To Deseq

str<-cond
str(colnames(rawcountPCAcol))
ExpDesign <- data.frame(row.names=colnames(rawcountPCAcol), condition = c(rep(cond, sum(str_count(colnames(rawcountPCAcol),fixed(str,ignore_case=TRUE )))),rep("ctrl", sum(str_count(colnames(rawcountPCAcol), pattern = "Luc")))))
#ExpDesign <- data.frame(condition=c("DL2","DL2","DL2","ctrl","ctrl","ctrl","DL2-TUTs","DL2-TUTs","DL2-TUTs"))

DesqDset =  DESeqDataSetFromMatrix(countData = rawcountPCAcol, colData=ExpDesign, design=~condition)

### Negative Binomial Wald Test
#tests<-estimateSizeFactors(DesqDset)
#tests<-estimateDispersions(tests)
#tests <- nbinomWaldTest(tests)
#res_Wald <- results(tests)



#DESeq? this function performs...:
#Estimation of size factors, estimation of dispersion and negative binomial GLM fitting and Wald statistics
dds<-DESeq(DesqDset)

#rlog? this function transforms the count data to the log2 scale in a way which minimizes differences
#between samples for rows with small counts, and which normalizes with respect to library size
rld<- rlog(dds, blind=FALSE) #blind=TRUE if NT

#vsd? This is a wrapper for the varianceStabilizingTransformation (VST) that provides much faster estimation of the dispersion
#trend used to determine the formula for the VST
vsd <- vst(dds, blind=FALSE) #blind=TRUE if NT

rld_counts<-assay(rld) #rld_counts is only used for the heatmap all genes

#We make a table with normalized counts
rld_counts_2 <- data.frame(IDs=rownames(rld_counts),rld_counts,stringsAsFactors = F)
write.csv(rld_counts_2, paste(cond,"_rlog_counts.csv",sep=""),quote=FALSE,row.names = F)

vsd_counts<-assay(vsd) #vsd_counts is not used in the whole function

str_2<-str_extract(colnames(rld_counts),regex("DL2|DL2-TUTs|Luc", ignore_case=TRUE))
str_2<-str_extract(colnames(rld_counts),regex("DL2-TUTs|DL2|Luc", ignore_case=TRUE))

#Plots - PCA PLOT################################################################################
setwd("C:/Users/juan_/OneDrive/Escritorio/PhD tough stuff/RNA-seq/paper dis3l2")
plt <- plotPCA(rld, returnData=TRUE)
percentVar <- round(100 * attr(plt, "percentVar"))
sample=colnames(table)
ggplt_PCA<-ggplot(plt, aes(PC1, PC2, color=sample, shape=condition)) +
  geom_point(size=3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=0, size = 14, face = "bold",hjust = 0.5,vjust = 1),axis.title.x=element_text(size = 14,vjust = 0.2), 
        axis.title.y = element_text(size = 14), text = element_text(face = "bold"), plot.title = element_text(hjust = 0.5, vjust = 1, size = 16),
        legend.text = element_text(size = 14),legend.title = element_text(size=14),axis.text.y = element_text(size = 14, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  labs(title = "Principal component analysis plot") +
  coord_fixed()
ggplt_PCA
save(ggplt_PCA, file = paste(cond,"_PCA_rld.rds", sep =""))

ggsave(paste(cond,"_PCA_rld.jpeg", sep =""),plot = ggplt_PCA, width=13, height=6)
#################################################################################################

write.table(assay(normTransform(dds)), paste(cond,"_normalized_counts.tsv", sep =""), quote=FALSE, sep='\t', row.names=FALSE)
#this table is not used in the whole function

#Results - de DESeq2 - extracts a result table from a DESeq analysis giving base means across samples, log2 fold changes,
#standard errors, test statistics, p-values and adjusted p-values
res<-results(dds)

#Plots - MAplot##################################################################################
setwd("C:/Users/juan_/OneDrive/Escritorio/PhD tough stuff/RNA-seq/paper dis3l2")
DESeq2::plotMA(res,ylim=c(-6,6), colNonSig = "gray55",
               colSig = "red3")#,main='DESeq2')
dev.copy(png, file = paste(cond,"_MA_plot.png", sep=""))
dev.off()
#################################################################################################

save(res, file = paste(cond,"_res_DESeq.rds", sep = ""))

# We generate files with base mean, log2 fold changes, p-values...
write.table(res, paste(cond,"_all_genes_pvalues.tsv", sep =""), quote=FALSE, sep='\t', col.names = NA)
write.table(res, paste(cond,"_all_genes_pvalues.csv", sep =""), quote=FALSE, sep=',', col.names = NA)

#### Heatmap of the 500 highest variance genes 

select <- order(rowVars(assay(rld)),decreasing=TRUE)[1:500]
log2.norm.counts <- assay(rld)[select,]
df <- as.data.frame(colData(dds)[,c("condition", "sizeFactor")])


plot_pheatmap<-pheatmap(log2.norm.counts, scale="row",
                        cluster_rows=T, show_rownames=F,
                        cluster_cols=T, annotation_col=df)
ggsave(paste(cond,"_Heatmap_500_var_rld.jpeg", sep =""), plot = plot_pheatmap, width=13, height=9)

save(plot_pheatmap, file = paste(cond,"_Heatmap_500_var_rld.rds", sep =""))


#### Heatmap of all genes

mat_scaled = t(apply(rld_counts, 1, scale)) # scaling applied by row: negative values correspond to ones that are lower than the row mean
mat_scaled <- na.omit(mat_scaled)
colnames(mat_scaled) = colnames(rld_counts)
type = str_2
ha = HeatmapAnnotation(df = data.frame(type = type))

ht1=Heatmap(mat_scaled, name = "expression", km = 2, col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
            top_annotation = ha, top_annotation_height = unit(4, "mm"), 
            show_row_names = F, show_column_names = T)

save(ht1, file = paste(cond,"_Heatmap_10_counts_rld.rds", sep =""))
draw(ht1)
dev.copy(png, file = paste(cond,"_Heatmap_10_counts_rld.png", sep=""))
dev.off()

#### Sample to sample heatmap


distsRL <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(distsRL)
rownames(sampleDistMatrix) <- rownames(colData(dds)[,c("condition", "sizeFactor")])
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)


sample_heatmap<-pheatmap(sampleDistMatrix,
                         clustering_distance_rows=distsRL,
                         clustering_distance_cols=distsRL,
                         col=colors)

ggsave(paste(cond,"_Sample_to_sample_heatmap_rld.jpeg", sep =""), plot = sample_heatmap, width=13, height=9 )
save(sample_heatmap, file = paste(cond,"_Sample_to_sample_heatmap_rld.rds", sep =""))


#### Mapping

p <- res
head(p)
p<-na.omit(p)
p<-p[p$padj<=0.05,]
prownames(p)<-sub("\\.[0-9]+", "\\1", rownames(p))
DE_genes_6 <- p  
DE_genes_3 <- p[p$log2FoldChange>=1|p$log2FoldChange<=-1,]  # Upregulated and downregulated genes only
DE_genes_4 <- p[p$log2FoldChange>=1,]  # Upregulated (enriched) genes only
DE_genes_2 <- p[p$log2FoldChange>=1.5,]  # Upregulated (enriched) genes only
DE_genes_1 <- p[p$log2FoldChange>=2,]  # Upregulated (enriched) genes only

head(DE_genes_2)

DE_genes_2 <- as.data.frame(DE_genes_2)
DE_genes_2$symbol <- mapIds(org.Hs.eg.db, 
                            keys=row.names(DE_genes_2), 
                            column="SYMBOL", 
                            keytype="ENSEMBL",
                            multiVals="first")

DE_genes_2$entrez <- mapIds(org.Hs.eg.db, 
                            keys=row.names(DE_genes_2), 
                            column="ENTREZID", 
                            keytype="ENSEMBL",
                            multiVals="first")

DE_genes_2$name =   mapIds(org.Hs.eg.db,
                           keys=row.names(DE_genes_2), 
                           column="GENENAME",
                           keytype="ENSEMBL",
                           multiVals="first")

DE_genes_2$UNIPROT =   mapIds(org.Hs.eg.db,
                              keys=row.names(DE_genes_2), 
                              column="UNIPROT",
                              keytype="ENSEMBL",
                              multiVals="first")
#From duplicated gene names, select the ones with highest overall expression values

write.csv(DE_genes_2, file = paste(cond,"_All_gene_names_DE.csv", sep =""))

save(DE_genes_2, file = paste(cond,"_All_gene_names_DE.rds", sep =""))

#### Kegg terms

#First of all, we redefine the universe
universe_DE=read.csv("DL2_all_genes_pvalues.csv", header=T, stringsAsFactors = F)
for (i in 1:length(universe_DE[,1])){
  universe_DE[i,1]=unlist(strsplit(universe_DE[i,1], split=".", fixed=T))[1]
}
#We can also redefine the universe, taking out the NA values, with... na.omit

#universeENTREZ_from_ENSMBL=mapIds(org.Hs.eg.db, keys=(universe_DE[,1]),column="ENTREZID", keytype = "ENSEMBL", multivals="first")
#universeSYMBOL_from_ENSMBL=mapIds(org.Hs.eg.db, keys=(universe_DE[,1]),column="SYMBOL", keytype = "ENSEMBL", multivals="first")

universeENTREZ_from_ENSMBL=mapIds(org.Hs.eg.db, keys=(rownames(universe_DE)),column="ENTREZID", keytype = "ENSEMBL", multivals="first")
universeSYMBOL_from_ENSMBL=mapIds(org.Hs.eg.db, keys=(rownames(universe_DE)),column="SYMBOL", keytype = "ENSEMBL", multivals="first")
universeUNIPROT_from_ENSMBL=mapIds(org.Hs.eg.db, keys=(rownames(universe_DE)),column="UNIPROT", keytype = "ENSEMBL", multivals="first")


#Aentrezgenes_enrichKeggres<-enrichKEGG(DE_genes_2[DE_genes_2$padj <= 0.05,]$entrez,organism = "hsa", keyType = "ncbi-geneid")
#Aentrezgenes_enrichKeggres<-as.data.frame(enrichKEGG(DE_genes_2$entrez,organism = "hsa", keyType = "ncbi-geneid", universe = universeENTREZ_from_ENSMBL ,pvalueCutoff = 0.2))
Aentrezgenes_enrichKeggres<-as.data.frame(enrichKEGG(DE_genes_2$entrez,organism = "hsa", keyType = "ncbi-geneid",universe = universeENTREZ_from_ENSMBL))
#So we have decided to play with the whole universe
#Aentrezgenes_enrichKeggres<-as.data.frame(Aentrezgenes_enrichKeggres)
write.csv(Aentrezgenes_enrichKeggres,file=paste(cond,"_enrichKegg_results_significant.csv", sep=""))
save(Aentrezgenes_enrichKeggres, file = paste(cond,"_enrichKegg_results_significant.rds", sep =""))

#RAW_FUNC_ENRICH_DATAFRAME_A<-Aentrezgenes_enrichKeggres

#write.csv(RAW_FUNC_ENRICH_DATAFRAME_A,file=paste(cond,"_enrichKegg_results_significant.csv", sep=""))

#save(RAW_FUNC_ENRICH_DATAFRAME_A, file = paste(cond,"_enrichKegg_results_significant.rds", sep =""))



# Marina's algorytm to group the GO terms

#Aentrezgenes_enrichGOres=enrichGO(DE_genes_2[DE_genes_2$padj <= 0.05,]$entrez,'org.Hs.eg.db',ont="BP")                     ## Biological Process all transcriptome universe
#if(length(as.data.frame(Aentrezgenes_enrichGOres))>0){Aentrezgenes_enrichGOres=clusterProfiler::simplify(Aentrezgenes_enrichGOres)}

#Aentrezgenes_enrichGOres<-as.data.frame(enrichGO(DE_genes_2$entrez,'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", universe = universeENTREZ_from_ENSMBL))                     ## Biological Process all transcriptome universe

Aentrezgenes_enrichGOres<-as.data.frame(enrichGO(DE_genes_2[DE_genes_2$padj <= 0.05,]$entrez,'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", universe = universeENTREZ_from_ENSMBL))                     ## Biological Process all transcriptome universe
#Aentrezgenes_enrichGOres$geneID<-mapIds(org.Hs.eg.db,keys = (Aentrezgenes_enrichGOres[,8]),column = "SYMBOL", keytype = "GENEID", multiVals = "first")
write.csv(Aentrezgenes_enrichGOres,file=paste(cond,"_BP_enrichGO_results_significant.csv",sep=""))
save(Aentrezgenes_enrichGOres, file = paste(cond,"_BP_enrichGO_results_significant.rds", sep =""))

#Aentrezgenes_enrichGOres<-as.data.frame(Aentrezgenes_enrichGOres)

#Aentrezgenes_enrichGOres_MF=enrichGO(DE_genes_2[DE_genes_2$padj <= 0.05,]$entrez,'org.Hs.eg.db',ont="MF")                  ## Molecular Function all transcriptome universe
#if(length(as.data.frame(Aentrezgenes_enrichGOres_MF))>0){Aentrezgenes_enrichGOres_MF=clusterProfiler::simplify(Aentrezgenes_enrichGOres_MF)}

#Aentrezgenes_enrichGOres_MF2=as.data.frame(enrichGO(DE_genes_2$entrez,'org.Hs.eg.db', keyType = "ENTREZID", ont="MF", universe = universeENTREZ_from_ENSMBL))

Aentrezgenes_enrichGOres_MF=as.data.frame(enrichGO(DE_genes_2[DE_genes_2$padj <= 0.05,]$entrez,'org.Hs.eg.db', keyType = "ENTREZID", ont="MF", universe = universeENTREZ_from_ENSMBL))
write.csv(Aentrezgenes_enrichGOres_MF,file=paste(cond,"_MF_enrichGO_results_significant.csv",sep=""))
save(Aentrezgenes_enrichGOres_MF, file = paste(cond,"_MF_enrichGO_results_significant.rds", sep =""))

#Aentrezgenes_enrichGOres_MF<-as.data.frame(Aentrezgenes_enrichGOres_MF)

GenesinGO<-{}
Gene.Name<-{}
Gene.Name_df<-{}
j<-{}
forlist=list(Aentrezgenes_enrichGOres, Aentrezgenes_enrichGOres_MF)
name_file<-c("BP", "MF")

for (count in 1:length(forlist)) {
  j= as.data.frame(forlist[count])
  ## Convert entrez to gene name
  if (length(j)>0) {
    GenesinGO<-unique(unlist(strsplit(j$geneID,"[////]|[^[:print:]]",fixed=FALSE))) # Get all the entrez
    
    
    Gene.Name <- mapIds(org.Hs.eg.db, 
                        keys=GenesinGO, 
                        column="GENENAME", 
                        keytype="ENTREZID",
                        multiVals="first")
    
    Gene.Name_df<-data.frame(entrez=names(Gene.Name), name=unname(Gene.Name)) # Make a DF with the entrez and Gene names
    
    j<-FindReplace(j, "geneID", Gene.Name_df, from = "entrez", to = "name",
                   exact = F, vector = FALSE) # Replace entrez by gene names in the GO terms DF
    save(j, file =paste(cond,"_",name_file[count],"_enrichGO_results_significant.rds", sep="") )
    write.csv(j,file=paste(cond,"_",name_file[count],"_enrichGO_results_significant.csv", sep=""))
    j<-{}
  }else{j<-{}}
  
}

Aentrezgenes_enrichGOres=as.data.frame(Aentrezgenes_enrichGOres)

destfile=paste(cond,"_BP_enrichGO_results_significant.rds", sep="")

if(file.exists(destfile)){
  load(destfile)
  RAW_FUNC_ENRICH_DATAFRAME_A=j
} else {RAW_FUNC_ENRICH_DATAFRAME_A<-{}}


RAW_FUNC_ENRICH_DATAFRAME_A=Aentrezgenes_enrichGOres


##################################################################
##################################################################DESeq2 all samples
##################################################################
setwd("C:/Users/juan_/OneDrive/Escritorio/Lab/R/RNAseq results/Raw counts")

x<-rawcountPCA_all
cond<-"ALL"

# Raw counts table

write.table(x, paste(cond,".tsv", sep =""), quote=FALSE, sep='\t', row.names=FALSE)

table<-read.csv(paste(cond,".tsv", sep =""),sep="\t",row.names="id")

colnames(table)<-colnames(x)[2:length(colnames(x))]
table <- na.omit(table)
rawcountPCAcol <- as.matrix(table)

### To Deseq

str<-cond

str(colnames(rawcountPCAcol))


ExpDesign <- data.frame(row.names=colnames(rawcountPCAcol), condition = c("DL2",       "DL2",       "DL2",       "DL2-TUTs",  "DL2-TUTs",  "DL2-TUTs", "ctrl",       "ctrl",       "ctrl"))


DesqDset =  DESeqDataSetFromMatrix(countData = rawcountPCAcol, colData=ExpDesign, design=~condition)

### Negative Binomial Wald Test
#tests<-estimateSizeFactors(DesqDset)
#tests<-estimateDispersions(tests)
#tests <- nbinomWaldTest(tests)
#res_Wald <- results(tests)



dds<-DESeq(DesqDset)
rld<- rlog(dds, blind=FALSE) #blind=TRUE if NT
vsd <- vst(dds, blind=FALSE) #blind=TRUE if NT
#dds <- estimateSizeFactors(DesqDset)
#dds <- estimateDispersionsGeneEst(dds)
#dispersions(dds) <- mcols(dds)$dispGeneEst
#dds<-nbinomWaldTest(dds)

rld_counts<-assay(rld)
vsd_counts<-assay(vsd)

str_2<-c("DL2",       "DL2",       "DL2",       "DL2-TUTs",       "DL2-TUTs",  "DL2-TUTs",    "ctrl",       "ctrl",       "ctrl")

# PCA plot
setwd("C:/Users/juan_/Desktop/PhD tough stuff/RNA-seq")

plt <- plotPCA(rld, returnData=TRUE)
percentVar <- round(100 * attr(plt, "percentVar"))
sample=colnames(table)
ggplt_PCA<-ggplot(plt, aes(PC1, PC2, color=sample, shape=condition)) +
  geom_point(size=3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=0, size = 14, face = "bold",hjust = 0.5,vjust = 1),axis.title.x=element_text(size = 14,vjust = 0.2), 
        axis.title.y = element_text(size = 14), text = element_text(face = "bold"), plot.title = element_text(hjust = 0.5, vjust = 1, size = 16),
        legend.text = element_text(size = 14),legend.title = element_text(size=14),axis.text.y = element_text(size = 14, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  labs(title = "Principal component analysis plot") +
  coord_fixed()
print(ggplt_PCA)

save(ggplt_PCA, file = paste(cond,"_PCA_rld.rds", sep =""))
ggsave(paste(cond,"_PCA_rld.jpeg", sep =""),plot = ggplt_PCA, width=13, height=6)

### https://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/

setwd("C:/Users/juan_/OneDrive/Escritorio/Lab/R/RNAseq results/Raw counts")

write.table(assay(normTransform(dds)), paste(cond,"_normalized_counts.tsv", sep =""), quote=FALSE, sep='\t', row.names=FALSE)



res<-results(dds)

plotMA(res)
dev.copy(png, file = paste(cond,"_MA_plot.png", sep=""))
dev.off()
#ggsave(paste(cond,"_MA_plot.jpeg", sep = ""))
#plotSparsity(dds, normalized=TRUE)
#ggsave(paste(cond,"_Sparsity_plot.jpeg", sep = ""))

save(res, file = paste(cond,"_res_DESeq.rds", sep = ""))

# We generate files with base mean, log2 fold changes, p-values...
write.table(res, paste(cond,"_all_genes_pvalues.tsv", sep =""), quote=FALSE, sep='\t', col.names = NA)
write.table(res, paste(cond,"_all_genes_pvalues.csv", sep =""), quote=FALSE, sep=',', col.names = NA)

#### Heatmap of all genes
setwd("C:/Users/juan_/Desktop/PhD tough stuff/RNA-seq")

mat_scaled = t(apply(rld_counts, 1, scale)) # scaling applied by row: negative values correspond to ones that are lower than the row mean
mat_scaled <- na.omit(mat_scaled)
colnames(mat_scaled) = colnames(rld_counts)
type = str_2
ha = HeatmapAnnotation(df = data.frame(type = type))

ht1=Heatmap(mat_scaled, name = "expression", km = 2, col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
            top_annotation = ha, 
            show_row_names = F, show_column_names = T)
draw(ht1)

dev.copy(png, file = paste(cond,"_Heatmap_10_counts_rld.png", sep=""))
dev.off()
#ggsave(paste(cond,"_Heatmap_10_counts_rld.jpeg", sep =""))
#save(ht1, file = paste(cond,"_Heatmap_10_counts_rld.rds", sep =""))


### Sample to sample heatmap


distsRL <- dist(t(assay(rld)))
#mat<- as.matrix(distsRL)
#rownames(mat) <- rownames(colData(dds)[,c("condition", "sizeFactor")])
#hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
#hc <- hclust(distsRL)

sampleDistMatrix <- as.matrix(distsRL)
rownames(sampleDistMatrix) <- rownames(colData(dds)[,c("condition", "sizeFactor")])
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)



#heatmap.2(mat, Rowv=as.dendrogram(hc),
#symm=TRUE, trace='none',
#col = rev(hmcol), margin=c(13, 13))

pheatmap(sampleDistMatrix,
         clustering_distance_rows=distsRL,
         clustering_distance_cols=distsRL,
         col=colors)

sample_heatmap<-pheatmap(sampleDistMatrix,
                         clustering_distance_rows=distsRL,
                         clustering_distance_cols=distsRL,
                         col=colors)

ggsave(paste(cond,"_Sample_to_sample_heatmap_rld.jpeg", sep =""), plot = sample_heatmap, width=13, height=9 )
save(sample_heatmap, file = paste(cond,"_Sample_to_sample_heatmap_rld.rds", sep =""))

