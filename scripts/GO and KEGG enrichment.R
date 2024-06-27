##############################################
#--------------------------------------------# 
#---        GO and KEGG enrichment        ---# 
#--------------------------------------------# 
##############################################

library(clusterProfiler)
library(GOSemSim)
library(ggplot2)

SemSimGO2 <- function(GO.df,ont,SimSemThreshold =0.7,method ="Lin"){
  library(GOSemSim)
  repeat{
    ss.selection <- data.frame(NA,NA,NA,NA)
    colnames(ss.selection) = colnames(GO.df)
    output <- ss.selection
    
    #prepares data and analizes the SemSim 
    if(ont=="BP"){
      if(!exists("ssDataBP")){ssDataBP <<- godata('org.Hs.eg.db', ont="BP", computeIC=T)}
      ss <- as.data.frame(termSim(GO.df[,1],GO.df[,1],semData=ssDataBP,method=method))
    }
    if(ont=="MF"){
      if(!exists("ssDataMF")){ssDataMF <<- godata('org.Hs.eg.db', ont="MF", computeIC=T)}
      ss <- as.data.frame(termSim(GO.df[,1],GO.df[,1],semData=ssDataMF,method=method))
    }
    
    #extracts similar GOs from analysis output according to threshold & selects the one with < pvalue
    for (c in 1:ncol(ss)){
      GOs <- rownames(ss[ss[,c]>SimSemThreshold,])
      f.df <- filter.df(GOs,GO.df,1)
      f.df[1,4] <- paste(unique(unlist(strsplit(f.df$geneID,split="/"))),collapse = "/")
      ss.selection <- rbind(ss.selection,f.df[1,])
    }
    
    #need to remove 1st line of NAs
    ss.selection <- ss.selection[-1,]
    
    for (s in unique(ss.selection$ID)){
      f.ss.selection <- filter.df(s,ss.selection,1)
      f.ss.selection[1,4] <- paste(unique(unlist(strsplit(f.ss.selection$geneID,split="/"))),collapse = "/")
      output <- rbind(output,f.ss.selection[1,])
    }
    
    #need to remove 1st line of NAs
    output <- output[-1,]
    
    if(nrow(output)==nrow(GO.df)){break} #if there is no change stops loop
    if(nrow(output)!=nrow(GO.df)){GO.df <- output} #if there is change let it loop
  }
  for(i in 1:nrow(output)){
    output[i,5] <- length(unique(unlist(strsplit(output$geneID[i],split="/"))))
  } #final step counts proteins in each select GOs
  output
}
filter.df <- function(filter,dataframe,col){
  d <- data.frame(dataframe[is.element(dataframe[,col],filter),],stringsAsFactors = F)
  colnames(d)[col] <- colnames(dataframe)[col]
  d
}

#Load gene list DIS3L2
setwd("C:/Users/juan_/OneDrive/Escritorio/Lab/R/RNAseq results")
DL2<-read.csv("DL2_All_gene_names_pvalues.csv",header=T, stringsAsFactors = F)
row.names(DL2)<-DL2$X
DL2<-DL2[,-1]
DL2<-na.omit(DL2)
DL2_significant<-DL2[DL2$padj<=0.05,]
DL2_significant_up<-DL2_significant[DL2_significant$log2FoldChange>0,]
DL2_significant_down<-DL2_significant[DL2_significant$log2FoldChange<0,]

#Load gene list DIS3L2+TUTs
setwd("C:/Users/juan_/OneDrive/Escritorio/Lab/R/RNAseq results")
DL2_TUTs<-read.csv("DL2_TUTs_All_gene_names_pvalues.csv",header=T, stringsAsFactors = F)
row.names(DL2_TUTs)<-DL2_TUTs$X
DL2_TUTs<-DL2_TUTs[,-1]
DL2_TUTs<-na.omit(DL2_TUTs)
DL2_TUTs_significant<-DL2_TUTs[DL2_TUTs$padj<=0.05,]
DL2_TUTs_significant_up<-DL2_TUTs_significant[DL2_TUTs_significant$log2FoldChange>0,]
DL2_TUTs_significant_down<-DL2_TUTs_significant[DL2_TUTs_significant$log2FoldChange<0,]

#GO and KEGG Enrichment DIS3L2 KD - with upregulated targets
setwd("C:/Users/juan_/OneDrive/Escritorio/Lab/R/RNAseq results/Enrichments/GO and Kegg if ALL UP genes are considered")

GO_KEGG_DIS3L2_denovo<-as.data.frame(enrichKEGG(DL2_significant_up$entrez,organism = "hsa", keyType = "ncbi-geneid"))
row.names(GO_KEGG_DIS3L2_denovo) <- NULL
write.csv(GO_KEGG_DIS3L2_denovo, file = "KEGG_all_upregulated_DIS3L2.csv")
GO_KEGG_DIS3L2_denovo <- read.csv(file =  "DL2_enrichKegg_results_significant.csv", stringsAsFactors = F)

GO_BP_DIS3L2_denovo<-as.data.frame(enrichGO(DL2_significant_up$entrez,'org.Hs.eg.db', keyType = "ENTREZID", ont="BP"))
row.names(GO_BP_DIS3L2_denovo) <- NULL
write.csv(GO_BP_DIS3L2_denovo, file = "GO_BP_all_upregulated_DIS3L2.csv")
GO_BP_DIS3L2_denovo <- read.csv(file =  "DL2_BP_enrichGO_results_significant.csv", stringsAsFactors = F)

GO_BP_DIS3L2_plot <- GO_BP_DIS3L2_denovo[c(2,1,14,43,7,11,12,17,18,36,39,78,79,80,81,95,96,144),]
GO_BP_DIS3L2_plot$ratio <- GO_BP_DIS3L2_plot$Count/2702
row.names(GO_BP_DIS3L2_plot) <- NULL
GO_BP_DIS3L2_plot$Description <- as.factor(GO_BP_DIS3L2_plot$Description)

setwd("C:/Users/juan_/Desktop/PhD tough stuff/RNA-seq/GO terms")
BP_DL2<-ggplot(data=GO_BP_DIS3L2_plot, aes(x=Description, y=ratio, fill = p.adjust)) +
  geom_bar(stat="identity") + theme_minimal() + coord_flip() +
  scale_fill_continuous(low="steelblue1", high="dodgerblue4", trans= "reverse") +
  scale_x_discrete(limits=rev(GO_BP_DIS3L2_plot$Description)) + 
  labs(title="Biological Processes GO Enrichment - DIS3L2 KD",x="GO terms",y="Ratio",fill="p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=0, size = 20, face = "bold",hjust = 0,vjust = 1),axis.title.x=element_text(size = 20,vjust = 0.2), 
        axis.title.y = element_blank(), text = element_text(face = "bold"), plot.title = element_text(hjust=0.5, vjust = 1, size = 28),
        legend.text = element_text(size = 20),legend.title = element_text(size=24),axis.text.y = element_text(size = 20, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA))
BP_DL2
ggsave("BP_DL2.tiff", width = 42, height = 25, units = "cm")

GO_KEGG_DIS3L2_plot <- GO_KEGG_DIS3L2_denovo
GO_KEGG_DIS3L2_plot$ratio <- GO_KEGG_DIS3L2_plot$Count/1295
row.names(GO_KEGG_DIS3L2_plot) <- NULL 
GO_KEGG_DIS3L2_plot$Description <- as.factor(GO_KEGG_DIS3L2_plot$Description)

KEGG_DL2<-ggplot(data=GO_KEGG_DIS3L2_plot, aes(x=Description, y=ratio, fill = p.adjust)) +
  geom_bar(stat="identity") + theme_minimal() + coord_flip() +
  scale_fill_continuous(low="chartreuse", high="chartreuse4", trans= "reverse") +
  scale_x_discrete(limits=rev(GO_KEGG_DIS3L2_plot$Description)) + 
  labs(title="KEGG enrichment analysis - DIS3L2 KD",x="GO terms",y="Ratio",fill="p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=0, size = 20, face = "bold",hjust = 0,vjust = 1),axis.title.x=element_text(size = 20,vjust = 0.2), 
        axis.title.y = element_blank(), text = element_text(face = "bold"), plot.title = element_text(hjust=0.5, vjust = 1, size = 28),
        legend.text = element_text(size = 20),legend.title = element_text(size=24),axis.text.y = element_text(size = 20, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA))
KEGG_DL2
ggsave("KEGG_DL2.tiff", width = 42, height = 25, units = "cm")

#GO and KEGG Enrichment DIS3L2 KD - with downregulated targets
setwd("C:/Users/juan_/OneDrive/Escritorio/Lab/R/RNAseq results/Enrichments/GO and Kegg if ALL DOWN genes are considered")

GO_KEGG_DIS3L2_denovo<-as.data.frame(enrichKEGG(DL2_significant_down$entrez,organism = "hsa", keyType = "ncbi-geneid"))
row.names(GO_KEGG_DIS3L2_denovo) <- NULL
write.csv(GO_KEGG_DIS3L2_denovo, file = "KEGG_all_dowregulated_DIS3L2.csv")
GO_KEGG_DIS3L2_denovo <- read.csv(file =  "KEGG_all_dowregulated_DIS3L2.csv", stringsAsFactors = F)

GO_BP_DIS3L2_denovo<-as.data.frame(enrichGO(DL2_significant_down$entrez,'org.Hs.eg.db', keyType = "ENTREZID", ont="BP"))
row.names(GO_BP_DIS3L2_denovo) <- NULL
write.csv(GO_BP_DIS3L2_denovo, file = "BP_all_dowregulated_DIS3L2.csv")
GO_BP_DIS3L2_denovo <- read.csv(file = "BP_all_dowregulated_DIS3L2.csv", stringsAsFactors = F)


GO_BP_DIS3L2_plot_down <- GO_BP_DIS3L2_denovo[c(1,3,6,8,15,18,24,28,29,65,73,109,146,155,170,247,270,330),]
GO_BP_DIS3L2_plot_down$ratio <- GO_BP_DIS3L2_plot_down$Count/3137
row.names(GO_BP_DIS3L2_plot_down) <- NULL
GO_BP_DIS3L2_plot_down$Description <- as.factor(GO_BP_DIS3L2_plot_down$Description)

setwd("C:/Users/juan_/Desktop/PhD tough stuff/RNA-seq/GO terms")
BP_DL2_DOWN<-ggplot(data=GO_BP_DIS3L2_plot_down, aes(x=Description, y=ratio, fill = p.adjust)) +
  geom_bar(stat="identity") + theme_minimal() + coord_flip() +
  scale_fill_continuous(low="firebrick1", high="firebrick4", trans= "reverse") +
  scale_x_discrete(limits=rev(GO_BP_DIS3L2_plot_down$Description)) + 
  labs(title="Biological Processes GO Enrichment - DIS3L2 KD",x="GO terms",y="Ratio",fill="p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=0, size = 20, face = "bold",hjust = 0,vjust = 1),axis.title.x=element_text(size = 20,vjust = 0.2), 
        axis.title.y = element_blank(), text = element_text(face = "bold"), plot.title = element_text(hjust=0.5, vjust = 1, size = 28),
        legend.text = element_text(size = 20),legend.title = element_text(size=24),axis.text.y = element_text(size = 20, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA))
BP_DL2_DOWN
ggsave("BP_DL2_DOWN.tiff", width = 42, height = 25, units = "cm")


GO_KEGG_DIS3L2_plot_down <- GO_KEGG_DIS3L2_denovo[c(2,4,7,8,9,10,11,15,18,25,26,27,28,29,31,36,44,47),]
GO_KEGG_DIS3L2_plot_down$ratio <- GO_KEGG_DIS3L2_plot_down$Count/1469
row.names(GO_KEGG_DIS3L2_plot_down) <- NULL 
GO_KEGG_DIS3L2_plot_down$Description <- as.factor(GO_KEGG_DIS3L2_plot_down$Description)

KEGG_DL2_DOWN<-ggplot(data=GO_KEGG_DIS3L2_plot_down, aes(x=Description, y=ratio, fill = p.adjust)) +
  geom_bar(stat="identity") + theme_minimal() + coord_flip() +
  scale_fill_continuous(low="goldenrod1", high="goldenrod4", trans= "reverse") +
  scale_x_discrete(limits=rev(GO_KEGG_DIS3L2_plot_down$Description)) + 
  labs(title="KEGG enrichment analysis - DIS3L2 KD",x="GO terms",y="Ratio",fill="p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=0, size = 20, face = "bold",hjust = 0,vjust = 1),axis.title.x=element_text(size = 20,vjust = 0.2), 
        axis.title.y = element_blank(), text = element_text(face = "bold"), plot.title = element_text(hjust=0.5, vjust = 1, size = 28),
        legend.text = element_text(size = 20),legend.title = element_text(size=24),axis.text.y = element_text(size = 20, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA))
KEGG_DL2_DOWN
ggsave("KEGG_DL2_DOWN.tiff", width = 42, height = 25, units = "cm")

#GO and KEGG Enrichment DIS3L2 + TUTs KD - with upregulated targets
GO_BP_DL2_TUTs_denovo<-as.data.frame(enrichGO(DL2_TUTs_significant_up$entrez,'org.Hs.eg.db', keyType = "ENTREZID", ont="BP"))
GO_KEGG_DL2_TUTs_denovo<-as.data.frame(enrichKEGG(DL2_TUTs_significant_up$entrez,organism = "hsa", keyType = "ncbi-geneid"))
row.names(GO_BP_DL2_TUTs_denovo) <- NULL
row.names(GO_KEGG_DL2_TUTs_denovo) <- NULL

setwd("C:/Users/juan_/OneDrive/Escritorio/Lab/R/RNAseq results/Enrichments/GO and Kegg if ALL UP genes are considered")
GO_BP_DIS3L2_TUTs_denovo <- read.csv(file =  "DL2_TUTs_BP_enrichGO_results_significant.csv", stringsAsFactors = F)
GO_KEGG_DIS3L2_TUTs_denovo <- read.csv(file =  "DL2_TUTs_enrichKegg_results_significant.csv", stringsAsFactors = F)

GO_BP_DL2_TUTs_plot <- GO_BP_DIS3L2_TUTs_denovo[c(1,6,8,10,12,16,20,30,34,40,41,44,70,83,89,92,93,98),]
GO_BP_DL2_TUTs_plot$ratio <- GO_BP_DL2_TUTs_plot$Count/2423
row.names(GO_BP_DL2_TUTs_plot) <- NULL
GO_BP_DL2_TUTs_plot$Description <- as.factor(GO_BP_DL2_TUTs_plot$Description)

setwd("C:/Users/juan_/Desktop/PhD tough stuff/RNA-seq/GO terms")
BP_DL2_TUTs<-ggplot(data=GO_BP_DL2_TUTs_plot, aes(x=Description, y=ratio, fill = p.adjust)) +
  geom_bar(stat="identity") + theme_minimal() + coord_flip() +
  scale_fill_continuous(low="steelblue1", high="dodgerblue4", trans= "reverse") +
  scale_x_discrete(limits=rev(GO_BP_DL2_TUTs_plot$Description)) + 
  labs(title="Biological Processes GO Enrichment - DIS3L2+TUTs KD",x="GO terms",y="Ratio",fill="p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=0, size = 20, face = "bold",hjust = 0,vjust = 1),axis.title.x=element_text(size = 20,vjust = 0.2), 
        axis.title.y = element_blank(), text = element_text(face = "bold"), plot.title = element_text(hjust=0.5, vjust = 1, size = 28),
        legend.text = element_text(size = 20),legend.title = element_text(size=24),axis.text.y = element_text(size = 20, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA))
BP_DL2_TUTs
ggsave("BP_DL2_TUTs.tiff", width = 42, height = 25, units = "cm")

GO_KEGG_DL2_TUTs_plot <- GO_KEGG_DIS3L2_TUTs_denovo[-c(19:24),]
GO_KEGG_DL2_TUTs_plot$ratio <- GO_KEGG_DL2_TUTs_plot$Count/1073
row.names(GO_KEGG_DL2_TUTs_plot) <- NULL 
GO_KEGG_DL2_TUTs_plot$Description <- as.factor(GO_KEGG_DL2_TUTs_plot$Description)

KEGG_DL2_TUTs <-ggplot(data=GO_KEGG_DL2_TUTs_plot, aes(x=Description, y=ratio, fill = p.adjust)) +
  geom_bar(stat="identity") + theme_minimal() + coord_flip() +
  scale_fill_continuous(low="chartreuse", high="chartreuse4", trans= "reverse") +
  scale_x_discrete(limits=rev(GO_KEGG_DL2_TUTs_plot$Description)) + 
  labs(title="KEGG enrichment analysis - DIS3L2+TUTs KD",x="GO terms",y="Ratio",fill="p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=0, size = 20, face = "bold",hjust = 0,vjust = 1),axis.title.x=element_text(size = 20,vjust = 0.2), 
        axis.title.y = element_blank(), text = element_text(face = "bold"), plot.title = element_text(hjust=0.5, vjust = 1, size = 28),
        legend.text = element_text(size = 20),legend.title = element_text(size=24),axis.text.y = element_text(size = 20, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA))
KEGG_DL2_TUTs
ggsave("KEGG_DL2_TUTs.tiff", width = 42, height = 25, units = "cm")

#GO and KEGG Enrichment DIS3L2 + TUTs KD - with down-regulated targets
setwd("C:/Users/juan_/OneDrive/Escritorio/Lab/R/RNAseq results/Enrichments/GO and Kegg if ALL DOWN genes are considered")
GO_KEGG_DL2_TUTs_denovo<-as.data.frame(enrichKEGG(DL2_TUTs_significant_down$entrez,organism = "hsa", keyType = "ncbi-geneid"))
row.names(GO_KEGG_DL2_TUTs_denovo) <- NULL
write.csv(GO_KEGG_DL2_TUTs_denovo, file = "KEGG_all_dowregulated_DL2_TUTs.csv")

GO_BP_DL2_TUTs_denovo<-as.data.frame(enrichGO(DL2_TUTs_significant_down$entrez,'org.Hs.eg.db', keyType = "ENTREZID", ont="BP"))
row.names(GO_BP_DL2_TUTs_denovo) <- NULL
write.csv(GO_BP_DL2_TUTs_denovo, file = "BP_all_dowregulated_DL2_TUTs.csv")

GO_BP_DIS3L2_TUTs_denovo <- read.csv(file = "BP_all_dowregulated_DL2_TUTs.csv", stringsAsFactors = F)
GO_KEGG_DIS3L2_TUTs_denovo <- read.csv(file =  "KEGG_all_dowregulated_DL2_TUTs.csv", stringsAsFactors = F)

setwd("C:/Users/juan_/Desktop/PhD tough stuff/RNA-seq/GO terms")

GO_BP_DL2_TUTs_plot_down <- GO_BP_DIS3L2_TUTs_denovo[c(3,4,5,9,13,17,22,27,35,37,39,48,69,71,99,113,132,146),]
GO_BP_DL2_TUTs_plot_down$ratio <- GO_BP_DL2_TUTs_plot_down$Count/2627
row.names(GO_BP_DL2_TUTs_plot_down) <- NULL
GO_BP_DL2_TUTs_plot_down$Description <- as.factor(GO_BP_DL2_TUTs_plot_down$Description)

BP_DL2_TUTs_DOWN<-ggplot(data=GO_BP_DL2_TUTs_plot_down, aes(x=Description, y=ratio, fill = p.adjust)) +
  geom_bar(stat="identity") + theme_minimal() + coord_flip() +
  scale_fill_continuous(low="firebrick1", high="firebrick4", trans= "reverse") +
  scale_x_discrete(limits=rev(GO_BP_DL2_TUTs_plot_down$Description)) + 
  labs(title="Biological Processes GO Enrichment - DIS3L2+TUTs KD",x="GO terms",y="Ratio",fill="p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=0, size = 20, face = "bold",hjust = 0,vjust = 1),axis.title.x=element_text(size = 20,vjust = 0.2), 
        axis.title.y = element_blank(), text = element_text(face = "bold"), plot.title = element_text(hjust=0.5, vjust = 1, size = 28),
        legend.text = element_text(size = 20),legend.title = element_text(size=24),axis.text.y = element_text(size = 20, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA))
BP_DL2_TUTs_DOWN
ggsave("BP_DL2_TUTs_DOWN.tiff", width = 42, height = 25, units = "cm")


GO_KEGG_DL2_TUTs_plot_down <- GO_KEGG_DIS3L2_TUTs_denovo[-c(5,20:24),]
GO_KEGG_DL2_TUTs_plot_down$ratio <- GO_KEGG_DL2_TUTs_plot_down$Count/1260
row.names(GO_KEGG_DL2_TUTs_plot_down) <- NULL 
GO_KEGG_DL2_TUTs_plot_down$Description <- as.factor(GO_KEGG_DL2_TUTs_plot_down$Description)

KEGG_DL2_TUTs_DOWN<-ggplot(data=GO_KEGG_DL2_TUTs_plot_down, aes(x=Description, y=ratio, fill = p.adjust)) +
  geom_bar(stat="identity") + theme_minimal() + coord_flip() +
  scale_fill_continuous(low="goldenrod1", high="goldenrod4", trans= "reverse") +
  scale_x_discrete(limits=rev(GO_KEGG_DL2_TUTs_plot_down$Description)) + 
  labs(title="KEGG enrichment analysis - DIS3L2+TUTs KD",x="GO terms",y="Ratio",fill="p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=0, size = 20, face = "bold",hjust = 0,vjust = 1),axis.title.x=element_text(size = 20,vjust = 0.2), 
        axis.title.y = element_blank(), text = element_text(face = "bold"), plot.title = element_text(hjust=0.5, vjust = 1, size = 28),
        legend.text = element_text(size = 20),legend.title = element_text(size=24),axis.text.y = element_text(size = 20, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA))
KEGG_DL2_TUTs_DOWN
ggsave("KEGG_DL2_TUTs_DOWN.tiff", width = 42, height = 25, units = "cm")
