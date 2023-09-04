
library(DESeq2)
library(biomaRt)

# function to add gene symbol (and description) to a DESeqResults
# return a dataframe
add_gene_symbol <- function(deseq_res){
  mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
  genes <- rownames(deseq_res)
  G_list <- getBM(filters = "ensembl_gene_id",
                  attributes = c("ensembl_gene_id", "mgi_symbol",
                                 "mgi_description"),
                  values = genes, mart = mart)
  
  genes <- as.data.frame(deseq_res)
  genes$ensembl_gene_id <- rownames(genes)
  
  Genes <- merge(genes, G_list, by="ensembl_gene_id")
  
  rownames(Genes) <- Genes$ensembl_gene_id
  Genes <- subset(Genes, select = -c(ensembl_gene_id))
  
  Genes <- Genes[,c("mgi_symbol", "mgi_description", "baseMean", 
                    "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  
  Genes <- Genes[order(Genes$log2FoldChange),]
  
  return(Genes)
}