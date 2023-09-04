
library(DESeq2) 
library(ggplot2) 

dds <- readRDS(file="TruSeq_dds.rds")

vst.dds <- vst(dds, blind = TRUE)

vst.dds$coated_with <- recode_factor(vst.dds$coated_with, 
                                     control = "0 µg BMP2 / 0 µg L51P (control)", 
                                     L51P = "2.5 µg L51P", 
                                     highBMP = "2.5 µg BMP2", 
                                     lowBMP = "0.25 µg BMP2", 
                                     lowBMPL51P = "0.25 µg BMP2 / 2.5 µg L51P")

data <- plotPCA(vst.dds, intgroup = c("coated_with", "timepoint"), ntop=500, 
                returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, shape=coated_with, color=timepoint)) + 
  geom_point() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_shape_manual(values=c(1, 4, 16, 0, 7)) + 
  labs(shape='coated with') 

ggsave("PCA_coatings_timepoint.png")
