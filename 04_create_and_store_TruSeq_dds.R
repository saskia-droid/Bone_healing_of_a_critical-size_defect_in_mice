
library(xlsx)
library(DESeq2)
library(gtools)

# import UMI count matrix:
data <- read.csv("../TruSeq_counts_matrix.txt", sep="\t", header=T, row.names=1) 

# sort columns in order that character strings containing embedded numbers 
# are sorted numerically rather than sorted by character value
data <- data[ , mixedsort(names(data))]

data <- subset(data, select = -c(X))

# create DESeq object

### change colnames 

metaData <- read.xlsx2("../RNAseq batch upload_Strunz.xlsx", 1, header=TRUE)

sample_names <- metaData$Name
rm(metaData)

sample_names <- gsub(" ", "", sample_names, fixed = TRUE)
sample_names <- gsub("+", "", sample_names, fixed = TRUE)

sample_names <- paste0(sample_names, rep(c("_1", "_2", "_3"), times=32))

colnames(data) <- sample_names

## transform the dataframe into a matrix
count_matrix <- as.matrix(data);
rm(data)

## create coldata matrix

operation <- c(rep(c('OVX', 'sham'), each=18), 
               rep(c('OVX', 'sham'), each=30))

medication <- c(rep(c('vehicle', 'ALN', 'vehicle', 'ALN'), each=9),
                rep(c('vehicle', 'ALN', 'vehicle', 'ALN'), each=15))

timepoint <- c(rep('6w', each=36), 
               rep('12w', each=60))

coated_with <- c(rep(rep(c('lowBMP', 'lowBMPL51P', 'highBMP'), each=3), times=4), 
                 rep(rep(c('control', 'L51P', 'lowBMP', 'lowBMPL51P', 'highBMP'), each=3), times=4))

coldata <- cbind(sample_names, operation, medication, timepoint, coated_with)
# put samples as rownames:
rownames(coldata) <- coldata[,1]
coldata <- coldata[, -1]

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ operation + medication + timepoint + coated_with)

# run DESeq analysis:
dds <- DESeq(dds)

saveRDS(dds, file="TruSeq_dds.rds")

