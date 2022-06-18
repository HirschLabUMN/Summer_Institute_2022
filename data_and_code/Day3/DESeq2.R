### Get all needed packages ###
# DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")

# gplots
if (!requireNamespace("gplots", quietly = TRUE))
    install.packages("gplots")

library(DESeq2)
library(gplots)

### Set working directory and read in matrix data ###
#setwd("~/Desktop/")
fpkm <- read.delim(file="Transcriptomics_Report_Expression_Matrix_2022.txt", head=T, sep="\t", row.names=1)

### Set up the data
# Create dataframe with tissues as conditions 
design = data.frame(row.names=colnames(fpkm), condition=c("18DAP_Pericarp", "12DAP_Endopsperm", "V1_4D_PE_Primary_root", "V1_4D_PE_Primary_root", "V1_4D_PE_Primary_root", "V7_Tip_of_transition_leaf", "V7_Tip_of_transition_leaf", "V7_Tip_of_transition_leaf", "V7_Bottom_of_transition_leaf", "V7_Bottom_of_transition_leaf", "V7_Bottom_of_transition_leaf", "V13_Immature_tassel", "V13_Immature_tassel", "V13_Immature_tassel", "V18_Meiotic_Tassel", "V18_Meiotic_Tassel", "V18_Immature_Cob", "V18_Immature_Cob", "V18_Immature_Cob", "R1_Silks", "R1_Silks", "R1_Silks", "R1_Anthers", "R1_Anthers", "R1_Anthers", "2DAP_Whole_seed", "2DAP_Whole_seed", "2DAP_Whole_seed", "18DAP_Pericarp", "18DAP_Pericarp", "22DAP_Whole_Seed", "22DAP_Whole_Seed", "22DAP_Whole_Seed", "22DAP_Endosperm", "22DAP_Endosperm", "22DAP_Endosperm", "V9_Eighth_Leaf", "V9_Eighth_Leaf", "V9_Eighth_Leaf", "V9_Eleventh_Leaf", "V9_Eleventh_Leaf", "V9_Eleventh_Leaf", "V9_Thirteenth_Leaf", "V9_Thirteenth_Leaf", "V9_Thirteenth_Leaf", "V9_Papery_Leaves", "V9_Papery_Leaves", "V9_Papery_Leaves", "V18_Meiotic_Tassel", "VT_Thirteenth_Leaf", "VT_Thirteenth_Leaf", "VT_Thirteenth_Leaf", "12DAP_Endopsperm", "12DAP_Endopsperm", "16DAP_Embryo", "16DAP_Embryo", "16DAP_Embryo"))
condition = design$condition

# Create count data set varialble within DESeq2
cds <- DESeqDataSetFromMatrix(countData = fpkm,
                              colData = design,
                              design = ~ condition)

# Normalize the data
cds <- estimateSizeFactors(cds)
sizeFactors(cds)

cds <- estimateDispersions(cds)

# Create object with normalized values
norm <- counts(cds, normalized = TRUE)

# Add tissue names to data frame and replace SRA numbers 
# create a vector called tissues - take tissues from condition factor in design df (same order as SRA numbers in norm), denote them as a character object using as.character
# use names function to geta set of names from tissues and assign row.names from design (SRA numbers that match each tissue/rep) - now each Tissue is assigned to an SRA number in a vector
# apply tissues as column names of norm 
tissues = as.character(design$condition) 
names(tissues) = row.names(design)
colnames(norm) = tissues

### Plot correlations of normalized expression levels ###
# Calculate correlations fromt he normalized matrix
cor <- cor(norm, use = "pairwise.complete.obs", method = "spearman")
cor_matrix <- as.matrix(cor,  center = FALSE)

# plot the heatmap
pdf("heatmap.pdf", pointsize=10, font="Helvetica")
heatmap.2(cor_matrix, 
          key = TRUE, 
          keysize = 0.8, 
          density.info = "non",
          trace = "none",
          dendrogram = "row",
          col = "redgreen",
          margin = c(10,10))
dev.off()

# plot the heirarchical tree
dist <- as.dist((1-cor)/2)
pdf(file="tree.pdf", pointsize=9)
plot(hclust(dist))
dev.off()

## Test for DE genes ##
test <- DESeq(cds)

## Write out the test results ##
res <- results(test)

## Pulling out a specific contrast ##
endo_seed <- results(test, contrast = c("condition", "22DAP_Endosperm", "22DAP_Whole_Seed"))
nrow(subset(endo_seed, padj<0.01))

write.table(endo_seed, file="EndospermVsSeed_DEseq.txt", sep="\t")


