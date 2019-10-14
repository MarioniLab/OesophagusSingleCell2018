# This script manibulates the data for online storage
# Read in SCE object
library(Matrix)
sce <- readRDS("Dropbox (Personal)/Oesophagus_single_cell/All_corrected_sce.rds")
Matrix::writeMM(obj = counts(sce), file = "Dropbox (Personal)/Oesophagus_single_cell/For_sharing/raw_counts.mtx")

# colData
cData <- colData(sce)

# Add other column
tissue_type_pub <- colData(sce)$Tissue
tissue_type_pub[tissue_type_pub == "GC"] <- "NG"
tissue_type_pub[tissue_type_pub == "D2"] <- "ND"
tissue_type_pub[tissue_type_pub == "GOJ"] <- "N_SCJ"
tissue_type_pub[tissue_type_pub == "SCJ"] <- "B_SCJ"
cData <- colData(sce)
cData$Tissue_in_publication <- tissue_type_pub
cData <- cData[,c("Barcode", "Patient", "Tissue", "Tissue_in_publication", "Tissue_cluster", "cell_type", "tissue_type", "include", "confidence")]
write.csv(x = cData, file = "Dropbox (Personal)/Oesophagus_single_cell/For_sharing/cell_metadata.csv", quote = FALSE)


# rowData
rData <- rowData(sce)
write.csv(x = rData, file = "Dropbox (Personal)/Oesophagus_single_cell/For_sharing/gene_metadata.csv", quote = FALSE)
