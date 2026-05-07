# Libraries
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(Seqinfo)
library(GenomeInfoDb)
library(BiocGenerics)
library(zellkonverter)

# the 10x hdf5 file contains both data types.
inputdata.10x <- Read10X_h5("./data/EGAD00001010038/P13_1mult/outs/filtered_feature_bc_matrix.h5")
# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
# Create Seurat object
obj.multi <- CreateSeuratObject(counts = rna_counts)
# Get % of mitochondrial genes
obj.multi[["percent.mt"]] <- PercentageFeatureSet(obj.multi, pattern = "^MT-")

# add the ATAC-seq assay
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
# Get gene annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# Change style to UCSC
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
# File with ATAC per fragment information file
frag.file <- "./data/EGAD00001010038/P13_1mult/outs/atac_fragments.tsv.gz"
# Add in ATAC-seq data as ChromatinAssay object
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
# Add the ATAC assay to the multiome object
obj.multi[["ATAC"]] <- chrom_assay
# Filter ATAC data based on QC metrics
obj.multi <- subset(
  x = obj.multi,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)
