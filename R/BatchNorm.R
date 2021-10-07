library(Matrix)
library(Seurat)
library(dplyr)
library(ggcyto)
library(RColorBrewer)
library(lisi)

colors.use <- c(RColorBrewer::brewer.pal(8, 'Set1')[c(1:5,7:8)],
                RColorBrewer::brewer.pal(7, 'Dark2'),
                'blue4', 'orangered', 'orangered4',
                c(RColorBrewer::brewer.pal(12, 'Paired'),
                  RColorBrewer::brewer.pal(6, 'Accent'))[c(1:10, 12, 15:18)])

# Function to z-score each PC's percentage of the total variance captured
#' Test PCA
#'
#' Identify important PCs from total principal components analysis (within a Seurat object).
#' TestPCA functions by generating a z-score corresponding to each respective PC's proportional contribution to the total variance.
#' Can be used similarly to the Seurat function ElbowPlot, which plots each successive PC by its standard deviation.
#' @param genes.use The vector of variable features used to construct the PCs
#' @param mtx.use The expression matrix used to construct the PCs.
#' At a minimum, this matrix must included the variable features included in "genes.use"
#' @return Returns a table containing the z-score of the cumulative percent of total variance for each PC
#' @export
TestPCA <- function(genes.use = object@assays$RNA@var.features,
                    mtx.use = object@assays$RNA@scale.data){
  data.use <- mtx.use[genes.use, ]
  pca.results <- svd(x = t(data.use))
  sdev <- pca.results$d/sqrt(max(1, ncol(data.use) - 1))
  PCVariance <- rbind(SD = sdev,
                      Proportion = (sdev^2)/sum(sdev^2),
                      Cumulative = cumsum(sdev^2)/sum(sdev^2))
  m <- mean(PCVariance['Proportion', ])
  s <- sd(PCVariance['Proportion', ])
  PCVariance <- rbind(PCVariance,
                      ZScore = (PCVariance['Proportion', ] - m)/s)
  return(PCVariance)
}


# Function to select mitochondrial percent by 5 Std Dev above median
#' Mitochondrial Transcripts Filter
#'
#' Automatically identify an appropriate threshold for subsetting a Seurat object by mitochondrial percentage (fraction of total transcripts beginning with "MT-".
#' MitoFilter first identifies the median mitochondrial percentage value five standard deviations above the median for a given sample.
#' MitoFilter then subsets a Seurat object, retaining only cells which have a lesser mitochondrial percentage than the cutoff (median +5 SD).
#' The object is returned containing only cells below threshold, as well as only cells with greter than 100 unique RNA transcripts.
#' @param obj The vector of variable features used to construct the PCs
#' @param mtx.use The Seurat object to subset.
#' This object must already contain a metadata feature "percent.mito" containing the total mitochondrial percentage, by cell.
#' This metadata feature can be added by the "PercentageFeatureSet" function.
#' @return Returns a filtered Seurat object with only cells passing the mitochondrial percentage filter, as well as containing greater than 100 unique RNA transcripts.
#' @examples
#' my_seurat_object <- MitoFilter(obj = my_seurat_object)
#' my_seurat_object <- my_seurat_object %>%
#'                     MitoFilter()
#' @export
MitoFilter <- function(obj = NULL){
  max.mito <<- 5*sd(obj$percent.mito) + median(obj$percent.mito)
  obj <- obj %>% subset(percent.mito < max.mito) %>% subset(nFeature_RNA > 100)
  return(obj)
  rm(max.mito)
}

# Function to correct feature names
#' Split "CD"
#'
#' Split feature names to extract "CDx" from ADT labels.
#' split.CD converts names in the format "CDx_xxxxxx" to "CDx", while eliminating whitespace and correcting capitalization.
#' A vector of feature names in the corrected format is returned.
#' @param x A vector of feature names to be corrected
#'
#' @return Returns a character vector of feature names in the desired format.
#' @examples
#' split.CD("CD19_TotalSeqC")
#' corrected.names <- split.CD(c("CD19_TotalSeqC", "CD45_TotalSeqC"))
#' @export
split.CD <- function(x = NULL){
  toupper(trimws(unlist(lapply(strsplit(x, "_"), '[[', 1))))
}

# Function to convert LISI score to "iLISI" score
#' Get iLISI
#'
#' Compute LISI scores (for each cell) for a Seurat object with a specified number of samples.
#' The median LISI score across all cells is then corrected for the total number of samples, and converted into a format where 0 = maximum integration and 1 = maximum segregation.
#' A single score on a scale of 0-1 is returned.
#' @param object A Seurat object with a UMAP computed, and sample identifiers stored as "orig.ident" in the object metadata slot.
#' @param nSamples The number of unique samples within the Seurat object (the number of samples which have been integrated in the UMAP).
#' @return Returns a numeric score from 0-1 (best-worst).
#' @examples
#' GetiLISI(object = my_seurat_object, nSamples = 3)
#' @export
GetiLISI <- function(object = NULL, nSamples = NULL) {
  coords <- object@reductions$umap@cell.embeddings
  sample_cats <- data.frame("orig.ident" = object@meta.data$orig.ident)
  res <- lisi::compute_lisi(coords, sample_cats, c('orig.ident'))
  sample.lisi <- median(res$orig.ident)
  adjusted.lisi <- 1 - (sample.lisi - 1) / (nSamples - 1)
  return(adjusted.lisi)
}

# Function to retrieve a CMS score given a single Seurat object
#' Get CMS score
#'
#' Compute CMS scores (for cells of a single sample) from a multi-sample Seurat object.
#' Requires a Seurat object with cell classifications stored as "Cell_Type" in the metadata and with sample ID stored as "orig.ident" in the metadata.
#' Also requires a dataframe containing the "reference" cell classifications generated from the single-sample dataset.
#' A single score on a scale of 0-1 is returned.
#' @param object A multi-sample Seurat object with sample identifiers stored as "orig.ident" and cell classifications as "Cell_Type" in the object metadata slot.
#' @param sample.ID The name of the single sample to be CMS-scored (must correspond to an ID stored in "orig.ident").
#' @param reference.ID A dataframe consisting of a single column, "Cell_Type", containing the cell classifications generated during a single-sample workflow (for use as a reference set without multi-sample batch effects).;
#' @return Returns a numeric score from 0-1 (best-worst).
#' @examples
#' GetCMS(object = my_seurat_object, sample.ID = "PBMC_01", reference.ID = PBMC01_SingleSample_RefID)
#' @export
GetCMS <- function(object = NULL, sample.ID = NULL, reference.ID = NULL){
  # Set object idents to sample ID, and subset the selected sample
  Idents(object) <- object[["orig.ident"]]
  id <- subset(object, idents = sample.ID)
  # Convert reference ID cell barcodes to match object barcodes format
  rownames(reference.ID) <- paste0(id, "_", rownames(reference.ID))
  # Select only cells shared between the current object and single-sample ID reference
  cells.use <- intersect(rownames(reference.ID), WhichCells(object))
  # Store ID as a vector
  test_ID <- as.character(object@meta.data[cells.use, 'Cell_Type'])
  reference.ID <- as.character(reference.ID[cells.use, 'Cell_Type'])
  # Create matrix with cell type columns
  Mismatch.mtx <- cbind(ref_ID, test_ID)
  rownames(Mismatch.mtx) <- cells.use
  # Find which test (group-workflow) IDs match reference (single-workflow)
  match <- Mismatch.mtx[, 'ref_ID'] == Mismatch.mtx[, 'test_ID']
  # Divide the number of matches by the total number of cells
  CMS <- 1 - sum(match, na.rm = T)/length(match)
  return(CMS)
}

