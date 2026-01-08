#' @title Dataset Registry
#' @description Registry of all available datasets organized by modality.
#' @name registry
#' @keywords internal
NULL

# Registry structure: list of lists
# Each modality has named entries with: package, loader/dataset, class, description

.registry <- list(
  
  # ============================================================================
  # Single-cell RNA-seq datasets
  # ============================================================================
  scrna = list(
    zeisel = list(
      package = "scRNAseq",
      loader = function() scRNAseq::ZeiselBrainData(),
      class = "SingleCellExperiment",
      description = "Mouse brain single-cell RNA-seq (3,005 cells, Zeisel et al.)"
    ),
    darmanis = list(
      package = "scRNAseq",
      loader = function() scRNAseq::DarmanisBrainData(),
      class = "SingleCellExperiment",
      description = "Human brain single-cell RNA-seq (466 cells, Darmanis et al.)"
    ),
    pbmc = list(
      package = "scRNAseq",
      loader = function() scRNAseq::MairPBMCData(),
      class = "SingleCellExperiment",
      description = "Human PBMC CITE-seq data (Mair et al.)"
    ),
    baron = list(
      package = "scRNAseq",
      loader = function() scRNAseq::BaronPancreasData("human"),
      class = "SingleCellExperiment",
      description = "Human pancreas single-cell RNA-seq (Baron et al.)"
    ),
    muraro = list(
      package = "scRNAseq",
      loader = function() scRNAseq::MuraroPancreasData(),
      class = "SingleCellExperiment",
      description = "Human pancreas single-cell RNA-seq (Muraro et al.)"
    ),
    campbell = list(
      package = "scRNAseq",
      loader = function() scRNAseq::CampbellBrainData(),
      class = "SingleCellExperiment",
      description = "Mouse hypothalamus single-cell RNA-seq (Campbell et al.)"
    )
  ),
  
  # ============================================================================
  # Bulk RNA-seq datasets
  # ============================================================================
  bulk = list(
    airway = list(
      package = "airway",
      loader = function() {
        .require_pkg("airway")
        data("airway", package = "airway", envir = environment())
        get("airway", envir = environment())
      },
      class = "RangedSummarizedExperiment",
      description = "Airway smooth muscle RNA-seq (8 samples, dexamethasone treatment)"
    ),
    tcga_brca = list(
      package = "curatedTCGAData",
      loader = function() {
        .require_pkg("curatedTCGAData")
        curatedTCGAData::curatedTCGAData(
          diseaseCode = "BRCA",
          assays = "RNASeq2GeneNorm",
          version = "2.0.1",
          dry.run = FALSE
        )
      },
      class = "MultiAssayExperiment",
      description = "TCGA Breast Cancer RNA-seq (curatedTCGAData)"
    ),
    tcga_luad = list(
      package = "curatedTCGAData",
      loader = function() {
        .require_pkg("curatedTCGAData")
        curatedTCGAData::curatedTCGAData(
          diseaseCode = "LUAD",
          assays = "RNASeq2GeneNorm",
          version = "2.0.1",
          dry.run = FALSE
        )
      },
      class = "MultiAssayExperiment",
      description = "TCGA Lung Adenocarcinoma RNA-seq (curatedTCGAData)"
    ),
    tcga_paad = list(
      package = "curatedTCGAData",
      loader = function() {
        .require_pkg("curatedTCGAData")
        curatedTCGAData::curatedTCGAData(
          diseaseCode = "PAAD",
          assays = "RNASeq2GeneNorm",
          version = "2.0.1",
          dry.run = FALSE
        )
      },
      class = "MultiAssayExperiment",
      description = "TCGA Pancreatic Cancer RNA-seq (curatedTCGAData)"
    )
  ),
  
  # ============================================================================
  # Spatial transcriptomics datasets
  # ============================================================================
  spatial = list(
    dlpfc = list(
      package = "STexampleData",
      loader = function() {
        .require_pkg("STexampleData")
        STexampleData::Visium_humanDLPFC()
      },
      class = "SpatialExperiment",
      description = "Human dorsolateral prefrontal cortex Visium data (Maynard et al.)"
    ),
    mouse_embryo = list(
      package = "STexampleData",
      loader = function() {
        .require_pkg("STexampleData")
        STexampleData::seqFISH_mouseEmbryo()
      },
      class = "SpatialExperiment",
      description = "Mouse embryo seqFISH data (Lohoff et al.)"
    ),
    brain_libd = list(
      package = "spatialLIBD",
      loader = function() {
        .require_pkg("spatialLIBD")
        spatialLIBD::fetch_data(type = "spe")
      },
      class = "SpatialExperiment",
      description = "LIBD Human DLPFC Visium (47,681 spots, spatialLIBD)"
    ),
    ovarian = list(
      package = "STexampleData",
      loader = function() {
        .require_pkg("STexampleData")
        STexampleData::Visium_mouseCoronal()
      },
      class = "SpatialExperiment",
      description = "Mouse brain coronal section Visium (10x Genomics)"
    )
  ),
  
  # ============================================================================
  # ATAC-seq datasets
  # ============================================================================
  atac = list(
    tcga_brca = list(
      package = "curatedTCGAData",
      loader = function() {
        .require_pkg("curatedTCGAData")
        # Note: TCGA ATAC data available through GDC, using methylation as proxy
        curatedTCGAData::curatedTCGAData(
          diseaseCode = "BRCA",
          assays = "Methylation",
          version = "2.0.1",
          dry.run = FALSE
        )
      },
      class = "MultiAssayExperiment",
      description = "TCGA Breast Cancer Methylation (as ATAC proxy)"
    )
  ),
  
  # ============================================================================
  # CyTOF / Mass Cytometry datasets
  # ============================================================================
  cytof = list(
    levine32 = list(
      package = "HDCytoData",
      loader = function() {
        .require_pkg("HDCytoData")
        HDCytoData::Levine_32dim_SE()
      },
      class = "SummarizedExperiment",
      description = "CyTOF bone marrow (32 markers, 265k cells, Levine et al.)"
    ),
    levine13 = list(
      package = "HDCytoData",
      loader = function() {
        .require_pkg("HDCytoData")
        HDCytoData::Levine_13dim_SE()
      },
      class = "SummarizedExperiment",
      description = "CyTOF bone marrow (13 markers, Levine et al.)"
    ),
    samusik = list(
      package = "HDCytoData",
      loader = function() {
        .require_pkg("HDCytoData")
        HDCytoData::Samusik_01_SE()
      },
      class = "SummarizedExperiment",
      description = "CyTOF mouse bone marrow (39 markers, 86k cells, Samusik et al.)"
    ),
    bodenmiller = list(
      package = "HDCytoData",
      loader = function() {
        .require_pkg("HDCytoData")
        HDCytoData::Bodenmiller_BCR_XL_SE()
      },
      class = "SummarizedExperiment",
      description = "CyTOF BCR-XL stimulation (16 samples, Bodenmiller et al.)"
    ),
    krieg = list(
      package = "HDCytoData",
      loader = function() {
        .require_pkg("HDCytoData")
        HDCytoData::Krieg_Anti_PD_1_SE()
      },
      class = "SummarizedExperiment",
      description = "CyTOF melanoma anti-PD-1 response (20 samples, Krieg et al.)"
    )
  ),
  
  # ============================================================================
  # ChIP-seq datasets
  # ============================================================================
  chipseq = list(
    tamoxifen = list(
      package = "DiffBind",
      loader = function() {
        .require_pkg("DiffBind")
        data("tamoxifen_counts", package = "DiffBind", envir = environment())
        dba_obj <- get("tamoxifen", envir = environment())
        # Convert to SummarizedExperiment
        DiffBind::dba(dba_obj, bSummarizedExperiment = TRUE)
      },
      class = "SummarizedExperiment",
      description = "ER ChIP-seq tamoxifen resistance (11 samples, 5 cell lines)"
    ),
    tamoxifen_peaks = list(
      package = "DiffBind",
      loader = function() {
        .require_pkg("DiffBind")
        data("tamoxifen_peaks", package = "DiffBind", envir = environment())
        get("tamoxifen", envir = environment())
      },
      class = "DBA",
      description = "ER ChIP-seq tamoxifen peaks only (DiffBind object)"
    )
  ),
  
  # ============================================================================
  # Genotype / Variant datasets
  # ============================================================================
  genotype = list(
    chr22 = list(
      package = "VariantAnnotation",
      loader = function() {
        .require_pkg("VariantAnnotation")
        fl <- system.file("extdata", "chr22.vcf.gz", package = "VariantAnnotation")
        VariantAnnotation::readVcf(fl, "hg19")
      },
      class = "CollapsedVCF",
      description = "1000 Genomes chr22 variants (10k variants, 5 samples)"
    ),
    structural = list(
      package = "VariantAnnotation",
      loader = function() {
        .require_pkg("VariantAnnotation")
        fl <- system.file("extdata", "structural.vcf", package = "VariantAnnotation")
        VariantAnnotation::readVcf(fl, "hg19")
      },
      class = "CollapsedVCF",
      description = "Structural variants example"
    )
  ),
  
  # ============================================================================
  # Perturbation / CRISPR screen datasets
  # ============================================================================
  perturb = list(
    # Note: Most perturb-seq data requires custom download
    # These are placeholder entries pointing to scRNAseq alternatives
    k562_crispr = list(
      package = "scRNAseq",
      loader = function() {
        .require_pkg("scRNAseq")
        # Placeholder - using a scRNAseq dataset
        # Real perturb-seq would come from scPerturb or similar
        message("Note: For full perturb-seq data, see scPerturb database")
        scRNAseq::LunSpikeInData()
      },
      class = "SingleCellExperiment",
      description = "Spike-in dataset (placeholder for perturb-seq)"
    )
  ),
  
  # ============================================================================
  # Multi-omics datasets
  # ============================================================================
  multiomics = list(
    tcga_brca_multi = list(
      package = "curatedTCGAData",
      loader = function() {
        .require_pkg("curatedTCGAData")
        curatedTCGAData::curatedTCGAData(
          diseaseCode = "BRCA",
          assays = c("RNASeq2GeneNorm", "Methylation", "RPPAArray"),
          version = "2.0.1",
          dry.run = FALSE
        )
      },
      class = "MultiAssayExperiment",
      description = "TCGA BRCA multi-omics (RNA-seq + Methylation + RPPA)"
    ),
    tcga_acc_multi = list(
      package = "curatedTCGAData",
      loader = function() {
        .require_pkg("curatedTCGAData")
        curatedTCGAData::curatedTCGAData(
          diseaseCode = "ACC",
          assays = c("RNASeq2GeneNorm", "RPPAArray"),
          version = "2.0.1",
          dry.run = FALSE
        )
      },
      class = "MultiAssayExperiment",
      description = "TCGA Adrenocortical Carcinoma (RNA-seq + RPPA)"
    )
  )
)

#' Get all available modalities
#' @return Character vector of modality names
#' @keywords internal
.get_modalities <- function() {
  names(.registry)
}

#' Get all datasets for a modality
#' @param modality Modality name
#' @return Named list of dataset entries
#' @keywords internal
.get_datasets <- function(modality) {
  if (!modality %in% names(.registry)) {
    stop(sprintf("Unknown modality: %s. Available: %s",
                 modality, paste(.get_modalities(), collapse = ", ")))
  }
  .registry[[modality]]
}

#' Get a specific dataset entry
#' @param modality Modality name
#' @param name Dataset name
#' @return Dataset entry list
#' @keywords internal
.get_entry <- function(modality, name) {
  datasets <- .get_datasets(modality)
  if (!name %in% names(datasets)) {
    stop(sprintf("Unknown dataset '%s' for modality '%s'. Available: %s",
                 name, modality, paste(names(datasets), collapse = ", ")))
  }
  datasets[[name]]
}
