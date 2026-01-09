.registry <- list(
  
  scrna = list(
    zeisel_2015 = function() scRNAseq::ZeiselBrainData(),
    darmanis_2015 = function() scRNAseq::DarmanisBrainData(),
    campbell_2017 = function() scRNAseq::CampbellBrainData(),
    tasic_2016 = function() scRNAseq::TasicBrainData(),
    romanov_2017 = function() scRNAseq::RomanovBrainData(),
    muraro_2016 = function() scRNAseq::MuraroPancreasData(),
    segerstolpe_2016 = function() scRNAseq::SegerstolpePancreasData(),
    lawlor_2017 = function() scRNAseq::LawlorPancreasData(),
    grun_2016 = function() scRNAseq::GrunPancreasData(),
    pbmc_2020 = function() scRNAseq::MairPBMCData(),
    kotliarov_2020 = function() scRNAseq::KotliarovPBMCData(),
    macosko_2015 = function() scRNAseq::MacoskoRetinaData(),
    zhong_2018 = function() scRNAseq::ZhongPrefrontalData(),
    bach_2017 = function() scRNAseq::BachMammaryData(),
    ernst_2019 = function() scRNAseq::ErnstSpermatogenesisData(),
    aztekin_2019 = function() scRNAseq::AztekinTailData()
  ),
  
  bulk = list(
    airway_2014 = function() { data("airway", package = "airway", envir = environment()); get("airway", envir = environment()) },
    pasilla_2010 = function() tidySummarizedExperiment::pasilla
  ),
  
  spatial = list(
    dlpfc_2021 = function() STexampleData::Visium_humanDLPFC(),
    mouse_coronal_2021 = function() STexampleData::Visium_mouseCoronal()
  ),
  
  cytof = list(
    levine32_2015 = function() HDCytoData::Levine_32dim_SE(),
    levine13_2015 = function() HDCytoData::Levine_13dim_SE(),
    samusik_2016 = function() HDCytoData::Samusik_01_SE(),
    bodenmiller_2012 = function() HDCytoData::Bodenmiller_BCR_XL_SE()
  )
  
  # chipseq = list(
  #   cstest_2008 = function() { data("cstest", package = "chipseq", envir = environment()); get("cstest", envir = environment()) }
  # ),
  # 
  # atac = list(
  #   buenrostro_2018 = function() {
  #     .require_pkg("scATAC.Explorer")
  #     scATAC.Explorer::queryATAC(accession = "GSE89362")[[1]]
  #   },
  #   satpathy_2019 = function() {
  #     .require_pkg("scATAC.Explorer")
  #     scATAC.Explorer::queryATAC(accession = "GSE129785")[[1]]
  #   }
  # ),
  # 
  # perturb = list(
  #   dixit_2016 = function() {
  #     .require_github("markowetzlab/SCperturb")
  #     list(
  #       counts = get("counts.dixit2016_K562_lowmoi", envir = asNamespace("SCperturb")),
  #       rowdata = get("rowmetadata.dixit2016_K562_lowmoi", envir = asNamespace("SCperturb")),
  #       coldata = get("colmetadata.dixit2016_K562_lowmoi", envir = asNamespace("SCperturb"))
  #     )
  #   },
  #   datlinger_2017 = function() {
  #     .require_github("markowetzlab/SCperturb")
  #     list(
  #       counts = get("counts.datlinger2017_stim", envir = asNamespace("SCperturb")),
  #       rowdata = get("rowmetadata.datlinger2017_stim", envir = asNamespace("SCperturb")),
  #       coldata = get("colmetadata.datlinger2017_stim", envir = asNamespace("SCperturb"))
  #     )
  #   }
  # )
  
)