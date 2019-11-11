
#' @title List of housekeeping genes
#'
#' @description The genes contained in the list are ubiquitously expressed
#' in most cells.
#'
#' @export
SCANDAL_HOUSEKEEPING_GENES_LIST <- c("ACTB", "B2M", "HNRPLL", "HPRT", "PSMB2", "PSMB4", "PPIA", "PRPS1", "PRPS1L1", "PRPS1L3", "PRPS2", "PRPSAP1", "PRPSAP2",
                                     "RPL10", "RPL10A", "RPL10L", "RPL11", "RPL12", "RPL13", "RPL14", "RPL15", "RPL17", "RPL18", "RPL19", "RPL21", "RPL22",
                                     "RPL22L1", "RPL23", "RPL24", "RPL26", "RPL27", "RPL28", "RPL29", "RPL3", "RPL30", "RPL32", "RPL34", "RPL35", "RPL36", "RPL37",
                                     "RPL38", "RPL39", "RPL39L", "RPL3L", "RPL4", "RPL41", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL7L1", "RPL8", "RPL9", "RPLP0",
                                     "RPLP1", "RPLP2", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS17", "RPS18", "RPS19", "RPS20",
                                     "RPS21", "RPS24", "RPS25", "RPS26", "RPS27", "RPS27A", "RPS27L", "RPS28", "RPS29", "RPS3", "RPS3A", "RPS4X", "RPS5", "RPS6",
                                     "RPS6KA1", "RPS6KA2", "RPS6KA3", "RPS6KA4", "RPS6KA5", "RPS6KA6", "RPS6KB1", "RPS6KB2", "RPS6KC1", "RPS6KL1", "RPS7", "RPS8",
                                     "RPS9", "RPSA", "TRPS1", "UBB")

#'
#' @title List of G1S cell-cycle markers
#'
#' @export
SCANDAL_G1S_MARKERS <- c("RRM2", "TYMS" , "UBE2T", "CDK1", "HMGB2", "MAD2L1", "PCNA", "UBE2C", "PBK", "TOP2A", "NUSAP1", "KIAA0101", "HIST1H4C",
                         "MLF1IP", "GMNN", "BIRC5", "FAM64A", "RNASEH2A", "MELK", "CENPK", "PTTG1", "TK1", "TPX2", "TMEM106C", "CDCA5", "CKS1B",
                         "CDC45", "MCM3", "CENPM", "AURKB", "PKMYT1", "KIF22","MCM4", "ASF1B", "GINS2", "MCM2", "NUF2", "CDKN3", "GGH", "NDC80",
                         "FEN1", "RRM1", "PRC1" , "DUT", "RAD51AP1", "CKS2", "MCM7", "CCNE2", "ZWINT")

#'
#' @title List of G2M cell-cycle markers
#'
#' @export
SCANDAL_G2M_MARKERS <- c("CCNB1", "UBE2C", "PTTG1", "CDC20", "CCNB2", "TOP2A", "FAM64A", "NUSAP1", "CDKN3", "PBK", "PLK1", "HMGB2", "TPX2", "BIRC5",
                         "MAD2L1", "PRC1", "NUF2", "UBE2T", "CDK1", "CKS2", "CCNA2", "CKAP2", "KNSTRN", "RACGAP1", "CDCA3", "TROAP", "KIF2C", "AURKA",
                         "CENPF", "KPNA2", "KIF20A", "ECT2", "BUB1", "CDCA8", "BUB1B", "TACC3", "NDC80", "TTK", "TUBA1C", "NCAPD2", "ARL6IP1", "KIF4A",
                         "CKAP2L", "MZT1", "KIFC1", "KIF22", "TYMS", "SPAG5", "ANP32E", "KIF11", "PSRC1", "TUBB4B", "SMC4", "MXD3", "CDC25B", "OIP5",
                         "GGH", "REEP4", "FOXM1", "TMPO", "GPSM2", "HMGB3", "ARHGAP11A", "RANGAP1", "H2AFZ")

#'
#' @title List of oligodendrocyte markers
#'
#' @export
SCANDAL_OLIGODENDROCYTE_MARKERS <- c("MOBP", "OPALIN", "CLDN11", "PLP1", "TF", "MBP", "MOG", "MAG")

#'
#' @title List of macrophage markers
#'
#' @export
SCANDAL_MACROPHAGE_MARKERS <- c("CD14", "AIF1", "CD163", "TYROBP", "CSF1R", "CD74", "CD68")

#'
#' @title List of AC programs markers from IDH-O dataset (Oligodendroglioma)
#'
#' @export
SCANDAL_IDH_O_AC_MARKERS <- c("APOE", "SPARCL1", "SPOCK1", "CRYAB", "ALDOC", "CLU", "EZR", "SORL1", "MLC1", "ABCA1",
                              "ATP1B2", "PAPLN", "CA12", "BBOX1", "RGMA", "AGT", "EEPD1", "CST3", "SSTR2", "SOX9",
                              "RND3", "EDNRB", "GABRB1", "PLTP", "JUNB", "DKK3", "ID4", "ADCYAP1R1", "GLUL", "EPAS1",
                              "PFKFB3", "ANLN", "HEPN1", "CPE", "RASL10A", "SEMA6A", "ZFP36L1", "HEY1", "PRLHR", "TACR1",
                              "JUN", "GADD45B", "SLC1A3", "CDC42EP4", "MMD2", "CPNE5", "CPVL", "RHOB", "NTRK2", "CBS")

#'
#' @title List of OC programs markers from IDH-O dataset (Oligodendroglioma)
#'
#' @export
SCANDAL_IDH_O_OC_MARKERS <- c("LMF1", "OLIG1", "SNX22", "POLR2F", "LPPR1", "GPR17", "DLL3", "ANGPTL2", "SOX8", "RPS2",
                              "FERMT1", "PHLDA1", "RPS23", "NEU4", "SLC1A1", "LIMA1", "ATCAY", "SERINC5", "CDH13", "CXADR",
                              "LHFPL3", "ARL4A", "SHD", "RPL31", "GAP43", "IFITM10", "SIRT2", "OMG", "RGMB", "HIPK2",
                              "APOD", "NPPA", "EEF1B2", "RPS17L", "FXYD6", "MYT1", "RGR", "OLIG2", "ZCCHC24", "MTSS1",
                              "GNB2L1", "C17orf76-AS1", "ACTG1", "EPN2", "PGRMC1", "TMSB10", "NAP1L1", "EEF2", "MIAT", "CDHR1")

#'
#' @title List of Stemness programs markers from IDH-O dataset (Oligodendroglioma)
#'
#' @export
SCANDAL_IDH_O_STEMNESS_MARKERS <- c("SOX4", "CCND2", "SOX11", "RBM6", "HNRNPH1", "HNRNPL", "PTMA", "TRA2A", "SET", "C6orf62",
                                    "PTPRS", "CHD7", "CD24", "H3F3B", "C14orf23", "NFIB", "SRGAP2C", "STMN2", "SOX2", "TFDP2",
                                    "CORO1C", "EIF4B", "FBLIM1", "SPDYE7P", "TCF4", "ORC6", "SPDYE1", "NCRUPAR", "BAZ2B", "NELL2",
                                    "OPHN1", "SPHKAP", "RAB42", "LOH12CR2", "ASCL1", "BOC", "ZBTB8A", "ZNF793", "TOX3", "EGFR",
                                    "PGM5P2", "EEF1A1", "MALAT1", "TATDN3", "CCL5", "EVI2A", "LYZ", "POU5F1", "FBXO27", "CAMK2N1")

#'
#' @title List of G1S (cell cycle) programs markers from IDH-O dataset (Oligodendroglioma)
#'
#' @export
SCANDAL_IDH_O_G1S_MARKERS <- c("MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
                               "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1",
                               "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2",
                               "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1")

#'
#' @title List of G2M (cell cycle) programs markers from IDH-O dataset (Oligodendroglioma)
#'
#' @export
SCANDAL_IDH_O_G2M_MARKERS <- c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2",
                               "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2",
                               "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "HJURP", "CDCA3",
                               "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8")

#'
#' @title List of AC programs markers from IDH-A dataset (Astrocytoma)
#'
#' @export
SCANDAL_IDH_A_AC_MARKERS <- c("APOE", "SPARCL1", "VIM", "ID4", "TIMP3", "EDNRB", "MLC1", "ID3", "CLU", "TNC",
                              "ZFP36L1", "ARHGEF26", "ATP1B2", "AGT", "RGMA", "JUN", "PFKFB3", "EZR", "SLC1A3", "ALDOC",
                              "JUNB", "ATP1A2", "DTNA", "ZFP36", "SOX9", "TRIL", "NDRG2", "NMB", "GFAP", "SLC1A2",
                              "RFX4", "MALAT1", "LRIG1", "FOS", "EGR1", "STK17B", "FOSB", "ATF3", "ABCA1", "ADCYAP1R1",
                              "GLUL", "IER2", "ZFP36L2", "ADHFE1", "MSI2", "CPE", "KLF6", "DOCK7", "IRF2BP2", "SPRY2")

#'
#' @title List of OC programs markers from IDH-A dataset (Astrocytoma)
#'
#' @export
SCANDAL_IDH_A_OC_MARKERS <- c("OLIG1", "NEU4", "GPR17", "SLC1A1", "ATCAY", "SIRT2", "APOD", "MYT1", "OLIG2", "TMEFF2",
                              "OMG", "ELMO1", "RTKN", "HIP1R", "TNR", "RPSA", "MEGF11", "EVI2A", "OPCML", "LHFPL3",
                              "RAB33A", "GRIA4", "SERINC5", "NXPH1", "BIN1", "BMP4", "EHD3", "GNAI1", "CSPG4", "DSCAM",
                              "GALNT13", "ZDHHC9", "ABCG1", "FKBP1A", "LRRN1", "ST8SIA3", "DNM3", "RAPGEF4", "CNP", "PDGFRA",
                              "PTGDS", "CHGA", "BCAS1", "PLXNB3", "NFASC", "SLC44A1", "GNG4", "PHLDB1", "CD82", "PRKCZ")

#'
#' @title List of Stemness programs markers from IDH-A dataset (Astrocytoma)
#'
#' @export
SCANDAL_IDH_A_STEMNESS_MARKERS <- c("SOX4", "DCX", "IGFBPL1", "SOX11", "TCF4", "NREP", "RND3", "CCND2", "MIAT", "CAMK2N1",
                                    "STMN4", "STMN1", "MYT1L", "HN1", "RNF122", "PROX1", "KLHDC8A", "ELAVL4", "NMNAT2", "TUBB",
                                    "ROBO1", "NELL2", "MLLT11", "CELF4", "POU3F2", "H3F3B", "ENC1", "GNG2", "ACOT7", "AKT3",
                                    "ARL4C", "FNBP1L", "VOPP1", "TOX3", "TUBB3", "SCG2", "TMSB15A", "TFDP2", "TMSB4X", "CDC42",
                                    "STMN2", "KCTD13", "RPH3A", "KIF5C", "NFIX", "CALM1", "TNPO2", "BOC", "KLHL13", "PGAP1")

#'
#' @title List of G1S (cell cycle) programs markers from IDH-A dataset (Astrocytoma)
#'
#' @export
SCANDAL_IDH_A_G1S_MARKERS <- c("MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
                               "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1",
                               "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2",
                               "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1")

#'
#' @title List of G2M (cell cycle) programs markers from IDH-A dataset (Astrocytoma)
#'
#' @export
SCANDAL_IDH_A_G2M_MARKERS <- c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2",
                               "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2",
                               "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "HJURP", "CDCA3",
                               "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8")




