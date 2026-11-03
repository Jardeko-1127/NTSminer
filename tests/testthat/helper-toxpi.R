.toxpi_fixture_epi <- function() {
  data.frame(
    SMILES = c("CCO", "CC", "CCC"),
    episuite_Name = c("ethanol", "ethane", "propane"),
    LogKow_est = c(-0.31, 1.81, 2.42),
    LogP_pred = c(-0.31, 1.81, 2.42),
    LogKoc_est = c(1.5, 2.1, 2.8),
    LogKoa_est = c(4.0, 5.1, 5.5),
    LogBCF = c(0.5, 1.2, 1.8),
    LogBAF = c(0.6, 1.4, 2.0),
    Biowin1 = c(0.8, 0.5, 0.4),
    Biowin2 = c(0.7, 0.4, 0.3),
    Biowin3 = c(4.0, 3.0, 2.0),
    Biowin4 = c(0.6, 0.45, 0.3),
    Biowin5 = c(0.85, 0.55, 0.45),
    Biowin6 = c(0.65, 0.4, 0.2),
    Biowin7 = c(0.4, 0.3, 0.2),
    Fish_LC50 = c(100, 10, 1),
    Daphnid_LC50 = c(80, 8, 0.8),
    Algae_LC50 = c(50, 5, 0.5),
    Fish_chV = c(10, 1, 0.1),
    Daphnid_chV = c(8, 0.8, 0.08),
    Algae_chV = c(5, 0.5, 0.05),
    stringsAsFactors = FALSE
  )
}

.toxpi_fixture_admet <- function() {
  data.frame(
    SMILES = c("CCO", "CC", "CCC"),
    LD50_Zhu = c(1.5, 0.5, -0.2),
    DILI = c(0.1, 0.4, 0.8),
    ClinTox = c(0.05, 0.3, 0.6),
    hERG = c(0.02, 0.25, 0.7),
    AMES = c(0.01, 0.2, 0.5),
    Carcinogens_Lagunin = c(0.02, 0.15, 0.4),
    `NR-AR` = c(0.1, 0.3, 0.6),
    `NR-AR-LBD` = c(0.05, 0.2, 0.5),
    `NR-AhR` = c(0.08, 0.25, 0.55),
    `NR-Aromatase` = c(0.03, 0.1, 0.2),
    `NR-ER` = c(0.04, 0.15, 0.35),
    `NR-ER-LBD` = c(0.06, 0.18, 0.4),
    `NR-PPAR-gamma` = c(0.02, 0.12, 0.3),
    `SR-ARE` = c(0.07, 0.22, 0.45),
    `SR-ATAD5` = c(0.05, 0.18, 0.38),
    `SR-HSE` = c(0.04, 0.16, 0.33),
    `SR-MMP` = c(0.09, 0.28, 0.5),
    `SR-p53` = c(0.03, 0.14, 0.29),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
