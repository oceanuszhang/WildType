"""Gene symbols from Ollie's real Embark Veterinary Report
(data/ollie_embark/Ollie_HealthReport.pdf, 255 conditions tested,
extracted 2026-08-15). Used ONLY as a cross-validation aid — to flag when
WildType's own from-scratch scan lands on a locus Embark independently
also tests, so the two can be compared on the same gene. WildType's
pipeline never reads Embark's report or its risk calls during its own
investigation; this list is checked afterward, for comparison purposes,
not fed into the agent's reasoning.

This is Ollie/dog-specific and only meaningful for this cross-check — the
species-agnostic pipeline itself has no dependency on this file.
"""
from __future__ import annotations

EMBARK_PANEL_GENES: frozenset[str] = frozenset(
    {
        "PRCD", "BEST1", "CUBN", "NHEJ1", "SLC37A2", "CNGB3", "OLFM3", "HSF4", "LAMB3",
        "ABCB1", "CLCN1", "CLN5", "CLN6", "CLN8", "ADAMTS17", "FAM20C", "FAM134B",
        "VPS13B", "SLC2A9", "XDH", "MANBA", "FGF4", "GPT", "APRT", "GDNF-AS", "SLC19A3",
        "NDRG1", "GFAP", "EDA", "RHO", "IGFBP5", "GP9", "MSTN", "SPTB", "FUCA1", "ITGB2",
        "FERMT3", "SERAC1", "YARS2", "PTPLA", "VLDLR", "ITGA10", "ADAMTS20", "DLX6", "C3",
        "NSDHL", "TPO", "SLC5A5", "TUBB1", "COLQ", "CHAT", "CHRNE", "LRIT3", "RPE65",
        "SLC3A1", "SLC7A9", "CNGA3", "MYO7A", "SOD1A", "SBF2", "MTRM13", "MIA3", "INPP5E",
        "RBM20", "PDK4", "TTN", "PRKG2", "FAM83H", "COL7A1", "LOXHD1", "EPS8L2", "SEL1L",
        "ADAMTS2", "ENAM", "BCAN", "DNM1", "F7", "F11", "COL4A4", "FAN1", "MFN2", "ITGA2B",
        "GALC", "G6PC", "AGL", "PFKM", "GLB1", "HEXA", "HEXB", "SLC4A3", "TTC8", "F8", "F9",
        "RAB24", "FAM83G", "DSG1", "SUV39H2", "VDR", "CAT", "FNIP2", "ALPL", "NIPAL4",
        "ASPRV1", "SLC27A4", "KRT10", "PNPLA1", "SLC25A12", "BIN1", "ACSL5", "LAMA3",
        "LGI2", "RAB3GAP1", "DIRAS1", "L2HGDH", "ATG4D", "RAPGEF6", "CAPN1", "ATP13A2",
        "ARHGEF10", "GJA9", "MKLN1", "TSEN54", "PLG", "SGCD", "SGCA", "KCNQ1", "LEPREL1",
        "CHST6", "RYR1", "MYH9", "CYB5R3", "RBP4", "NAGLU", "SGSH", "ARSB", "GUSB", "DMD",
        "ADAMTSL2", "HCRTR2", "NEB", "SPTBN2", "ATF2", "LAMP3", "VPS11", "TECPR2", "PPT1",
        "CTSD", "TPP1", "MFSD8", "ARSG", "SLC45A2", "COL9A2", "SLC13A1", "COL1A2",
        "SERPINH1", "COL1A1", "P2Y12", "KRT16", "PIGN", "AMHR2", "POU1F1", "TMEM16F",
        "PKD1", "GAA", "KLKB1", "NME5", "CCDC39", "AGXT", "ADAMTS10", "SAG", "IFT122",
        "BBS2", "CNGA1", "PDE6B", "RPGRIP1", "CNGB1", "FAM161A", "GH1", "NPHS1", "PDP1",
        "PKLR", "AKNA", "FLCN", "SIX6", "PRKDC", "RAG1", "PLP1", "MTBP", "COL11A2", "PKP1",
        "SCN8A", "KCNJ10", "ATP1B2", "ABCA4", "ALDH5A1", "RASGRP1", "COL6A3", "COL6A1",
        "PTPRQ", "VWF", "COL4A5", "MTM1", "RPGR", "IL2RG",
    }
)
