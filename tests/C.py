REQUIREMENTS = """\
import os
from pasmopy import Individualization

from . import __path__
"""


INDIVIDUALIZATION = """\
incorporating_gene_expression_levels = Individualization(
    parameters=C.NAMES,
    species=V.NAMES,
    transcriptomic_data="https://raw.githubusercontent.com/pasmopy/breast_cancer/master/transcriptomic_data/TPM_RLE_postComBat_BRCA_BREAST.tar.xz",
    gene_expression={
        "ErbB1": ["EGFR"],
        "ErbB2": ["ERBB2"],
        "ErbB3": ["ERBB3"],
        "ErbB4": ["ERBB4"],
        "Grb2": ["GRB2"],
        "Shc": ["SHC1", "SHC2", "SHC3", "SHC4"],
        "RasGAP": ["RASA1", "RASA2", "RASA3"],
        "PI3K": ["PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"],
        "PTEN": ["PTEN"],
        "SOS": ["SOS1", "SOS2"],
        "Gab1": ["GAB1"],
        "RasGDP": ["HRAS", "KRAS", "NRAS"],
        "Raf": ["ARAF", "BRAF", "RAF1"],
        "MEK": ["MAP2K1", "MAP2K2"],
        "ERK": ["MAPK1", "MAPK3"],
        "Akt": ["AKT1", "AKT2"],
        "PTP1B": ["PTPN1"],
        "GSK3b": ["GSK3B"],
        "DUSP": ["DUSP5", "DUSP6", "DUSP7"],
        "cMyc": ["MYC"],
    },
    read_csv_kws={"index_col": "Description"},
)
"""


INCORPORATION = """\
        x[C.V291] = incorporating_gene_expression_levels.as_reaction_rate(
            __path__[0].split(os.sep)[-1], x, "V291", "DUSP"
        )
        x[C.V310] = incorporating_gene_expression_levels.as_reaction_rate(
            __path__[0].split(os.sep)[-1], x, "V310", "cMyc"
        )
        y0 = incorporating_gene_expression_levels.as_initial_conditions(
            __path__[0].split(os.sep)[-1], x, y0
        )
"""
