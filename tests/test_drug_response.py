import os
import shutil
from typing import Final

import pandas as pd

from pasmopy.validation import CancerCellLineEncyclopedia

ErbB_expression_ratio: Final[pd.DataFrame] = pd.read_csv(
    "https://raw.githubusercontent.com/pasmopy/breast_cancer/master/drug_response/data/ErbB_expression_ratio.csv",
    index_col=0,
)

ccle: Final = CancerCellLineEncyclopedia()


def test_dose_response_curve():
    for drug in ["Erlotinib", "Lapatinib"]:
        assert (
            ccle.save_dose_response_curve(
                ErbB_expression_ratio,
                {"value": ["high", "low"]},
                drug,
                labels=["EGFR high", "EGFR low"],
            )
            is None
        )

    for drug in ["Erlotinib", "Lapatinib"]:
        assert os.path.isfile(os.path.join("dose_response", "EGFR", f"{drug}.pdf"))


def test_activity_area():
    for drug in ["Erlotinib", "Lapatinib"]:
        assert (
            ccle.save_activity_area(
                ErbB_expression_ratio,
                {"value": ["high", "low"]},
                drug,
                labels=["EGFR high", "EGFR low"],
            )
            is None
        )
    for drug in ["Erlotinib", "Lapatinib"]:
        assert os.path.isfile(os.path.join("activity_area", "EGFR", f"{drug}.pdf"))


def test_cleanup():
    for dirname in ["dose_response", "activity_area"]:
        if os.path.isdir(dirname):
            shutil.rmtree(dirname)
