import os
from dataclasses import dataclass
from typing import List, NamedTuple, Optional

import pandas as pd


@dataclass
class TCGA(object):
    """
    The Cancer Genome Atlas (TCGA).
    """

    pass


class DrugResponse(NamedTuple):
    ccle_cell_line_name: str
    primary_cell_line_name: str
    compound: str
    target: str
    doses: List[float]
    activity_data: List[float]
    activity_sd: List[float]
    num_data: int
    fit_type: str
    ec50: float
    ic50: float
    amax: float
    act_area: float


@dataclass
class CCLE(object):
    """
    Cancer Cell Line Encyclopedia (CCLE).

    Attributes
    ----------
    gene_expression : str (filepath or URL)
        RNA-seq dataset containing the TPM normalized values.

    drug_response : str (filepath or URL)
        Pharmacologic profiles for 24 anticancer drugs across 504 cell lines.
    """

    gene_expression: str = (
        "https://data.broadinstitute.org/ccle/CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz"
    )
    drug_response: str = (
        "https://data.broadinstitute.org/ccle_legacy_data/pharmacological_profiling/"
        "CCLE_NP24.2009_Drug_data_2015.02.24.csv"
    )

    def __post_init__(self):
        self._gene_expression_data: pd.DataFrame = pd.read_table(self.gene_expression, index_col=0)
        self._drug_response_data: pd.DataFrame = pd.read_csv(self.drug_response)

    @property
    def gene_expression_data(self) -> pd.DataFrame:
        return self._gene_expression_data

    @property
    def drug_response_data(self) -> pd.DataFrame:
        return self._drug_response_data

    def save_transcriptomic_data(self):
        os.makedirs("transcriptomic_data", exist_ok=True)
        self.gene_expression_data.to_csv("transcriptomic_data/CCLE_tpm_values.csv")

    def get_drug_response(
        self,
        cell_line: Optional[List[str]] = None,
        compound: Optional[List[str]] = None,
    ) -> List[DrugResponse]:
        """
        Extract drug-response data.

        Parameters
        ----------
        cell_line : list of strings, optional
            List of CCLE Cell Line Names.
        compound : list of strings, optional
            List of drug names.

        Returns
        -------
        drug_response_info : list of DrugResponse
            Drug response information.
        """

        if cell_line is None:
            cell_line = list(self.drug_response_data.loc[:, "CCLE Cell Line Name"])
        if compound is None:
            compound = list(self.drug_response_data.loc[:, "Compound"])
        df = self.drug_response_data[
            (self.drug_response_data["CCLE Cell Line Name"].isin(cell_line))
            & (self.drug_response_data["Compound"].isin(compound))
        ]
        drug_response_info = []
        for i in df.index:
            drug_response_info.append(
                DrugResponse(
                    df.at[i, "CCLE Cell Line Name"],
                    df.at[i, "Primary Cell Line Name"],
                    df.at[i, "Compound"],
                    df.at[i, "Target"],
                    [float(val) for val in df.at[i, "Doses (uM)"].split(",")],
                    [float(val) for val in df.at[i, "Activity Data (median)"].split(",")],
                    [float(val) for val in df.at[i, "Activity SD"].split(",")],
                    df.at[i, "Num Data"],
                    df.at[i, "FitType"],
                    df.at[i, "EC50 (uM)"],
                    df.at[i, "IC50 (uM)"],
                    df.at[i, "Amax"],
                    df.at[i, "ActArea"],
                )
            )
        return drug_response_info
