from dataclasses import dataclass, field
from typing import Dict, List, Optional

import pandas as pd


@dataclass
class Individualization(object):
    """
    Individualize a mechanistic model by incorporating gene expression levels.

    Attributes
    ----------
    parameters : List[str]
        List of model parameters.

    species : List[str]
        List of model species.

    tpm_values : str
        Path to (1) RLE-normalized and (2) post ComBat TPM values (CSV-formatted).

    structure : Dict[str, List[str]]
        Pairs of proteins and their related genes.

    read_csv_kwargs : dict, optional
        Arguments to pandas.read_csv.

    prefix : str (default: "w_")
        Prefix of weighting factors on gene expression levels.
    """

    parameters: List[str]
    species: List[str]
    tpm_values: str
    structure: Dict[str, List[str]]
    read_csv_kwargs: Optional[dict] = field(default=None)
    prefix: str = field(default="w_", init=False)

    def __post_init__(self) -> None:
        kwargs = self.read_csv_kwargs
        if kwargs is None:
            kwargs = {}
        kwargs.setdefault("index_col", 2)
        self._tpm_rle_postComBat: pd.DataFrame = pd.read_csv(self.tpm_values, **kwargs)
        del kwargs

    @property
    def tpm_rle_postComBat(self):
        return self._tpm_rle_postComBat

    def _get_tpm(self, gene: str, id: str) -> float:
        return self.tpm_rle_postComBat.at[gene, id]

    def _calculate_weighted_sum(
        self,
        id: str,
        x: List[float],
    ) -> Dict[str, float]:
        """
        Incorporate gene expression levels in the model.

        Parameters
        ----------
        id : str
            CCLE_ID or TCGA_ID.

        x : List[float]
            Model parameters.

        Returns
        -------
        weighted_sum : Dict[str, float]
            Estimated protein levels after incorporating transcriptomic data.
        """
        weighted_sum = dict.fromkeys(self.structure, 0.0)
        for (protein, genes) in self.structure.items():
            for gene in genes:
                weighted_sum[protein] += x[
                    self.parameters.index(self.prefix + gene)
                ] * self._get_tpm(gene, id)
        return weighted_sum

    def as_reaction_rate(
        self,
        id: str,
        x: List[float],
        param_name: str,
        protein: str,
    ) -> float:
        """
        Gene expression levels are incorporated as a reaction rate.

        Parameters
        ----------
        id : str
            CCLE_ID or TCGA_ID.

        x : List[float]
            List of parameter values.

        param_name : str
            Parameter incorporating gene_expression_data.

        protein: str
            Protein involved in the reaction.

        Returns
        -------
        param_value : float
        """
        weighted_sum = self._calculate_weighted_sum(id, x)
        param_value = x[self.parameters.index(param_name)]
        param_value *= weighted_sum[protein]
        return param_value

    def as_initial_conditions(
        self,
        id: str,
        x: List[float],
        y0: List[float],
    ) -> List[float]:
        """
        Gene expression levels are incorporated as initial conditions.

        Parameters
        ----------
        id : str
            CCLE_ID or TCGA_ID.

        x : List[float]
            List of parameter values.

        y0 : List[float]
            List of initial values.

        Returns
        -------
        y0 (individualized) : List[float]
            Cell-line- or patient-specific initial conditions.
        """
        weighted_sum = self._calculate_weighted_sum(id, x)
        for protein in self.structure.keys():
            y0[self.species.index(protein)] *= weighted_sum[protein]
        return y0
