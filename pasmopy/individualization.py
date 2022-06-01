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

    transcriptomic_data : str
        Path to normalized gene expression data (CSV-formatted),
        e.g., (1) RLE-normalized and (2) post-ComBat TPM values.
        Below is an example of data table.

        =========== ======== ======== ======== ===
        Description patient1 patient2 patient3 ...
        =========== ======== ======== ======== ===
        gene1       value1,1 value1,2 value1,3 ...
        gene2       value2,1 value2,2 value2,3 ...
        gene3       value3,1 value3,2 value3,3 ...
        ...         ...      ...      ...      ...
        =========== ======== ======== ======== ===

    gene_expression : Dict[str, List[str]]
        Pairs of proteins and their related genes.

    read_csv_kws : dict, optional
        Keyword arguments to pass to ``pandas.read_csv``.

    prefix : str (default: "w_")
        Prefix of weighting factors on gene expression levels.

    Examples
    --------
    ``search_param.py``

    .. code-block:: python

        import os
        import numpy as np
        from pasmopy import Individualization
        from . import __path__
        from .name2idx import C, V
        from .ode import initial_values, param_values

        incorporating_gene_expression_levels = Individualization(
            parameters=C.NAMES,
            species=V.NAMES,
            transcriptomic_data=os.path.join("transcriptomic_data", "TPM_RLE_postComBat_BRCA_BREAST.csv"),
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
            read_csv_kws={"index_col": "Description"}
        )

        ...

        def update(self, indiv):
            x = param_values()
            y0 = initial_values()
            for i, j in enumerate(self.idx_params):
                x[j] = indiv[i]
            for i, j in enumerate(self.idx_initials):
                y0[j] = indiv[i + len(self.idx_params)]
            # As maximal transcription rate
            x[C.V291] = incorporating_gene_expression_levels.as_reaction_rate(
                __path__[0].split(os.sep)[-1], x, "V291", "DUSP"
            )
            x[C.V310] = incorporating_gene_expression_levels.as_reaction_rate(
                __path__[0].split(os.sep)[-1], x, "V310", "cMyc"
            )
            # As initial conditions
            y0 = incorporating_gene_expression_levels.as_initial_conditions(
                __path__[0].split(os.sep)[-1], x, y0
            )

            ...
    """

    parameters: List[str]
    species: List[str]
    transcriptomic_data: str
    gene_expression: Dict[str, List[str]]
    read_csv_kws: Optional[dict] = field(default=None)
    prefix: str = field(default="w_", init=False)

    def __post_init__(self) -> None:
        kwargs = self.read_csv_kws
        if kwargs is None:
            kwargs = {}
        self._expression_level: pd.DataFrame = pd.read_csv(self.transcriptomic_data, **kwargs)

    @property
    def expression_level(self) -> pd.DataFrame:
        return self._expression_level

    def _calculate_weighted_sum(
        self,
        id: str,
        x: List[float],
    ) -> Dict[str, float]:
        """
        Incorporate gene expression levels in the model.

        Returns
        -------
        weighted_sum : Dict[str, float]
            Estimated protein levels after incorporating transcriptomic data.
        """
        weighted_sum = dict.fromkeys(self.gene_expression, 0.0)
        for (protein, genes) in self.gene_expression.items():
            for gene in genes:
                weighted_sum[protein] += (
                    x[self.parameters.index(self.prefix + gene)]
                    * self.expression_level.at[gene, id]
                )
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
            Name of the parameter incorporating gene_expression_data.

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
        for protein in self.gene_expression.keys():
            y0[self.species.index(protein)] *= weighted_sum[protein]
        return y0
