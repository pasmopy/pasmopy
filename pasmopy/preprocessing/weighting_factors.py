import os
from dataclasses import dataclass, field
from typing import Dict, List

from biomass.model_object import ModelObject


@dataclass
class WeightingFactors(object):
    """
    Prepare for adding information about gene expression data to model.

    Attributes
    ----------
    model : :py:class:`biomass.model_object.ModelObject`
        BioMASS model object.
    gene_expression : dict
        Pairs of proteins and their related genes.
    weighting_factors : list of strings
        List of weighting factors.
    prefix : str (default: "w_")
        Prefix of weighting factors on gene expression levels.
    indentation : str (default: 4 spaces)
        How many spaces as indentation.

    Examples
    --------
    >>> import erbb_network
    >>> model = Model(erbb_network.__package__).create()
    >>> gene_expression = {
    ...     "ErbB1": ["EGFR"],
    ...     "ErbB2": ["ERBB2"],
    ...     "ErbB3": ["ERBB3"],
    ...     "ErbB4": ["ERBB4"],
    ...     "Grb2": ["GRB2"],
    ...     "Shc": ["SHC1", "SHC2", "SHC3", "SHC4"],
    ...     "RasGAP": ["RASA1", "RASA2", "RASA3"],
    ...     "PI3K": ["PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"],
    ...     "PTEN": ["PTEN"],
    ...     "SOS": ["SOS1", "SOS2"],
    ...     "Gab1": ["GAB1"],
    ...     "RasGDP": ["HRAS", "KRAS", "NRAS"],
    ...     "Raf": ["ARAF", "BRAF", "RAF1"],
    ...     "MEK": ["MAP2K1", "MAP2K2"],
    ...     "ERK": ["MAPK1", "MAPK3"],
    ...     "Akt": ["AKT1", "AKT2"],
    ...     "PTP1B": ["PTPN1"],
    ...     "GSK3b": ["GSK3B"],
    ...     "DUSP": ["DUSP5", "DUSP6", "DUSP7"],
    ...     "cMyc": ["MYC"],
    ... }
    >>> weighting_factors = WeightingFactors(model, gene_expression)
    >>> weighting_factors.add_to_params()
    >>> weighting_factors.set_search_bounds()
    """

    model: ModelObject
    gene_expression: Dict[str, List[str]]
    weighting_factors: List[str] = field(default_factory=list, init=False)
    prefix: str = field(default="w_", init=False)
    indentation: str = field(default=" " * 4, init=False)

    def add_to_params(self) -> None:
        """
        Add weighting factors to model parameters.
        """
        for genes in self.gene_expression.values():
            for gene in genes:
                if self.prefix + gene not in self.model.parameters:
                    self.weighting_factors.append(self.prefix + gene)
        if self.weighting_factors:
            with open(
                os.path.join(self.model.path, "name2idx", "parameters.py"),
                mode="r",
                encoding="utf-8",
            ) as f:
                lines = f.readlines()
            for line_num, line in enumerate(lines):
                if line.startswith("NUM: int"):
                    lines[line_num] = "NAMES.extend(\n"
                    lines[line_num] += (
                        f"{self.indentation}[\n"
                        + f'{2 * self.indentation}"'
                        + f'",\n{2 * self.indentation}"'.join(self.weighting_factors)
                        + f'",\n{self.indentation}]\n'
                    )
                    lines[line_num] += ")\n\nNUM: int = len(NAMES)\n"
            with open(
                os.path.join(self.model.path, "name2idx", "parameters.py"),
                mode="w",
                encoding="utf-8",
            ) as f:
                f.writelines(lines)

    def set_search_bounds(self, lb: float = 0.01, ub: float = 100.0) -> None:
        """
        Set search bounds for weighting factors.

        Parameters
        ----------
        lb : float (default: 0.01)
            Lower bound.
        ub : float (default: 100.0)
            Upper bound.
        """
        if self.weighting_factors:
            search_bounds = [
                f"search_rgn[:, C.{wf}] = [{lb}, {ub}]" for wf in self.weighting_factors
            ]
            with open(
                os.path.join(self.model.path, "search_param.py"),
                mode="r",
                encoding="utf-8",
            ) as f:
                lines = f.readlines()
            for line_num, line in enumerate(lines):
                if line.startswith(f"{2 * self.indentation}]") and lines[line_num + 3].startswith(
                    f"{2 * self.indentation}self.idx_initials = "
                ):
                    lines[line_num] = (
                        f"{3 * self.indentation}"
                        + f",\n{3 * self.indentation}".join(
                            map(lambda _s: "C." + _s, self.weighting_factors)
                        )
                        + f",\n{2 * self.indentation}]\n"
                    )
                elif line.startswith(f"{2 * self.indentation}search_rgn = convert_scale("):
                    lines[line_num] = 2 * self.indentation + f"\n{2 * self.indentation}".join(
                        search_bounds
                    )
                    lines[line_num] += f"\n\n{2 * self.indentation}search_rgn = convert_scale(\n"
            with open(
                os.path.join(self.model.path, "search_param.py"),
                mode="w",
                encoding="utf-8",
            ) as f:
                f.writelines(lines)
