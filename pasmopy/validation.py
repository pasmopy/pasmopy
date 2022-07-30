"""
Validate model-based predictions using drug-response data from cancer cell lines.
"""

import os
import ssl
from dataclasses import dataclass, field
from typing import Dict, List, NamedTuple, NoReturn, Optional
from urllib.request import urlopen

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import brunnermunzel


class DrugResponse(NamedTuple):
    ccle_cell_line_name: str
    primary_cell_line_name: str
    compound: str
    target: str
    doses: np.ndarray
    activity_data: np.ndarray
    activity_sd: np.ndarray
    num_data: int
    fit_type: str
    ec50: float
    ic50: float
    amax: float
    act_area: float


@dataclass
class CancerCellLineEncyclopedia(object):
    """
    Cancer Cell Line Encyclopedia (CCLE)
    https://portals.broadinstitute.org/ccle

    Attributes
    ----------
    drug_alias : dict
        Other drug names.
    _drug_response_data : ``pandas.DataFrame``
        Pharmacologic profiles for 24 anticancer drugs across 504 cell lines.
    """

    drug_alias: Dict[str, str] = field(
        default_factory=lambda: {"AZD6244": "Selumetinib", "ZD-6474": "Vandetanib"},
        init=False,
    )

    def __post_init__(self):
        url = "https://data.broadinstitute.org/ccle_legacy_data/pharmacological_profiling/CCLE_NP24.2009_Drug_data_2015.02.24.csv"
        context = ssl.create_default_context()
        context.set_ciphers("DEFAULT")
        result = urlopen(url, context=context)
        self._drug_response_data: pd.DataFrame = pd.read_csv(result)

    @property
    def drug_response_data(self) -> pd.DataFrame:
        return self._drug_response_data

    def _convert_drug_name(self, name: str) -> str:
        if name in self.drug_alias.keys():
            return self.drug_alias[name]
        else:
            return name

    def _drug2target(self, drug: str) -> str:
        target = list(
            self.drug_response_data[self.drug_response_data["Compound"] == drug]["Target"]
        )
        return target[0]

    def _extract_drug_response(
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
                    np.array([float(val) for val in df.at[i, "Doses (uM)"].split(",")]),
                    np.array(
                        [float(val) for val in df.at[i, "Activity Data (median)"].split(",")]
                    ),
                    np.array([float(val) for val in df.at[i, "Activity SD"].split(",")]),
                    df.at[i, "Num Data"],
                    df.at[i, "FitType"],
                    df.at[i, "EC50 (uM)"],
                    df.at[i, "IC50 (uM)"],
                    df.at[i, "Amax"],
                    df.at[i, "ActArea"],
                )
            )
        return drug_response_info

    def _check_args(self, drug: str) -> Optional[NoReturn]:
        if drug not in set(self.drug_response_data["Compound"]):
            raise ValueError(
                f"{drug} doesn't exist in the CCLE drug response data.\n"
                f"Should be one of {', '.join(list(set(self.drug_response_data['Compound'])))}"
            )

    @staticmethod
    def _plot_dose_response_curve(
        x: List[float],
        population: List[DrugResponse],
        color: str,
        label: str,
        show_individual: bool,
    ):
        if show_individual:
            for i, _ in enumerate(population):
                plt.plot(
                    x,
                    np.interp(x, population[i].doses, population[i].activity_data) + 100,
                    "o",
                    color=color,
                )
        y_mean = np.mean(
            [
                (np.interp(x, population[i].doses, population[i].activity_data) + 100)
                for i, _ in enumerate(population)
            ],
            axis=0,
        )
        y_err = np.std(
            [
                (np.interp(x, population[i].doses, population[i].activity_data) + 100)
                for i, _ in enumerate(population)
            ],
            axis=0,
            ddof=1,
        )
        plt.plot(
            x,
            y_mean,
            "-",
            color=color,
            label=label,
        )
        plt.fill_between(
            x,
            y_mean - y_err,
            y_mean + y_err,
            lw=0,
            color=color,
            alpha=0.1,
        )

    def _get_drug_responses(
        self, expression_ratio: pd.DataFrame, classifier: Dict[str, List[str]], drug: str
    ) -> list:
        population = [[] for _ in range(2)]
        for column, values in classifier.items():
            for i, value in enumerate(values):
                for j, level in enumerate(list(expression_ratio.loc[:, column])):
                    if level == value:
                        population[i].append(expression_ratio.index[j])
        drug_responses = [[] for _ in range(2)]
        for i in range(2):
            drug_responses[i] = self._extract_drug_response(
                population[i], [drug] * len(population[i])
            )

        return drug_responses

    def save_dose_response_curve(
        self,
        expression_ratio: pd.DataFrame,
        classifier: Dict[str, List[str]],
        drug: str,
        *,
        labels: List[str],
        config: Optional[dict] = None,
        show_individual: bool = False,
    ) -> None:
        """
        Save dose-response curves.

        Examples
        --------
        >>> from pasmopy.validation import CancerCellLineEncyclopedia
        >>> ErbB_expression_ratio = pd.read_csv(
        ...     "https://raw.githubusercontent.com/pasmopy/breast_cancer/master/drug_response/data/ErbB_expression_ratio.csv",
        ...     index_col=0,
        ... )
        >>> ccle = CancerCellLineEncyclopedia()
        >>> for drug in ["Erlotinib", "Lapatinib"]:
        ...     ccle.save_dose_response_curve(
        ...         ErbB_expression_ratio,
        ...         {"value": ["high", "low"]},
        ...         drug,
        ...         labels=["EGFR high", "EGFR low"],
        ...     )
        """
        self._check_args(drug)

        os.makedirs(os.path.join("dose_response", f"{self._drug2target(drug)}"), exist_ok=True)

        drug_responses = self._get_drug_responses(expression_ratio, classifier, drug)
        assert len(drug_responses) == 2

        if config is None:
            config = {}
        config.setdefault("figure.figsize", (7, 5))
        config.setdefault("font.size", 24)
        config.setdefault("font.family", "Arial")
        config.setdefault("mathtext.it", "Arial:italic")
        config.setdefault("axes.linewidth", 2.4)
        config.setdefault("xtick.major.width", 2.4)
        config.setdefault("ytick.major.width", 2.4)
        config.setdefault("lines.linewidth", 3)
        config.setdefault("lines.markersize", 1)
        config.setdefault("savefig.bbox", "tight")
        config.setdefault("savefig.format", "pdf")

        plt.rcParams.update(config)

        plt.gca().spines["right"].set_visible(False)
        plt.gca().spines["top"].set_visible(False)

        DOSE = [2.50e-03, 8.00e-03, 2.50e-02, 8.00e-02, 2.50e-01, 8.00e-01, 2.53e00, 8.00e00]

        for i in range(2):
            self._plot_dose_response_curve(
                DOSE,
                drug_responses[i],
                color="darkmagenta" if i == 0 else "goldenrod",
                label=labels[i],
                show_individual=show_individual,
            )
        plt.xscale("log")
        plt.xlabel("Concentration (μM)", fontsize=28)
        plt.title(f"{self._convert_drug_name(drug)}", fontsize=36)

        plt.ylabel("Relative viability (%)", fontsize=28)
        plt.yticks([0, 25, 50, 75, 100])

        plt.legend(loc="lower left", frameon=False, labelspacing=1)
        # plt.show()
        plt.savefig(
            os.path.join(
                "dose_response",
                f"{self._drug2target(drug)}",
                f"{self._convert_drug_name(drug)}",
            ),
        )
        plt.close()

    def save_activity_area(
        self,
        expression_ratio: pd.DataFrame,
        classifier: Dict[str, List[str]],
        drug: str,
        *,
        labels: List[str],
        config: Optional[dict] = None,
    ) -> None:
        """
        Save ActArea.

        Examples
        --------
        >>> from pasmopy.validation import CancerCellLineEncyclopedia
        >>> ErbB_expression_ratio = pd.read_csv(
        ...     "https://raw.githubusercontent.com/pasmopy/breast_cancer/master/drug_response/data/ErbB_expression_ratio.csv",
        ...     index_col=0,
        ... )
        >>> ccle = CancerCellLineEncyclopedia()
        >>> for drug in ["Erlotinib", "Lapatinib"]:
        ...     ccle.save_activity_area(
        ...         ErbB_expression_ratio,
        ...         {"value": ["high", "low"]},
        ...         drug,
        ...         labels=["EGFR high", "EGFR low"],
        ...     )
        """
        self._check_args(drug)

        os.makedirs(
            os.path.join(
                "activity_area",
                f"{self._drug2target(drug)}",
            ),
            exist_ok=True,
        )

        drug_responses = self._get_drug_responses(expression_ratio, classifier, drug)
        assert len(drug_responses) == 2

        activity_area = np.array(
            [
                [drug_responses[0][i].act_area for i, _ in enumerate(drug_responses[0])],
                [drug_responses[1][i].act_area for i, _ in enumerate(drug_responses[1])],
            ],
            dtype=object,
        )

        if config is None:
            config = {}
        config.setdefault("figure.figsize", (3, 5))
        config.setdefault("font.family", "Arial")
        config.setdefault("mathtext.it", "Arial:italic")
        config.setdefault("font.size", 24)
        config.setdefault("axes.linewidth", 2.4)
        config.setdefault("xtick.major.width", 2.4)
        config.setdefault("ytick.major.width", 2.4)
        config.setdefault("lines.linewidth", 1.8)
        config.setdefault("lines.markersize", 11)
        config.setdefault("savefig.bbox", "tight")
        config.setdefault("savefig.format", "pdf")

        plt.rcParams.update(config)

        plt.gca().spines["right"].set_visible(False)
        plt.gca().spines["top"].set_visible(False)

        p_value = brunnermunzel(activity_area[0], activity_area[1]).pvalue

        sns.boxplot(
            data=(activity_area[0], activity_area[1]),
            palette=sns.color_palette(["darkmagenta", "goldenrod"]),
        )
        sns.swarmplot(data=(activity_area[0], activity_area[1]), color=".48")
        plt.xticks([0, 1], [s.replace(" ", "\n") for s in labels])
        plt.ylabel("Activity area", fontsize=28)
        plt.suptitle(
            (r"$\it{p}$" + " = {:.2e}".format(p_value)) if p_value < 0.05 else "n.s.",
            fontsize=22,
        )
        # plt.show()
        plt.savefig(
            os.path.join(
                "activity_area",
                f"{self._drug2target(drug)}",
                f"{self._convert_drug_name(drug)}",
            ),
        )
        plt.close()
