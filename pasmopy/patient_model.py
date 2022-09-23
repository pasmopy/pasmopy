import csv
import multiprocessing
import os
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Literal, Optional, Union

import numpy as np
import pandas as pd
import seaborn as sns
from biomass import create_model, run_analysis, run_simulation
from biomass.model_object import ModelObject
from scipy.integrate import simpson
from tqdm import tqdm


@dataclass
class InSilico(object):
    """
    Patient-specific in silico models.

    Attributes
    ----------
    path_to_models : str
        Path (dot-separated) to the directory containing patient-specific models.

    patients : list of strings
        List of patients' names or identifiers.
    """

    path_to_models: str
    patients: List[str]

    def __post_init__(self) -> None:
        """
        Check for duplicates in self.patients.
        """
        duplicate = [patient for patient in set(self.patients) if self.patients.count(patient) > 1]
        if duplicate:
            raise NameError(f"Duplicate patient: {', '.join(duplicate)}")

    def parallel_execute(
        self,
        func: Callable[[str], None],
        n_proc: int,
        method: Literal["spawn", "fork", "forkserver"],
        progress: bool,
    ) -> None:
        """
        Execute multiple models in parallel.

        Parameters
        ----------
        func : Callable
            Function executing a single patient-specific model.

        n_proc : int
            The number of worker processes to use.

        method : Literal["spawn", "fork", "forkserver"]
            Start method in ``multiprocessing``.

        progress : bool
            Whether the progress bar is animating or not.
        """

        ctx = multiprocessing.get_context(method)
        p = ctx.Pool(processes=n_proc)
        try:
            with tqdm(total=len(self.patients), disable=not progress) as t:
                for _ in p.imap_unordered(func, self.patients):
                    t.update(1)
        finally:
            p.close()

    @staticmethod
    def _check_ctx(context: str) -> None:
        """
        Check whether the context is appropriate.
        """
        if context not in (contexts := ["spawn", "fork", "forkserver"]):
            raise ValueError("context must be one of '{}'.".format("', '".join(contexts)))


@dataclass
class PatientModelSimulations(InSilico):
    """
    Run simulations of patient-specific models.

    Attributes
    ----------
    biomass_kws : dict, optional
        Keyword arguments to pass to ``biomass.run_simulation``.

    response_characteristics : dict[str, Callable[[1d-array], int ot float]]
        A dictionary containing functions to extract dynamic response characteristics
        from time-course simulations.
    """

    biomass_kws: Optional[dict] = field(default=None)
    response_characteristics: Dict[str, Callable[[np.ndarray], Union[int, float]]] = field(
        default_factory=lambda: dict(
            max=np.max,
            AUC=simpson,
        ),
        init=False,
    )

    def _run_single_patient(self, patient: str) -> None:
        """
        Run a single patient-specifc model simulation.
        """

        kwargs = self.biomass_kws
        if kwargs is None:
            kwargs = {}
        kwargs.setdefault("viz_type", "average")
        kwargs.setdefault("stdev", True)

        model = create_model(".".join([self.path_to_models, patient.strip()]))
        run_simulation(model, **kwargs)

    def run(
        self,
        n_proc: Optional[int] = None,
        context: Literal["spawn", "fork", "forkserver"] = "spawn",
        progress: bool = True,
    ) -> None:
        """
        Run simulations of multiple patient-specific models in parallel.

        Parameters
        ----------
        n_proc : int, optional
            The number of worker processes to use.

        context : Literal["spawn", "fork", "forkserver"] (default: "spawn")
            The context used for starting the worker processes.

        progress : bool (default: :obj:`True`)
            If :obj:`True`, the progress indicator will be shown.
        """
        if n_proc is None:
            n_proc = multiprocessing.cpu_count() - 1
        self._check_ctx(context)
        self.parallel_execute(self._run_single_patient, n_proc, context, progress)

    @staticmethod
    def _cleanup_csv(dirname: str) -> None:
        """
        Delete CSV files in folder.

        Parameters
        ----------
        dirname: str
            Name of the directory.
        """
        files = os.listdir(dirname)
        for file in files:
            if file.endswith(".csv"):
                os.remove(os.path.join(dirname, f"{file}"))

    @staticmethod
    def _normalize(
        data: np.ndarray,
        patient_specific: ModelObject,
        obs_name: str,
        normalization: dict,
    ) -> np.ndarray:
        """
        Return normalized simulation results.

        Parameters
        ----------
        data : numpy.ndarray
            Raw simulation results.
        patient_specific : biomass.model_object.ModelObject
            Patient-specific model object.
        obs_name : str
            Observable name.
        normalization : dict
            Normalization condition.

        Returns
        -------
        data : numpy.ndarray
            Normalized simulation results.
        """
        if not normalization[obs_name]["condition"]:
            normalization[obs_name]["condition"] = patient_specific.problem.conditions
        for i in range(data.shape[0]):
            if not np.isnan(data[i]).all() and not np.all(data[i] == 0.0):
                data[i] /= (
                    data[i][
                        [
                            patient_specific.problem.conditions.index(c)
                            for c in normalization[obs_name]["condition"]
                        ],
                        normalization[obs_name]["timepoint"],
                    ]
                    if normalization[obs_name]["timepoint"] is not None
                    else np.nanmax(
                        data[i][
                            [
                                patient_specific.problem.conditions.index(c)
                                for c in normalization[obs_name]["condition"]
                            ],
                        ]
                    )
                )
        data = np.nanmean(data, axis=0)
        norm_max: float = np.max(
            data[
                [
                    patient_specific.problem.conditions.index(c)
                    for c in normalization[obs_name]["condition"]
                ],
            ]
        )
        if normalization[obs_name]["timepoint"] is None and norm_max != 0.0:
            data /= norm_max
        return data

    def _extract(
        self,
        dynamical_features: Dict[str, Dict[str, List[str]]],
        normalization: dict,
        progress: bool,
    ) -> None:
        """
        Extract response characteristics from patient-specific signaling dynamics.
        """
        os.makedirs("classification", exist_ok=True)
        self._cleanup_csv("classification")
        for obs_name, conditions_and_metrics in dynamical_features.items():
            with open(
                os.path.join("classification", f"{obs_name}.csv"),
                "w",
                newline="",
            ) as f:
                writer = csv.writer(f)
                header = ["Sample"]
                for condition, metrics in conditions_and_metrics.items():
                    for metric in metrics:
                        header.append(f"{condition}_{metric}")
                writer.writerow(header)

                for patient in tqdm(self.patients, disable=not progress):
                    patient_specific = create_model(
                        ".".join([self.path_to_models, patient.strip()])
                    )
                    all_data = np.load(
                        os.path.join(
                            patient_specific.path,
                            "simulation_data",
                            "simulations_all.npy",
                        )
                    )
                    data = np.array(all_data[patient_specific.observables.index(obs_name)])
                    if obs_name in normalization.keys():
                        data = self._normalize(data, patient_specific, obs_name, normalization)
                    patient_specific_characteristics = [patient]
                    for h in header[1:]:
                        condition, metric = h.split("_")
                        patient_specific_characteristics.append(
                            str(
                                self.response_characteristics[metric](
                                    data[patient_specific.problem.conditions.index(condition)],
                                )
                            )
                        )
                    writer = csv.writer(f, lineterminator="\n")
                    writer.writerow(patient_specific_characteristics)

    def subtyping(
        self,
        fname: Optional[str],
        dynamical_features: Dict[str, Dict[str, List[str]]],
        normalization: Optional[dict] = None,
        progress: bool = True,
        *,
        clustermap_kws: Optional[dict] = None,
    ):
        """
        Classify patients based on dynamic characteristics extracted from simulation results.

        Parameters
        ----------
        fname : str, path-like or :obj:`None`
            The clustermap is saved as fname if it is not :obj:`None`.

        dynamical_features : Dict[str, Dict[str, List[str]]]
            ``{"observable": {"condition": ["metric", ...], ...}, ...}``.
            Characteristics in the signaling dynamics used for classification.

        normalization : dict, optional (default: :obj:`None`)
            * 'timepoint' : Optional[int]
                The time point at which simulated values are normalized.
                If :obj:`None`, the maximum value will be used for normalization.
            * 'condition' : list of strings
                The experimental conditions to use for normalization.
                If empty, all conditions defined in ``sim.conditions`` will be used.

        progress : bool (default: :obj:`True`)
            If :obj:`True`, the progress indicator will be shown.

        clustermap_kws : dict, optional
            Keyword arguments to pass to ``seaborn.clustermap()``.

        Examples
        --------
        Subtype classification

        >>> with open("models/breast/sample_names.txt", mode="r") as f:
        ...    TCGA_ID = f.read().splitlines()
        >>> from pasmopy import PatientModelSimulations
        >>> simulations = PatientModelSimulations("models.breast", TCGA_ID)
        >>> simulations.subtyping(
        ...    "subtype_classification.pdf",
        ...    {
        ...        "Phosphorylated_Akt": {"EGF": ["max"], "HRG": ["max"]},
        ...        "Phosphorylated_ERK": {"EGF": ["max"], "HRG": ["max"]},
        ...        "Phosphorylated_c-Myc": {"EGF": ["max"], "HRG": ["max"]},
        ...    },
        ...    {
        ...        "Phosphorylated_Akt": {"timepoint": None, "condition": ["EGF", "HRG"]},
        ...        "Phosphorylated_ERK": {"timepoint": None, "condition": ["EGF", "HRG"]},
        ...        "Phosphorylated_c-Myc": {"timepoint": None, "condition": ["EGF", "HRG"]},
        ...    },
        ...    clustermap_kws={"figsize": (9, 12)}
        ... )

        Add new characteristics

        >>> import numpy as np
        >>> def get_droprate(time_course: np.ndarray) -> float:
        ...     return - (time_course[-1] - np.max(time_course)) / (len(time_course) - np.argmax(time_course))
        >>> simulations.response_characteristics["droprate"] = get_droprate
        """
        if normalization is None:
            normalization = {}
        # seaborn clustermap
        if clustermap_kws is None:
            clustermap_kws = {}
        clustermap_kws.setdefault("z_score", 1)
        clustermap_kws.setdefault("cmap", "RdBu_r")
        clustermap_kws.setdefault("center", 0)
        # extract response characteristics
        self._extract(dynamical_features, normalization, progress)
        if fname is not None:
            characteristics: List[pd.DataFrame] = []
            files = os.listdir("classification")
            for file in files:
                observable, ext = os.path.splitext(file)
                if ext == ".csv":
                    df = pd.read_csv(os.path.join("classification", file), index_col="Sample")
                    characteristics.append(
                        df.rename(columns=lambda s: observable.replace("_", " ") + "_" + s)
                    )
            all_info = pd.concat(characteristics, axis=1)
            all_info.index.name = ""
            fig = sns.clustermap(all_info, **clustermap_kws)
            fig.savefig(fname)


@dataclass
class PatientModelAnalyses(InSilico):
    """
    Run analyses of patient-specific models.

    Attributes
    ----------
    biomass_kws : dict, optional
        Keyword arguments to pass to ``biomass.run_analysis``.
    """

    biomass_kws: Optional[dict] = field(default=None)

    def _run_single_patient(self, patient: str) -> None:
        """
        Run a single patient-specifc model analysis.
        """

        kwargs = self.biomass_kws
        if kwargs is None:
            kwargs = {}
        kwargs.setdefault("target", "initial_condition")
        kwargs.setdefault("metric", "integral")
        kwargs.setdefault("style", "heatmap")
        kwargs.setdefault("options", None)

        model = create_model(".".join([self.path_to_models, patient.strip()]))
        run_analysis(model, **kwargs)

    def run(
        self,
        n_proc: Optional[int] = None,
        context: Literal["spawn", "fork", "forkserver"] = "spawn",
        progress: bool = True,
    ) -> None:
        """
        Run analyses of multiple patient-specific models in parallel.

        Parameters
        ----------
        n_proc : int, optional
            The number of worker processes to use.

        context : Literal["spawn", "fork", "forkserver"] (default: "spawn")
            The context used for starting the worker processes.

        progress : bool (default: :obj:`True`)
            If :obj:`True`, the progress indicator will be shown.
        """
        if n_proc is None:
            n_proc = multiprocessing.cpu_count() - 1
        self._check_ctx(context)
        self.parallel_execute(self._run_single_patient, n_proc, context, progress)
