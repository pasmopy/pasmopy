import csv
import multiprocessing
import os
import platform
from dataclasses import dataclass, field
from importlib import import_module
from typing import Any, Callable, Dict, List, Optional

import numpy as np
from biomass import ModelObject, run_analysis, run_simulation
from scipy.integrate import simps
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

    def import_model_package(self, patient: str) -> Any:
        """
        Import biomass-formatted model package.

        Parameters
        ----------
        patient : str
            Name (ID) of each patient.

        Returns
        -------
        biomass_model : module
            Patient-specific BioMASS model.
        """
        try:
            biomass_model = import_module(".".join([self.path_to_models, patient]))
            return biomass_model
        except ImportError:
            print(f"cannot import {patient.strip()} from {self.path_to_models}.")

    def parallel_execute(
        self,
        func: Callable[[str], None],
        n_proc: int,
    ) -> None:
        """
        Execute multiple models in parallel.

        Parameters
        ----------
        func : Callable
            Function executing a single patient-specific model.

        n_proc : int
            The number of worker processes to use.
        """

        if platform.system() == "Darwin":
            # fork() has always been unsafe on Mac
            # spawn* functions should be instead
            ctx = multiprocessing.get_context("spawn")
            p = ctx.Pool(processes=n_proc)
        else:
            p = multiprocessing.Pool(processes=n_proc)

        with tqdm(total=len(self.patients)) as t:
            for _ in p.imap_unordered(func, self.patients):
                t.update(1)
        p.close()


@dataclass
class PatientModelSimulations(InSilico):
    """
    Run simulations of patient-specific models.

    Attributes
    ----------
    biomass_kwargs : dict, optional
        Arguments to biomass.run_simulation.
    """

    biomass_kwargs: Optional[dict] = field(default=None)

    def _run_single_patient(self, patient: str) -> None:
        """
        Run a single patient-specifc model simulation.

        Parameters
        ----------
        patient : str
            Name (ID) of each patient.
        """

        kwargs = self.biomass_kwargs
        if kwargs is None:
            kwargs = {}
        kwargs.setdefault("viz_type", "average")
        kwargs.setdefault("show_all", False)
        kwargs.setdefault("stdev", True)
        kwargs.setdefault("save_format", "pdf")
        kwargs.setdefault("param_range", None)

        biomass_model = self.import_model_package(patient.strip())
        run_simulation(ModelObject(biomass_model.create()), **kwargs)

    def run(self, n_proc: int = multiprocessing.cpu_count() - 1) -> None:
        """
        Run simulations of multiple patient-specific models in parallel.

        Parameters
        ----------
        n_proc : int (default: multiprocessing.cpu_count() - 1)
            The number of worker processes to use.
        """

        self.parallel_execute(self._run_single_patient, n_proc)

    @staticmethod
    def _calc_response_characteristics(
        time_course_data: np.ndarray,
        metric: str,
    ) -> str:
        if metric.lower() == "max":
            response_characteristics = np.max(time_course_data)
        elif metric.lower() == "auc":
            response_characteristics = simps(time_course_data)
        elif metric.lower == "droprate":
            response_characteristics = (np.max(time_course_data) - time_course_data[-1]) / (
                len(time_course_data) - np.argmax(time_course_data)
            )
        else:
            raise ValueError("Available metrics are: 'max', 'AUC', 'droprate'.")

        return str(response_characteristics)

    def extract(
        self,
        dynamic_characteristics: Dict[str, Dict[str, List[str]]],
        normalization: bool = True,
    ) -> None:
        """
        Extract response characteristics from patient-specific signaling dynamics.

        Parameters
        ----------
        dynamic_characteristics : Dict[str, Dict[str, List[str]]]
            {"observable": {"condition": ["metric", ...], ...}, ...}.
            Characteristics in the signaling dynamics used for classification.
            'metric' must be one of 'max', 'AUC', 'droprate'.
        normalization : bool (default: True)
            Whether to perform max-normalization.

        Examples
        --------
        >>> with open ("models/breast/sample_names.txt", mode="r") as f:
                TCGA_ID = f.read().splitlines()
        >>> from dyaus import PatientModelSimulations
        >>> simulations = PatientModelSimulations("models.breast", TCGA_ID)
        >>> simulations.extract(
                {
                    "Phosphorylated_Akt": {"EGF": ["max"], "HRG": ["max"]},
                    "Phosphorylated_ERK": {"EGF": ["max"], "HRG": ["max"]},
                    "Phosphorylated_c-Myc": {"EGF": ["max"], "HRG": ["max"]},
                }
            )
        """
        os.makedirs("classification", exist_ok=True)
        for obs_name, conditions_and_metrics in dynamic_characteristics.items():
            with open(os.path.join("classification", f"{obs_name}.csv"), "w", newline="") as f:
                writer = csv.writer(f)
                header = ["Sample"]
                for condition, metrics in conditions_and_metrics.items():
                    for metric in metrics:
                        header.append(f"{condition}_{metric}")
                writer.writerow(header)

                for patient in self.patients:
                    biomass_model = self.import_model_package(patient.strip())
                    patient_specific = biomass_model.create()
                    all_data = np.load(
                        os.path.join(
                            patient_specific.path,
                            "simulation_data",
                            "simulations_all.npy",
                        )
                    )
                    data = np.array(all_data[patient_specific.obs.index(obs_name)])
                    if normalization:
                        for i in range(data.shape[0]):
                            if not np.isnan(data[i]).all():
                                data[i] /= np.nanmax(data[i])
                        data = np.nanmean(data, axis=0)
                        data /= np.max(data)
                    patient_specific_characteristics = [patient]
                    for h in header[1:]:
                        condition, metric = h.split("_")
                        patient_specific_characteristics.append(
                            self._calc_response_characteristics(
                                data[:, patient_specific.sim.conditions.index(condition)],
                                metric,
                            )
                        )
                    writer = csv.writer(f, lineterminator="\n")
                    writer.writerow(patient_specific_characteristics)


@dataclass
class PatientModelAnalyses(InSilico):
    """
    Run analyses of patient-specific models.

    Attributes
    ----------
    biomass_kwargs : dict, optional
        Arguments to biomass.run_analysis.
    """

    biomass_kwargs: Optional[dict] = field(default=None)

    def _run_single_patient(self, patient: str) -> None:
        """
        Run a single patient-specifc model analysis.

        Parameters
        ----------
        patient : str
            Name (ID) of each patient.
        """

        kwargs = self.biomass_kwargs
        if kwargs is None:
            kwargs = {}
        kwargs.setdefault("target", "initial_condition")
        kwargs.setdefault("metric", "integral")
        kwargs.setdefault("style", "heatmap")
        kwargs.setdefault("options", None)

        biomass_model = self.import_model_package(patient.strip())
        run_analysis(ModelObject(biomass_model.create()), **kwargs)

    def run(self, n_proc: int = multiprocessing.cpu_count() - 1) -> None:
        """
        Run analyses of multiple patient-specific models in parallel.

        Parameters
        ----------
        n_proc : int (default: multiprocessing.cpu_count() - 1)
            The number of worker processes to use.
        """

        self.parallel_execute(self._run_single_patient, n_proc)
