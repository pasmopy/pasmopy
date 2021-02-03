import multiprocessing
import os
import platform
from dataclasses import dataclass, field
from typing import List, Optional

import biomass
from tqdm import tqdm


@dataclass
class PatientSpecificModel(object):
    """
    Run simulations and analyses of patient-specific models.

    Attributes
    ----------
    path_to_models : str
        Path to the directory containing patient-specific models.

    patients : list of strings
        List of patients' names or indices.

    options : dict, optional
        Arguments of biomass.run_simulation.

    """

    path_to_models: str
    patients: List[str]
    options: Optional[dict] = field(default=None)

    def _run_single_patient(self, index: int) -> None:
        """
        Parameters
        ----------
        idx : int
            Index of each patient.

        """
        options = self.options
        if options is None:
            options = {}
        options.setdefault("viz_type", "average")
        options.setdefault("show_all", False)
        options.setdefault("stdev", True)
        options.setdefault("save_format", "pdf")
        options.setdefault("param_range", None)

        try:
            exec(
                f"from {self.path_to_models.lstrip(f'.{os.sep}').replace(os.sep, '.')} "
                f"import {self.patients[index].strip()}",
            )
        except ImportError as e:
            print(f"Cannot import {self.patients[index].strip()} from {self.path_to_models}", e)

        model = eval(f"{self.patients[index].strip()}.create()")
        biomass.run_simulation(model, **options)
        print(f"[{self.patients[index].strip()}] finished.")

    def run(self, n_proc: int = multiprocessing.cpu_count() - 1) -> None:
        """
        Run simulations of multiple patient-specific models in parallel.

        Parameters
        ----------
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
            for _ in p.imap_unordered(
                self._run_single_patient,
                range(len(self.patients)),
            ):
                t.update(1)
        p.close()
