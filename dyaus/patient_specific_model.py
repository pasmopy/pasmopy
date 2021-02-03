import multiprocessing
import os
from dataclasses import dataclass
from typing import List, Optional

import biomass
from tqdm import tqdm


@dataclass
class PatientSpecificModel(object):
    """
    Run simulations and analyses of patient-specific models.

    Attributes
    ----------
    patients : list of strings
        List of patients' names or indices.

    """

    patients: List[str]

    def run(
        self,
        path_to_models: str,
        options: Optional[dict] = None,
    ) -> None:
        """
        Run simulations of patient-specific models.

        Parameters
        ----------
        path_to_models : str
            Path to the directory containing patient-specific models.

        options : dict, optional
            Arguments of biomass.run_simulation.

        """
        if options is None:
            options = {}
        options.setdefault("viz_type", "average")
        options.setdefault("show_all", False)
        options.setdefault("stdev", True)
        options.setdefault("save_format", "pdf")
        for patient in tqdm(self.patients):
            exec(
                f"from {path_to_models.lstrip(f'.{os.sep}').replace(os.sep, '.')} "
                f"import {patient}"
            )
            exec(f"model = {patient}.create()")
            exec("biomass.run_simulation(model, **options)")

    def update_out(self, path_to_out: str):
        pass
