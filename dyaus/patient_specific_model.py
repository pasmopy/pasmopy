import os
from dataclasses import dataclass
from importlib import import_module
from typing import List, Optional

import biomass
from tqdm import tqdm


@dataclass
class PatientSpecificModel(object):
    patients: List[str]

    def run(self, path_to_models: str, options: Optional[dict] = None):
        if options is None:
            options = {}
        options.setdefault("viz_type", "average")
        options.setdefault("show_all", False)
        options.setdefault("stdev", True)
        options.setdefault("save_format", None)
        for patient in tqdm(self.patients):
            biological_system = import_module(
                os.path.join(
                    path_to_models,
                    patient,
                )
            )
            model = biological_system.create()
            biomass.run_simulation(model, **options)

    def update_out(self, path_to_out: str):
        pass
