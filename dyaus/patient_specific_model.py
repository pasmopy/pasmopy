import sys
from dataclasses import dataclass
from typing import List, Optional

import biomass


@dataclass
class PatientSpecificModel(object):
    patients: List[str]

    def run(self, options: Optional[dict] = None):
        if options is None:
            options = {}
        options.setdefault("viz_type", "average")
        options.setdefault("show_all", False)
        options.setdefault("stdev", True)
        options.setdefault("save_format", None)
        for idx, patient in enumerate(self.patients):
            exec(f"from models import {patient.strip()}")
            exec(f"model = {patient.strip()}.create()")
            biomass.run_simulation(model, **options)
            sys.stdout.write(
                "\r{} is done ({:d} / {:d}).\n".format(
                    patient.strip(),
                    idx,
                    len(self.patients),
                )
            )

    def update_out(self, path_to_out: str):
        pass
