import os
import shutil

from biomass import models

from pasmopy import PatientModelSimulations


class BiomassExamples(PatientModelSimulations):
    pass


def test_biomass_examples():
    model_list = [
        "circadian_clock",
        "insulin_signaling",
        "mapk_cascade",
        "nfkb_pathway",
        "tgfb_smad",
    ]
    biomass_examples = BiomassExamples(models.__package__, model_list)
    assert biomass_examples.run() is None


def test_cleanup():
    if os.path.isdir("biomass"):
        shutil.rmtree("biomass")
