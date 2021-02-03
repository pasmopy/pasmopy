from ..patient_specific_model import PatientSpecificModel


class BioMassExamples(PatientSpecificModel):
    pass


def test_run_with_biomass_examples():
    model_list = [
        "circadian_clock",
        "mapk_cascade",
        "nfkb_pathway",
        "tgfb_smad",
    ]
    biomass_examples = BioMassExamples(
        "biomass.models", model_list, options={"viz_type": "original"}
    )
    assert biomass_examples.run() is None
