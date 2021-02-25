from dyaus import PatientModelSimulations


class BioMassExamples(PatientModelSimulations):
    pass


def test_run_with_biomass_examples():
    model_list = [
        "circadian_clock",
        "mapk_cascade",
        "nfkb_pathway",
        "tgfb_smad",
    ]
    biomass_examples = BioMassExamples(
        "biomass.models",
        model_list,
        biomass_options={"viz_type": "original"},
    )
    assert biomass_examples.run() is None
