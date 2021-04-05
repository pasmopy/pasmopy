from dyaus import PatientModelSimulations


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
    biomass_examples = BiomassExamples(
        "biomass.models",
        model_list,
        biomass_kwargs={"viz_type": "original"},
    )
    assert biomass_examples.run() is None
