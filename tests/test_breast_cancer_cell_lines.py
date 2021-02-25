from dyaus import PatientModelAnalyses, PatientModelSimulations

cell_lines = ["MCF7", "BT474", "SKBR3", "MDAMB231"]
for i, name in enumerate(cell_lines):
    cell_lines[i] = name + "_BREAST"


def test_run_simulations():
    pass
    """
    breast_cancer_simulations = PatientModelSimulations("models.breast", cell_lines)
    assert breast_cancer_simulations.run() is None
    """


def test_run_analyses():
    pass
    """
    breast_cancer_analyses = PatientModelAnalyses("models.breast", cell_lines)
    assert breast_cancer_analyses.run() is None
    """
