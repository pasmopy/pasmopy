import os
import shutil
import time
from typing import List, Optional

from pasmopy import PatientModelAnalyses, PatientModelSimulations, Text2Model, create_model
from pasmopy.preprocessing import WeightingFactors

from .C import INCORPORATION, INDIVIDUALIZATION, REQUIREMENTS

try:
    import tests.models.breast
except ImportError:
    print("can't import 'breast' from 'tests.models'.")


PATH_TO_MODELS: str = os.path.join("tests", "models", "breast")


def path_to_patient(patient_id: str) -> str:
    return os.path.join(PATH_TO_MODELS, patient_id)


with open(path_to_patient("sample_names.txt"), mode="r") as f:
    TCGA_ID = f.read().splitlines()


with open(path_to_patient("selected_tnbc.txt"), mode="r") as f:
    TNBC_ID = f.read().splitlines()


def test_model_construction():
    # building a model
    modeldir = os.path.join("tests", "models", "erbb_network")
    if os.path.isdir(modeldir):
        shutil.rmtree(modeldir)
    Text2Model(os.path.join("tests", "models", "erbb_network.txt")).convert()

    model = create_model("tests.models.erbb_network")
    # add weighting factors
    gene_expression = {
        "ErbB1": ["EGFR"],
        "ErbB2": ["ERBB2"],
        "ErbB3": ["ERBB3"],
        "ErbB4": ["ERBB4"],
        "Grb2": ["GRB2"],
        "Shc": ["SHC1", "SHC2", "SHC3", "SHC4"],
        "RasGAP": ["RASA1", "RASA2", "RASA3"],
        "PI3K": ["PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"],
        "PTEN": ["PTEN"],
        "SOS": ["SOS1", "SOS2"],
        "Gab1": ["GAB1"],
        "RasGDP": ["HRAS", "KRAS", "NRAS"],
        "Raf": ["ARAF", "BRAF", "RAF1"],
        "MEK": ["MAP2K1", "MAP2K2"],
        "ERK": ["MAPK1", "MAPK3"],
        "Akt": ["AKT1", "AKT2"],
        "PTP1B": ["PTPN1"],
        "GSK3b": ["GSK3B"],
        "DUSP": ["DUSP5", "DUSP6", "DUSP7"],
        "cMyc": ["MYC"],
    }

    weighting_factors = WeightingFactors(model, gene_expression)
    weighting_factors.add_to_params()
    weighting_factors.set_search_bounds()
    # Edit observable.py, update normalization
    with open(
        os.path.join("tests", "models", "erbb_network", "observable.py"),
        mode="r",
        encoding="utf-8",
    ) as f:
        lines = f.readlines()
    for line_num, line in enumerate(lines):
        if line.startswith(f"{' ' * 8}self.normalization: dict = {{}}"):
            lines[line_num] = (
                f"{' ' * 8}self.normalization: dict = {{}}\n"
                + f"{' ' * 8}for obs_name in self.obs_names:\n"
                + f"{' ' * 12}self.normalization[obs_name] = "
                + f"{{'timepoint': None, 'condition': ['EGF', 'HRG']}}\n"
            )
    with open(
        os.path.join("tests", "models", "erbb_network", "observable.py"),
        mode="w",
        encoding="utf-8",
    ) as f:
        f.writelines(lines)
    # Edit search_param.py
    with open(
        os.path.join("tests", "models", "erbb_network", "search_param.py"),
        mode="r",
        encoding="utf-8",
    ) as f:
        lines = f.readlines()
    for line_num, line in enumerate(lines):
        if line.startswith("from .name2idx import C, V"):
            lines[line_num - 1] = f"{REQUIREMENTS}"
        elif line.startswith("class SearchParam(object):"):
            lines[line_num - 1] = f"\n{INDIVIDUALIZATION}\n\n"
        elif line.startswith(f"{' ' * 8}self.idx_initials = []"):
            lines[line_num] = f"{' ' * 8}self.idx_initials = [V.PIP2]\n"
        elif line.startswith(f"{' ' * 8}# parameter constraints"):
            lines[line_num - 1] = f"\n{INCORPORATION}\n"
    with open(
        os.path.join("tests", "models", "erbb_network", "search_param.py"),
        mode="w",
        encoding="utf-8",
    ) as f:
        f.writelines(lines)
    tcgamodel = os.path.join(PATH_TO_MODELS, "TCGA_3C_AALK_01A")
    if os.path.isdir(tcgamodel):
        shutil.rmtree(tcgamodel)
    shutil.move(modeldir, tcgamodel)
    if os.path.isdir(path_to_patient("TCGA_3C_AALK_01A")):
        from tests.models.breast import TCGA_3C_AALK_01A

        model = create_model(TCGA_3C_AALK_01A.__package__)
        # 220 parameters to be estimated & initial amount of PIP2.
        assert len(model.problem.idx_params) + len(model.problem.idx_initials) == 221
    else:
        raise FileNotFoundError("TCGA_3C_AALK_01A/ does not exist in tests/models/breast/.")


def test_patient_model_simulations(
    exec_model: bool = False,
    dynamical_features: Optional[List[str]] = None,
):
    # Initialization
    for patient in TCGA_ID:
        if patient in os.listdir(PATH_TO_MODELS) and patient != "TCGA_3C_AALK_01A":
            shutil.rmtree(path_to_patient(f"{patient}"))
    if os.path.isdir(os.path.join(path_to_patient("TCGA_3C_AALK_01A"), "out")):
        shutil.rmtree(os.path.join(path_to_patient("TCGA_3C_AALK_01A"), "out"))
    # Set optimized parameter sets
    breast_cancer_models: List[str] = []
    for f in os.listdir(PATH_TO_MODELS):
        if os.path.isdir(path_to_patient(f)) and (f.startswith("TCGA_") or f.endswith("_BREAST")):
            breast_cancer_models.append(f)
    for model in breast_cancer_models:
        if os.path.isdir(os.path.join(path_to_patient(f"{model}"), "out")):
            shutil.rmtree(os.path.join(path_to_patient(f"{model}"), "out"))
        shutil.copytree(
            os.path.join("tests", "models", "out"),
            os.path.join(path_to_patient(f"{model}"), "out"),
        )
    # Create patient-specific models
    for patient in TNBC_ID:
        shutil.copytree(path_to_patient("TCGA_3C_AALK_01A"), path_to_patient(f"{patient}"))
    # Execute patient-specific models
    simulations = PatientModelSimulations(
        tests.models.breast.__package__, TNBC_ID if exec_model else [TNBC_ID[0]]
    )
    start = time.time()
    assert simulations.run() is None
    elapsed = time.time() - start
    print(f"Computation time for simulating {len(TNBC_ID)} patients: {elapsed/60:.1f} [min]")
    if exec_model:
        # Extract response characteristics and visualize patient classification
        if dynamical_features is None:
            dynamical_features = ["AUC"]
        simulations.subtyping(
            "subtype_classification.pdf",
            {
                "Phosphorylated_Akt": {"EGF": dynamical_features, "HRG": dynamical_features},
                "Phosphorylated_ERK": {"EGF": dynamical_features, "HRG": dynamical_features},
                "Phosphorylated_c-Myc": {"EGF": dynamical_features, "HRG": dynamical_features},
            },
            {
                "Phosphorylated_Akt": {"timepoint": None, "condition": ["EGF", "HRG"]},
                "Phosphorylated_ERK": {"timepoint": None, "condition": ["EGF", "HRG"]},
                "Phosphorylated_c-Myc": {"timepoint": None, "condition": ["EGF", "HRG"]},
            },
            progress=False,
        )
        obs_names = ["Phosphorylated_Akt", "Phosphorylated_ERK", "Phosphorylated_c-Myc"]
        for observable in obs_names:
            assert os.path.isfile(os.path.join("classification", f"{observable}.csv"))
        assert os.path.isfile("subtype_classification.pdf")


def test_patient_model_analyses(exec_model: bool = False):
    for patient in TNBC_ID:
        assert os.path.isdir(os.path.join(PATH_TO_MODELS, patient))
    if exec_model:
        patients = TNBC_ID
        analyses = PatientModelAnalyses(
            tests.models.breast.__package__,
            patients,
            biomass_kws={
                "metric": "maximum",
                "style": "heatmap",
                "options": {"excluded_initials": ["PIP2"]},
            },
        )
        assert analyses.run() is None
        for patient in patients:
            assert os.path.isfile(
                os.path.join(
                    PATH_TO_MODELS,
                    patient,
                    "sensitivity_coefficients",
                    "initial_condition",
                    "maximum.npy",
                )
            )


def test_cleanup_models():
    # patients
    for patient in TCGA_ID:
        if patient in os.listdir(PATH_TO_MODELS):
            shutil.rmtree(path_to_patient(f"{patient}"))
    # patient classification
    if os.path.isdir("classification"):
        shutil.rmtree("classification")
    if os.path.isfile("subtype_classification.pdf"):
        os.remove("subtype_classification.pdf")
