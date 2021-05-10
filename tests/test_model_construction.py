import os
import shutil

import numpy as np
from biomass import run_simulation
from pasmopy import Text2Model


def test_preprocessing():
    for model in [
        "michaelis_menten",
        "Kholodenko_JBC_1999",
        "michaelis_menten_jl",
        "Kholodenko_JBC_1999_jl",
    ]:
        if os.path.isdir(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                model,
            )
        ):
            shutil.rmtree(
                os.path.join(
                    os.path.dirname(__file__),
                    "text_files",
                    model,
                )
            )


def test_text2model():
    for model in ["michaelis_menten", "Kholodenko_JBC_1999"]:
        if os.path.isdir(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                model,
            )
        ):
            shutil.rmtree(
                os.path.join(
                    os.path.dirname(__file__),
                    "text_files",
                    model,
                )
            )
        for lang in ["python", "julia"]:
            if model == "michaelis_menten":
                mm_kinetics = Text2Model(
                    os.path.join(
                        os.path.dirname(__file__),
                        "text_files",
                        f"{model}.txt",
                    ),
                    lang=lang,
                )
                mm_kinetics.register_word("is_dissociated", "releases")
                mm_kinetics.to_biomass_model()
            elif model == "Kholodenko_JBC_1999":
                mapk_cascade = Text2Model(
                    os.path.join(
                        os.path.dirname(__file__),
                        "text_files",
                        f"{model}.txt",
                    ),
                    lang=lang,
                )
                mapk_cascade.to_biomass_model()


def test_run_simulation():
    try:
        from .text_files import Kholodenko_JBC_1999, michaelis_menten

        for model in ["michaelis_menten", "Kholodenko_JBC_1999"]:
            run_simulation(eval(f"{model}.create()"), viz_type="original")
            simulated_values = np.load(
                os.path.join(
                    os.path.dirname(__file__),
                    "text_files",
                    model,
                    "simulation_data",
                    "simulations_original.npy",
                )
            )
            assert np.isfinite(simulated_values).all()
    except ImportError as e:
        print(e)


def test_text2markdown():
    for model in ["michaelis_menten", "Kholodenko_JBC_1999"]:
        if model == "michaelis_menten":
            mm_kinetics = Text2Model(
                os.path.join(
                    os.path.dirname(__file__),
                    "text_files",
                    f"{model}.txt",
                )
            )
            mm_kinetics.register_word("is_dissociated", "releases")
            mm_kinetics.to_markdown(n_reaction=2)
        elif model == "Kholodenko_JBC_1999":
            mapk_cascade = Text2Model(
                os.path.join(
                    os.path.dirname(__file__),
                    "text_files",
                    f"{model}.txt",
                )
            )
            mapk_cascade.to_markdown(n_reaction=25)
        assert os.path.isfile(os.path.join("markdown", model, "rate_equation.md"))
        assert os.path.isfile(os.path.join("markdown", model, "differential_equation.md"))
        shutil.rmtree("markdown")


def test_julia_models():
    necessities = [
        os.path.join("name2idx", "parameters.jl"),
        os.path.join("name2idx", "species.jl"),
        "set_model.jl",
        "observable.jl",
        "simulation.jl",
        "experimental_data.jl",
        "set_search_param.jl",
        "fitness.jl",
    ]
    for model in ["michaelis_menten", "Kholodenko_JBC_1999"]:
        for file in necessities:
            assert os.path.isfile(
                os.path.join(
                    os.path.dirname(__file__),
                    "text_files",
                    f"{model}_jl",
                    file,
                )
            )


def test_cleanup():
    for model in [
        "michaelis_menten",
        "Kholodenko_JBC_1999",
        "michaelis_menten_jl",
        "Kholodenko_JBC_1999_jl",
    ]:
        assert os.path.isdir(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                model,
            )
        )
        shutil.rmtree(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                model,
            )
        )
