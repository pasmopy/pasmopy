import os
import shutil
from typing import Optional

from pasmopy import Model, ScipyDifferentialEvolution, run_simulation
from .models import Nakakuki_Cell_2010


model = Model(Nakakuki_Cell_2010.__package__).create()
optimizer = ScipyDifferentialEvolution(model)

def objective(x):
    '''An objective function to be minimized.'''
    return optimizer.get_obj_val(x)

def test_parameter_estimation(options: Optional[dict] = None):

    if options is None:
        options = {}
    options.setdefault("maxiter", 10)
    options.setdefault("popsize", 2)
    options.setdefault("workers", -1)

    files = [
        "best_fitness.npy",
        "count_num.npy",
        "fit_param{}.npy".format(options["maxiter"]),
        "generation.npy",
        "optimization.log",
    ]
    for i in range(1, 4):
        res = optimizer.minimize(objective, i, optimizer_options=options)
        optimizer.save_param(res, i)
        assert run_simulation(model, viz_type=str(i)) is None
        for fname in files:
            assert os.path.isfile(os.path.join(model.path, "out", str(i), fname))


def test_cleanup():
    for res in ["out", "simulation_data", "figure", "optimization_results"]:
        if os.path.isdir(os.path.join(model.path, res)):
            shutil.rmtree(os.path.join(model.path, res))