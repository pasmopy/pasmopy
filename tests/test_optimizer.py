import os
import shutil
from typing import Optional

from pasmopy import Model, OptimizationResults, ScipyDifferentialEvolution, run_simulation

from .models import Nakakuki_Cell_2010

model = Model(Nakakuki_Cell_2010.__package__).create()


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
        optimizer = ScipyDifferentialEvolution(model, i)
        def objective(x):
            """An objective function to be minimized."""
            return optimizer.get_obj_val(x)
        res = optimizer.minimize(objective, optimizer_options=options)
        optimizer.save_param(res)
        assert run_simulation(model, viz_type=str(i)) is None
        for fname in files:
            assert os.path.isfile(os.path.join(model.path, "out", str(i), fname))
    output = OptimizationResults(model)
    output.savefig(figsize=(16, 5), boxplot_kws={"orient": "v"})
    output.trace_obj()
    for file in ["optimized_params.csv", "estimated_parameter_sets.pdf", "obj_func_traces.pdf"]:
        assert os.path.isfile(os.path.join(model.path, "optimization_results", file))


def test_cleanup():
    for res in ["out", "simulation_data", "figure", "optimization_results"]:
        if os.path.isdir(os.path.join(model.path, res)):
            shutil.rmtree(os.path.join(model.path, res))
