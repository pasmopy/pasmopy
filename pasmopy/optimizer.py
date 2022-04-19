import os
import sys
import warnings
from dataclasses import dataclass
from typing import Optional

import numpy as np
from biomass.exec_model import ModelObject
from scipy.optimize import differential_evolution, OptimizeResult


class _Logger(object):
    """
    Duplicate stdout to optimization.log.
    """

    def __init__(self, model_path: str, x_id: int, dirname: str):
        self.terminal = sys.stdout
        self.log = open(
            os.path.join(model_path, dirname, f"{x_id}", "optimization.log"),
            mode="w",
            encoding="utf-8",
        )

    def write(self, message: str):
        self.terminal.write(message)
        self.log.write(message)


@dataclass
class ScipyDifferentialEvolution(object):
    """
    Use ``scipy.optimize.differential_evolution`` to estimate kinetic parameters.

    Attributes
    ----------
    model : :class:`biomass.exec_model.ModelObject`
        BioMASS model object.

    Examples
    --------
    >>> from pasmopy import Model, ScipyDifferentialEvolution, run_simulation
    >>> import your_model
    >>> model = Model(your_model.__package__).create()
    >>> optimizer = ScipyDifferentialEvolution(model)
    >>> optimizer.minimize(1, optimizer_options={"workers", -1})
    >>> run_simulation(model, viz_type="1")
    """

    model: ModelObject

    def __post_init__(self):
        self.default_stdout = sys.stdout

    def _gene2objval(self, genes: np.ndarray) -> float:
        """
        Objective function: from gene vector to objective value.
        """
        param_values = self.model.problem.gene2val(genes)
        obj_val = self.model.problem.objective(param_values)

        return obj_val

    def _import_solution(self, res: OptimizeResult, x_id: int, dirname: str) -> None:
        """
        Import the solution of the optimization to the model.
        The solution vector `x` will be saved to `path_to_model`/`dirname`/`x_id`/.
        Use ``pasmopy.run_simulation`` to visualize the optimization result.

        Parameters
        ----------
        x : Union[np.ndarray, List[float]]
            The solution of the optimization.
        x_id : int (default: 0)
            Index of the parameter set.
        dirname : str
            Where to save optimization result
        """

        param_values = self.model.problem.gene2val(res.x)
        best_fitness: float = self._gene2objval(res.x)
        n_iter: int = 0
        with open(
            os.path.join(self.model.path, dirname, f"{x_id:d}", "optimization.log"),
            mode="r",
            encoding="utf-8",
        ) as f:
            log_file = f.readlines()
        for message in log_file:
            if len(message.strip()) > 0:
                n_iter += 1
        np.save(
            os.path.join(self.model.path, dirname, f"{x_id:d}", "best_fitness"),
            best_fitness,
        )
        np.save(
            os.path.join(self.model.path, dirname, f"{x_id:d}", "count_num"),
            n_iter,
        )
        np.save(
            os.path.join(self.model.path, dirname, f"{x_id:d}", "generation"),
            n_iter,
        )
        np.save(
            os.path.join(self.model.path, dirname, f"{x_id:d}", f"fit_param{n_iter:d}"),
            param_values,
        )

    def minimize(
        self,
        x_id: int,
        *,
        dirname: str = "out",
        optimizer_options: Optional[dict] = None,
    ) -> None:
        """
        Run ``scipy.optimize.differential_evolution``.
        The optimization result will be saved in ``model.path/dirname/x_id``.

        Parameters
        ----------
        x_id: int
            Index  of parameter set to estimate.
        dirname : str (default: "out")
            Where to save optimization result.
        optimizer_options : dict, optional
            Keyword arguments to pass to ``scipy.optimize.differential_evolution``.
        """
        if os.path.isdir(os.path.join(self.model.path, dirname, f"{x_id:d}")):
            raise ValueError(
                f"{dirname}{os.sep}{x_id:d} already exists in {self.model.path}. "
                "Use another parameter id."
            )
        else:
            os.makedirs(os.path.join(self.model.path, dirname, f"{x_id:d}"))

        if optimizer_options is None:
            optimizer_options = {}
        optimizer_options.setdefault("strategy", "best2bin")
        optimizer_options.setdefault("maxiter", 50)
        optimizer_options.setdefault("popsize", 3)
        optimizer_options.setdefault("tol", 1e-4)
        optimizer_options.setdefault("mutation", 0.1)
        optimizer_options.setdefault("recombination", 0.5)
        optimizer_options.setdefault("disp", True)
        optimizer_options.setdefault("polish", False)
        optimizer_options.setdefault("workers", 1)

        try:
            sys.stdout = _Logger(self.model.path, x_id, dirname)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = differential_evolution(
                    self._gene2objval,
                    [(0, 1) for _ in range(len(self.model.problem.bounds))],
                    **optimizer_options,
                )
            self._import_solution(res, x_id, dirname)
        finally:
            sys.stdout = self.default_stdout

