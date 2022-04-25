import os
import shutil
import sys
import warnings
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
from biomass.exec_model import ModelObject, ExecModel
from scipy.optimize import differential_evolution, OptimizeResult


DIRNAME = "_tmp"

class _Logger(object):
    """
    Duplicate stdout to optimization.log.
    """

    def __init__(self, model_path: str, disp_here: bool):
        self.disp_here = disp_here
        self.terminal = sys.stdout
        self.log_file = open(
            os.path.join(model_path, "out", DIRNAME, "optimization.log"),
            mode="w",
            encoding="utf-8",
        )

    def write(self, message: str):
        if self.disp_here:
            self.terminal.write(message)
        self.log_file.write(message)


@dataclass
class ScipyDifferentialEvolution(ExecModel):
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
    >>> def objective(x):
    ...     '''An objective function to be minimized.'''
    ...     return optimizer.get_obj_val(x)
    >>> param_idx = 1
    >>> res = optimizer.minimize(objective, optimizer_options={"workers": -1})
    >>> optimizer.save_param(res, param_idx)
    >>> run_simulation(model, viz_type=str(param_idx))
    """

    model: ModelObject

    def __post_init__(self) -> None:
        self.default_stdout = sys.stdout

    def _get_n_iter(self, x_id: int) -> int:

        savedir = os.path.join(self.model.path, "out", f"{x_id:d}")
        n_iter = 0
        with open(
            os.path.join(savedir, "optimization.log"),
            mode="r",
            encoding="utf-8",
        ) as f:
            log_file = f.readlines()
        for message in log_file:
            if len(message.strip()) > 0:
                n_iter += 1
        return n_iter

    def save_param(self, res: OptimizeResult, x_id: int) -> None:
        """
        Import the solution of the optimization to the model.
        The solution vector `x` will be saved to `path_to_model`/out/`x_id`/.
        Use ``pasmopy.run_simulation`` to visualize the optimization result.

        Parameters
        ----------
        res : OptimizeResult
            The optimization result.
        x_id : int
            Index of the parameter set.
        """
        if os.path.isdir(os.path.join(self.model.path, "out", f"{x_id:d}")):
            raise ValueError(
                f"out{os.sep}{x_id:d} already exists in {self.model.path}. "
                "Use another parameter id."
            )
        else:
            os.makedirs(os.path.join(self.model.path, "out", f"{x_id:d}"))
        savedir = os.path.join(self.model.path, "out", f"{x_id:d}")
        shutil.move(os.path.join(self.model.path, "out", DIRNAME, "optimization.log"), savedir)

        param_values = self.model.problem.gene2val(res.x)
        best_fitness: float = self.get_obj_val(res.x)
        n_iter = self._get_n_iter(x_id)
        np.save(os.path.join(savedir, "best_fitness"), best_fitness)
        np.save(os.path.join(savedir, "count_num"), n_iter)
        np.save(os.path.join(savedir, "generation"), n_iter)
        np.save(os.path.join(savedir, f"fit_param{n_iter:d}"), param_values)

    def minimize(
        self,
        objective: Callable[[np.ndarray], float],
        *,
        optimizer_options: Optional[dict] = None,
        disp_here: bool = True,
    ) -> OptimizeResult:
        """
        Run ``scipy.optimize.differential_evolution``.
        The optimization result will be saved in ``model.path/_tmp/x_id``.

        Parameters
        ----------
        objective: Callable
            An objective function to be minimized. Define it using ``get_obj_val()``.
            Please see example below.
        x_id: int
            Index  of parameter set to estimate.
        optimizer_options : dict, optional
            Keyword arguments to pass to ``scipy.optimize.differential_evolution``.
        disp_here: bool (default: False)
            Whether to show the evaluated *objective* at every iteration.

        Returns
        -------
        res : OptimizeResult
            The optimization result.
        """
        os.makedirs(os.path.join(self.model.path, "out", DIRNAME), exist_ok=True)

        if optimizer_options is None:
            optimizer_options = {}
        optimizer_options.setdefault("strategy", "best1bin")
        optimizer_options.setdefault("maxiter", 50)
        optimizer_options.setdefault("popsize", 3)
        optimizer_options.setdefault("tol", 1e-4)
        optimizer_options.setdefault("mutation", 0.1)
        optimizer_options.setdefault("recombination", 0.5)
        optimizer_options.setdefault("disp", True)
        optimizer_options.setdefault("polish", False)
        optimizer_options.setdefault("workers", 1)

        if not optimizer_options["disp"]:
            raise ValueError(
                "Set optimizer_options['disp'] to True. "
                "If you don't want to see the evaluated objective function at every iteration, "
                "set the keyword argument `disp_here` to False."
            )

        try:
            sys.stdout = _Logger(self.model.path, disp_here)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = differential_evolution(
                    objective,
                    [(0, 1) for _ in range(len(self.model.problem.bounds))],
                    **optimizer_options,
                )
            return res
        finally:
            sys.stdout = self.default_stdout
