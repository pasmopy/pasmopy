"""
:py:class:`Text2Model` is a useful class to build an ordinary differential equation (ODE) model from a text file describing biochemical systems.
It was originally developed in the pasmopy project [1]_ but migrated to the `biomass <https://github.com/biomass-dev/biomass>`_ framework.

* **How to use:** https://pasmopy.readthedocs.io/en/latest/model_development.html
* **API reference:** https://biomass-core.readthedocs.io/en/latest/api/construction.html

References
----------
.. [1] Imoto, H., Yamashiro, S. & Okada, M.
    A text-based computational framework for patient -specific modeling for classification of cancers.
    iScience 25, 103944 (2022).
"""
from typing import Final

from biomass import Text2Model, create_model, run_simulation

MICHAELIS_MENTEN: Final = """\n
E + S ⇄ ES | kf=0.003, kr=0.001 | E=100, S=50
ES → E + P | kf=0.002

@obs Substrate: u[S]
@obs E_free: u[E]
@obs E_total: u[E] + u[ES]
@obs Product: u[P]
@obs Complex: u[ES]

@sim tspan: [0, 100]
"""


def run_example_model() -> None:
    """
    A simple Michaelis-Menten two-step enzyme catalysis model.
    See README at https://github.com/pasmopy/pasmopy.
    """
    with open("michaelis_menten.txt", mode="w") as f:
        f.writelines(MICHAELIS_MENTEN)
    description = Text2Model("michaelis_menten.txt")
    description.convert()
    model = create_model("michaelis_menten")
    run_simulation(model)
