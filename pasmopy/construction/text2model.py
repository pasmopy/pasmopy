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

try:
    from biomass import Text2Model
except ImportError:
    print("Use biomass version 0.8 or more.")