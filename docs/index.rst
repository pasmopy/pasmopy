=============================================
Pasmopy – Patient-Specific Modeling in Python
=============================================

.. image:: https://raw.githubusercontent.com/pasmopy/pasmopy/master/docs/_static/img/overview.png

|Actions Status| |Documentation Status| |PyPI version| |License| |Downloads| |Python versions| |Code quality| |Pre commit| |Code style|

**Pasmopy** is a scalable toolkit for classifying cancer subtypes based on intracellular signaling dynamics generated from kinetic modeling. It is compatible with `biomass <https://github.com/biomass-dev/biomass>`_ and offers the following features:

* Construction of mechanistic models from text
* Personalization of the model using transcriptome data
* Prediction of patient outcome based on *in silico* signaling dynamics
* Sensitivity analysis for prediction of potential drug targets

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   about
   installation
   model_development
   modules/index

.. |Actions Status| image:: https://github.com/pasmopy/pasmopy/workflows/Tests/badge.svg
   :target: https://github.com/pasmopy/pasmopy/actions
   :alt: Actions Status

.. |Documentation Status| image:: https://img.shields.io/readthedocs/pasmopy/latest.svg?logo=read%20the%20docs&logoColor=white&&label=Docs&version=latest
   :target: https://pasmopy.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. |PyPI version| image:: https://img.shields.io/pypi/v/pasmopy.svg?logo=PyPI&logoColor=white
   :target: https://pypi.python.org/pypi/pasmopy/
   :alt: PyPI version

.. |License| image:: https://img.shields.io/badge/License-Apache%202.0-green.svg
   :target: https://opensource.org/licenses/Apache-2.0
   :alt: License

.. |Downloads| image:: https://pepy.tech/badge/pasmopy
   :target: https://pepy.tech/project/pasmopy
   :alt: Downloads

.. |Python versions| image:: https://img.shields.io/pypi/pyversions/pasmopy.svg?logo=Python&logoColor=white
   :target: https://pypi.python.org/pypi/pasmopy/
   :alt: Python versions

.. |Code quality| image:: https://img.shields.io/lgtm/grade/python/g/pasmopy/pasmopy.svg?logo=lgtm&logoWidth=18
   :target: https://lgtm.com/projects/g/pasmopy/pasmopy/context:python
   :alt: Code quality: Python

.. |Pre commit| image:: https://results.pre-commit.ci/badge/github/pasmopy/pasmopy/master.svg
   :target: https://results.pre-commit.ci/latest/github/pasmopy/pasmopy/master
   :alt: pre-commit.ci status

.. |Code style| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black
   :alt: Code style: black