<br>
<p align="center">
    <a href="https://pasmopy.readthedocs.io/en/latest">
        <img src="https://raw.githubusercontent.com/pasmopy/pasmopy/master/docs/_static/img/pasmopy-project-logo.png" width="400">
    </a>
</p>

[![PyPI version](https://img.shields.io/pypi/v/pasmopy.svg?logo=PyPI&logoColor=white)](https://pypi.python.org/pypi/pasmopy)
[![Actions Status](https://github.com/pasmopy/pasmopy/workflows/Tests/badge.svg)](https://github.com/pasmopy/pasmopy/actions)
[![Documentation Status](https://img.shields.io/readthedocs/pasmopy/latest.svg?logo=read%20the%20docs&logoColor=white&&label=Docs&version=latest)](https://pasmopy.readthedocs.io/en/latest/?badge=latest)
[![License](https://img.shields.io/badge/License-Apache%202.0-green.svg?logo=apache)](https://opensource.org/licenses/Apache-2.0)
[![Downloads](https://pepy.tech/badge/pasmopy)](https://pepy.tech/project/pasmopy)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/pasmopy.svg?logo=Python&logoColor=white)](https://pypi.python.org/pypi/pasmopy)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/pasmopy/pasmopy/master.svg)](https://results.pre-commit.ci/latest/github/pasmopy/pasmopy/master)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)
[![iScience Paper](https://img.shields.io/badge/DOI-10.1016%2Fj.isci.2022.103944-blue)](https://doi.org/10.1016/j.isci.2022.103944)

**Pasmopy** is a scalable toolkit to identify prognostic factors for cancers based on intracellular signaling dynamics generated from personalized kinetic models. It is compatible with [biomass](https://github.com/biomass-dev/biomass) and offers the following features:

- Construction of mechanistic models from text
- Personalization of the model using transcriptome data
- Prediction of patient outcome based on _in silico_ signaling dynamics
- Sensitivity analysis for prediction of potential drug targets

## Documentation

Online documentation is available at https://pasmopy.readthedocs.io.

You can also build the documentation locally by running the following commands:

```shell
$ cd docs
$ make html
```

The site will live in `_build/html/index.html`.

## Installation

The latest stable release (and required dependencies) can be installed from [PyPI](https://pypi.python.org/pypi/pasmopy):

```
$ pip install pasmopy
```

Pasmopy requires Python 3.8+ to run.

## Example

### Building mathematical models of biochemical systems from text

This example shows you how to build a simple Michaelis-Menten two-step enzyme catalysis model with Pasmopy.

> E + S ⇄ ES → E + P

_An enzyme, E, binding to a substrate, S, to form a complex, ES, which in turn releases a product, P, regenerating the original enzyme._

```python
import os
from pasmopy import Text2Model, create_model, run_simulation

# Prepare a text file describing the biochemical reactions (e.g., `michaelis_menten.txt`)
reactions = """\
E + S ⇄ ES | kf=0.003, kr=0.001 | E=100, S=50
ES → E + P | kf=0.002
"""

observables = """
@obs Substrate: u[S]
@obs E_free: u[E]
@obs E_total: u[E] + u[ES]
@obs Product: u[P]
@obs Complex: u[ES]
"""

simulation_condition = """
@sim tspan: [0, 100]
"""

with open("michaelis_menten.txt", mode="w") as f:
   f.writelines(reactions)
   f.writelines(observables)
   f.writelines(simulation_condition)

# Convert the text into an executable model
description = Text2Model("michaelis_menten.txt")
description.convert()
assert os.path.isdir("michaelis_menten")

# Run simulation
model = create_model("michaelis_menten")
run_simulation(model)
```

[![michaelis_menten](https://raw.githubusercontent.com/pasmopy/pasmopy/master/docs/_static/img/michaelis_menten_sim.png)](https://pasmopy.readthedocs.io/en/latest/model_development.html#michaelis-menten-enzyme-kinetics)

For more examples, please refer to the [Documentation](https://pasmopy.readthedocs.io/en/latest/).

### Personalized signaling models for cancer patient stratification

Using Pasmopy, we built a mechanistic model of ErbB receptor signaling network, trained with protein quantification data obtained from cultured cell lines, and performed _in silico_ simulation of the pathway activities on breast cancer patients using The Cancer Genome Atlas (TCGA) transcriptome datasets. The temporal activation dynamics of Akt, extracellular signal-regulated kinase (ERK), and c-Myc in each patient were able to accurately predict the difference in prognosis and sensitivity to kinase inhibitors in triple-negative breast cancer (TNBC).

For further details, please see https://pasmopy.readthedocs.io/en/latest/personalized_model.html.

## References

- Imoto, H., Yamashiro, S. & Okada, M. A text-based computational framework for patient -specific modeling for classification of cancers. _iScience_ **25**, 103944 (2022). https://doi.org/10.1016/j.isci.2022.103944

- Imoto, H., Yamashiro, S., Murakami, K. & Okada, M. Protocol for stratification of triple-negative breast cancer patients using _in silico_ signaling dynamics. _STAR Protocols_ **3**, 101619 (2022). https://doi.org/10.1016/j.xpro.2022.101619

## Author

[Hiroaki Imoto](https://github.com/himoto)

## License

[Apache License 2.0](https://github.com/pasmopy/pasmopy/blob/master/LICENSE)
