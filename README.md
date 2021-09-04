# Pasmopy – Patient-Specific Modeling in Python

![overview](https://raw.githubusercontent.com/pasmopy/pasmopy/master/docs/_static/img/overview.png)

[![Actions Status](https://github.com/pasmopy/pasmopy/workflows/Tests/badge.svg)](https://github.com/pasmopy/pasmopy/actions)
[![License](https://img.shields.io/badge/License-Apache%202.0-green.svg?logo=apache)](https://opensource.org/licenses/Apache-2.0)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

**Pasmopy** is a scalable toolkit to identify prognostic factors for cancers based on intracellular signaling dynamics generated from personalized kinetic models. It is compatible with [biomass](https://github.com/biomass-dev/biomass) and offers the following features:

- Construction of mechanistic models from text
- Personalization of the model using transcriptome data
- Prediction of patient outcome through classification based on signaling dynamics
- Sensitivity analysis for prediction of potential drug targets

## Installation

The latest stable release (and required dependencies) can be installed from PyPI:

```bash
$ pip install pasmopy
```

Pasmopy requires Python 3.7+ to run.

## Example

### Building mathematical models of biochemical systems from text

This example shows you how to build a simple Michaelis-Menten two-step enzyme catalysis model with Pasmopy.

> E + S ⇄ ES → E + P

_An enzyme, E, binding to a substrate, S, to form a complex, ES, which in turn releases a product, P, regenerating the original enzyme._

1. Prepare a text file describing biochemical reactions (`michaelis_menten.txt`)

   ```
   E binds S --> ES | kf=0.003, kr=0.001 | E=100, S=50
   ES dissociates to E and P | kf=0.002, kr=0

   @obs Substrate: u[S]
   @obs E_free: u[E]
   @obs E_total: u[E] + u[ES]
   @obs Product: u[P]
   @obs Complex: u[ES]

   @sim tspan: [0, 100]
   ```

1. Convert text into an executable model

   ```python
   from pasmopy import Text2Model

   Text2Model("michaelis_menten.txt").convert()
   ```

1. Run simulation with [`biomass`](https://github.com/biomass-dev/biomass)

   ```python
   from biomass import Model, run_simulation
   import michaelis_menten

   model = Model(michaelis_menten.__package__).create()
   run_simulation(model)
   ```

   ![michaelis_menten](https://raw.githubusercontent.com/pasmopy/pasmopy/master/docs/_static/img/michaelis_menten_sim.png)

For more examples, please refer to the Documentation.

## Author

[Hiroaki Imoto](https://github.com/himoto)

## License

[Apache License 2.0](https://github.com/pasmopy/pasmopy/blob/master/LICENSE)
