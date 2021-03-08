# Dyaus – Dynamics-driven automatic subtyping

[![Actions Status](https://github.com/okadalabipr/dyaus/workflows/Tests/badge.svg)](https://github.com/okadalabipr/dyaus/actions)
[![License](https://img.shields.io/badge/License-Apache%202.0-brightgreen.svg)](https://opensource.org/licenses/Apache-2.0)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

**Dyaus** is a scalable framework for classifying cancer subtypes based on intracellular signaling dynamics generated from kinetic modeling.

![overview](resources/images/overview.png)

<!--
![overview](https://raw.githubusercontent.com/okadalabipr/dyaus/master/resources/images/overview.png)
-->

## Features

- Data integration
- Model construction
- Parameter estimation
- Personalized predictions of patient outcomes
- Cancer subtype classification

## Requirements

| Language      | Dependent packages                                 |
| ------------- | -------------------------------------------------- |
| Python >= 3.7 | See [requirements.txt](requirements.txt)           |
| Julia >= 1.5  | [BioMASS.jl](https://github.com/himoto/BioMASS.jl) |
| R             | [TODO] Write dependent packages here.              |

## Workflow for classifying breast cancer subtypes

- Integrate TCGA and CCLE data

  [TODO] Write analysis procedure here.

- Build an executable model of the ErbB signaling network

  ```python
  from biomass import Text2Model

  Text2Model("models/erbb_network.txt").to_biomass()
  ```

- Train model parameters against time-course datasets obtained from breast cancer cell lines

  ```shell
  $ cd training
  $ mkdir errout
  $ sh optimize_parallel.sh
  ```

- Run patient-specific models

  ```python
  from dyaus import PatientModelSimulations

  with open ("models/breast/sample_names.txt", mode="r") as f:
      TCGA_ID = f.read().splitlines()

  simulations = PatientModelSimulations("models.breast", TCGA_ID)

  simulations.run()
  ```

- Classify cancer subtypes based on the ErbB signaling dynamics

  [TODO] Write analysis procedure here.

## Installation

```
$ git clone https://github.com/okadalabipr/dyaus.git
```

## Author

- Hiroaki Imoto
  > - Building a mechanistic model of the ErbB signaling network
  > - Parameter estimation using quantitative experimental measurements
- Sawa Yamashiro
  > - Integration of transcriptomic data from TCGA and CCLE
  > - Cancer subtype classification based on inferred dynamic features

## License

[Apache-2.0 License](LICENSE)
