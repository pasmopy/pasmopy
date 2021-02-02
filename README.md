# Dyaus

[![License](https://img.shields.io/badge/License-Apache%202.0-brightgreen.svg)](https://opensource.org/licenses/Apache-2.0)

Dyaus (**Dy**namics-driven **au**tomatic **s**ubtyping) is a scalable framework for classifying cancer subtypes based on intracellular signaling dynamics generated from kinetic modeling.

![logo](resource/images/logo.png)

## Features

- Model construction
- Data integration
- Parameter estimation
- Personalized predictions of patient outcomes
- Cancer subtype classification

## Requirements

| Language      | Dependent packages                                 |
| ------------- | -------------------------------------------------- |
| Python >= 3.7 | [biomass](https://github.com/okadalabipr/biomass)  |
| Julia >= 1.5  | [BioMASS.jl](https://github.com/himoto/BioMASS.jl) |
| R             | [TODO] Write dependent packages here.              |

## Workflow

- Build an executable model of the ErbB signaling network

  ```python
  from biomass import TextToModel

  TextToModel("erbb_network.txt").to_biomass()
  ```

- Integrate TCGA and CCLE data

  [TODO] Write analysis procedure here.

- Estimate model parameters from experimental data

  ```bash
  $ cd training
  $ sh optimize_parallel.sh
  ```

- Run patient-specific models

  ```python
  from dyaus import PatientSpecificModel

  with open ("TCGA_breast.txt", mode="r") as f:
      TCGA_ID = f.readlines()

  PatientSpecificModel(TCGA_ID).run()
  ```

- Classify cancer subtypes based on the ErbB signaling dynamics

  [TODO] Write analysis procedure here.

## Installation

```
$ git clone https://github.com/okadalabipr/dyaus.git
```

## Author

- [Hiroaki Imoto](https://github.com/himoto)
- Sawa Yamashiro

## License

[Apache-2.0 License](LICENSE)
