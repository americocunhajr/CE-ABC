## Cross-Entropy Approximate Bayesian Computation

**CE-ABC: Cross-Entropy Approximate Bayesian Computation** is a Matlab package that implements a framework for uncertainty quantification in mechanistic epidemic models defined by ordinary differential equations (ODEs). This package combines the cross-entropy method for optimization and approximate Bayesian computation for statistical inference. With straightforward adaptations, the CE-ABC strategy can be applied to various other systems, including mechanical, electrical, and coupled systems.

### Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Usage](#usage)
- [Documentation](#documentation)
- [Reproducibility](#reproducibility)
- [Authors](#authors)
- [Citing CE-ABC](#citing-ce-abc)
- [License](#license)
- [Institutional support](#institutional-support)
- [Funding](#funding)

### Overview
**CE-ABC** addresses model calibration and uncertainty quantification in mechanistic models, primarily for epidemic modeling. The package integrates the cross-entropy method, which is a powerful optimization technique, with approximate Bayesian computation, a statistical inference method. This combination allows for efficient and accurate calibration and uncertainty quantification in ODE-based models.

For more details, refer to the following paper:
- **A. Cunha Jr, D. A. W. Barton, and T. G. Ritto**, *Uncertainty quantification in mechanistic epidemic models via cross-entropy approximate Bayesian computation*, Nonlinear Dynamics, vol. 111, pp. 9649–9679, 2023. <a href="https://doi.org/10.1007/s11071-023-08327-8" target="_blank">DOI</a> 

Preprint available <a href="https://arxiv.org/abs/2207.12111" target="_blank">here</a>.

### Features
- Combines cross-entropy method for optimization with approximate Bayesian computation for statistical inference
- Applicable to mechanistic models defined by ODEs
- Flexible framework for various systems (mechanical, electrical, coupled, etc.)
- Numerically robust and efficient implementation
- Educational style for intuitive use
- Includes example scripts for representative benchmark tests

### Usage
To get started with **CE-ABC**, follow these steps:
1. Clone the repository:
   ```bash
   git clone https://github.com/americocunhajr/CE-ABC.git
   ```
2. Navigate to the code directory:
   ```bash
   cd CE-ABC/CE-ABC-1.0
   ```
3. For a deterministic simulation with SEIRpAHD model, execute:
   ```bash
   Main_IVP_SEIRpAHD
   ```
4. For a stochastic simulation with SEIRpAHD model, execute:
   ```bash
   Main_CE_ABC_SEIRpAHD
   ```
5. For a stochastic simulation with SEIRpAHDbeta model, execute:
   ```bash
   Main_CE_ABC_SEIRpAHDbeta
   ```
5. To plot Rio de Janeiro COVID-19 data, execute:
   ```bash
   Main_COVID19RJ_Data_plot
   ```

### Documentation
**CE-ABC** routines are well-commented to explain their functionality. Each routine includes a description of its purpose and a list of inputs and outputs. Examples with representative benchmark tests are provided to illustrate the code's functionality.

### Reproducibility
Simulations done with **CE-ABC** are fully reproducible, as can be seen on this <a href="https://codeocean.com/capsule/5200426/tree/v2" target="_blank">CodeOcean capsule</a>.

### Authors
- Americo Cunha Jr
- David A. W. Barton
- Thiago G. Ritto

### Citing CE-ABC
If you use **CE-ABC** in your research, please cite the following publication:
- *A. Cunha Jr, D. A. W. Barton, and T. G. Ritto, Uncertainty quantification in mechanistic epidemic models via cross-entropy approximate Bayesian computation, Nonlinear Dynamics, vol. 111, pp. 9649–9679, 2023 https://doi.org/10.1007/s11071-023-08327-8*

```
@article{CunhaJr2023p9649,
   author  = {A {Cunha~Jr} and D. A. W. Barton and T. G. Ritto},
   title   = {Uncertainty quantification in mechanistic epidemic models via cross-entropy approximate Bayesian computation},
   journal = {Nonlinear Dynamics},
   year    = {2023},
   volume  = {111},
   pages   = {9649–9679},
   doi    = {10.1007/s11071-023-08327-8},
}
```

### License
**CE-ABC** is released under the MIT license. See the LICENSE file for details. All new contributions must be made under the MIT license.

<img src="logo/mit_license_red.png" width="10%"> 

### Institutional support

<img src="logo/logo_uerj_color.jpeg" width="10%"> &nbsp; &nbsp; <img src="logo/logo_bristol.png" width="32%"> &nbsp; &nbsp; <img src="logo/logo_ufrj.png" width="8%">

### Funding

<img src="logo/cnpq.png" width="20%"> &nbsp; &nbsp; <img src="logo/capes.png" width="10%">  &nbsp; &nbsp; &nbsp; <img src="logo/faperj.jpg" width="20%"> &nbsp; &nbsp; &nbsp; <img src="logo/epsrc.png" width="18%">
