# Predictor-driven Birth-Death Phylodynamics using Generalized Linear Models (Publication in preparation)

This repository contains code, and analysis XMLs accompanying the manuscript:

**Valenzuela Agüí C., Allen B. J., Stolz U., Vaughan T. G., and Stadler T.**  
*Predictor-driven Birth-Death Phylodynamics using Generalized Linear Models*

## Overview

Generalized linear models (GLMs) are used in phylodynamics to link external predictors to parameters of the phylodynamics model. In this manuscript, we provide two example applications of a flexible GLM framework for birth–death models implemented in **BEAST2**. The `GLMPrior` BEAST2 package is available at https://github.com/cecivale/GLMPrior.

This repository includes two applications from the manuscript:

-  Early spread of SARS-CoV-2 in Europe
-  Fossil sampling in non-avian dinosaurs


## Repository structure

```text
.
├── README.md
├── epiglm_sarscov2_application/
│   ├── config/
│   ├── resources/beast_xml/
│   └── workflow/
└── macroglm_dinosaur_application/
    ├── R/
    ├── predictors/
    └── XML/
```


