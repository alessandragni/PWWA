# Nonparametric Estimation of the Patient-Weighted While-Alive Estimand

**Supplementary information & reproducible research code**

This repository contains all the scripts 
required to reproduce the tables and figures of the manuscript:

> **“Nonparametric Estimation of the Patient-Weighted While-Alive Estimand”**  
> **Authors:** Alessandra Ragni, Torben Martinussen, Thomas Scheike

[arXiv preprint](https://arxiv.org/abs/2412.03246)

📧 **Contact:** alessandra.ragni@polimi.it

---

## Overview

The purpose of this repository is to ensure full reproducibility of the results presented in the manuscript.

It includes:
- Real-data case studies
- Extensive simulation studies
- Scripts for generating all tables and figures in the paper and web appendices

All analyses are implemented in **R**.

The repository is organized as follows:


```
.
├── data/ 
├── casestudies/ 
│ ├── HFAction.R
│ ├── ColorectalCancer.R
│ └── results/
├── simulations/ 
│ ├── intermediate_results/
│ ├── results/
│ ├── *.R
│ ├── *.sh
│ ├── *.qsub
| └── *.cpp
└── README.md
└── README.txt
```

The working directory should be set to the **main repository folder**.

---

## `./data/`

- **`hhosp.rda`**

Original data from the **HF-ACTION randomized controlled trial**, used in **Section 6** of the manuscript:

> O’Connor, C. M., et al. (2009).  
> *Efficacy and safety of exercise training in patients with chronic heart failure.*  
> **JAMA**, 301(14), 1439–1450.

This dataset **cannot be shared** due to privacy and confidentiality restrictions.  
It was kindly provided by the BioLINCC of the National Heart, Lung, and Blood Institute.


- **`colorectal.rda`**

Original data from a colorectal cancer trial, used in **Web Appendix D**:

> Ducreux, M., et al. (2011).  
> *Sequential versus combination chemotherapy for advanced colorectal cancer.*  
> **The Lancet Oncology**, 12(11), 1032–1044.

---

## `./casestudies/`

- **`HFAction.R`**  
  Reproduces the analysis presented in **Section 6** of the manuscript.

- **`ColorectalCancer.R`**  
  Reproduces the analysis presented in **Web Appendix D**.

- `./casestudies/results/`

  Contains:
  - PDF plots
  - Text file containing the body of **Table 3**

Running both scripts reproduces all results for Section 6 and Web Appendix D.

---

## 🔬 Simulation Studies

All simulation code is located in `./simulations/`.

Due to computational cost, simulations save intermediate results that are later used to generate tables.


### Web Appendix C

- **`1_SimTabC1.R`**  
  Performs simulations for **Web Table C.1**.  
  Uses functions defined in `utilsTabC1.R` and `functionsTabC1.cpp`.

- **`2_SimTabC2.R`**  
  Performs simulations for **Web Table C.2**.  
  Uses functions defined in `utilsTabC2.R`.

Results are saved in:
```
./intermediate_results/TabC1TabC2/
```

Default number of replications: `nsim = 1000` (can be reduced for faster execution).


### Tables 1, 2 and Web Appendix E.1 and E.2

- **`3_SimTab1TabE1.R`**  
  Simulations for **Table 1** and **Web Table E.1**

- **`4_SimTab2TabE2.R`**  
  Simulations for **Table 2** and **Web Table E.2**

Results are saved in:
```
./intermediate_results/Tab1TabE1/
./intermediate_results/Tab2TabE2/
```

Default number of replications: `nsim = 5000` (can be reduced for faster execution).



### Batch Submission Scripts

The following scripts are provided for running simulations on a Linux cluster using **PBS/qsub**:

- **`1_and_2_Sims.qsub`**  
  Runs simulations for Tables C.1 and C.2 in parallel.

- **`3_SimTab1TabE1.sh`**  
  Submits multiple jobs for Table 1 and Web Table E.1.

- **`4_SimTab2TabE2.sh`**  
  Submits multiple jobs for Table 2 and Web Table E.2.

⚠️ You may need to adapt the following lines depending on your system setup:

```bash
source $HOME/spack-1.0/share/spack/setup-env.sh
spack load r
```

### Generating Tables
All scripts used to generate the LaTeX tables appearing in the manuscript and web appendices are located in:
```
./simulations/results/
```

This folder contains:
- `.txt` files with the **LaTeX table bodies**
- R scripts used to generate these files from the simulation results


  - **`1_PrintTabC1.R`**  
    Generates the body of **Web Table C.1** (`TableC1body.txt`)

  - **`2_PrintTabC2.R`**  
    Generates the body of **Web Table C.2** (`TableC2body.txt`)

  - **`3_PrintTab1TabE1.R`**  
    Generates the bodies of **Table 1** and **Web Table E.1**  
    (`Table1Latex.txt`, `TableE1Latex.txt`)

  - **`4_PrintTab2TabE2.R`**  
    Generates the bodies of **Table 2** and **Web Table E.2**  
    (`Table2Latex.txt`, `TableE2Latex.txt`)


## 💻 Software Environment

### Local Machine (macOS)

### Linux HPC Cluster (Simulation Studies)

