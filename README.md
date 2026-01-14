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
  Contains all results (PDF plots and **Table 3**) for Section 6 and Web Appendix D, produced running both scripts.


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

Default number of replications: `nsim = 1000` (can be reduced for faster execution).
Results are saved in:
```
./intermediate_results/TabC1TabC2/
```


### Tables 1, 2 and Web Appendix E.1 and E.2

- **`3_SimTab1TabE1.R`**  
  Simulations for **Table 1** and **Web Table E.1**

- **`4_SimTab2TabE2.R`**  
  Simulations for **Table 2** and **Web Table E.2**

Default number of replications: `nsim = 5000` (can be reduced for faster execution).
Results are saved, respectively, in:
```
./intermediate_results/Tab1TabE1/
./intermediate_results/Tab2TabE2/
```


### Batch Submission Scripts

The following scripts are provided for running simulations on a Linux cluster using **PBS/qsub**:

- **`1_and_2_Sims.qsub`**  
  Runs simulations for Tables C.1 and C.2 in parallel.

- **`3_SimTab1TabE1.sh`**  
  Submits multiple jobs for Table 1 and Web Table E.1.

- **`4_SimTab2TabE2.sh`**  
  Submits multiple jobs for Table 2 and Web Table E.2.

You may need to adapt the following lines depending on your system setup:

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


## Software Environment

### Local Machine (macOS)

```
R version 4.5.2 (2025-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS Sequoia 15.1.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Rome
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] xtable_1.8-4    latex2exp_0.9.8 mets_1.3.9     

loaded via a namespace (and not attached):
 [1] digest_0.6.39          numDeriv_2016.8-1.1    codetools_0.2-20       RcppArmadillo_15.2.3-1
 [5] Matrix_1.7-4           lattice_0.22-7         timereg_2.0.7          splines_4.5.2         
 [9] lava_1.8.2             parallel_4.5.2         mvtnorm_1.3-3          parallelly_1.46.1     
[13] grid_4.5.2             future_1.68.0          compiler_4.5.2         globals_0.18.0        
[17] tools_4.5.2            future.apply_1.20.1    listenv_0.10.0         Rcpp_1.1.1            
[21] survival_3.8-3        
```

### Linux HPC Cluster (Simulation Studies)
```
R version 4.5.1 (2025-06-13)
Platform: x86_64-pc-linux-gnu
Running under: Rocky Linux 9.6 (Blue Onyx)

Matrix products: default
BLAS/LAPACK: /u/ragni/spack-1.0/opt/spack/linux-zen3/openblas-0.3.29-2jjevuoakaxplu6eyjesbl2tuc4gw2kr/lib/libopenblas-r0.3.29.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/Rome
tzcode source: internal

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
[1] xtable_1.8-4           RcppArmadillo_15.2.3-1 doMC_1.3.8            
[4] iterators_1.0.14       foreach_1.5.2          mets_1.3.9            
[7] Rcpp_1.1.1             timereg_2.0.7          survival_3.8-3        

loaded via a namespace (and not attached):
 [1] digest_0.6.39       numDeriv_2016.8-1.1 codetools_0.2-20   
 [4] Matrix_1.7-4        lattice_0.22-7      splines_4.5.1      
 [7] lava_1.8.2          mvtnorm_1.3-3       parallelly_1.46.1  
[10] future_1.68.0       grid_4.5.1          compiler_4.5.1     
[13] globals_0.18.0      tools_4.5.1         future.apply_1.20.1
[16] listenv_0.10.0   
```
