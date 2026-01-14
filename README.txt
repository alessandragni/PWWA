Supplementary information / reproducible research files for the manuscript 
Title: "Nonparametric estimation of the Patient Weighted While-Alive Estimand"

Authors: Ragni, A., Martinussen, T., Scheike, T.
In case of questions or comments please contact alessandra.ragni@polimi.it



The code was written/evaluated in R with the following software versions:

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




Simulations were run on a Linux server with software versions:

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



This folder contains the following data and files that can be used to reproduce 
all analysis and figures of the manuscript.
The working directory shold be set in this main folder.
It contains four subfolders containing the following files:

./data/:

    hhosp.rda
    An file containing the original data of the HF-ACTION study
    used for analysis in Section 6
    from O’Connor, Christopher M., et al. "Efficacy and safety of exercise training 
    in patients with chronic heart failure: HF-ACTION randomized controlled trial." 
    Jama 301.14 (2009): 1439-1450.
    Data cannot be shared due to privacy and confidentiality considerations! 
    For our paper, they were kindly made available by the BioLINCC of the National Heart, Lung, 
    and Blood Institute.
    
    colorectal.rda
    An file containing the original data of the colorectal cancer study 
    used for analysis in Web Appendix D
    from Ducreux, Michel, et al. "Sequential versus combination chemotherapy for the treatment 
    of advanced colorectal cancer (FFCD 2000–05): an open-label, randomised, phase 3 trial." 
    The lancet oncology 12.11 (2011): 1032-1044
    

./casestudies/:

    ColorectalCancer.R
    An R script that contains the code of the analysis reported in Web Appendix D of the paper.
    
    HFAction.R
    An R script that contains the code of the analysis reported in Section 6 of the paper. 
    
    ./results/
    A subfolder containing the plots saved in .pdf and the body of Table 3 saved in .txt  
    the result of evaluating both the ColorectalCancer.R and the HFAction.R scripts
    from above, allowing to reproduce all the results of Section 6 and Web Appendix D of the manuscript.
    
    
./simulations/

    1_SimTabC1.R
    An R script that performs the simulations reported Web Appendix C (Web Table C.1). 
    It makes use of the functions defined in utilsTabC1.R and functionsTabC1.cpp.
    Results of the simulation are saved into the folder ./intermediate_results/TabC1TabC2/ 
    since computation of the whole simulation can take several hours.
    For faster computation the number of replications within each simulation scenario could be reduced from 1000 
    by changing the nsim variable.
        
    2_SimTabC2.R
    An R script that performs the simulations reported Web Appendix C (Web Table C.2). 
    It makes use of the functions defined in utilsTabC2.R.
    Results of the simulation are saved into the folder ./intermediate_results/TabC1TabC2/ 
    since computation of the whole simulation can take several hours.
    For faster computation the number of replications within each simulation scenario could be reduced from 1000 
    by changing the nsim variable.
    
    1_and_2_Sims.qsub
    A PBS batch job submission script written in Bash. 
    It runs two R simulation scripts in parallel on a Linux cluster, through qsub 1_and_2_Sims.qsub
    Note that it might be necessary to change these two lines
    source $HOME/spack-1.0/share/spack/setup-env.sh
    spack load r
    
    3_SimTab1TabE1.R
    An R script that performs the simulations reported in the first part of the simulation study (Tables 1 and Web Table E.1). 
    It makes use of the functions defined in utils.R.
    Results of the simulation are saved into the folder ./intermediate_results/Tab1TabE1/ 
    since computation of the whole simulation can take several hours.
    For faster computation the number of replications within each simulation scenario could be reduced from 5000 
    by changing the nsim variable.
    
    3_SimTab1TabE1.sh
    A Bash submission script that iterates over parameter values and submits multiple PBS batch jobs using qsub. 
    Each job runs an R simulation script with different (varz, cens) parameter combinations.
    It can be run through
    chmod +x 3_SimTab1TabE1.sh
    ./3_SimTab1TabE1.sh
    Note that it might be necessary to change these two lines depending on how R was installed
    source $HOME/spack-1.0/share/spack/setup-env.sh
    spack load r
    
    4_SimTab2TabE2.R
    An R script that performs the simulations reported in the second part of the simulation study (Tables 2 and Web Table E.2). 
    It makes use of the functions defined in utils.R.
    Results of the simulation are saved into the folder ./intermediate_results/Tab2TabE2/ 
    since computation of the whole simulation can take several hours.
    For faster computation the number of replications within each simulation scenario could be reduced from 5000 
    by changing the nsim variable.
    
    4_SimTab2TabE2.sh
    A Bash submission script that iterates over parameter values and submits multiple PBS batch jobs using qsub. 
    Each job runs an R simulation script with different (type, dep, varz, scaled) parameter combinations.
    It can be run through
    chmod +x 4_SimTab2TabE2.sh
    ./4_SimTab2TabE2.sh
    Note that it might be necessary to change these two lines
    source $HOME/spack-1.0/share/spack/setup-env.sh
    spack load r
    
    
    ./intermediate_results/
    A subfolder containing the intermediate results of the simulations performed by the above four R scripts.
    Files are grouped into subdirectories based on the tables for which they are subsequently used:
        ./Tab1TabE1/  for Tables 1 and E.1
        ./Tab2TabE2/  for Tables 2 and E.2
        ./TaC1TabC2/  for Tables C.1 and C.2
    
    
    ./results/
    A subfolder containing the .txt files containing the body of the LaTex tables 
    and the R scripts for generating them:
        1_PrintTabC1.R for generating TableC1body.txt
        2_PrintTabC2.R for generating TableC2body.txt
        3_PrintTab1TabE1.R for generating Table1Latex.txt and TableE1Latex.txt
        4_PrintTab2TabE2.R for generating Table2Latex.txt and TableE2Latex.txt
        
        
        
    
