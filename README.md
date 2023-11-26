## Code and instructions to reproduce results from "Selective inference for k-means clustering"

Chen YT and Gao LL. (2023+) [Testing for a difference in means of a single feature after clustering](https://arxiv.org/abs/). arXiv:TBD [statME].

### Step-by-step instructions for reproducing figures

N.B.: Please double check your working directories and modify the relevant lines for sourcing data and saving outputs in the code accordingly (e.g., if you are running them on a computing cluster).

1. Download this repository and install relevant packages.
```
git clone https://github.com/yiqunchen/CADET-experiments.git
cd CADET-experiments
Rscript ./install_packages.R
```
2. Generate Figure 1 (an example to illustrate the inflated Type I error rate of the naive approach, and the selective Type I error rate control of our proposal).
```
Rscript ./Gen_Figure_1.R
```
2. Generate Figure 2 (an example to illustrate the intuition for the perturbation x'(j,phi)).
```
Rscript ./Gen_Figure_2.R
```
3. Generate data used on Figures 3 and 4 (simulated data sets to assess empirical Type I error and power). This step can be computationally intensive!
``` 
Rscript ./Gen_Data_Figure_3.R
Rscript ./Gen_Data_Figure_4.R
```
4. Generate Figure 3 (demonstrate selective Type I error control ). 
```
Rscript ./Plot_Figure_3.R
```
5. Generate Figures 4, 5, and Supplementary Figure 1 (conditional power, detection probability, and the regression-spline-based estimate of power, respectively). 
```
Rscript ./Plot_Figure_4.R 
```
6. Generate Figures 6 and 7 (analysis of the single-cell data). You need to install the relevant single-cell packages in `install_packages.R`, and download the [Bone Marrow single cell data set](https://figshare.com/ndownloader/files/34701967) from the Tabular Sapiens consortium.
```
Rscript ./Gen_Figure_6.R
Rscript ./Gen_Figure_7.R
```

### Reference: 
- Chen YT and Gao LL. (2023+) Testing for a difference in means of a single feature after clustering. arXiv preprint. 




