# Redeconve
Deconvolution of spatial transcriptomics at single-cell resolution.


* This package is still under development. Now the basic functions are complete and we will include more functions and optimization in the future.
* To see codes analyzing the paper, please refer to https://codeocean.com/capsule/5481250/tree (not published yet).


### History
---------------------



### Installation
---------------------
Please use the following codes to install *Redeconve*:
```{r}
# install.packages("devtools")
devtools::install_github("ZxZhou4150/Redeconve", build_vignettes = F)
```
You may also set `build_vignettes` as `T`, but it may take some time to build the vignette (1 hour or so). If you don't want to build vignette, you can refer to "Redeconve manual" here, which is a simplified manual.

```{r}
browseVignettes("Redeconve")
```
