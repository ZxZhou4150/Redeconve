# Redeconve
Deconvolution of spatial transcriptomics at single-cell resolution.


* This package is still under development. Now the basic functions are complete and we will include more functions and refinement in the future.
* To see codes analyzing the paper, please refer to https://codeocean.com/capsule/1351962/tree/v1.




---------------------
### History
08/21/2023: The first official version (v1.1.0) of Redeconve is released on Github.

12/01/2023: The article was published in the journal nature communications (https://www.nature.com/articles/s41467-023-43600-9).

04/01/2024: Version v1.1.1 is released with minor bugfixes and new functions about co-localization.

06/07/2024: Version v1.1.2 is released. In this version, all non-R codes are deleted, making this package no longer needs compilation.


---------------------
### Installation
Please use the following codes to install *Redeconve*:
```{r}
# install.packages("devtools")
devtools::install_github("ZxZhou4150/Redeconve", build_vignettes = F)
```
You may also set `build_vignettes` as `T`, but it may take some time to build the vignette (tens of minutes). 


```{r}
browseVignettes("Redeconve")
```


<font color=Red>**!Note on 06/07/2024!**</font>
Setting `build_vignettes` as `T` may encounter this error: 
```
Error in graph.adjacency.dense(adjmatrix, mode = mode, weighted = weighted, :
At vendor/cigraph/src/constructors/adjacency.c:535 : Adjacency matrix should be symmetric to produce an undirected graph. Invalid value
```
This is not solved yet (actually the passed adjacency matrix IS symmetric). So please don't build the vignette at this time. See the next section for a rendered manual.


---------------------
### Manual

If you didn't build the vignette, you can refer [**here**](Redeconve manual.ipynb) for jupyter notebook and [**<u>here</u>**](https://zxzhou4150.github.io/Redeconve%20manual.html) for the rendered html page. 