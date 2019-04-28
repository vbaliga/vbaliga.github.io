---
layout: post
tags: ["R","packages","package-installation","package-loading"]
title: Check if a package is installed [and install if not] in R
---

Say you have an R script that you share with others. You may not be sure that each user has installed all the packages the script will require. Using `install.packages()` would be unnessary for users who already have the packages and simply need to load them.

Here’s some code that provides an easy way to check whether specific packages are in the default Library. If they are, they’re simply loaded. If any packages are missing, they’re installed (with dependencies) into the default Library and are then loaded.

(This is a re-post of an entry that appeared on my old blog - see [here](https://www.vikram-baliga.com/blog/2015/7/19/a-hassle-free-way-to-verify-that-r-packages-are-installed-and-loaded)).
<!---more--->

## Install | install & load packages

``` r
## If a package is installed, it will be loaded. If any 
## are not, the missing package(s) will be installed 
## from CRAN and then loaded.

## First specify the packages of interest
packages = c("tidyverse", "geomorph",
             "phytools", "viridis")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
```


Since I already have all these packages installed, I see the following messages:

    ## Loading required package: tidyverse

    ## -- Attaching packages ------------------------------------- tidyverse 1.2.1.9000 --

    ## v ggplot2 3.1.1       v purrr   0.3.2  
    ## v tibble  2.1.1       v dplyr   0.8.0.1
    ## v tidyr   0.8.3       v stringr 1.4.0  
    ## v readr   1.3.1       v forcats 0.4.0

    ## -- Conflicts --------------------------------------------- tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    ## Loading required package: geomorph

    ## Loading required package: RRPP

    ## Loading required package: rgl

    ## Loading required package: phytools

    ## Loading required package: ape

    ## Loading required package: maps

    ## 
    ## Attaching package: 'maps'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     map

    ## Loading required package: viridis

    ## Loading required package: viridisLite

The logic of the `package.check()` function basically goes: 
* Using `lapply()` to the list of `packages`  
* If a package is not installed, install it  
* Otherwise, load it  

You can then use `search()` to determine whether all the packages have
loaded.

``` r
search()
```

    ##  [1] ".GlobalEnv"          "package:viridis"     "package:viridisLite"
    ##  [4] "package:phytools"    "package:maps"        "package:ape"        
    ##  [7] "package:geomorph"    "package:rgl"         "package:RRPP"       
    ## [10] "package:forcats"     "package:stringr"     "package:dplyr"      
    ## [13] "package:purrr"       "package:readr"       "package:tidyr"      
    ## [16] "package:tibble"      "package:ggplot2"     "package:tidyverse"  
    ## [19] "package:stats"       "package:graphics"    "package:grDevices"  
    ## [22] "package:utils"       "package:datasets"    "package:methods"    
    ## [25] "Autoloads"           "package:base"


That’s all!  
🐢
