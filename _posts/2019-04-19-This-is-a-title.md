---
layout: post
title: this is a title
---

This is an R Markdown format used for publishing markdown documents to GitHub. When you click the **Knit** button all R code chunks are run and a markdown file (.md) suitable for publishing to GitHub is generated.
<!---more--->

Including Code
--------------

You can include R code in the document as follows:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

Including Plots
---------------

You can also embed plots, for example:
``` r
plot(pressure, pch=19, col="black", xlim=c(0,800), 
     ylim=c(0,800), tck=0.02, bty="n")
```
![pressure](/images/pressure-1.png "presssssure")

Note that the `echo = TRUE` parameter was added to the code chunk to allow printing of the R code that generated the plot.

Try out one more syntatical thing:
``` r
SEX<-as.factor(c(rep("FEMALE",44),rep("MALE",56)))
```
