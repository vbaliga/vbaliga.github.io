---
layout: post
tags: ["R","PCA","ancestral-states","phylogenetic-PCA","variance"]
title: PCA and inferring ancestral states - some observations
---  

<meta name="description" content="How is ancestral state estimation affected by phylogenetic vs. regular PCA?">

<p>
There is a healthy debate on the circumstances under which phylogenetic principal component analysis (pPCA; Revell 2009) should be used. The technique can be valuable because it provides a rotation of multivariate data that accounts for the effects of phylogeny. 
<img src="https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-08-30/ancestral_states_vs_pca-1.png" alt="ancestral states via PCA techniques" style="float:right;width:200px;height:200px;margin-left:30px;">
But unlike 'vanilla' PCA, pPCA results in species' scores being correlated across axes. Summed eigenvalues also don’t match the total variance in the original data. Accordingly, some researchers prefer to use vanilla PCA even if phylogenetic signal in the data is strong.
</p>

I started to wonder how ancestral state estimation could be affected by choice of PCA technique. Surely, inferring ancestral states directly from scores from vanilla PC axes would lead to biased results? What are the best practices to avoid this kind of bias? Does the software we use already account for all this?

In this post, I'll show how <b>inferring ancestral states from scores on vanilla PC axes indeed leads to biased estimates. But if ancestral states are estimated beforehand and then rotated during PCA (as is done in the `geomorph` package), this problem goes away. Inferring ancestral states from phylogenetic PCA scores, however, seems to be fine.</b>

<!---more--->


## Simulate some data

We’ll first simulate some data using functions in `phytools`. We'll also use the `geiger` and `geomorph` packages later on.

``` r
# First load | install&load packages we'll need
packages = c("phytools", "geiger", "geomorph")
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

For simplicity, we’ll just simulate one tree and one data set under Brownian motion.

``` r
# Now simulate a tree with 300 tips
tree <- pbtree(n = 300)

# Now use phytools::fastBM() to simulate data under brownian motion
# We'll generate 6 traits total with varying sigma^2 values
bm_sig_small <- fastBM(tree, sig2 = 0.01, nsim = 2)
bm_sig_med <- fastBM(tree, sig2 = 0.1, nsim = 2)
bm_sig_large <- fastBM(tree, sig2 = 1, nsim = 2)

# sort them by tip name so that species are matched
tmp1 <- bm_sig_small[order(rownames(bm_sig_small)), ]
tmp2 <- bm_sig_med[order(rownames(bm_sig_med)), ]
tmp3 <- bm_sig_large[order(rownames(bm_sig_large)), ]

# put it together
data <- data.frame(tmp1, tmp2, tmp3)
colnames(data) <- paste(rep("trait", 6), 1:6, sep = "")

# now use geiger::treedata() to combine the phylogeny and data
# this ensures data are sorted and matched with the phylogeny
# making subsequent analyses a bit easier
tree_data <- treedata(tree, data, sort = TRUE)
```

## Infer ancestral states from the original data

We’ll now infer ancestral states from the original raw data. This will be done on a per-trait basis in the scale of each original trait.

Although ancestral state estimation can be tricky in and of itself, we’ll pretend that these estimates are the ‘actual’ values of ancestors. These values will be what we’d hope to find after PCA techniques have been used.

``` r
# Ancestral states from raw data via phytools::fastAnc()
real_states <-  matrix(nrow = length(tree_data$phy$tip.label) - 1,
                       ncol = ncol(tree_data$data))
for (i in 1:dim(tree_data$data)[2]) {
  real_states[, i] <- fastAnc(tree_data$phy, tree_data$data[, i])
}
colnames(real_states) <- paste(rep("a", 6), 1:6, "_real", sep = "")
head(real_states)
```

    ##          a1_real     a2_real    a3_real     a4_real   a5_real  a6_real
    ## [1,] -0.01477320 -0.04941062 -0.1980860  0.34403583 0.5868818 1.081041
    ## [2,]  0.03831289 -0.02263758 -0.4644985 -0.09940381 1.2831621 2.166376
    ## [3,]  0.13964368  0.02538198 -0.4718416 -0.28313480 2.6832781 2.256867
    ## [4,]  0.06422408 -0.03913391 -0.6715974 -0.82400049 3.2936564 2.265131
    ## [5,] -0.01869110 -0.05138655 -0.1784241  0.37676289 0.5354944 1.000941
    ## [6,] -0.01186845 -0.09610738 -0.2148796  0.43092911 0.4930362 0.705014

## Vanilla PCA & ancestral states

Time for regular PCA. We’ll first compute the PC axes. Then we’ll use the scores on these axes to estimate ancestral characters. This should lead to faulty estimates, and it will be cool to see how faulty these estimates are and whether the faultiness is proportional to underlying trait variance.

``` r
# we'll first run prcomp to get an object we can work with
reg_pca <- prcomp(tree_data$data)
# use the covariance matrix instead
reg_pca_eig <- eigen(cov(tree_data$data))
# replace the scores in the object prcomp() 
# with data %*% eigenvectors
reg_pca$x <- tree_data$data %*% reg_pca_eig$vectors

# visualize the first two axes
phylomorphospace(tree_data$phy, reg_pca$x[, 1:2])
```

![](https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-08-30/vanilla_PCA-1.png)<!-- -->

``` r
# Now ancestral states
reg_states <- matrix(nrow = length(tree_data$phy$tip.label) - 1,
                     ncol = ncol(tree_data$data))
for (i in 1:dim(tree_data$data)[2]) {
  reg_states[, i] <- fastAnc(tree_data$phy, reg_pca$x[, i])
}
colnames(reg_states) <- paste(rep("a", 6), 1:6, "_reg", sep = "")
```

An important thing to remember is that the ancestral states we just inferred are defined via PCA axes and not the original raw data. We’ll need to transform them (through eigenvectors) back to the raw data scale via linear algebra.

``` r
as.matrix(reg_states) %*% reg_pca_eig$vectors -> reg_states_transformed
```

## Phylogenetic PCA & ancestral states

And here’s phylogenetic PCA. Because we used the covariance matrix above, we’ll do so here as well. We’ll estimate phylogenetic signal beforehand to check whether `method='lamba'` should be used instead. Since the traits seem to evolve via BM (which should be the case since we simulated them under BM\!), we can opt to use `method='BM'`.

``` r
# Estimate phylogenetic signal to verify that BM should be used
physignal(tree_data$data, tree_data$phy, 
          iter = 10000, print.progress = FALSE)
```

    ## 
    ## Call:
    ## physignal(A = tree_data$data, phy = tree_data$phy, iter = 10000,  
    ##     print.progress = FALSE) 
    ## 
    ## 
    ## 
    ## Observed Phylogenetic Signal (K): 0.91927
    ## 
    ## P-value: 1e-04
    ## 
    ## Based on 10001 random permutations

``` r
# Use phytools::phyl.pca() with method="BM" and covariance mode
phy_pca <-phyl.pca(tree_data$phy,
                   as.data.frame(tree_data$data),
                   method = "BM",
                   mode = "cov")


# visualize
phylomorphospace(tree_data$phy, phy_pca$S[, 1:2])
```

![](https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-08-30/phylo_PCA-1.png)<!-- -->

``` r
# And now ancestral states from phyloPCA
phy_states_PC <-matrix(nrow = length(tree_data$phy$tip.label) - 1,
                       ncol = ncol(tree_data$data))
for (i in 1:dim(tree_data$data)[2]) {
  phy_states_PC[, i] <- fastAnc(tree_data$phy, phy_pca$S[, i])
}
colnames(phy_states_PC) <- paste(rep("a", 6), 1:6, "_phy", sep = "")
```

Again, the ancestral states we just inferred are defined via pPCA axes and not the original raw data. We’ll need to transform them (via eigenvectors) back to the raw data scale. This can be tricky if not done carefully.

We’ll basically re-perform the pPCA step-by-step and then supply the ancestral states instead of the pPCA scores.

``` r
# Most of this is taken directly from phytools::phyl.pca()
as.data.frame(tree_data$data) -> Y
Y <- as.matrix(Y)
n <- nrow(Y)
m <- ncol(Y)
C <- vcv.phylo(tree_data$phy)[rownames(Y), rownames(Y)]
temp <- phyl.vcv(Y, C, 1)
V <- temp$R
a <- t(temp$alpha)
C <- temp$C
invC <- solve(C)
es = eigen(V)
result <- list()
result$Eval <- diag(es$values)
result$Evec <- es$vectors
dimnames(result$Eval) <- list(paste("PC", 1:ncol(Y), sep = ""),
                              paste("PC", 1:ncol(Y), sep = ""))
dimnames(result$Evec) <- list(colnames(Y), paste("PC", 1:ncol(Y),
                                                 sep = ""))
A <- matrix(rep(a, n), n, m, byrow = T)
# remove the first row because number of 
# ancestral states is 1 fewer than scores
A[-1, ] -> AncA

# transform the ancestral states back to the raw data scale
(phy_states_PC %*% solve(result$Evec)) + AncA -> phy_states_transformed
```

## Compare the results

``` r
par(mfrow = c(3, 4))

plot_fun <- function(i){
    plot(real_states[,i], reg_states_transformed[,i],
         pch = 19, tck = 0.02, bty = "n",  
         col = rgb(0, 0, 0, alpha = 0.2),
         main = paste("trait ",i," ancestral states \nvia regular PCA",sep=""),
         xlab = paste("trait ",i," expected value", sep=""),
         ylab = paste("trait ",i," obtained value", set=""),
         asp = 1
         );abline(a=0,b=1)
    plot(real_states[,i], phy_states_transformed[,i],
         pch = 19, tck = 0.02, bty = "n", 
         col = rgb(0, 0, 0, alpha = 0.2),
         main = paste("trait ",i," ancestral states \nvia phylo PCA",sep=""),
         xlab = paste("trait ",i," expected value", sep=""),
         ylab = paste("trait ",i," obtained value", set=""),
         asp = 1
         );abline(a=0,b=1)
}

for (i in 1:ncol(data)){
plot_fun(i)
}
```
![](https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-08-30/ancestral_states_vs_pca-1.png)<!-- -->

[Click here](https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-08-30/ancestral_states_vs_pca-1.png) to see the full-sized image.

As a reminder, traits 1 and 2 in the first row of plots had the smallest original variance (simulated under `sig2 = 0.01`). In each successive row, sigma^2 increases: `sig2 = 0.1` and `sig2 = 1`, respectively.

Please note that in these visualizations, the orientation of the relationship between variables doesn’t matter, but rather it’s the strength of the covariance we’re interested in. This is because the direction of the first prinicpal component is arbitrary which can have a cascading effect on other PCs. It just so happens that in this example, the directions of all PC axes were similar to those of the original data, but YMMV.

In any case, it’s easy to see that ancestral states derived from scores on vanilla PC axes are faulty, whereas those from pPCA axes look just fine.

## Rotate the ancestral states too\!

Of course, it may be more reasonable to perform ancestral character estimation on the raw data and then rotate these ancestral values according to the PCA (and then throw them in the visualization). `gm.prcomp()` and the related function `plotGMphylomorphospace()` both from `geomorph` take these steps. Let’s try it out:

``` r
GM_states <- plotGMPhyloMorphoSpace(tree_data$phy,tree_data$data)
```

![](https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-08-30/plots_GM-1.png)<!-- -->

So if everything went according to plan, the ancestral states from this function should be identical to those in `real_states`.

``` r
head(GM_states)
```

    ##          trait1      trait2     trait3      trait4    trait5   trait6
    ## 301 -0.01477320 -0.04941062 -0.1980860  0.34403583 0.5868818 1.081041
    ## 302  0.03831289 -0.02263758 -0.4644985 -0.09940381 1.2831621 2.166376
    ## 303  0.13964368  0.02538198 -0.4718416 -0.28313480 2.6832781 2.256867
    ## 304  0.06422408 -0.03913391 -0.6715974 -0.82400049 3.2936564 2.265131
    ## 305 -0.01869110 -0.05138655 -0.1784241  0.37676289 0.5354944 1.000941
    ## 306 -0.01186845 -0.09610738 -0.2148796  0.43092911 0.4930362 0.705014

``` r
head(real_states)
```

    ##          a1_real     a2_real    a3_real     a4_real   a5_real  a6_real
    ## [1,] -0.01477320 -0.04941062 -0.1980860  0.34403583 0.5868818 1.081041
    ## [2,]  0.03831289 -0.02263758 -0.4644985 -0.09940381 1.2831621 2.166376
    ## [3,]  0.13964368  0.02538198 -0.4718416 -0.28313480 2.6832781 2.256867
    ## [4,]  0.06422408 -0.03913391 -0.6715974 -0.82400049 3.2936564 2.265131
    ## [5,] -0.01869110 -0.05138655 -0.1784241  0.37676289 0.5354944 1.000941
    ## [6,] -0.01186845 -0.09610738 -0.2148796  0.43092911 0.4930362 0.705014

``` r
# commenting out the line below because it's a large matrix
# but all the values will be 0 (or very close approximations)
# real_states - GM_states 
```
Yup! We get basically identical values. 


## Observations and recommendations

  - **As underlying trait variance increased, ancestral character estimation from scores on PC axes was closer to the ‘real’ values** (but was never perfect). I am uncertain if this pattern is generalizable. It may be because of how I combined these traits together – traits with higher variance may simply have larger influence on PC axes. Therefore, low-variance traits would get washed out and direct inference of their values would be increasingly nonsensical. Running further simulations could confirm if this holds up. In any case, hopefully no one is inferring ancestral states directly from scores on PC axes in published research.

  - Should regular PCA be prefered (for whatever reason), **first perform ancestral character estimation on the raw data and then rotate these ancestral values according to the PCA**. Again, `geomorph::gm.prcomp()` and `geomorph::plotGMphylomorphospace()` already implement these steps.

  - **pPCA produced the correct estimates when we inferred ancestral states from scores on pPC axes**. But, we should keep in mind that the estimates had to be back-transformed to get them to the original raw variable scale. Hopefully this is a step that people remember to take\!


🐢

