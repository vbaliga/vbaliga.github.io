---
layout: post
tags: ["R","phytools","simmap","stochastic-character","parallel","make-simmap"]
title: Run phytools' make.simmap() in parallel
---  
<meta name="description" content="Here I provide code to run in parallel the `make.simmap()` function from phytools. It’s a Windows-friendly approach and similar to my code from another blog post, I make use of `parLapply()`.">
<p>
<img src="https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-05-19/simmap_parallel-1.png" alt="simmap parallel anoledata" style="float:right;width:200px;height:200px;margin-left:30px;">
In macroevolutionary studies, we often use stochastic character mapping to infer how a discrete trait may have evolved.
</p> 

I am very grateful that the [phytools](https://github.com/liamrevell/phytools) package allows easy implementation of character mapping via the `make.simmap()` function. Of course, this method uses a Markovian process where we sample character histories in proportion to their posterior probabilities under a given model. So we need to simulate many, many (hundreds, thousands...) of potential histories to get meaningful results. 

As with any other algorithm that we'd like to run repeatedly, it makes sense to see if parallelization can help us.

**Here I provide code to run `make.simmap()` in parallel.** It’s a Windows-friendly approach and similar to my code from [another blog post](https://vbaliga.github.io/parallel-processing-for-mcmcglmm-in-r-windows-friendly/), I make use of `parLapply()`.

<!---more--->

## Set up data

I’ll use the tree and states included in the `anoletree` dataset provided by [phytools](https://github.com/liamrevell/phytools).

``` r
# First load | install&load packages we'll need
packages = c("phytools", "parallel", "viridis")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Import the anoletree dataset from phytools
data(anoletree)
tree <- anoletree

# Use phytools::getStates() to determine the states at the tips
states <- getStates(anoletree, type="tips")

# We'll simply use an ER model
fitER <- fitMk(anoletree, states, model = "ER")

# It takes some work to re-format the Q matrix for use in
# make.simmap()
Q <- matrix(NA, length(fitER$states),
            length(fitER$states))
Q[] <- c(0, fitER$rates)[fitER$index.matrix + 1]
diag(Q) <- 0
diag(Q) <- -rowSums(Q)
colnames(Q) <- rownames(Q) <- fitER$states
```
This tree (`tree`), tip data (`states`), and transition matrix (`Q`) should give us everything we need to run through an example comparison.

## Run make.simmap() in parallel

Similar to my code from [this post](https://vbaliga.github.io/parallel-processing-for-mcmcglmm-in-r-windows-friendly/), **we will now use `parLapply()` from the [parallel](https://www.rdocumentation.org/packages/parallel/versions/3.6.0) package to execute runs of `make.simmap()` in parallel.** 

To show how much time can be saved (using these parameters and with my laptop's specs) I will time it using `Sys.time()` and then compare its timing to a 'vanilla' version run in series later on in this post. 

``` r
# Now import all these objects on to parallel clusters
# And run make.simmap() in parallel
t0 <- Sys.time()
setCores <-
  round(detectCores() * 0.85) #use 85% of available logical processors
cl <- parallel::makeCluster(getOption("cl.cores", setCores))
cl.pkg <- parallel::clusterEvalQ(cl, library(phytools))
parallel::clusterExport(cl, "tree")
parallel::clusterExport(cl, "states")
parallel::clusterExport(cl, "Q")
NRUNS <- 10   #10 'runs' of nsim = 50 -> 500 total simmaps
ER_Mkparallel <- parallel::parLapply(cl = cl, 1:NRUNS, function(i) {
  make.simmap(tree, states, Q = Q,
              nsim = 50)
})

# Super important! Make sure you stop the cluster!
# Otherwise, cores on your computer will still be dedicated 
# to running (nonexistant) processes
stopCluster(cl)

# Finally, measure the time elapsed
t1 <- Sys.time()
t1 - t0
```

    ## Time difference of 44.06479 secs

For reference, I ran this on my laptop, which is hilariously less powerful than my desktop computer. Some relevant specs: i7 processor w/ 4 cores @ ~2 GHz each; 8 GB RAM; Windows 10 Pro 64-bit; R 3.5.2 in RStudio 1.2.1335.

If you missed it above, please make sure you stop the virtual cluster using `stopCluster(cl)` at this point!

**One more important bit:** we'll now collect all products from these parallel runs into a single object which will be of class `multiSimmap`. 

``` r
# Collect all the stochastic maps
ER_Mk <- do.call(c, ER_Mkparallel)
class(ER_Mk) <- c("multiSimmap", class(ER_Mk))
# Summarize them
ER_Mk_summary <- summary(ER_Mk, plot = FALSE)
```

Then just for fun, we'll do some basic plotting of the results to showcase our simulated histories in a couple ways.

``` r
# Just for fun, I'll plot some of the results
# I'm a big fan of viridis colors
cols <- viridis(length(levels(as.factor(states))))
names(cols) <- levels(as.factor(states))

# Two example plots
# First is a summary
plot(
  ER_Mk_summary,
  colors = cols,
  type = "fan",
  fsize = 0.8,
  cex = c(0.5, 0.5)
)
```

![](https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-05-19/simmap_parallel-1.png)<!-- -->

``` r
# And now just one of the maps
plotSimmap(ER_Mk[[1]], colors = cols)
```

![](https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-05-19/simmap_parallel-2.png)<!-- -->

## So...how much time does this save us?

So it would be good to compare how fast this batch of runs completed with how it would fare in normal circumstances. Of course, since this is a stochastic process, we can't exactly generate the same batch of simmapped trees again (well, I wonder if `set.seed` does anything here -- something to check later!). 

In any case, we'll just map 500 trees based on the same data & tree using the usual process and see how long it takes. 

``` r
# Now compare with running in series
t2 <- Sys.time()
ER_Mkseries <- make.simmap(tree, states, Q = Q, nsim = 500)
```

    ## make.simmap is sampling character histories conditioned on the transition matrix
    ## 
    ## Q =
    ##             CG          GB          TC          TG          Tr          Tw
    ## CG -0.11570735  0.02314147  0.02314147  0.02314147  0.02314147  0.02314147
    ## GB  0.02314147 -0.11570735  0.02314147  0.02314147  0.02314147  0.02314147
    ## TC  0.02314147  0.02314147 -0.11570735  0.02314147  0.02314147  0.02314147
    ## TG  0.02314147  0.02314147  0.02314147 -0.11570735  0.02314147  0.02314147
    ## Tr  0.02314147  0.02314147  0.02314147  0.02314147 -0.11570735  0.02314147
    ## Tw  0.02314147  0.02314147  0.02314147  0.02314147  0.02314147 -0.11570735
    ## (specified by the user);
    ## and (mean) root node prior probabilities
    ## pi =
    ##        CG        GB        TC        TG        Tr        Tw 
    ## 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667

    ## Done.

``` r
t3 <- Sys.time()
t3 - t2
```

    ## Time difference of 1.532948 mins

``` r
# Express these time differences as a ratio
as.numeric(difftime(t3,t2,units='secs'))/as.numeric(difftime(t1,t0,units='secs'))
```

    ## [1] 2.08731

So parallelization is roughly twice as fast, at least for these specifications and on my laptop. I assume computational time varies based on the underlying data as well as the speed & number of your processors -- it would be cool to map this out as a function of tree size...etc later on!

🐢
