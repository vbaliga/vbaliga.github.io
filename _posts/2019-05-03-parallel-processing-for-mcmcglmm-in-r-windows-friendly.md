---
layout: post
tags: ["R","MCMCglmm","parallel","parallel-processing","Windows"]
title: Parallel processing for MCMCglmm in R (Windows-friendly)
---
<meta name="description" content="I set up a virtual cluster and then use the parallel::parLapply() function to run iterations of MCMCglmm() in parallel for computers running Windows.">

Lately, I have been using the [MCMCglmm](https://cran.r-project.org/web/packages/MCMCglmm/index.html) package to run linear
mixed-models in a Bayesian framework. The documentation is generally very good but there seems to be relatively little support for using parallel processing (here: using multiple cores on your machine) to efficiently run large volumes of mcmc runs. This is especially true for Windows users, who cannot use functions like `parallel::mclapply()`.

I’m happy to share that I have worked out a solution using the [parallel](https://www.rdocumentation.org/packages/parallel/versions/3.5.1) package. Basically, I set up a virtual cluster and then use the `parallel::parLapply()` function to run iterations of `MCMCglmm()` in parallel.

<small>(This is a re-post of an entry that appeared on my old blog - see **[here](https://www.vikram-baliga.com/blog/2018/9/30/parallel-processing-for-mcmcglmm-in-r-windows-friendly)**).</small>
<!---more--->

## Data

I’ll use “Example 2” from the [MCMCglmm() function help](https://www.rdocumentation.org/packages/MCMCglmm/versions/2.26/topics/MCMCglmm). You can skip ahead to the next section if instead you’d like to tailor this to your own data & analysis.

First load (or install\&load) the `MCMCglmm` and `parallel` packages:

``` r
packages = c("MCMCglmm", "parallel")

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

    ## Loading required package: MCMCglmm

    ## Loading required package: Matrix

    ## Loading required package: coda

    ## Loading required package: ape

    ## Loading required package: parallel

With the packages loaded, we’ll prep our data set. Lifting this directly from the `MCMCglmm()` help page:

``` r
data(bird.families)
phylo.effect <- rbv(bird.families, 1, nodes = "TIPS")
phenotype <- phylo.effect + rnorm(dim(phylo.effect)[1], 0, 1)

# simulate phylogenetic and residual effects
# with unit variance
test.data <- data.frame(phenotype = phenotype,
                        taxon = row.names(phenotype))
Ainv <- inverseA(bird.families)$Ainv

# inverse matrix of shared phyloegnetic history
prior <- list(R = list(V = 1, nu = 0.002), 
              G = list(G1 = list(V = 1, nu = 0.002)))

model2 <- MCMCglmm(
  phenotype ~ 1,
  random =  ~ taxon,
  ginverse = list(taxon = Ainv),
  data = test.data,
  prior = prior,
  verbose = FALSE,
  nitt = 1300,
  burnin = 300,
  thin = 1
)
summary(model2)
```

    ## 
    ##  Iterations = 301:1300
    ##  Thinning interval  = 1
    ##  Sample size  = 1000 
    ## 
    ##  DIC: 468.589 
    ## 
    ##  G-structure:  ~taxon
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## taxon    0.4607  0.04937    1.424    4.904
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units     1.568   0.8002    2.134    19.49
    ## 
    ##  Location effects: phenotype ~ 1 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp pMCMC
    ## (Intercept)   -0.1574  -0.6576   0.3099      812 0.488

Of course, the example provided sets `nitt` to only 1300, yielding an ESS of only ~800 for the fixed effect. I am guessing this is intended to make sure the example is quick to run.

Boosting this to `nitt = 100000`, `burnin = 10000`, and `thin = 10` gives a more healthy ESS > 8000. But please note that this will take a lot longer to finish (I’ll leave it up to you to use the `Sys.time()` function to time it yourself).

## Run MCMC chains in parallel

Whenever conducting MCMC-based analyses, it’s advisable to conduct multiple runs (different chains) and then assess convergence. I’ll leave the convergence assessments for another day (but here’s [a good StackExchange post](https://stats.stackexchange.com/questions/507/what-is-the-best-method-for-checking-convergence-in-mcmc)). For now we’ll just conduct 10 runs of this model, each using `nitt = 100000`, using parallel processing. 

***PLEASE NOTE**: I am setting this up to use only 80% of your machine’s total logical processors. You can certainly harness all of your CPUs if you’d like, although I advise against doing so if any of your MCMC runs take more than a few minutes. It also doesn’t make sense to set the number of logical processors to be greater than the number of runs (chains), but more on that later. Anyway, treat your silicon well!*

``` r
# use detectCores() by itself if you want all CPUs
setCores <- round(detectCores() * 0.8)

# make the cluster
cl <- makeCluster(getOption("cl.cores", setCores))
  # EDIT ON 2020-07-27: I have been informed that Mac users 
  # may have better luck using:
  # cl <- parallel::makeCluster(getOption("cl.cores", setCores), 
  #                             setup_strategy = "sequential")
  # This is due to an apparent issue in RStudio. 
  # See this stackoverflow page for details:
  # https://stackoverflow.com/questions/61700586/r-makecluster-command-used-to-work-but-now-fails-in-rstudio

# load the MCMCglmm package within the cluster
cl.pkg <- clusterEvalQ(cl, library(MCMCglmm))

# import each object that's necessary to run the function
clusterExport(cl, "prior")
clusterExport(cl, "test.data")
clusterExport(cl, "Ainv")

# use parLapply() to execute 10 runs of MCMCglmm()
# each with nitt=100000
model2_10runs <- parLapply(cl = cl, 1:10, function(i) {
  MCMCglmm(
    phenotype ~ 1,
    random =  ~ taxon,
    ginverse = list(taxon = Ainv),
    data = test.data,
    prior = prior,
    verbose = FALSE,
    nitt = 100000,
    burnin = 10000,
    thin = 10
  )
})

# once it's finished, use stopCluster() to stop running
# the parallel cluster
stopCluster(cl)
```

The `model2_10runs` object is a list that contains each of the 10 mcmc models. You can perform all the usual summarization, plotting…etc, but just be sure to specify models within the list, e.g.:
`summary(model2_10runs[[3]])` to summarize the third model out of the 10

As I mentioned above, we’ll leave convergence and other fun topics like autocorrelation for another day.

That’s all!

🐢
