---
layout: post
tags: ["R","sample-variance","sample-size","variance"]
title: Dangers of sample variance at small sample size
---
Sample variance gives an unbiased estimate of the true population variance, but that doesn’t mean it’s necessarily a reliable estimate of population variance. Here, I show that sample variance itself has high variance at low sample sizes. 
<img style="float: left;" src="https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-05-07/sample_variance_vs_sample_size.png" style="width:250px; height:250px">
<!---more--->

First, we’ll create a normally-distributed parent population with a known mean, variance, and sample size. This represents a natural
population of something we’d like to study but for sake of time, money, or feasibility, we cannot measure everything. Our goal is to figure out how reliable smaller samples are with respect to estimates of variance. We’ll take increasingly larger samples from this population and see how sample variance fares.

``` r
mean = 0
SD = 20 # Therefore population variance should be ~ 400. 
# We'll set population size low-ish for sake of 
# computational time
popsize = 1000 

set.seed(100) # reproducibility
# generate the parent population
pop <- rnorm(popsize, mean, SD)

# Determine the true population variance.
# It may be different from SD^2, since we are simulating
# from a normal distribution
var(pop)
```

    ## [1] 424.8446

``` r
# Now create a sequence from 1 to popsize
# in increments of 1
Ns <- seq(1, popsize, 1)
# Within a sample size, we'll create 1000 replicates
# to help us generalize our findings
reps = 1000

# Create a function that gives an unbiased estimate of 
# population variance, given a sample. This follows from
# the definition of sample variance
var.p <- function(x) {
  var(x) * (length(x) - 1) / length(x)
}
```

# How does sample variance ‘behave’?

Using our sequence of increasing sample size (`Ns`), we’ll now create a matrix of variances. Each row number will correspond to its sample size. E.g. all values in row \[50,\] are variances from random samples of n = 50 taken from the parent population. Therefore, samples in row \[1000,\] should be identical and equal to the parent population’s variance, since we are drawing all 1000 samples from the parent population.

This process is repeated 1000 (`reps`) times for each sample size.

``` r
# This may take some time.
mymat = matrix(nrow = length(Ns), ncol = reps)
for (i in 1:dim(mymat)[1])
{
  for (j in 1:dim(mymat)[2])
  {
    mymat[i, j] = var.p(sample(pop, Ns[i]))
  }
}
rownames(mymat) <- seq(1, length(Ns))

# By definition, all the values in row [1,] will be "NA", 
# since variance cannot be computed for N = 1. 
# So we'll just remove the row.
mymat[-1, ] -> varmat 
```

It’s always good to visualize data. We’ll first plot these raw estimtates of variance.

``` r
# Sample size will be on the x-axis
# and sample variance will be on the y.
# For a given sample size, 1000 reps were performed.
plot(
  rep(2, ncol(varmat)),
  mymat[2, ],
  ylim = c(0, max(varmat)),
  xlim = c(0, nrow(mymat)),
  pch = 19,
  col = rgb(0, 0, 0, alpha = 0.2),
  xlab = 'sample size',
  ylab = 'sample variance',
  tck = 0.02,
  bty = "n"
)
for (i in 3:nrow(mymat)) {
  points(rep(i, ncol(varmat)),
         mymat[i, ],
         pch = 19,
         col = rgb(0, 0, 0, alpha = 0.2))
}

# True population variance is ~ 400.
abline(h = var(pop),
       col = rgb(0, 0, 1, alpha = 0.5),
       lwd = 3)

# Compute the mean of sample variance at each sample size
# and add it to the plot.
lines(2:popsize, rowMeans(varmat),
      col = 'orange', lwd = 3)

# Add a legend
legend(500, 1500, 
       legend=c("True population variance",
                "Means of sample variance"),
       col=c(rgb(0, 0, 1, alpha = 0.5), "orange"), 
       lty=1, lwd=3, box.lty=0)
```

![](https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-05-07/sample_variance_vs_sample_size.png)<!-- -->

We can see some pretty crazy trends. 
* **The variation in sample variance is tremendous at small sample sizes**
* **The mean of sample variance (orange line) deviates strongly at small sample sizes but otherwise generally captures the true population variance (blue).**

Let’s figure out at what point the means of sample variances seem to become unreliable. Since we know this happens at small sample sizes,
we’ll just plot cases where sample size varies from 1 to 250 to get a more refined view of the data.

``` r
# Mean variance of each sample size
rowMeans(varmat) -> varmeanz

# Plot of mean variance of each sample size
plot(
  varmeanz[1:(length(pop) / 4)],
  pch = 19,
  col = rgb(0, 0, 0, alpha = 0.2),
  xlab = 'sample size',
  ylab = 'means of sample variance',
  ylim = c(0, 500),
  tck = 0.02,
  bty = "n"
)

# Add a red line for the population variance
abline(h = var(pop),
       col = rgb(0, 0, 1, alpha = 0.5),
       lwd = 3)

# At what sample size do we get at least 95% of the 
# population variance?
abline(v = min(which(varmeanz > 0.95 * var(pop))),
       col = rgb(1, 0, 0, alpha = 0.8))
text(
  x = min(which(
    varmeanz > 0.95 * var(pop)
    )) + (0.05 * length(varmeanz)),
  y = 0.2 * var(pop),
  pos = 4,
  paste("95% of true variance at",
        min(which(
          varmeanz > 0.95 * var(pop)
        )),
        "samples")
)
```

![](https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-05-07/mean_of_sample_variance.png)<!-- -->

The vertical red line shows the minimum sample size after which 95% of true variance is achieved.

# Can we find general patterns?

At what point is sample size large enough to trust its estimation of the true variance? The answer, of course, likely depends on the parent
population’s actual variance.

Let’s create a few other examples and see if we can find common patterns. We’ll fix population means at 0, population sizes to be 1000
but vary standard deviations (and therefore variance) widely.

``` r
mean = 0
reps = 1000

# Specify our SDs and set popsize to 1000 in each case.
params <- expand.grid(SD = c(0.1, 0.5, 1, 2, 5,
                             10, 50, 100),
                      popsize = 1000)
paramslist <- split(params, seq(nrow(params)))

# This function takes all the steps we did in the 
# previous analysis and function-izes it.
varSamplr <- function(SD, popsize) {
  pop <- rnorm(popsize, mean, SD)
  Ns <- seq(1, popsize, 1)
  
  mymat = matrix(nrow = length(Ns), ncol = reps)
  for (i in 1:dim(mymat)[1])
  {
    for (j in 1:dim(mymat)[2])
    {
      mymat[i, j] = var.p(sample(pop, Ns[i]))
    }
  }
  rownames(mymat) <- seq(1, length(Ns))
  mymat[-1, ] -> varmat
  rowMeans(varmat) -> varmeanz
  
  plot(
    varmeanz[1:(length(pop) / 4)],
    pch = 19,
    col = rgb(0, 0, 0, alpha = 0.2),
    tck = 0.02,
    bty = "n",
    ylab = "means of sample variance",
    xlab = "sample size",
    main = paste(
      "N = ",
      popsize,
      "; ",
      "SD = ",
      SD,
      "\n95% of true variance at ",
      min(which(varmeanz > 0.95 * var(pop))),
      " samples",
      sep = ""
    )
  )
  abline(h = var(pop),
         col = rgb(0, 0, 1, alpha = 0.5),
         lwd = 3)
  abline(v = min(which(varmeanz > 0.95 * var(pop))),
         col = rgb(1, 0, 0, alpha = 0.8),)
}

# Organize so we can multi-plot
par(mfrow = c(4, 2))

# Run. This will take some time.
for (i in 1:nrow(params)) {
  print(i)
  varSamplr(params[i, 1], params[i, 2])
}
```

![](https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-05-07/mean_of_sample_variance_vs_pop_variance.png)<!-- -->

Pretty interesting\! 
> **Although the standard deviation varies widely across these data sets (from 0.1 to 100), taking samples of size 1 through \~ 20 severely underestimates the true population variance. So we’re seeing that samples of \< 2% of the true population size are relatively unreliable.**

# Does population size matter?

One more thing I’d like to determine is if this is a consequence of
fixing the population size at 1000. So, we’ll repeat this but instead of
varying standard deviation, we’ll vary population size.

``` r
mean = 0
reps = 1000

# This time SD is set to 1 and popsize varies
params <- expand.grid(SD = 1,
                      popsize = c(50, 100, 200, 500,
                                  1000, 2000, 5000,
                                  10000))
paramslist <- split(params, seq(nrow(params)))

# Organize so we can multi-plot
par(mfrow = c(4, 2))

# Run. This will take some time.
for (i in 1:nrow(params)) {
  print(i)
  varSamplr(params[i, 1], params[i, 2])
}
```

![](https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-05-07/sample_variance_vs_pop_size.png)<!-- -->

So it seems that no matter the population size, sample variance hits 95% of population variance after sample sizes \> 20. To put it another way, 
> **if we want to get trustworthy estimates of population variance, our sample sizes should generally be \> 20**. 

I don’t think I would have predicted that sample size should be greater than a particular number rather than as a proportion of population size. I’m sure this is covered by theory - something for me to look in to\!

