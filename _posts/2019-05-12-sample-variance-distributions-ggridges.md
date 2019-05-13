---
layout: post
tags: ["R","ggridges","sample-variance","sample-size","variance"]
title: Sample variance at small sample sizes II - distributions
---  

<p>
<img src="https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-05-12/sample_variance_distributions_ggridges.png" alt="sample variance ggridges" style="float:right;width:200px;height:200px;margin-left:30px;">
In my <b><a href="https://vbaliga.github.io/dangers-sample-variance-small-sample-size/">previous post</a></b>, I showed that although sample variance on average gives an unbiased estimate of population variance, it is highly unreliable at extremely small sample sizes. 

This time, I will focus more closely on the <i>distributions</i> of sample variance. <b>Does sample size affect the distributions of sample variance?</b> And how might this inform how we determine which sample sizes are too small? I'll use one of my favorite new(ish) packages, <b><a href="https://github.com/clauswilke/ggridges">ggridges</a></b>, to plot the sets of distributions from one example simulation.

</p>
<!---more--->

## ‘Behavior’ of sample variance at small sample sizes.

I’ll first re-create part of what I showed in my previous post. We’ll simulate a dataset for a ‘parent’ population and then take many random samples of increasingly larger sample sizes to get a sense of how sample variance behaves.

To cut down on time, I will only take samples at particular sample sizes, based on the sample sizes which seemed interesting (to me\!) in the previous post.

``` r
mean = 0
SD = 20
popsize = 1000

set.seed(123) # to get the same parent pop as last time
pop <- rnorm(popsize, mean, SD)

# verify that we get the same variance as last time
var(pop)
```

    ## [1] 393.3836

Checks out - great!

``` r
# pick specific Ns this time
Ns <- c(2, 3, 5, 10, 15, 30, 45, 
        90, 180, 250, 500, 750)
reps = 1000

var.p <- function(x) {
  var(x) * (length(x) - 1) / length(x)
}

# Please note that I don't use set.seed() here.
# So these samples will not be identical to those we
# got in mymat & varmat last time.
varmat = matrix(nrow = length(Ns), ncol = reps)
for (i in 1:dim(varmat)[1])
{
  for (j in 1:dim(varmat)[2])
  {
    varmat[i, j] = var(sample(pop, Ns[i]))
  }
}
rownames(varmat) <- Ns

plot(
  rep(2, ncol(varmat)),
  varmat[1,],
  ylim = c(0, max(varmat)),
  xlim = c(0, (max(Ns) + 50)),
  pch = 19,
  col = rgb(0, 0, 0, alpha = 0.2),
  xlab = 'sample size',
  ylab = 'sample variance',
  tck = 0.02,
  bty = "n"
)
for (i in 2:length(Ns)) {
  points(rep(Ns[i], ncol(varmat)),
         varmat[i,],
         pch = 19,
         col = rgb(0, 0, 0, alpha = 0.2))
}
abline(h = var.p(pop),
       col = rgb(0, 0, 1, alpha = 0.5),
       lwd = 3)
points(Ns, rowMeans(varmat),
       col = 'orange', pch = 19)
legend(300, 3000, 
       legend=c("True population variance",
                "Means of sample variance"),
       col=c(rgb(0, 0, 1, alpha = 0.5), "orange"), 
       lty=1, lwd=3, box.lty=0)
```

![](https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-05-12/sample_variance_sample_sizes_all.png)<!-- -->

This looks generally similar to the first figure from my previous post. The parent population is identical to the one I used previously, as I called `set.seed(123)` prior to each simulation. Because I didn't want to bother with seeding within `for()` loops I did not use `set.seed()` prior to sampling from the parent, which means that `varmat` will be different each time.

It may be hard to see what’s going on at the smallest sample sizes, so here’s the data at sample size \<= 15:

![](https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-05-12/sample_variance_sample_sizes_smallSS.png)<!-- -->

Of course, it’s (hopefully) very likely that no published study would try to say anything conclusive about population variance based on a sample size of 2 or 3.

## Distributions of sample variance at each sample size

We’ll now take a look at how the distributions of sample variance at each sample size (i.e. each vertical strip) vary. This will involve plotting via `ggplot2` along with the `geom_density_ridges2()` function from `ggridges` to get some nice visualizations of the distributions.

``` r
# First load | install&load packages we'll need
packages = c("ggplot2", "ggridges", "tidyr", 
             "forcats", "dplyr")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Re-organize our data in tidy format
df <- tidyr::gather(as.data.frame(t(varmat)))
# Sorting data so it appears on the y-axis in the 
# correct order is tricky in ggplot2.
# We'll use dplyr::mutate() with forcats::fct_relevel()
# to re-organize the data prior to plotting
df %>% mutate(key = fct_relevel(key, as.character(Ns))) -> df

# Use ggplot() with ggridges::geom_density_ridges2()
# The coord_flip() argument flips the axes so they are
# in the same orientation as in the previous figure.
p <- ggplot(df, aes(x = value, y = key)) +
  geom_vline(
    xintercept = var.p(pop),
    col = rgb(0, 0, 1, alpha = 0.5),
    lwd = 1
  ) +
  geom_density_ridges2(
    rel_min_height = 0.001,
    scale = 2,
    fill = rgb(0, 0, 0, alpha = 0.75)
  ) +
  coord_flip() +
  ylab("sample size") + xlab("sample variance") +
  
  theme_ridges()
p
```
See [this post](https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/) for info on the first few lines of the code chunk.

![](https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-05-12/sample_variance_distributions_ggridges.png)<!-- -->

***At small sample sizes, we see extremely skewed distributions of sample variance.***

As in the other plots, the horizontal blue line shows where the true population variance is. Please also note that the sample sizes are not evenly distributed along the x-axis.

For sample size 15 or below (among our cherry-picked examples), we’re seeing extremely long right-tailed distributions. The shapes of the distributions indicate that median and/or mode might strongly differ from the mean sample variance. Let’s take a look.

We’ll just focus on median vs. mean:

``` r
medianz <- apply(varmat, 1, median)
meanz <- apply(varmat, 1, mean)

plot(
  Ns,
  meanz,
  pch = 19,
  col = "#42AB5D",
  xlab = "sample size",
  ylab = "value",
  ylim = c(0, (max(meanz) + 100)),
  xlim = c(0, (max(Ns) + 50)),
  tck = 0.02,
  bty = "n"
)
points(
  Ns,
  medianz,
  pch = 19,
  col = "#DD3497",
  ylab = "value",
  xlab = "sample size"
)
abline(h = var.p(pop),
       col = rgb(0, 0, 1, alpha = 0.5),
       lwd = 3)
legend(
  300,
  200,
  legend = c("Means of sample variance",
             "Medians of sample variance"),
  col = c("#42AB5D", "#DD3497"),
  pch = 19,
  box.lty = 0
)
```

![](https://github.com/vbaliga/vbaliga.github.io/raw/master/images/2019-05-12/variance_mean_median.png)<!-- -->

***So at small sample sizes, the means of sample variance are close to but overshoot our population variance, whereas the medians sharply underestimate population variance***.

But at sample sizes over 45, means and medians of sample variance are nearly identical to the true population variance. It is roughly around that sample size where we see the sample variances start to show normal distributions.

So we're getting closer to understanding why small sample sizes poorly describe our toy example's population variance. The distributions of sample variance show crazy skew when generated from small samples. In my previous post, I learned that higher population sizes seem to need higher minimum sample sizes to capture population variance, but that was all based on using the `changepoint` algorithm and looking for breakpoints in the variance of sample variance. It might be a good subject to revisit by just plotting distributions as I did here. How does increasing population size affect the behavior of sample variance distributions? 

🐢
