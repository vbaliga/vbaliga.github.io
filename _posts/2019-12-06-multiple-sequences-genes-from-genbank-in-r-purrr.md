---
layout: post
tags: ["R","GenBank","purrr","tidyverse","ape","DNA","tree-inference","phylogeny","FASTA"]
title: Pulling multiple sequences for multiple genes from GenBank in R
---  

<meta name="description" content="Batch download GenBank sequences via the ape and rentrez packages">

Building off [my previous
post](https://vbaliga.github.io/download-genbank-dna-protein-sequences-in-R/),
I have now devised a way to not only batch download GenBank sequences
for a given gene, but also across multiple genes. This post will give a
worked-out example using the sets of genes I used to build a [phylogeny
of 220 birds](https://www.vikram-baliga.com/data-software-code#trees) as
part of [Baliga et al. (2019)](https://doi.org/10.1126/sciadv.aaw6670).

One wrinkle [I had
encountered](https://twitter.com/labrichthys/status/1202273293127966720)
was that the GenBank API enforces limits to batch requests (as it
should\!). Luckily, `purrr::slowly()` provides a way around this by introducing
pauses between successive calls.


<blockquote class="twitter-tweet" data-lang="en" data-theme="dark"><p lang="en" dir="ltr">resolved my API access issue: purrr::slowly() did the trick! introduced a 5-sec delay between requests. thanks again to <a href="https://twitter.com/jamie_lendrum?ref_src=twsrc%5Etfw">@jamie_lendrum</a> and <a href="https://twitter.com/BenjaminWolfe?ref_src=twsrc%5Etfw">@BenjaminWolfe</a> <a href="https://twitter.com/hashtag/rstats?src=hash&amp;ref_src=twsrc%5Etfw">#rstats</a> <a href="https://twitter.com/hashtag/tidyverse?src=hash&amp;ref_src=twsrc%5Etfw">#tidyverse</a> <a href="https://t.co/nTucwjICDk">pic.twitter.com/nTucwjICDk</a></p>&mdash; vikram baliga 🐢 (@labrichthys) <a href="https://twitter.com/labrichthys/status/1202730287983034368?ref_src=twsrc%5Etfw">December 5, 2019</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script> 



The code, which available at
[vbaliga/genbank\_downloadR](https://github.com/vbaliga/genbank_downloadR),
currently works but is also in a developmental stage – I hope to soon add in
functionality for getting proteins and other sequence types. A stable
version of the code that works for DNA sequences can be found [in this
release](https://github.com/vbaliga/genbank_downloadR/releases/tag/0.1.0),
as well as in the text of this post.

<!---more--->

## Get the `tidyverse` and `ape` packages

We’ll be using functions from `tidyverse` and `ape`, so we’ll first get
these packages and load them.

As [explained
here](https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/)
the code block below first checks to see if you have these packages. If
you do, they’re simply loaded. If not, the missing packages are
installed from CRAN and then loaded.

``` r
# First load | install&load packages we'll need
packages <- c("tidyverse", "ape")
package.check <- lapply(packages,
                        FUN = function(x) {
                          if (!require(x, character.only = TRUE)) {
                            install.packages(x, dependencies = TRUE)
                            library(x, character.only = TRUE)
                            }
                          }
                        )
```

## Import and organize lists of accession IDs

We’ll use the list of accession IDs from [Baliga et
al. (2019)](https://doi.org/10.1126/sciadv.aaw6670) as an example. I have made
the accession IDs available through [vbaliga/genbank\_downloadR](https://github.com/vbaliga/genbank_downloadR) 
as a `CSV` file that can be downloaded remotely (see [here](https://github.com/vbaliga/genbank_downloadR/raw/master/Baliga_et_al_2019_SciAdv_all_genes.csv)).

You shouldn’t need to run `library(tidyverse)` or `library(ape)` after
running the chunk above, but please make sure you have `tidyverse` and
`ape` loaded or the stuff that follows won’t work.

``` r
## The imported csv file should have: 
##  - column 1 = species names; example file uses "binomial" as this column name
##  - each successive column represents one gene
##  - each row is one species
accessions_batch <- read_csv("https://github.com/vbaliga/genbank_downloadR/raw/master/Baliga_et_al_2019_SciAdv_all_genes.csv")
```

Which produces a `tibble` that contains:

    ## # A tibble: 222 x 13
    ##    binomial CYTB  COI   ND1   ND2   `12S` `16S` FGB   MUSK  ODC   RAG1 
    ##    <chr>    <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr>
    ##  1 Accipit~ AY98~ AY66~ <NA>  AY98~ <NA>  KM04~ AY98~ <NA>  <NA>  <NA> 
    ##  2 "Accipi~ EU58~ AY66~ KP85~ EU58~ KF78~ <NA>  DQ88~ <NA>  DQ88~ DQ88~
    ##  3 Accipit~ U833~ HM03~ KM87~ KX08~ <NA>  KM04~ <NA>  <NA>  <NA>  <NA> 
    ##  4 Actitis~ AY89~ EU52~ AY89~ AY89~ DQ67~ DQ67~ AY69~ <NA>  <NA>  KC96~
    ##  5 Aechmop~ AY56~ AY66~ <NA>  <NA>  DQ67~ AF33~ <NA>  <NA>  <NA>  <NA> 
    ##  6 "Aegoli~ EU07~ AY66~ <NA>  EU60~ U837~ <NA>  <NA>  <NA>  <NA>  EU34~
    ##  7 Aeronau~ EU16~ DQ43~ EU16~ EU16~ EU16~ <NA>  <NA>  <NA>  GU16~ <NA> 
    ##  8 Agaporn~ AF00~ <NA>  <NA>  EU32~ EU19~ EU19~ GQ39~ <NA>  GQ39~ GQ39~
    ##  9 Aix spo~ EU58~ AY66~ <NA>  EU58~ HM06~ <NA>  <NA>  <NA>  <NA>  <NA> 
    ## 10 Alector~ FJ43~ KT18~ <NA>  EU84~ KC74~ KC98~ DQ30~ <NA>  <NA>  <NA> 
    ## # ... with 212 more rows, and 2 more variables: TGFB2 <chr>, ZENK <chr>

There are 1488 GenBank accession IDs in the file, each of which
corresponds to a gene sequence from a bird or outgroup taxon. These
accession IDs were determined beforehand as being appropriate for the
study.

The `CSV` (and corresponding `tibble` after it’s imported into R) is
organized to have column 1 as species’ names (`binomial`) and successive
columns are each devoted to one type of gene. Rows correspond to
species. Should you like to use your own sets of accessions, please
format accordingly.

We’ll now make a list. Each element will correspond to a gene, within
which will be its associated set of accession IDs. All `NA`s will be
removed.

``` r
## Now exclude the column with species names and make a list of each gene set
accessions_list <- 
    accessions_batch %>% 
  ## exclude species names
    dplyr::select(-binomial) %>% 
  ## coerce to list
    as.list() %>% 
  ## remove NAs
    lapply(function(x) x[!is.na(x)])

str(accessions_list)
```

    ## List of 12
    ##  $ CYTB : chr [1:209] "AY987308.1" "EU583321.1" "U83305.1" "AY894231.1" ...
    ##  $ COI  : chr [1:202] "AY666504.1" "AY666498.1" "HM033205.1" "EU525288.1" ...
    ##  $ ND1  : chr [1:56] "KP857864.1" "KM875999.1" "AY894282.1" "EU166921.1" ...
    ##  $ ND2  : chr [1:189] "AY987130.1" "EU583260.1" "KX083645.1" "AY894180.1" ...
    ##  $ 12S  : chr [1:156] "KF781309.1" "DQ674579.1" "DQ674555.1" "U83759.1" ...
    ##  $ 16S  : chr [1:107] "KM042918.1" "KM042935.1" "DQ674617.1" "AF339361.1" ...
    ##  $ FGB  : chr [1:124] "AY987211.1" "DQ881938.1" "AY695182.1" "GQ395348.1" ...
    ##  $ MUSK : chr [1:60] "EU739751.1" "EU739761.1" "EU739817.1" "EU739743.1" ...
    ##  $ ODC  : chr [1:79] "DQ881709.1" "GU166927.1" "GQ395349.1" "KF589090.1" ...
    ##  $ RAG1 : chr [1:113] "DQ881796.1" "KC969114.1" "EU348863.1" "GQ395351.1" ...
    ##  $ TGFB2: chr [1:91] "KP857808.1" "KP857790.1" "EU600970.1" "GQ395362.1" ...
    ##  $ ZENK : chr [1:62] "AF490154.1" "AF490161.1" "AF492508.1" "GQ395352.1" ...

## Pull multiple sequences for multiple genes `slowly()`

We’re now ready for batch downloading, again via `ape::read.GenBank()`.
As I mentioned above, the issue here is that feeding a constant set of
requests to the GenBank servers won’t work. After a few hundred, you’ll
run into a “`429 Too Many Requests`” error.

A few folks [clued me
in](https://twitter.com/labrichthys/status/1202273293127966720) to
`purrr::slowly`, which as it turns out is perfect for this use case:

<blockquote class="twitter-tweet" data-lang="en" data-theme="dark"><p lang="en" dir="ltr">purrr::slowly() is invaluable for this.</p>&mdash; Jamie Lendrum (@jamie_lendrum) <a href="https://twitter.com/jamie_lendrum/status/1202278201952784384?ref_src=twsrc%5Etfw">December 4, 2019</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script> 



So, the code block below sets up a `rate_delay()` of 2 seconds, which is
then used by `slowly()` to enforce a 2-second delay between successive
calls of `ape::read.GenBank()`. Ultimately, this is all used by
`purrr::map()` to apply this delayed downloading method (`slow_pulls()`)
to each of the genes of interest (listed in `accessions_list`).

The last part of the pipe simply renames the sequences according to
species’ binomial names instead of by accession IDs.

``` r
## We'll use purrr::slowly() along with map().
## Use rate_delay() to make R momentarily pause before running the next line.
## This prevents overloading the API.

## I've found a 2-sec pause seems to work; YMMV.
sleep_timer <- rate_delay(2)

## Now use slowly() to take ape::read.GenBank() and wait 2 secs
slow_pulls <- slowly(~ read.GenBank(.x, species.names = TRUE), 
                     rate = sleep_timer,
                     quiet = FALSE)

## Now map this function, using accessions_list
sequences_set <- 
    map(accessions_list, slow_pulls) %>%
  ## And rename according to species instead of accession IDs
    lapply(function(x) {
      attr(x, "species") -> attr(x, "names")
      return(x)
    })
```

Which reports to us each time the call completes and the rate delay kicks in:

    ## Retrying in 2 seconds.
    ## Retrying in 2 seconds.
    ## Retrying in 2 seconds.
    ## Retrying in 2 seconds.
    ## Retrying in 2 seconds.
    ## Retrying in 2 seconds.
    ## Retrying in 2 seconds.
    ## Retrying in 2 seconds.
    ## Retrying in 2 seconds.
    ## Retrying in 2 seconds.
    ## Retrying in 2 seconds.

`sequences_set` now has everything we need.

Home stretch: time to export to `fasta` files. We’ll extract the gene
names and use them to automate the naming of exported files.

``` r
## now use species names, as found within the GenBank entry
gene_names <- 
  colnames(accessions_batch)[-1] %>%
  paste0("_seqs.fasta")

## Write each set of sequences to FASTA files
invisible(lapply(1:length(sequences_set),
                 function(x) {
                   write.dna(sequences_set[[x]],
                             gene_names[[x]],
                             format = "fasta")
                 }))
```

That’s all\!

🐢
