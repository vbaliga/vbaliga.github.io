---
layout: post
tags: ["R","GenBank","rentrez","ape","DNA","protein","tree-inference","phylogeny","FASTA"]
title: GenBank in R - Download DNA or protein sequences using the ape and rentrez packages
---  

<meta name="description" content="Batch download GenBank sequences via the ape and rentrez packages">

I realized the other day that using `ape::read.Genbank()` does not work for downloading protein sequences in batch from Genbank.

<img src="https://raw.githubusercontent.com/vbaliga/vbaliga.github.io/master/images/2019-11-26/genbank_seqs.png" alt="unaligned sequences from genbank" style="float:right;width:200px;height:97px;margin-left:30px;"> 

This post will cover how to use the `rentrez` package to download protein sequences from GenBank while also recapping how `read.Genbank()` can do a similar thing for a set of DNA seqs. I’ll actually start with the DNA example because I suspect it’s the more common use case.  

Both pipelines generate `FASTA` files, which can then be used for multiple sequence alignment etc etc.

<!---more--->


## Get the `ape` and `rentrez` packages

We’ll be using functions from `ape` and `rentrez`, so we’ll first get
these packages and load them.

As [explained here](https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/) the code below first checks to see if you have these packages. If you do, they’re simply loaded. If not, the missing packages are installed from CRAN and then loaded.

``` r
# First load | install&load packages we'll need
packages = c("ape", "rentrez")
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

    ## Loading required package: ape
    
    ## Loading required package: rentrez

## Downloading DNA sequences from GenBank

OK, DNA sequences first. We’ll use the list of COI accession IDs from [Baliga and Law (2016)](http://dx.doi.org/10.1016/j.ympev.2015.09.006) as an example.

You shouldn’t need to run `library(ape)` after running the chunk above, but please make sure you have `ape` loaded or the stuff that follows won’t work.

``` r
## import the sequence list from a CSV file
## NO "strings as factor"
coi <- read.table("https://raw.githubusercontent.com/vbaliga/vbaliga.github.io/master/images/2019-11-26/COI_BaligaLaw2016.csv",
                  quote="\"", stringsAsFactors=FALSE)

## convert to character lists
as.list(coi)$V1 -> coi_list

## see how it's formatted
str(coi_list)
```

    ##  chr [1:282] "AY662765" "KF434769.1" "GQ341585.1" "EF609278.1" ...

There are 282 GenBank accession IDs in the file, each of which corresponds to a COI sequence from a species of wrasse or parrotfish or outgroup taxon. These accession IDs were determined beforehand as being appropriate for the study.

The code is a bit old-school; there’s a more elegant way to do import a character list via `tidyverse`, but I’m opting to keep it this way to reduce the number of dependencies in this post.

Anyway, we’re ready for `ape::read.Genbank()`:

``` r
## use read.GenBank to acquire the sequences
coi_gen <- read.GenBank(coi_list, species.names = TRUE)

## but note that these are named according to the sequence ID
coi_gen
```

    ## 282 DNA sequences in binary format stored in a list.
    ## 
    ## Mean sequence length: 647.986 
    ##    Shortest sequence: 526 
    ##     Longest sequence: 712 
    ## 
    ## Labels:
    ## AY662765
    ## KF434769.1
    ## GQ341585.1
    ## EF609278.1
    ## AY662764
    ## DQ119198
    ## ...
    ## 
    ## Base composition:
    ##     a     c     g     t 
    ## 0.229 0.291 0.187 0.292 
    ## (Total: 182.73 kb)

So the `read.Genbank()` function is super handy for automating the download and organization of sequences, but note from the code above that the names given to the sequences follow the GenBank accession IDs.

In the the block below, we’ll extract the species’ names and use them instead.

``` r
## now use species names, as found within the GenBank entry
names_coi <- data.frame(species = attr(coi_gen, "species"),
                        accs = names(coi_gen))

## take a look at how this is organized
head(names_coi)
```

    ##                        species       accs
    ## 1          Abudefduf_saxatilis   AY662765
    ## 2         Abudefduf_vaigiensis KF434769.1
    ## 3        Acantholabrus_palloni GQ341585.1
    ## 4           Achoerodus_viridis EF609278.1
    ## 5 Amblyglyphidodon_leucogaster   AY662764
    ## 6      Amphilophus_citrinellus   DQ119198

``` r
## add these names 
names(coi_gen) <- attr(coi_gen, "species")

## take a peak at the names to verify
coi_gen
```

    ## 282 DNA sequences in binary format stored in a list.
    ## 
    ## Mean sequence length: 647.986 
    ##    Shortest sequence: 526 
    ##     Longest sequence: 712 
    ## 
    ## Labels:
    ## Abudefduf_saxatilis
    ## Abudefduf_vaigiensis
    ## Acantholabrus_palloni
    ## Achoerodus_viridis
    ## Amblyglyphidodon_leucogaster
    ## Amphilophus_citrinellus
    ## ...
    ## 
    ## Base composition:
    ##     a     c     g     t 
    ## 0.229 0.291 0.187 0.292 
    ## (Total: 182.73 kb)

**Noice\!**

We’re now ready for export, which we can do via `ape::write.dna()`.

``` r
## export the file to your working directory
## commenting out the line below so this post knits properly
## be sure to un-comment it so it actually runs:

#write.dna(coi_gen, "my_COI_seqs.fasta", format = "fasta")
```

You should now have an exported `.fasta` file that is ready for alignment or other analyses.

If you import it into a multiple sequence alignment program like [MEGA](https://www.megasoftware.net), you should see something like this:

![](https://raw.githubusercontent.com/vbaliga/vbaliga.github.io/master/images/2019-11-26/genbank_seqs.png)

### Using your own accession IDs

If you want to switch to your own set of accession IDs, you can supply your own `.csv` file or a character vector.

Here’s some examples with a character vector:

``` r
## make a character vector of accession IDs
my_coi_IDs <- c("AY662765",
                "KF434769.1",
                "GQ341585.1",
                "EF609278.1",
                "AY662764")

## one way
#coi_gen_chars <- read.GenBank(my_coi_IDs, species.names = TRUE)
#coi_gen_chars

## or just dump the character vector right in to read.GenBank()
coi_gen_chars_two <- read.GenBank(c("AY662765",
                                    "KF434769.1",
                                    "GQ341585.1",
                                    "EF609278.1",
                                    "AY662764"),
                                  species.names = TRUE) 
coi_gen_chars_two
```

    ## 5 DNA sequences in binary format stored in a list.
    ## 
    ## Mean sequence length: 624.4 
    ##    Shortest sequence: 587 
    ##     Longest sequence: 678 
    ## 
    ## Labels:
    ## AY662765
    ## KF434769.1
    ## GQ341585.1
    ## EF609278.1
    ## AY662764
    ## 
    ## Base composition:
    ##     a     c     g     t 
    ## 0.233 0.291 0.184 0.293 
    ## (Total: 3.12 kb)

## Downloading protein sequences from GenBank

We’ll do a similar song and dance here, but note that we can’t use `read.Genbank()` for this use case.

Instead, `rentrez::entrez_fetch()` does the trick. Again, please note that if you ran the first chunk in this post successfully, you shouldn’t need to run `library(rentrez)`. But do make sure that package is loaded\!

I don’t really work with proteins myself, so I’m going to grab a couple random accession IDs just for sake of example.

We’ll jump right into it:

``` r
## make a character vector of accession IDs
my_protein_IDs <- c("AAF73255.1",
                    "P04573.3",
                    "ACF39884.1"
                    )

## OR if you want to import a list of accession IDs from a CSV
# my_protein_IDs <- read.table("./my_protein_accessions.csv",
#                              quote="\"", stringsAsFactors=FALSE)
# ## and convert to character list
# as.list(my_protein_IDs)$V1 -> prot_list

## use entrez_fetch() to get the protein FASTAs
## set `db = "protein"` to specify the database type as protein
prot_gen <- entrez_fetch(id = my_protein_IDs,
                         db = "protein", 
                         rettype = "fasta")

## the resultant object isn't as neat as what read.GenBank() constructs
## it's basically just a character. meh.
str(prot_gen)
```

    ##  chr ">AAF73255.1 alcohol dehydrogenase class 3 [Branchiostoma lanceolatum]\nMAETAGKPITCRAAVAWEAKKPLVIETIEVAPPRAHEVRI"| __truncated__

``` r
## to write to file:
#write(prot_gen, file="my_proteins.fasta")
```

There you go\! I’m not as big a fan of `rentrez::entrez_fetch()` as I am of `ape::read.GenBank` because of how the objects are handled, but the ability to grab sequences from virtually any of the GenBank databases is handy.

That’s all\!

🐢
