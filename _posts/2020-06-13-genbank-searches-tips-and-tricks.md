---
layout: post
tags: ["GenBank","DNA","tree-inference","phylogeny","FASTA"]
title: Searching for sequences on GenBank - a few tips and tricks
--- 

<meta name="description" content="Some things I've learned from performing searches for DNA 
sequences on GenBank">

This will be a "low-tech" post, i.e. one that doesn't showcase code. 

Over the past few years, I've done a lot of searches for genetic sequences on GenBank with the 
endgame of building phylogenetic trees. Here's a list of things that have & haven't worked
for me when it has come to the task of finding specific genes for specific taxa.

A few good ways to set up the search term for nuclear genes are (ordered by specificity):
- (Genus_species) NOT (whole genome) NOT predicted 
- (Genus_species) NOT (whole genome) NOT predicted NOT mitochondri* 
- (Genus_species) NOT (whole genome) NOT predicted NOT mitochondri* 
AND (gene1 OR gene2)

Note 1: I always avoid including sequences that have the word "predicted" in the title; I'd 
rather rely on sequences that have established identity. That said, including `NOT predicted`
can be dangerous as the word "predicted" may appear e.g. in the title of the corresponding 
paper, but the sequence itself may be known with more certainity.

Note 2: For similar reasons, `NOT (whole genome)` may fare better than `NOT genome`

Note 3: `OR` needs to be nested within parenthetical statements. Taking the third bullet point
as an example, a search without the the parentheses surrounding the `OR` statement such as 
`(Genus_species) NOT (whole genome) NOT predicted NOT mitochondri* AND gene1 OR gene2` would 
amount to searching for `(Genus_species) NOT (whole genome) NOT predicted NOT mitochondri* AND gene1` 
OR `gene2`. Therefore, all cases of gene2 for every species ever would appear in the search results!

Perhaps this list of notes will grow more someday, but that's all for now 

🐢
