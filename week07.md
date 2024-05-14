---
layout: default
---

<a name="top"></a>


# Week 7 lecture and tutorial
1. [Protein domains](#domains)
2. [Protein secondary structure prediction](#structure)
3. [Conserved domain prediction and visualization with R and CD-search](pred)


## <ins>**Protein domains**</ins> <a name="domains"></a>

## <ins>**Protein secondary structure prediction**</ins> <a name="structure"></a>

## <ins>**Conserved domain prediction and visualization with R and CD-search**</ins> <a name="pred"></a>

To create the figure above, we will need to generate three things:
1. A multiple sequences alignment of the protein sequences
2. A phylogenetic tree inferred from the protein alignment
3. A table of conserved domains

We will generate these files using commands we've used before (mafft and IQ-tree) as well as a web-based conserved domain prediction tool

### 1. A multiple sequences alignment of the protein sequences
Unaligned sequences for the species are located at: `/shared/forsythe/BB485/Week07/input_files/`. Choose one file (at random) for your analysis. Check with your classmates to make sure you haven't choosen the same one.

Run mafft from the command line to create a multiple sequence alignment. Hint: you may need to switch into your phylogenetics conda environment.

### 2. A phylogenetic tree inferred from the protein alignment

Next, use IQtree to infer a phylogeny from your alignment. Hint: we ran iqtree with a 'system call' from inside of a python script. Since we only need to infer one tree (as opposed to ~8000), you can call iqtree directly from the command line. There should be clues in your python script on the syntax for running iqtree from the command line.

### 3. A table of conserved domains

To predict domains across multiple species, we will use the web-based tool [CDsearch](https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi).

Input your multiple sequence alignment in


[Back to Top](#top)
