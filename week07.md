---
layout: default
---

<a name="top"></a>


# Week 7 lecture and tutorial
1. [Protein secondary structure](#structure)
2. [Protein domains](#domains)
3. [Conserved domain prediction and visualization with R and CD-search](pred)

## <ins>**Protein secondary structure**</ins> <a name="structure"></a>

![prot01](/Images/Week07/prot01.png)

![prot02](/Images/Week07/prot02.png)


## <ins>**Protein domains**</ins> <a name="domains"></a>

![prot03](/Images/Week07/prot03.png)

"Conserved domains" refers to the idea that domains can be conserve in homologous proteins between species. We can use phylogentics/comparative genomics to get a sense of how the domain composition has evoloved over time.

![prot04](/Images/Week07/prot04.png)

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

Input your multiple sequence alignment into the input box and download the table. The easiest way to 'download' the output is to create a new tsv file on the HPC and then copy and paste from the web browser into your file.

Edit the table file so that it looks like this:

```bash
Query   Hit_type        PSSM-ID From    To      E-Value Bitscore        Accession       Short_name      Incomplete      Superfamily
A_ang_AANG005961        specific        438889  113     158     1.68696e-14     66.4916 cd22117 F-box_FBXL4      -      cl45894
```
- note that you need to manually do the following:
  - change `Hit type` to `Hit_type`
  - change `Short name` to `Short_name`
  - Remove `Q#2 - >` etc...

## Creating a tree/domain figure in R
Next we're going to create an R script to generate our figure. To do this we'll need to:
1. Setup a conda environment for running R
2. Create an R script that:
  - Loads needed packages
  - Reads in needed files
  - Creates a figure


### 1. Setting up a conda environment for running R

To aid in creating an R environment, you can use the following job script:
```bash
#!/bin/bash
#SBATCH --job-name=install
#SBATCH --ntasks-per-node=2
#SBATCH --time=24:0:0
#SBATCH --output=install.out
#SBATCH --error=install.err
#SBATCH --mail-user=evan.forsythe@osucascades.edu
#SBATCH --mail-type=END


#Bash commands to run within job
conda create -n Renv r-essentials r-base
y
conda activate Renv
conda install conda-forge::r-biocmanager
y
conda install bioconda::bioconductor-ggtree
y
conda install conda-forge::r-ape
y
conda install bioconda::r-seqinr
y
conda install -c conda-forge r-dplyr
y
```

### 2. Create an R script

Use a text editor to create an R script with
```bash
code Domain_evo.R
```

Add the following blocks and test to see if they're working as you go

- Set working directory:
```R
###Set the working directory (typically the directory where this script is stored)
#Add the full path to your desired working directory to the quotes
setwd("<the full path to the directory in which this scipt lives>")
```

- Load required packages:
```R
#Load packages by 'librarying them
library("BiocManager")
library("ggtree")
```

- Read in your phylogenetic tree:
```R
tree<-read.tree("path-to-your-newick-tree-file")
#Print the content of the variable "tree"
tree
```

- Read in your tsv file and store as R dataframe:
```R
#Read the tsv file into R
domain_df<-read.table(file = "<name of the domain tsv file you created>", header = TRUE, sep = "\t")

#Clean up this dataframe a bit
names(domain_df)[1]<-"Newick_label"
```

- Add a column to your dataframe:
```R
###Add a column that gives the length of each sequence
#Read in seq file (in a different format)
seqs2<-seqinr::read.fasta(file = <path-to-alignment file>, seqtype = "AA")

#Create a df of sequence lengths and join it to the domain data
domain_dat_full<-right_join(domain_df, data.frame(Newick_label=names(seqs2), Seq_ln=getLength(seqs2)), by = "Newick_label")
```
- Do some reformatting of the dataframe:
```R
#Change the classes in the dataframe
domain_dat_full[,1]<-paste(domain_dat_full[,1])
domain_dat_full[,4]<-as.numeric(paste(domain_dat_full[,4]))
domain_dat_full[,5]<-as.numeric(paste(domain_dat_full[,5]))
domain_dat_full[,6]<-as.numeric(paste(domain_dat_full[,6]))
domain_dat_full[,7]<-as.numeric(paste(domain_dat_full[,7]))
domain_dat_full[,8]<-paste(domain_dat_full[,8])
domain_dat_full[,9]<-paste(domain_dat_full[,9])
domain_dat_full[,10]<-paste(domain_dat_full[,10])
domain_dat_full[,11]<-paste(domain_dat_full[,11])
domain_dat_full[,12]<-as.numeric(paste(domain_dat_full[,12]))
#Make a new column that's the same as newick labels
domain_dat_full[,13]<-paste(domain_dat_full[,1])
names(domain_dat_full)[13]<-"TipLabels"
```

- Use ggtree to create the plot

```R
### Begin creating the tree/domain plot using ggplot
#Make a ggtree object 
p1<-ggtree(tree, branch.length ='none', ladderize = TRUE)

#Add tip names in as a facet
p2<-facet_plot(p1, panel='tip_labels',
               data=domain_dat_full, geom=geom_text, 
               mapping=aes(x=0, label= TipLabels), size=3)

#add seq length line
p3<-facet_plot(p2, panel = "domains", data = domain_dat_full, geom= geom_segment, 
               mapping = aes(x=0, xend=Seq_ln, y=y, yend=y), size=0.5, color='black')

#Add domains
p4<-facet_plot(p3, panel = "domains", data = domain_dat_full, geom=geom_segment, 
               aes(x=From, xend=To, y=y, yend=y, col=Short_name), size=3) +
  theme(legend.position = "right")

#Plot the final plot
pdf("output.pdf", width=10, height=10)
p4
dev.off()
```



[Back to Top](#top)