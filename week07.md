---
layout: default
---

<a name="top"></a>


# Week 7 lecture and tutorial
1. [Protein secondary structure](#structure)
2. [Protein domains](#domains)
3. [Conserved domain prediction and visualization with R and CD-search](#pred)
    - **A.** [Code update](#update)
5. [Project write-up assignement](#write)
    

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

- Input your **unaligned sequences** into the input box and download the table.
- Note: you'll need to download your results to your personal computer using the download button (see picture below). You'll do work in R on your computer, so it's ok to leave this file on your local computer to use in the R step.

![cdsearch](/Images/Week07/cdsearch.png)
- Note: be sure to download the 'standard' format file


The downloaded file will come with some uneeded metadata at the top. You can remove the unneeded lines at the top that look like this:
```bash
#Batch CD-search tool	NIH/NLM/NCBI
#cdsid	QM3-qcdsearch-65209C2A234ECDF-1A43972510CBD590
#datatype	hitsConcise Results
#status	0
#Start time	2024-05-18T15:16:16	Run time	0:00:00:02
#status	success
```

Edit the remaining table file so that it looks like this:

```bash
Query	Hit_type	PSSM-ID	From	To	E-Value	Bitscore	Accession	Short_name	Incomplete	Superfamily
A_tha_AT3G07310.1	specific	461674	174	244	8.71002e-19	79.8209	pfam05542	DUF760	 - 	cl09377
A_tha_AT5G48590.1	specific	461674	162	240	1.71747e-16	73.2725	pfam05542	DUF760	 - 	cl09377
A_tha_AT5G48590.1	superfamily	471877	254	321	0.0046136	35.5229	cl09377	DUF760 superfamily	 - 	 - 
E_sal_Thhalv10020977m	specific	461674	174	244	6.28311e-21	85.5989	pfam05542	DUF760	 - 	cl09377
E_sal_Thhalv10004513m	specific	461674	170	240	1.80105e-14	67.8797	pfam05542	DUF760	 - 	cl09377
O_sat_LOC_Os11g26890.1	specific	461674	171	252	3.73386e-13	64.4129	pfam05542	DUF760	 - 	cl09377
S_lyc_Solyc09g057580.3.1	specific	461674	179	249	4.22176e-18	78.2801	pfam05542	DUF760	 - 	cl09377
S_lyc_Solyc06g007560.2.1	specific	461674	120	187	2.8519e-12	61.3313	pfam05542	DUF760	 - 	cl09377
S_lyc_Solyc06g007550.1.1	specific	461674	108	183	4.76798e-14	65.9537	pfam05542	DUF760	 - 	cl09377
S_pol_Spipo9G0068100	specific	461674	275	342	2.39807e-13	65.1833	pfam05542	DUF760	 - 	cl09377
```
- note that you need to manually do the following:
  - change `Hit type` to `Hit_type`
  - change `Short name` to `Short_name`
  - Remove `Q#2 - >` etc...
  - Note: if you have rows that seem formatted 'weird' in any way, you can remove those rows.

## Creating a tree/domain figure in R
1. Download your sequences, tree, and domain table to your laptop.
  - You can download files using the scp command (described in Week 2).
  - If scp is causing problems, you can also create empty text files and copy and paste the contents of your HPC files into them (this way is not preferred but will work in a pinch). 
2. Create an R script file and paste the R code from below to generate a plot. 

Here is an R script that you can use to create a domain evolution figure. You'll need to change information about the names of files. I recommend using R-studio on your local computer to run the R code.

## R code:
NOTE: this code includes extra packages, that could cause issues when attempting to install them. See the block of code further down on this page for an alternative R script, which only includes the minimal set of packages.
```R
###Set the working directory (typically the directory where this script is stored)
#Add the full path to your desired working directory to the quotes
setwd("<location-of-your-R-script>")


#The BiocManager package is needed in order to install some other packages
#This is an if-statement that asks whether the package is already installed and installs if not
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
library("BiocManager")

#Install and load treeio if it isn't already
if (!requireNamespace("treeio", quietly = TRUE)){
  BiocManager::install("treeio")

}
library(treeio)

#Install and load ggtree if it isn't already
if (!requireNamespace("ggtree", quietly = TRUE)){
  BiocManager::install("ggtree")
}
library(ggtree)

###Install and load several packages that are installed with the standard base install function
#Make a list of package names that we'll need
package_list<-c("ape", "ips", "Biostrings", "phytools", "seqinr", "dplyr", "ggplot2")

#Loop to check if package is installed and loaded. If not, install/load
#If you get a warning saying "there is no package called <XYZ>", run the loop again
for(k in 1:length(package_list)){
  
  if (!require(package_list[k], character.only = TRUE)) {
    install.packages(package_list[k], dependencies = TRUE)
    library(package_list[k], character.only=TRUE)
  }
}

#Read in the tree file
tree<-read.tree("path-to-your-newick-tree-file")

#Print the content of the variable "tree"
tree

#Read in the domain table and store as R dataframe 
domain_df<-read.table(file = "<name of the domain tsv file you created>", header = TRUE, sep = "\t")

#Clean up this dataframe a bit
names(domain_df)[1]<-"Newick_label"

###Add a column that gives the length of each sequence
#Read in seq file (note: this should not be the aligned sequences)
seqs<-seqinr::read.fasta(file = <path-to-unaligned-seqs-file>, seqtype = "AA")

#Create a df of sequence lengths and join it to the domain data
domain_dat_full<-left_join(domain_df, data.frame(Newick_label=names(seqs), Seq_ln=getLength(seqs)), by = "Newick_label")

#check out the dataframe
domain_dat_full

# Do some reformatting of the dataframe
#Change the classes in the dataframe so that R can recognize the numbers as number and the strings as strings.
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
p4

# You can now use the "Export" button to export a pdf or image file of your tree and domains figure.
# Note: part of the sequence IDs may be cutoff. You can stretch the image window or change the dimensions of the image/pdf size to help fix this a bit, but it's ok if some of the IDs are slightly cutoff for this assignment. We could use adobe illustrator to manually do some post-processing and fix this issue if we were preparing a figure for a publication.
```

### <ins>**Code update**</ins> <a name="update"></a>

UPDATE (may 22nd): Below is a block of code that seems less likely to cause installation problems. Note, ggtree is still causing issues. If you're unable to install ggtree, we may need to take a pass on the final ggtree step. If you are unable, to install ggtree, send me an email telling me the name of your input file (e.g. OG0009631.fasta). I'll plot the tree for you and send it to you and you can use that plot to answer the assignment questions.

Here is an updated block of code to give one last try.
```R
#Add the full path to your desired working directory to the quotes
setwd("path-to-folder-where-this-script-lives")

#Make a list of package names that we'll need
package_list<-c("ape", "dplyr", "seqinr", "ggplot2", "ggtree")

#Loop to check if package is installed and loaded. If not, install/load
#If you get a warning saying "there is no package called <XYZ>", run the loop again
for(k in 1:length(package_list)){
  
  if (!require(package_list[k], character.only = TRUE)) {
    install.packages(package_list[k], dependencies = TRUE)
    library(package_list[k], character.only=TRUE)
  }
}

# Read in the tree file
tree <- ape::read.tree("name-of-your-tree.treefile")

# Read in the domain table and store it as an R dataframe
domain_df <- read.table(file = "name-of-your-domains-table.tsv", header = TRUE, sep = "\t")

# Clean up this dataframe a bit
names(domain_df)[1] <- "Newick_label"

# Read in unaligned sequences
seqs <- seqinr::read.fasta(file = "name-of-your-seq-file", seqtype = "AA")

# Create a dataframe of sequence lengths and join it with the domain data
domain_dat_full <- left_join(domain_df, 
                             data.frame(Newick_label = names(seqs), Seq_ln = seqinr::getLength(seqs)), 
                             by = "Newick_label")

# Check out the dataframe
print(head(domain_dat_full))

# Do some reformatting of the dataframe
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


# Check the modified dataframe
print(head(domain_dat_full))

### NOTE: If you can't succesfully install ggtree, you can skip the steps below. 
### Instead, write you domain_dat_full dataframe as a tsv file using:
### write.table(domain_dat_full, file = "new_file_name.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
### Note that this should be similar to your orignal file but it should inlcude an extra column, "Seq_ln"

### Plotting the tree and domains using ggtree and ggplot2

# Create a ggtree object
p1 <- ggtree(tree, branch.length = 'none', ladderize = TRUE)

# Add tip names as a facet
p2 <- facet_plot(p1, panel = 'tip_labels', data = domain_dat_full, geom = geom_text, 
                 mapping = aes(x = 0, label = TipLabels), size = 3)

# Add sequence length line
p3 <- facet_plot(p2, panel = "domains", data = domain_dat_full, geom = geom_segment, 
                 mapping = aes(x = 0, xend = Seq_ln, y = y, yend = y), size = 0.5, color = 'black')

# Add domains
p4 <- facet_plot(p3, panel = "domains", data = domain_dat_full, geom = geom_segment, 
                 aes(x = From, xend = To, y = y, yend = y, col = Short_name), size = 3) +
  theme(legend.position = "right")

# Plot the final plot
print(p4)
```



<!---
### NOTE TO SELF: I had a hard time installing all the needed R packages on the HPC. Pivoting to having the students use R studio on their own computers.

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
-->


## <ins>**Project write-up assignement**</ins> <a name="write"></a>

**Provide the following in your canvas submission:**

- Submit your domain evolution pdf/image.
- Describe one change that has ocurred in the evolution of your protein of interest. Write 2-3 sentances describing the evolution of this change. See the example questions below to get ideas. You don't need to answer all of the questions. You can describe any pattern that you notice. 
  - Is there a domain that is present/absent is a subset of species?
  - Has a domain undergone a tandem duplication, meaning there are two copies of the domain next to eachother?
  - Did this change occur in just one species or did it occur more anciently in the ancestor of several species?
    - Note that the beginning of the sequence ID's indicate the speices. (e.g. A_tha = Arabidopsis thaliana).
    - Has there been gene duplication of this gene?



[Back to Top](#top)
