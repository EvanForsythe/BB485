---
layout: default
---

<a name="top"></a>


# Week 7 lecture and tutorial
1. [Protein secondary structure](#structure)
2. [Protein domains](#domains)
3. [Computational prediction of protein structure](#comp)
	- **A.** [Machine learning and pattern recognition in biology](#ml)
	- **B.** [Alphafold prediction of protein structure](#alphafold)
3. [Tutorial assignment](#tut)
4. [Conserved domain prediction and visualization with R and CD-search](#pred)
5. [Project write-up assignement](#write)
    

## <ins>**Protein secondary structure**</ins> <a name="structure"></a>

![prot01](/Images/Week07/prot01.png)

![prot02](/Images/Week07/prot02.png)

## <ins>**Protein domains**</ins> <a name="domains"></a>

![prot03](/Images/Week07/prot03.png)


## <ins>**Computational prediction of protein structure**</ins> <a name="comp"></a>

## <ins>**Machine learning and pattern recognition in biology**</ins> <a name="ml"></a>

# Introduction to Artificial Intelligence in Biology

Artificial Intelligence (AI) refers to the ability of machines to perform tasks that typically require human intelligence. This includes learning from data, recognizing patterns, making decisions, and adapting over time. One of the most powerful branches of AI is **machine learning**, where algorithms improve automatically through experience and data exposure.


## Applications of AI in Biology

Some major applications of AI in biological research include:

- **Genomics and Transcriptomics**: AI algorithms help identify gene variants, predict gene expression levels, and interpret the regulatory impact of mutations.
- **Drug Discovery**: Machine learning models can predict how molecules interact with proteins, speeding up the drug development process.
- **Medical Imaging**: AI can detect diseases from medical images such as MRIs, CT scans, and histology slides with high accuracy.
- **Epidemiology**: AI helps model the spread of diseases and predict future outbreaks using vast epidemiological data.
- **Ecology and Evolution**: AI techniques are used to infer phylogenetic relationships, model species distributions, and simulate evolutionary scenarios.

## How Machine Learning Recognizes Patterns

A central strength of machine learning is its ability to detect patterns in large, complex datasets—patterns that may be too subtle or multidimensional for humans to spot.

### Example: Gene Expression Classification

Imagine you have thousands of gene expression profiles from healthy and cancerous tissues. Each profile includes the expression levels of 20,000 genes.

A machine learning model (e.g. **neural network**) can be trained on a subset of this data labeled as "healthy" or "cancer." The algorithm learns patterns in how gene expression levels differ between the two classes.

Once trained, the model can:
- Predict whether a new, unlabeled expression profile comes from a healthy or cancerous sample.
- Identify which genes are most important for distinguishing between the two states.


## Introduction to pattern recognition in bioinformatics research: **Principal Component Analysis (PCA)** and **K-means clustering**.

#### Principal Component Analysis (PCA)

**PCA** is a statistical method used to reduce the dimensionality of large datasets while preserving as much variance as possible. In gene expression analysis:

- Each sample may have expression data for thousands of genes.
- PCA transforms the data into a smaller number of **principal components**—new axes that capture the most variance.
- This allows researchers to **visualize high-dimensional data** and sometimes observe natural groupings or outliers.

PCA is useful for **exploration and visualization**, but does not perform classification or prediction.

![pca](/Images/Week07/pca.png)

#### K-means Clustering

**K-means clustering** is an algorithm that partitions data into **K distinct groups** based on similarity.

- It starts by randomly assigning K "centroids" and then iteratively assigns each data point to the nearest centroid.
- The centroids are updated to reflect the mean of all points in each cluster, and the process repeats until convergence.
- In gene expression studies, K-means can be used to group samples or genes that show similar patterns of expression.

Like PCA, K-means is **unsupervised**—it does not use known labels (like “cancer” or “healthy”)—but it **actively groups data** based on internal structure.

#### Comparison to Machine Learning

| Feature                   | PCA (Statistical Method)                  | K-means Clustering                         | Machine Learning (e.g., SVM, Neural Network)      |
|---------------------------|--------------------------------------------|---------------------------------------------|----------------------------------------------------|
| Type                      | Unsupervised                              | Unsupervised                                | Often supervised (requires labeled data)           |
| Goal                      | Dimensionality reduction, visualization   | Grouping samples by similarity              | Prediction, classification, or regression          |
| Uses Labels               | No                                         | No                                          | Yes (for supervised learning)                      |
| Output                    | Principal components                      | Cluster assignments                         | Classification labels or probability scores        |
| Interpretability          | High (components can be examined)         | Medium (cluster centers and members)        | Variable (some models are "black boxes")           |
| Flexibility/Accuracy      | Good for simple patterns                  | Good for obvious groupings                  | Better for complex or nonlinear relationships      |


![sup_vs_unsup](/Images/Week07/sup_vs_unsup.png)

## <ins>**Alphafold prediction of protein structure**</ins> <a name="alphafold"></a>

Alphafold2 is a supervised machine learning program that predicts protein structure from primary amino acid sequence. It is trained (i.e. supervised) using known protein structures to identify patterns connecting AA sequences and folded structure.

![alpha](/Images/Week07/alpha.png)


## <ins>**Tutorial assignment**</ins> <a name="tut"></a>

1. Complete the [online tutorial from alphafold](https://www.ebi.ac.uk/training/online/courses/alphafold/) to learn about the alphafold program.

2. Install the following R packages on your conda environment with the following commands:

```bash
conda install conda-forge::r-ape -y
conda install conda-forge::r-dplyr -y
conda install bioconda::r-seqinr -y
conda install conda-forge::r-ggplot2 -y
conda install bioconda::bioconductor-ggtree -y
```

I would recommend using a job submission to install these. I'm not certain they'll work, but give them a try and we'll trouble shoot in class if it doesn't work.

<br>

## <ins>**Conserved domain prediction and visualization with R and CD-search**</ins> <a name="pred"></a>

"Conserved domains" refers to the idea that domains can be conserve in homologous proteins between species. We can use phylogentics/comparative genomics to get a sense of how the domain composition has evoloved over time. We will work on an analysis of conserved domain on Thursday.

![prot04](/Images/Week07/prot04.png)

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

## Manuscript to find figures:
Here is my [MANUSCRIPT](https://www.sciencedirect.com/science/article/pii/S0021925822000497?via%3Dihub)


## <ins>**Project write-up assignment**</ins> <a name="write"></a>

**Provide the following in your canvas submission:**

- Submit your domain evolution pdf/image.
- Describe one change that has ocurred in the evolution of your protein of interest. Write 2-3 sentances describing the evolution of this change. See the example questions below to get ideas. You don't need to answer all of the questions. You can describe any pattern that you notice. 
  - Is there a domain that is present/absent in a subset of species?
  - Has a domain undergone a tandem duplication, meaning there are two copies of the domain next to eachother?
  - Did this change occur in just one species or did it occur more anciently in the ancestor of several species?
    - Note that the beginning of the sequence ID's indicate the speices. (e.g. A_tha = Arabidopsis thaliana).
    - Has there been gene duplication of this gene?



[Back to Top](#top)
