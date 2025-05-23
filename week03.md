---
layout: default
---

<a name="top"></a>

# Week 3 lecture and tutorial
1. [Conceptual background on phylogenetics](#Homology)
   - **A.** [Intro to phylogenetic trees](#trees)
   - **B.** [Homology](#Homology)
   - **C.** [Phylogenetic inference](#phylogenies)
   - **D.** [Parsimony](#parsimony)
   - **E.** [Inferring difficult relationships](#difficult)
   - **F.** [Multiple sequence alignments](#MSA)
   - **G.** [Maximum likelihood](#ML)
3. [Setting up a computing environment with conda](#conda)
4. [Preparing data for a phylogenetic analysis](#analysis)
   - **A.** [Homolog protein sequences](#seqs)
   - **B.** [Multiple sequence alignment with MAFFT](#mafft)
5. [Lecture 3 Tutorial Assignment](#tut_assign)
6. [Maximum Likelihood Phylogenetic Inference with RAxML](#raxml)
7. [Project write-up](#write)

## <ins>**Intro to phylogenetic trees**</ins> <a name="trees"></a>
This week we will learn how to perform phylogenetic analyses. Phylogenetic trees are diagrams that describe our best hypothesis for the relationship between different species/groups/populations/taxa that diverged from a common ancestor. 

![phy01](/Images/Week03/phy01.png)

Phylogenetic trees can be visualized in many different ways. Phylogentic tree infromation is stored and processed on computers in paranthetical text format **(newick format)**. Newick format is generally not very 'human readable', so we can use different software to depict the same information in the human-readable tree diagram we're used to seeing. 

![phy04](/Images/Week03/phy04.png)

Below are some of the terms used in phylogenetics. We will employ these terms as we analyze and interpret phylogenies. 

![phy05](/Images/Week03/phy05.png)

## <ins>**Homology: features shared due to common ancestry**</ins> <a name="Homology"></a>

The concept of homology is extremely important in evolutionary biology. As we'll see, phylogenetics relies on the ability to compare homologous structures/characters/traits across multiple species.

![phy02](/Images/Week03/phy02.png)

Phylogenetics relies specificaly on homologous structures/characters/traits that have changed during the evolution of the species being studied. A trait that hasn't changed doesn't provide any information about the relationships among species.

![phy03](/Images/Week03/phy03.png)

Synapomorphies are particularily important for phylogenetic analyses because they can help define a clade (i.e. they are 'phylogenetically informative'). 

**Synapomorphy:** A derived trait that is shared by two or more species (i.e. a clade). 

![phy06](/Images/Week03/phy06.png)

## <ins>**Phylogenetic inference**</ins> <a name="phylogenies"></a>
The phylogenitic trees that we're used to seeing are the output of complex statistical analyses. We can never know the true relationships of species (without a time-machine), so a phylogeny is our 'best guess' based on stastical analysis of data. The process of using statistics to make a best guess is called. <ins>phylogenetic inference</ins>.

![phy07](/Images/Week03/phy07.png)

Note that including an 'outgroup' in the analysis can be very helpful in thinking about which traits are ancestral vs derived.

- **Outgroup:** a taxon that we can safely assume is distantly related to all other species included in a phylogenetic analysis. Including an outgroup, allows us to 'root' a phylogenetic tree.

## <ins>**Parsimony: the simplest explanation is the best explanation**</ins> <a name="parsimony"></a>

Parsimony is a framework/worldview used in phylogenetics to 'choose' between alternative phylogenetic hypotheses. The idea is that the tree that requires the fewest changes/mutations to explain the data is the most parsimonious (i.e. the 'best'). 

![phy08](/Images/Week03/phy08.png)

Parsimony can be used in the 'tree search' phase of phylogenetic inference. 

![phy09](/Images/Week03/phy09.png)

Parsimony is a very simple and elegant framework. However, it may be an oversimplification of how evolution works in real life. For example, parsimony is prone to produce error when homoplasy has occured.

![phy10](/Images/Week03/phy10.png)

## <ins>**Inferring difficult relationships**</ins> <a name="difficult"></a>

Below we will discuss a famous phylogenetic relationship that has been difficult to resolve. We will learn about different statistical frameworks and how they can lead to conflicting results. 

![phy11](/Images/Week03/phy11.png)

The phylogenetic placement of whales has been notoriously tricky. Different researchers have argued about whether whales are closely or distantly related to hippos, leading this debate to be referred to as the 'whippo' debate. 

![phy12](/Images/Week03/phy12.png)

Early analyses based on parsimony-based analyses of a small number of morphological traits suggested that whales are located outside of the clade that contains hippos. However, bringing in more 'modern' analytical approaches may question this finding (discussed below). 

## <ins>**Multiple sequence alignment**</ins> <a name="MSA"></a>

The vast majority of modern phylogenetic analyses use molecular sequences as the 'trait' to study in phylogenetic analyses.

- Advantages of 'molecular phylogenetics':
   - The sequence of DNA (A, T, G, C) represents distinct character states, which are less prone to differences in interpretation/subjectivity. 
   - Each DNA nucleotide can be considered as a 'trait'. This means that the DNA sequence of a single gene contains thousands of traits.
      - Phylogenetic statistics (and all statistics) gain statistical power from larger sample size.
   - Theoretical work has yielded a detailed understanding of how DNA mutates and evolves in nature, which can help us form realistic models and incorperate them into phylogenetic inference.

![phy13](/Images/Week03/phy13.png)

DNA sequences can be homologous. The idea is that if a group of extant species have a given gene, the common ancestor must of had that gene as well.

![phy14](/Images/Week03/phy14.png)

To use DNA sequences for phylogenetic analyses, it is critical that we be able to 'line the sequences up' such that you could draw a vertical line and the lined up nucleotides would be the homologous site in the DNA sequence.

![phy15](/Images/Week03/phy15.png)

Lining up DNA (or protein) sequences across several species is called <ins>multiple sequence alignment</ins>. 

Inferring a MSA is also a statistical inference based on a set of rules/scores (i.e. a statisticasl algorithm). 

![phy16](/Images/Week03/phy16.png)

The concept of parsimony can be applied to DNA sequences in phylogenetic inference. 

![phy17](/Images/Week03/phy17.png)

## <ins>**Maximum likelihood**</ins> <a name="ML"></a>

Maximum Likelihood (ML) is a statistical framework that a assesses a given tree in terms of how 'likely' it would be to observe the data (the DNA multiple sequnce alignment) if that particular tree were true.

![phy18](/Images/Week03/phy18.png)

The maximum likelihood framework incorperates more complex models of sequence evolution than parsimony. For example, ML can include models that specifically account for homoplasy; therefore, ML is less prone to errors when homoplasy has occured in the data. 

![phy19](/Images/Week03/phy19.png)

Algorithms that perform 'tree search' analyses are often called 'hill-climbing' algorithms. The reason for this terminology has to do with our theoretical idea of 'tree space' (i.e. all of the possible trees for a given group of species).

![phy20](/Images/Week03/phy20.png)

![phy21](/Images/Week03/phy21.png)

![phy22](/Images/Week03/phy22.png)

![phy23](/Images/Week03/phy23.png)

## <ins>**Setting up a computing environment with conda**</ins> <a name="conda"></a>

A **computing environment** refers to the hardware and software available to run analyses. We don't have much control over the hardware, but we can manage the software in our environment by installing software.

## Checking your PATH

Your PATH is an important component of your computing environment. The PATH refers to the folders that your computer 'looks' in to find programs to run. For example, you can only run a given unix program (i.e. a command) from the command line if the folder that contains that program is in your PATH.

You can view your current path at any time by running:
```bash
echo $PATH
```
- Your path should contain several different full paths to folders where programs are stored. When you run a command, your computer searches in these folders (in the same order they're printed) until it finds a command with the name you entered.

You can also look up the location of a specific command you're running. For example, to find the location of the `head` command run:
```bash
which head
```
- The path that gets listed will be listed in your PATH.

**Issues with computing environment (i.e. software installation and PATH) are the largest source of frustration for new bioinfromaticians. The good news is that these are very easy issues to fix if you understand your PATH!**

## Creating **virtual environments** with **conda**

A virtual environment is a temporary computing environment that you can move in and out of very easily. You can setup multiple virtual environments on the same computer and move between them for different project that require different environments.

`conda` is a system for setting up virtual environments (called 'conda environments'). Conda provides the following benefits:
- The ability to download and install software with a single command.
- The ability to install different software in different environments. This includes installing different <ins> versions </ins> of software, to ensure that all the software packages in a given conda environment are compatible with eachother.

## Step-by-step tutorial for setting up a conda envirnonment
1. First, print your PATH so that you have a before/after comparison point.
```bash
echo $PATH
```
- Copy and paste your PATH into a text editor so we can look at it later. 
2. Load the anaconda module:
```bash
module load python/anaconda/3.11
```

3. Activate conda:
- this step has been a little weird on the HPC. One of the two commands should work.
```bash
conda activate
```
OR
```bash
source activate
```
- You'll know the above commands worked if you see `(base)` next to your command line prompt. This indicates that you're 'in' a conda environment. If you're seeing `(base)`, then that means that you don't need to repeat the above steps.

4. Create a new conda environment:
```bash
conda create -n <your-env-name> python=3.12
```
- `-n <your-env-name>` specifies the name of the environment. Name this for the specific task/analyses you're setting the env up for.
- `python=3.12` is used to specify the specific version of python to use. 3.12 is the latest stable release at the time of this tutorial. If we experience unexpected compatibility issues, we can create a different environment with older versions.
- Note, you only need to create an environment once. Once you've created an environment, you just need to activate it each time you log on to conda. 

5. Activate the conda environment:
```bash
conda activate <your-env-name>
```
- Note: while you only need to create an environment once, you need to activate it each time you'd like to use it.
  - you can get a list of your previously created environments with: `conda info --envs`
  - Note that environment names don't auto-fill when you tab.

6. Install software (The "-y" flag in each command should automatically answer "yes" to any questions):
- python packages we'll use frequently
  - `conda install -y pandas`
  - `conda install -y numpy`
  - `conda install -y conda-forge::biopython`
- mulitple sequence alignment software
  - `conda install -y bioconda::mafft`
- phylogenetics software
  - `conda install -y -c bioconda raxml`
  - `conda install -y bioconda::newick_utils`
 
7. Checking the status of the environment:
- Run `echo $PATH` to see how it's changed from the beginning.
- start typing mafft and use [TAB] to see if it autofills.
- If you want to exit an environment, you can run `conda deactivate` at any time.


## <ins>**Performing a phylogenetic analysis**</ins> <a name="analysis"></a>
We are going to perform a phylogenetic analysis using the command line programs, MAFFT (multiple sequence alignment) and raxML (phylogenetic inference). These programs impliment algorithms for performing the statistical calculations we discussed earlier. Software that uses complex algorithms often exists as a program that we call from the command line (as opposed to running within python for example). The rationale for this is that it makes it easier to submit a 'run' as a unix command from within a job submission script.

## <ins>**Homolog protein sequences**</ins> <a name="seqs"></a>

Phylogenetic analyses require that we obtain a protein homolog from each species we'd like to study. If there is exactly one homolog in each species, we call them "orthologs". However, if there has been gene duplication and multiple copies are present in some species, we call them "paralogs". "Homologs" is a catch all that includes any related protein. 

![atp_synth](/Images/Week03/atp_synth.png)

For this example, I  used BLAST to identify the protein sequence for a subunit of the ATP synthase enzyme complex. This complex is essential to life, so we expect it to be present in all organisms.

You can find a copy of the file at: `/shared/forsythe/BB485/Week03/ATP_synthase_subunitC.fasta`

Copy this file to your scratch directory and use `code` to open to file to look at it.

## <ins>**Multiple sequence alignment with MAFFT**</ins> <a name="mafft"></a>

MAFFT is very easy to run from the command line. It does not require many arguments or user-defined parameters. To run MAFFT, use the following command:

```bash
mafft <input-seq-file>.fasta
```
- Notice that the command above prints the alignment to your screen (i.e. the "standard out"). We typically want to save the outout as a file. You can easily do that using the "redirect" symbol (i.e. `>`).

```bash
mafft [input-seq-file] > [name-of-new-file]
```

Check your alignment file out using command like `more`, and `wc -l`.

## Viewing alignments

Viewing the alignment in a more human-readable format is a little more tricky. There are several GUI applications you can install on your personal computer (e.g. Mega or Geneious), but it's difficult to view alignments from the command line. 

Here is a [web-based tool from the National Center of Biotechnology Information](https://www.ncbi.nlm.nih.gov/projects/msaviewer/).

![aln_view](/Images/Week03/aln_view.png)


## <ins>**Lecture 3 Tutorial Assignment**<ins> <a name="tut_assign"></a>
 
<div style="border: 5px solid black; padding: 10px; margin: 10px 0;">
   <strong>Tutorial Assignment:</strong>
   <br />
   For the ATP sythase protein sequence file:
   <ol>
      <li>Perform multiple sequence alignment and save the output as a new fasta file</li>
      <li>How many characters (including gaps) is each sequence after alignment? Note: they should all me the same length after alignment.</li>
      <li>Which species seems to have the most gaps in it's alignment?</li>
   </ol>
   Submit the answers to the questions above on Canvas
</div>

<br />
<br />



## <ins>**Maximum Likelihood Phylogenetic Inference with RAxML**</ins> <a name="raxml"></a>

RaxML is a maximum likelihood algorithm for inferring phylogenies. Maximum likelihood is a complex statistical framework, meaning there are more arguments. I will provide a brief description of the most relevant ones. 

```bash
raxmlHPC-PTHREADS-SSE3 -s <name-of-alignment-file> -n <some-text-for-output> -m PROTGAMMALGF -p 12345 -x 12345 -f a -# 100 -T 2
```

- `raxmlHPC-PTHREADS-SSE3` is the version of the algorithm we're running.
- `-s` is the name of the input file. **This should be a MSA.**
- `-n` is used to specify what we want our output files to be named. E.g. `TEST`.
- `-m` is the model of evolution used when calculating likelihood scores. This model is what let's ML use a more 'realistic' framework (e.g. allowing for some homoplasy).
- `-p` and `-x` give raxml a 'random' number to help it decide where to begin the 'tree search'. (i.e. which tree to start with during the hill-climbing algorithm).
- `-f` tells raxml what mode/task to run.
- `-#` tells raxml how many bootstrap replicated to run.
- `-T` tells raxml how many parallel threads (i.e. CPUs) to use. You can change this to a larger number and run the command from inside of a job submission.




Use `head` and `wc -l` to check out the file that starts in `RAxML_bipartitions`. This file is in **newick** format and it contains the phylogeny that you inferred with the maximum likelihood algorithm.

## Viewing trees
Newick format trees are definitely not human-readable so we'll need some outside software to be able to view and interpret the trees.

For a quick view, there is a cool command-line tool for viewing the tree. We installed this tool as part of the 'newick utilities' package. To view your tree run:

```bash
nw_display <name-of-newick-file>
```

This gives a quick view but it is not suitible for a publication-quality figure (or an assignment-quality) figure. 

## Creating high-quality figures displaying phylogenies
### Web-based visualization
Similar to viewing alignments, there are several desktop apps you can install on your computer (e.g. Mega or FigTree). To avoid installing apps on your computer, you can use [this web-based tool from NCBI.](https://www.ncbi.nlm.nih.gov/projects/treeview/).


- **NOTE:** make sure that when you paste the contents of your newick file, it is all on one line. If you need to, use a text editor to get rid of weird spaces/new-lines that may have resulted from copy-and-pasting.
- **NOTE:** You'll need to root the tree with your outgroup (Arabidopsis thaliana) and change the layout to 'rectangle cladeogram'. 

### Rooted vs unrooted trees

![rooting](/Images/Week03/root.png)

### Python-based figure generation
We can use the Phylo module (part of Biopython) to work with phylogenies in python. This is a good option because there are fewer manual point-and-click steps needed. This allows us to reproduce out analyses easily!

Here is an example of python code that can be used to create a pdf file displaying a tree. There are many options for customizing your figures (see biopython help page for more info). 

```python
from Bio import Phylo
import matplotlib.pyplot as plt

# Read the Newick file into a tree object
newick_file = "RAxML_bipartitions.TEST"
tree = Phylo.read(newick_file, 'newick')

# Root the tree by the taxon "Arabidopsis_thaliana_ATP_synthase_subunit_C_family_protein"
root_taxon = "Arabidopsis_thaliana_ATP_synthase_subunit_C_family_protein"
tree.root_with_outgroup({'name': root_taxon})

# Set up figure and axis
fig, ax = plt.subplots(figsize=(10, 10))

# Plot the tree without X and Y axis
Phylo.draw(tree, axes=ax)
ax.set_xticks([])
ax.set_yticks([])

# Remove the box around the outside of the tree
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

# Save the plot as a PDF
output_pdf = "TEST_phylogenetic_tree.pdf"
plt.savefig(output_pdf)
```

You may need to install matplotlib module as follows:
```bash
conda install -y conda-forge::matplotlib
```

You'll need to modify the lines that specify:
1. The name of the input newick file
2. The name of the outgroup taxon used for rooting
3. The name of the pdf file being created


You can view the pdf in VScode, but first you'll need to install a pdf extension:

![VScode](/Images/Week03/VScode.png)


## <ins>**Weekly project write-up assignment**</ins> <a name="write"></a>
For the weekly project writeup, you will perform a phylogenetic analysis on of the same sequence but include "all domains of life". 

Steps:
1. Go to this website: `https://www.shoot.bio/`
2. Paste in the Arabidopsis thaliana ATP synthase protein sequence.
   - You can find this in the original fasta file from the shared directory. It should include the ">Arabidopsis_thaliana..." line plus the protein sequence directly after that.
4. Search "All domains of life"
5. Download the sequence files ("Export Sequences" at the top).
   - You'll need to download the file to your computer. Then you can open in any text editor and then use copy/paste to create a new fasta file on the HPC. The end goal is that you want a fasta file with all the (unaligned) sequences in your HPC so you can work with it for downstream steps.
7. Use mafft to make a multiple sequence alignment
8. Use the multiple sequence alignment to infer a ML tree
   - This analysis could take much longer than our example from class. You'll need to use a job submission to run your raxml command. I recommend 12 threads and 12 CPUs. Please ask me if you have questions about how to do that. 
9. Once you have a Newick file, use python to create a PDF of that file.
   - You can use to python code from the tutorial. Note that you don't need to run that python code from inside of a job.
   - Note that we probably don't want to root with Arabidopsis thaliana. However, we don't have an a priori outgroup to use for rooting in this case, so there's no obvious choice for an outgroup. Instead, root your tree by the branch that looks the longest when you visualize the tree. Using the longest branch as the root is a decent assumption when you don't have a clear outgroup.
11. Submit a writeup assignment that includes the following items:
- A PDF of the tree
- A screenshot of the MSA alignment (visualized on the NCBI sequence viewer)
- Write a "Results section" (expected length about 2-3 sentences per question) that answers the following questions based on the tree:
   - Do the plant species on the tree form a monophyletic clade? If so, what is the bootstrap support for the branch leading to the clade? If not, describe the topology of the tree and how it does not place plants as a clade.
   - Do the animal species on the tree from a monophyletic clade? If so, what is the bootstrap support for the branch leading to the clade? If not, describe the topology of the tree and how it does not place animals as a clade.
   - Has there been gene duplication of the ATP synthase protein sequence?


[Back to Top](#top)



