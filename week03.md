---
layout: default
---

# Week 3 lecture and tutorial
1. [Conceptual background on phylogenetics](#Homology)
   - **A.** [Homology](#Homology)
   - **B.** [Phylogenetic inference](#phylogenies)
   - **C.** [Parsimony](#parsimony)
   - **D.** [Multiple sequence alignments](#MSA)
   - **E.** [Maximum likelihood](#ML)
2. [Setting up a computing environment with conda](#conda)

## <ins>**Homology: features shared due to common ancestry**<ins> <a name="Homology"></a>

## <ins>**Phylogenetic inference**<ins> <a name="phylogenies"></a>

## <ins>**Parsimony: the simplest explaination is the best explanation**<ins> <a name="parsimony"></a>

## <ins>**Multiple sequence alignment**<ins> <a name="MSA"></a>

## <ins>**Maximum likelihood**<ins> <a name="ML"></a>

## <ins>**Setting up a computing environment with conda**<ins> <a name="conda"></a>

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
- The ability to install different software in different environments. This includes installing different <ins> versions <ins> of software, to ensure that all the software packages in a given conda environment are compatible with eachother.

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

6. Install software (answer "y" when it asks a question):
- python packages we'll use frequently
  - `conda install pandas`
  - `conda install numpy`
  - `conda install Bio`
- mulitple sequence alignment software
  - `conda install mafft`
- phylogenetics software
  - `conda install raxml`
  - `conda install newick_utilities`
 
7. Checking the status of the environment:
- Run `echo $PATH` to see how it's changed from the beginning.
- start typing mafft and use [TAB] to see if it autofills.
- If you want to exit an environment, you can run `conda deactivate` at any time.




  


