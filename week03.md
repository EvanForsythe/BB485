---
layout: default
---

# Week 3 lecture and tutorial
1. [Homology](#Homology)
2. [Multiple sequence alignments](#MSA)
3. [Phylogenetic inference](#phylogenies)
4. [Setting up a computing environment with conda](#conda)

## <ins>**Homology: features shared due to common ancestry**<ins> <a name="Homology"></a>


## <ins>**Multiple sequence alignment**<ins> <a name="MSA"></a>


## <ins>**Phylogenetic inference**<ins> <a name="phylogenies"></a>


## <ins>**Setting up a computing environment with conda**<ins> <a name="conda"></a>

A **computing environment** refers to the hardware and software available to run analyses. We don't have much control over the hardware, but we can manage the software in our environment by installing software.

## Checking your PATH
Your PATH is an important component of your computing environment. The PATH refers to the folders that your computer 'looks' in to find programs to run. For example, you can only run a give unix program (i.e. a command) from the command line if the folder that contains that program is in your path.

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

**Issues with computing environment (i.e. software installation and PATH) are the largest source of frustration for new bioinfromaticians. The good news is that these are very easy issues to fix if you understand your PATH**





