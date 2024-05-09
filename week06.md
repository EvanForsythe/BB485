---
layout: default
---

<a name="top"></a>

# Week 6 lecture and tutorial
1. [Comparative genomics and genome evolution](#comparative)
	- **A.** [Large chromosomal mutations](#inversions)
	- **B.** [Whole genome duplication](#duplication)
2. [Emergent properties revealed by whole-genome analyses](#emergent)
	- **A.** [Phylogenetic incongruence](#incongruence)
	- **B.** [Hybridization and introgression](#int)
3. [Phylogenomic analyses](#phylogenomic)
4. [Bioinformatic pipelines](#pipelines)
5. [Collaborative code development with git and github](#git)
	- **A.** [Setting up your first github repo](#repo)
	- **B.** [Cloning your remote repo to a local machine](#clone)
	- **C.** [Configuration of your github credentials](#config)
6. [Tutorial assignment](#tut)
7. [Developing a phylogenomics pipeline](#pipe)
8. [Weekly write-up assignment](#writeup)



Developing a phylogenomics pipeline

## <ins>**Comparative genomics and genome evolution**</ins> <a name="comparative"></a>

![comp00](/Images/Week06/comp00.png)

### <ins>**Large chromosomal mutations**</ins> <a name="inversions"></a>

![comp04](/Images/Week06/comp04.png)

Comparing the 'synteny' between two genomes can reveal large mutations that have occured.

![comp05](/Images/Week06/comp05.png)

### <ins>**Whole genome duplication**</ins> <a name="duplication"></a>

![comp01](/Images/Week06/comp01.png)

Comparing genomes can help us understand other types of dramatic mutations, such as **polyploidy (i.e. whole genome duplication)**.

![comp02](/Images/Week06/comp02.png)

![comp03](/Images/Week06/comp03.png)

## <ins>**Emergent properties revealed by whole-genome analyses**</ins> <a name="emergent"></a>

A genome-wide view can reveal things that we could not detect when looking at a single region of the genome.

![comp07](/Images/Week06/comp07.png)

![covid](/Images/Week06/covid.png)

### <ins>**Phylogenetic incongruence**</ins> <a name="incongruence"></a>

Phylogenetic incongruence refers to when 'gene trees' do not match eachother. This makes it difficult to know when tree is 'correct'. BUT, the incongruence itself can also be an important clue about important processes that occured in the past!

![comp06](/Images/Week06/comp06.png)

### <ins>**Hybridization and introgression**</ins> <a name="int"></a>

![comp08](/Images/Week06/comp08.png)

![comp09](/Images/Week06/comp09.png)

## <ins>**Phylogenomic analyses**</ins> <a name="phylogenomic"></a>

![comp10](/Images/Week06/comp10.png)

![comp11](/Images/Week06/comp11.png)

![comp12](/Images/Week06/comp12.png)

![comp13](/Images/Week06/comp13.png)

## <ins>**Bioinformatic pipelines**</ins> <a name="pipelines"></a>

A bioinformatic 'pipeline' is code that performs several steps, using the output from one step as the input for the next step.

![comp14](/Images/Week06/comp14.png)

Pipelines can automate extremely complicated tasks, allowing us to do 'high-throughput' analyses and detect the emergent properties of genomes!

![comp15](/Images/Week06/comp15.png)

## <ins>**Collaborative code development with git and github**</ins> <a name="git"></a>

Git and GitHub are widely-used tools for managing the difficult task of 'version control' when developing code.

![comp16](/Images/Week06/comp16.png)

![comp17](/Images/Week06/comp17.png)

![comp18](/Images/Week06/comp18.png)

### <ins>**Setting up your first github repo**</ins> <a name="repo"></a>

![comp19](/Images/Week06/comp19.png)

![comp20](/Images/Week06/comp20.png)

![comp21](/Images/Week06/comp21.png)

![comp22](/Images/Week06/comp22.png)

### <ins>**Cloning your remote repo to a local machine**</ins> <a name="clone"></a>

![comp23](/Images/Week06/comp23.png)

1. Copy the URL from github
2. In your command line on the HPC, navigate to the folder in which you want to create the repo
3. Clone with repo with: `git clone <paste-url-here>`
4. Use the following to check what branch you're on in the repo (this also confirms you're in a repo): `git branch -a`
5. You can now make edits to files or create files
6. To 'push' the edits to the remote repo do the following steps in order:
- Add the files with: `git add *`
- Make a commit with: `git commit -m "a quick message describing what you just did"`
- Push to the remote repo with: `git push origin main`

### <ins>**Configuration of your github credentials**</ins> <a name="config"></a>

1. **Open Terminal or Command Prompt:**
   Depending on your operating system, open Terminal (on macOS and Linux) or Command Prompt (on Windows).

2. **Enter Git configuration commands:**
   Use the following commands to configure Git with your access token and GitHub username:

   ```bash
   git config --global credential.helper store
   git config --global user.name "Your GitHub Username"
   git config --global user.email "your_email@example.com"
   ```

   Replace `"Your GitHub Username"` with your GitHub username and `"your_email@example.com"` with the email associated with your GitHub account.

3. **Store the Access Token:**
   Run the following command to store your access token:

   ```bash
   git credential approve
   ```

   This will prompt you to enter your GitHub username and password. Enter your username and paste the access token you copied earlier when prompted for the password.

4. **Verify Configuration:**
   You can verify that your configuration is set up correctly by running:

   ```bash
   git config --list
   ```

## <ins>**Tutorial assignment**</ins> <a name="tut"></a>
- Clone your newly made github repo to the HPC
- Create a new text file to develop into a python script.
	- Should be named something like: `phylogenomic_analyses.py`
- Outline the major steps of your phylogenomic pipeline (see the flow chart above)
	- Do this outlining inside of your new .py script by adding a series of comments inside of the script in order.
 		- e.g. `#Create a list of unaligned files that we'll need to align`
- `git add`, `git commit`, and `git push` your repo to the remote repository.
- Submit the url of your git account on canvas as the tutorial assignment.    


## <ins>**Developing a phylogenomics pipeline**</ins> <a name="pipe"></a>
Here are the steps we'll need to accomplish:
- Import needed modules
```python
# Import needed modules
import os
import sys
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
```

- Run mafft using a "system call"
```python
aln_cmd = 'mafft --auto --quiet '+file+' > '+new_file_path
print(aln_cmd)
os.system(aln_cmd)
```


- Read in the trees and test the topology
```python
    #Read in the tree and store as phylo object
    temp_tree = Phylo.read(tree, "newick")

    #Loop through the tips in the tree to find which one contains Es (the outgroup)
    for tip in temp_tree.get_terminals():
        if "Es_" in tip.name:
            es_tip = tip
            #Stope the loop once we found the correct tip
            break
    
    #Root the tree by the outgroup taxon
    temp_tree.root_with_outgroup(es_tip)
    
    #Get a list of all terminal (aka tips) branches
    all_terminal_branches = temp_tree.get_terminals()
    
    #Loop through the branches and store the names of the tips of each
    for t in all_terminal_branches:
        if "Bs_" in t.name:
            Bs_temp=t 
        elif "Cr_" in t.name:
            Cr_temp=t
        elif "At_" in t.name:
            At_temp=t
        else:
            out_temp=t
        
    #Make lists of pairs of branches, so that we can ask which is monophyletic
    P1_and_P2=[Bs_temp, Cr_temp]
    P1_and_P3=[Bs_temp, At_temp]
    P2_and_P3=[Cr_temp, At_temp]
    

    #Use series of if/else statemetns to ask which pair in monophyletic
    if bool(temp_tree.is_monophyletic(P1_and_P2)):
        topo_str = "12top"
    elif bool(temp_tree.is_monophyletic(P1_and_P3)):
        topo_str = "13top"
    elif bool(temp_tree.is_monophyletic(P2_and_P3)):
        topo_str = "23top"
    else:
        topo_str = "Unknown"
    
```

- Create a figure
Use google to figure out how to impliment this in python

## <ins>**Weekly write-up assignment**</ins> <a name="writeup"></a>
- Use your github repo to develop a phylogenomics pipeline that performs the tasks above
- Add a final step to your pipeline that generates a pdf figure displaying the topology counts accross all of your trees. Choose whichever type of plot you think would work best (e.g. pie chart, bar chart, etc..) and google how to generate one in python.
- Add the pdf to your git repo and push to github
- Answer the following questions and submit on Cannas:
	- <ins> questions TBD <ins/> 

[Back to Top](#top)
