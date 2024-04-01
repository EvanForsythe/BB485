---
layout: default
---
# Week 1 lecture and tutorial

## **Course structure**

### **Overview**:
BB485 is an advanced, upper-level course designed to give students first-hand experience performing real bioinformatics analyses akin to academic/biotech-industry research. 

This class will be structured similar to a graduate-level course that you would encounter in biology-related masters or PhD program. This will take some adjustment but will ultimately be extremely valuable in preparing students for the 'real world' (or at least one sector of the 'real world')

Compared to introductory biology courses, you can expect this course to:
- **Place a strong emphasis on projects.** Academic/industry research is driven by results. Similarly, this course will mainly be project-based. 
- **Place a less emphasis on exams and quizzes.** We will not have exams in this course. Instead, student learning is focused on developing critical-thinking, problem-solving, and trouble-shooting skills in an applied context. 
- **Final projects will constitute real, publishable research.** By the end of the term, you will contribute to real research. You'll be conducting novel analyses and asking questions that your instructor doesn't already know the answer to (nobody know the answer. That's scientific inquiry!).

### **Grading structure:**

![Points structure of BB 485](/Images/Week01/points_piechart.png)

**Tutorial Assingments:** Each Tuesday will begin with an interactive tutorial, which will include a small assignment. Students will most often complete tutorial assignments by the end of class.

**Project Write-ups:** We will spend Thursdays working on a much more in depth analysis of real biological datasets. These analyses will require much more time, trouble shooting, and computational power to accomplish. These analyses will constitute the weekly project write-ups. These weekly projects are the major emphasis of the course because they are your chance to apply your bioinformatics skillset to real-world bioinformatic analyses. 

![Assignment timeline](/Images/Week01/assignment_timeline.png)

**Project makeup work:** You will be allowed to revise and resubmit two of the six weekly project write-ups after grading. This means that you'll have a chance to improve your score on write-ups, which will have a large positive impact on your grade in the class.

![what is bioinformatics](/Images/Week01/bioinformatics.png)

## What is bioinformatics?
Bioinformatics refers to using computers to learn about biological organisms. This often involves computational analyses of 'biological data'. Biological data can be a lot of different things (from measuring the metabolic rate of a desert mouse to the length of a frog's whiskers). In this course, biological data typically relates to DNA, RNA, or protein (sequence, expression level, folding, etc...). Bioinformatics is a set of tools/skills that allow us to glean biological insights from biological data.

**What's the difference between "Bioinformatics" and "Computational Biology"?** Nothing! Some definitions may disagree, but your instructor uses the words interchangably.

**What's the best way to learn bioinformatics?** By doing it! Let's dive in....

## **Connecting to Jupyter Hub course workspace**

The link to the lab Jupyter Hub is on Canvas under Class Resources.

When you select "Start My Server" you'll get a choice of different server size options. Choose "Standard" unless instructed otherwise.

![Jupyter Hub map](/Images/Week01/JH_nav.png)

## Linux command line navigation

![Anatomy of a unix command](/Images/Week01/unix_command.png)

# Introduction to Basic Linux/Unix Commands

## A few example commands:

### 1. **`pwd`**
   - **Description:** Print the current working directory.
   - **Usage:** `pwd`

### 2. **`ls`**
   - **Description:** List directory contents.
   - **Usage:** `ls <options> <directory to view the contents of>`
   - **Options:**
     - `-l`: Long format, displaying detailed information.
     - `-a`: Include hidden files.
     
### 3. **`cd`**
   - **Description:** Change directory.
   - **Usage:** `cd <full or relative path to directory to change into>`
   - **Example:** 
     ```
     cd Documents
     ```

### 4. **`mkdir`**
   - **Description:** Create a new directory.
   - **Usage:** `mkdir <name of directory you'd like to create>`
   - **Example:** 
     ```
     mkdir Projects
     ```

### 5. **`history`**
   - **Description:** Prints previous commands you've used.
   - **Usage:** `history`
   - **Example:** 


### 6. **`cp`**
   - **Description:** Copy files and directories.
   - **Usage:** `cp <options> <source> <destination>`
   - **Example:** 
     ```
     cp file1.txt file2.txt
     ```

### 7. **`mv`**
   - **Description:** Move or rename files and directories.
   - **Usage:** `mv <source> <destination>`
   - **Example:** 
     ```
     mv file1.txt newdir/
     ```

### 8. **`rm`**
   - **Description:** Remove files and directories.
   - **Usage:** `rm <options> <file or directory to remove>`
   - **Options:**
     - `-r`: Recursively (used to remove directories).
   - **Example:** 
     ```
     rm file.txt
     ```

### 9. **`cat`**
   - **Description:** Display file content. 'cat' is short for 'concatenate' because it can also be used to combine multiple text files.
   - **Usage:** `cat <file to print to screen>`
   - **Example:** 
     ```
     cat example.txt
     ```

### 10. **`grep`**
   - **Description:** Search text patterns within files.
   - **Usage:** `grep <options> "pattern" <file to search in>`
   - **Options:**
     - `-i`: Case-insensitive search.
     - `-r`: Recursively search directories.
     - There are many more options. Grep is a powerful tool.
   - **Example:** 
     ```
     grep "pattern" file.txt
     ```

## Python and the BioPython module

# Introduction to Biopython

## Overview
Biopython is a Python library designed to enable bioinformatics tasks such as sequence analysis, molecular biology, and bioinformatics data manipulation. It provides a wide range of functionalities to work with biological data efficiently. Below are some common commands and examples to get you started with Biopython.

## BioPython commands:

### 1. **`Seq`**
   - **Description:** Represents a biological sequence (DNA, RNA, or protein).
   - **Example:** 
     ```python
     from Bio.Seq import Seq
     my_seq = Seq("AGTACACTGGT")
     ```

### 2. **Reverse Complement**
   - **Description:** Obtain the reverse complement of a DNA sequence.
   - **Example:** 
     ```python
     from Bio.Seq import Seq
     my_seq = Seq("AGTACACTGGT")
     reverse_complement = my_seq.reverse_complement()
     ```

### 3. **Translation**
   - **Description:** Translate a DNA sequence into a protein sequence.
   - **Example:** 
     ```python
     from Bio.Seq import Seq
     my_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
     protein_seq = my_dna.translate()
     ```

### 4. **`SeqRecord`**
   - **Description:** Represents a sequence with associated metadata (e.g., ID, description).
   - **Example:** 
     ```python
     from Bio.SeqRecord import SeqRecord
     from Bio.Seq import Seq
     record = SeqRecord(Seq("AGTACACTGGT"), id="1", description="Example sequence")
     ```

### 5. **`SeqIO`**
   - **Description:** Input/output operations for sequence files in various formats (FASTA, GenBank, etc.).
   - **Example:** 
     ```python
     from Bio import SeqIO
     for record in SeqIO.parse("sequences.fasta", "fasta"):
         print(record.id, len(record))
     ```

### 6. **`Align`**
   - **Description:** Perform sequence alignment.
   - **Example:** 
     ```python
     from Bio import Align
     aligner = Align.PairwiseAligner()
     alignments = aligner.align("ACGT", "ACG")
     ```

### 7. **`BLAST`**
   - **Description:** Interface to the NCBI BLAST suite for sequence similarity searching.
   - **Example:** 
     ```python
     from Bio.Blast import NCBIWWW
     result_handle = NCBIWWW.qblast("blastn", "nt", "ACGT")
     ```

### 8. **`Bio.Entrez`**
   - **Description:** Access NCBI databases including PubMed and GenBank.
   - **Example:** 
     ```python
     from Bio import Entrez
     Entrez.email = "your@email.com"
     handle = Entrez.efetch(db="nucleotide", id="71066805", rettype="gb", retmode="text")
     ```

## Data visualization in python

## Overview of programs/servers/software/databases/repositories we'll use this term
