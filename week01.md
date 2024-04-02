---
layout: default
---
# Week 1 lecture and tutorial


## <ins>**Course structure**<ins>

### **Course Overview**:
BB485 is an advanced, upper-level course designed to give students first-hand experience performing real bioinformatics analyses akin to academic/biotech-industry research. 

This class will be structured similar to a graduate-level course that you would encounter in biology-related masters or PhD program. This format will take some adjustment, but will ultimately be great preparation for the 'real world' (or at least one sector of the 'real world')

Compared to introductory biology courses, you can expect this course to:
- **Place a strong emphasis on projects.** Academic/industry research is driven by results. Similarly, this course will mainly be project-based. 
- **Place a less emphasis on exams and quizzes.** We will not have exams in this course. Instead, student learning is focused on developing critical-thinking, problem-solving, and trouble-shooting skills in an applied context. 
- **Final projects will constitute real, publishable research.** By the end of the term, you will contribute to real research. You'll be conducting novel analyses and asking questions that your instructor doesn't already know the answer to (nobody know the answer. That's scientific inquiry!).

### **Grading structure:**

![Points structure of BB 485](/Images/Week01/points_piechart.png)

**Tutorial Assignments:** Each Tuesday will begin with an interactive tutorial, which will include a small assignment. Students will most often complete tutorial assignments by the end of class.

**Project Write-ups:** We will spend Thursdays working on a much more in depth analysis of real biological datasets. These analyses will require much more time, trouble shooting, and computational power to accomplish. These analyses will constitute the weekly project write-ups. These weekly projects are the major emphasis of the course because they are your chance to apply your bioinformatics skillset to real-world bioinformatic analyses. 

![Assignment timeline](/Images/Week01/assignment_timeline.png)

**Project makeup work:** You will be allowed to revise and resubmit two of the six weekly project write-ups after grading. This means that you'll have a chance to improve your score on write-ups, which will have a large positive impact on your grade in the class.



## <ins>**Bioinformatics overview**<ins>

![what is bioinformatics](/Images/Week01/bioinformatics.png)

**What is bioinformatics?** Bioinformatics refers to using computers to learn about biological organisms. This often involves computational analyses of 'biological data'. Biological data can be a lot of different things (from measuring the metabolic rate of a desert mouse to the length of a frog's whiskers). In this course, biological data typically relates to DNA, RNA, or protein (sequence, expression level, folding, etc...). Bioinformatics is a set of tools/skills that allow us to glean biological insights from biological data.

**What's the difference between "Bioinformatics" and "Computational Biology"?** Nothing! Some definitions may disagree, but your instructor uses the words interchangably.

**What's the best way to learn bioinformatics?** <ins>By doing it!<ins> Let's dive in....


## <ins>**Connecting to Jupyter Hub course workspace**<ins>

The link to the lab Jupyter Hub is on Canvas under Class Resources.

When you select "Start My Server" you'll get a choice of different server size options. Choose "Standard" unless instructed otherwise.

![Jupyter Hub map](/Images/Week01/JH_nav.png)


## <ins>**Linux command line navigation**<ins>

![Anatomy of a unix command](/Images/Week01/unix_command.png)

## <ins>**Introduction to Basic Linux/Unix Commands**<ins>

## A few example commands:
<!-- TODO: Add more, head, tail, wc -->

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

<div style="border: 1px solid black; padding: 10px; margin: 10px 0;">
   <strong>Task 1:</strong>
   <ol>
      <li>Log in to the jupyter hub (JH)</li>
      <li>Open a command line terminal in the JH</li>
      <li>Use the cd command to navigate the file system to find the example data for today's lecture (it will be in a shared folder that we all have access to).</li>
      <li>Find the full path to the location of a fasta file.</li>
      <li>Create a text file in your home directory to use for today's notes. Paste the full path to the fasta file into your notes file (we'll need the path later).</li>
   </ol>
</div>


## <ins>**Python Refresher**<ins>

<!-- TODO: Add slides and examples for strings, lists, dictionaries, strings, for loops, functions -->

<!-- TODO: Add a simple python task -->


<div style="border: 1px solid black; padding: 10px; margin: 10px 0;">
   <strong>Task 2:</strong>
   <ol>
      <li>Open a new python notebook. Name the notebook and save in your home directory.</li>
      <li>Get python to say "hello world"</li>
      <li>Create a variable that stores the full path to the fasta file (from the last task) as a string object</li>
      <li>Use the open command to create a file handle for 'reading' in the text file.</li>
      <li>Loop through the lines of the fasta file and add each sequence to a python dictionary</li>
      <li>Print the names of all the sequences in the dictionary (i.e. the 'keys' of the dictionary.</li>
   </ol>
</div>

## **Python notebooks versus python scripts**
Python notebooks are a great tool for developing new python code. However, when we're using our python code to accomplish a bioinformatics task, we typically need to be able to run that code from the command line. The best way to run python code from the command line is to put the code inside of a python script.

**Python scripts:** A python script is a text file that contains python code. The script can be named anything but should end in ".py" (e.g. my_first_script.py). You can run all of the code within a python script from the command line by running:
```python <name of python script>```

**Developing python scripts:** To develop a python script you can start writing your python code in a text file OR you can write your python code in a Jupyter notebook and then convert the notebook to a python script from the command line with: ```jupyter nbconvert --to script [YOUR_NOTEBOOK].ipynb```

<div style="border: 1px solid black; padding: 10px; margin: 10px 0;">
   <strong>Task 3:</strong>
   <ol>
      <li>Create a python script that does what your python notebook (from the last task) does.</li>
      <li>Run the script from the command line.</li>
   </ol>
</div>



## <ins>**The BioPython module**<ins>
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

## <ins>**Data visualization in python**<ins>

## <ins>**Overview of programs/servers/software/databases/repositories we'll use this term**<ins>



