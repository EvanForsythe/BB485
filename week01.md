---
layout: default
---
# Week 1 lecture and tutorial

## <ins>**Course structure**<ins>

### **Course Overview**:
BB485 is an advanced, upper-level course designed to give students first-hand experience performing real bioinformatics analyses akin to academic/biotech industry research. 

This class will be structured similar to a graduate-level course that you would encounter in biology-related masters or PhD program. This format will take some adjustment, but will ultimately be great preparation for the 'real world' (or at least one sector of the 'real world')

Compared to introductory biology courses, you can expect this course to:
- **Place a strong emphasis on projects.** Academic/industry research is driven by results. Similarly, this course will mainly be project-based. 
- **Place a less of an emphasis on exams and quizzes.** We will not have exams in this course. Instead, student learning is focused on developing critical-thinking, problem-solving, and trouble-shooting skills in an applied context. 
- **Include a final project that constitutes real, publishable research.** By the end of the term, you will contribute to real research. You'll be conducting novel analyses and asking questions that your instructor doesn't already know the answer to (nobody knows the answer! That's scientific inquiry!).

### **Grading structure:**

![Points structure of BB 485](/Images/Week01/points_piechart.png)

**Tutorial Assignments:** Each Tuesday will begin with an interactive tutorial, which will include a small assignment. Students will most often complete their tutorial assignments by the end of class.

**Project Write-ups:** We will spend Thursdays working on a much more in depth analysis of real biological datasets. These analyses will require much more time, trouble shooting, and computational power to accomplish. These analyses will constitute the weekly project write-ups. These weekly projects are the major emphasis of the course because they are your chance to apply your bioinformatics skillset to real-world bioinformatic analyses. 

![Assignment timeline](/Images/Week01/assignment_timeline.png)

**Project makeup allowance:** You will be allowed to revise and resubmit two of the six weekly project write-ups after grading. This means that you'll have a chance to improve your score on write-ups, which will have a large positive impact on your grade in the class.

<br />
<br />

## <ins>**Overview of bioinformatics**<ins>

![what is bioinformatics](/Images/Week01/bioinformatics.png)

<br />
<br />

**What is bioinformatics?** Bioinformatics refers to the use of computers to learn about biological organisms. This often involves computational analyses of 'biological data'. Biological data can be a lot of different things (from measuring the metabolic rate of a desert mouse to the length of a frog's whiskers). In this course, biological data typically relates to DNA, RNA, or protein (sequence, expression level, folding, etc...). Bioinformatics is a set of tools/skills that allow us to glean biological insights from biological data.

<br />
<br />

## Biological data
![computers](/Images/Week01/computers.png)

<br />
<br />

**What's the difference between "Bioinformatics" and "Computational Biology"?** Nothing! Some definitions may disagree, but your instructor uses the words interchangably.

**What's the best way to learn bioinformatics?** By doing it! Let's dive in....

<br />
<br />

## <ins>**Connecting to Jupyter Hub course workspace**<ins>

We will use remote servers in this course, which we will connect to via the internet.

![servers](/Images/Week01/servers.png)

## How to access the class Jupyter Hub:

The link to the lab Jupyter Hub is on Canvas under 'Class Resources'.

When you select "Start My Server" you'll get a choice of different server size options. Choose "Standard" unless instructed otherwise.

![Jupyter Hub map](/Images/Week01/JH_nav.png)

<br />
<br />

## <ins>**Unix/Linux command line navigation**<ins>

![Unix and linux](/Images/Week01/unix_linux.png)

I will sometimes use the terms linux/unix interchangeably. However, they are different things. Here is my understanding of the differences:
- **Linux:** An operating system (similar to mac or windows). Most super computers/clusters run the linux operating system.
- **Unix:** A coding language (similar to python or R). When running commands from the command line, you are technically running unix code.

### Here is the anatomy of a typical shell command:

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

### 5. **`more`**
   - **Description:** View the contents of a file one screen at a time.
   - **Usage:** `more <file>`
   - **Example:** 
     ```
     more example.txt
     ```

### 6. **`head`**
   - **Description:** Display the first few lines of a file.
   - **Usage:** `head <options> <file>`
   - **Options:**
     - `-n`: Specify the number of lines to display.
   - **Example:** 
     ```
     head -n 10 file.txt
     ```

### 7. **`tail`**
   - **Description:** Display the last few lines of a file.
   - **Usage:** `tail <options> <file>`
   - **Options:**
     - `-n`: Specify the number of lines to display.
   - **Example:** 
     ```
     tail -n 20 file.txt
     ```

### 8. **`wc`**
   - **Description:** Count the number of lines, words, and characters in a file.
   - **Usage:** `wc <options> <file>`
   - **Options:**
     - `-l`: Count lines.
     - `-w`: Count words.
     - `-c`: Count characters.
   - **Example:** 
     ```
     wc -l file.txt
     ```

### 9. **`cp`**
   - **Description:** Copy files and directories.
   - **Usage:** `cp <options> <source> <destination>`
   - **Example:** 
     ```
     cp file1.txt file2.txt
     ```

### 10. **`mv`**
   - **Description:** Move or rename files and directories.
   - **Usage:** `mv <source> <destination>`
   - **Example:** 
     ```
     mv file1.txt newdir/
     ```

### 11. **`rm`**
   - **Description:** Remove files and directories.
   - **Usage:** `rm <options> <file or directory to remove>`
   - **Options:**
     - `-r`: Recursively (used to remove directories).
   - **Example:** 
     ```
     rm file.txt
     ```
   - **IMPORTANT NOTE:** rm perminantly deletes files. Double-check what you're doing before you use rm. 

### 12. **`cat`**
   - **Description:** Display file content. 'cat' is short for 'concatenate' because it can also be used to combine multiple text files.
   - **Usage:** `cat <file to print to screen>`
   - **Example:** 
     ```
     cat example.txt
     ```

### 13. **`grep`**
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
<br />
<br />

## Navigating a linux file system and visualizing a 'file tree'

![file navigation](/Images/Week01/file_nav.png)

<br />
<br />

<div style="border: 1px solid black; padding: 10px; margin: 10px 0;">
   <strong>Task 1: Navigating the file system.</strong>
   <ol>
      <li>Log in to the jupyter hub (JH) and open a command line terminal</li>
      <li>Use the cd command to navigate the file system to find the example data for today's lecture (it will be in a shared folder that we all have access to).</li>
      <li>Find the full path to the location of a DNA fasta file.</li>
      <li>Create a text file in your home directory to use for today's notes. Paste the full path to the fasta file into your notes file (we'll need the path later).</li>
   </ol>
</div>

<br />
<br />

## <ins>**Python Refresher**<ins>

![Python in bioinformatics](/Images/Week01/python.png)

## Short list of common python commands

### 1. **`range`**
   - **Description:** Generate a sequence of numbers.
   - **Usage:** 
     - `range(<stop position number>)` Generates numbers from 0 to stop-1.
     - `range(<start position number>, <stop position number>)`Generates numbers from start to stop-1.
     - `range(start position number>, <stop position number>, <increment/step size>)`Generates numbers from start to stop-1 with step increments.
   - **Example:** 
     ```python
     for i in range(0, 50, 3):
        print(i)
     ```

### 2. **`print`**
   - **Description:** Print messages or variables to the console.
   - **Usage:** `print(<message or variable>)`
   - **Example:** 
     ```python
     print("Hello, world!")
     ```

   - **`print` with f-strings:** 
     - **Description:** Print strings with variables embedded using f-strings.
     - **Usage:** `print(f"<string> {<variable>} <string>")`
     - **Example:** 
       ```python
       name = "Alice"
       age = 30
       print(f"My name is {name} and I am {age} years old.")
       ```

### 3. **`for`**
   - **Description:** Iterate over a sequence (list, tuple, string, etc.).
   - **Usage:** `for <variable> in <sequence>:`
   - **Example:** 
     ```python
     fruits = ["apple", "banana", "cherry"]
     for fruit in fruits:
         print(fruit)
     ```

### 4. **`len`**
   - **Description:** Get the length of a sequence (list, tuple, string, etc.).
   - **Usage:** `len(<sequence>)`
   - **Example:** 
     ```python
     fruits = ["apple", "banana", "cherry"]
     print(len(fruits))
     ```

### 5. **`replace method`**
   - **Description:** replace is a method that can be applied to string objects. Note that the syntax for 'methods' is that they come after the object they're applied to.
   - **Usage:** `my_string.replace(<string to replace>, <replace with>)`
   - **Example:** 
     ```python
     dna_seq = "ATG GTGCCAGA GTAAGC G"
     dna_seq.replace(" ", "")
     print(dna_seq)
     ```

<br />
<br />

<div style="border: 1px solid black; padding: 10px; margin: 10px 0;">
   <strong>Task 2: Working with string objects. </strong>
   <ol>
      <li>Open a new python notebook. Name the notebook and save in your home directory.</li>
      <li>Store this DNA sequences as a variable: "ATGGAGGACCCTTTGTTGACTCAGAGTGAGCACATCGTCGATGACGTTACAATCCATGGC
GATTCTTCTTCAAATGAAGAGCACATCGTCGACGTTACAACCAATGGCAATCCTTCATCA
GCTGATGAGAAAAGACCGCATGAGGGTGTCCAATGGAGTGATATATTTACATTTACCACT
GTTTGTATTTTGGTCGAGTTTGTTGTGGCTTTAGTCCAGATTGTTGCCGCCATTGTTGTT
CTGACTCTGGCAAAAGATGAACAACCTCCACAAAAAATGTTTCCTACACTGATCCTCAGT
TATACCGGTTGCTGTATTGCCACACTCCCTATTCTAGGTTTGCGTTTCTGGCATTCTTAC
CGAAGTGTTAGCACAGAGACAAGAATCTACGAGGTGGTGGACATTTTGAAAAAGATGCTT
GAATATTTCTTCGTGGGTTGGGTTGTGGTGCTTTTATGGCATCTTATCAACAACTCATCA
TCTATAGATAATACTACGCAACAGTTCTGGTTATGTATGACTTTCCTTGCTATCAGCTGC
ATTCTACATGTTCTTCGTAATCTCCCCTGTGCGGGAGTTTGTTTTCTGTATCCTATGATA
CTATATCTTTCCCAATCGATAGACTTCGTTGGTGACATTACTGATGAGATAAATTTGACT
ACGTCTATAATCCTATTATGCTTTGGAATCTTTGCCTGCATTATCTGTGGTTGTTGTTCC
AGATGCTTATGCAGATAA"</li>
      <li>Loop through the sequence and print the 3-letter sequence of each codon.</li>
   </ol>
</div>

<br />
<br />

## Reading in a fasta file using a file handle

Below is a block of python code that can be used to read in a fasta file and create a dictionary object (object named "seq_dict")

```python
#Create a file handle
seq_handle = open(seq_file_path, “r”)

#Create an empty dictionary
seq_dict = {}
#Loop through the line in the file
for line in seq_handle:
    if line.startswith(“>”):
        id_temp = line.strip() #Removes “\n”
        id_clean = id_temp.replace(“>”, “”) #Removes “>” by replacing it with nothing.
        #Add the item to the dictionary
        seq_dict[id_clean]=“” # id_clean is the key, the value is an empty string (for now)
    else:
        seq_line = line.strip() #Removes “\n”
        #append this line to the dictionary value, using the key (which is still “id_clean” from the previous line)
        seq_dict[id_clean] += seq_line
```

<br />
<br />

<div style="border: 1px solid black; padding: 10px; margin: 10px 0;">
   <strong>Task 3: reading in a fasta file and storing sequences in a dictionary object.</strong>
   <ol>
      <li>Open a new python notebook. Name the notebook and save in your home directory.</li>
      <li>Create a variable that stores the full path to the fasta file (from the last task) as a string object</li>
      <li>Use the block of python code above to read in fasta file and create a dictionary. </li>
      <li>Print the names of all the sequences in the dictionary (i.e. the 'keys' of the dictionary.</li>
       <li>Loop through each of the sequences and print a statement like "The sequence [sequence ID] is [length of sequence] long." for each sequence in the dictionary.</li>
   </ol>
</div>

<br />
<br />


## **Python notebooks versus python scripts**
Python notebooks are a great tool for developing new python code. However, when we're using our python code to accomplish a bioinformatics task, we typically need to be able to run that code from the command line. The best way to run python code from the command line is to put the code inside of a python script.

**Python scripts:** A python script is a text file that contains python code. The script can be named anything but should end in ".py" (e.g. my_first_script.py). You can run all of the code within a python script from the command line by running:
```python <name of python script>```

**Developing python scripts:** To develop a python script you can start writing your python code in a text file OR you can write your python code in a Jupyter notebook and then convert the notebook to a python script from the command line with: ```jupyter nbconvert --to script [YOUR_NOTEBOOK].ipynb```

<br />
<br />

<div style="border: 1px solid black; padding: 10px; margin: 10px 0;">
   <strong>Task 4: running python code from a python script.</strong>
   <ol>
      <li>Create a python script that does what your python notebook (from the last task) does.</li>
      <li>Run the script from the command line.</li>
   </ol>
</div>

<br />
<br />

**Making python scripts 'executable':** It is sometimes cleaner/easier to run your script without the need to put type "python" before the name of your script. You can accomplish this by making your script 'executable'.

- Add the 'shebang' line to the top of your script. This must be the first line in your script: ```#!/usr/bin/python3```.
- Next, you'll need to change the 'permissions' on your python script file.
   - To see the permission on your python file, use ```ls -l```
   - To change the permissions to make you file executable: ```chmod +x <name of script>```
   - Run ```ls -l``` again to see how it changed.
   - Now you can run your python script like this: ```./<name of script>```

<br />
<br />
 
<div style="border: 5px solid black; padding: 10px; margin: 10px 0;">
   <strong>Tutorial Assignment:</strong>
   In a text file, address each of the questions below. Submit your responses on Canvas.
   <ol>
      <li>We you able to successfully complete each of the four tasks above? If not, what issues did you encounter? Work with me to find solutions to any issues by the end of the week.</li>
      <li>Please install VScode (https://code.visualstudio.com/download) on your laptop. We will need to use it next week. Were you able to install it? If not, that's OK , we can talk about solutions/work-arounds.</li>
      <li>Please install Slack (google 'slack download') on your laptop. We will need to use it next week. Were you able to install it? If not, that's OK, we can talk about solutions/work-arounds.</li>
   </ol>
</div>

<br />
<br />

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

### 7. **`SeqIO.to_dict`**
   - **Description:** The SeqIO.to_dict function is used in conjunction with the SeqIO.parse function to read in a fasta file and store as a dictionary.
   - **Example:** 
     ```python
     from Bio import SeqIO
     input_file = open("input.fasta", "r")
     my_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
     ```

<br />
<br />

## Python Modules

**BioPython is a python module.** In python, modules are a collection of python functions. The functions that are automatically loaded every time you use python are called 'base python' function. When we need to do more specialized things in python, we can explicitly 'import' a module with specialized python functions. In addition to BioPython, some exmaples of common modules we'll use this term are, NumPy (common math/numerical functions), pandas (dataframes for working with data tables), and matplotlib (plotting figures).

## Importing modules
To import a module, you typically add the following to the top of you python script/notebook: `import Bio`

Sometimes is cleanest to import specific functions as well so that you can easily call them by name later. Like this:
```python
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
```
- *Full dislosure: I don't have a great system for when I use from/import versus import. TBH, I just follow the conventions that it seems like most people are using on the internet for a particular module.*
- *On a related note, people seem to frequently give nicknames to some modules using import/as. I also tend to follow those conventions as well. Example: `import pandas as pd`*

<br />
<br />

<div style="border: 1px solid black; padding: 10px; margin: 10px 0;">
   <strong>Task 5: working with sequences using BioPython functions</strong>
   <ol>
      <li>Use your python notebook you started earlier this week to work with BioPython functions</li>
      <li>Import the following functions: Seq, SeqRecord, SeqIO</li>
      <li>Convert your string object into a BioPython Seq object.</li>
      <li>Convert your BioPython Seq object into a BioPython SeqRecord object.</li>
      <li>Add a sequence ID to your BioPython SeqRecord object</li>
      <li>Read in the fasta file from earlier and store as a dictionary. This time do it using BioPython SeqIO.parse and SeqIO.to_dict functions.</li>
      <li>Loop through the first ten sequences in the dictionary and translate the each sequence to protein sequence.</li>
   </ol>
</div>

<br />
<br />


## <ins>**Data visualization in python**<ins>
Creating figures/graphs/tables is an extremely important part of bioinformatics. The 'biological insights' (i.e. the end goal of bioinformatic analyses) usually come from being able to visualize data in an informative way that reveals patterns/trends in the data.

## Plotting data from data tables (i.e. pandas dataFrames)

#### Step 1: Import necessary modules

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
```

#### Step 2: Create a DataFrame

```python
# Generate 20 random age measurements between 50 and 70
random_ages = np.random.randint(20, 71, size=50)

data = {'Age': random_ages}
df = pd.DataFrame(data)
```

#### Step 3: Plot the Histogram

Now, let's plot the histogram based on the "Age" column in our DataFrame:

```python
plt.hist(df['Age'], bins=10, color='skyblue', edgecolor='black')
plt.title('Age Distribution')
plt.xlabel('Age')
plt.ylabel('Frequency')
plt.show()
```

In the `plt.hist` function:
- `df['Age']` specifies the data to plot.
- `bins=10` defines how many bins (segments) the data should be divided into.
- `color` and `edgecolor` properties are used to customize the appearance of the histogram.

#### Here is what mine looked like:

![histogram](/Images/Week01/hist.png)

<br />
<br />

## Creating a DataFrame describing the length of sequences in our fasta file

#### Step 1: Make an empty dataframe
```python
#Make an empty dataframe
df = pd.DataFrame()

#Add columns
df["ID"] = [] #This is an empty list for data that you want to add
df["Length"] = [] #This is an empty list for data that you want to add
```
#### Step 2: Loop through our dictionary and add rows to the dataframe
```python
#Loop through the dictionary and add info about each sequence to the dataframe
for key in seq_dict.keys():
    id_temp = key
    seq_len_temp = len(seq_dict[key])
    
    #Use loc to add to the dataframe
    df.loc[len(df.index)] = [id_temp, seq_len_temp]
```

#### Step 3 (optional): Write the data in the dataframe to a text file (.csv)
```python
#Save the dataframe as a csv file
df.to_csv("<name of file to create>.csv")
```

<br />
<br />

<div style="border: 1px solid black; padding: 10px; margin: 10px 0;">
   <strong>Task 6: Plot a histogram of gene lengths</strong>
   <ol>
      <li>Use the code above to create a dataframe of sequence lengths and plot a histogram of the lengths.</li>
      <li>Do the same thing again BUT this time only include the genes in the chloroplast genome (ID starts with "ATC")</li>
      <li>Does it look like the distributions differ?</li>
   
   </ol>
</div>

<br />
<br />


