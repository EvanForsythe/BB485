---
layout: default
---

## Course structure

![Points structure of BB 485](/Images/Week01/points_piechart.png)

## What is bioinformatics?

![what is bioinformatics](/Images/Week01/bioinformatics.png)

## Connecting to Jupyter Hub course workspace

The link to the lab Jupyter Hub is on Canvas under Class Resources.

When you select "Start My Server" you'll get a choice of different server size options. Choose "Standard" unless instructed otherwise.

![Jupyter Hub map](/Images/Week01/JH_nav.png)

## Linux command line navigation

![Anatomy of a unix command](/Images/Week01/unix_command.png)

# Introduction to Basic Linux/Unix Commands

## A few example commands

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

## Data visualization in python

## Overview of programs/servers/software/databases/repositories we'll use this term


Here is an exmaple of a code block in the tutorial. 

```python{linenos=True}
def hello_world():
    print("Hello, world!")
```

