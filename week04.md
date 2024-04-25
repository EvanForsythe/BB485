---
layout: default
---

<a name="top"></a>

# Week 4 lecture and tutorial
1. [Overview of DNA sequencing](#DNAseq)
   - **A.** [Sanger sequencing](#sanger)
   - **B.** [Shotgun sequencing](#shotgun)
2. [Next-generation sequencing technologies](#next_gen)
   - **A.** [Short-read vs long-read sequencing](#versus)
3. [Genome assembly](#assembly)
4. [Tutorial assignment](#assign)
5. [More details about genome assembly](#details)
   - **A.** [De novo assembly versus reference-guided assembly](#versus)
   - **B.** [Paired-end sequencing](#paired)
   - **C.** [contigs and scaffolds](#scaff)
   - **D.** [Assembly quality statistics](#stats)
   - **E.** [fasta format versus fastq format](#formats)
6. [De novo genome assembly with SPAdes](#spades)
7. [Using a reference genome to assess assembly quality](#quality)
8. [Weekly project write-up assignment](#writeup)


## <ins>**Overview of DNA sequencing**<ins> <a name="DNAseq"></a>
DNA sequencing is a process used to determine the precise order of nucleotides within a DNA molecule. Nucleotides are the building blocks of DNA, and they consist of four types: adenine (A), thymine (T), cytosine (C), and guanine (G). By sequencing DNA, scientists can decipher the genetic information encoded within an organism's genome, including genes, regulatory sequences, and other genomic features. This information is crucial for understanding various biological processes, such as inheritance, evolution, and disease. DNA sequencing techniques have evolved significantly over the years, becoming faster, more accurate, and more cost-effective, thus enabling a wide range of applications in fields such as medicine, agriculture, forensics, and basic research.

- DNA sequencing technologies have been among the fastest growing technological sectors in the history of technology!!

![assemb11](/Images/Week04/assemb11.png)

- The ability to generate sequencing data is vastly outpacing computing technology to store/process the data.
   - Because of this, bioinformatics is an extremely important field!!!

### <ins>**Sanger sequencing**<ins> <a name="sanger"></a>

When researchers want to obtain the sequence of a single gene/region, they routinely apply a method known as 'cloning and sequencing'. 

![assemb01](/Images/Week04/assemb01.png)

Sequencing a single region of the genome often uses a sequenceing technology called **Sanger sequencing.** This method of sequencing was invented in 1977. It has been refined over the years but relies on the same elegant concept as the original invention.

![assemb04](/Images/Week04/assemb04.png)

Here is a video that nicely illustrates how Sanger sequencing works. 

[![Video 1 Thumbnail](https://img.youtube.com/vi/l0JVVPt4vNw/0.jpg)](https://www.youtube.com/watch?v=l0JVVPt4vNw)

### <ins>**Shotgun sequencing**<ins> <a name="shotgun"></a>

Sanger sequencing technology is not scalable to the scale needed to sequence a full genome. Full-genome sequencing has been made possible by 'next-generation' sequencing technologies. Next-gen sequencing relies on **shotgun sequencing**, which involves breaking long DNA molecules into smaller pieces.

![assemb02](/Images/Week04/assemb02.png)

## <ins>**Next-generation sequencing technologies**<ins> <a name="next_gen"></a>

There are several different DNA sequencing technologies available. Illumina is probably the most widely-used (for now), but there are strengths and weaknesses to each method (see below)

[![Video 2 Thumbnail](https://img.youtube.com/vi/CZeN-IgjYCo/0.jpg)](https://www.youtube.com/watch?v=CZeN-IgjYCo)

### <ins>**Short-read vs long-read sequencing**<ins> <a name="versus"></a>

![assemb03](/Images/Week04/assemb03.png)

## <ins>**Genome assembly**<ins> <a name="assembly"></a>
Below are a series of slides describing genome assembly from the [Bioinformatic.ca workshop](https://bioinformaticsdotca.github.io/BiCG_2019).



![assemb05](/Images/Week04/assemb05.png)
![assemb06](/Images/Week04/assemb06.png)
![assemb07](/Images/Week04/assemb07.png)
![assemb08](/Images/Week04/assemb08.png)
![assemb09](/Images/Week04/assemb09.png)
![assemb10](/Images/Week04/assemb10.png)

The algorithms used in genome assembly a very cool. However, since this course is focused on <ins>applying</ins> bioinformatics, we will not do a deep dive into the algorithms. If you're curious, below is a good video that walks through how de bruijn graphs and K-mers are used in genome assembly.

### De Bruijn graphs
[![Video Thumbnail 3](https://img.youtube.com/vi/TNYZZKrjCSk/0.jpg)](https://www.youtube.com/watch?v=TNYZZKrjCSk)


## <ins>**Tutorial assignment**<ins> <a name="assign"></a>
See canvas for this week's tutorial as part of that assingment, you'll need to install, Spades, an assembly program on the HPC using conda.

Below are example commands I used to setup my environment.
```bash
conda create -n assemble python=3.12
conda activate assemble
conda install bioconda::spades
```

## <ins>**More details about genome assembly**<ins> <a name="details"></a>

### <ins>**De novo assembly versus reference-guided assembly**<ins> <a name="versus"></a>

![ref](/Images/Week04/ref.png)

### <ins>**Paired-end sequencing**<ins> <a name="paired"></a>

![paired](/Images/Week04/paired.png)


[![Video Thumbnail](https://img.youtube.com/vi/WTbnk91e2WU/0.jpg)](https://www.youtube.com/watch?v=WTbnk91e2WU)

### <ins>**contigs and scaffolds**<ins> <a name="scaff"></a>

![scaff](/Images/Week04/scaffold.png)

### <ins>**Assembly quality statistics**<ins> <a name="stats"></a>

![quality](/Images/Week04/stats.png)

### <ins>**fasta format versus fastq format**<ins> <a name="formats"></a>
In addition to fasta files, you may also encounter 'fastq' files. Here is a description of the differences. 

### FASTA Files
A FASTA file contains biological sequences, such as DNA, RNA, or protein sequences. Each sequence entry consists of two parts:

1. **Header line**: Starts with a greater than symbol (`>`), followed by a unique identifier for the sequence and optional description.
2. **Sequence data**: A series of letters representing the nucleotide (A, T, G, C, U) or amino acid residues (A, R, N, etc.).

#### Example

```
>Sequence_1
ATCGATCGATCG
>Sequence_2
GCTAGCTAGCTA
```

### FASTQ Files
A FASTQ file also contains biological sequences, but with additional quality score information for each base in the sequence. Each sequence entry consists of four parts:

1. **Header line**: Same as FASTA, starts with a greater than symbol (`>`), followed by a unique identifier and optional description.
2. **Sequence data**: Same as FASTA, a series of letters representing nucleotide or amino acid residues.
3. **Separator line**: Starts with a plus symbol (`+`), optionally followed by the same identifier as the header line.
4. **Quality scores**: A series of ASCII characters representing the quality scores for each base in the sequence.

#### Example

```
@Sequence_1
ATCGATCGATCG
+Sequence_1
HHHHHHHHHHHH
@Sequence_2
GCTAGCTAGCTA
+Sequence_2
BBBBBBBBBBBB
```

### Key Differences

1. **Information**: FASTA files provide sequence data and optional annotations, while FASTQ files include sequence data along with quality scores.
2. **Format**: FASTA files are simpler and contain only sequences and headers, while FASTQ files are more complex, containing sequence data, headers, quality scores, and separator lines.
3. **Use Cases**: FASTA files are often used for storing and exchanging sequences, while FASTQ files are preferred for downstream analysis, particularly in next-generation sequencing (NGS) applications where quality scores are crucial for error correction and variant calling.



## <ins>**De novo genome assembly with SPAdes**<ins> <a name="spades"></a>

SPAdes is a popular genome assembler used for illumina sequencing data. [Here](https://github.com/ablab/spades) is the github page, which includes a quick start guide to using the software.

1. You should have already created a conda environment for doing assembly work and installed spades with: `conda install bioconda::spades`. Check that `spades.py` is in your PATH. 
2. Create a directory in which to work on this week's project.
3. Run a test run of spades directly from the command line with `spades.py --test`. This should run a test run of spades on a 'toy dataset' provided with spades. It will create an output folder called `spades_test`.
4. Once this is working we can create a job submission script to run a full assembly of paired-end illumina sequencing reads from Staphylococcus aureus (bacterium that cases Staph infections), which I downloaded from [this database](https://gage.cbcb.umd.edu/data/). Follow the instructions below to run spades inside of a job submission script.

Spades command (to go inside of a job submission script):
```bash
spades.py -1 /shared/forsythe/BB485/Week04/Staphylococcus_aureus_data/Illumina_reads/frag_1.fastq.gz -2 /shared/forsythe/BB485/Week04/Staphylococcus_aureus_data/Illumina_reads/frag_2.fastq.gz -o OUT/ -t 16
```
- `spades.py` is the program (which should be somewhere in your PATH)
- `-1` is the full path to the file containing the reads from **one of the paired ends**. Note that it is in fastq format and the file is compressed with gzip. Spades can work directly from the zipped version of the file, so there's no need to unzip it.
- `-2` is the full path to the file containing the reads from the **other paired end**. Note that it is in fastq format and the file is compressed with gzip. Spades can work directly from the zipped version of the file, so there's no need to unzip it.
- `-o` is used to tell spades where to put the output files. You don't need to make this folder in advance. Spades will do it automatically.
- `-t` tells spades how many thread (aka cores, aka CPUs) to use. **Be sure to set this to the same number near the top of your job submission script**.

### ALTERNATIVELY:
To test your spades program:
```bash
/shared/forsythe/BB485/Week04/Spades/SPAdes-3.15.5-Linux/bin/spades.py --test
```

```bash
/shared/forsythe/BB485/Week04/Spades/SPAdes-3.15.5-Linux/bin/spades.py -1 /shared/forsythe/BB485/Week04/Staphylococcus_aureus_data/Illumina_reads/frag_1.fastq.gz -2 /shared/forsythe/BB485/Week04/Staphylococcus_aureus_data/Illumina_reads/frag_2.fastq.gz -o OUT/ -t 16
```

### SPAdes output:
- the most important output file that is generated will be called `contigs.fasta`. You can count the number of contigs your assembly achieved with `grep ">" contigs.fasta`.

## <ins>**Using a reference genome to assess assembly quality**<ins> <a name="quality"></a>

Because of the importance of Staph infection in medicine, large efforts have been made by reserchers to assemble of 'finished' version of the genome, in which all the contigs have been scaffolded together such that number of contigs is equal to the number of chromosomes.

I downloaded a finished version of the genome to use as a reference genome to assess the quality of our assembly. The file is located at:
```bash
/shared/forsythe/BB485/Week04/Staphylococcus_aureus_data/Reference_genome/GCA_037039335.1_ASM3703933v1_genomic.fa
```
We will use a program called QUAST (Quality Assessment Tool for Genome Assemblies) to compare our assembly to the reference genome to assess how well our assembly worked.

- **Note:** I had trouble installing QUAST with conda, so instead we'll need to go a different route to use the QUAST program:
   - We are going to run quast by providing the full path to where a copy of the program lives on the HPC: ```/shared/forsythe/BB485/Week04/Quast/quast-5.2.0/quast.py```

 To run QUAST use (you don't need to run it in a job because it's fast):
 ```bash
/shared/forsythe/BB485/Week04/Quast/quast-5.2.0/quast.py <path/to/your/assembled_genome.fasta> -r /shared/forsythe/BB485/Week04/Staphylococcus_aureus_data/Reference_genome/GCA_037039335.1_ASM3703933v1_genomic.fa -o <path/to/your/out/dir/>
```
### QUAST OUTPUT:
- All the relevant info you'll need is in: `report.txt`

## <ins>**Weekly project write-up assignment**<ins> <a name="writeup"></a>
- Perform spades genome assembly with and without the `--careful` flag.
   - This means that you'll need to run spades two seperate times (as two seperate job submissions).
- Use QUAST to generate assembly statistics.
- Use the info contained in report.txt to fill in the following table (create it in excel, etc...):

| Method               | Number of Contigs | Length of Largest Contig | N50   | Genome Fraction (%) |
|----------------------|-------------------|--------------------------|-------|---------------------|
| Without Careful Argument |                   |                          |       |                     |
| With Careful Argument    |                   |                          |       |                     |


- Answer the following questions:
   - TBD

[Back to Top](#top)

