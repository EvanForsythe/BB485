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
   - **A.** [de novo assembly](#de_novo)
   - **B.** [reference-guided assembly](#ref)
4. [Tutorial assignment](#assign)


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


[Back to Top](#top)

