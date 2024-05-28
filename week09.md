---
layout: default
---

<a name="top"></a>


# Week 9 lecture and tutorial
1. [Introduction to Squeakuences](#intro)
2. [Squeakuences bioinformatics tool development](#tool)
3. [Testing and benchmarking Squeakuences](#test)
5. [Tutorial assignment](#tut)
    

## <ins>**Introduction to Squeakuences**</ins> <a name="intro"></a>
Linnea is a recent OSU-Cascades Biology/Computer Science graduatae and is the lead developer on Squeakuences. Below is an introduction to the tool.

![squeak01](/Images/Week09/squeak01.png)

![squeak01](/Images/Week09/squeak02.png)

![squeak01](/Images/Week09/squeak03.png)

![squeak01](/Images/Week09/squeak04.png)

![squeak01](/Images/Week09/squeak05.png)

![squeak01](/Images/Week09/squeak06.png)

![squeak01](/Images/Week09/squeak07.png)

## <ins>**Squeakuences bioinformatics tool development**</ins> <a name="tool"></a>

### Planning a manuscript

**Draft Abstract**

## <ins>**Testing and benchmarking Squeakuences**</ins> <a name="test"></a>

### Obtaining squeakuences

We will each use git/github to 'clone' squeakuences to the HPC so that we can run it.

### Obtaining test datasets

Squeakuences is designed to work with **fasta files**. Fasta files are used to store a wide variety of biological sequence information. We would like Squeakuences to work on any fasta file from any source. In order to identify unexpected problems/bugs caused by specific data types/sources, we need to test it as widely as possible. 

- NCBI genome sequences: `https://www.ncbi.nlm.nih.gov/genome/`
- Phytozome plant genome/proteome database: `https://phytozome-next.jgi.doe.gov/`
- NCBI raw sequencing data, Sequence Read Archive (SRA): `https://www.ncbi.nlm.nih.gov/sra/`

Datasets of specific interest:
- Parasitoid wasps:
    - `https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=63433`
    - `https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=7425`
- Devaleraea mollis seaweed sequencing:
    - Genome assembly: `https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_032361265.1/`
    - Raw sequencing reads: `https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR25035936&display=metadata`
- Seaweed genome sequencing:
    - Genome assembly: `https://figshare.com/s/517f180deb56a5b31b14`
    - Raw sequencing reads: `https://www.ncbi.nlm.nih.gov/nuccore?term=PRJNA597239`

## <ins>**Tutorial assignment**</ins> <a name="tool"></a>

### Benchmarking table

Download three genomic datasets and apply sqeakuences to each. Fill in the following information.

|             Stat         | Dataset 1 | Dataset 2 | Dataset 3 |
|--------------------------|---------|---------|---------|
| Taxon name        |         |         |         |
| Data Source (URL)        |         |         |         |
| Data type (prot/RNA/DNA)        |         |         |         |
| Processing Time (hours) |         |         |         |
| Memory (MB)              |         |         |         |
| Starting File Size (MB)  |         |         |         |
| Ending File Size (MB)    |         |         |         |
| Number of sequences cleaned |         |         |         |
| Example seq ID (before)  |         |         |         |
| Example seq ID (after)   |         |         |         |





[Back to Top](#top)
