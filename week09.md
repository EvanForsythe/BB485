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

**Journal submission guidelines:**

- (Bioinformatics journal submission guidelines)[https://academic.oup.com/bioinformatics/pages/instructions_for_authors]

**Draft Abstract**

## <ins>**Testing and benchmarking Squeakuences**</ins> <a name="test"></a>

### Obtaining squeakuences

We will each use git/github to 'clone' squeakuences to the HPC so that we can run it. Go to the [squeakuences github page](https://github.com/EvanForsythe/Squeakuences) and then use `git clone` to clone it to the HPC.

### Obtaining test datasets

Squeakuences is designed to process **fasta files**. Fasta files are used to store a wide variety of biological sequence information. We would like Squeakuences to work on any fasta file from any source. In order to identify unexpected problems/bugs caused by specific data types/sources, we need to test it as widely as possible. 

- NCBI genome sequences: `https://www.ncbi.nlm.nih.gov/genome/`
- Phytozome plant genome/proteome database (account login required): `https://phytozome-next.jgi.doe.gov/`
- NCBI raw sequencing data, Sequence Read Archive (SRA): `https://www.ncbi.nlm.nih.gov/sra/`

Datasets of specific interest:
- Parasitoid wasps:
    - `https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=63433`
    - `https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=7425`
- Devaleraea mollis seaweed sequencing:
    - Genome assembly: `https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_032361265.1/`
    - Raw sequencing reads: `https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR25035936&display=metadata`
- Sargassum fusiforme seaweed genome sequencing:
    - Genome assembly: `https://figshare.com/s/517f180deb56a5b31b14`
    - Raw sequencing reads: `https://www.ncbi.nlm.nih.gov/nuccore?term=PRJNA597239`

To download a dataset directly to the HPC, you can run the following command from the HPC:
```bash
wget <URL of download link>
```
- Note that this doesn't work on every database. Unfortunately, some databases require that we download to our local machine and the upload to the HPC (using the `scp` command). 

## <ins>**Tutorial assignment**</ins> <a name="tut"></a>

### Benchmarking table

1. Download three genomic datasets and apply sqeakuences to each. Fill in the following information.
- Try to include three different database sources or data types.

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

2. Create an 'Issue' on github suggesting a way to improve Squeakuences. If you experience any bugs, you'll definitely want to report those. If you don't see any bugs, use the Issue to suggest an improvement in the usability of Squeakuences. 

[Back to Top](#top)
