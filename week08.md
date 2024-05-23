---
layout: default
---

<a name="top"></a>

# Week 8 lecture and tutorial
1. [The genome and the transcriptome](#transcriptome)
2. [Layers of the transcriptome](#layers)
3. [Gene expression](#gene)
4. [Tools for studying gene expression](#study)
	- **A.** [Microarrays](#micro)
	- **B.** [RNAseq](#RNAseq)
5. [RNAseq workflow](#workflow)
	- **A.** [Library prep](#library)
	- **B.** [In silico analysis](#silico)
6. [Downstream analyses of RNAseq data](#down)
	- **A.** [Differential expression](#diff)
	- **B.** [Gene co-expression networks](#net)
7. [Tutorial assignment](#tut)
8. [Performing a transcriptome analysis with hisat2 on example data](#hisat)
9. [Analyzing a full genome dataset](#full)




## <ins>**The genome and the transcriptome**</ins> <a name="transcriptome"></a>
- Gen**ome** = all the genes in a cell/tissue/organism
- Transcrip**ome** = all the transcripts (RNAs) in a cell/tissue/organism
![trans01](/Images/Week08/trans01.png)

The transcriptome is highly dynamic. The ability to <ins>profile</ins> the transcriptome gives us extremely important clues about the cellular function of each gene.

## <ins>**Layers of the transcriptome**</ins> <a name="layers"></a>
The full transcriptome includes all types of RNAs present in a cell. However, in many cases, the term 'transcriptome' tends to refer to the full set of mRNAs (protein-coding transcripts) in the cell.

![trans02](/Images/Week08/trans02.png)

## <ins>**Gene expression**</ins> <a name="gene"></a>

The transcriptome essentially gives us a readout of gene expression for all the genes in a genome in a given cell/tissue/organism. 

- The transcriptome differs across developmental stages. A different set of gene is expressed in the embryo than in an adult.
	- Identifying which genes are expressed when gives us a headstart on identifying genes involved in development.

![trans03](/Images/Week08/trans03.png)

Observing differences in the transcriptome is akin to observing differences in gene regulation.

![trans04](/Images/Week08/trans04.png)

Gene expression also differs between cell types and tissues in an organism.
- The image below is a readout from a plant genome/transcriptome database called [The Arabidopsis Information Resource (TAIR)](https://www.arabidopsis.org/)
	- Try searching on gene on TAIR:
 		- Search `TIC110` (or any gene of interest)
   		- Scroll down to `Expression`
		- Set Data Source to `Developmental_map` 

![trans05](/Images/Week08/trans05.png)

Gene expression also changes in response to environmental stimulants.

![trans06](/Images/Week08/trans06.png)

## <ins>**Tools for studying gene expression**</ins> <a name="study"></a>

Advances in tools/technology used to study transcriptomes have revolutionized molecular biology/medical research.

Below is a 'heat-map' showing transcriptome expression levels inferred using two different methods.

![trans10](/Images/Week08/trans10.png)

### <ins>**Microarrays**</ins> <a name="micro"></a>

![trans07](/Images/Week08/trans07.png)

![trans08](/Images/Week08/trans08.png)

![trans09](/Images/Week08/trans09.png)

Because of the drawbacks to microarrays listed above, RNA seq has become more popular in recent years.

![trans11](/Images/Week08/trans11.png)

### <ins>**RNAseq**</ins> <a name="RNAseq"></a>

![trans12](/Images/Week08/trans12.png)

## <ins>**RNAseq workflow**</ins> <a name="workflow"></a>

![trans13](/Images/Week08/trans13.png)

### <ins>**Library prep**</ins> <a name="library"></a>

Library prep involves synthesizing cDNA (complementary DNA)

![trans14](/Images/Week08/trans14.png)

### <ins>**In silico analysis**</ins> <a name="silico"></a>

![trans15](/Images/Week08/trans15.png)

## <ins>**Downstream analyses of RNAseq data**</ins> <a name="down"></a>

![trans18](/Images/Week08/trans18.png)

![trans19](/Images/Week08/trans19.png)

### <ins>**Differential expression**</ins> <a name="diff"></a>

Differential expression analyses can point to genes that are upregulated or downregulated in cells/tissues/environments/etc

![trans16](/Images/Week08/trans16.png)

### <ins>**Gene co-expression networks**</ins> <a name="net"></a>

Co-expression networks identify pairs of genes that show correlated expression across a variety of cells/tissues/environments/etc. Co-expression is a good indicator that two proteins interact with eachother.

![trans17](/Images/Week08/trans17.png)

## <ins>**Tutorial assignment**</ins> <a name="tut"></a>

For the tutorial assignment, setup a new conda environment, which we'll use to perform RNA seq data on Thursday

Install the following packages:
## NEW VERSION:
- `conda install bioconda::hisat2`
- `conda install bioconda::stringtie`
- `conda install bioconda::samtools`


## OLD VERSION (dont install these afterall)
- tophat2
- cufflinks
- DESeq2

## <ins>**Performing a transcriptome analysis with hisat2**</ins> <a name="hisat"></a>
hisat2 (hierarchical indexing for spliced alignment of transcripts) is an algorithm for aligning whole-transcriptome sequencing reads to a reference genome.

## Installing hisat2
- You should already have hisat2 installed in your most recent conda environment. Check that it is installed by running `hisat2`. You should see the abreviated user manual.
- [here](https://daehwankimlab.github.io/hisat2/manual/) is the full-length user manual.




SAM (Sequence Alignment/Map) format is a text-based format for storing sequence alignments against a reference genome.



## Structure of a SAM File

A SAM file has two main sections:
1. **Header Section**: Meta-information about the alignments.
2. **Alignment Section**: Contains read alignment information.

### Header Section
- Lines start with `@`
- Common headers:
  - `@HD`: Header
  - `@SQ`: Reference sequences
  - `@RG`: Read groups
  - `@PG`: Programs

### Alignment Section
- One line per read alignment
- 11 mandatory fields:
  1. `QNAME`: Query name
  2. `FLAG`: Bitwise flag
  3. `RNAME`: Reference name
  4. `POS`: Position
  5. `MAPQ`: Mapping quality
  6. `CIGAR`: CIGAR string
  7. `RNEXT`: Mate reference name
  8. `PNEXT`: Mate position
  9. `TLEN`: Template length
  10. `SEQ`: Sequence
  11. `QUAL`: Quality

Optional fields (e.g., `NM`, `MD`, `AS`, `XS`)

## <ins>**Analyzing a full genome dataset**</ins> <a name="full"></a>

Full datasets from a RNAseq experiment on a human sample are located at this path:
```bash
/shared/forsythe/BB485/Week08/hg38/
```
- `GCA_000001405.28_GRCh38.p13_genomic.fna` is the genome assembly in fasta format
- `SRR28621297.fasta` is the RNA seq reads file
- `GCA_000001405.28_GRCh38.p13_genomic.gff` is the gene annotation file




[Back to Top](#top)
