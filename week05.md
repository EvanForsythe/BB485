---
layout: default
---

<a name="top"></a>

# Week 4 lecture and tutorial
1. [Overview of genome annotation](#overview)
   - **A.** [What is genome annotation?](#what)
   - **B.** [Functional elements in the genome](#elements)
   - **C.** [Does "junk DNA" exist?](#junk)
2. [How annotation information is stored](#stored)
   - **A.** [GFF3 format](#gff3)
   - **B.** [BED format](#bed)
4. [Working with annotation data using base python](#python)
5. [Tools for working with annotation data with BioPython](#biopython)
6. [Web-based genome browsers](#web)
7. [Week 5 tutorial assignment](#tutorial)


## <ins>**Overview of genome annotation**<ins> <a name="overview"></a>
### <ins>**What is genome annotation?**<ins> <a name="what"></a>
### <ins>**Functional elements in the genome**<ins> <a name="elements"></a>
### <ins>**Does "junk DNA" exist?**<ins> <a name="junk"></a>
## <ins>**How annotation information is stored**<ins> <a name="stored"></a>
Storing genomic annotation data for a given genome assembly is akin to creating a road map to the location of a region of interest within the genome. Consistent with the idea of a map, genome annotations are stored as a set of <ins>**Genomic Coordinates**</ins>.

Genomic coordinates much always store four basic pieces of information (see below). This is the core information needed to reliably locate a specific region of the genome. In addition to this basic coordinate information, there is often (always) further information that describes things like what the feature is, where it came from, how it was predicted, and more.

<ins>**The four basic variables of genomic coordinates:**</ins>
- **Chromosome:** This represents the specific chromosome (e.g., chr3, chrY) where a feature is located. <ins>Note:</ins> this often refers to a contig, rather than a chromosome, depending on whether the assembly is of 'chromosome-level' vs 'contig-level' assembly quality. 
- **Start Position:** The beginning position of a feature along the chromosome in a genomic coordinate system. It is typically represented as a numerical value denoting the base pair position, with the first base pair being numbered 0 (or 1, depending on whether files are 0-indexed vs 1-indexed).
- **Stop Position:** The ending position of a feature along the chromosome in a genomic coordinate system. It represents the last base pair included in the feature, and it is typically represented as a numerical value.
- **Strand:** In molecular biology, DNA is composed of two complementary strands. The strand indicates the orientation of the DNA molecule and can be either positive (+) or negative (-). In genomic coordinate data files, it denotes which strand of the DNA the feature is located on.


### <ins>**GFF3 format**<ins> <a name="gff3"></a>

GFF3 files are tab-seperated value (.tsv) table in which each row represents a different 'feature' in the genome. The columns provide different information needed to access/sort/interpret the annotation data for the genomic feature.

<ins>**Column headers:**</ins>
- **Seqid:** the name of the sequence ID (e.g. chromosome number) that the feature is associated with.
- **Source:** the original database that the data came from (this type of information is called ‘metadata’)
- **Type:** the type of feature that the row is referencing
- **Start:** the position along the sequence where the feature starts
- **End:** the position along the sequence where the feature ends
- **Strand:** Whether the feature is shown on the strand contained in the reference sequence
- **Attribute:** This field can store misc. information about the feature (e.g. gene names, etc..)



### <ins>**BED format**<ins> <a name="bed"></a>

<ins>**Column headers:**</ins>
- **chrom:** The name of the chromosome (e.g. chr3, chrY) or scaffold.
- **chromStart:** The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
- **chromEnd:** The “non-inclusive” ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature, so it behaves much like ranges and substrings in python. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
- **name:** Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
- **score:** A score between 0 and 1000.
- **strand:** Defines the strand - either ’+’ or ’-’.
- **thickStart:** The starting position at which the feature is drawn thickly (for example, the start codon in gene displays).
- **thickEnd:** The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
- **itemRgb:** An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to ”On”, this RBG value will determine the display color of the data contained in this BED line.
- **blockCount:** The number of blocks (exons) in the BED line.
- **blockSizes:** A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount. Exons sizes are stored here if this is gene data.
- **blockStarts:** A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount. Exons start positions are stored.




## <ins>**Working with annotation data using base python**<ins> <a name="pythob"></a>
## <ins>**Tools for working with annotation data with BioPython**<ins> <a name="biopython"></a>
## <ins>**Web-based genome browsers**<ins> <a name="web"></a>
## <ins>**Week 5 tutorial assignment**<ins> <a name="tutorial"></a>


[Back to Top](#top)
