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
Storing genomic annotation data for a given genome assembly is akin to creating a 'road map' to the location of a region of interest within the genome. Consistent with the idea of a 'map', genome annotations are stored as a set of <ins>**Genomic Coordinates**</ins>.

Genomic coordinates much always store four basic pieces of information (see below). This is the core information needed to reliably locate a specific region of the genome. In addition to this basic coordinate information, there is often (always) further information that describes things like: what the feature is, where it came from, how it was predicted, and more.

<ins>**The four basic variables of genomic coordinates:**</ins>
- **Chromosome:** This represents the specific chromosome (e.g., chr3, chrY) where a feature is located. <ins>Note:</ins> this often refers to a contig, rather than a chromosome, depending on whether the assembly is of 'chromosome-level' vs 'contig-level' assembly quality. 
- **Start Position:** The beginning position of a feature along the chromosome in a genomic coordinate system. It is typically represented as a numerical value denoting the base pair position, with the first base pair being numbered 0 (or 1, depending on whether files are 0-indexed vs 1-indexed).
- **Stop Position:** The ending position of a feature along the chromosome in a genomic coordinate system. It represents the last base pair included in the feature, and it is typically represented as a numerical value.
- **Strand:** In molecular biology, DNA is composed of two complementary strands. The strand indicates the orientation of the DNA molecule and can be either positive (+) or negative (-). In genomic coordinate data files, it denotes which strand of the DNA the feature is located on.

###0-based vs 1-baded indexing in common annotation file formats

| Format            | Position system          |
|-------------------|--------------------------|
| GFF/GTF           | 1-based                  |
| BLAST results     | 1-based                  |
| BLAT results      | 1-based                  |
| maf               | 0-based                  |
| bam               | 0-based                  |
| sam               | 1-based                  |
| BED               | 0-based start, 1-based end |



### <ins>**GFF3 format**<ins> <a name="gff3"></a>

GFF3 files are tab-seperated value (.tsv) table in which each row represents a different 'feature' in the genome. The columns provide different information needed to access/sort/interpret the annotation data for the genomic feature.

<ins>**Column headers:**</ins>
| Field    | Description                                                                                           | Example |
|----------|-------------------------------------------------------------------------------------------------------|---------|
| Seqid    | The name of the sequence ID (e.g. chromosome number) that the feature is associated with.            | Chr5 |
| Source   | The original database that the data came from (this type of information is called ‘metadata’)         |phytozomev12|
| Type     | The type of feature that the row is referencing                                                       | CDS |
| Start    | The position along the sequence where the feature starts                                             | 995 |
| End      | The position along the sequence where the feature ends                                               | 5156|
| Strand   | Whether the feature is shown on the strand contained in the reference sequence                        | + |
| Attribute| This field can store misc. information about the feature (e.g. gene names, etc..)                   | ID=AT5G01010.2.Araport11.447;Name=AT5G01010.2;pacid=37423170;longest=1;Parent=AT5G01010.Araport11.447
![image](https://github.com/EvanForsythe/BB485/assets/40008668/d6b76d43-b61f-4ec2-9471-2782bf09c247)
 |


### <ins>**GTF format**<ins> <a name="gtf"></a>
A GTF file is very similar to a GFF file, but with a few different specifications. The first eight GTF fields are the same as GFF. The major difference is that, in GTF format, attribute field has been expanded to allow for a list of information within a single column. Each attribute consists of a type/value pair. Attributes must end in a semi-colon, and be separated from any following attribute by exactly one space. The attribute list must begin with the two mandatory attributes:
1. gene id value - A globally unique identifier for the genomic source of the sequence.
2. transcript id value - A globally unique identifier for the predicted transcript.

| Variable   | Description                                                                                         | Example            | Notes   |
|------------|-----------------------------------------------------------------------------------------------------|--------------------|---------|
| seqname    | The name of the sequence. Must be a chromosome or scaffold.                                        | chr1, chrX, scaffold123 |         |
| source     | The program or database that generated or annotated this feature.                                    | Ensembl, GENCODE, Gnomon |         |
| feature    | The type of feature, such as "exon", "CDS", "start_codon", or "stop_codon".                         | exon, CDS, start_codon |         |
| start      | The starting position of the feature in the sequence. The first base is numbered 1.                | 1001, 5000, 23789 |         |
| end        | The ending position of the feature in the sequence, inclusive.                                       | 2000, 6000, 24567 |         |
| score      | A score associated with the feature, typically ranging from 0 to 1000. If not applicable, a dot "." is used. | 100, 500, . | Optional |
| strand     | The strand of the feature, represented as "+" (forward) or "-" (reverse). Unknown or undefined strands may be represented as ".". | +, -, . | Optional |
| frame      | For features that represent coding sequences (CDS), the reading frame of the first base. It's a number between 0 and 2, representing the remainder of the division of the feature's start position by 3. For non-coding features, it's typically ".". | 0, 1, 2, . | Optional |
| attributes | Additional information about the feature, often stored as key-value pairs. These can include gene IDs, transcript IDs, gene names, and other relevant information. | gene_id "ENSG00000112345"; transcript_id "ENST00000123456"; gene_name "TP53" | Optional |


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

### <ins>**SAM/BAM format**<ins> <a name="sam"></a>

| Field | Type   | Brief description                     | Example      |
|-------|--------|--------------------------------------|--------------|
| QNAME | String | Query template NAME                  | "query1"     |
| FLAG  | Int    | bitwise FLAG                         | 16           |
| RNAME | String | Reference sequence NAME              | "chr1"       |
| POS   | Int    | 1-based leftmost mapping POSition   | 1001         |
| MAPQ  | Int    | MAPping Quality                      | 30           |
| CIGAR | String | CIGAR string                         | "50M"        |
| RNEXT | String | Ref. name of the mate read           | "chr2"       |
| PNEXT | Int    | Position of the mate/next read       | 2000         |
| TLEN  | Int    | observed Template LENgth             | 300          |
| SEQ   | String | segment SEQuence                     | "ATCGGTA"    |
| QUAL  | String | ASCII of Phred-scaled base QUALity+33| "!&#$^@&"   |


## <ins>**Working with annotation data using base python**<ins> <a name="pythob"></a>
## <ins>**Tools for working with annotation data with BioPython**<ins> <a name="biopython"></a>
## <ins>**Web-based genome browsers**<ins> <a name="web"></a>
## <ins>**Week 5 tutorial assignment**<ins> <a name="tutorial"></a>


[Back to Top](#top)
