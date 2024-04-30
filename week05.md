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

## <ins>**Working with annotation data using base python**<ins> <a name="pythob"></a>
## <ins>**Tools for working with annotation data with BioPython**<ins> <a name="biopython"></a>
## <ins>**Web-based genome browsers**<ins> <a name="web"></a>
## <ins>**Week 5 tutorial assignment**<ins> <a name="tutorial"></a>


[Back to Top](#top)
