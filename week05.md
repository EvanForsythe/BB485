---
layout: default
---

<a name="top"></a>

# Week 5 lecture and tutorial
1. [Overview of genome annotation](#overview)
2. [How annotation information is stored](#stored)
   - **A.** [GFF3 format](#gff3)
   - **B.** [GTF format](#gtf)
   - **C.** [BED format](#bed)
   - **D.** [SAM/BAM format](#sam)
3. [Working with annotation data in python](#python)
4. [Week 5 tutorial assignment](#tutorial)
5. [Performing genome annotation of a prokaryotic genome](#annot)
   - **A.** [Using Prokka annotation software](#prokka)
   - **B.** [Visualizing annotations with Circos plots](#circos)
   - **C.** [Using python to extract sequences of interest](#seqs)
6. [Week 5 write-up assignment](#writeup)

## <ins>**Overview of genome annotation**<ins> <a name="overview"></a>

![annot01](/Images/Week05/annot01.png)
![annot02](/Images/Week05/annot02.png)
![annot03](/Images/Week05/annot03.png)

## <ins>**How annotation information is stored**<ins> <a name="stored"></a>
Storing genomic annotation data for a given genome assembly is akin to creating a 'road map' to the location of a region of interest within the genome. Consistent with the idea of a 'map', genome annotations are stored as a set of <ins>**Genomic Coordinates**</ins>.

Genomic coordinates, at the most basic level, will always store four basic pieces of information (see below). This is the core information needed to reliably locate a specific region of the genome. In addition to this basic coordinate information, there is often further information that describes things like: what the feature is, where it came from, how it was predicted, and more.

<ins>**The four basic variables of genomic coordinates:**</ins>
- **Chromosome:** This represents the specific chromosome (e.g., chr3, chrY) where a feature is located. <ins>Note:</ins> this often refers to a contig, rather than a chromosome, depending on whether the assembly is of 'chromosome-level' vs 'contig-level' assembly quality. 
- **Start Position:** The beginning position of a feature along the chromosome in a genomic coordinate system. It is typically represented as a numerical value denoting the base pair position, with the first base pair being numbered 0 (or 1, depending on whether files are 0-indexed vs 1-indexed).
- **Stop Position:** The ending position of a feature along the chromosome in a genomic coordinate system. It represents the last base pair included in the feature, and it is typically represented as a numerical value.
- **Strand:** In molecular biology, DNA is composed of two complementary strands. The strand indicates the orientation of the DNA molecule and can be either positive (+) or negative (-). In genomic coordinate data files, it denotes which strand of the DNA the feature is located on.

### 0-based vs 1-based indexing in common annotation file formats

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

![annot04](/Images/Week05/annot04.png)
![annot05](/Images/Week05/annot05.png)
![annot06](/Images/Week05/annot06.png)


GFF3 files are tab-seperated value (.tsv) table in which each row represents a different 'feature' in the genome. The columns provide different information needed to access/sort/interpret the annotation data for the genomic feature.

| Field    | Description                                                                                           | Example |
|----------|-------------------------------------------------------------------------------------------------------|---------|
| Seqid    | The name of the sequence ID (e.g. chromosome number) that the feature is associated with.            | Chr5 |
| Source   | The original database that the data came from (this type of information is called ‘metadata’)         |phytozomev12|
| Type     | The type of feature that the row is referencing                                                       | CDS |
| Start    | The position along the sequence where the feature starts                                             | 995 |
| End      | The position along the sequence where the feature ends                                               | 5156|
| Strand   | Whether the feature is shown on the strand contained in the reference sequence                        | + |
| Attribute| This field can store misc. information about the feature (e.g. gene names, etc..)                   | ID=AT5G01010.2.Araport11.447; Name=AT5G01010.2; pacid=37423170; longest=1; Parent=AT5G01010.Araport11.447|


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

BED files are often used to store annotation data used in genome browsers, like the [UCSC Genome Browser](https://genome.ucsc.edu/)


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

SAM (Sequence Alignment/Map) and BAM (Binary Alignment/Map) files are commonly used in genome annotation and represent formats for storing sequence alignment data, such as reads aligned to a reference genome. These are not considered genome annotation per se becuase they describe a sequencing experiment, which is different from the information we've defined as an 'annotation'. However, information about how sequencing reads map to a genome is often used as a very important layer of functional annotation. We will discuss this type of quantitative read mapping when we discuss **<ins>Transcriptomics and RNA sequencing</ins>** in a couple weeks.


Here is basic information contained in SAM/BAM files:

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


## <ins>**Working with annotation data using base python**<ins> <a name="python"></a>
Genome annotation data are stored in table format. Python has excellent tools forworking with table-shaped data by storing it as a pandas DataFrame. Let's practice creating a python script for reading in an annotation file and using it to extract important genomic information.

```python3
#Import needed modules
import pandas as pd

# Create a string object that is a full path to a tsv file
file_path = <full path to the tsv file>

# Read the tsv file in as a dataframe
annot_df = pd.read_csv(file_path, delimiter= "\t")

# Print the dataframe to get a look at it
print(annot_df)
```

To practice extracting information from a dataframe, loop through the rows in the dataframe and print the start and end coordinates.
```python3
#Loop through the rows of the dataframe
for i in annot_df.index:
    
    #Get the type of feature
    temp_type = annot_df.iloc[i]['type']
    
    #Get the start position
    temp_start = annot_df.iloc[i]['start']

    #Get the stop position
    temp_stop = annot_df.iloc[i]['end']

    #Create a print statement
    print(f"Type: {temp_type}, Start: {temp_start}, End: {temp_stop}")

```

Remember: an annotation file can't tell us anything without the corresponding genome assembly file to which it pertains. If we want to extract sequence information, we need to read the genome assembly file into python as well.

```python3
# Read in the DNA sequence associated with the annotations

#Get the full path to the DNA sequence
seq_file_path = "/shared/forsythe/BB485/Week05/A_thaliana_chr5_short.fa"

seq_file_handle = open(seq_file_path, "r")

#Create an empty dictionary
seq_dict = {}

#Loop through the line in the file
for line in seq_file_handle:
    if line.startswith(">"):
        id_temp = line.strip() #Removes "\n"
        id_clean = id_temp.replace(">", "") #Removes ">" by replacing with nothing.
        
        #Add the item to the dictionary
        seq_dict[id_clean]="" # id_clean is the key, the value is an empty string (for now)
    else:
        seq_line = line.strip() #Removes "\n"
        
        #append this line to the dictionary value, using the key (which is still "id_clean" from the previous line)
        seq_dict[id_clean] += seq_line

#This is printing the first 100 nucleotides of chromosome 5. You can alter this to print from any start to any end coordinate.
print(seq_dict["Chr5"][0:100])


```
## <ins>**Week 5 tutorial assignment**<ins> <a name="tutorial"></a>
- Modify the python code we've been working on to:
   - Loop through genomic features
   - IF the feature is of type "gene"
      - print the DNA sequence from start to end
   - Paste your final python code and the first lines of output into the canvas assignment.

- **HINTS:**
   - Perform the following steps in order (these should each correspond to a block of code from the tutorial, but you may need to switch the order compared to how they're shown on the tutorial):
      - Read in the tsv file and store as dataframe
      - Read in the CSV file and store as dictionary
      - Loop through the tsv file (add an if statement inside of the for loop and print the DNA sequence from inside the if statement).
         - Here's what my print statement looked like: `print(seq_dict["Chr5"][temp_start:temp_stop])` 


## <ins>**Performing genome annotation of a prokaryotic genome**<ins> <a name="annot"></a>

Annotating genomes is a complex process that can make use of many different types of information. We can annotate genomes 'de novo' by searching for patterns/hallmarks of the genomic features we're expecting to find. For example, we can use our knowledge of the genetic code to create a (relatively) simple bioinformatics program that scans DNA sequences to identify <ins>**Open Reading Frames (ORFs)**</ins>.

![orf](/Images/Week05/orf.png)

In addition to the 'de novo' approach to genome annotation, it is very powerful to use existing annotations from a previously annotated model organism (e.g. E. coli) as a starting point and search for similar features in an unannoted genome of interest. This type of annotation is referred to as <ins>**homology-based genome annotation**</ins>.

### <ins>**Using Prokka annotation software**<ins> <a name="prokka"></a>

[Prokka](https://github.com/tseemann/prokka) is a software package for prokaryotic genome annotation. It performs homology-based genome annotation by comparing our unannotated genome of interest (e.g. the contigs from a genome assembly) to curated databases of annotations from other prokaryotes. 

![prokka](/Images/Week05/prokka.png)

To run prokka, we need to install it. We will try two different options. First, our typical method:
```bash
conda install bioconda::prokka
```

If this doesn't work, we'll use the module load system:
```bash
module load prokka
```

Once installed, you can run prokka from <ins>**within a job submission script**</ins> by adding the following lines to your script.
```bash
#load prokka module
module load prokka

#run prokka on assembly dataset
prokka <full/path/to/assembly file>
```
- Note that your assembly file from last week should be named `contigs.fasta`. You can use either of the assemblies that you generated.
- Note: you should only need 8 CPUs to run this job and it should finish in a couple minutes.


This will output several files in a folder that starts with: `PROKKA_`
- The file ending in `.tbl` is a descriptive table that tells us about each annotation
- The file ending in `.txt` gives a summary of what annotations you found in the genome
- The file ending in `.gff` is a gff format file of the annotations

<ins>**IMPORTANT NOTE:**</ins> If you couldn't get prokka to run, you can use the output of the prokka run I did. These files can be found at:
`/shared/forsythe/BB485/Week05/Staphylococcus_aureus_assembly/PROKKA_05022024/`

### <ins>**Visualizing annotations with Circos plots**<ins> <a name="circos"></a>

![circos](/Images/Week05/circos.png)

### <ins>**Using python to extract sequences of interest**<ins> <a name="seqs"></a>
We can modify our python script from earlier in the week, to print the sequence of any genomic feature.

Here is a block of python code you can use to read in a gff file:
```python
import pandas as pd

# Define the columns of the GFF3 file
columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

# Add the string for the full path to your annotation file from prokka
gff_path = <path/to/file>

# Read the GFF3 file into a DataFrame
df_temp = pd.read_csv(gff_path, sep="\t", comment="#", header=None, names=columns)

df = df_temp.dropna(subset=['attributes'])

# Display the DataFrame
print(df.head())
```

For the assignment, you'll be instructed to write about a random feature in the genome. Once you've read in your gff table as a dataframe, you can show just one row (at random) from it with:

```python
#Get one random row
random_row = df.sample(n=1)
#Print that row
print(random_row)
```

You can then extract the values from the random row with:

```python
# Extracting values from specific rows
seqid_value = random_row.at[random_row.index[0], 'seqid']
type_value = random_row.at[random_row.index[0], 'type']
start_value = random_row.at[random_row.index[0], 'start']
end_value = random_row.at[random_row.index[0], 'end']
attributes_value = random_row.at[random_row.index[0], 'attributes']

print("seqid:", seqid_value)
print("type:", type_value)
print("start:", start_value)
print("end:", end_value)
print("attributes:", attributes_value)
```

Once you have the above information, you can use you python code from earlier this week to extract a DNA sequence corresponding to your random feature.

For the assignment below, combine blocks of python code to create a script that does the following:
- Read in the annotation gff file
- Get a random annotation from the larger dataframe
- Extract coordinates and other relevant info for that annotation
- Read the a fasta file in as a dictionary:
   - `/shared/forsythe/BB485/Week05/Staphylococcus_aureus_assembly/filtered_contigs.fasta`
- Once you have the fastas file stored as a dictionary object, you can extract the DNA sequence corresponding to your feature by using the `start_value` and  `end_value` variables. Save this DNA string as a new variable `my_dna_seq` (or something similar).
- Write a python expression that calculates the GC content of your sequence. In plain english, the way to calculate GC content is:
   - GC content = (number of G's + number of C's) divided by the length of the DNA.
- HINT: debug all of this and then once it's running, run one final run of the script to generate 'your' random gene to use to fill in the second table in the assignment. 


## <ins>**Week 5 write-up assignment**<ins> <a name="writeup"></a>
1. Perform genome annotation with prokka on the genome assembly from last week's assignment. The file should be named `contigs.fasta`

2. Fill in the following table with your annotation results (you don't need python for this):

   | Annot type | Count    |
   |------------|----------|
   | contigs    |          |
   | bases      |          |
   |   CDS      |          |
   |    tRNA    |          |

3. Fill in the following table about a random feature from your annotation file (use python to extract info from the larger gff and fasta files):

   | Info                  |   Value      |
   |-----------------------|--------------|
   | Gene description      |              |
   | Gene location (start) |              |
   | Gene location (end)   |              |
   | Gene length           |              |
   | GC content            |              |

   -Note that the gene description will be somewhere inside of the 'attribute' field of the dataframe. The infomation varies from gene to gene.

4. Answer the following questions about your random feature.
   - What type of feature is your feature? (protein? rRNA? tRNA?)
   - Try googling the name/description of your feature. You're looking for any info you can find about what the gene does in the cell.
      - Try searching in standard google
      - Try searching in google scholar
      - Try searching including the term 'bacteria' or 'prokaryote'
   - Write a one paragraph description of what your gene is and what it does in prokaryotic cells. The information will vary depending on which random gene you pull, but try to dig up as many interesting/informative tidbits as you can! Include sources (wikipedia is ok, but try to find one academic publication)

5. Submit your responses for the weekly write-up assignment.


[Back to Top](#top)
