---
layout: default
---

<a name="top"></a>


# Week 9B lecture and tutorial (seaweed antimicrobial peptides)
1. [Intro to algae by Dr. Fox](#algae)
2. [Intro to anti-microbial peptides by Dr. Cotten](#amp)
3. [Predicting AMPs in genome sequences](#pred)
4. [Tutorial assignment](#tut)
    

## <ins>**Intro to algae by Dr. Fox**</ins> <a name="algae"></a>

Devaleraea mollis (pacific dulse) is a red algae that grows along the Oregon coast.

Types of bioactives (secondary metabolites):
- Allelopathy effectors
- Pigments (for photosythesis)
- Anti-microbial peptides (AMPs)

Questions:
- How do growth conditions impact the 'bioactive profile'?
    - Which bioactives are most important to monitor?

![algae](/Images/Week09/algae.png)

## <ins>**Intro to anti-microbial peptides by Dr. Cotten**</ins> <a name="amp"></a>

![amps](/Images/Week09/amps.png)


## <ins>**Predicting AMPs in genome sequences**</ins> <a name="pred"></a>

### Sequence characteristics of AMPs
- Length of peptides?
- Chemical properties of peptides?
- Structural properties of peptides?
- **Located within larger proteins vs stand-alone ORFs?**
    - What is an 'encrypted' peptide (sensu Torres et al)?     

### Programs for predicting AMPs from sequence
#### Web-based
- [iAMPpred](http://cabgrid.res.in:8080/amppred/index.html)

- [AxPEP](https://app.cbbio.online/ampep/home)
    - Can input either proteome sequences or genome (DNA) sequences

#### Command line
- [amPEPpy](https://github.com/tlawrence3/amPEPpy)
    - Written in python
    - Has multi-threading capablities

#### Additional tools
- [ORFFinder](https://pypi.org/project/orffinder/)
    - A python tool for screening DNA sequences for ORFs

- [marcel](https://github.com/BigDataBiology/macrel)
    - pipeline for prediction from metagenomes

### References:
- `https://www.sciencedirect.com/science/article/abs/pii/S2211926423004150?via%3Dihub`
- `https://www.nature.com/articles/s41551-021-00801-1`
- `https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8181880/`

## <ins>**Tutorial assignment**</ins> <a name="tut"></a>

- This week's tutorial assignment is to choose your final individual project and begin outlining the specific aims of your project. The project should be related to either Squeakuences (discussed last lecture) or antimicrobial peptide prediction (discussed this lecture). 

- For this assignment, begin describing your vision for a project. We will work together next week to refine the scope of the project.


1. What is the overarching goal of your project?

2. Create a bulleted list of the steps you'll take to accomplish your goals.

3. What type(s) of input data will you need? In what file format(s) do you expect these data to be?

4. What types of output (including file format) do you expect to produce?

5. What tools/software do you expect to use for your project (e.g. python, HPC job submission, specific software packages?)

6. How will you summarize/visualize/interpret your results? 


[Back to Top](#top)
