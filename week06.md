---
layout: default
---

<a name="top"></a>

# Week 6 lecture and tutorial
1. [Comparative genomics and genome evolution](#comparative)
	- **A.** [Large chromosomal mutations](#inversions)
	- **B.** [Whole genome duplication](#duplication)
2. [Emergent properties revealed by whole-genome analyses](#emergent)
	- **A.** [Phylogenetic incongruence](#incongruence)
	- **B.** [Hybridization and introgression](#int)
3. [Phylogenomic analyses](#phylogenomic)
4. [Bioinformatic pipelines](#pipelines)
5. [Collaborative code development with git and github](#git)
	- **A.** [Setting up your first github repo](#repo)
	- **B.** [Cloning your remote repo to a local machine](#clone)
	- **C.** [Configuration of your github credentials](#config)
6. [Tutorial assignment](#tut)


## <ins>**Comparative genomics and genome evolution**</ins> <a name="comparative"></a>

![comp00](/Images/Week06/comp00.png)

### <ins>**Large chromosomal mutations**</ins> <a name="inversions"></a>

![comp04](/Images/Week06/comp04.png)

![comp05](/Images/Week06/comp05.png)

### <ins>**Whole genome duplication**</ins> <a name="duplication"></a>

![comp01](/Images/Week06/comp01.png)

![comp02](/Images/Week06/comp02.png)

![comp03](/Images/Week06/comp03.png)

## <ins>**Emergent properties revealed by whole-genome analyses**</ins> <a name="emergent"></a>

![comp07](/Images/Week06/comp07.png)

![covid](/Images/Week06/covid.png)

### <ins>**Phylogenetic incongruence**</ins> <a name="incongruence"></a>

![comp06](/Images/Week06/comp06.png)

### <ins>**Hybridization and introgression**</ins> <a name="int"></a>

![comp08](/Images/Week06/comp08.png)

![comp09](/Images/Week06/comp09.png)

## <ins>**Phylogenomic analyses**</ins> <a name="phylogenomic"></a>

![comp10](/Images/Week06/comp10.png)

![comp11](/Images/Week06/comp11.png)

![comp12](/Images/Week06/comp12.png)

![comp13](/Images/Week06/comp13.png)

## <ins>**Bioinformatic pipelines**</ins> <a name="pipelines"></a>

![comp14](/Images/Week06/comp14.png)

![comp15](/Images/Week06/comp15.png)

## <ins>**Collaborative code development with git and github**</ins> <a name="git"></a>

![comp16](/Images/Week06/comp16.png)

![comp17](/Images/Week06/comp17.png)

![comp18](/Images/Week06/comp18.png)

### <ins>**Setting up your first github repo**</ins> <a name="repo"></a>

![comp19](/Images/Week06/comp19.png)

![comp20](/Images/Week06/comp20.png)

![comp21](/Images/Week06/comp21.png)

![comp22](/Images/Week06/comp22.png)

### <ins>**Cloning your remote repo to a local machine**</ins> <a name="clone"></a>

![comp23](/Images/Week06/comp23.png)

1. Copy the URL from github
2. In your command line on the HPC, navigate to the folder in which you want to create the repo
3. Clone with repo with: `git clone <paste-url-here>`
4. Use the following to check what branch you're on in the repo (this also confirms you're in a repo): `git branch -a`
5. You can now make edits to files or create files
6. To 'push' the edits to the remote repo do the following steps in order:
- Add the files with: `git add *`
- Make a commit with: `git commit -m "a quick message describing what you just did"`
- Push to the remote repo with: `git push origin main`

### <ins>**Configuration of your github credentials**</ins> <a name="config"></a>

1. **Open Terminal or Command Prompt:**
   Depending on your operating system, open Terminal (on macOS and Linux) or Command Prompt (on Windows).

2. **Enter Git configuration commands:**
   Use the following commands to configure Git with your access token and GitHub username:

   ```bash
   git config --global credential.helper store
   git config --global user.name "Your GitHub Username"
   git config --global user.email "your_email@example.com"
   ```

   Replace `"Your GitHub Username"` with your GitHub username and `"your_email@example.com"` with the email associated with your GitHub account.

3. **Store the Access Token:**
   Run the following command to store your access token:

   ```bash
   git credential approve
   ```

   This will prompt you to enter your GitHub username and password. Enter your username and paste the access token you copied earlier when prompted for the password.

4. **Verify Configuration:**
   You can verify that your configuration is set up correctly by running:

   ```bash
   git config --list
   ```

## <ins>**Tutorial assignment**</ins> <a name="tut"></a>
- Clone your newly made github repo to the HPC
- Create a new text file to develop into a python script.
	- Should be named something like: `phylogenomic_analyses.py`
- Outline the major steps of your phylogenomic pipeline (see the flow chart above)
	- Do this outlining inside of your new .py script by adding a series of comments inside of the script in order.
 		- e.g. `#Create a list of unaligned files that we'll need to align`
- `git add`, `git commit`, and `git push` your repo to the remote repository.
- Submit the url of your git account on canvas as the tutorial assignment.    


[Back to Top](#top)
