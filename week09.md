---
layout: default
---

<a name="top"></a>

# Week 9 lecture and tutorial
1. [Guest lecture introduction](#intro)
2. [ERCgo background](#ERCgo)
3. [Tutorial assignment](#tut)
4. [Weekly project write-up assignment](#write)

## <ins>**Guest lecture introduction**</ins> <a name="intro"></a>
This week's lecture will be led by Linnea Lane, a Research Technician in the Forsythe Lab. Linnea graduated from OSU-Cascades in 2022 with the Biology degree and Computer Science minor. Linnea created a new python-based software package, which we will learn about this week!

## <ins>**ERCgo background**</ins> <a name="ERCgo"></a>

![ercgo1](/Images/Week09/Slide1.png)
![ercgo2](/Images/Week09/Slide2.png)
![ercgo3](/Images/Week09/Slide3.png)
![ercgo4](/Images/Week09/Slide4.png)
![ercgo5](/Images/Week09/Slide5.png)
![ercgo6](/Images/Week09/Slide6.png)
![ercgo7](/Images/Week09/Slide7.png)
![ercgo8](/Images/Week09/Slide8.png)
![ercgo9](/Images/Week09/Slide9.png)
![ercgo10](/Images/Week09/Slide10.png)
![ercgo11](/Images/Week09/Slide11.png)
![ercgo12](/Images/Week09/Slide12.png)
![ercgo13](/Images/Week09/Slide13.png)
![ercgo14](/Images/Week09/Slide14.png)
![ercgo15](/Images/Week09/Slide15.png)
![ercgo16](/Images/Week09/Slide16.png)
![ercgo17](/Images/Week09/Slide17.png)
![ercgo18](/Images/Week09/Slide18.png)
![ercgo19](/Images/Week09/Slide19.png)
![ercgo20](/Images/Week09/Slide20.png)
![ercgo21](/Images/Week09/Slide21.png)
![ercgo22](/Images/Week09/Slide22.png)
![ercgo23](/Images/Week09/Slide23.png)



## <ins>**Tutorial assignment**</ins> <a name="tut"></a>

Do the following before class on Thursday:
- Clone the ERCgo git repository (to your scratch dir on the HPC). Use this `git clone https://github.com/linxlane/ERCgo.git`
- Make sure that you're in a specific branch of the ERCgo repo (the BB485 branch). To see what branch you're in: `git branch -a`. To change to the needed branch: `git checkout BB485`
- Use the included .yml file to create a new conda environment for running ERCgo. Use this command: `conda env create -f environment.yml -y`
- Respond "Done!" on canvas once you've finished the above steps.

## <ins>**Weekly project write-up assignment**</ins> <a name="write"></a>

For your project write-up this week, please do the tasks below and answer the associated questions. Submit your writeup on Canvas.

1. Run ERCgo on the provided interaction dataset and answer the following questions:
	- 1A: Do we see a general correlation between the confidence value for interactions and the GO overlap score that ERCgo calculates? What pieces of qualitative/quantitative evidence help you assess this?
	- 1B: Where do the Clp interactions fall in your plot? How does this compare to the ERC-based network that Linnea showed?
	- 1C: Do the 'observed' interactions that exist in our network tend to have a higher degree of GO-overlap than 'non-hits'? What evidence are you using to assess this? 
2. ERCgo should provide a tsv file that includes all of the GO-overlap scores for the network. Open this file in excel (or google sheets) and sort it so the highest GO overlap scores are at the top.
	- 2A: For the top five interactions, look up information on gene A and gene B involved in each interaction. Look up functional information by searching the ATXGXXXXXX ID at [The Arabidopsis Information Resource website](https://www.arabidopsis.org/)

| Gene A ID | Gene A TAIR Description | Gene B ID | Gene B TAIR Description | GO overlap score | Does the overlap score make sense? |
|-----------|--------------------------|-----------|--------------------------|-----------------|-----------------------------------|
|           |                          |           |                          |                 |                                   |
|           |                          |           |                          |                 |                                   |
|           |                          |           |                          |                 |                                   |
|           |                          |           |                          |                 |                                   |
|           |                          |           |                          |                 |                                   |

3. You are the first 'beta testers' of ERCgo. To help us improve ERCgo and get practice using github features, I would like each student to create an "Issue" on the ERCgo github page. As you run ERCgo, keep an eye out for bugs or even just wonky implimentation/features/output. Create a github Issue that describes it and suggests an improvement. Alternatively, if you do not see any bugs, then create an issue that described a cool new feature that you would like to see in ERCgo.
	- 3A: Respond "Done!" to let me know you created an issue. I will look on the ERCgo git repo to check out your issue. 
	
	
[Back to Top](#top)
