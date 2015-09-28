# ROOTFIT


ROOT-FIT is descriptive model of fitted growth functions to root system architecture dynamics. ROOT-FIT was developed by Magdalena Julkowska and Christa Testerink. The R version was made by Guillaume Lobet and Magdalena Julkowska.


## Running ROOT-FIT

To launch ROOT-FIT, you need to have R installed on your computer. This is the only prerequisite.

To run ROOT-FIT: use the following command:


	library(shiny)
	shiny::runGitHub("guillaumelobet/rootfit", "guillaumelobet") 
	
	
PS: an internet connection is needed, as the most recent version of ROOT-FIT will be download from the server.


## Using ROOT-FIT

![](https://raw.githubusercontent.com/guillaumelobet/rootfit/master/www/rootfit_ui.png)

1. Load your .csv datafile
2. Choose the fitting type
3. Set the analysis options (these will be updated when the cvs file is loaded)
	- Time of treatment: day 0 of your experiment (by default the smallest day value in the datafile)
	- Reference genotype: used for the computation of the relative growth factors
	- Reference treatment: used for the computation of the relative growth factors
	- Find best fit on # samples: number of samples to used to find the best fitting function
4. Choose the genotypes to plot (these will be updated when the cvs file is loaded)
5. Unleash ROOT-FIT	
6. Inspect your data and download the plots
7. Download the results

## Data structure 

The input datafile for ROOT-FIT  needs to be in a cvs form, with the following columns:

	Plate.name | Media | Genotype | Plant.id | Age | MR.path.length | Number.LR.MR | Total.root.size
	
- **Plate.name** : the name of the image
- **Media** : the media type, or treatment
- **Genotype** : the genotype
- **Plant.id** : the identifier of the plant (needs to be unique by image / genotype / treatment)
- **Age** : age of the plant, in days
- **MR.path.length** : the length of the primary root
- **Number.LR.MR** : the number of lateral roots
- **Total.root.size** : the total root system size (primary + laterals)