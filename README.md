<h1>slidinghet</h1>
The R script to perform the sliding window analysis performed in Velo-Antón et. al. 2017 in prep

In this project you can find an example file (acanthodactylusDataset.csv) and the RScript to run the analysis (slidingHet.R). A short description of the script can be found below as well as an example on how to use to reproduce the analysis. The script is provided as is.

<h3>Dataset</h3>
The dataset acanthodactylusDataset.csv contains the genotypes and 16S sequence for all 130 individuals used in Velo-Antón <em>et. al.</em> 2018 <em>in prep</em>

<h3>Functions</h3>

This script contains a series of functions meant to be used together.

He: Given a data.frame containing a column representing sample CODE and a series of columns representing alleles per locus, calculates expected heterozygosity of the given dataset.


nDiv: Given a data.frame containing a column representing sample CODE and a column with sequences, calculates the nucleotide diversity of the given dataset.


dataIntegrity: A preliminary function meant to test the integrity of the dataset, it checks if the necessary columns are present (CODE,SPECIES,LATITUDE,LONGITUDE). It also checks if the number of allele columns are pair (Diploid individuals)


linearCoordinates: This function takes as input the data data.frame, and transforms one of the two axis using a linear model to summarize the spatial distributions of samples into a single axis. Parameters: (data, transformAxis = "Longitude"(default) or "Latitude")


slidingWindow: This function performs the sliding window analysis, it takes the min and max transformed variable values to initiate the sliding window and per window estimates the the chosen statistics. Parameters: (lineardata,window_width= 0.5(default), slide_by = 1(default), seq_included = TRUE(default) or FALSE). window_width corresponds to the range of half the window (total range = 2*window_width); slide_by manages how much the window center moves per iteration.


dataVisualizer: A preliminary function for data visualization, it requires the package rworldmap. The idea of this function is to draw a map that discriminates how the slidingWindow split the data. Parameters: (lineardata, slidingWindowData)


<h3>Example:</h3>

```
source("a_script2.R")

input <- read.csv("acanthosDataset.csv", header = T,stringsAsFactors = F)



##Sample filtering to mimic manuscript dataset
##Remember that samples away from the coast were excluded from the sliding window analysis
input <- input[!(input$CODE %in% c("A765","6443","9999","13266","13267","13268","7245","13265")),]


##################################
####Check for data consistency####
##################################

dataIntegrity(input)                  
##Checks for header integrity. Its just a indication, even if 
##it passes, errors might still exist inside your data!!!

##################################
#######Linearize coordinates######
##################################

linearize <- linearCoordinates(input)
##linearize[[1]] contains the input dataframe with a new transformed column
##linearize[[2]] contains the coefficients of the linear model


#################################
#####Sliding Window analysis#####
#################################


sWindows <- slidingWindow(linearize)  
##sWindows[[1]] contains the sample window allocation
##sWindows[[2]] contains the results per window
##Two files are automatically generated at the working directory
###sampleGroups.csv which contains sWindows[[1]]
###perWindowResults.csv which contains sWindows[[2]]


################################
#####Results visualization######
################################

dataVisualizer (linearize,sWindows)  
##Again this graph should be used just as a preliminary representation
##The function requires ALOT of work before being usefull in any other way


```
