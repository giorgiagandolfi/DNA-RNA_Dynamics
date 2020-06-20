# Report DNA/RNA Dynamics
## Analysis of a DNA methylation dataset performed with R software and R Packages. 
### Step 1: Load raw data with minfi and create an object called RGset storing the RGChannelSet object.
First of all, it is important to remove all the varibale from R workspace and set the working directory containing the raw data of DNA methylation dataset. 

`rm(list=ls())`

`getwd()`

`setwd()`

Then let's load minfi package, previously installed from Bioconductor. Minfi is a flexible and comprehensive tool for the analysis of Infifium DNA methylation microarrays. Minfi loading is performed using the suppressMessages function to avoid displaying messages. 

`suppressMessages(library(minfi))`

Let's now import the raw data from the folder /Input_Data. In pariticulary, we are intereted to the SampleSheet csv file: we can load this file using the function `read.metharray.sheet` and then assign it to a object of class RGChannelSet using the function `read.metharray.exp`.

`SampleSheet <- read.table("Input_Data/Samplesheet_report_2020.csv",sep=",",header=T)`

`baseDir <- ("Input_data")`

`targets <- read.metharray.sheet(baseDir)`

`RGset <- read.metharray.exp(targets = targets)`

`save(RGset,file="RGset.RData")`



