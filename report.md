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

### Step 2: Create the dataframes Red and Green to store the red and green fluorescences respectively.
Let's create the two dataframes for storing red and green fluorescences using the functions getGreen and getRed, provided by Minfi. 
Creation of Red dataframe:

`Red <- data.frame(getRed(RGset))`

`dim(Red)`

`[1] 622399      8`

Creation of Green dataframe:

`Green <- data.frame(getGreen(RGset))`

`dim(Green)`

`[1] 622399      8`

### Step 3: Find Red and Green fluorescence associated to the address 39802405. 
In order to find red and green fluorescences for the specfic address the following command is performed in Red and Green dataframes respectively:

`Red[rownames(Red)=="39802405",]`

`Green[rownames(Green)=="39802405",]`

Sample | Red flourescence | Green flourescence | Type
------------ | ------------- | ------------ | -------------
X5775278051_R01C01 | 6382 | 7885 | II
X5775278051_R04C02 | 7313 | 9600 | II
X5775278078_R02C01 | 4963 | 6844 | II
X5775278078_R05C01 | 5786 | 7405 | II
X5775278078_R05C02 | 5747 | 7203 | II
X5930514034_R01C02 | 6152 | 7582 | II
X5930514035_R04C02 | 5715 | 10793 | II
X5930514035_R06C02 | 6268 | 10337 | II

Since is a the address is associated with a type II probe, no color is specified in the table above. It is possible to check if the address is associated with a type I or type II probes in the following way: Let's first load the cleaned Illumina450Manifest in R workspace:

`load("/Users/giorg/Desktop/DNA_RNA DYNAMICS/MODULE_2/Lesson 2-20200514/Illumina450Manifest_clean.RData")`

Then the assigned address is searched in the Illumina450Manifest both for AddressA_ID and AddressB_ID:

`Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="39802405",]`

`Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID=="39802405",]`

In paritcular, no AddressB_ID is found with the selected code; otherwise, AddressA_ID with code 39802405 is associted to the probe  cg01086462, a type II probe localized on the forward strand on chromosome Y. 

### Step 4: Create the object MSet.raw 
The function MSet.raw allows to extract methylated and unmethylated signals from the raw data in RGset. 

`MSet.raw <- preprocessRaw(RGset)`

`MSet.raw`

It is possible to sae MSet.raw object as RData fiel:

`save(MSet.raw,file="MSet_raw.RData")`

### Step 5.1: Perform the quality check with QCplot.
The quality control is first perfomed with the function getQC, that estimates sample-specific quality checks for methylation data 
The result of getQC function is a dataframe with two columns referring to the chipwide medians of the mthylated and unmethilated channels.

`qc <- getQC(MSet.raw)`

Then it is possible to plot the results.

`plotQC(qc)`

![GitHub Logo](/images/logo.png)



