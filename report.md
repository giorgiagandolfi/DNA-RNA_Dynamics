# Report DNA/RNA Dynamics
## Analysis of a DNA methylation dataset performed with R software and R Packages. 
### Step 1: Load raw data with minfi and create an object called RGset storing the RGChannelSet object.
First of all, it is important to remove all the varibale from R workspace and set the working directory containing the raw data of DNA methylation dataset. 

```
rm(list=ls())
getwd()
setwd()
```

Then let's load minfi package, previously installed from Bioconductor. Minfi is a flexible and comprehensive tool for the analysis of Infifium DNA methylation microarrays. Minfi loading is performed using the suppressMessages function to avoid displaying messages. 

`suppressMessages(library(minfi))`

Let's now import the raw data from the folder /Input_Data. In pariticulary, we are intereted to the SampleSheet csv file: we can load this file using the function `read.metharray.sheet` and then assign it to a object of class RGChannelSet using the function `read.metharray.exp`.

```
SampleSheet <- read.table("Input_Data/Samplesheet_report_2020.csv",sep=",",header=T)
baseDir <- ("Input_data")
targets <- read.metharray.sheet(baseDir)
RGset <- read.metharray.exp(targets = targets)
save(RGset,file="RGset.RData")
```

### Step 2: Create the dataframes Red and Green to store the red and green fluorescences respectively.
Let's create the two dataframes for storing red and green fluorescences using the functions getGreen and getRed, provided by Minfi. 
Creation of Red dataframe:

```
Red <- data.frame(getRed(RGset))
dim(Red)
[1] 622399      8
```
Creation of Green dataframe:
```
Green <- data.frame(getGreen(RGset))
dim(Green)
[1] 622399      8
```

### Step 3: Find Red and Green fluorescence associated to the address 39802405. 
In order to find red and green fluorescences for the specfic address the following command is performed in Red and Green dataframes respectively:
```
Red[rownames(Red)=="39802405",]
Green[rownames(Green)=="39802405",]
```

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
```
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="39802405",]
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID=="39802405",]
```

In paritcular, no AddressB_ID is found with the selected code; otherwise, AddressA_ID with code 39802405 is associted to the probe  cg01086462, a type II probe localized on the forward strand on chromosome Y. 

### Step 4: Create the object MSet.raw 
The function MSet.raw allows to extract methylated and unmethylated signals from the raw data in RGset. 
```
MSet.raw <- preprocessRaw(RGset)
MSet.raw
```

It is possible to sae MSet.raw object as RData fiel:

`save(MSet.raw,file="MSet_raw.RData")`

### Step 5.1: Perform the quality check with QCplot.
The quality control is first perfomed with the function getQC, that estimates sample-specific quality checks for methylation data 
The result of getQC function is a dataframe with two columns referring to the chipwide medians of the mthylated and unmethilated channels.

`qc <- getQC(MSet.raw)`

Then it is possible to plot the results.

`plotQC(qc)`

![QCplot](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/qcplot.png)

As it can be seen, all eight samples have a good quality, since they are plotted on the right side of the line.  
### Step 5.2: Check the intensity of negative controls using minfi.
Negative controls are sample-dependent controls that by hybridizing with the sample DNA allow to evalute the performance across samples. Generally, the signal intensities of negative controls vary between 100 and 1000 units and higher values lead to poor DNA template quality. Let's evaluate the singal intenisties of negative controls in the dataset taking advantage of the controlStripPlot function (must specify the type of controls probes).

`controlStripPlot(RGset, controls="NEGATIVE")`

mettere immagine
As it can be seen from the plot, all negative controls for both red and green flourescence are good, since their values are below l000 (log2(1000) = 10). 
### Step 5.3: Calculate detection pValues and count how many probes have a dection pValue higher than 0.05.
Let's calculate the detection pValue with the function detectionP that identifies failed positions defined as both the methylated and unmethylated channel reporting background singals. 

`detP <- detectionP(RGset)`

In order to evaluate how many probes have a pValue higher than the selected threshold (i.e. 0.05), a new variable is generated by filtering the detectionP result.

`failed <- detP>0.05`

In the summary, the boolean "TRUE" refers to the number of probes which have pValue higher than the threshold, while the boolean "FALSE" refers to those probes whose pValue is lower than or equal to 0.05. The table below summurizes the results.

Sample | Failed positions |
------------ | ------------- |
X5775278051_R01C01 | 247 |
X5775278051_R04C02 | 210 |
X5775278078_R02C01 | 264 |
X5775278078_R05C01 | 413 |
X5775278078_R05C02 | 385 |
X5930514034_R01C02 | 91 |
X5930514035_R04C02 | 46 |
X5930514035_R06C02 | 115 |

Only 1771 among 3884096 probes have no significant pValue (0.045% of the total number of probes).

### Step 6: Calculate raw beta and M values and plot the densities of mean methylation values, dividing the samples in DS and WT.
Let's divide the samples in Wild Type and with Down Sydrome, thus according to the field Group, first in the SampleSheet file then in the MSet.raw file.
```
wt <- SampleSheet[SampleSheet$Group=="WT", "Basename"]
ds <- SampleSheet[SampleSheet$Group=="DS", "Basename"]
wtset<-MSet.raw[,colnames(MSet.raw) %in% wt]
dsset<-MSet.raw[,colnames(MSet.raw) %in% ds]
```

Now raw M and beta values are extracted for both groups using getM and getBeta functions, respectively. 

Then, the mean of each M and beta values of probes are calculated and plotted according to their density distributions.

```
mean_of_beta_wtset <- apply(beta_wtset,1,mean,na.rm=T)
mean_of_beta_dsset <- apply(beta_dsset,1,mean,na.rm=T)

d_mean_of_beta_wtset <- density(mean_of_beta_wtset)
d_mean_of_beta_dsset <- density(mean_of_beta_dsset)

mean_of_M_wtset <- apply(M_wtset,1,mean,na.rm=T)
mean_of_M_dsset <- apply(M_dsset,1,mean,na.rm=T)

d_mean_of_M_wtset <- density(mean_of_M_wtset)
d_mean_of_M_dsset <- density(mean_of_M_dsset)
```

The denisty distributions of beta and M values are almost overlapping in the two groups, as we can see from the plots below. 
mettere plots

### Step 7: Normalize the data using the function preprocessSWAN to each student and compare raw data and normalized data.
Normalization procedure aims to resolving the systematic errors and bias introduced by the microarray experimental platform. In this step Subset-quantile Within Array Normalization (SWAN) method is performed using preprocessSWAN function. The effects of the normalization procedure are then evaluated comparing densities of mean and standard deviation beta values. Boxplot representation for beta values will be used. The comparison will consider the chemistry types probe as well. 
Raw data:

```
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
dfI <- droplevels(dfI)

dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
dfII <- droplevels(dfII)`

beta <- getBeta(MSet.raw)
beta_I <- beta[rownames(beta) %in% dfI$IlmnID,]
beta_II <- beta[rownames(beta) %in% dfII$IlmnID,]

mean_of_beta_MSet_I <- apply(beta_I,1,mean)
mean_of_beta_MSet_II <- apply(beta_II,1,mean)
d_mean_of_beta_MSet_I <- density(mean_of_beta_MSet_I, na.rm=T)
d_mean_of_beta_MSet_II <- density(mean_of_beta_MSet_II, na.rm=T)

sd_of_beta_I <- apply(beta_I,1,sd,na.rm=T)
sd_of_beta_II <- apply(beta_II,1,sd,na.rm=T)
d_sd_of_beta_I <- density(sd_of_beta_I)
d_sd_of_beta_II <- density(sd_of_beta_II)
```
Normalized data:

```
preprocessSWAN_results <- preprocessSWAN(RGset)
beta_preprocessSWAN <- getBeta(preprocessSWAN_results)
beta_preprocessSWAN_I <-beta_preprocessSWAN[rownames(beta_preprocessSWAN) %in% dfI$IlmnID,]
beta_preprocessSWAN_II <- beta_preprocessSWAN[rownames(beta_preprocessSWAN) %in% dfII$IlmnID,]

mean_of_beta_preprocessSWAN_I <- apply(beta_preprocessSWAN_I,1,mean)
mean_of_beta_preprocessSWAN_II <- apply(beta_preprocessSWAN_II,1,mean)
d_mean_of_beta_preprocessSWAN_I <- density(mean_of_beta_preprocessSWAN_I,na.rm=T)
d_mean_of_beta_preprocessSWAN_II <- density(mean_of_beta_preprocessSWAN_II,na.rm=T)

sd_of_beta_preprocessSWAN_I <- apply(beta_preprocessSWAN_I,1,sd)
sd_of_beta_preprocessSWAN_II <- apply(beta_preprocessSWAN_II,1,sd)
d_sd_of_beta_preprocessSWAN_I <- density(sd_of_beta_preprocessSWAN_I,na.rm=T)
d_sd_of_beta_preprocessSWAN_II <- density(sd_of_beta_preprocessSWAN_II,na.rm=T)
```
mettere plot+legenda
mettere commento
### Step 8: Perform a PCA on the beta matrix generated in step 7.






