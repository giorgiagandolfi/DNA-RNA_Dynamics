# Report DNA/RNA Dynamics
## Analysis of a DNA methylation dataset performed with R software and R Packages. 
### Step 1: Load raw data with minfi and create an object called RGset storing the RGChannelSet object.
First of all, it is important to remove all the variables from R workspace and set the working directory containing the raw data of DNA methylation dataset. 

```
rm(list=ls())
getwd()
setwd()
```

Then let's load **minfi** package, previously installed from Bioconductor. Minfi is a flexible and comprehensive tool for the analysis of Infifium DNA methylation microarrays. Minfi loading is performed using the suppressMessages function to avoid displaying messages. 

`suppressMessages(library(minfi))`

Let's now import the raw data from the folder /Input_Data. In pariticulary, we are intereted to the SampleSheet csv file: we can load this file using **read.metharray.sheet** function and then assign it to a object of class RGChannelSet using **read.metharray.exp** function.

```
SampleSheet <- read.table("Input_Data/Samplesheet_report_2020.csv",sep=",",header=T)
baseDir <- ("Input_data")
targets <- read.metharray.sheet(baseDir)
RGset <- read.metharray.exp(targets = targets)
save(RGset,file="RGset.RData")
```

### Step 2: Create the dataframes Red and Green to store the red and green fluorescences respectively.
Let's create the two dataframes for storing red and green fluorescences using **getGreen** and **getRed** functions, provided by Minfi. 

Generation of Red dataframe:

```
Red <- data.frame(getRed(RGset))
dim(Red)
[1] 622399      8
```
Generation of Green dataframe:
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

Since the address is associated with a type II probe, no color is specified in the table above. It is possible to check if the address is associated with a type I or type II probes in the following way: Let's first load the cleaned Illumina450Manifest in R workspace:

`load("/Users/giorg/Desktop/DNA_RNA DYNAMICS/MODULE_2/Lesson 2-20200514/Illumina450Manifest_clean.RData")`

Then the assigned address is searched in the Illumina450Manifest both for AddressA_ID and AddressB_ID:
```
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="39802405",]
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID=="39802405",]
```

In paritcular, no AddressB_ID is found with the selected address; otherwise, AddressA_ID with code 39802405 is associted to the probe  cg01086462, a type II probe localized on the forward strand on chromosome Y. 

### Step 4: Create the object MSet.raw 
**MSet.raw** function allows to extract methylated and unmethylated signals from the raw data in RGChannelSet.
```
MSet.raw <- preprocessRaw(RGset)
MSet.raw
```

It is possible to sae MSet.raw object as RData file:

`save(MSet.raw,file="MSet_raw.RData")`

### Step 5.1: Perform the quality check with QCplot.
The quality control is first perfomed with **getQC** function, that estimates sample-specific quality checks for methylation data 
The result of getQC function is a dataframe with two columns referring to the chipwide medians of the mthylated and unmethilated channels.

```
qc <- getQC(MSet.raw)
qc
DataFrame with 8 rows and 2 columns
                       mMed      uMed
                  <numeric> <numeric>
5775278051_R01C01   11.7616   11.8222
5775278051_R04C02   12.0427   12.0668
5775278078_R02C01   11.5774   11.6170
5775278078_R05C01   11.7645   11.7444
5775278078_R05C02   11.7288   11.7241
5930514034_R01C02   11.5038   11.6416
5930514035_R04C02   11.7211   11.7661
5930514035_R06C02   11.9436   11.9035
```

Then it is possible to plot the results.

`plotQC(qc)`

![QCplot](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/qcplot.png)

As it can be seen, all eight samples have a good quality, since they are plotted on the right side of the line.  
### Step 5.2: Check the intensity of negative controls using minfi.
**Negative controls** are sample-dependent controls that, by hybridizing with the sample DNA, allow to evalute the performance across samples. Generally, the signal intensities of negative controls vary between 100 and 1000 units and higher values lead to poor DNA template quality. Let's evaluate the singal intenisties of negative controls in the dataset taking advantage of **controlStripPlot** function (must specify the type of controls probes).

`controlStripPlot(RGset, controls="NEGATIVE")`

![control negative](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/controls_neg.png)

As it can be seen from the plot, all negative controls for both red and green flourescence are good, since their values are below l000 (log2(1000) = 10). 
### Step 5.3: Calculate detection pValues and count how many probes have a dection pValue higher than 0.05.
Let's calculate the detection pValue with **detectionP** function that identifies failed positions, such as those methylated and unmethylated channels reporting background singals. 

`detP <- detectionP(RGset)`

In order to evaluate how many probes have a pValue higher than the selected threshold (i.e. 0.05), detectionP results are filtered up.

`failed <- detP>0.05`

In the summary, the boolean "TRUE" refers to the number of probes which have pValue higher than the threshold, while the boolean "FALSE" refers to those probes whose pValue is lower than or equal to 0.05. The table below summurizes the results according each samples.

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
Let's divide the samples in Wild Type (WT) and with Down Sydrome (DS), thus according to the field Group, first in the SampleSheet file then in the MSet.raw file.
```
wt <- SampleSheet[SampleSheet$Group=="WT", "Basename"]
ds <- SampleSheet[SampleSheet$Group=="DS", "Basename"]
wtset<-MSet.raw[,colnames(MSet.raw) %in% wt]
dsset<-MSet.raw[,colnames(MSet.raw) %in% ds]
```

Now raw M and beta values are extracted for both groups using getM and getBeta functions, respectively. 

Then, the means of each M and beta values of probes are calculated and plotted according to their density distributions.

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
Beta values and M values can be plotted both separately and in overlapping graphs.

```
#single plots
par(mfrow=c(2,2))
plot(d_mean_of_beta_wtset,main="Density of Beta Values wt",col="orange")
plot(d_mean_of_beta_dsset,main="Density of Beta Values ds",col="orange")
plot(d_mean_of_M_wtset,main="Density of M Values wt",col="purple")
plot(d_mean_of_M_dsset,main="Density of M Values ds",col="purple")
#overlapping plots
par(mfrow=c(1,2))
plot(d_mean_of_beta_wtset,main="Density of Beta Values",col="orange")
lines(d_mean_of_beta_dsset,col="green")
plot(d_mean_of_M_wtset,main="Density of M Values",col="purple")
lines(d_mean_of_M_dsset, col="green")
```

The denisty distributions of beta and M values are almost overlapping in the two groups, as we can see from the plots below. 

![plots](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/beta_m_ds_vs_wt.png)
![overlapping plots](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/beta_and_M_overlapping.png)

### Step 7: Normalize the data using the function preprocessSWAN to each student and compare raw data and normalized data.
Normalization procedure aims to resolving the systematic errors and bias introduced by the microarray experimental platform. In this step Subset-quantile Within Array Normalization (SWAN) method is performed using **preprocessSWAN** function. The effects of the normalization procedure are then evaluated comparing densities of mean and standard deviation beta values. Boxplot representation for beta values will be used. The comparison will consider the chemistry types probe as well. 
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
![plots](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/preprocessSWAN_plots.png)
It is possible to observe that normalized beta and M values are not so different from those obtained from the raw data, despite an higher peak in 
### Step 8: Perform a PCA on the beta matrix generated in step 7.
Principal Component Analysis (**PCA**) is an exploratory tool to find predominant gene expression pattern along the identified componetents as well as to detect possible outliers and batch effects. PCA reduces the high-dimensionality space by finding the gratest variances in the data. In R enviroment, PCA is performed using **prcomp()** function, that takes as argument the matrix of normalized beta values. Then a **scree plot** showing the cumulative variance explained by each principal component is generated.

```
pca_results <- prcomp(t(beta_preprocessSWAN), scale=T)
plot(pca_results, col="pink",main="Scree plot")
```

![scree plot](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/Scree%20plot.png)

It is possible to visualize the results in the plot by colouring the dots according to a particular phenotype, i.e. the group (either WT or DS). 

```

![pca plot](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/pca_group.png)

### Step 9: Using the matrix of normalized beta values generated in step 7, identify differentially methylated probes between group DS and group WT using the Mann-Whitney test. 
Let's identify deifferentialy methylated probes between group DS and group WT using a non-parametric test called Mann-Witheny or Wilcoxon rank-sum test. It is the equivalent of the most-known Student's t-test for unpaired data when no assumptions regarding normal distribution of data is made. We will use wilcox.test function for each row of the dataframe of preprocessSWAN beta values; thus an *ad hoc* function is created and applied to preprocessSWAN beta values.

```
My_mannwhitney_function <- function(x) {
  wilcox <- wilcox.test(x~ pheno$Group)
  return(wilcox$p.value)
} 
pValues_wilcox <- apply(beta_preprocessSWAN,1, My_mannwhitney_function)
```

We can calculate how many probes are differntialy methylated filtering the resulting pValues according to a threshold of 0.05. 

```
final_wilcox<-data.frame(beta_preprocessSWAN,pValues_wilcox)
final_wilcox <- final_wilcox[order(final_wilcox$pValues_wilcox),]
finaL_wilcox_0.05<-final_wilcox[final_wilcox$pValues_wilcox<=0.05,]
dim(finaL_wilcox_0.05)
[1] 22351     9
```
22351 CpG probes are differentially methylated according to the Mann-Whitney test.

### Step 10: Apply multiple test correction and set a significant threshold of 0.05. How many probes do you identify as differentially methylated considering nominal pValues? How many after Bonferroni correction? How many after BH correction?
When the number of probes is in the order of thousands, a correction of the signficance threshold is required: it is called multiple test correction. As we can see, by applying a threshold of 0.05 to the raw pValues of Mann-Whtiney test, a high number of false positives have been detected. In this step we are going to apply multiple test correction, such as Bonferroni correction and Benjamini-Hochberg correction, and to compare the results with the raw test.

```
#create a vector storing the raw pValues
raw_pValues<-final_wilcox[,9]

#apply BH and Bonferroni corrections
corrected_pValues_BH <- p.adjust(raw_pValues,"BH")
corrected_pValues_Bonf <- p.adjust(raw_pValues,"bonferroni")

#create a dataframe storing raw and corrected pValues for each probe
final_wilcox_corrected<-data.frame(final_wilcox, corrected_pValues_BH, corrected_pValues_Bonf)
```
Let's now see how many probes are differentially methylated after Bonferroni end BH correction.

```
dim(final_wilcox[final_wilcox$corrected_pValues_Bonf<=0.05,])
[1] 0 9
dim(final_wilcox[final_wilcox$corrected_pValues_BH<=0.05,])
[1] 0 9
```

None false postives are detected when multiple test correction is applied. A boxplot of the raw and corrected pValues shows their distributions.

![boxplot](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/boxplot%20pink.png)

### Step 11: Produce an heatmap of the top 100 differentially mehtylated probes.
Heatpmap is a data visualization technique showing the level of differentially expressed genes in a matrix with different colours. Heatmaps are always coupled with hierarchical clustering, whose aim is to link genes or samples with similar profiles to form
a dendrogram. Different linkage methods can be used to calculate the distance among clusters: (1) single linkage: (2) complete linkage and (3) average linkage. They are all implemented in heatmap.2 function. Let's first generate the input for the heatmaps and the color bar according to the phenotype group.

```
install.packages("gplots")
library(gplots)
input_heatmap=as.matrix(final_wilcox_corrected[1:100,1:8])
pheno$Group
colorbar<-c("red", "red", "blue","blue","blue","red", "red", "blue")

#complete linkage (default)
heatmap.2(input_heatmap,col=cm.colors(100),Rowv=T,Colv=T,dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F)
#single linkage
heatmap.2(input_heatmap,col=cm.colors(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'single'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F)
#average linkage
heatmap.2(input_heatmap,col=cm.colors(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'average'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F)
```
#### Complete linkage
![complete linkage](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/complete_link_ord.png)

#### Single linkage
![single linkage](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/single_link.png)

#### Average linkage
![average linkage](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/average_link.png)

A dendrogram generated using complete linkage method is charatcterized by compact and well-defined clusters both for the probes and for the samples,on the right side and on the upper side of the heatmap, respectively. While the dendrograms generated with single linkage method show the chaining effect, as we can see expecially for the one of the probes. However, all heatmaps show a clear division between WT samples (blue) and DS samples (red), as well as a distiction between hypermethylated (highly expressed genes) and ipomethylated (poorly expressed genes) CpG islands. 

### Step 12: Produce a volcano plot and a Manhattan plot of the results of differential methylation analysis.
A volcano plot is a type of scatterplot that displays the statistical significance against the fold change. Since in this case, the pValues have been generated using a non-parametric test, the typical volcano shape of the plot is missed. The fold change is calculated as the difference between the average beta values of WT samples and the average beta values of DS samples, in the following way.

```
beta_group <- final_wilcox_corrected[,1:8]
beta_groupDS <- beta_group[,SampleSheet$Group=="DS"]
mean_beta_groupDS <- apply(beta_groupDS,1,mean)
beta_groupWT <- beta_group[,SampleSheet$Group=="WT"]
mean_beta_groupWT <- apply(beta_groupWT,1,mean)
delta <- mean_beta_groupWT - mean_beta_groupDS
```

A data frame is created storing in one column, the delta values, and on the other, the -log10 of pValues; these values are then plotted using the plot function with some options.

```
toVolcPlot <- data.frame(delta, -log10(final_wilcox_corrected$pValues_wilcox))
plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5)
```

![volcano plot](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/volcano_plot.png)

Together with the volcano plot, the so-called Manhattan plot can be generated using the **gap** package. A Manhattan plot is a type of scatter plot useful for the visualization of large number of data, that are plotted according to the genomic coordinates (X-axis) and the negative log10 of the pValues. For the same reason of volcano plot, the typical profile of skyscrapers charterizing the Manhattan plots is lost. 
After downloading the required package, the dataframe of CpG probes with associted pValues has to be annotated by extracting genomic information from the Illumina450Manifest (cleaned version). The merge function is used to perform this step. However, since the CpG probes are stored in the rownames and not in the columns (as required by merge() funtion), a new dataframe storing the CpG probe IDs on the column is created.

```
final_wilcox_corrected_col <- data.frame(rownames(final_wilcox_corrected), final_wilcox_corrected)
head(final_wilcox_corrected_col)
colnames(final_wilcox_corrected_col) [1] <- "IlmnID"
final_wilcox_corrected_col_annotated <- merge(final_wilcox_corrected_col, Illumina450Manifest_clean, by="IlmnID")
dim(final_wilcox_corrected_col_annotated)
[1] 485512     44
```

Then the Manhattan plot is generated using **mhtplot** function.

```
input_Manhattan <- data.frame(final_wilcox_corrected_col_annotated$CHR, final_wilcox_corrected_col_annotated$MAPINFO, final_wilcox_corrected_col_annotated$pValues_wilcox)
dim(input_Manhattan)

levels(input_Manhattan$final_wilcox_corrected_col_annotated.CHR)

input_Manhattan$final_wilcox_corrected_col_annotated.CHR <- factor(input_Manhattan$final_wilcox_corrected_col_annotated.CHR,levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))
levels(input_Manhattan$final_wilcox_corrected_col_annotated.CHR)

palette <-c("red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","brown","red","blue","green")
palette

mhtplot(input_Manhattan, control=mht.control(colors=palette))
axis(2,cex=0.5)
```





