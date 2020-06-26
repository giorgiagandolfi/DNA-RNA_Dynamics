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
targets

Sample_Name Group Age      Slide  Array
1        1020    DS  29 5775278051 R01C01
2        1036    DS  34 5775278051 R04C02
3        3038    WT  46 5775278078 R02C01
4        3042    WT  32 5775278078 R05C01
5        3052    WT  31 5775278078 R05C02
6        1016    DS  43 5930514034 R01C02
7        1029    DS  32 5930514035 R04C02
8        3029    WT  35 5930514035 R06C02
                      Basename
1 Input_data/5775278051_R01C01
2 Input_data/5775278051_R04C02
3 Input_data/5775278078_R02C01
4 Input_data/5775278078_R05C01
5 Input_data/5775278078_R05C02
6 Input_data/5930514034_R01C02
7 Input_data/5930514035_R04C02
8 Input_data/5930514035_R06C02

RGset <- read.metharray.exp(targets = targets)
save(RGset,file="RGset.RData")
RGset
str(RGset)

class: RGChannelSet 
dim: 622399 8 
metadata(0):
assays(2): Green Red
rownames(622399): 10600313 10600322 ... 74810490 74810492
rowData names(0):
colnames(8): 5775278051_R01C01 5775278051_R04C02 ... 5930514035_R04C02 5930514035_R06C02
colData names(7): Sample_Name Group ... Basename filenames
Annotation
  array: IlluminaHumanMethylation450k
  annotation: ilmn12.hg19
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
       IlmnID       Name AddressA_ID
32 cg01086462 cg01086462    39802405
                                     AlleleA_ProbeSeq AddressB_ID
32 AATAAACCATCCCTAATTAACCCTACCTACTTAATATACACCATACCTAC            
   AlleleB_ProbeSeq Infinium_Design_Type Next_Base Color_Channel
32                                    II                        
                                                                                                               Forward_Sequence
32 CGCCTGCTTTAAATAGACCATCCCTGATTGACCCTGCCTGCTTAATGTACACCATACCTG[CG]TTATCGCACCTCCTTTAAATACACAACACTTGATTGACCCTGCCTGCTTTATGTAGACAC
   Genome_Build CHR MAPINFO
32           37   Y 7429349
                                            SourceSeq Chromosome_36
32 CGCAGGTATGGTGTACATTAAGCAGGCAGGGTCAATCAGGGATGGTCTAT             Y
   Coordinate_36 Strand Probe_SNPs Probe_SNPs_10 Random_Loci
32       7489349      R                                   NA
   Methyl27_Loci UCSC_RefGene_Name UCSC_RefGene_Accession
32            NA                                         
   UCSC_RefGene_Group UCSC_CpG_Islands_Name
32                     chrY:7428179-7428424
   Relation_to_UCSC_CpG_Island Phantom DMR Enhancer HMM_Island
32                     S_Shore                   NA           
   Regulatory_Feature_Name Regulatory_Feature_Group DHS
32                                                   NA

Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID=="39802405",]
 [1] IlmnID                      Name                       
 [3] AddressA_ID                 AlleleA_ProbeSeq           
 [5] AddressB_ID                 AlleleB_ProbeSeq           
 [7] Infinium_Design_Type        Next_Base                  
 [9] Color_Channel               Forward_Sequence           
[11] Genome_Build                CHR                        
[13] MAPINFO                     SourceSeq                  
[15] Chromosome_36               Coordinate_36              
[17] Strand                      Probe_SNPs                 
[19] Probe_SNPs_10               Random_Loci                
[21] Methyl27_Loci               UCSC_RefGene_Name          
[23] UCSC_RefGene_Accession      UCSC_RefGene_Group         
[25] UCSC_CpG_Islands_Name       Relation_to_UCSC_CpG_Island
[27] Phantom                     DMR                        
[29] Enhancer                    HMM_Island                 
[31] Regulatory_Feature_Name     Regulatory_Feature_Group   
[33] DHS                        
<0 rows> (or 0-length row.names)

```

In paritcular, no AddressB_ID is found with the selected address; otherwise, AddressA_ID with code 39802405 is associted to the probe  cg01086462, a type II probe localized on the forward strand on chromosome Y. Since it is a Type II probe, no color channel is associated to it (the emission flourescence depends on the methylation status of the target CpG).

### Step 4: Create the object MSet.raw 
**MSet.raw** function allows to extract methylated and unmethylated signals from the raw data in RGChannelSet.
```
MSet.raw <- preprocessRaw(RGset)
MSet.raw

class: MethylSet 
dim: 485512 8 
metadata(0):
assays(2): Meth Unmeth
rownames(485512): cg00050873 cg00212031 ... ch.22.47579720R
  ch.22.48274842R
rowData names(0):
colnames(8): 5775278051_R01C01 5775278051_R04C02 ...
  5930514035_R04C02 5930514035_R06C02
colData names(7): Sample_Name Group ... Basename filenames
Annotation
  array: IlluminaHumanMethylation450k
  annotation: ilmn12.hg19
Preprocessing
  Method: Raw (no normalization or bg correction)
  minfi version: 1.34.0
  Manifest version: 0.4.0
```

It is possible to sae MSet.raw object as RData file:

`save(MSet.raw,file="MSet_raw.RData")`

### Step 5.1: Perform the quality check with QCplot.
The quality control is first perfomed with **getQC** function, that estimates sample-specific quality checks for methylation data. 
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

```
table(failed)
failed
  FALSE    TRUE 
3882325    1771 
```

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
It is possible to observe that normalized beta and M values are not so different from those of raw data, despite an higher peak in the distribution of type II standrd deviation and a up-shifted median of beta values in the boxplots.
### Step 8: Perform a PCA on the beta matrix generated in step 7.
Principal Component Analysis (**PCA**) is an exploratory tool to find predominant gene expression pattern along the identified componetents as well as to detect possible outliers and batch effects. PCA reduces the high-dimensionality space by finding the gratest variances in the data. In R enviroment, PCA is performed using **prcomp()** function, that takes as argument the matrix of normalized beta values. Then a **scree plot** showing the cumulative variance explained by each principal component is generated.

```
pca_results <- prcomp(t(beta_preprocessSWAN), scale=T)
print(summary(pca_results))
Importance of components:
                            PC1      PC2      PC3       PC4       PC5       PC6       PC7       PC8
Standard deviation     461.1061 269.7775 233.1360 203.65763 189.44997 185.91455 183.92640 4.542e-12
Proportion of Variance   0.4379   0.1499   0.1119   0.08543   0.07392   0.07119   0.06968 0.000e+00
Cumulative Proportion    0.4379   0.5878   0.6998   0.78521   0.85913   0.93032   1.00000 1.000e+00
plot(pca_results, col="pink",main="Scree plot")
pca_results$x
                         PC1        PC2        PC3        PC4        PC5         PC6         PC7
5775278051_R01C01 -490.39517  273.66838  -46.53628   11.80547 -132.53730    9.079777 -340.629412
5775278051_R04C02 -580.39316  320.07385  103.43103   29.52610  120.37819   42.834122  288.382072
5775278078_R02C01  -55.67497 -195.75628  -13.04294 -164.58776 -376.02822  -92.090704  158.506414
5775278078_R05C01 -213.07150 -284.22106  -35.89420   39.24178  211.49176 -343.417885  -52.793862
5775278078_R05C02 -197.45533 -391.14454 -148.14290  124.19658   58.87745  320.053770  -11.075602
5930514034_R01C02  589.54629  228.41669 -418.95788 -114.00845   81.24534   -2.236843   39.704411
5930514035_R04C02  590.14977   82.13768  200.98379  375.69546  -87.35420  -34.241623    7.064962
5930514035_R06C02  357.29407  -33.17472  358.15937 -301.86919  123.92698  100.019385  -89.158983
                            PC8
5775278051_R01C01  5.818104e-12
5775278051_R04C02  1.962048e-12
5775278078_R02C01  1.660567e-12
5775278078_R05C01  1.868676e-12
5775278078_R05C02 -1.768683e-11
5930514034_R01C02  6.477131e-12
5930514035_R04C02  5.833994e-12
5930514035_R06C02 -3.655050e-12
```

![scree plot](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/Scree%20plot.png)

It is possible to visualize the results in the plot by colouring the dots according to a particular phenotype, i.e. the group (either WT or DS). As it can be observed, WT and DS groups are well-separated and a fitting line can be plotted to better evaluate the separation.

```
pheno <- read.csv("/Input_data/Samplesheet_report_2020.csv",header=T, stringsAsFactors=T)
palette(c("red","green"))
plot(pca_results$x[,1], pca_results$x[,2],cex=2.5,pch=17,col=pheno$Group,xlab="PC1",ylab="PC2",xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2],labels=pheno$Array,cex=1.0,pos=3)
legend("bottomright",legend=levels(pheno$Group),col=c(1:nlevels(pheno$Group)),pch=17)
```

![pca plot](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/pca_group.png)

### Step 9: Using the matrix of normalized beta values generated in step 7, identify differentially methylated probes between group DS and group WT using the Mann-Whitney test. 
Let's identify differentialy methylated probes between group DS and group WT using a non-parametric test called **Mann-Witheny** or **Wilcoxon rank-sum test**. It is the equivalent of the most-known Student's t-test for unpaired data when no assumptions regarding normal distribution of data are made. We will use **wilcox.test** function for each row of the dataframe of preprocessSWAN beta values; thus an *ad hoc* function is created and applied to preprocessSWAN beta values.

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
When the number of probes is in the order of thousands, a correction of the signficance threshold is required: it is called **multiple test correction**. As we can see, by applying a threshold of 0.05 to the raw pValues of Mann-Whitney test, a high number of false positives has been detected. In this step we are going to apply multiple test correction, such as **Bonferroni correction** and **Benjamini-Hochberg correction**, and to compare the results with the ones from the raw test. Both multiple test corrections have been performed using **p.adjust** function.

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
**Heatpmap** is a data visualization technique showing the level of differentially expressed genes in a matrix with different colours. Heatmaps are always coupled with **hierarchical clustering**, whose aim is to link genes or samples with similar profiles to generate a dendrogram. Different linkage methods can be used to calculate the distance among clusters: (1) **single linkage**, (2) **complete linkage** and (3) **average linkage**. They are all implemented in **heatmap.2** function. Let's first generate the input for the heatmaps and the color bar according to the phenotype group.

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
![single linkage](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/single_link_or.png)

#### Average linkage
![average linkage](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/average_link_ord.png)

Dendrogram generated using complete linkage method is charatcterized by compact and well-defined clusters both for the probes and for the samples. While dendrograms generated with single linkage method show the chaining effect, expecially in the probes' one. However, all heatmaps exhibit a clear division between WT samples (blue) and DS samples (red), as well as a distiction between hypermethylated (high expression level) and ipomethylated (low expression level) CpG islands. 

### Step 12: Produce a volcano plot and a Manhattan plot of the results of differential methylation analysis.
A **volcano plot** is a type of scatterplot that displays the statistical significance against the fold change; it can be coupled with Gene Ontology enrichment. The fold change is calculated as the difference between the average beta values of WT samples and the average beta values of DS samples, in the following way.

```
beta_group <- final_wilcox_corrected[,1:8]
beta_groupDS <- beta_group[,SampleSheet$Group=="DS"]
mean_beta_groupDS <- apply(beta_groupDS,1,mean)
beta_groupWT <- beta_group[,SampleSheet$Group=="WT"]
mean_beta_groupWT <- apply(beta_groupWT,1,mean)
delta <- mean_beta_groupWT - mean_beta_groupDS
```

A data frame is created storing in one column, the delta values, and on the other, the -log10 of pValues; these values are then plotted using **plot** function with some options.

```
toVolcPlot <- data.frame(delta, -log10(final_wilcox_corrected$pValues_wilcox))
plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5)
```

![volcano plot](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/volcano_plot.png)

Since in this case the pValues have been generated using a non-parametric test, the typical "volcano" shape of the plot is less detectable than in a volcano plot generated using pValues from parametric tests. 

Together with the volcano plot, a **Manhattan plot** can be generated using the **gap** package. A Manhattan plot is a type of scatter plot useful for the visualization of large number of data, that are plotted according to the genomic coordinates (X-axis) and the negative log10 of the pValues (Y-axis). For the same reason of volcano plot, the typical profile of skyscrapers charterizing the Manhattan plots is lost. 
After downloading the required package, the dataframe of CpG probes with associted pValues has to be annotated by extracting genomic information from the Illumina450Manifest (cleaned version). The **merge** function is used to perform this step. However, since CpG probes are stored in the rownames and not in the columns (as required by merge() funtion), a new dataframe storing the CpG probe IDs on the column is generated.

```
final_wilcox_corrected_col <- data.frame(rownames(final_wilcox_corrected), final_wilcox_corrected)
head(final_wilcox_corrected_col)
colnames(final_wilcox_corrected_col) [1] <- "IlmnID"
final_wilcox_corrected_col_annotated <- merge(final_wilcox_corrected_col, Illumina450Manifest_clean, by="IlmnID")
dim(final_wilcox_corrected_col_annotated)
[1] 485512     44
```

Then the Manhattan plot is created using **mhtplot** function.

```
input_Manhattan <- data.frame(final_wilcox_corrected_col_annotated$CHR, final_wilcox_corrected_col_annotated$MAPINFO, final_wilcox_corrected_col_annotated$pValues_wilcox)
dim(input_Manhattan)

levels(input_Manhattan$final_wilcox_corrected_col_annotated.CHR)

input_Manhattan$final_wilcox_corrected_col_annotated.CHR <- factor(input_Manhattan$final_wilcox_corrected_col_annotated.CHR,levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))
levels(input_Manhattan$final_wilcox_corrected_col_annotated.CHR)

palette <-rainbow(24)
palette

mhtplot(input_Manhattan, control=mht.control(colors=palette))
axis(2,cex=0.5)
```

![manhattan plot](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/manhattan_plot.png)

### Optional: As DS is caused by the trisomy of chromosome 21, try also to plot the density of the methylation values of the probes mapping on chromosome 21. Do you see a very clear difference between the samples? How many differentially methylated probes do you find on chromosome 21?

```
chr21 <- Illumina450Manifest_clean[Illumina450Manifest_clean$CHR==21,]
chr21<-droplevels(chr21)

#beta_wtset and beta_dsset have been already calculated in step 6
beta_wtset_21 <- beta_wtset[rownames(beta_wtset) %in% chr21$IlmnID,]
beta_dsset_21 <- beta_dsset[rownames(beta_dsset) %in% chr21$IlmnID,]
M_wtset_21 <- M_wtset[rownames(M_wtset) %in% chr21$IlmnID,]
M_dsset_21 <- M_dsset[rownames(M_dsset) %in% chr21$IlmnID,]
mean_of_beta_wtset_21 <- apply(beta_wtset_21,1,mean,na.rm=T)
mean_of_beta_dsset_21 <- apply(beta_dsset_21,1,mean,na.rm=T)

d_mean_of_beta_wtset_21 <- density(mean_of_beta_wtset_21)
d_mean_of_beta_dsset_21 <- density(mean_of_beta_dsset_21)

mean_of_M_wtset_21 <- apply(M_wtset_21,1,mean,na.rm=T)
mean_of_M_dsset_21 <- apply(M_dsset_21,1,mean,na.rm=T)
d_mean_of_M_wtset_21 <- density(mean_of_M_wtset_21)
d_mean_of_M_dsset_21 <- density(mean_of_M_dsset_21)
```
Then the plots are generated as overlapping plots.

```
par(mfrow=c(1,2))
plot(d_mean_of_beta_wtset_21,main="Density of Beta Values",col="orange")
lines(d_mean_of_beta_dsset_21,col="purple")
plot(d_mean_of_M_wtset_21,main="Density of M Values",col="orange")
lines(d_mean_of_M_dsset_21,col="purple")
```
![optional plots](https://github.com/giorgiagandolfi/DNA-RNA_Dynamics/blob/master/optional_chr21.png)

No great differences in the distribution of beta and M values between the two groups can be appreciated.

In order to find how many differentially methylated probes are in chromosome 21, from the vector storing all pValues lower than the significance threshold, the ones associated with CpG probes in chromosme 21 are extracted. 

```
all_chromosomes<-data.frame(Illumina450Manifest_clean$IlmnID,Illumina450Manifest_clean$CHR)
only_chr21<-ID_CHR[ID_CHR$Illumina450Manifest_clean.CHR==21,]
extract_pVal_chr21<-rownames(final_wilcox)%in%cpg21$Illumina450Manifest_clean.IlmnID
finaL_wilcox_0.05_chr21<-final_wilcox[final_wilcox$pValues_wilcox<=0.05 & extract_pVal_chr21,]
dim(finaL_wilcox_0.05_chr21)
[1] 298   9
```




