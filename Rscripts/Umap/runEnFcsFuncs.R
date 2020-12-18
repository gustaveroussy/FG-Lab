#Evan Newell 2015; Ming 2016; Etienne Becht 2016; Evan 2017, 2018

## New Names.csv file sytem:  using y to indicate channels to use still works but also indicates to use default cytof transform
# Now if "f" is indicated default (old version) fluor logicle params will be used (not recommended), "a" instructs auto-logicle transform (recommended for fluor data)
# in addition, "c" also means cytof transform, "l" linear normalized to max = 4.5, "n" linear normalized min=0, max=4.5, "as5" ArcSinh x/5 (for cytof), "as150" ArcSinh x/150 (for fluor)
# Additional transform options can be implemented with this sytem.

# Startup Settings:
if(!require(devtools)){
  install.packages("devtools")
}
library(devtools) # A workaround to Rtools incompatibility with more recent versions of R.
assignInNamespace("version_info", c(devtools:::version_info, list("3.5" = list(version_min = "3.3.0", version_max = "99.99.99", path = "bin"))), "devtools")

# Parameters are split into basic and advanced.
## Basic Parameters:

##### Choose functions to run-------------------------------------------------------------------
RunAlgorithms = T   # change to F if you just want to make meaning plots or heatplots using output files (no need to rerun tSNE etc.)
RunMeaningPlot = F
RunFileHeatplots = F
    #IF TRUE  ... In progress
    DoClusters=F  # makes heatmaps to describe phenograph clusters - median marker intensity and frequencies (lin and log) within samples
    clustParam="Phenograph"  # e.g., "Phenograph" or "FlowSOM"
Run2DPlot = F
Run3DPlot = F

##### Main Parameters --------------------------------------------------------------------------
sourceFcsFolder = "fcs"      # Folder with samples (even 1 sample ok)
ceil = 200000                  # Number of events to take PER SAMPLE (unless sample has less)
FNnames="names2.csv" # Parameters to analyze in csv file
OutputSuffix = "p30_md0.01_k30out"         # Adds this to all new directories and file names


##### Algorithms -------------------------------------------------------------------------------
DotSNE = T
    tSNEperpelxity = 30      #default = 30; increased perplexity => increased spread
DoUMAP = T
    min_dist = 0.01           #default = 0.2
    n_neighbors = 15         #default = 15
Do3DtSNE = F                 #T for running 3D tSNE
Do3DUMAP = F
Dophenograph = T             #clusters cells using Rphenograpy. 
    kValue = 30              #default = 30; k is low => more clusters and inversely

###### Meaningplot Parameters ------------------------------------------------------------------
Xaxis = "UMAP1"   ## examples: tSNE1 or UMAP1  Add * to color by an additional round of tSNE or UMAP
Yaxis = "UMAP2"  ## examples: tSNE2 or UMAP2
prefix = paste0(sourceFcsFolder,OutputSuffix,"_UMAP") #edit that last part to give this plot a name

#### 2D Plot Parameters ------------------------------------------------------------------------
TwoDColorBy = "Phenograph"   # e.g. "FlowSOM", "Phenograph", use "InFile" to color by sample file
TwoDparN1 = "UMAP1"          # e.g., UMAP1 or tSNE1
TwoDparN2 = "UMAP2"          # e.g., UMAP2 or tSNE
TwoDlabelClusters = T
TwoDOutputSuffix = "tSNE"

#### 3D Plot Parameters ------------------------------------------------------------------------ 
Color3Dby = "Phenograph"     # e.g. "FlowSOM", "Phenograph", use "InFile" to color by sample file
parN1 = "TD_UMAP1"           # e.g., TD_UMAP1 or ThreeDtSNE1
parN2 = "TD_UMAP2"           # e.g., TD_UMAP2 or ThreeDtSNE2
parN3 = "TD_UMAP3"           # e.g., TD_UMAP3 or ThreeDtSNE3
labelClusters = T
TDOutputSuffix = "3D-tSNE"



#######  Advanced #######  #######  Advanced #######  

#### General Settings ---------------------------------------------------------------------------
useCSV = T   # T if using .csv files instead of fcs files - will make csv output also
TransformOutput = F  # F Puts derived data on linear scale 1-10k - best for viewing in FJ10 otherwise T is setup for FJ9. Talk to Evan about implications for one-sense

#### Additional algorithms and settings ---------------------------------------------------------
DoOneSENSE = F  #needs a modified names.csv file with column names for each category
DoOneSUMAP = F # umap version of One-SENSE 
DoFlowSOM = F 
MaxClusters = 30
DoIsomap = F
DoDiffMap = F
DoPCA = F

#### Meaningplot Parameters ---------------------------------------------------------------------
MeaningTopPercentile = 1
MeaningBotPercentile = 0
DoOneForEach = F #T => will do meaning plots for each individual file, F => for the concat file
palette=c("black","blue","green","yellow","red")
color.scale.type="relative"  # choose "relative" or "absolute"
resolution=150
cex=0.2
pch=16
sourceFcsForMeaningPlot <- paste0(sourceFcsFolder,"_",OutputSuffix)

##### FileHeatplot Settings ---------------------------------------------------------------------
fileHPoutputsuffix = paste0(OutputSuffix,"HP") # Edit last part for unique suffix
HPpalette=c("black","blue","lightblue","green","yellow","darkorange","darkred")

##### 2D Plot Settings --------------------------------------------------------------------------
TwoDNcolors = 0  # set to zero to have this automatically set to maximum number
TwoDpch = 16
TwoDcex = 0.25
TwoDresoution = 600
if (!require(RColorBrewer)) {
  install.packages("RColorBrewer") # Install Rtsne package from CRAN
  library(RColorBrewer) # Load package
}
library(RColorBrewer) 
# code for making random colors, good for clustering
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
paletteFor2D <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#same as for meaning plots
#paletteFor2D <- c("black","blue","green","yellow","red")
TwoDColorList = paletteFor2D # can choose colors if numbers match otherwise it will interpolate
# ColorList = c("colorX", "colorY", ...) put n color names if n clusters in order to manually define the colors of the clusters
# Color palette in R: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf)

#### 3D Plot Settings --------------------------------------------------------------------------
Ncolors = 0  # set to zero to have this automatically set to maximum number
# code for making random colors, good for clustering
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
paletteFor3D <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ColorList = paletteFor3D # can choose colors if numbers match otherwise it will interpolate
# ColorList = c("colorX", "colorY", ...) put n color names if n clusters in order to manually define the colors of the clusters
# Color palette in R: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf)





###################################### EXECUTE FUNCTIONS ###################################################

##### ##### No edit below here ## No edit below here ## No edit below here #### ####

source('ENfcsFuncs.R')

if (RunAlgorithms) 
  FCStSNEone(LoaderPATH = sourceFcsFolder, 
           useCSV = useCSV, # T if using .csv files instead of fcs files - will make csv output also
           TransformOutput = TransformOutput, # F Puts derived data on linear scale 1-10k - best for viewing in FJ10 otherwise T is setup for FJ9. Talk to Evan about implications for one-sense
           ceil = ceil, #number of events to take per sample (unless sample has less)
           FNnames=FNnames,#Parameters to analyze in csv file
           OutputSuffix = OutputSuffix, # Adds this to all new directories and file names
           DotSNE = DotSNE,
           tSNEperpelxity = tSNEperpelxity, #default = 30; increased perplexity => increased spread
           DoOneSENSE = DoOneSENSE, #needs a modified names.csv file with column names for each category
           Dophenograph = Dophenograph, #clusters cells using Rphenograpy
           kValue = kValue, #default = 30; k is low => more clusters and inversely
           DoFlowSOM = DoFlowSOM,
           MaxClusters = MaxClusters,
           DoIsomap = DoIsomap,
           DoDiffMap = DoDiffMap,
           DoPCA = DoPCA, 
           Do3DtSNE = Do3DtSNE, #T for running 3D tSNE
           DoUMAP = DoUMAP, #still in prep- will add many parameters for this
           Do3DUMAP = Do3DUMAP,
           min_dist = min_dist,
           n_neighbors = n_neighbors,
           DoOneSUMAP = DoOneSUMAP) 


## Still in prep:
if (RunMeaningPlot)
   meaningPlot(LoaderPATH =sourceFcsForMeaningPlot,
            useCSV = useCSV, # T if using .csv files instead of fcs files - will make csv output also
            FNnames = FNnames,
            ceil = ceil,
            TransformOutput = TransformOutput,
            MeaningTopPercentile = MeaningTopPercentile,
            MeaningBotPercentile = MeaningBotPercentile,
            PC1 = Xaxis,
            PC2 = Yaxis,
            DoOneForEach = DoOneForEach,
            prefix = prefix,
            palette=palette,
            color.scale.type=color.scale.type,
            resolution=resolution,
            cex=cex,
            pch=pch)

if(RunFileHeatplots)
  fileHeatplot(LoaderPATH =sourceFcsForMeaningPlot,
                          ceil = ceil,
                          TransformOutput = TransformOutput,
                          FNnames = FNnames,
                          OutputSuffix = fileHPoutputsuffix,
                          palette=HPpalette,
                         DoClusters=DoClusters,
                         clustParam=clustParam)

if(Run2DPlot)
  TwoDPlot (LoaderPATH =sourceFcsForMeaningPlot,
              FNnames = FNnames,
              ceil = ceil, 
              OutputSuffix = TwoDOutputSuffix,
              parN1 = TwoDparN1,
              parN2 = TwoDparN2, 
              Ncolors = TwoDNcolors,
              colparN = TwoDColorBy,
              labelClusters = TwoDlabelClusters,
              pch = TwoDpch,
              cex = TwoDcex,
              resolution = TwoDresoution,
              ColorList = TwoDColorList)  # can use "InFile" to color by sample file

if(Run3DPlot)
  ThreeDPlot (LoaderPATH =sourceFcsForMeaningPlot,
              FNnames = FNnames,
              ceil = ceil, 
              OutputSuffix = TDOutputSuffix,
              parN1 = parN1,
              parN2 = parN2, 
              parN3= parN3, 
              labelClusters = labelClusters,
              Ncolors = Ncolors,
              colparN = Color3Dby,
              ColorList = ColorList)  # can use "InFile" to color by sample file

print("Completed without errors!")

#gatingDataPoints (cropping = F, useUMAP = T)

##### Useful script for making names.csv file:
# if (!require(flowCore)) { 
#   source("http://bioconductor.org/biocLite.R")
#   biocLite("flowCore")
# } 
# 
# inFileN <- "test.fcs"        #FCS file to make names file for 
# outNamesCsv <- "names.csv"   #Output csv file name 
# 
# FF <- read.FCS(inFileN)
# colNames <- FF@parameters$desc
# empties <- which(is.na(colNames) | colNames== " ")
# colNames[empties] <- FF@parameters$name[empties]
# write.csv(colNames, outNamesCsv, row.names = F)
