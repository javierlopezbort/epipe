############   EPIPE CONFIGURATION FILE   ####################
##############################################################

## Project configuration:
project_name<-'PIK3CA-METH3'
author <- "Marina Vilardell"
description <- 'Dna methylation analysis in PIK3CA mutations'


################################################################################

## Sample sheet path:
data_paths<-c(example_EPICv2=system.file("extdata/EPICv2/samplesheet_EPICv2.rds", package = "epipe")) # Substitute system.file for the path to your samplesheet

arraytype <- "EPICv2"

## IDATS folder path:
idats_folder <- system.file("extdata/EPICv2/idats/", package = "epipe")
# NULL if the Basename column of the sample sheet correctly contains the path to idats green and red files

################################################################################

## Folder paths:
results_folder <- "results/"
analysis_folder <- "analysis/"

################################################################################

## Metadata:
idcol= "Sample_Name" #column with unique sample ids.
sampGroups='Sample_Group' # Column name indicating Sample group or the variable of interest


## Normalization:
# Choose a method for normalization:
# noob, ssnoob, swan, funn, noob_pq, noob_swan, quantile
norm_function <-'noob'


## Filters:
remove_sex <- TRUE
sex_prediction<-TRUE
#crossreactive=TRUE
#snps=TRUE


################################################################################

### Color palettes:

pal_discrete =  c("#1B9E77","#7570B3","#D95F02", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666","#FF0010")


################################################################################

## Deconvolution:

# Deconvolution will predict the following cell types proportion:
# Leukocyte (leuk.prop), neutrophil (Neu), Monocytes (Mono), CD8T, CD4T, Natural Killer (NK), Bcell


################################################################################

## Variables for correlation analysis:
# By default: All variables
# variables='ALL'

# But you can specify the variables that you are interested in:
variables=c('Sample_Group','Sentrix_ID','Condition','predictedAge','predictedSex','CD4T','NK','Bcell','Mono','Neu','leuk.prop')


## PCA top variable:
topN <- c(100,1000,5000,10000)

# Variables to select for the Bplots
bplots_var<-c('Sample_Group','Sentrix_ID','Sentrix_Position','Condition','predictedSex','predictedAge',
              'CD4T','NK','Bcell','Mono','Neu','leuk.prop')


################################################################################

## Model:
group_var = "Sample_Group" # variable of interest

# Define the contrast interested in.
#Contrasts <- 'Control-Treated'
Contrasts <- NULL

# Covariates:
# covs<-c('Neu','predictedSex','age')
covs=NULL


########################################################################################

## Differential methylation analysis:

# Parameters for DMP:
mDiffDMP = 0.01
adjp.valueDMP =0.05
p.valueDMP = 0.05

# Parameters for DMR:
min.cpgDMR = 3
fdrDMR = 0.01
mdiffDMR = 0.05



################################################################################################

# Report options:

report_colab <- tibble::tibble(
  report_name = "report_colab",
  qc_include_analyst = "false",
  qc_include_collab = "true",
  corrplot_include = "false",
  parameter_tunning_plots = "false",
  code_fold = TRUE,
  show_warning = FALSE,
  shows_message = FALSE

)

report_analyst <- tibble::tibble(
  report_name = "report_analyst",
  qc_include_analyst = "true",
  qc_include_collab = "false",
  corrplot_include = "true",
  parameter_tunning_plots = "true",
  code_fold = FALSE,
  show_warning = TRUE,
  shows_message = TRUE

)

