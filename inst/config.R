############   EPIPE CONFIGURATION FILE   ####################
##############################################################

## Project configuration:
project_name<-'PIK3CA-METH3'
author <- "Marina Vilardell"
description <- 'Dna methylation analysis in PIK3CA mutations'


################################################################################

## Sample sheet path:
data_paths<-c(example_EPICv2=system.file("extdata/EPICv2/samplesheet_EPICv2.rds", package = "epipe")) # Substitute system.file for the path to your samplesheet
data_names <- c(names(data_paths))

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
sampGroups='Condition' # Column name indicating Sample group or the variable of interest


## Normalization:
# Choose a method for normalization:
# noob, ssnoob, swan, funn, noob_pq, noob_swan, quantile
norm_function <-'quantile'


## Filters:
remove_sex <- TRUE
sex_prediction<-TRUE
#crossreactive=TRUE
#snps=TRUE


################################################################################
### Color palettes:

pal_discrete =  c("#1B9E77","#7570B3","#D95F02", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666","#FF0010")

################################################################################
### Deconvolution:

# Deconvolution will predict the following cell types proportion:
# Leukocyte (leuk.prop), neutrophil (Neu), Monocytes (Mono), CD8T, CD4T, Natural Killer (NK), Bcell


#######################################################################################################
### PCA:

## PCA top variable:
topN <- c(100,1000,5000,10000)

# Variables to put in the correlation analysis:
# By default: All variables
# variables='ALL'
# But you can specify the variables that you are interested in.

variables=c('Sentrix_ID','mutation','condition','age','predictedSex','purity','CD4T','NK','Bcell','Mono','Neu','leuk.prop')


# Variables to select for the Bplots
bplots_var<-c('Sentrix_ID','Sentrix_position','mutation','condition','purity','predictedSex','predictedAge',
              'CD4T','NK','Bcell','Mono','Neu','leuk.prop','age','sex')


########################################################################################
#### Model:

group_var = "condition" # variable of interest that we want to model.

# Define the contrast we are interested in. If NULL, the function will automatically generate them.
#Contrasts<-NULL
Contrasts <- 'lymphatic-control'

# Covariates:
# covs<-c('Neu','predictedSex','age')
# covs<-c('predictedSex','age','mutation')
# covs=NULL
covs='predictedSex'
# covs=c('Gestational_obesity_group','Pregestational.obesity_group','Sex')


########################################################################################
##### Differential methylation analysis:

# Parameters for DMP:
mDiffDMP = 0.01
adjp.valueDMP =0.05
p.valueDMP = 0.05

# Parameters for DMR:
min.cpgDMR = 3
fdrDMR = 0.01
mdiffDMR = 0.05


# DO NOT MODIFY:
values <- data.table(cbind(data_names,data_paths,arraytype,project_name,mDiffDMP,adjp.valueDMP,p.valueDMP,fdrDMR,mdiffDMR,min.cpgDMR))
values$mDiffDMP=as.numeric(values$mDiffDMP)
values$p.valueDMP=as.numeric(values$p.valueDMP)
values$adjp.valueDMP=as.numeric(values$adjp.valueDMP)
values$min.cpgDMR=as.numeric(values$min.cpgDMR)
values$fdrDMR=as.numeric(values$fdrDMR)
values$mdiffDMR=as.numeric(values$mdiffDMR)

values <- values[,.(norm=rlang::syms(norm_function)),by=data_names][values,on=.(data_names)]





### Advanced functions:


# models (In case a more specific model is desired). Overrides group_var & covs:

# A) Specific model for each sample_sheet:
# models<- c(
#   m1 = "~ 0 + Condition",
#   m2 = "~ 0 + Type"
# )
# values <- values$models<-models

# B) Multiple models for each sample_sheet:
# models<- c(
#   m1 = "~ 0 + Type * Condition",
#   m2 = "~ 0 + Condition",
#   m3 = "~ 0 + Type"
# )
# values <- values[,.(models=models),by=data_names][values,on=.(data_names)]

# Singular is one against all the rest of variables in Sample_Group
# Pairwise is all against all
# gr: groups together multiple values making the mean.
# Subset: Type == Control, subset the sample_sheet so only a subgroup is left


### Raw data paths (Genomics Unit):
path_to_Samples_on_array_excel = "W://GENOMICS_UNIT/SampleDB/Samples_on_array_last.xls"
path_to_SDNA.xls = "W://GENOMICS_UNIT/DBB/SDNA.xls"
path_to_storage = "/ijc/LABS/GENOMICS/RAW/Arrays/"
path_to_LTS = "/ijc/LABS/GENOMICS/LTS"



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

report_opts <- rbind(report_colab,report_analyst)

report_parameters <- tibble::tibble(
  report_opts,
  path_to_Samples_on_array_excel = path_to_Samples_on_array_excel,
  path_to_SDNA.xls = path_to_SDNA.xls,
  path_to_storage = path_to_storage,
  path_to_LTS = path_to_LTS
)

