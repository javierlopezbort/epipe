# EPIPE: (E)pigenetic (PIPE)line

EPIPE is a bioinformatics pipeline designed for the analysis of methylation data obtained from microarrays. It implements the whole methylation array workflow steps from processing raw data to pathway annotation of the differentially methylation positions (DMP) and regions (DMR) found.

The pipeline is set up to use `targets` to manage the execution of each step in the analysis process.

## Installation

You can install the development version of epipe with:

``` r
# install.packages("devtools")
devtools::install_github("ijcBIT/epipe")
```

## Configuration

### Input data requirements

EPIPE requires the following input files:

-   **Idats/**: A folder containing all the raw idat files for the methylation analysis.

-   **Samplesheet.rds**: A sample sheet in RDS format.

Before running the pipeline, ensure that the sample sheet adheres to the following guidelines:

-   Must be saved as an RDS object.

-   Include a Basename column that specifies the path to the idats files (e.g., inst/extdata/EPICv2/idats/207505140099_R08C01).

-   All idats files must be located within the same folder.

-   Avoid spaces in column names.

-   A `Sample_Name` column should be included, containing unique identifiers for each sample.

-   Must include a `Sentrix_ID` column.

For an example, refer to the provided sample sheet:

``` r
data(samplesheeet)
```

### Configuration file: config.R

EPIPE uses a `config.R` file to define key parameters for the pipeline. You can customize the configuration file according to your analysis needs.

Access the configuration file via the following:

``` r
# Path to the config.R file
config_file_path <- system.file("config.R", package = "epipe")

# Open the file in RStudio
file.edit(config_file_path)
```
Once you have the file save it and make a copy.


### Targets.R file

When you install or load the EPIPE package with `library(epipe)`, the \_targets.R file will be automatically copied from the package to your current working directory (if it is not already present). This allows you to immediately start running your bioinformatics pipeline without needing to manually copy the file.

However, you can also manually copy the \_targets.R file to your working directory, by using the following R command:

``` r
file.copy(system.file("_targets.R", package = "epipe"),"_targets.R")
```

## Running the pipeline

Once the `config.R` file is configured, you can run the pipeline as follows:

``` r
library(epipe)
targets::tar_make()
```

This will start the pipeline based on the targets defined in the config.R file.

## Output

The primary outputs are HTML reports designed for different audiences:

1.  **Collaborator Report**: A concise overview highlighting key findings and summaries.
2.  **Analyst Report**: A detailed report containing comprehensive data, visualizations, and in-depth analysis.

#### Report Contents

Reports include:

-   Quality control plots.
-   Methylation Values: Raw and normalized beta values in tabular formats.
-   Exploratory analysis: PCA plots, heatmaps, and sex-specific visualizations.
-   DMPs: Lists, summaries, and plots (e.g., volcano and Manhattan plots).
-   DMRs: Detailed region analysis results and summaries.
-   Pathway enrichment: Results for DMPs and DMRs.


## Examples

To help you get started with the pipeline, we provide some example data in the inst/exdata/ folder of the epipe package.

-   Idats folder:

``` r
system.file("extdata/EPICv2/idats/", package = "epipe")
```

-   Sample sheet

``` r
system.file("extdata/EPICv2/samplesheet_EPICv2.rds", package = "epipe")
```

These files serve as input data for the default pipeline.

Use the default config.R file and run the pipeline in R:

``` r
library(epipe)
targets:tar_make()
```

