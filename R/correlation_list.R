#' Correlation Analysis
#'
#' @description This function performs correlation analysis on a Rgset object. It identifies correlated variables, and generates plots.
#'
#' @param clean_object A Rgset object.
#' @param path The path to save plots.
#' @param variables A vector of variable names to include in the analysis, or "ALL" to include all variables.
#' @param sg A grouping variable to color the correlation plot.
#'
#' @return A list of correlated variable pairs.
#'
#' @import ggplot2
#' @import GGally
#' @importFrom stats cor.test t.test aov chisq.test prop.test
#' @import data.table
#' @export

# path<-custom_paths$qc_folder
# x<-correlation_analysis(clean_EPICv2,path=custom_paths$qc_folder)


correlation_analysis<-function(clean_object,path,variables,sg){
  data_prepared<-preparation(clean_object)
  dataset_df<-dataset_var(data_prepared,variables)
  listed_variables<-correlated_list(dataset_df)
  correlation_plot(listed_variables,dataset_df,sg,path)
  return(listed_variables)
}

## EXAMPLES:

#x<-correlation_analysis(clean_EPICv2)
#
# data_prepared<-preparation(clean_EPICv2)
# listed_variables<-correlated_list(data_prepared)
# correlation_plot(listed_variables,data_prepared)



#############################################################################################

#' Sample sheet check and preparation
#'
#' @description Check and prepare the sample sheet for correlation analysis by filtering and processing relevant metadata variables.
#'
#' @param object RgSet object that includes metadata
#'
#' @import S4Vectors
#' @import Biobase
#' @import data.table
#'
#' @return A data.table containing the processed and relevant variables for correlation analysis.

# Punt de partida-->clean object (includes sex variable)

preparation<-function(object){

  # Just select the variables that may be correlated (not all the variables that we have)
  ###### Ex: Sample_name is irrelevant

  #metadata<-colData(object)
  metadata = object@colData

  df<-data.frame(metadata)
  ss_new<-get_category(df)

  variables<-names(ss_new)
  category<-attributes(ss_new)$category

  # Select variables that are not 'ids'
  non_ids_variables <- variables[category != "ids"]

  # Subset the dataframe with selected variables
  ss_new_subset <- ss_new[, non_ids_variables, with = FALSE]

  ##################################################################

  # Transform character variables into a factor:
  for (col in names(ss_new_subset)) {
    if (is.character(ss_new_subset[[col]])) {
      ss_new_subset[[col]] <- as.factor(ss_new_subset[[col]])
    }
  }

  ###############################################################

  # Can happen that some variables have more levels than what they truly have (ex: Type variable: covs, case,control)-->covs is not a level


  # Iterate over each column in the data.table
  for (col in names(ss_new_subset)) {
    # Check if the column is a factor and if "covs" is one of its levels
    if (is.factor(ss_new_subset[[col]]) && "covs" %in% levels(ss_new_subset[[col]])) {
      # Remove the level "covs" from the factor variable
      ss_new_subset[[col]] <- droplevels(ss_new_subset[[col]], exclude = "covs")
    }
    if (is.factor(ss_new_subset[[col]]) && "batch" %in% levels(ss_new_subset[[col]])) {
      # Remove the level "covs" from the factor variable
      ss_new_subset[[col]] <- droplevels(ss_new_subset[[col]], exclude = "batch")
    }
  }

  return(ss_new_subset)

}



#############################################################################################################################

#' Select Variables from the Dataset
#'
#' @description This function selects specified variables from a dataset. If no variables are specified, the full dataset is returned.
#'
#' @param dataset A data.table or data.frame containing the dataset.
#' @param variables A vector of variable names to select. If set to `'ALL'`, the full dataset is returned.
#'
#' @return A subset of the dataset containing only the specified variables, or the full dataset if `'ALL'` is specified.
#' @import data.table

dataset_var<-function(dataset,variables=NULL){

  if (length(variables) == 1 && variables == 'ALL') {
    # Do nothing, use the full dataset
  } else {
    dataset <- dataset[, ..variables]
  }

  return(dataset)

}



#####################################################################################################################

#' List of correlated variables
#'
#' @description Obtain a list of variables that are correlated or its design is not well balanced, depending on each type of variables(numerical/categorical)
#'
#' @param dataset A data.frame or data.table containing the data to be analyzed.
#'
#' @return A list where each element contains a pair of correlated variables and their correlation score or p-value.
#' @import stats

correlated_list<-function(dataset){

  cor_data<-list()
  p.values<-NULL
  # Loop over each column in the data.table
  for (i in 1:length(names(dataset))) {
    variable1 <- names(dataset)[i]  # Get the name of the first variable

    # Start the inner loop from the next column to enforce the order
    for (j in (i + 1):length(names(dataset))) {
      variable2 <- names(dataset)[j]  # Get the name of the second variable


      ##################################################################
      # # If both variables are numeric--> Pearson correlation test
      if (is.numeric(dataset[[variable1]]) & is.numeric(dataset[[variable2]]) & variable1!=variable2){

        cor_score<-cor.test(x=dataset[[variable1]],y=dataset[[variable2]])$estimate
        #print(cor_score)
        if (abs(cor_score)>0.5){
          cor_data <- append(cor_data, list(c(variable1, variable2,cor_score)))
          #cor_data<-c(cor_data,paste0(variable1,' and ', variable2,' are correlated'))
          #p.values<-c(p.values,cor_score)
        }

      }

      ###########################################################
      # #If variable 1 is numerical and variable 2 is caregorical
      if (is.numeric(dataset[[variable1]]) & is.factor(dataset[[variable2]])){

        # If categorical variable has just one group---> NO fit a regression model
        # If categorical variable has two groups --> t.test
        # If categorical variable has more than 2 groups --> anova test
        if (length(levels(dataset[[variable2]]))==2){
          p.value<-t.test(dataset[[variable1]]~dataset[[variable2]])$p.value
          if (p.value<0.05){
            cor_data <- append(cor_data, list(c(variable1, variable2,p.value)))
            #cor_data<-c(cor_data,paste0(variable1,'and', variable2, 'design is not balanced'))
            #p.values<-c(p.values,p.value)
          }

        }else if(length(levels(dataset[[variable2]]))>2){
          anova<-aov(dataset[[variable1]]~dataset[[variable2]])
          p.value<-summary(anova)[[1]]$`Pr(>F)`[1]
          if (p.value<0.05){
            cor_data <- append(cor_data, list(c(variable1, variable2,p.value)))
            #cor_data<-c(cor_data,paste0(variable1,'and', variable2, 'design is not balanced'))
            #p.values<-c(p.values,p.value)
          }

        }


      }

      #############################################################
      # #If variable 1 is categorical and variable 2 is numerical

      if (is.numeric(dataset[[variable2]]) & is.factor(dataset[[variable1]])){

        # If categorical variable has just one group---> NO fit a regression model
        # If categorical variable has two groups --> t.test
        # If categorical variable has more than 2 groups --> anova test
        if (length(levels(dataset[[variable1]]))==2){
          p.value<-t.test(dataset[[variable2]]~dataset[[variable1]])$p.value
          if (p.value<0.05){
            cor_data <- append(cor_data, list(c(variable1, variable2,p.value)))
            #cor_data<-c(cor_data,paste0(variable1,'and', variable2, 'design is not balanced'))
            #p.values<-c(p.values,p.value)
          }


        }else if (length(levels(dataset[[variable1]]))>2){
          anova<-aov(dataset[[variable2]]~dataset[[variable1]])
          p.value<-summary(anova)[[1]]$`Pr(>F)`[1]
          if (p.value<0.05){
            cor_data <- append(cor_data, list(c(variable1, variable2,p.value)))
            #cor_data<-c(cor_data,paste0(variable1,'and', variable2, 'design is not balanced'))
            #p.values<-c(p.values,p.value)
          }

        }


      }


      ##########################################################################
      #If both variables are categorical---> chisq.test (multinomial) or test.equal proportions (binomial)


       if (is.factor(dataset[[variable1]]) & is.factor(dataset[[variable2]]) & variable1!=variable2){


         # if both variables have 2 levels --> test equal proportions

         if (length(levels(dataset[[variable1]]))==2 & length(levels(dataset[[variable2]]))==2){
           contingency_table<-table(dataset[[variable1]],dataset[[variable2]])
           p.value<-prop.test(contingency_table)$p.value
           if (p.value<0.05){
             cor_data <- append(cor_data, list(c(variable1, variable2,p.value)))
                 #p.values<-c(p.values,p.value)
               }

         }


         # if at least one variable has more than 2 levels --> chisq.test

         if (length(levels(dataset[[variable1]])) > 2 | length(levels(dataset[[variable2]])) > 2) {

           # Perform chi-square test of independence
           contingency_table <- table(dataset[[variable1]], dataset[[variable2]])
           chisq_test_result <- chisq.test(contingency_table)
           p.value <- chisq_test_result$p.value
           if (p.value < 0.05) {
             cor_data <- append(cor_data, list(c(variable1, variable2,p.value)))
             #p.values <- c(p.values, p.value)
           }

         }




       }


    }

  }

  return(cor_data)

}



# How to use:
#coorrelated_l<-correlated_list(tips)
