# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
# library(epipe)
library(data.table)
library(quarto)
# Load other packages as needed. # nolint
library(crew)
library(crew.cluster)

# Libraries from Illumina
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

max_ncores=RcppParallel::defaultNumThreads() # Not in use but may be useful

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

source("config.R") # Source other scripts as needed. # nolint

# Set the number of workers proportional to the number of rows in values.
# workers are downgraded if idle for 3 seconds.
# controller_local <- crew::crew_controller_local(
#   name = "my_controller",
#   workers =   6,#max_ncores -2, # nrow(values)*3,
#   seconds_idle = 3
# )

# Controller to run with slurm that converts each task to a slurm job
# Needs configuration not working.
# 1. save a singularity image on the HPC
# 2. Find a way to load the singularity image and run the commands inside
# controller_hpc <- crew.cluster::crew_controller_slurm(
#   name = "minastirith_controller",
#   workers = 4,
#   seconds_idle = Inf,
#   script_lines = c("module load singularity","singularity shell /mnt/beegfs/idevillasante/apps/rocker/images/meth-dev.sif"),
#   slurm_log_output = "log_folder/",
#   slurm_memory_gigabytes_per_cpu = 8,
#   host = "minastirith")

# async_controller <- crew_controller_async(workers = 12, processes = 4, host="minastirith")

# controller_hpc$launcher$script("minastirith_controller")

# Set target options:
tar_option_set(
  packages = c("tibble","data.table","S4Vectors"), # packages that your targets need to run
  # format = "rds", # default storage format
  format= "qs", # faster than default rds
  # controller = controller_local,
  storage = "worker",
  retrieval = "worker"
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multiprocess")

# tar_make_future() configuration (okay to leave alone):
future::plan(future.callr::callr)




# Replace the target list below with your own:
top_betas_N_target <- tar_target(top_betas_N,topN)

targets <- tarchetypes::tar_map(
  values = values,
  names = data_names, #"data_source", # Select columns from `values` for target names.
  tar_target(vals,values[data_names,,on="data_names"]),
  tar_target(sample_names_config,idcol),
  tar_target(sample_groups_config,sampGroups),
  #tar_target(slices, values, pattern = slice(data_names)),
  tar_target(custom_paths,make_results_dirs(subf=data_names, results_folder = results_folder, analysis_folder = analysis_folder)),
  tar_target(samplesheet_path, data_paths, format = "file"),
  tar_target(samplesheet, readRDS(samplesheet_path)),
  tar_target(ss, get_category(samplesheet)),

  # Read idats:
  tar_target(rgSet, read_idats(ss,idats_folder=idat_folder,arraytype = arraytype,idcol=idcol, author= author, description = description)), # makes rgChannelset object using minfi and annotates according to arraytype

  # Qc report:
  tar_target(QC_plotss, qc(rgSet,
                          sampGroups = sampGroups, #names(ss)[attributes(ss)$category == "mgroups"][1],
                          sampNames=idcol,
                          idcol=idcol,
                          pal=pal_discrete,
                          qc_folder = custom_paths[["qc_folder"]]
                          )
             ,packages = "minfi",deployment = "worker" ),

  # Filters:  -probes: pval<0.01, -samples: 10% probes fail, Plots: Sample p.values barplot (colMeans)
  tar_target(filtered, filter(
    rgSet=rgSet,
    sg= sg, #names(ss)[attributes(ss)$category == "mgroups"][1],
    qc_folder = custom_paths[["qc_folder"]],
    save_barplot=T)),

  #Tumor purity
  #tar_target(purity, purify(rgSet = rgSet, arraytype = arraytype),error = "continue"),

  #Get raw beta values (before normalization)
  tar_target(raw_beta, minfi::getBeta(filtered)),
  tar_target(save_raw_beta,write.csv(raw_beta, file=file.path(custom_paths$betas_folder,'raw_beta.csv'))),

  #Normalization
  tar_target(normalize, norm(filtered)),
  tar_target(dplot_normalize,denplot(normalize,
                                ss,
                                sampGroups = sampGroups, #names(ss)[attributes(ss)$category == "mgroups"][1],
                                path=custom_paths[["qc_folder"]],
                                norm_method = norm_function)),
  tar_target(normalize_all, normalization_all_methods(filtered,
                                                      ss,
                                                      sampGroups = sampGroups, #names(ss)[attributes(ss)$category == "mgroups"][1],
                                                      path=custom_paths[["qc_folder"]],
                                                      norm_method = norm_function)),

  # Get beta values after normalization
  tar_target(norm_beta, minfi::getBeta(normalize)),
  tar_target(save_norm_beta,write.csv(norm_beta, file=file.path(custom_paths$betas_folder,'norm_beta.csv'))),


  #Preprocessing (probe removal)
  tar_target(clean, prep(normalize, remove_sex = remove_sex, arraytype = arraytype,sexplot_folder= custom_paths[["sexplot_folder"]],predict_sex=sex_prediction)),


  #Predict age
  tar_target(clean_age,ageprediction(clean),error='continue'),

  #Deconvolution
  tar_target(clean_all_deconv,celldeconvolution(rgSet,clean_age,arraytype = arraytype),error = 'continue'),

  # Correlation analysis
  tar_target(cor_analysis, correlation_analysis(clean_all_deconv,
                                                path=custom_paths[["qc_folder"]],variables=variables,sg=sg)),

  # Create a data frame with all the variables predicted
  tar_target(ss_clean,{colData(clean_age)}),
  tar_target(ss_clean_allvariables,
             savecoldata(clean_all_deconv,
                         dir = custom_paths[["ss_clean_path"]], file = "ss_clean.csv",
                         quote = F,sep = ","),
             deployment = "main"),


  tar_target(ss_clean_dt,data.table::as.data.table(data.frame(ss_clean)),deployment = "main"),
  tar_target(plotvars,bplots_var,
             packages = "SummarizedExperiment"),


  # T.test and boxplot for deconvoluted data
  tar_target(boxplot,ttest_boxplot(data=ss_clean_allvariables,
                                   sampGroup = sampGroups,
                                   path = custom_paths[["qc_folder"]])),

  # Distribution plots:
  tar_target(distribution_pl,distribution_plots(data=ss_clean_allvariables,
                                                variable_interest = sampGroups,
                                                path = custom_paths[["qc_folder"]])),

  #Get annotation
  tar_target(ann, minfi::getAnnotation(clean_age)),
  # tar_target(betas_path, file.path(custom_paths$temp,"betas.bk"), format = "file"),

  #Get beta values and save to dis
  tar_target(betas, betas_disk(clean_age,backingfile=file.path(custom_paths$temp,"betas")),
  packages = c("minfi","bigstatsr")),
  # tar_target(betas, minfi::getBeta(clean)),
  tar_target(betas_idx, {b<-minfi::getBeta(clean_age);list(rn=rownames(b),cn=colnames(b))}),

  #Select top beta values
  tar_target(top,top_beta(betas,betas_idx,n=top_betas_N),pattern=map(top_betas_N)),


  # PCA
  tar_target(pca, pca_res(top,
                          ss,
                          sampGroups = sampGroups,
                          filename=paste0(data_names,"_pc_plot",NROW(top)),
                          path=custom_paths[["bplots_folder"]]),pattern=map(top)),
  tar_target(pca_corrplot,corpca(beta_top100 = top,                              # Correlation plot
                                 metadata=ss_clean_allvariables,
                                 path=paste0(custom_paths$corrplot_folder,"/",NROW(top),"/"),
                                 filename=paste0(data_names,"_pca_corrplot",NROW(top),".png"),
                                 title=paste0("PC1-6 correlations with ",data_names," clinical vars(top ",NROW(top)," )")
                                 ),
             pattern=map(top)),
  tar_target(bplots, bplot(pca,                                                  # Bi plots for PCA components
                           ss=ss_clean_allvariables,
                           colgroup=plotvars,
                           s=sg,
                           pal=pal_discrete,
                           alfa=0.7,
                           idcol=idcol,
                           folder = paste0(custom_paths$bplots_folder,"/",NROW(pca$rotation),"/")),
             packages = c("ggfortify","ggrepel","gplots","ggplot2"),
             pattern=map(pca)
  ),

  # Heatmaps
  tar_target(heatmap, heatmap_top(top,
                                  ss,
                                  sampGroups = sampGroups,
                                  path=custom_paths[["heatmap_folder"]],
                                  filename=paste0(data_names,"_heatmap_",NROW(top))),pattern = map(top)),


  #Create the Models
  tar_target(model_covs,covs),
  tar_target(model, mod(object = betas, betas_idx = betas_idx, group_var = group_var, contrasts = Contrasts, metadata = ss_clean_allvariables,covs=covs),
                        packages=c("limma","stats","bigstatsr")
  ),
  tar_target(dmps, DMPextr( fit= model,
                            betas_idx=betas_idx,
                            ContrastsDM = colnames(model$contrasts),
                            beta_normalized = betas,
                            p.value = 1,
                            mDiff = 0,
                            ann = ann,
                            writeOut = F),
             deployment = "worker"
  ),
  tar_target(save_dmps,
             write.table(dmps,
                         file.path(
                           custom_paths$dmp_folder,
                           paste0(as.character(quote(norm)),"_dmps.txt")
                           )
                         )
             ),

  # volcano plot
  tar_target(volcano,volcanoplot(dmps,path = custom_paths[["dmpplots_folder"]])),

  # manhattan plot
  tar_target(manhattan,manhattanplot(dmps,path = custom_paths[["dmpplots_folder"]])),

  # tar_target(dmp_battery,priority = 1,                                           # DMPs distribution along params.
  #            apply_filter_dmps(
  #              dmps = dmps,path=file.path(paste0(custom_paths$dmp_folder,data_names))),
  #            error ="continue",deployment = "worker",memory = "transient"),

  tar_target(dmps_f , filter_dmps(dmps, adj.p.value=adjp.valueDMP,p.value = p.valueDMP, mDiff = mDiffDMP)),      # Choose filter for DMPs
  tar_target(save_dmps_f,
             write.table(dmps_f,
                         file.path(
                           custom_paths$dmp_folder,
                           paste0(as.character(quote(norm)),"_dmps_filtered.txt")
                         )
             )
  ),
  tar_target(dmp_genes,{
    if(nrow(dmps_f)>0){
      dt <- data.table::setDT(dmps_f)
      dt_g <- dt[,.(Gene=unlist(lapply(strsplit(UCSC_RefGene_Name, ";"),'['))),by=c("diff_meanMeth","Contrast")]
      dt_g$meandiff<-dt_g$diff_meanMeth
      dt_g
    }else{
      dmp_genes<-data.table::data.table(meandiff=numeric(0),Contrast=character(0),Gene=character(0))
      }
  }),
  tar_target(dmp_pathways, get_pathways(dmp_genes,res.folder =paste0(custom_paths$pathway_folder,"/DMPS"),savefile=TRUE)),
  tar_target(dmps_summary,                                                       # Summary statistics for DMPs
            summary_dmps(dmps_f, dir = custom_paths$dmp_folder,name=data_names,write=TRUE),error ="continue"),
  tar_target(dmpplots, plotDMP(dmps_f,path=custom_paths$dmpplots_folder),error ="continue"),   # Barplots hipo/hyper, genomic region, CpG islands.

  # DMRS

  # tar_target(
  #   dmrs,
  #   { nrow(dmps_summary)
  #     load('dmrs_obj.RData')
  #     return(dmrs_obj)  # Return the dmrs object
  #   }
  # ),

  tar_target(dmrs,                                                              # Finds DMRs with dmrcate can be relaxed here and filter by HMFDR later
             find_dmrs(object=dmps_f,model=model,
                       fdr = fdrDMR,bcutoff = 0.05, min.cpg=min.cpgDMR),
             deployment = "worker"),
  #tar_target(dmrs_battery,  priority = 1, error ="continue",                                       # DMRs distribution along params
             # apply_filter_dmrs(
             #   dmrs = dmrs,path=paste0(custom_paths$dmrs_folder,data_names))),
  tar_target(dmrs_f, filter_dmrs(dmrs,p.value = 'FDR', mDiff = mdiffDMR, min.cpg=min.cpgDMR)),
  tar_target(save_dmrs,                                                          # Saves DMRs
             writexl::write_xlsx(
               dmrs, paste0(custom_paths$dmrs_folder,"_",data_names,".xlsx"))),
  tar_target(sumaries,summarize(dmrs = dmrs_f, dmps = dmps_f,path = paste0(custom_paths$results_folder,"/",data_names,"/"))),
  tar_target(dmr_genes,{
    if(nrow(dmrs_f)>0){
      dt <- data.table::setDT(dmrs_f)
      dt_g <- dt[,.(Gene=unlist(lapply(strsplit(overlapping.genes, ";"),'['))),by=c("meandiff","Contrast")]

      dt_g
    }else{
      dt_g<-data.table::data.table(meandiff=numeric(0),Contrast=character(0),Gene=character(0))
    }
  }),

  tar_target(dmr_pathways, get_pathways(dmr_genes,res.folder =paste0(custom_paths$pathway_folder,"/DMRs/"),savefile=TRUE)),
  tar_target(dmrs_summary,                                                       # Summary stats for DMRs
              summary_dmrs(
                dmrs,path=file.path(custom_paths$dmrs_folder,paste0("full_dmrs_summary",data_names,".csv"))),
              error = "continue"),













  #####################################################################
  tar_target(valss,{
    nrow(dmrs_summary)
    valueslist = makelist(vals)
    return(valueslist)
  }),
  tar_target(ep,tibble::tibble(report_parameters,
                               values_row=jsonlite::toJSON(valss),
                               paths=jsonlite::toJSON(custom_paths),
                               output_file = file.path(custom_paths$report_folder,paste0(report_parameters$report_name,".html"))

  # ))
  )),
  tar_quarto_rep(report,
                 path = "_qmd/report.qmd",
                 execute_params = ep,
                 priority = 0
  )
)
# quarto::quarto_render(input = "report.qmd",execute_params = tar_read("report_params_ex_EPIC")[1,])

combined <- tarchetypes::tar_combine(
  combined_summary,
  targets[["sumaries"]],
  command = dplyr::bind_rows(!!!.x,.id = "experiment")
)

combined_ss <- tarchetypes::tar_combine(
  ss_comb,
  targets[["ss_clean_dt"]],
  command = dplyr::bind_rows(!!!.x,.id = "experiment")
)

save_combined<-tar_target(
  save_combined_summary,{
    data.table::setorder(combined_summary,experiment,Contrast)
    data.table::fwrite(combined_summary,paste0(results_folder,"/full_summary.csv"))
  }
)


list(
  top_betas_N_target,
  targets,combined,save_combined,combined_ss)
