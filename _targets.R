####### Epigenetics PIPEline Main Workflow #######

# # Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(data.table)
library(quarto)
library(SummarizedExperiment)
suppressPackageStartupMessages(library(qs2))
library(crew)


# library(crew.cluster)
#
# # Libraries from Illumina
# library(IlluminaHumanMethylationEPICv2manifest)
# library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
#
# max_ncores=RcppParallel::defaultNumThreads() # Not in use but may be useful

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

# Source scripts files
#tar_source(system.file(package = "epipe"))


# Source config.R file
#source(system.file("config.R", package = "epipe"))
source("config.R")


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
  tar_target(custom_paths,epipe::make_results_dirs(subf=data_names, results_folder = results_folder, analysis_folder = analysis_folder)),
  tar_target(samplesheet_path, data_paths, format = "file"),
  tar_target(samplesheet, readRDS(samplesheet_path)),
  tar_target(ss, epipe::get_category(samplesheet)),

  # Read idats:
  tar_target(rgSet, epipe::read_idats(ss,arraytype = arraytype,idats_folder=idats_folder,idcol=idcol, author= author, description = description)), # makes rgChannelset object using minfi and annotates according to arraytype

  # Qc report:
  tar_target(QC_plotss, epipe::qc(rgSet,
                          sampGroups = sampGroups, #names(ss)[attributes(ss)$category == "mgroups"][1],
                          idcol=idcol,
                          pal=pal_discrete,
                          qc_folder = custom_paths[["qc_folder"]]
                          )
             ,packages = "minfi",deployment = "worker" ),

  # Filters:  -probes: pval<0.01, -samples: 10% probes fail, Plots: Sample p.values barplot (colMeans)
  tar_target(filtered, epipe::filter(
    rgSet=rgSet,
    sampGroups= sampGroups, #names(ss)[attributes(ss)$category == "mgroups"][1],
    sampNames = idcol,
    qc_folder = custom_paths[["qc_folder"]],
    save_barplot=T)),

  #Tumor purity
  #tar_target(purity, purify(rgSet = rgSet, arraytype = arraytype),error = "continue"),

  #Get raw beta values (before normalization)
  tar_target(raw_beta, minfi::getBeta(filtered)),
  tar_target(save_raw_beta,write.csv(raw_beta, file=file.path(custom_paths$betas_folder,'raw_beta.csv'))),

  #Normalization
  tar_target(normalize, epipe::norm(filtered)),
  tar_target(dplot_normalize,epipe::denplot(normalize,
                                ss,
                                sampGroups = sampGroups, #names(ss)[attributes(ss)$category == "mgroups"][1],
                                path=custom_paths[["qc_folder"]],
                                norm_method = norm_function)),
  tar_target(normalize_all, epipe::normalization_all_methods(filtered,
                                                      ss,
                                                      sampGroups = sampGroups, #names(ss)[attributes(ss)$category == "mgroups"][1],
                                                      path=custom_paths[["qc_folder"]],
                                                      norm_method = norm_function)),

  # Get beta values after normalization
  tar_target(norm_beta, minfi::getBeta(normalize)),
  tar_target(save_norm_beta,write.csv(norm_beta, file=file.path(custom_paths$betas_folder,'norm_beta.csv'))),


  #Preprocessing (probe removal)
  tar_target(clean, epipe::prep(normalize, remove_sex = remove_sex, arraytype = arraytype,sexplot_folder= custom_paths[["sexplot_folder"]],predict_sex=sex_prediction)),


  #Predict age
  tar_target(clean_age,epipe::ageprediction(clean),error='continue'),

  #Deconvolution
  tar_target(clean_all_deconv,epipe::celldeconvolution(rgSet,clean_age,arraytype = arraytype),error = 'continue'),

  # Correlation analysis
  tar_target(cor_analysis, epipe::correlation_analysis(clean_all_deconv,
                                                path=custom_paths[["qc_folder"]],variables=variables,sampGroups=sampGroups)),

  # Create a data frame with all the variables predicted
  tar_target(ss_clean,{colData(clean_age)}),
  tar_target(ss_clean_allvariables,
             epipe::savecoldata(clean_all_deconv,
                         dir = custom_paths[["ss_clean_path"]], file = "ss_clean.csv",
                         quote = F,sep = ","),
             deployment = "main"),


  tar_target(ss_clean_dt,data.table::as.data.table(data.frame(ss_clean)),deployment = "main"),
  tar_target(plotvars,bplots_var,
             packages = "SummarizedExperiment"),


  # T.test and boxplot for deconvoluted data
  tar_target(boxplot,epipe::ttest_boxplot(data=ss_clean_allvariables,
                                   sampGroup = sampGroups,pal=pal_discrete,
                                   path = custom_paths[["qc_folder"]])),

  # Distribution plots:
  tar_target(distribution_pl,epipe::distribution_plots(data=ss_clean_allvariables,
                                                variable_interest = sampGroups,
                                                path = custom_paths[["qc_folder"]])),

  #Get annotation
  tar_target(ann, minfi::getAnnotation(clean_age)),
  # tar_target(betas_path, file.path(custom_paths$temp,"betas.bk"), format = "file"),

  #Get beta values and save to dis
  tar_target(betas, epipe::betas_disk(clean_age,backingfile=file.path(custom_paths$temp,"betas")),
  packages = c("minfi","bigstatsr")),
  # tar_target(betas, minfi::getBeta(clean)),
  tar_target(betas_idx, {b<-minfi::getBeta(clean_age);list(rn=rownames(b),cn=colnames(b))}),

  #Select top beta values
  tar_target(top,epipe::top_beta(betas,betas_idx,n=top_betas_N),pattern=map(top_betas_N)),


  # PCA
  tar_target(pca, epipe::pca_res(top,
                          ss,
                          sampGroups = sampGroups,
                          filename=paste0(data_names,"_pc_plot",NROW(top)),
                          path=custom_paths[["bplots_folder"]]),pattern=map(top)),
  tar_target(pca_corrplot,epipe::corpca(beta_top100 = top,                              # Correlation plot
                                 metadata=ss_clean_allvariables,
                                 path=paste0(custom_paths$corrplot_folder,"/",NROW(top),"/"),
                                 filename=paste0(data_names,"_pca_corrplot",NROW(top),".png"),
                                 title=paste0("PC1-6 correlations with ",data_names," clinical vars(top ",NROW(top)," )")
                                 ),
             pattern=map(top)),
  tar_target(bplots, epipe::bplot(pca,                                                  # Bi plots for PCA components
                           ss=ss_clean_allvariables,
                           colgroup=plotvars,
                           s=sampGroups,
                           pal=pal_discrete,
                           alfa=0.7,
                           idcol=idcol,
                           folder = paste0(custom_paths$bplots_folder,"/",NROW(pca$rotation),"/")),
             packages = c("ggfortify","ggrepel","gplots","ggplot2"),
             pattern=map(pca)
  ),

  # Heatmaps
  tar_target(heatmap, epipe::heatmap_top(top,
                                  ss,
                                  sampGroups = sampGroups,
                                  path=custom_paths[["heatmap_folder"]],
                                  filename=paste0(data_names,"_heatmap_",NROW(top))),pattern = map(top)),


  #Create the Models
  tar_target(model_covs,covs),
  tar_target(model, epipe::mod(object = betas, betas_idx = betas_idx, group_var = group_var, contrasts = Contrasts, metadata = ss_clean_allvariables,covs=covs),
                        packages=c("limma","stats","bigstatsr")
  ),
  tar_target(dmps, epipe::DMPextr( fit= model,
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
  tar_target(volcano,epipe::volcanoplot(dmps,path = custom_paths[["dmpplots_folder"]])),

  # manhattan plot
  tar_target(manhattan,epipe::manhattanplot(dmps,path = custom_paths[["dmpplots_folder"]])),

  # tar_target(dmp_battery,priority = 1,                                           # DMPs distribution along params.
  #            apply_filter_dmps(
  #              dmps = dmps,path=file.path(paste0(custom_paths$dmp_folder,data_names))),
  #            error ="continue",deployment = "worker",memory = "transient"),

  tar_target(dmps_f , epipe::filter_dmps(dmps, adj.p.value=adjp.valueDMP,p.value = p.valueDMP, mDiff = mDiffDMP)),      # Choose filter for DMPs
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
  tar_target(dmp_pathways, epipe::get_pathways(dmp_genes,res.folder =paste0(custom_paths$pathway_folder,"/DMPS"),savefile=TRUE)),
  tar_target(dmps_summary,                                                       # Summary statistics for DMPs
            epipe::summary_dmps(dmps_f, dir = custom_paths$dmp_folder,name=data_names,write=TRUE),error ="continue"),
  tar_target(dmpplots, epipe::plotDMP(dmps_f,path=custom_paths$dmpplots_folder),error ="continue"),   # Barplots hipo/hyper, genomic region, CpG islands.

  # DMRS

  # tar_target(
  #   dmrs,
  #   { nrow(dmps_summary)
  #     load('dmrs_obj.RData')
  #     return(dmrs_obj)  # Return the dmrs object
  #   }
  # ),

  tar_target(dmrs,                                                              # Finds DMRs with dmrcate can be relaxed here and filter by HMFDR later
             epipe::find_dmrs(object=dmps_f,model=model,
                       fdr = fdrDMR,bcutoff = 0.05, min.cpg=min.cpgDMR),deployment = "worker"),




  #tar_target(dmrs_battery,  priority = 1, error ="continue",                                       # DMRs distribution along params
             # apply_filter_dmrs(
             #   dmrs = dmrs,path=paste0(custom_paths$dmrs_folder,data_names))),
  tar_target(dmrs_f, epipe::filter_dmrs(dmrs,p.value = 'FDR', mDiff = mdiffDMR, min.cpg=min.cpgDMR)),
  tar_target(save_dmrs,                                                          # Saves DMRs
             writexl::write_xlsx(
               dmrs, paste0(custom_paths$dmrs_folder,"_",data_names,".xlsx"))),
  tar_target(sumaries,epipe::summarize(dmrs = dmrs_f, dmps = dmps_f,path = paste0(custom_paths$results_folder,"/",data_names,"/"))),
  tar_target(dmr_genes,{
    if(nrow(dmrs_f)>0){
      dt <- data.table::setDT(dmrs_f)
      dt_g <- dt[,.(Gene=unlist(lapply(strsplit(overlapping.genes, ";"),'['))),by=c("meandiff","Contrast")]

      dt_g
    }else{
      dt_g<-data.table::data.table(meandiff=numeric(0),Contrast=character(0),Gene=character(0))
    }
  }),

  tar_target(dmr_pathways, epipe::get_pathways(dmr_genes,res.folder =paste0(custom_paths$pathway_folder,"/DMRs/"),savefile=TRUE)),
  tar_target(dmrs_summary,                                                       # Summary stats for DMRs
              epipe::summary_dmrs(
                dmrs,path=file.path(custom_paths$dmrs_folder,paste0("full_dmrs_summary",data_names,".csv"))),
              error = "continue"),













  #####################################################################
  tar_target(valss,{
    nrow(dmrs_summary)
    valueslist = epipe::makelist(vals)
    return(valueslist)
  }),
  tar_target(ep,tibble::tibble(report_parameters,
                               values_row=jsonlite::toJSON(valss),
                               paths=jsonlite::toJSON(custom_paths),
                               output_file = file.path(custom_paths$report_folder,paste0(report_parameters$report_name,".html"))

  # ))
  )),
  tar_quarto_rep(report,
                 path = system.file("_qmd/report.qmd", package = "epipe"),
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
