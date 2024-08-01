

# Find dmrs--->Rscript:find_dmrs2
dmrs_obj<-find_dmrs(object=dmps_f,model=model,
          fdr = fdrDMR,bcutoff = 0.05, min.cpg=min.cpgDMR)


#dmrs_f2<-filter_dmrs(dmrs,p.value = 'FDR', mDiff = mdiffDMR, min.cpg=min.cpgDMR)

save(dmrs_obj,file='dmrs_obj.RData')
