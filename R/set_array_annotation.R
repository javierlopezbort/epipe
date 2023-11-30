# Helper function to set array-specific annotation
set_array_annotation <- function(mSetSqn, arraytype) {
  if (arraytype == "IlluminaHumanMethylation450k") {
    mSetSqn@annotation <- c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19")
  } else if (arraytype == "IlluminaHumanMethylationEPIC") {
    mSetSqn@annotation <- c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b2.hg19")
  } else if (arraytype == "IlluminaHumanMethylationEPICv2") {
    mSetSqn@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
  } else if (arraytype == "IlluminaHumanMethylation27k") {
    mSetSqn@annotation <- c(array = "IlluminaHumanMethylation27k", annotation = "ilmn12.hg19")
  } else if (arraytype == "HorvathMammalMethylChip40") {
    mSetSqn@annotation <- c(array = "HorvathMammalMethylChip40", annotation = "test.unknown")
  } else {
    mSetSqn@annotation <- c(array = "chip.unknown", annotation = "test.unknown")
  }
  return(mSetSqn)
}
