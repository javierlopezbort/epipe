# Install tidyr if you haven't already
# install.packages("tidyr")

library(dplyr)
library(tidyr)

# Separate the Basename column into sentrix_ID and Sentrix_position
ss_clean_2 <- ss_clean_new %>%
  separate(Basename, into = c("Sentrix_ID", "Sentrix_position"), sep = "_")

# View the first few rows of the updated tibble
head(ss_clean_2)


ss_clean_2 <- ss_clean_2 %>%
  rename(Basename = Basename_1)


ss_clean_2 <- subset(ss_clean_2, select = -c(ANA_dom1, ANA_dom2, ANA_dom3, ANA_wt1, ANA_wt2,ANA_wt3,ANA_ven,ANA_lymph,ANA_ven_lymph,ANA_mut,ANA_dom,ANA_WT,ANA_SANDRA))

# Keep specific rows
ss_clean_2 <- ss_clean_2 %>%
  slice(c(1, 2, 3, 10,16,21,25,26,27))

ss_clean_2$Sentrix_ID<-as.factor(ss_clean_2$Sentrix_ID)
ss_clean_2$mutation<-as.factor(ss_clean_2$mutation)
ss_clean_2$condition<-as.factor(ss_clean_2$condition)
ss_clean_2$age<-as.factor(ss_clean_2$age)

saveRDS(ss_clean_2,file='samplesheet.rds')
