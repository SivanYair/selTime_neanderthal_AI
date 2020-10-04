# Sample size files

# This is for number of individuals sampled, not number of genomes sampled
# population names need to exactly match those that were previously specified, but the order doesn't matter
# Sample sizes are saved for python with feather; it will be turned into a dictionary there
# Sample sizes*2 (for genomes) are saved as an RDS for use in R later

library(feather)

# OPTION 1
option=1
sampleSizes_indiv=c(Africa=100,Neanderthal=1,East_Asia=100,EurUP=6,Europe=100,WHG=50,EF=145,Steppe=20)
sampleSizes_indiv_python=data.frame(pop=names(sampleSizes_indiv),size_indiv=sampleSizes_indiv,row.names=NULL)
write_feather(sampleSizes_indiv_python, paste0("~/slim_sims/specification_files/option",option,"_fullPop_indivSampleSizes.feather"))
sampleSizes_genomes_R=sampleSizes_indiv*2
saveRDS(sampleSizes_genomes_R,paste0("~/slim_sims/specification_files/option",option,"_fullPop_genomeSampleSizes.RDS"))

# OPTION 2
