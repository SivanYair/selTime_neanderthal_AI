# Take feather data frame objects written by python,
# reformat them if they need to be vectors, and then
# convert them to RDS objects for long term storage and access in R

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("Please supply version letter, sample size option, scenario name, and run", call.=FALSE)
}

versionLetter=args[1]
sampleSize_option=args[2]
scenario=args[3]
run=args[4]

library(feather)

# input path
feather_file_path=paste0("selection_output/version",versionLetter,"/ftr_files/")
# output path
R_file_path=gsub("ftr_","R",feather_file_path)
# specifier for file name
file_prefix=paste0('ss',sampleSize_option,'_',scenario,'_run',run)

# frequency data frame
freqs = data.frame(read_feather(paste0(feather_file_path,file_prefix,"_freqs.ftr")))
saveRDS(freqs,paste0(R_file_path,file_prefix,"_freqs.RDS"))

# positions
pos = read_feather(paste0(feather_file_path,file_prefix,"_positions.ftr"))$pos
saveRDS(pos,paste0(R_file_path,file_prefix,"_positions.RDS"))

# sample sizes
samplesizes = unlist(read_feather(paste0(feather_file_path,file_prefix,"_sampleSizes.ftr")))
saveRDS(samplesizes,paste0(R_file_path,file_prefix,"_sampleSizes.RDS"))

# sel site index
ss_index = read_feather(paste0(feather_file_path,file_prefix,"_selSiteIndex.ftr"))$ss_index
saveRDS(ss_index,paste0(R_file_path,file_prefix,"_selSiteIndex.RDS"))

# sel site sampled freqs
sampledFreqs = unlist(data.frame(read_feather(paste0(feather_file_path,file_prefix,"_sampled_finFreqs.ftr"))))
saveRDS(sampledFreqs,paste0(R_file_path,file_prefix,"_sampled_finFreqs.RDS"))

# sel site true freqs
trueFreqs = unlist(data.frame(read_feather(paste0(feather_file_path,file_prefix,"_true_finFreqs.ftr"))))
saveRDS(trueFreqs,paste0(R_file_path,file_prefix,"_true_finFreqs.RDS"))
