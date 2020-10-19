# working directory: selTime_neanderthal_AI/method

# before running this script, need to run calcNeutralF.R


## arguments to supply on command line:
## (1) full path to file that describes a scenario that we simulated
## (2) sample size option: a number corresponding to the sample sizes used to obtain allele frequencies
##     options are defined in selTime_neanderthal_AI/specification_files/sampleSizes_allOptions.R
## (3) run number of the simulation

# load libraries
library(stringr)
library(MASS)
library(parallel)
library(gtools) # gives permutations() function
library(purrr)


args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("Please supply args file (full path), sample size option, and run number to analyze", call.=FALSE)
}
cat(args, sep = "\n")

sim_args_file=args[1]
sampleSize_option=args[2]
run=args[3]
print(paste("run:",run))

# obtain info about the simulations that produced the data we analyze
# "VersionLetter" and "scenario" are the only objects that we use later on, but we list
# the other contents of the file here as a reminder of how we describe a scenario
sim_args=readLines(sim_args_file)
versionLetter=sim_args[1] # letter describing the version of demography being used
true_t_btwn=as.numeric(sim_args[2])
true_sel_coef=as.numeric(sim_args[3])
true_fin_freq=as.numeric(sim_args[4])
true_selPops=strsplit(sim_args[5],split=" ")[[1]]
scenario=tail(strsplit(gsub(".txt","",sim_args_file),split="/")[[1]],n=1) # name of the scenario that we simulated

### path to data files ###
# path to data for selection simulations
input_path_target=paste0('../selection_simulations/selection_output/version',versionLetter,'/Rfiles/')
# prefix of filenames for data produced by selection simulations
input_prefix=paste0('ss',sampleSize_option,'_',scenario,'_run',run)
# prefix of filenames for the output of the method
general_prefix=paste0('ss',sampleSize_option,'_',scenario,'_run',run)
output_path=paste0('method_output/version',versionLetter,'/')

# file to save composite likelihood results
CL_file=paste0(output_path,'CL_files/',general_prefix,'.RDS')

# Load info on demography associated with simulations
# Loads pops (a vector of full/non-partitioned population names, in correct order),
#       div_times_from_present, nested_pops, nest_index, migs, t_sampled
load(paste0('~/slim_sims/specification_files/version',versionLetter,'_demography.RData'))

# parameters necessary to generate F matrices under selection
g=as.numeric(migs[migs$donor=="Neanderthal","rate"])
t_int=as.integer(migs[migs$donor=="Neanderthal","end_from_present"])
rec_rate=1e-8
npos=2e6 # number of sites simulated in bp; change if simulating larger locus


#### LOAD DATA PRODUCED AHEAD OF SCRIPT (PROCESSED SIMULATION OUTPUT)
positions=readRDS(paste0(input_path_target,input_prefix,'_positions.RDS') )
freqs_notRandomized=readRDS(paste0(input_path_target,input_prefix,'_freqs.RDS') )
sampleSizes=readRDS(paste0(input_path_target,input_prefix,'_sampleSizes.RDS') )
sel_site_index=readRDS(paste0(input_path_target,input_prefix,'_selSiteIndex.RDS') )
fin_freqs=readRDS(paste0(input_path_target,input_prefix,'_sampled_finFreqs.RDS') )

# name the posiiton corresponding to the selected site / partition site
# need to add 1 to sel_site_index (originally saved in python) to make indexing consistent with R, not python
names(positions)[sel_site_index+1]="sel_site"


#### NAMING OF POPULATIONS, USED IN VARIOUS SCRIPTS
whole_pops=pops # reassign pops to whole_pops, because we will start using partition names
ancient=c('EF','Steppe','EurUP','WHG') # ancient, whole/full population names
pops=colnames(freqs_notRandomized) # names of partitioned populations
## official populations used in analysis
# remove AA and Aa partition pops of a full population if it has Neanderthal allele at less than 5% frequency
investigate=pops[grepl('_AA',pops) | grepl('_Aa',pops) ]
for(pop in investigate){
  if(fin_freqs[gsub('_..$','',pop)]<0.05){
    pops=setdiff(pops,pop)
  }
}
freqs_notRandomized=freqs_notRandomized[,pops]
sampleSizes=sampleSizes[pops]

pops_ancient=pops[sapply(pops,function(p) any(sapply(ancient,function(a) grepl(a,p) ))  )] # names of partitioned ancient populations


# Include _AA or _aa populations if their information is not listed and
#         their corresponding _Aa population exists
#         Note: These additional populations will be excluded when calculating
#         composite likelihoods, since they have sample size = 0 at all sites,
#         but we need them when generating selection matrices, because _Aa
#         population's elements are a linear combination of _AA and _aa.
add_pops=c()
for(pop in ancient) { # only ancient pops can have _Aa partitions
  if(paste0(pop,'_Aa') %in% pops_ancient) {
    if( !(paste0(pop,'_aa') %in% pops_ancient)  ) {
      add_pops=c(add_pops,paste0(pop,'_aa'))
    }

    if( !(paste0(pop,'_AA') %in% pops_ancient) ) {
      add_pops=c(add_pops,paste0(pop,'_AA'))
    }
  }
}
## freq = NA and sample size = 0 (mainly used so that we can calculate Aa expectation, and for inference)
pops_ancient = c(pops_ancient,add_pops)
pops=c(pops,add_pops)
numPops=length(pops)

# PART POPS THAT EXIST CHANGE BASED ON THE SCENARIO WE ARE CONSIDERING

# name the "source" population of introgression and the population that can
# never be selected ("neverSel") in our models because it does not contain introgressed ancestry
source="Neanderthal"
neverSel="Africa"

presentDay_pops=c('Africa',grep('East_Asia',pops,value=T),grep('Europe',pops,value=T))

# putative selected populations
potential_selPops=sel_pops=pops[grep('_AA',pops)]


# randomize allele frequencies: get dim numPops x numLoci
allFreqs = apply(freqs_notRandomized, 1, function(my.freqs) {
  if(runif(1)<0.5){my.freqs<-1-my.freqs} ; my.freqs
})


# below script filters sites in the data set, loads neutral F matrix,
# and creates objects for relevant admixture graph parameters
source('prep_data_perRun.R')
print("finished data prep")

# Define functions that generate F matrix predictions under the null and selection models (first script) and that calculate composite likelihoods (second script)
source('genSelMatrices_fxns.R')
source('calcCompLike_fxns.R')


# Run script that generates selection matrices under different proposed
# parameter combinations of s and t_b
source('genSelMatrices_exec.R')

# Run script that calculates & saves composite likelihoods of s & t_b based on their predicted F^(S) or F^(N)
source('calcCompLike_exec.R')

# Run script that identifies estimates for t_b and s and plots their profile composite likelihood surfaces
source('plotCompLikeResults.R')
