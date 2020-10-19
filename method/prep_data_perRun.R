
### SITE FILTERING IN DATA SET

# remove monomorphic sites
polymorphic = apply(allFreqs,2,function(locus_freqs) !(all(locus_freqs==0) | all(locus_freqs==1)) )
allFreqs=allFreqs[,polymorphic]
positions=positions[polymorphic]

# remove cases in which all present day pops have freq = 0
good = apply(allFreqs[presentDay_pops,],2,function(site_freqs) (any(site_freqs>0) & any(site_freqs<1)) | (any(site_freqs==0) & any(site_freqs==1)) )
allFreqs=allFreqs[,good]
positions=positions[good]

# remove sites with private alleles at less than 1% frequency
# need to get frequency in whole pop to evaluate this
# first find out which partPops belong to whole pops
partPops=lapply(whole_pops,function(pop) grep(pop,real_pops,value=T))
names(partPops)=whole_pops
# then merge partPops frequencies into whole pops
allFreqs_whole = apply(allFreqs,2,function(site_freqs){
  sapply(whole_pops,function(pop){
    sum(site_freqs[partPops[[pop]]]*sampleSizes[partPops[[pop]]]) / sum(sampleSizes[partPops[[pop]]]) # get weighted average of frequency i.e. whole population frequencies
  })
})
#function to identify sites to keep If it returns F, exclude the site
keep_noLowFreqPrivate=apply(allFreqs_whole,2,function(row){
  if(sum(row>0,na.rm=T)>1 & sum(row<1,na.rm=T)>1){
    return(T)
  } else{
    filt_row = row[!is.na(row)]
    private_pop=which(filt_row>0 & filt_row<1)
    if(length(private_pop)==0){return(T)}
    if((filt_row[private_pop]<0.01 & sum(filt_row[-private_pop])==0) | (filt_row[private_pop]> 0.99 & sum(filt_row[-private_pop])==1)){
      return(F)
    } else{
      return(T)
    }
  }
} )
allFreqs=allFreqs[,keep_noLowFreqPrivate]
positions=positions[keep_noLowFreqPrivate]

# now, subset number of segregating sites to be more realistic:
# want it to match what we observe in reality:  ~12k sites
sel_site_pos_index=which(positions==positions["sel_site"]) # track position of sel site to make sure we include it
choose_index=sample(1:length(positions),size=12e3)
if(!(sel_site_pos_index %in% choose_index)) choose_index=c(choose_index,sel_site_pos_index) # add position of sel site in case it was excluded
choose_index=choose_index[order(choose_index)] # put indices of positions vector back in order
allFreqs=allFreqs[,choose_index]
positions=positions[choose_index]
rm(sel_site_pos_index) # remove this object, it no longer points to the index in the positions vector that corresponds to the selected site

### DEMOGRAPHY
EA="East_Asia"
Eur="Europe"
## define nestedness of populations and their divergence times
## this will be used to get t_div, t_anc, and possible selection scenarios (later, with t_b)
### populations with Neanderthal introgression
#### first, identify all partitioned populations corresponding to a certain 'whole' population
EA_pops=grep(EA,pops,value=T) # East Asian
Eur_pops=grep(Eur,pops,value=T) # European
Steppe_pops=grep("Steppe",pops,value=T) # Steppe
EF_pops=grep("EF",pops,value=T) # Early Farmer
WHG_pops=grep("WHG",pops,value=T) # Western Hunter Gatherer
EurUP_pops=grep("EurUP",pops,value=T) # West Eurasian Upper Paleolithic
#### then list these as the introgressed populations:
intro_pops = c(EurUP_pops,EA_pops,Eur_pops,Steppe_pops,EF_pops,WHG_pops)
### below list defines order of population splits <- each index corresponds to a different divergence time
## reduce nested_pops & corresponding div_times_from_present saved from version A demography to only include introgressed populations
## by doing -c(1,2); merge the last set of populations into one
nested_pops = lapply(nested_pops[-c(1,2)],function(pop_set) unlist(lapply(pop_set,function(pop) grep(pop,pops,value=T) ))  )
# need to merge the last two entries of nested pops, because they actually represent the same split from each other
# before, divergence represented when a given population AND the remaining diverge from the previous population
# now it represents when a given population diverges from the previous populations
nested_pops[[length(nested_pops)-1]]=c(nested_pops[[length(nested_pops)-1]],nested_pops[[length(nested_pops)]])
nested_pops=nested_pops[1:(length(nested_pops)-1)]
### assign nest index for each partitioned pop (which nest does it belong to?)
nest_index = sapply(intro_pops, function(pop) which(sapply(nested_pops, function(nest) pop %in% nest )) )
### following vector has same length as nested_pops, and describes the time that
### the populations in that nest split from the populations listed after that index
### times described are time between present and divergence
div_times = div_times_from_present[-c(1,2,3)]
### conversely, anc_times describes the amount of shared history a population in one index
### shares with the populations listed in the subsequent indices
### Note: anc_times is used in calcTotAddF functions, to figure out which pops share ancestral sweep
anc_times = t_int - div_times

### for ease later, create list that defines order of population splits
### specifically for potentially selected populations (keep _AA)
nested_selPops=lapply(nested_pops,function(nest) grep('_AA',nest,value=T) )
selPop_nest_index=sapply(potential_selPops, function(pop) which(sapply(nested_selPops, function(nest) pop %in% nest )) )


## DIVERGENCE TIME:
# t_div = number of generations from present until the pair of pops diverged
source('make_t_div.R')

# t_anc = number of generation from when 2 pops join until Neanderthal introgression
t_anc = t_int - t_div

## SAMPLING TIME:
# t_sampled = vector describing how long ago (from the present looking pastwards) individuals were sampled
# such that the time between when the pop was sampled to the time of Neanderthal introgression is t_int - t_sampled
# whole_pop_ID has already been defined: vector of length pops that specifies full pop that it came from
t_sampled=rep(0,length(pops)) #time that Neander + Africa were sampled doesn't matter, we'll never use it
names(t_sampled)=pops
t_sampled[whole_pop_ID %in% names(t_sampled_from_present)] = t_sampled_from_present[whole_pop_ID[whole_pop_ID %in% names(t_sampled_from_present)]]


### NEUTRAL F ESTIMATE: modify to include partitioned populations
F_estimate = readRDS("NeutralF_versionA_ss1.RDS")
# modify F_estimate by adding estimates for the subpops -- make them the same as the whole pop they came from
whole_pop_ID = gsub("_..\\b","",pops) # whole pops that each subpop (in order of object pops) belong to (\\b is for boundary, . is for any character)
for(i in 1:numPops) {
  F_estimate=cbind(F_estimate,F_estimate[,whole_pop_ID[i]])
  colnames(F_estimate)=c(colnames(F_estimate)[-ncol(F_estimate)],pops[i])
  F_estimate=rbind(F_estimate,F_estimate[whole_pop_ID[i],])
  rownames(F_estimate)=colnames(F_estimate)
}
# now remove the whole pop that the subpops were copied from
F_estimate=F_estimate[pops,pops]

### MIGRATION MATRIX: currently assumes all designated sel pops don't change across models
## source script that defines:
## p: matrix of admixture proportions (before being changed for different selection scenarios)
## mig_times: vector describing migration times (in generations) from Europe into ancient populations
source('make_admixture_objects.R')

## effective population size
## Note: we modify it for Neanderthals, but models only need this value for Neanderthal-introgressed populations
Ne = rep(10000,numPops)
names(Ne) = pops
Ne["Neander"] = 3000
