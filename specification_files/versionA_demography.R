options(stringsAsFactors = F)

### below list defines order of population splits <- each index corresponds to a different divergence time
# NOTE: the nestedness must correspond to order of divergence
#       populations split off from the previous population that split, so this
#       must be considered when ordering them
#       only last index may contain multiple populations (otherwise need to have a set of pointers from one pop to another)
# first population MUST be Neanderthals and second population MUST be ancestor of MH
nested_pops = list('Neanderthal','Africa','East_Asia','EurUP','Europe',list('Steppe','EF','WHG'))
# note: Eur has to split from EA separately from ancient pops to represent some shared time btwn Eur + ancient pops
### assign nest index for each partitioned pop (which nest does it belong to?)
pops=unlist(nested_pops) # index of pops here corresponds to index in slim script
nest_index = sapply(pops, function(pop) which(sapply(nested_pops, function(nest) pop %in% nest )) )
### following vector has same length as nested_pops, and describes the time that
### the populations in that nest split from the populations listed after that index
### times described are time between present and divergence
### each time corresponds to when given pop that matches index AND the rest split off from the previous pop index
div_times_from_present = c(NA,16000,2500,1724,1552,1379)

### migration information: recipient, donor, per-generation rate, start gen, end gen
### note that the rate should account for the fact that contribution will be diluted in later generations

### based on migration rates, WHG needs to be first migrating pop,  then EF, then Steppe
migs=as.data.frame(do.call('rbind',list(
  c(recipient='East_Asia',donor='Neanderthal',rate='0.02',start_from_present=2071, end_from_present=2070),# first pop that represents Eurasians is the one that gets introgressed
  c(recipient='Europe',donor='WHG',rate='0.999',start_from_present=227, end_from_present=226),
  c(recipient='Europe',donor='EF',rate='0.8125',start_from_present=225, end_from_present=224),
  c(recipient='Europe',donor='Steppe',rate='0.36',start_from_present=156, end_from_present=155)
)))

t_sampled_from_present=c(Neanderthal=1550,EurUP=1164,WHG=302,EF=246,Steppe=167)

# save data to load in slim script-making R scripts
save.image('~/slim_sims/specification_files/versionA_demography.RData')

# need to save pops to python readable object? perhaps just run this version once, then load the data..
library(feather)
write_feather(as.data.frame(pops), "~/slim_sims/specification_files/versionA_pops.feather")
