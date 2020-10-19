# This script is sourced in prep_data_perRun.R
# It defines the following objects:
## p: matrix of admixture proportions (before being changed for different selection scenarios)
## mig_times: vector describing migration times (in generations) from ancient populations into present day Europeans


# p = matrix defining proportion of the population that comes from another
# Row: recipient pop, Column: donor pop
numPops=length(pops)
p = matrix(data=rep(0,numPops^2),nrow=numPops, ncol=numPops, dimnames=list(pops,pops))
diag(p)=1 #start by assuming no migration
# add in cases for migration, if the Europe_AA population exists:
migration = paste0(Eur,'_AA') %in% pops & any(c('EF_AA','Steppe_AA','WHG_AA') %in% pops) # variable telling us whether we should account for migration
if(migration){
  # MIGRATION RATES
  if(any(grepl('EF_AA',pops))){
    p[paste0(Eur,'_AA'),'EF_AA']=0.52 #0.5176239857
  }
  if(any(grepl('Steppe_AA',pops))){
    p[paste0(Eur,'_AA'),'Steppe_AA']=0.36 #0.3585322047
  }
  if(any(grepl('WHG_AA',pops))){
    p[paste0(Eur,'_AA'),'WHG_AA']=0.12
  }

  # MIGRATION TIMES
  # consistent with simulations (version A demography)
  mig_times = matrix(data=rep(NA,numPops^2),nrow=numPops, ncol=numPops, dimnames=list(pops,pops))
  if('EF_AA' %in% pops){
    mig_times[paste0(Eur,'_AA'),'EF_AA']=as.integer(migs$start_from_present[migs$donor=='EF'])
  }
  if('Steppe_AA' %in% pops){
    mig_times[paste0(Eur,'_AA'),'Steppe_AA']=as.integer(migs$start_from_present[migs$donor=='Steppe']) 
  }
  if('WHG_AA' %in% pops){
    mig_times[paste0(Eur,'_AA'),'WHG_AA']=as.integer(migs$start_from_present[migs$donor=='WHG']) 
  }

}

# Following adjusts migration to self if migration rates were modified
for (i in 1:nrow(p)) {
  p[i,i]=1-sum(p[i,-i])
}


rm(migs) # removing original migs matrix from version A demography to avoid accidentally using it later
