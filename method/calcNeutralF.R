# working directory: selTime_neanderthal_AI/method

args=commandArgs(trailingOnly = T)
if(length(args)!=3) stop("Error: calcNeutralF.R needs 3 arguments: (1) demography version (letter) ; (2) sample size option (number) ; (3) number of neutral sim runs")
demography_version=args[1]
sampleSize_option=args[2]
numRuns=args[3]


# input and output path
path="" # add file path to neutral allele frequency matrices here

# get alt allele freq matrix: dim numLoci x numPops
allFreqs=do.call("rbind",lapply(1:numRuns, function(run) readRDS(paste0(path,"ss",sampleSize_option,"_run",run,"_freqs.RDS")) ))
sampleSizes=readRDS(paste0("../specification_files/option",sampleSize_option,"_fullPop_genomeSampleSizes.RDS")) # these are for number of chromosomes sampled (2*number of indiv)
# make population order of sampleSizes correspond to order of allFreqs
sampleSizes=sampleSizes[colnames(allFreqs)]

rownames(allFreqs)=NULL


## Remove monomorphic sites: all positions in which global pop frequencies = 0 or all = 1
keep_rows=apply(allFreqs,1, function(row) !(all(row==0,na.rm=T) | all(row==1,na.rm=T)))
allFreqs=allFreqs[keep_rows,]

# randomize allele frequencies: get matrix with dim numPops x numLoci
rand.freqs = apply(allFreqs, 1, function(my.freqs) {
  if(runif(1)<0.5){my.freqs<-1-my.freqs} ; my.freqs
})

pops = rownames(rand.freqs)

# get sum pi for all population pairs (just one matrix)
get_sum_pi = Vectorize(function(pop1,pop2){
  if(pop1==pop2){
    sampleSizes[pop1]/(sampleSizes[pop1]-1)*sum(rand.freqs[pop1,]*(1-rand.freqs[pop2,]) + rand.freqs[pop2,]*(1-rand.freqs[pop1,]) )
  } else {
    sum(rand.freqs[pop1,]*(1-rand.freqs[pop2,]) + rand.freqs[pop2,]*(1-rand.freqs[pop1,]) )
  }
})

keep_index = apply(rand.freqs[c('Africa','Europe','East_Asia'),],2,function(site_freqs) (any(site_freqs>0) & any(site_freqs<1)) | (any(site_freqs==0) & any(site_freqs==1)) ) &
              apply(rand.freqs,2,function(site_freqs){
                  if(sum(site_freqs>0,na.rm=T)>1 & sum(site_freqs<1,na.rm=T)>1){
                    return(T)
                  } else{
                    filt_site_freqs = site_freqs[!is.na(site_freqs)]
                    private_pop=which(filt_site_freqs>0 & filt_site_freqs<1)
                    if(length(private_pop)==0){return(T)}
                    if((filt_site_freqs[private_pop]<0.01 & sum(filt_site_freqs[-private_pop])==0) | (filt_site_freqs[private_pop]> 0.99 & sum(filt_site_freqs[-private_pop])==1)){
                      return(F)
                    } else{
                      return(T)
                    }
                  }
                } )
rand.freqs=rand.freqs[,keep_index]
sum_pi_matrix=outer(pops,pops,get_sum_pi)
dimnames(sum_pi_matrix)=list(pops,pops)

# get coancestry
F_estimate = 1-sum_pi_matrix/sum_pi_matrix['Neanderthal','Africa']
saveRDS(F_estimate,paste0(path,'NeutralF.RDS'))
