
# general objects needed to calculate CLs

# distance of each site to sel site
distances=abs(positions-positions["sel_site"])
# distance bin of each site: tells us the index of the list of F^(S) or F^(N) for
# each genetic distance that the site of interest corresponds to
distBins = as.numeric(cut(distances,my.seq))

# MVN parameters
# for mean centering matrices and frequencies
M = length(real_pops) #number of pops for this observation
Tmatrix = matrix(data = rep(-1 / M, (M - 1) * M), nrow = M - 1, ncol = M)
diag(Tmatrix) = (M - 1) / M
k=nrow(Tmatrix)
mu = as.matrix(rep(0, k))

# mean-cenetered frequencies
freqs_MC = apply(allFreqs,2, function(site_freqs) Tmatrix %*% site_freqs)
epsilons = colMeans(allFreqs,na.rm=T)

# sample error to add along diagonal of predicted F^(N) or F^(S)
sampleError=diag((1/sampleSizes),nrow = M, ncol = M)

## COMPOSITE LIKELIHOOD FOR NULL  MODEL
# returns a single value
Null_CL = sum(sapply(1 : length(distances), function(site) { # nrow(distances) = nrow(sampleSizes) -- just means to iterate over every position
  return(calcLikelihood_bin_null(site, F_Null, distBins, epsilons, freqs_MC, sampleError, Tmatrix,mu,k))
}))

print('finished null model CL calculation')

## COMPOSITE LIKELIHOODS FOR MODELS UNDER SELECTION
# saved as a matrix: rows correspond to proposed values of s, and columns 
#                    correspond to proposed values of t_b. Each element is 
#                    the composite likelihood for that parameter combination

SelModel_CLs = sapply(1:length(times_btwn), function(t_btwn){
  sapply(1 : length(sels), function(s) {
    # list of F^(S) matrices for each genetic distance bin under this combination of t_b and s
    Fs_allBins = selMatrices_all[[t_btwn]][[s]] 
    if( 'null' %in% Fs_allBins) {return(Null_CL) }
    if( 'invalid' %in% Fs_allBins) {return('invalid')}
    # vector of probability of vector of population allele frequencies at each neutral site
    likelihood_by_site = sapply(1 : length(distances), function(site) { # 1:length(distances) just means to iterate over every position
      calcLikelihood_bin_selection(site, Fs_allBins, distBins, epsilons, freqs_MC, sampleError, Tmatrix,mu,k)
    })
    return(sum(likelihood_by_site))
  })
})
dimnames(SelModel_CLs)=list(sels,times_btwn)

print('finished selection model CL calculations')

# need to save the results to specific files that we can reload to plot them

all_results=append(list(Null_CL),list(SelModel_CLs))
names(all_results)[1] = 'Null'
names(all_results)[2] = 'Selection'

saveRDS(all_results,CL_file)
