# Functions for generating matrices for probabilities of coalescing under 
# null (F^(N)) and selection (F^(S)) models

# How it works:
# Step 1: Define parameters of each function
# Step 2: Calculate F^(N) for range of representative genetic distances between a neutral and selected site
#         This returns F_Null, a list of matrices in which each element corresponds to a given genetic distance
# Step 3: Calculate F^(S) for each combination of t_b, s, and r (genetic distance). 
#         This returns selMatrices_all, a nested list of matrices:
#           - the first index corresponds to indices of times_btwn (vector of proposed t_b)
#           - the second index corresponds to indices of sels (vector of proposed s)
#           - the third index corresponds to indices of mid_genetic_distances (vector of representative genetic distances from the selected site)
#         For example, selMatrices_all[[5]][[3]][[900]] is a matrix for probabilities of coalescing for each population pair
#         for t_b=times_btwn[5], s=sels[3], and r=mid_genetic_distances[900]

#### DEFINE PARAMETERS OF FUNCTIONS F^(N) AND F^(S) (except admixture graph parameters that were previously defined) ####

##bin genetic distances
# creating absolute genetic distance bins based on min and max distance
# each neutral site can be from the selected site, as long as
# we simulate 2 Mb of sequence with the selected site at the center of the locus
# we want 1000 bins per genetic distance of 1 cM (that would be 1000 bins in a 1 Mb locus with uniform rec rate of 1e-8)
numBins = round((npos/2)*1e-8*1e5) #(=1e5 bins per genetic distance of 1 M = 1e3 bins per genetic distance of 1 cM)
my.seq = seq(-1, (npos/2 + 1), length.out = (numBins + 1))
midDistances = sapply(1:numBins, function(i) mean(c(my.seq[i], my.seq[i+1]))) # representative physical distances for each bin
mid_genetic_distances=rec_rate*midDistances # representative genetic distances for each bin

#proposed selection coefficients
sels = c(0.002, 0.005,0.01,0.015,0.02,0.025,0.03,0.05)

# proposed times between introgression and selection onset
times_btwn = seq(0,1700,by=100)

# x_s in models
sweep_finFreq = mean(fin_freqs[gsub("_..\\b","",potential_selPops)])

## Generate coancestry matrix under null model
F_Null=calc_F.N_all(F_est=F_estimate)

## below are used in function to generate F^(S)
# Neanderthal-admixed populations
intro_pops=setdiff(pops,c(source,neverSel,grep("_Aa",pops,value=T))) 
# probability of coalescing given both lineages did not descend from Neanderthals
.both_MH=outer(intro_pops,intro_pops, Vectorize(function(pop1,pop2) (F_estimate[pop1,pop2]-g^2*F_estimate[source,source])/(1-g)^2, USE.NAMES=F) )
dimnames(.both_MH)=list(intro_pops,intro_pops)

## Generate selection matrices for bins of genetic distances away from selected site
## for each combination of t_b and s
selMatrices_all = lapply(times_btwn,function(t_btwn){
  lapply(sels,function(S){
    calc_F.S_all(selPops=sel_pops, s=S,t_b=t_btwn,x=sweep_finFreq,migs=p)
  })
})
