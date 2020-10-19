calcLikelihood_bin_null = function(Site, F_N, distBins, epsilons, freqs_MC, sampleError,Tmatrix,mu,k) {
  # Calculates log-likelihood of data at a given position for null model
  #
  # Args:
  #	site: index of "positions" that the neutral site corresponds to
  # F_N: the list of F^(N) matrices for each genetic distance supplied by mid_genetic_distances
  # distBins: same as distBins defined in parent environment
  # epsilons: mean (proxy for ancestral) allele frequency at each position
  # freqs_MC: mean-centered allele frequencies at each position (i.e. original freq - epsilon)
  # sampleError, Tmatrix, mu, k: same as these objects defined in parent environment
  #
  # Returns:
  #	log-likelihood of s=0 at a given position
  #
  # Note: see multivariate normal probability density function to understand below
  #       likelihood calculation
  
  bin = distBins[Site] # index of mid_genetic_distances / of list of F^(N) that tells us which genetic distance bin to get predicted matrix
  my.e = epsilons[Site]*(1 - epsilons[Site]) # x_a*(1-x_a) in covariance calculation
  my.x = as.matrix(freqs_MC[ , Site]) # observed population allele frequencies
  
  #Add sample error to F^(N) matrix, then mean-center the matrix
  Fs_MC=(Tmatrix %*% (F_N[[bin]] + sampleError) %*% t(Tmatrix))
  
  # determinant and inverse of mean-cenetered and sample-size-corrected matrix
  det_F_N=det(Fs_MC)
  inv_F_N=ginv(Fs_MC)
  
  likelihood = 1 / (sqrt((2 * pi)^k * (det_F_N * my.e^k))) * 
    exp(-1 / 2 * t(my.x - mu) %*% (inv_F_N / my.e) %*% (my.x - mu))
  
  return(log(likelihood))

}

calcLikelihood_bin_selection = function(Site, F_S, distBins, epsilons, freqs_MC, sampleError,Tmatrix,mu,k) {
  # Calculates log-likelihood of data at a given position for selection model
  # under a certain combination of t_b and s
  #
  # Args:
  #	site: index of "positions" that the neutral site corresponds to
  # F_S: the list of F^(S) matrices for each genetic distance supplied by mid_genetic_distances
  # distBins: same as distBins defined in parent environment
  # epsilons: mean (proxy for ancestral) allele frequency at each position
  # freqs_MC: mean-centered allele frequencies at each position (i.e. original freq - epsilon)
  # sampleError, Tmatrix, mu, k: same as these objects defined in parent environment
  #
  # Returns:
  #	log-likelihood of s=0 at a given position
  #
  # Note: see multivariate normal probability density function to understand below
  #       likelihood calculation
  
  bin = distBins[Site] # index of mid_genetic_distances / of list of F^(N) that tells us which genetic distance bin to get predicted matrix
  my.e = epsilons[Site]*(1 - epsilons[Site]) # x_a*(1-x_a) in covariance calculation
  my.x = as.matrix(freqs_MC[ , Site]) # observed population allele frequencies
  
  #Add sample error to F^(S) matrix, then mean-center the matrix
  Fs_MC=(Tmatrix %*% (F_S[[bin]] + sampleError) %*% t(Tmatrix))
  
  # determinant and inverse of mean-cenetered and sample-size-corrected matrix
  det_F_S=det(Fs_MC)
  inv_F_S=ginv(Fs_MC)
  
  likelihood = 1 / (sqrt((2 * pi)^k * (det_F_S * my.e^k))) * 
    exp(-1 / 2 * t(my.x - mu) %*% (inv_F_S / my.e) %*% (my.x - mu))
  
  return(log(likelihood))
}