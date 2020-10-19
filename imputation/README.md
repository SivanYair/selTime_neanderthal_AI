We show how to impute genotypes in our ancient samples. We provide genotype likelihood
files from the BNC2 region for our example.


Information about contents of this directory:

* **impute_ancient_fullPop.txt**: bash script that imputes genotypes for each ancient population

* **readcounts_byRegion_v4**: directory that contains a vcf for each region with each ancient
sample's read counts and genotype likelihoods at each site in the region.

* **pop-samples**: directory that contains text files listing the IDs of ancient samples that
belong to a certain population

* **1000g_refPanel**: directory that contains a vcf for 1000 Genomes European and East Asian samples used in the reference panel for imputation. 

* Before running impute_ancient_fullPop.txt you need to add a directory here called "plink.GRCh37.map", which contains genetic maps downloaded from [this Beagle page](http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/).

The regions imputed correspond to those listed in ../region-info/inference-ranges-merged-2e-2. This file lists the start and end positions of the window containing 2cM around the candidate region of adaptive introgression. They are
"merged" because some of these windows for analysis overlapped.
