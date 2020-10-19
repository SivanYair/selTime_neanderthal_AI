This folder contains files specifying windows to analyze for a region with signals
of adaptive introgression.

**Format**: each file contains tab separated columns in the following order:
  1. chromosome
  2. left position
  3. right position
  4. region label

**partition-ranges-merged.txt**: lists windows for the 36 regions of adaptive introgression.

**inference-ranges-merged-2e-2.txt**: lists windows that we run the method on for a region, i.e. each region of adaptive introgression listed in the previous file,
plus 2cM surrounding it (1 cM on each side). There are fewer regions listed in this file because we merged regions with overlapping windows. We still run the method separately on each, we just merge them for efficiency when imputing or calculating allele frequencies.  
