#!/usr/bin/env python
# coding: utf-8

import msprime, pyslim
import numpy as np
import matplotlib.pyplot as plt
import tskit
import pandas as pd
import argparse

# start by getting/defining command line arguments
parser=argparse.ArgumentParser(description="get population allele frequencies from neutral slim tree sequences")
parser.add_argument("--demographyVersion",required=True,type=str,help="This is the letter that specifies which set of demographic parameters were used")
parser.add_argument("--sampleSize",required=True,type=str,help="This is the number that specifies which set of sample sizes to use")
parser.add_argument("--run",required=True,type=str,help="This is the number that specifies which simulation run should be analyzed")
parser.add_argument("--scenarioName",required=True,type=str,help="This is the file name (just prefix) of the selection scenario used to generate simulations")
args=parser.parse_args()

# Load the .trees file
trees_path="selection_output/version"+args.demographyVersion+"/"
run_specifier=args.scenarioName+"_run"+args.run
ts = pyslim.load(trees_path+'orig_trees/'+run_specifier+".trees")
#ts = pyslim.load('/home/sivan/Documents/Ch1/scripts/slim_sims/practice/output/selectionPractice.trees')
#ts = pyslim.load('/home/sivan/Documents/Ch1/scripts/slim_sims/practice/output/selectionTest3.trees')

# Specify the paths for population names and sample sizes, which will be used later
spec_path="../specification_files/"
pops_file=spec_path+"version"+args.demographyVersion+"_pops.feather"
sampleSizes_file=spec_path+"option"+args.sampleSize+"_fullPop_indivSampleSizes.feather"


recap = ts.recapitate(recombination_rate=1e-8, Ne=1e4)
#recap.dump(trees_path+'recapitated/'+run_specifier+".trees")

# Remember individuals from the past (aka sampled indiv in the past)
# find out what these times are:
# for t in np.unique(mutated.individual_times):
#  print(f"There are {np.sum(mutated.individual_times == t)} individuals from time {t}.")
# Order of sampling (from present to past) is:
# 1. modern, 2. Steppe, 3. EF, 4. WHG, 5. UP, 6. sim start

# individual metadata
tables=recap.dump_tables()
ind_md=list(pyslim.extract_individual_metadata(tables)) # each index is one individual's metadata

# sel site index before we add additional mutations: (this will be updated when we add mutations)
selSiteIndex=0

# population ID and name ####
# since population IDs in slim have different indexing than python lists,
# define slim population IDs by name in pop_ID dictionary
# then later we can use this dictionary to find the name associated with that ID
# and store population level information according to these names
# (we want to avoid any sort of indexing strategy because of the risk of inconsistency)

# load population IDs and sample sizes from R to python
# this allows for consistency between what's specified in multiple scripts
# pops start as a data frame, which we then make a series by indexing column name pops, which we then make a dictionary
# keys correspond to python indexing (add 1 to get slim and R indexing)
pops=list(pd.read_feather(pops_file)['pops'])
#pops = ['Neanderthal','Africa','UP','East_Asia','Europe','WHG','EF','Steppe']
#sampleSize_indiv={'Neanderthal':1,'Africa':100,'UP':10,'East_Asia':100,'Europe':100,'WHG':40,'EF':80,'Steppe':10}

# sample size (number of *individuals* sampled)
# load a data frame with 2 columns: 'pop' and 'size_indiv'
sampleSize_indiv_df=pd.read_feather(sampleSizes_file)
# turn the data frame into a dictionary with keys=popName and values=sampleSize
sampleSize_indiv=dict(zip(sampleSize_indiv_df['pop'],sampleSize_indiv_df['size_indiv']))

# ancient populations -- specifies which samples should be taken from paste
ancient_pops=['Neanderthal','EurUP','WHG','EF','Steppe']

# list of populations, in which each element is a list of individuals belonging to that population
# create empty list of length populations
pop_individuals=[ [] for i in range(len(pops)) ]

for t in np.unique(recap.individual_times)[:-1]: #[:-1] removes last item, which is the first gen in slim used for recapitation
    if t==0: # don't record ancient individuals in present day
        for indiv in recap.individuals_alive_at(t):
            popID=ind_md[indiv].population
            if pops[popID-1] not in ancient_pops: #popID-1 is the python index for the population
                pop_individuals[popID-1].append(indiv)
    else:
        for indiv in recap.individuals_alive_at(t):
            popID=ind_md[indiv].population
            pop_individuals[popID-1].append(indiv)

# turn pop_individuals into a dictionary with keys=pops and values=list of individual IDs
pop_individuals=dict(zip(pops,pop_individuals))

# take indiv samples from all possible individuals in each population
# dictionary: key = pop, value = list of individuals to sample
pop_indiv_samples={}
for pop in pops:
    pop_indiv_samples[pop]=np.random.choice(pop_individuals[pop],size=int(sampleSize_indiv[pop]),replace=False)

# each individual is two nodes (aka genomes); the nodes are what we sample in a tree
# so need to get list of nodes belonging to each population's sample
pop_node_samples={}
for pop in pops:
    if pop=='Africa' or pop=='Neanderthal':
        pop_node_samples[pop]=[]
        for indiv in pop_indiv_samples[pop]:
            pop_node_samples[pop].extend(recap.individual(indiv).nodes)
    elif pop in ancient_pops: # then assign nodes to dif genotype categories
        pop_node_samples[pop+'_AA']=[]
        pop_node_samples[pop+'_Aa']=[]
        pop_node_samples[pop+'_aa']=[]
        for indiv in pop_indiv_samples[pop]:
            ind_nodes=recap.individual(indiv).nodes
            # take sum of allele states at the two nodes (diploid genotype) to figure out indiv's pop assignment
            indiv_GT=sum(list(recap.variants(samples=ind_nodes))[selSiteIndex].genotypes)
            if indiv_GT==0: # no Neanderthal alleles
                pop_node_samples[pop+'_aa'].extend(ind_nodes)
            elif indiv_GT==1: # 1 Neanderthal allele
                pop_node_samples[pop+'_Aa'].extend(ind_nodes)
            else: # 2 Neanderthal alleles
                pop_node_samples[pop+'_AA'].extend(ind_nodes)
    else: #assign nodes to dif haplotype categories
        pop_node_samples[pop+'_AA']=[]
        pop_node_samples[pop+'_aa']=[]
        all_pop_node_samples=[]
        for indiv in pop_indiv_samples[pop]:
            all_pop_node_samples.extend(recap.individual(indiv).nodes)
        # in the order of nodes for each pop, get ancestral/derived state at sel site
        node_selSite_ancestry=list(recap.variants(samples=all_pop_node_samples))[selSiteIndex].genotypes
        # assign nodes to different partitioned populations according to ancestry at sel site
        # (subsetting by T/F of another vector only works when you make it a numpy array)
        pop_node_samples[pop+'_AA'].extend(np.array(all_pop_node_samples)[node_selSite_ancestry==1])
        pop_node_samples[pop+'_aa'].extend(np.array(all_pop_node_samples)[node_selSite_ancestry==0])

# clean up pop_node_samples dict by removing key,value pairs with empty values
# i.e. remove partitioned populations without any samples
pop_node_samples={k: v for k, v in pop_node_samples.items() if len(v) is not 0}

# get set of (partitioned) pops that we're working with
partPops=list(pop_node_samples.keys())

# SIMPLIFY ACCORDING TO NODES OF INTEREST

# first collapse all node samples into a single list, but track which index belongs to which partitioned population
all_pop_node_samples_recap=[]
for v in pop_node_samples.values():
	all_pop_node_samples_recap.extend(v)


tree_justNodes=recap.simplify(samples=all_pop_node_samples_recap,map_nodes=True,reduce_to_site_topology=False, filter_populations=False,
                          filter_individuals=False, filter_sites=False)
# with map_nodes=True, it returns an array with [0] tree sequence [2] numpy array that indicates the new sample ID of the node samples that previously existed
# below are the new node sample IDs corresponding to each pop
pop_node_samples_simple={k: tree_justNodes[1][v] for k,v in pop_node_samples.items()}



# ADD MUTATIONS
# lose info about ancestral samples if you do  below:
#mutated = msprime.mutate(ts, rate=1e-7, keep=True)
#mutated = pyslim.SlimTreeSequence(msprime.mutate(tree_justNodes[0], rate=1e-10, keep=True))
mutated = pyslim.SlimTreeSequence(msprime.mutate(tree_justNodes[0], rate=1e-7, keep=True))
#mutated.dump(trees_path+'mutated/'+run_specifier+".trees")


# get population SAMPLE allele frequencies
# the trick to getting genotypes at each site from each population is to subset samples from each pop one at a time,
# then get its genotype matrix (really just every sample's allelic state at each position)
# genotypes is list of pops, in which each item is list of lists: top level: site; lower level: sample states at a site
# each genotype matrix returned is the same length (same number of sites)

frequencies={}
for pop in partPops:
    # list of each site's allele freq in that population
    # iterate over variant generator to get population genotypes at each site (they say it's better to do it sequentially rather than all at once)
    frequencies[pop]=[sum(var.genotypes)/len(var.genotypes) for var in mutated.variants(samples=pop_node_samples_simple[pop]) ]

# get true selSiteIndex now that there are mutations
# need to find index of selected site
# because we want to know which position in a haplotype corresponds to the one m2 (selected) mutation
# this is the only site with mutation metadata, so can find its index this way:
is_SelSite=[]
for mut in pyslim.extract_mutation_metadata(mutated.dump_tables()):
    is_SelSite.append(len(mut))

selSiteIndex=is_SelSite.index(1) # sel site is the only one with mutation length 1, everything else is 0

# get sample sizes for each partitioned population
# when we save a named pandas series (soon), it loads in R as a named vector
sampleSizes_partPops_dict={k:len(v) for k,v in pop_node_samples_simple.items()}
sampleSizes_partPops=pd.Series(sampleSizes_partPops_dict)

# get full population (not partitioned) SAMPLE final frequencies at the selected site
selSite_finalFreqs_sampled={} # dictionary of each full population's Neanderthal allele frequency at sel site
for full_pop in pops:
	subPops=[pop for pop in partPops if full_pop in pop ]
	proportions=sampleSizes_partPops[subPops] / sampleSizes_partPops[subPops].sum()
	selSite_finalFreqs_sampled[full_pop]=sum([frequencies[subpop][selSiteIndex]*proportions[subpop] for subpop in subPops])
# the alternative to above is to just count the number of indiv in each partPop, but eh it's too late I already wrote this

# get full population (not partitioned) TRUE final frequencies at the selected site
fullPop_all_nodes={}
for pop in pops:
    fullPop_all_nodes[pop]=[]
    for indiv in pop_individuals[pop]:
        fullPop_all_nodes[pop].extend(recap.individual(indiv).nodes)

selSite_finalFreqs_true={} # dictionary of each full population's Neanderthal allele frequency at sel site
for full_pop in pops:
    selSite_info=list(recap.variants(samples=fullPop_all_nodes[full_pop]))[0]
    selSite_finalFreqs_true[full_pop]=sum(selSite_info.genotypes)/len(selSite_info.genotypes)


# get positions, in which index corresponds to frequencies in each pop
# Note: positions refer to base pair coordinates
# Note: don't really need positions for neutrality

# prep all data to be a pandas data frame
freq_df=pd.DataFrame(frequencies) # first get pandas data frame

positions=[site_md.position for site_md in list(mutated.sites())]
positions=pd.DataFrame(positions,columns=["pos"])

sampleSizes_partPops=pd.DataFrame([sampleSizes_partPops_dict.values()],columns=list(sampleSizes_partPops_dict.keys())).astype(int)
selSiteIndex=pd.DataFrame([selSiteIndex],columns=["ss_index"])

selSite_finalFreqs_sampled=pd.DataFrame([selSite_finalFreqs_sampled.values()],columns=list(selSite_finalFreqs_sampled.keys()))
selSite_finalFreqs_true=pd.DataFrame([selSite_finalFreqs_true.values()],columns=list(selSite_finalFreqs_true.keys()))


freq_df.to_feather(trees_path+'ftr_files/ss'+args.sampleSize+'_'+run_specifier+'_freqs.ftr')
positions.to_feather(trees_path+'ftr_files/ss'+args.sampleSize+'_'+run_specifier+'_positions.ftr')
sampleSizes_partPops.to_feather(trees_path+'ftr_files/ss'+args.sampleSize+'_'+run_specifier+'_sampleSizes.ftr')
selSite_finalFreqs_sampled.to_feather(trees_path+'ftr_files/ss'+args.sampleSize+'_'+run_specifier+'_sampled_finFreqs.ftr')
selSite_finalFreqs_true.to_feather(trees_path+'ftr_files/ss'+args.sampleSize+'_'+run_specifier+'_true_finFreqs.ftr')
selSiteIndex.to_feather(trees_path+'ftr_files/ss'+args.sampleSize+'_'+run_specifier+'_selSiteIndex.ftr')
