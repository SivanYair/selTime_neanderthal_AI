# Script to write slim script
# Required input is vector describing populations, splits, sample sizes
# I'll need to make some of this python readable so that it can be used there

# note: simulating 2 Mb of sequence with recombination rate of 1e-8 (to get 2e-2 Morgans of sequence) is hardcoded in script

args_file=commandArgs(trailingOnly = T)
if(length(args_file)!=1) stop('Need to supply file with arguments for make_selection_slim_script.R, with following 5 lines: demography version, time btwn, sel coef, fin freq, selPop Names (space delimited)')


args=readLines(args_file)

versionLetter=args[1] # letter describing the version of demography being used
t_btwn=as.numeric(args[2])
sel_coef=2*as.numeric(args[3]) # I define selection with 1,1+s,1+2s while slim defines it as 1,1+0.5s,1+s
fin_freq=as.numeric(args[4])
selPops=strsplit(args[5],split=" ")[[1]]

options(stringsAsFactors = F)

scenario=tail(strsplit(gsub(".txt","",args_file),split="/")[[1]],n=1)

filename=paste0('slim_scripts/version',versionLetter,'/',scenario,'.slim')

# File to output tree sequences is expected to be specified on slim command
# SLiM command should include -d trees_file="'${trees_file}'"

load(paste0('../specification_files/version',versionLetter,'_demography.RData'))

# define admixed pops from pops object defined in version demography
admixed_pops=setdiff(pops,c('Neanderthal','Africa'))
presentDay_admixed=c('Europe','East_Asia')

# make sure selected populations provided make sense according to demography

# get div times from beginning of simulation, assuming that
# the oldest divergence occurs when the simulation begins
# i.e. MH/Neanderthal split begins the simulation
# sim_duration is the full duration of simulation: starts at first pop split and lasts length of time specified by oldest divergence
sim_duration=2+div_times_from_present[2]
div_times_from_sim_start=sapply(div_times_from_present,function(t) sim_duration-t)
t_sampled_from_sim_start=sapply(t_sampled_from_present,function(t) sim_duration-t)
migs$start_from_sim_start=sapply(migs$start_from_present,function(t) sim_duration-as.integer(t))
migs$end_from_sim_start=sapply(migs$end_from_present,function(t) sim_duration-as.integer(t))


initialize_params=paste("initialize() {",
                        "\tinitializeTreeSeq();",
                        "\tinitializeMutationRate(0);",
                        '\tinitializeMutationType("m1", 0.5, "f", 0.0);',
                        '\tinitializeMutationType("m2", 0.5, "f", 0.0); // sweep mutation (additive selective effect), initially without a selective advantage',
                        '\tinitializeGenomicElementType("g1", m1, 1.0);',
                        '\tinitializeGenomicElement(g1, 0, 2e6);', # CHANGE AFTER TESTING
                        '\tinitializeRecombinationRate(1e-8);',
                        '}',
                        sep='\n')
write(file=filename,initialize_params)

# create first population (to be Neanderthals) & save run identifier for restarting based on derived allele freq trajectory
make_first_pop=paste('1 early() {',
                     "\t// save this run's identifier, used to save and restore",
                     '\tdefineConstant("simID", getSeed());',
                     '\tsim.addSubpop("p1", 10000);',
                     '}',
                     sep='\n')
write(file=filename,make_first_pop,append=T)

# create second population:ancestor of MH -- this specifies first population split
make_second_pop=paste('2 early() {',
                      '\tsim.addSubpopSplit("p2",10000,p1); // p2 currently represents ancestor of all MH, p1=Neanderthal',
                      '\tp1.setSubpopulationSize(2500);',
                      '}',
                      sep='\n')
write(file=filename,make_second_pop,append=T)

# create mutation in Neanderthals
introduce_mutation=paste('// introduce derived mutation in Neanderthals',
                         '3 late() {',
                         '\ttarget = sample(p1.genomes, 1);',
                         '\ttarget.addNewDrawnMutation(m2, 1e6);', # CHANGE AFTER TESTING
                         '\tsim.mutationsOfType(m2).setSelectionCoeff(0.05);',
                         '\t// save the state of the simulation',
                         '\tsim.treeSeqOutput("/tmp/slim_" + simID + ".trees");',
                         '}',
                         sep='\n')
write(file=filename,introduce_mutation,append=T)

# check that the mutation in Neanderthals is established (i.e. will reach fixation) before admixture (stop checking 2 gens before admixture)
establish_mutation=paste(paste0('3:',migs$start_from_sim_start[migs$donor=='Neanderthal']-2,' late() {'),
                         '\tmut = sim.mutationsOfType(m2);',
                         '\tif (size(mut) == 1)',
                         '\t{',
                         '\t\tif (sim.mutationFrequencies(p1, mut) > 0.1)',
                         '\t\t{',
                         '\t\t\tcatn("ESTABLISHED in Neanderthals");',
                         '\t\t\tsim.deregisterScriptBlock(self);',
                         '\t\t}',
                         '\t}',
                         '\telse',
                         '\t{',
                         '\t\tcatn("LOST FROM NEANDERTHALS – RESTARTING");',
                         '\t\t// go back to generation 3 & start a newly seeded run (but preserving tree sequence)',
                         '\t\tsim.readFromPopulationFile("/tmp/slim_" + simID + ".trees");',
                         '\t\tsetSeed(getSeed() + 1);',
                         '\t\tsim.mutationsOfType(m2).setSelectionCoeff(0.05); // we need the reminder because it seems like the simulation forgets these parameters',
                         '\t}',
                         '}',
                         sep='\n')
write(file=filename,establish_mutation,append=T)

# in the generation before admixture, check that Neanderthals fixed for derived allele, and change selection coefficient to be neutral
prep_for_admixture=paste(paste0(migs$start_from_sim_start[migs$donor=='Neanderthal']-1,' late() {'),
                         '\tmut = sim.mutationsOfType(m2);',
                         '\tif (size(mut) == 1)',
                         '\t{',
                         '\t\tif (sim.mutationFrequencies(p1, mut) == 1.0)',
                         '\t\t{',
                         '\t\t\tcatn("FIXED in Neanderthals");',
                         '\t\t\tsim.mutationsOfType(m2).setSelectionCoeff(0.0);',
                         '\t\t}',
                         '\t\telse',
                         '\t\t{',
                         '\t\t\tcatn("NOT FIXED IN NEANDERTHALS – RESTARTING");',
                         '\t\t\t// go back to generation 3 & start a newly seeded run (but preserving tree sequence)',
                         '\t\t\tsim.readFromPopulationFile("/tmp/slim_" + simID + ".trees");',
                         '\t\t\tsetSeed(getSeed() + 1);',
                         '\t\t\tsim.mutationsOfType(m2).setSelectionCoeff(0.05); // we need the reminder because it seems like the simulation forgets these parameters',
                         '\t\t}',
                         '\t}',
                         '\telse',
                         '\t{',
                         '\t\tcatn("NOT FIXED IN NEANDERTHALS – RESTARTING");',
                         '\t\t// go back to generation 3 & start a newly seeded run (but preserving tree sequence)',
                         '\t\tsim.readFromPopulationFile("/tmp/slim_" + simID + ".trees");',
                         '\t\tsetSeed(getSeed() + 1);',
                         '\t\tsim.mutationsOfType(m2).setSelectionCoeff(0.05); // we need the reminder because it seems like the simulation forgets these parameters',
                         '\t}',
                         '}',
                         sep='\n')
write(file=filename,prep_for_admixture,append=T)


# now add selection, divergence, migration, sampling among remaining populations

# following vector describes what to write for different generations
# we can later order this vector according to those generations, and then write the to the slim script
generation_instructions=c()
generation=c() # corresponding generation to generation_instructions

# work through each set of divergence times
for(i in 3:length(nested_pops)){ #pop_set can be 1 or more populations diverging
  pop_set=nested_pops[[i]]
  timing=div_times_from_sim_start[i]
  actual_info=paste(sapply(1:length(pop_set),function(p){
    code=paste0('\tsim.addSubpopSplit("p', (i-1)+p,'",10000,p',i-1,');') #(i-1) is previous pop, i is start of current pop set, but index starts from 1
    comment=ifelse(i==length(nested_pops),
                   paste0('// p',(i-1)+p,'=',pop_set[p],', p',i-1,'=',nested_pops[[i-1]]),
                   paste0('// p',i,' currently represents ancestor of ',paste(unlist(nested_pops[i:length(nested_pops)]),collapse='+'),', p',i-1,'=',nested_pops[[i-1]]) )
    paste(code,comment)
  }),collapse='\n')
  generation_chunk=paste(
    paste(timing,'late() {'),
    actual_info,
    '}',
    sep='\n'
  )
  generation_instructions=c(generation_instructions,generation_chunk)
  generation=c(generation,timing)
}

# work through each set of ancient population sampling times
for(pop in names(t_sampled_from_sim_start)){
  timing=t_sampled_from_sim_start[pop]
  sim_popID=paste0('p',which(pops==pop))
  code=paste0('\tsim.treeSeqRememberIndividuals(',sim_popID,'.individuals);')
  comment=paste('// record samples of all',pop,'individuals alive at this time')
  generation_chunk=paste(
    paste(timing,'early() {'),
    paste(code,comment),
    '}',
    sep='\n'
  )
  generation_instructions=c(generation_instructions,generation_chunk)
  generation=c(generation,timing)
}

# work through all sets of migration
for(i in 1:nrow(migs)){
  # recipient & donor pop ID
  recip_popID=paste0('p',which(pops==migs$recipient[i]))
  donor_popID=paste0('p',which(pops==migs$donor[i]))

  # print the chunks corresponding to the start and end of migration separately
  ## start chunk:
  start_timing=migs$start_from_sim_start[i]
  code=paste0('\t',recip_popID,'.setMigrationRates(',donor_popID,',',migs$rate[i],');')
  comment=ifelse(nest_index[migs$recipient[i]]!=length(nested_pops) & start_timing<div_times_from_sim_start[nest_index[migs$recipient[i]]+1],
                 paste('// ancestor of',paste(unlist(nested_pops[nest_index[migs$recipient[i]]:length(nested_pops)]),collapse='+'),'gets migrants from',migs$donor[i]),
                 paste('//',migs$recipient[i],'gets migrants from',migs$donor[i]))

  # make the chunk to start migration
  # if Neanderthals are the donor, save the state of the simulation
  if(migs$donor[i]=='Neanderthal'){
    start_chunk=paste(
      paste(start_timing, 'late() {'),
      paste(code,comment),

      # include lines to save simulation state (we will revert to this if selected Neanderthal allele doesn't persist in MH pop)
      '\t// save the state of the simulation',
	    '\tsim.treeSeqOutput("/tmp/slim_" + simID + ".trees");',

      '}',
      sep='\n'
    )
  } else{
    start_chunk=paste(
      paste(start_timing, 'late() {'),
      paste(code,comment),
      '}',
      sep='\n'
    )
  }

  generation_instructions=c(generation_instructions,start_chunk)
  generation=c(generation,start_timing)

  ##end chunk:
  end_timing=migs$end_from_sim_start[i]
  code=paste0('\t',recip_popID,'.setMigrationRates(',donor_popID,',0);')
  comment=ifelse(nest_index[migs$recipient[i]]!=length(nested_pops) & start_timing<div_times_from_sim_start[nest_index[migs$recipient[i]]+1],
                 paste('// ancestor of',paste(unlist(nested_pops[nest_index[migs$recipient[i]]:length(nested_pops)]),collapse='+'),'no longer gets migrants from',migs$donor[i]),
                 paste('//',migs$recipient[i],'no longer gets migrants from',migs$donor[i]))
  # if Neanderthals are the donor, print allele freqs
  if(migs$donor[i]=='Neanderthal'){
    end_chunk=paste(
      paste(end_timing, 'late() {'),
      paste(code,comment),

      # include lines to print allele freqs
      '\t// frequencies after migration',
      '\tmut=sim.mutationsOfType(m2);',
      '\tcatn("ALELLE FREQS JUST AFTER ADMIXTURE: " + sim.generation);',
      '\tpopFreqs=NULL;',
      '\tfor(pop in sim.subpopulations)',
      '\t{',
      '\t\tpopFreqs=c(popFreqs,sim.mutationFrequencies(pop,mut));',
      '\t}',
      '\tcatn(popFreqs,sep=",");',

      '}',
      sep='\n')


  } else{
    end_chunk=paste(
      paste(end_timing, 'late() {'),
      paste(code,comment),
      '}',
      sep='\n')
  }


  generation_instructions=c(generation_instructions,end_chunk)
  generation=c(generation,end_timing)

}

# work through selection chunks ###
#    if the timing of one of these chunks coincides with an already written chunk, we insert this text in there

# step 1: make mutation advantageous in selected populations
# a) get time frame of selection
selection_onset_from_sim_start = migs$start_from_sim_start[migs$donor=='Neanderthal']+t_btwn #time selection begins (with respect to start of simulation)
# use expected selection time for final freq to stop selection
g=as.numeric(migs$rate[migs$donor=='Neanderthal'])
sel_time = round((1/(sel_coef/2))*log((fin_freq*(1-g))/(g*(1-fin_freq)))) #sel_coef/2 because that's the true sel coef that we use in theory
selection_end_from_sim_start = selection_onset_from_sim_start + round(sel_time)
if(selection_end_from_sim_start>sim_duration) selection_end_from_sim_start=sim_duration-1

# does selection occur in all populations?
# if so, can define global sel coef
# if not, need to use fitness callbacks
use_globalSelCoef = all(admixed_pops %in% selPops)

if(use_globalSelCoef){
  # introduce the selective advantage, first checking if the allele is still segregating in all selected populations
  # if it is, record its allele frequencies, make it advantageous, and save the state of the simulation

  introduce_selection=paste('// start selection in admixed populations',
                            paste0(selection_onset_from_sim_start,' late() {'),
                            '\tmut = sim.mutationsOfType(m2);',
                            '\t// check if all selected populations are still segregating for Neanderthal mutation',
                            '\tpopFreqs=NULL;',
                            '\tfor(pop in setDifference(sim.subpopulations,c(p1,p2)))',
                            '\t{',
                            '\t\tpopFreqs=c(popFreqs,sim.mutationFrequencies(pop,mut));',
                            '\t}',
                            '\tif ( all(popFreqs>0) ) ',
                            '\t{',
                            '\t\tcatn("STARTING SELECTION in MH populations");',
                            paste0('\t\tsim.mutationsOfType(m2).setSelectionCoeff(', sel_coef, ');'),
                            # choosing not to save the state of the simulation, because sometimes allele freqs are too low they'll never reach high enough freq before selection ends, and when we keep restarting
                            # '\t\t// save the state of the simulation',
                            # '\t\tsim.treeSeqOutput("/tmp/slim_" + simID + ".trees");',
                            '\t\t// record what the initial allele frequencies are in each population',
                            '\t\tpopFreqsAll=NULL;',
                            '\t\tfor(pop in sim.subpopulations )',
                            '\t\t{',
                            '\t\t\tpopFreqsAll=c(popFreqsAll,sim.mutationFrequencies(pop,mut));',
                            '\t\t}',
                            '\t\tcat("ALLELE FREQS AT SEL START (" + sim.generation + "): ");',
                            '\t\tcatn(popFreqsAll,sep=", ");',
                            '\t}',
                            '\telse',
                            '\t{',
                            '\t\tcatn("NOT SEGREGATING in all selected pops – RESTARTING");',
                            '\t\t// go back to generation of admixture & start a newly seeded run (but preserving tree sequence)',
                            '\t\tsim.readFromPopulationFile("/tmp/slim_" + simID + ".trees");',
                            '\t\tsetSeed(getSeed() + 1);',
                            '\t\tsim.mutationsOfType(m2).setSelectionCoeff(0.0); // we need the reminder because it seems like the simulation forgets these parameters',
                            '\t}',
                            '}',
                            sep='\n')
  #if(selection_onset_from_sim_start %in% generation) stop('overlapping generation blocks, maybe bad?')
  generation_instructions=c(generation_instructions,introduce_selection)
  generation=c(generation,selection_onset_from_sim_start)

  # record frequency trajectory and monitor selection
  record_freq_traj = paste('// record frequency trajectory in sel pops & monitor that selected allele is not lost from any sel pops',
                           paste0(selection_onset_from_sim_start+1,':',selection_end_from_sim_start-1,' late() {'),
                           '\tmut = sim.mutationsOfType(m2);',
                           '\tpopFreqs=NULL;',
                           '\tfor(pop in setDifference(sim.subpopulations,c(p1,p2)))',
                           '\t{',
                           '\t\tpopFreqs=c(popFreqs,sim.mutationFrequencies(pop,mut));',
                           '\t}',
                           '\tselpopIDs_string="p"+3:length(sim.subpopulations);',
                           '\tcatn(c("gen:"+sim.generation,selpopIDs_string+":"+popFreqs),sep=",");',
                           '\tif ( !all(popFreqs>0) ) ',
                           '\t{',
                           '\t\tcatn("Neanderthal allele will NOT PERSIST in selected populations");',
                           '\t\t// go back to generation in which selection started & begin a newly seeded run (but preserving tree sequence)',
                           '\t\tsim.readFromPopulationFile("/tmp/slim_" + simID + ".trees");',
                           '\t\tsetSeed(getSeed() + 1);',
                           '\t\tsim.mutationsOfType(m2).setSelectionCoeff(0.0); // we need the reminder because it seems like the simulation forgets these parameters',
                           '\t}',
                           '}',
                           sep='\n')
  generation_instructions=c(generation_instructions,record_freq_traj)
  generation=c(generation,selection_onset_from_sim_start+1)

  # before finishing sweep, make sure that there's at least one population with the allele at greater than 20% frequency
  finish_selection = paste('// finish selection if final frequencies are appropriate',
                           paste0(selection_end_from_sim_start,' late() {'),
                           '\tmut=sim.mutationsOfType(m2);',
                           '\tpopFreqs=NULL;',
                           '\tfor(pop in setDifference(sim.subpopulations,c(p1,p2)))',
                           '\t{',
                           '\t\tpopFreqs=c(popFreqs,sim.mutationFrequencies(pop,mut));',
                           '\t}',
                           '\tif ( any(popFreqs>0.2 ) )',
                           '\t{',
                           '\t\tsim.mutationsOfType(m2).setSelectionCoeff(0.0);',
                           '\t\tpopFreqsAll=NULL;',
                           '\t\tfor(pop in sim.subpopulations )',
                           '\t\t{',
                           '\t\t\tpopFreqsAll=c(popFreqsAll,sim.mutationFrequencies(pop,mut));',
                           '\t\t}',
                           '\t\tcat("ALLELE FREQS AT SEL FINISH (" + sim.generation + "): ");',
                           '\t\tcatn(popFreqsAll,sep=", ");',
                           '\t}',
                           '\telse',
                           '\t{',
                           '\t\tcatn("Neanderthal allele DID NOT reach high enough freq (>20%) in any selected populations");',
                           '\t\t// go back to generation selection started & begin a newly seeded run (but preserving tree sequence)',
                           '\t\tsim.readFromPopulationFile("/tmp/slim_" + simID + ".trees");',
                           '\t\tsetSeed(getSeed() + 1);',
                           paste0('\t\tsim.mutationsOfType(m2).setSelectionCoeff(', sel_coef, '); // we need the reminder because it seems like the simulation forgets these parameters'),
                           '\t}',
                           '}',
                           sep='\n')

  #if(selection_end_from_sim_start %in% generation) stop('overlapping generation blocks, maybe bad?')
  generation_instructions=c(generation_instructions,finish_selection)
  generation=c(generation,selection_end_from_sim_start)
} else {
  # introduce fitness callback instead -- this should be the only difference from setting up mutations
  # fitness callback should be for all generations that selection exists
  # should look like: 14006:14480 fitness(m2,p2) { return 1.0 <-- neutral; otherwise 1+s  } <-- I can choose whether I want callback for sel or nonsel pops (which is easier / has fewer callbacks?)

  # find out which subpopulations experience selection (i.e. which are specified to be selected, and what populations exist at the time of selection)
  # get effective selPops: those that exist for at least some of the duration of the sweep (those that don't exist haven't been created yet)
  # if selection starts when not all populations have diverged, determine if there should be selection in an ancestral population
  # note that this is complicated to account for populations that may have been selected against, and so aren't specified as selected pops
  # which populations are selected, and are they represented in ancestral selection group?
  # (i.e. have they been defined by time selection starts and ends?)

  # we have 3 categories of selected populations:
  # 1) effective_selPops: populations that are directly selected and remain so
  # 2) selPops_thenAgainst: populations that are directly selected and later selected against
  # 3) selAgainst_after: populations whose ancestors were directly selected, but this population must be sel against

  notDefined_selPops_end = setdiff(selPops,unlist(nested_pops[which(div_times_from_sim_start<selection_end_from_sim_start)]))
  notDefined_selPops_start = setdiff(selPops,unlist(nested_pops[which(div_times_from_sim_start<selection_onset_from_sim_start)]))

  # if some selected populations experienced selection in a proxy ancestral population:
  if(length(notDefined_selPops_start)>0){
    # need to show selection in that ancestor: determine proxy pop
    proxyPop = nested_pops[[min(nest_index[notDefined_selPops_start])-1]]
    # keep adding to proxyPop as long as it's represented by another ancestral population
    while( div_times_from_sim_start[nest_index[tail(proxyPop,1)]] > selection_onset_from_sim_start ) {
      proxyPop=c(proxyPop,nested_pops[[nest_index[tail(proxyPop,1)]-1]])
    }

    effective_selPops = setdiff(selPops,notDefined_selPops_end) # populations that are selected and remain selected (not sel against later)
    selPops_thenAgainst=setdiff(proxyPop,effective_selPops) # populations that are selected and later sel against

    # is there selection in an ancestral pop, but one of its descendants isn't seg for Neander allele?
    # make that descendant population a "selAgainst_after" pop
    selAgainst_after=setdiff( setdiff(admixed_pops, unlist(nested_pops[which(div_times_from_sim_start<selection_end_from_sim_start)])) ,notDefined_selPops_end)

  } else { # if selected populations can experience selection in their own population
    effective_selPops=selPops
    selPops_thenAgainst=NULL
    selAgainst_after=NULL
  }

  # in order to introduce and monitor selected allele frequencies, we need to know which directly selected
  # populations are introduced at each time period
  # which directly selected populations exist at time selection starts?
  # create list of length time periods to consider:
  #   each element is a vector of existing sel pops
  direct_selPops = c(effective_selPops,selPops_thenAgainst)
  # identify the different time periods in which selected populations might changed, based on when population's split,
  # but ignore time periods corresponding to pop splits that happened before selection
  time_periods=unique(setdiff(nest_index[direct_selPops][order(nest_index[direct_selPops])], which(div_times_from_sim_start < selection_onset_from_sim_start)))
  selPop_exist_times=c() # record start times of a group of existing pops
  selPop_existing=list() # record group of pops existing at each sel time

  # if no selection is ancestral, i.e. starts independently in each sel pop
  if(length(time_periods)==0){
    selPop_existing[[1]] = direct_selPops
    selPop_exist_times[1] = selection_onset_from_sim_start
    selPop_exist_times[2] = selection_end_from_sim_start

    # IDs of subpopulations that are selected -- NOTE: we can provide fitness callbacks for populations that don't exist yet
    effective_selPop_IDs = paste0('p',which(pops %in% effective_selPops))

    # introduce the selective advantage
    # make fitness callback for populations that are directly selected & remain selected
    fitness_callback_effective=paste('// fitness of selected populations for sweep duration',
                                     paste(sapply(effective_selPop_IDs, function(p) paste0(selection_onset_from_sim_start,':',selection_end_from_sim_start,' fitness(m2,',p,') { return 1.0 + ', sel_coef, ' ; }')  ),
                                           collapse='\n'),
                                     sep='\n')
    generation_instructions=c(generation_instructions,fitness_callback_effective)
    generation=c(generation,selection_onset_from_sim_start)

  } else {
    # get the directly selected populations in the time period at which selection began
    # this time period is not accounted for in vector "time_periods" because it precedes
    # the first pop split that occurs during selection
    addPopGroup=intersect(direct_selPops,unlist(nested_pops[1:(time_periods[1]-1)]))
    selPop_exist_times=c(selPop_exist_times,selection_onset_from_sim_start)
    selPop_existing=append(selPop_existing,list(addPopGroup))

    # now consider changes in selected populations corresponding to population splits
    for(i in 1:length(time_periods)){ # need to consider which pops to add, and which pops to remove (maybe later selected against?)
      # before adding existingPops to selPop_existing list,
      # check if there were any changes in who's still considered a selected pop
      # (perhaps with a new pop divergence, one of them is relieved from being ancestral proxy pop)
      remove_pops = intersect(selPops_thenAgainst,selPop_existing[[i]]) # selPop_existing[[i]] actually does correspond to the previous index
      addPopGroup= setdiff( c(selPop_existing[[i]], nested_pops[[time_periods[i]]] ), remove_pops)

      selPop_exist_times=c(selPop_exist_times,div_times_from_sim_start[time_periods[i]])
      selPop_existing=append(selPop_existing,list(addPopGroup))

    }

    selPop_exist_times=c(selPop_exist_times,selection_end_from_sim_start) # add ending time for easy reference

    # IDs of subpopulations that are selected -- NOTE: we can provide fitness callbacks for populations that don't exist yet
    effective_selPop_IDs = paste0('p',which(pops %in% effective_selPops))

    # introduce the selective advantage
    # make fitness callback for populations that are directly selected & remain selected
    fitness_callback_effective=paste('// fitness of selected populations for sweep duration',
                                     paste(sapply(effective_selPop_IDs, function(p) paste0(selection_onset_from_sim_start,':',selection_end_from_sim_start,' fitness(m2,',p,') { return 1.0 + ', sel_coef, ' ; }')  ),
                                           collapse='\n'),
                                     sep='\n')
    generation_instructions=c(generation_instructions,fitness_callback_effective)
    generation=c(generation,selection_onset_from_sim_start)

    # make fitness callback(s) for populations that are directly selected and then selected against when they diverge from the rest
    for(pop in selPops_thenAgainst){
      start_deleterious=div_times_from_sim_start[nest_index[pop]+1]+1 # generation in which mutation switches from being beneficial to deleterious in this population

      stop_bene=ifelse(start_deleterious > selection_end_from_sim_start, selection_end_from_sim_start, start_deleterious )
      p=paste0('p',which(pops==pop))
      fitness_callback=paste(paste('// fitness of',pop,'for duration of sweep that it represents directly selected ancestral pop'),
                             paste0(selection_onset_from_sim_start,':',stop_bene,' fitness(m2,',p,') { if(homozygous) return 1.0 + ', sel_coef, ' ; else return 1.0 + mut.mutationType.dominanceCoeff*',sel_coef,'; }'),
                             paste('// fitness of',pop,'(sel against) after it diverges from directly selected ancestral pop'),
                             paste0(start_deleterious+1,':',sim_duration,' fitness(m2,',p,') { return 1.0 - 0.01 ; }'), # allow it to always be deleterious, despite homozygous/heterozygous
                             sep='\n')
      generation_instructions=c(generation_instructions,fitness_callback)
      generation=c(generation,selection_onset_from_sim_start)
    }

    # make fitness callback(s) for populations whose ancestors are directly selected, but they experience selection against when they diverge
    for(pop in selAgainst_after){
      start_deleterious=div_times_from_sim_start[nest_index[pop]+1]+1 # generation in which mutation switches from being  neutral to deleterious in this population
      p=paste0('p',which(pops==pop))
      fitness_callback=paste(paste('// fitness of',pop,'(sel against) after it diverges from directly selected ancestral pop'),
                             paste0(start_deleterious,':',sim_duration,' fitness(m2,',p,') { return 1.0 - 0.01 ; }'),
                             sep='\n')
      generation_instructions=c(generation_instructions,fitness_callback)
      generation=c(generation,selection_onset_from_sim_start)
    }

  }




  # don't allow selection to start unless allele is still segregating in all selected populations
  # if it is, record its allele frequencies, make it advantageous, and save the state of the simulation
  introduce_selection=paste('// check if we can start selection in admixed populations',
                            paste0(selection_onset_from_sim_start,' late() {'),
                            '\tmut = sim.mutationsOfType(m2);',
                            '\t// check if all selected populations are still segregating for Neanderthal mutation',
                            '\tpopFreqs=NULL;',
                            paste0('\tfor(pop in c(',paste(paste0('p',which(pops %in% selPop_existing[[1]])),collapse=','),') )'),
                            '\t{',
                            '\t\tpopFreqs=c(popFreqs,sim.mutationFrequencies(pop,mut));',
                            '\t}',
                            '\tif ( all(popFreqs>0) ) ',
                            '\t{',
                            '\t\tcatn(" STARTING SELECTION in MH populations");',
                            # choosing not to save the state of the simulation, because sometimes allele freqs are too low they'll never reach high enough freq before selection ends, and when we keep restarting
                            # '\t\t// save the state of the simulation',
                            # '\t\tsim.treeSeqOutput("/tmp/slim_" + simID + ".trees");',
                            '\t\t// record what the initial allele frequencies are in each population',
                            '\t\tpopFreqsAll=NULL;',
                            '\t\tfor(pop in sim.subpopulations )',
                            '\t\t{',
                            '\t\t\tpopFreqsAll=c(popFreqsAll,sim.mutationFrequencies(pop,mut));',
                            '\t\t}',
                            '\t\tcat("ALLELE FREQS AT SEL START (" + sim.generation + "): ");',
                            '\t\tcatn(popFreqsAll,sep=", ");',
                            '\t}',
                            '\telse',
                            '\t{',
                            '\t\tcatn("NOT SEGREGATING in all selected pops – RESTARTING");',
                            '\t\t// go back to generation of admixture & start a newly seeded run (but preserving tree sequence)',
                            '\t\tsim.readFromPopulationFile("/tmp/slim_" + simID + ".trees");',
                            '\t\tsetSeed(getSeed() + 1);',
                            '\t\tsim.mutationsOfType(m2).setSelectionCoeff(0.0); // we need the reminder because it seems like the simulation forgets these parameters',
                            '\t}',
                            '}',
                            sep='\n')

  #if(selection_onset_from_sim_start %in% generation) stop('overlapping generation blocks, maybe bad?')
  generation_instructions=c(generation_instructions,introduce_selection)
  generation=c(generation,selection_onset_from_sim_start)

  # record frequency trajectory and monitor selection
  record_freq_traj = paste('// record frequency trajectory in sel pops & monitor that selected allele is not lost from any sel pops',
                           paste( sapply(1:length(selPop_existing),function(i){
                             paste(paste0(selPop_exist_times[i]+1,':',selPop_exist_times[i+1]-1,' late() {'),
                                   '\tmut = sim.mutationsOfType(m2);',
                                   '\tpopFreqs=NULL;',
                                   paste0('\tfor(pop in c(',paste(paste0('p',which(pops %in% selPop_existing[[i]])),collapse=','),') )'),
                                   '\t{',
                                   '\t\tpopFreqs=c(popFreqs,sim.mutationFrequencies(pop,mut));',
                                   '\t}',
                                   paste0('\tselpopIDs_string=c(',paste(paste0('"p',which(pops %in% selPop_existing[[i]]),'"'),collapse=','),');'),
                                   '\tcatn(c("gen:"+sim.generation,selpopIDs_string+":"+popFreqs),sep=",");',
                                   '\tif ( !all(popFreqs>0) ) ',
                                   '\t{',
                                   '\t\tcatn("Neanderthal allele will NOT PERSIST in selected populations");',
                                   '\t\t// go back to generation in which selection started & begin a newly seeded run (but preserving tree sequence)',
                                   '\t\tsim.readFromPopulationFile("/tmp/slim_" + simID + ".trees");',
                                   '\t\tsetSeed(getSeed() + 1);',
                                   '\t\tsim.mutationsOfType(m2).setSelectionCoeff(0.0); // we need the reminder because it seems like the simulation forgets these parameters',
                                   '\t}',
                                   '}',
                                   sep='\n')
                           }), collapse='\n'),
                           sep='\n')
  generation_instructions=c(generation_instructions,record_freq_traj)
  generation=c(generation,selection_onset_from_sim_start+1)

  # before finishing sweep, make sure that there's at least one population with the allele at greater than 20% frequency
  finish_selection = paste('// finish selection if final frequencies are appropriate',
                           paste0(selection_end_from_sim_start,' late() {'),
                           '\tmut = sim.mutationsOfType(m2);',
                           '\tpopFreqs=NULL;',
                           paste0('\tfor(pop in c(',paste(paste0('p',which(pops %in% selPop_existing[[length(selPop_existing)]])),collapse=','),') )'),
                           '\t{',
                           '\t\tpopFreqs=c(popFreqs,sim.mutationFrequencies(pop,mut));',
                           '\t}',
                           '\tif ( !any(popFreqs>0.2 ) )',
                           '\t{',
                           '\t\tcatn("Neanderthal allele DID NOT reach high enough freq (>20%) in any selected populations");',
                           '\t\t// go back to generation selection started & begin a newly seeded run (but preserving tree sequence)',
                           '\t\tsim.readFromPopulationFile("/tmp/slim_" + simID + ".trees");',
                           '\t\tsetSeed(getSeed() + 1);',
                           '\t\tsim.mutationsOfType(m2).setSelectionCoeff(0.0); // we need the reminder because it seems like the simulation forgets these parameters',
                           '\t}',
                           '\telse',
                           '\t{',
                           '\t\tpopFreqsAll=NULL;',
                           '\t\tfor(pop in sim.subpopulations )',
                           '\t\t{',
                           '\t\t\tpopFreqsAll=c(popFreqsAll,sim.mutationFrequencies(pop,mut));',
                           '\t\t}',
                           '\t\tcat("ALLELE FREQS AT SEL FINISH (" + sim.generation + "): ");',
                           '\t\tcatn(popFreqsAll,sep=", ");',
                           '\t}',
                           '}',
                           sep='\n')
  #if(selection_end_from_sim_start %in% generation) stop('overlapping generation blocks, maybe bad?')
  generation_instructions=c(generation_instructions,finish_selection)
  generation=c(generation,selection_end_from_sim_start)

  # # find out if a population experienced selected in the ancestor, but then selection against Neander allele
  # # these are the introgressed populations not included in effective_selPops nor the non-effective sel pops
  # intro_pops=setdiff(pops,c('Neanderthal','Africa'))
  # setdiff(selPops,effective_selPops) # <-- these are the populations whose ancestors were selected
  #
  # # see which population(s) selection starts in:
  # orig_selPops=setdiff(unlist(nested_pops[div_times_from_sim_start<=selection_onset_from_sim_start]),'Africa')
  #
  # # identify descendants of those originally selected populations that shouldnt have Neanderthal allele at high freq:
  # selAgainst_pops=setdiff(unlist(nested_pops[div_times_from_sim_start>selection_onset_from_sim_start]),selPops)
  #
  #
  # # these are the populations whose ancestors may have been selected:
  # selAgainst_pops = setdiff(intro_pops,selPops)



}







#####   WRITE TO SLIM SCRIPT #######

# and write to the slim script according to timing (ascending order)
write(file=filename,generation_instructions[order(generation)],append=T,sep='\n\n')

# add the final chunk to finish simulation and record tree sequence to file
finish_chunk=paste(
  paste(sim_duration,'late() {'),
  '//make sure allele freqs are high enough in at least one presentDay, selected population before finishing',
  '\tmut = sim.mutationsOfType(m2);',
  '\tpopFreqs=NULL;',
  paste0('\tfor(pop in c(',paste(paste0('p',which(pops %in% unique(c(intersect(selPops,presentDay_admixed)),'Europe'))),collapse=','),') )'),
  '\t{',
  '\t\tpopFreqs=c(popFreqs,sim.mutationFrequencies(pop,mut));',
  '\t}',
  '\tif ( !any(popFreqs>0.2 ) )',
  '\t{',
  '\t\tcatn("AT SIM END, Neanderthal allele DID NOT reach high enough freq (>20%) in any selected populations -- RESTARTING SELECTION");',
  '\t\t// go back to generation selection started & begin a newly seeded run (but preserving tree sequence)',
  '\t\tsim.readFromPopulationFile("/tmp/slim_" + simID + ".trees");',
  '\t\tsetSeed(getSeed() + 1);',
  '\t\tsim.mutationsOfType(m2).setSelectionCoeff(0.0); // we need the reminder because it seems like the simulation forgets these parameters',
  '\t}',
  '\telse',
  '\t{',
  '\t//record final frequencies as sanity check',
  '\t\tcatn("FINISHED! ALELLE FREQS AT FINAL GEN: " + sim.generation);',
  '\t\tpopFreqsAll=NULL;',
  '\t\tfor(pop in sim.subpopulations)',
  '\t\t{',
  '\t\t\tpopFreqsAll=c(popFreqsAll,sim.mutationFrequencies(pop,mut));',
  '\t\t}',
  '\t\tcatn(popFreqsAll,sep=",");',
  '',
  '\t\tsim.treeSeqOutput(trees_file);', # trees_file should be specified in slim command
  '\t\tsim.simulationFinished();',
  '\t}',
  '}',
  sep='\n'
)
write(file=filename,finish_chunk,append=T)
