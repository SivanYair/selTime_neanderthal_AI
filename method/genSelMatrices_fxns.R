
## These are a set of functions that are used in genSelMatrices_exec.R
## They generate predicted selection matrices for each model according to the
## distance bin away from selected site, neutral F matrix, admixture graph parameters,
## set of putative selected populations, final frequency of putative selected lineage,
## and the parameters we aim to infer, s (selection coefficient) and t_b (time between
## admixture and onset of selection)

### GENERAL FUNCTION TERMS USED LATER ###
# probability a pair of lineages coalesce, given that neither of them migrate back to Neanderthals
both_MH=function(pop1,pop2) (F_estimate[pop1,pop2]-g^2*F_estimate[source,source])/(1-g)^2

# Rx = prob of not recombining out from time sampled to time of sweep completion <-- varies among populations due to different sampling times (t_sampled)
Rx=function(pop) {
  if(t_int-t_b-unname(t_sampled[pop]) >= t_s ){
    exp(-r*(t_int-t_b-t_s-unname(t_sampled[pop])))
  } else{
    1
  }
} 

# probability that an lineage *sampled from the AA partition of a selected population* is linked to the selected variant at the time the sweep finishes (SF)
AA_selBG_SF=function(pop) Rx(pop) + (1-Rx(pop))*fin_freq
# probability that an lineage *sampled from the aa partition of a selected population* is linked to the selected variant at the time the sweep finishes (SF)
aa_selBG_SF=function(pop) (1-Rx(pop))*fin_freq

# probability that an lineage does not recombines out of its background from the time it was sampled to the time of admixture
p_nr_wholeTime=function(pop) exp(-r*(t_int-t_sampled[pop]))

# probability that a coalescence or recombination event occurs within a given amount of time,
# along with the conditional probability that the first event is coalescence or recombination
# returns vector with names: p_event('w' in math), p_coal, p_rec
p_events=function(pop1,pop2,time,bg_freq){ 
  ## args:
  ## pop1 and pop2: the two populations that the pair of lineages were sampled from
  ## time: amount of time that the event can occur
  ## bg_freq: frequency of the background that both lineages are associated with
  w=unname(1-exp(-time*(1/(2*mean(c(Ne[pop1],Ne[pop2]))*bg_freq) + 2*r)))
  p_coal=1/(1+4*mean(c(Ne[pop1],Ne[pop2]))*bg_freq*r)
  p_rec=(4*mean(c(Ne[pop1],Ne[pop2]))*bg_freq*r/(1+4*mean(c(Ne[pop1],Ne[pop2]))*bg_freq*r))
  return(c(w=w,p_coal=p_coal,p_rec=p_rec))
}

# probability a pair of lineages coalesce, given that they are linked to the selected variant
# at the time that the sweep finishes
P=function(pop1,pop2){
  # t_shared_g = the amount of time the selected variant is at freq g in the ancestral population of pop1 and pop2
  # t_nonShared_g = the amount of time the selected variant is at freq g when pop1 and pop2 are isolated/diverged from each other

  # If sweep starts in ancestor of pair of populations (lineages are in same population for whole time selected variant is at freq g):
  t_shared_g=ifelse((t_anc[pop1,pop2]-t_b)>0, t_b, t_anc[pop1,pop2])
  t_nonShared_g=t_b-t_shared_g # t_nonShared_g=0 if sweep starts in ancestor of pair of populations

  # nr_selBG = probability of never recombining out of selected variant background between
  #          sweep finish and time lineages are in same population
  #          this involves two time periods: during sweep & during time t_nonShared_g

  nr_selBG = Ry_selBG*exp(-r*t_nonShared_g)

  # this is the final term to return
  nr_selBG^2*(P.g(pop1,pop2,t_shared_g)) +
    2*nr_selBG*(1-nr_selBG)*(Rz*g*F_estimate[source,source] + (1-Rz)*F_estimate[pop1,pop2]) +
    (1-nr_selBG)^2*F_estimate[pop1,pop2]

}

# probability a pair of lineages coalesce, given that at the time that the sweep finishes,
# one is linked to the selected variant while the other is not
Q=function(pop1,pop2){
  Ry_selBG*Rz*Ry_nonSelBG*Rz*0 + #neither rec out of bg
    (1-Ry_selBG*Rz)*(1-g)*Ry_nonSelBG*Rz*both_MH(pop1,pop2) + #linked recombines out and know it's type MH, unlinked does not
    Ry_selBG*Rz*(1-Ry_nonSelBG*Rz)*g*F_estimate[source,source] + # linked does not rec out (type source), unlinked does and know it's type source
    (1-Ry_selBG*Rz)*(1-Ry_nonSelBG*Rz)*F_estimate[pop1,pop2] #both rec out
}

# probability a pair of lineages coalesce, given that at the time that the sweep finishes,
# both are not linked to the selected variant
U=function(pop1,pop2){
  Ry_nonSelBG*Rz^2*both_MH(pop1,pop2) + (1-Ry_nonSelBG*Rz)^2*F_estimate[pop1,pop2] + 2*Ry_nonSelBG*Rz*(1-Ry_nonSelBG*Rz)*(1-g)*both_MH(pop1,pop2)
}

# probability a pair of lineages coalesce, if they both recombine out at least once
# between time of sampling and the sweep finishing
p_coal_anyBG=function(pop1,pop2){
  fin_freq^2*P(pop1,pop2) + (1-fin_freq)^2*U(pop1,pop2) + 2*fin_freq*(1-fin_freq)*Q(pop1,pop2)
}

# used for relationships with non-selected pops:
# probability a pair of lineages coalesce, given that they are segregating in the same population
# and are linked to the selected variant when it's at freq g,
# for the given amount of time the lineage is at this frequency before admixture
P.g=function(pop1,pop2,time){
  # t_shared_g = the amount of time the selected variant is at freq g in the ancestral population of pop1 and pop2
  t_shared_g=time

  p_events_b4_admix=p_events(pop1,pop2,t_shared_g,g)

  # this is the final term to return
  (1-p_events_b4_admix['w'])*F_estimate[source,source] + p_events_b4_admix['w']*(p_events_b4_admix['p_coal'] +
                                                                                   p_events_b4_admix['p_rec']*(Rz*g*F_estimate[source,source] + (1-Rz)*F_estimate[pop1,pop2]) )

}


### SELECTION MODEL ###
calc_F.S_all = function(selPops, s,t_b,x,migs) {
  # Computes list of F^(S) for every genetic distance bin under the parameter combination defined by function input
  #
  # Args:
  # selPops: putative selected populations
  # s: selection coefficient
  # t_b: time between admixture and selection onset
  # x: final frequency of the selected variant
  # migs: matrix of proportion ancestry that column pop j contributes to row pop i, among Neanderthal-admixed populations
  
  # Returns list of length numBins because F^(S) will differ
  # with distance of the neutral site to the selected site

  # duration of sweep phase
  if(x==1){
    t_s = (1/s)*log((2*mean(Ne[selPops])-1)*(1-g)/g)
  } else{
    t_s = (1/s)*log((x*(1-g))/(g*(1-x)))
  }
  
  # parameter combinations that do not make sense
  if(t_s > t_int){return('null')} # if sweep phase duration exceeds time until introgression, note that the null model applies
  if(t_s > t_int-t_b){return('invalid')} # if sweep phase duration exceeds allotted time before present, don't calculate F^(S)

  # determine who the selected populations are, based on time between sweep and selection:
  # determine which pops have selection ancestrally shared selection, which of those may have subsequently lost the selected lineage, and which pops have independent selection
  shared_time = anc_times - t_b
  shared_time[shared_time<0]=0
  if(all(shared_time==0)){
    selPops_anc_AA=c()
    selPops_ind_AA=selPops
    selPops_thenAgainst_aa=c()
  } else {
    # see who shares the sweep ancestrally
    anc_nest_index = min(which(shared_time>0))
    selPops_anc_whole = unique(gsub('_..$','',unlist(nested_pops[anc_nest_index:length(nested_pops)])))
    if(any(sapply(selPops_anc_whole,function(pop) grepl(pop,selPops)))){
      selPops_anc_AA = paste0(selPops_anc_whole[sapply(selPops_anc_whole,function(whole_pop) paste0(whole_pop,'_AA') %in% pops )],'_AA')
      sel_against_whole=selPops_anc_whole[sapply(selPops_anc_whole,function(whole_pop) !(paste0(whole_pop,'_AA') %in% pops) )]
      if(length(sel_against_whole) > 0){
        selPops_thenAgainst_aa = paste0(sel_against_whole,'_aa')
        selPops_thenAgainst_aa = unname(sapply(selPops_thenAgainst_aa,function(pop){
          if(pop %in% pops){
            return(pop)
          } else{
            return(gsub("_aa",'',pop))
          }
        }))
      } else{
        selPops_thenAgainst_aa=c()
      }
    } else{
      selPops_anc_AA=c()
      selPops_thenAgainst_aa=c()
    }
    selPops_ind_AA=setdiff(selPops,c(selPops_anc_AA,selPops_thenAgainst_aa))
  }
  
  # remove a selected population if it was sampled after selection started
  selPops_ind_AA=selPops_ind_AA[t_sampled[selPops_ind_AA] < t_int - t_b ]
  
  # remove non-existent populations
  selPops_anc_aa = gsub('AA','aa',selPops_anc_AA)
  selPops_anc_aa=selPops_anc_aa[selPops_anc_aa %in% pops] 
  selPops_ind_aa = gsub('AA','aa',selPops_ind_AA)
  selPops_ind_aa=selPops_ind_aa[selPops_ind_aa %in% pops] 
  
  # are there any aa partitions of sel pops to keep track of?
  # need to check or else we get an error with rbind() later
  selPops_aa_exist=ifelse(length(c(selPops_anc_aa,selPops_ind_aa))==0,F,T)
  
  intro_sel_pops = c(selPops_anc_AA,selPops_ind_AA,selPops_ind_aa,selPops_anc_aa,selPops_thenAgainst_aa)
  n_intro_sel_pops=length(intro_sel_pops)
  
  # partitioned populations sampled from populations that were introgressed but did not experience selection
  noSel_AA=pops[ (grepl('_AA',pops)) & !(pops %in% c(source,neverSel,selPops_anc_AA,selPops_ind_AA,selPops_thenAgainst_aa) )]
  noSel_aa=pops[ ((grepl('_aa',pops)) & !(pops %in% c(source,neverSel,selPops_anc_aa,selPops_ind_aa,selPops_thenAgainst_aa) )) | (pops%in%intro_sel_pops & !(grepl("_..$",pops)) & !(pops %in% selPops_thenAgainst_aa))]
  
  calc_F.S = function(r) {
    # Computes F^(S) for a single genetic distance bin and parameter combination provided in containing function
    #
    # Args:
    # r: recombination rate (genetic distance) between neutral and selected site
    
    #### step 1: define a matrix of within & between pop coancestry for _AA and _aa pops
    
    # partitioned populations corresponding to selected populations
    # assumes there's always an 'aa' population corresponding to 'AA' population
    
    # probability of not recombining out of background during neutral phase I
    Rz=exp(-r*(t_b))
    
    # probability of never recombining onto other background during sweep phase
    if(x==1){
      # probabilities of *not* recombining onto
      # the other background; depends on initial background
      # options: linked to the sel lineage, or not linked to the sel lineage
      Ry_selBG=(x*(1-g)/g)^(-r/s)
      Ry_nonSelBG=(2*mean(Ne[selPops]))^(-r/s) 
    } else{
      # probabilities of *not* recombining onto
      # the other background; depends on initial background
      # options: linked to the sel lineage, or not linked to the sel lineage
      Ry_selBG=(x*(1-g)/g)^(-r/s)
      Ry_nonSelBG=(1/(1-x))^(-r/s)
    }
    
    # make a vector of these probabilities because sometimes they need to be modified
    Ry_selBG = rep(Ry_selBG,n_intro_sel_pops)
    names(Ry_selBG) = intro_sel_pops
    Ry_nonSelBG = rep(Ry_nonSelBG,n_intro_sel_pops)
    names(Ry_nonSelBG) = intro_sel_pops
    

    if(any(anc_times[nest_index[selPops_thenAgainst_aa]] < t_b+t_s)  ){ # this returns FALSE if there are no selAgainst_pops
      # need to modify Ry_selBG and Ry_nonSelBG for these populations if at least one
      # of them split off during the sweep (they have a different final freq so a 
      # different probability of staying linked)
      
      change_pops = selPops_thenAgainst_aa[anc_times[nest_index[selPops_thenAgainst_aa]] < t_b+t_s]
      
      Ry_selBG_selAgainst_fxn = Vectorize(function(pop){
        new_ts = anc_times[nest_index[pop]] - t_b
        new_x = g*exp(s*new_ts) / (1-g + g*exp(s*new_ts))
        return((new_x*(1-g)/g)^(-r/s))
      },USE.NAMES = F)
      
      Ry_nonSelBG_selAgainst_fxn = Vectorize(function(pop){
        new_ts = anc_times[nest_index[pop]] - t_b
        new_x = g*exp(s*new_ts) / (1-g + g*exp(s*new_ts))
        return((1/(1-new_x))^(-r/s))
      },USE.NAMES = F)
      
      Ry_selBG[change_pops] = Ry_selBG_selAgainst_fxn(change_pops)
      Ry_nonSelBG[change_pops] = Ry_nonSelBG_selAgainst_fxn(change_pops)
      
    }
    
    # if sweep ends before sampling, we lower final frequency in that population
    if(any(t_int-t_b-t_s < t_sampled[intro_sel_pops])){
      # need to modify Ry_selBG and Ry_nonSelBG for these populations if at least one
      # of them split off during the sweep (they have a different final freq so a 
      # different probability of staying linked)
      
      change_pops = intro_sel_pops[t_int-t_b-t_s < t_sampled[intro_sel_pops]]
      
      Ry_selBG_midSweep_fxn = Vectorize(function(pop){
        new_ts = t_int - t_b - t_sampled[pop]
        new_x = g*exp(s*new_ts) / (1-g + g*exp(s*new_ts))
        return((new_x*(1-g)/g)^(-r/s))
      },USE.NAMES = F)
      
      Ry_nonSelBG_midSweep_fxn = Vectorize(function(pop){
        new_ts = t_int - t_b - t_sampled[pop]
        new_x = g*exp(s*new_ts) / (1-g + g*exp(s*new_ts))
        return((1/(1-new_x))^(-r/s))
      },USE.NAMES = F)
      
      Ry_selBG[change_pops] = Ry_selBG_midSweep_fxn(change_pops)
      Ry_nonSelBG[change_pops] = Ry_nonSelBG_midSweep_fxn(change_pops)
      
    }
    
    # any function called here should use the environment of the parent frame
    environment(Rx)=environment()
    environment(p_events)=environment()
    environment(P.g)=environment()
    environment(p_nr_wholeTime)=environment()
    
    # probability of never recombining out of background during neutral phase II
    .Rx=sapply(intro_sel_pops,Rx) ; names(.Rx)=intro_sel_pops
    # probability of never recombining out of background during entire time from sampling to introgression
    .p_nr_wholeTime=sapply(intro_pops,p_nr_wholeTime) ; names(.p_nr_wholeTime)=intro_pops
    
    # all of the below matrices are symmetric
    .P=outer(intro_sel_pops,intro_sel_pops,Vectorize( function(pop1,pop2){
      # t_shared_g = the amount of time the selected variant is at freq g in the ancestral population of pop1 and pop2
      # t_nonShared_g = the amount of time the selected variant is at freq g when pop1 and pop2 are isolated/diverged from each other
      
      # If sweep starts in ancestor of pair of populations (lineages are in same population for whole time selected variant is at freq g):
      t_shared_g=ifelse((t_anc[pop1,pop2]-t_b)>0, t_b, t_anc[pop1,pop2])
      t_nonShared_g=t_b-t_shared_g # t_nonShared_g=0 if sweep starts in ancestor of pair of populations
      
      # nr_selBG = probability of never recombining out of selected variant background between
      #          sweep finish and time lineages are in same population
      #          this involves two time periods: during sweep & during time t_nonShared_g
      
      nr_selBG_pop1 = Ry_selBG[pop1]*exp(-r*t_nonShared_g)
      nr_selBG_pop2 = Ry_selBG[pop2]*exp(-r*t_nonShared_g)
      
      # this is the final term to return
      nr_selBG_pop1*nr_selBG_pop2*(P.g(pop1,pop2,t_shared_g)) +
        (nr_selBG_pop1*(1-nr_selBG_pop2) + (1-nr_selBG_pop1)*nr_selBG_pop2)*(Rz*g*F_estimate[source,source] + (1-Rz)*F_estimate[pop1,pop2]) +
        (1-nr_selBG_pop1)*(1-nr_selBG_pop2)*F_estimate[pop1,pop2]
      
    }, USE.NAMES=F ))
    .Q=outer(intro_sel_pops,intro_sel_pops,Vectorize( function(pop1,pop2){
      # need to average over the two possibilities of which lineage recombined out and which stayed linked
      
      # when pop1 stays linked
      if_pop1_linked = Ry_selBG[pop1]*Rz*Ry_nonSelBG[pop2]*Rz*0 + #neither rec out of bg
        (1-Ry_selBG[pop1]*Rz)*(1-g)*Ry_nonSelBG[pop2]*Rz*.both_MH[pop1,pop2] + #linked recs out and know it's type MH, unlinked does not
        Ry_selBG[pop1]*Rz*(1-Ry_nonSelBG[pop2]*Rz)*g*F_estimate[source,source] + # linked does not rec out (type source), unlinked does and know it's type source
        (1-Ry_selBG[pop1]*Rz)*(1-Ry_nonSelBG[pop2]*Rz)*F_estimate[pop1,pop2] #both rec out
      # when pop2 stays linked
      if_pop2_linked = Ry_selBG[pop2]*Rz*Ry_nonSelBG[pop1]*Rz*0 + #neither rec out of bg
        (1-Ry_selBG[pop2]*Rz)*(1-g)*Ry_nonSelBG[pop1]*Rz*.both_MH[pop1,pop2] + #linked recs out and know it's type MH, unlinked does not
        Ry_selBG[pop2]*Rz*(1-Ry_nonSelBG[pop1]*Rz)*g*F_estimate[source,source] + # linked does not rec out (type source), unlinked does and know it's type source
        (1-Ry_selBG[pop2]*Rz)*(1-Ry_nonSelBG[pop1]*Rz)*F_estimate[pop1,pop2] #both rec out
      
      return(mean(c(if_pop1_linked,if_pop2_linked)))
      
    }, USE.NAMES=F ))
    .U=outer(intro_sel_pops,intro_sel_pops,Vectorize(function(pop1,pop2){
      Ry_nonSelBG[pop1]*Ry_nonSelBG[pop2]*Rz^2*.both_MH[pop1,pop2] + (1-Ry_nonSelBG[pop1]*Rz)*(1-Ry_nonSelBG[pop2]*Rz)*F_estimate[pop1,pop2] + 
        (Ry_nonSelBG[pop1]*Rz*(1-Ry_nonSelBG[pop2]*Rz) + Ry_nonSelBG[pop2]*Rz*(1-Ry_nonSelBG[pop1]*Rz) )*(1-g)*.both_MH[pop1,pop2]
    }, USE.NAMES=F ))
    .p_coal_anyBG=x^2*.P + (1-x)^2*.U + 2*x*(1-x)*.Q
    
    
    # probability selPop_thenAgainst_aa pop doesn't recombine out of its non-Neanderthal background
    # during the time the lineage is at frequency x within its own population
    # Note: we don't need to specify these times between populations, because that will be
    # calculated for each population pair
    # Also, these values can't be negative if this category of populations is specified in the first place
    Rx_selAgainst=sapply(selPops_thenAgainst_aa, function(pop){
      p_no_rec=exp(-r*(anc_times[nest_index[pop]] - t_b - t_s))
      if(p_no_rec>1) p_no_rec=1 # if there actually isn't any time that the lineage is at freq x in a population (happens when sweep starts in ancestral pop but finishes after a population diverged)
      return(p_no_rec)
    })
    
    # probability selPop_thenAgainst_aa lineage does not recombine out of 
    # non-Neanderthal background before it's at frequency x in these populations
    R_nonX_selAgainst=sapply(selPops_thenAgainst_aa, function(pop){
      p_no_rec=exp(-r*(t_int-anc_times[nest_index[pop]]-t_sampled[pop]))
      return(p_no_rec)
    },USE.NAMES = F)
    
    
    # start getting probabilities of coalescing under selection
    F.S_noMig = F_estimate # define a matrix of these probabilities under no migration
    
    ### Relationships between Neanderthal-admixed populations and source (Neanderthals), neverSel (Africans)
    
    # probability of definitely being linked to Neanderthal lineage at time of the sweep completion
    neaBG_selPops_AA_fxn = Vectorize(function(pop){
      (.Rx[pop] + (1-.Rx[pop])*x)
    },USE.NAMES = F)
    neaBG_selPops_AA = neaBG_selPops_AA_fxn(c(selPops_anc_AA,selPops_ind_AA))
    
    neaBG_selPops_aa_fxn = Vectorize(function(pop){
      (1-.Rx[pop])*x
    },USE.NAMES = F)
    neaBG_selPops_aa = neaBG_selPops_aa_fxn(c(selPops_anc_aa,selPops_ind_aa))
    
    neaBG_selPops_thenAgainst_aa_fxn = Vectorize(function(pop){
      (1-R_nonX_selAgainst[pop])*Rx_selAgainst[pop]*fin_freqs[gsub('_aa','',pop)] + (1-Rx_selAgainst[pop])*x
    },USE.NAMES = F)
    neaBG_selPops_thenAgainst_aa = neaBG_selPops_thenAgainst_aa_fxn(selPops_thenAgainst_aa)
    
    # combine info from all pops
    neaBG = unlist(c(neaBG_selPops_AA,neaBG_selPops_aa,neaBG_selPops_thenAgainst_aa))  #unlist is necessary if there are no specified pops in a category (AA or aa)
    neaBG = neaBG[intro_sel_pops]
    
    # probability of definitely being linked to MH (modern human / non-Neanderthal) lineage at time of the sweep completion
    mhBG_selPops_AA_fxn = Vectorize(function(pop){
      (1-.Rx[pop])*(1-x)
    },USE.NAMES = F)
    mhBG_selPops_AA = mhBG_selPops_AA_fxn(c(selPops_anc_AA,selPops_ind_AA))
    
    mhBG_selPops_aa_fxn = Vectorize(function(pop){
      .Rx[pop] + (1-.Rx[pop])*(1-x)
    },USE.NAMES = F)
    mhBG_selPops_aa = mhBG_selPops_aa_fxn(c(selPops_anc_aa,selPops_ind_aa))
    
    mhBG_selPops_thenAgainst_aa_fxn = Vectorize(function(pop){
      # the following is split up into probabilities of being on mh or nea background before joining other pops
      (R_nonX_selAgainst[pop] + (1-R_nonX_selAgainst[pop])*(1-fin_freqs[gsub('_aa','',pop)]))*(Rx_selAgainst[pop] + (1-Rx_selAgainst[pop])*(1-x)) + 
        (1-R_nonX_selAgainst[pop])*fin_freqs[gsub('_aa','',pop)]*(1-Rx_selAgainst[pop])*(1-x)
    },USE.NAMES = F)
    mhBG_selPops_thenAgainst_aa = mhBG_selPops_thenAgainst_aa_fxn(selPops_thenAgainst_aa)
    
    # combine info from all pops
    mhBG = unlist(c(mhBG_selPops_AA,mhBG_selPops_aa,mhBG_selPops_thenAgainst_aa))  #unlist is necessary if there are no specified pops in a category (AA or aa)
    mhBG = mhBG[intro_sel_pops]
    
    # sanity check that it sums to 1: 
    # all(mhBG + neaBG == 1)
    
    # btwn selected, Neanderthal-admixed pops and source (Neanderthal) pop
    F.S_noMig[source,intro_sel_pops] = F.S_noMig[intro_sel_pops,source] = neaBG*Ry_selBG*Rz*F_estimate[source,source] + (mhBG*(1-Ry_nonSelBG*Rz) + neaBG*(1-Ry_selBG*Rz))*F_estimate[source,intro_sel_pops]
    # btwn selected, Neanderthal-admixed pops and non-Neanderthal-admixed pop
    F.S_noMig[neverSel,intro_sel_pops] = F.S_noMig[intro_sel_pops,neverSel] = (mhBG*Ry_nonSelBG*Rz/(1-g) + mhBG*(1-Ry_nonSelBG*Rz) + neaBG*(1-Ry_selBG*Rz))*F_estimate[neverSel,intro_sel_pops]
    
    # btwn non-selected, Neanderthal-admixed pops and source (Neanderthal) pop
    F.S_noMig[source,noSel_AA] = F.S_noMig[noSel_AA,source] = .p_nr_wholeTime[noSel_AA]*F_estimate[source,source] + (1-.p_nr_wholeTime[noSel_AA])*F_estimate[source,noSel_AA]
    F.S_noMig[source,noSel_aa] = F.S_noMig[noSel_aa,source] = (1-.p_nr_wholeTime[noSel_aa])*F_estimate[source,noSel_aa]
    
    # btwn non-selected, Neanderthal-admixed pops and and non-Neanderthal-admixed pop
    F.S_noMig[neverSel,noSel_AA] = F.S_noMig[noSel_AA,neverSel] = (1-.p_nr_wholeTime[noSel_AA])*F_estimate[neverSel,noSel_AA]
    F.S_noMig[neverSel,noSel_aa] = F.S_noMig[noSel_aa,neverSel] = (.p_nr_wholeTime[noSel_aa]/(1-g) + (1-.p_nr_wholeTime[noSel_aa]))*F_estimate[neverSel,noSel_aa]
    
    
    ### Relationships between Neanderthal-admixed selected and non-selected populations
    
    # rows: sel pops, cols: non-sel pops
    # I want to multiply each row of "specific" by "general" to get a matrix with dimensions  number of sel pops x number of non-sel pops 
    # "specific" refers to terms specific to a given non-selected population
    
    # defining two cases:
    # case1: probability lineage sampled from non-selected (noSel) pop coalesces
    #        lineage sampled from selected pop, when noSel lineage does *not*
    #        recombine out of its original background 
    # case2: probability lineage sampled from non-selected (noSel) pop coalesces
    #        lineage sampled from selected pop, when noSel lineage *does*
    #        recombine out of its original background 
    # case 1 and case 2 are defined separately for AA and aa partitions because
    # lineages sampled from each begin on different ancestry backgrounds
    
    general2 = neaBG*Ry_selBG*Rz*g*F_estimate[source,source] # this line applies to both AA and aa partitions
    
    if(length(noSel_aa)>0){
      # case 1
      # noSel_aa lineage doesn't recombine out of MH background: lineages can coalesce if they both did not descend from Neanderthals
      specific1_aa = .p_nr_wholeTime[noSel_aa] * .both_MH[noSel_aa,intro_sel_pops]
      general1_aa = mhBG*(Ry_nonSelBG*Rz + (1-Ry_nonSelBG*Rz)*(1-g)) + neaBG*(1-Ry_selBG*Rz)*(1-g)
      case1_aa = t(specific1_aa)*general1_aa 
      
      # case 2
      # noSel_aa lineage does recombine out of MH background: now we condition on what happens
      # in lineage sampled from selected population to get probability that they coalesce
      specific2a = (1-.p_nr_wholeTime[noSel_aa])
      specific2b = t(F_estimate[noSel_aa,intro_sel_pops])*(mhBG*(1-Ry_nonSelBG*Rz) + neaBG*(1-Ry_selBG*Rz))
      specific2c = t(.both_MH[noSel_aa,intro_sel_pops])*mhBG*Ry_nonSelBG*Rz*(1-g)

      if(length(noSel_aa)>1){
        case2_aa = ( specific2b + specific2c + general2 ) %*% diag(specific2a)
      } else {
        case2_aa =  specific2a * ( specific2b + specific2c + general2 )
      }
      
      F.S_selPops_noSelaa = case1_aa + case2_aa
      
      # assignment
      F.S_noMig[intro_sel_pops,noSel_aa] = F.S_selPops_noSelaa
      F.S_noMig[noSel_aa,intro_sel_pops] = t(F.S_selPops_noSelaa)
      
    }
    
    if(length(noSel_AA)>0){
      # noSel_AA lineage doesn't recombine out of Neanderthal background: lineages can coalesce if they both descended from Neanderthals
      general1_AA = F_estimate[source,source]*(mhBG*(1-Ry_nonSelBG*Rz)*g + neaBG*(Ry_selBG*Rz + (1-Ry_selBG*Rz)*g))
      case1_AA = general1_AA %*% t(.p_nr_wholeTime[noSel_AA])
      
      # case 2
      # noSel_AA lineage does recombine out of Neanderthal background: now we condition on what happens
      # in lineage sampled from selected population to get probability that they coalesce
      specific2a = (1-.p_nr_wholeTime[noSel_AA])
      specific2b = t(F_estimate[noSel_AA,intro_sel_pops])*(mhBG*(1-Ry_nonSelBG*Rz) + neaBG*(1-Ry_selBG*Rz))
      specific2c = t(.both_MH[noSel_AA,intro_sel_pops])*mhBG*Ry_nonSelBG*Rz*(1-g)
      
      if(length(noSel_AA)>1){
        case2_AA = ( specific2b + specific2c + general2 ) %*% diag(specific2a)
      } else {
        case2_AA =  t(specific2a * ( specific2b + specific2c + general2 ))
      }
      
      F.S_selPops_noSelAA = (case1_AA + case2_AA)
      
      # assignment
      F.S_noMig[intro_sel_pops,noSel_AA] = F.S_selPops_noSelAA
      F.S_noMig[noSel_AA,intro_sel_pops] = t(F.S_selPops_noSelAA)
      
      
    }
    
    ### Relationships within non-selected populations
    
    # both lineages sampled from noSel_aa populations
    never_rec_aa = matrix(.p_nr_wholeTime[noSel_aa],ncol=1) 
    F.S_noMig[noSel_aa,noSel_aa] = (never_rec_aa %*% t(never_rec_aa)) * .both_MH[noSel_aa,noSel_aa] +
      ((1-never_rec_aa) %*% t(1-never_rec_aa)) * F_estimate[noSel_aa,noSel_aa] +
      (never_rec_aa %*% t(1-never_rec_aa) + (1-never_rec_aa) %*% t(never_rec_aa) )*(1-g)*.both_MH[noSel_aa,noSel_aa]
    
    # both lineages sampled from noSel_AA populations
    never_rec_AA = matrix(.p_nr_wholeTime[noSel_AA],ncol=1) 
    F.S_noMig[noSel_AA,noSel_AA] = (never_rec_AA %*% t(never_rec_AA)) * F_estimate[source,source] +
      ((1-never_rec_AA) %*% t(1-never_rec_AA)) * F_estimate[noSel_AA,noSel_AA] +
      (never_rec_AA %*% t(1-never_rec_AA) + (1-never_rec_AA) %*% t(never_rec_AA) )*g*F_estimate[source,source]
    
    
    # one lineage sampled from noSel_aa (rows) and the other lineage sampled from noSel_AA (columns)
    
    # neither recombine out
    noSel_btwn = (never_rec_aa %*% t(never_rec_AA) * 0) + 
      # aa don't recombine out, AA do
      ((never_rec_aa %*% t(1-never_rec_AA)) *(1-g)*.both_MH[noSel_aa,noSel_AA] ) +
      # aa recombine out, AA don't
      ((never_rec_aa %*% t(1-never_rec_AA) ) *g*F_estimate[source,source]  ) +
      # both recombine out
      ((1-never_rec_aa) %*% t(1-never_rec_AA) ) * F_estimate[noSel_aa,noSel_AA]
    
    F.S_noMig[noSel_aa,noSel_AA] = noSel_btwn
    F.S_noMig[noSel_AA,noSel_aa] = t(noSel_btwn)
    
    
    ### Relationships within selected populations
    
    # t_shared_x represents the amount of time the selected lineage is at frequency x 
    # when populations i & j share a common ancestor -- this is the amount of time
    # between the sweep finish and their divergence
    t_shared_x=t_int-t_b-t_s-t_div[intro_sel_pops,intro_sel_pops]
    # remove t_sampled from t_shared_x for all part pops that actually belong to the same pop
    same_pops=which(t_div[intro_sel_pops,intro_sel_pops]==0,arr.ind = T)
    same_pops=cbind(rownames(same_pops),same_pops)
    dimnames(same_pops)[[2]][1] = 'ref_pop'
    # change t_shared_x among populations that are actually partitions of the same population
    # to account for sampling time (t_div matrix is measured from the present)
    t_shared_x[t_div[intro_sel_pops,intro_sel_pops]==0] = apply(same_pops,1,function(index_set)  t_shared_x[as.integer(index_set['row']),as.integer(index_set['col'])] - t_sampled[index_set['ref_pop']]  )
    # change t_shared_x for within selPops_thenAgainst_aa, as they'll only have the lineage
    # at frequency x when they merge back with selected populations
    # anc_times[nest_index[pop]] specifies the amount of shared time that population has with the populations it diverged from when it split off
    diag(t_shared_x)[selPops_thenAgainst_aa] = unlist(sapply(selPops_thenAgainst_aa,function(pop) anc_times[nest_index[pop]]-t_b-t_s ))
    # any negative values should be made zero (meaning the sweep  didn't complete in their shared ancestor)
    t_shared_x[t_shared_x<0]=0
    
    # empty matrix to fill in relationships among selected pops
    start_matrix=matrix(NA,nrow=n_intro_sel_pops,ncol=n_intro_sel_pops,dimnames=list(intro_sel_pops,intro_sel_pops))
    
    # pairs of populations whose relationships we will fill into the matrix
    pop_pair_matrix = combn(intro_sel_pops,m=2) 
    # order pairs to correspond to upper triangle of matrix when written as a vector, reading row by row (excludes diagonal)
    pairs = split(pop_pair_matrix, rep(1:ncol(pop_pair_matrix), each = nrow(pop_pair_matrix))) 
    
    ## PART 1: get values corresponding to the competing Poisson processes of
    #          coalescence and recombination during neutral phase II
    #          (recall that lineages can coalesce due to association with the 
    #           same background, when they're in the same population, during neutral phase II)
    
    # get a list in which each element contains results for each population pair
    # these results are a named vector of the probability of an event occurring, the probability that the first event is coalescence, and the probability that the first event is recombination
    
    # we will get a different matrix to represent each of these values for each population pair
    # population pairs are in the order in which you fill the lower triangle of a matrix
    # subsetting with lower.tri() gives vector of values column by column e.g. c(col1[2:n],col2[3:n],col3[4:n],...) 
    
    # WHEN BOTH LINEAGES ARE LINKED TO SEL ALLELE AT THE SWEEP FINISH
    # get values for upper triangle of matrix for probabilities of coalescing/recombining/either of these events occurring
    p_events_b4_SF_selBG_upper = lapply(pairs, function(pop_pair) p_events(pop1=pop_pair[1],pop2=pop_pair[2],time=t_shared_x[pop_pair[1],pop_pair[2]],bg_freq=x) )
    # get within population values
    p_events_b4_SF_selBG_diag = lapply(intro_sel_pops, function(pop) p_events(pop1=pop,pop2=pop,time=t_shared_x[pop,pop],bg_freq=x) )
    
    w_selBG = start_matrix
    w_selBG[lower.tri(w_selBG)] = sapply(p_events_b4_SF_selBG_upper, function(pair_vec) pair_vec['w'] )
    diag(w_selBG) = sapply(p_events_b4_SF_selBG_diag, function(pair_vec) pair_vec['w'] )
    w_selBG[upper.tri(w_selBG)] = t(w_selBG)[upper.tri(w_selBG)]
    
    p_coal_selBG = start_matrix
    p_coal_selBG[lower.tri(p_coal_selBG)] = sapply(p_events_b4_SF_selBG_upper, function(pair_vec) pair_vec['p_coal'] )
    diag(p_coal_selBG) = sapply(p_events_b4_SF_selBG_diag, function(pair_vec) pair_vec['p_coal'] )
    p_coal_selBG[upper.tri(p_coal_selBG)] = t(p_coal_selBG)[upper.tri(p_coal_selBG)]
    
    p_rec_selBG = start_matrix
    p_rec_selBG[lower.tri(p_rec_selBG)] = sapply(p_events_b4_SF_selBG_upper, function(pair_vec) pair_vec['p_rec'] )
    diag(p_rec_selBG) = sapply(p_events_b4_SF_selBG_diag, function(pair_vec) pair_vec['p_rec'] )
    p_rec_selBG[upper.tri(p_rec_selBG)] = t(p_rec_selBG)[upper.tri(p_rec_selBG)]
    
    
    # WHEN BOTH lineageS ARE LINKED TO NON SEL lineage
    # get values for upper triangle
    p_events_b4_SF_nonSelBG_upper = lapply(pairs, function(pop_pair) p_events(pop1=pop_pair[1],pop2=pop_pair[2],time=t_shared_x[pop_pair[1],pop_pair[2]],bg_freq=1-x) )
    # get within population values
    p_events_b4_SF_nonSelBG_diag = lapply(intro_sel_pops, function(pop) p_events(pop1=pop,pop2=pop,time=t_shared_x[pop,pop],bg_freq=1-x) )
    
    w_nonSelBG = start_matrix
    w_nonSelBG[lower.tri(w_nonSelBG)] = sapply(p_events_b4_SF_nonSelBG_upper, function(pair_vec) pair_vec['w'] )
    diag(w_nonSelBG) = sapply(p_events_b4_SF_nonSelBG_diag, function(pair_vec) pair_vec['w'] )
    w_nonSelBG[upper.tri(w_nonSelBG)] = t(w_nonSelBG)[upper.tri(w_nonSelBG)]
    
    # if x=1, then freq of non sel bg = 0 --> get p_coal of 1/(2*N*0) = NaN
    # since there's no change lineages will be linked to nonSel BG in this case, 
    # anything related to nonSel BG won't contribute, but to avoid NaNs we set 
    # the prob of an event equal to 0
    # note that p_coal vs p_rec are totally wrong, but once again those won't be evaluated
    if(x==1) w_nonSelBG[]=0 
    
    p_coal_nonSelBG = start_matrix
    p_coal_nonSelBG[lower.tri(p_coal_nonSelBG)] = sapply(p_events_b4_SF_nonSelBG_upper, function(pair_vec) pair_vec['p_coal'] )
    diag(p_coal_nonSelBG) = sapply(p_events_b4_SF_nonSelBG_diag, function(pair_vec) pair_vec['p_coal'] )
    p_coal_nonSelBG[upper.tri(p_coal_nonSelBG)] = t(p_coal_nonSelBG)[upper.tri(p_coal_nonSelBG)]
    
    p_rec_nonSelBG = start_matrix
    p_rec_nonSelBG[lower.tri(p_rec_nonSelBG)] = sapply(p_events_b4_SF_nonSelBG_upper, function(pair_vec) pair_vec['p_rec'] )
    diag(p_rec_nonSelBG) = sapply(p_events_b4_SF_nonSelBG_diag, function(pair_vec) pair_vec['p_rec'] )
    p_rec_nonSelBG[upper.tri(p_rec_nonSelBG)] = t(p_rec_nonSelBG)[upper.tri(p_rec_nonSelBG)]
    
    ## PART 2: get probabilities of lineages being on a certain background
    #          with each transition between phases of the selected allele frequency
    #          trajectory and transitions between when lineages belong to a
    #          shared ancestral population
    
    # probability of not recombining out of background while lineages are in the same pop, but before the sweep has finished
    Rx2 = exp(-r * t_shared_x)
    
    # probability lineage sampled in ROW POP does not recombine out of background while in different pop from lineage sampled in COLUMN POP
    # row: population of interest, column: relationship with other pop
    # note: this matrix must be filled out
    # following is the same as: Rx1 = t(sapply(selPops_AA,function(pop) .Rx[pop] / Rx2[pop,] ))
    Rx1 = diag(.Rx) %*% (1/Rx2) 
    dimnames(Rx1)=list(intro_sel_pops,intro_sel_pops)
    
    # define transpose of Rx1. When you multiply it with Rx1 (element by element), each element
    # equals the probability that both lineages didn't recombine out by the time they're in the same population
    t_Rx1 = t(Rx1)
    
    # for all populations, can describe probability of lineage being on sel vs non sel background by time
    # it is in the same population as the other lineage
    # following matrix describes probability lineage from ROW POP is linked to sel lineage at time
    # it begins to be in the same population as COL POP (looking pastwards)
    
    # make separate matrices for different categories of populations
    # then we can bind those together
    
    linked_selBG_selPops_AA = t(sapply(c(selPops_anc_AA,selPops_ind_AA),function(pop){
      Rx1[pop,] + (1-Rx1[pop,])*x
    }))
    
    if(selPops_aa_exist){
      linked_selBG_selPops_aa = t(sapply(c(selPops_anc_aa,selPops_ind_aa),function(pop){
        (1-Rx1[pop,])*x
      }))
    } else{
      linked_selBG_selPops_aa = matrix(nrow=0,ncol=n_intro_sel_pops,dimnames=list(NULL,intro_sel_pops))
    }
    
    
    if(length(selPops_thenAgainst_aa)>0){
      linked_selBG_selPops_thenAgainst_aa = t(sapply(selPops_thenAgainst_aa,function(pop){
        # how much time is the lineage at freq x before the populations join?
        # how much time is the lineage at freq non-x (i.e. when removed by sel/drift) before the populations join?
        
        #R_nonX_selAgainst covers the amount of time the lineage is at low freq after selection
        #--by definition this population has no shared history with others during this phase
        
        # similarly, t_shared_x describes amount of time they share at freq x
        # perhaps lineage is at freq x in selAgainst pop before they join, though, so
        # need to account for that amount of time
        
        # amount of time lineage is at freq x in focal pop
        total_x_pop=ifelse(anc_times[nest_index[pop]] - t_b >= t_s, anc_times[nest_index[pop]] - t_b - t_s, 0)
        # amount of time lineage is at freq x when not in common ancestor with other pops
        t_nonShared_x = total_x_pop - t_shared_x[pop,] # when total_x_pop=0, t_shared_x = 0 
        # prob of not recombining out in above time
        p_noRec_nonShared_x = exp(-r*t_nonShared_x)
        
        # first term below doesn't need conditioning on what happens during nonX phase, 
        # because any of the scenarios can still lead to recombination during t_nonShared_x
        (1-p_noRec_nonShared_x)*x + (1-R_nonX_selAgainst[pop])*fin_freqs[gsub('_aa','',pop)]*p_noRec_nonShared_x
        
      }))
      
      # results from all populations
      linked_selBG = rbind(linked_selBG_selPops_AA,linked_selBG_selPops_aa,linked_selBG_selPops_thenAgainst_aa)
      
    } else{
      # results from all populations (and there's no selPops_thenAgainst_aa to include)
      linked_selBG = rbind(linked_selBG_selPops_AA,linked_selBG_selPops_aa)
    }
  
    # final product: probability that lineage from ROW POP is linked to
    # selected or non-selected background at the time it begins to be in the same
    # population as COL POP (looking pastwards)
    linked_selBG = linked_selBG[intro_sel_pops,intro_sel_pops]
    linked_nonSelBG = 1-linked_selBG
    
    # probability that both, none, or one of the lineages are linked to the selected lineage by
    # the time the two lineages are segregating in the same ancestral population
    both_selBG = linked_selBG*t(linked_selBG)
    both_nonSelBG = linked_nonSelBG*t(linked_nonSelBG)
    dif_BG = linked_selBG*t(linked_nonSelBG) + linked_nonSelBG*t(linked_selBG)
    
    ## PART 3: probabilities of coalescing between each pair of selected populations
    #          conditioned on whether by the time the lineages are in the same 
    #          population, that they're both on the selected, non-selected, or 
    #          different backgrounds
    
    
    # probability of coalescing if both lineages are on sel background when they're in the same population
    if_both_selBG=(1-w_selBG)*.P + w_selBG*(p_coal_selBG + p_rec_selBG*(Rx2*(x*.P + (1-x)*.Q) + (1-Rx2)*(.p_coal_anyBG)))
    
    # probability of coalescing if both lineages are on non-sel background when they're in the same population
    if_both_nonSelBG=(1-w_nonSelBG)*.U + w_nonSelBG*(p_coal_nonSelBG + p_rec_nonSelBG*(Rx2*(x*.Q + (1-x)*.U) + (1-Rx2)*.p_coal_anyBG ))
    
    # probability of coalescing if lineages are on different backgrounds when they're in the same population
    if_difBG=Rx2^2*.Q + Rx2*(1-Rx2)*(x*(.P+.Q)+(1-x)*(.U+.Q)) +(1-Rx2)^2*.p_coal_anyBG
    
    F.S_noMig[intro_sel_pops,intro_sel_pops] = both_selBG*if_both_selBG + both_nonSelBG*if_both_nonSelBG + dif_BG*if_difBG
    
    # Make final Fsel matrix
    selMatrix = F.S_noMig
    
    # incorporate migration if Europe_AA exists (we ask whether this is the case in make_admixture_objects.R)
    # in the future, when there are more populations that receive migrants, we can make this a vectorized function
    if(migration){
      i = paste0(Eur,'_AA') # recipient population
      donors=names(which(migs[i,]!=0))
      
      # then calculate F sel matrix elements accordingly
      aa = gsub("AA","aa",i)
      if(!(aa %in% pops)){aa=i} # if there isn't an aa version, just make Europe like itself
      t_mig=mean(mig_times[i,],na.rm=T)
      # probability does not recombine out of current pop before migration into it
      u=exp(-r*t_mig*(1-x)) 
      
      # relationships between pop i and given pop j (vector elements represent each pop j)
      other_pops = setdiff(pops,i)
      
      # this function is called in below function, and later for within pop values
      # gets probability of coalescing with pop j weighted by probability of migration to donor d
      coal_donor=Vectorize(function(d,J){
        migs[i,d]*F.S_noMig[J,d] #prob of migrating to other pop * prob coalesce btwn j & that donor pop
      },vectorize.args = "d")
      
      # probability of coalescence with pop j, weighted by migration possibilities
      coal_pop=Vectorize(function(j){
        unlist(u*sum(coal_donor(donors,J=j)) + (1-u)*(F.S_noMig[aa,j]))
      })
      
      selMatrix[i,other_pops] = selMatrix[other_pops,i] = coal_pop(other_pops)
      
      # NOW DO WITHIN POPS
      # use permutations() to generate all possibilities for where each lineage originates from
      # permutations(n,r) is a fxn that generates all permutations of n elements taken r at a time
      # n is size of options (one number input)
      # repeats are allowed because it's possible to sample 2 lineages from the same population
      # output: each row specifies a different permutation, but just gives index of donors to represent population
      # we sample 2 lineages, each of which could come from any of the donor/recipient populations
      options=permutations(n=length(donors),r=2,repeats.allowed = T)
      # get probability of coalescing for each permutation, and add them together
      prob_coal = Vectorize(function(combo){
        pop1=donors[options[combo, 1]]
        pop2=donors[options[combo, 2]]
        if(pop1==i & pop2==i){
          migs[i,pop1]*migs[i,pop2]*F.S_noMig[i,i]
        } else if(any(c(pop1,pop2)==i)){ # one population is i, the other is not
          migs[i,pop1]*migs[i,pop2]*F.S_noMig[pop1,pop2]
        } else{
          return(migs[i,pop1]*migs[i,pop2]*F.S_noMig[pop1,pop2])
        }
      })
      # get prob coalescing between lineage that migrates and lineage that doesnt
      # (this includes probability of lineage migrating to a certain donor pop in the first place)
      
      both_linked=sum(prob_coal(1:nrow(options)))
      selMatrix[i,i]=unlist(u^2*both_linked + (1-u)^2*F.S_noMig[aa,aa] +
                              2*u*(1-u)*(sum(coal_donor(donors,J=aa))))
    }
    
    
    ### now modify the _Aa subpops of populations that have them
    Aa_pops = grep("_Aa",pops,value=T)
    
    if(length(Aa_pops)>0){
      # function to get between pop values, which is called within the next function
      Aa_btwn = Vectorize(function(m,AA_pop,aa_pop){
        # Args: m = population that isn't Aa
        if(m %in% Aa_pops){
          m_AA = gsub("_Aa","_AA",m)
          m_aa = gsub("_Aa","_aa",m)
          return(mean(c(selMatrix[AA_pop,m_AA], selMatrix[AA_pop,m_aa], selMatrix[aa_pop,m_AA], selMatrix[aa_pop,m_aa] ) ))
        } else{
          return(0.5*selMatrix[AA_pop,m] + 0.5*selMatrix[aa_pop,m])
        }
      },vectorize.args = "m")
      
      # function to get values for within/btwn given Aa pop and all other pops that aren't Aa pops
      Aa_values = Vectorize(function(Aa){
        AA = gsub("_Aa","_AA",Aa) # selected subpop
        aa = gsub("_Aa","_aa",Aa) # neutral subpop
        
        btwn_values = Aa_btwn(setdiff(pops,Aa),AA_pop=AA,aa_pop=aa)
        within_value =  0.25*selMatrix[aa,aa] + 0.25*selMatrix[AA,AA] + 0.5*selMatrix[AA,aa]
        names(within_value)=Aa
        
        return(c(btwn_values,within_value)[pops])
      })
      
      Aa_relationships=Aa_values(Aa_pops)
      
      selMatrix[pops,Aa_pops] = Aa_relationships
      selMatrix[Aa_pops,pops] = t(Aa_relationships)
    }
    
  
    return(selMatrix[real_pops,real_pops])
  }
  
  # now that we've defined the function to get F^(S) for all recombination rates
  # under the combination of parameters defined above, we can define those matrices:
  
  F.S_all_rec = purrr::map(mid_genetic_distances, calc_F.S)
  return(F.S_all_rec)
}


### NULL MODEL THAT ACCOUNTS FOR PARTITIONING ###
calc_F.N = function(r,F_estimate) {
  # Computes F^(N) for a single genetic distance bin
  #
  # Args:
  # r: recombination rate (genetic distance) between neutral and selected site
  #
  # Returns:
  #	F matrix under the null model, F^(N), for a single genetic distance bin
  
  #### step 1: define a matrix of within & between pop coancestry for _AA and _aa pops
  
  #non-African modern humans,non of whom experience selection in the null model
  mh_null_AA=pops[ (grepl('_AA',pops)) & !(pops %in% c(source,neverSel) )]
  mh_null_aa=pops[  !(grepl('_Aa',pops)) & !(pops %in% c(source,neverSel,mh_null_AA) )]
  
  F_null=F_estimate # start with baseline of neutral F, from which to modify
  
  p_r=exp(-r*t_int) # probability of not recombining out of background from time sampled to time of introgression
  
  for(i in mh_null_AA) {
    
    # within
    w_ii=1-exp(-(t_int-t_sampled[i])*((1/(2*Ne[i]*g))+2*r))
    term1= (1/(1+4*Ne[i]*g*r)) + (4*Ne[i]*g*r/(1+4*Ne[i]*g*r))*((1-exp(-r*(t_int-t_sampled[i])))*F_estimate[i,i] + exp(-r*(t_int-t_sampled[i]))*g*F_estimate[source,source])
    F_null[i,i]=w_ii*term1 + (1-w_ii)*F_estimate[source,source]
    
    #btwn mh_null_AA pops
    for(j in mh_null_AA){
      if(i!=j){
        w_ij=1-exp(-t_anc[i,j]*((1/(2*Ne[i]*g))+2*r))
        term1= (1/(1+4*mean(Ne[i],Ne[j])*g*r)) + (4*mean(Ne[i],Ne[j])*g*r/(1+4*mean(Ne[i],Ne[j])*g*r))*((1-exp(-r*t_anc[i,j]))*F_estimate[i,j] + exp(-r*t_anc[i,j])*g*F_estimate[source,source])
        F_null[i,j]=F_null[j,i]=exp(-r*(t_div[i,j]-t_sampled[i]))*exp(-r*(t_div[i,j]-t_sampled[j]))*(w_ii*term1 +(1-w_ii)*F_estimate[source,source]) +
          (1-exp(-r*(t_div[i,j]-t_sampled[i])))*(1-exp(-r*(t_div[i,j]-t_sampled[j])))*F_estimate[i,j] +
          (exp(-r*(t_div[i,j]-t_sampled[i]))*(1-exp(-r*(t_div[i,j]-t_sampled[j]))) + (1-exp(-r*(t_div[i,j]-t_sampled[i])))*exp(-r*(t_div[i,j]-t_sampled[j])))*(exp(-r*t_anc[i,j])*g*F_estimate[source,source] + (1-exp(-r*t_anc[i,j]))*F_estimate[i,j])
      }
    }
    
    #btwn mh_null_AA and mh_null_aa
    for (l in mh_null_aa) {
      
      both_MH=(F_estimate[i,l]-g^2*F_estimate[source,source])/(1-g)^2 # prob coalesce given two lineages are type MH
      F_null[i,l]=F_null[l,i]= (1-p_r)^2*F_estimate[i,l] + p_r*(1-p_r)*(g*F_estimate[source,source] + (1-g)*both_MH) # option without considering a being linked to non sel background
      
    }
    
    #btwn mh_null_AA and source population
    F_null[i,source]=F_null[source,i]=(1-p_r)*F_estimate[i,source] + p_r*F_estimate[source,source]
    
    #btwn mh_null_AA and never possibly selected pops
    for (k in neverSel) {
      F_null[i,k]=F_null[k,i]=(1-p_r)*F_estimate[i,k]
    }
    
  }
  
  
  
  # also need to modify _aa pops since they are less similar to Neanderthals
  for(i in mh_null_aa) {
    # already modified relationships with mh_null_AA pops
    
    # with aa pops, including self:
    for (j in mh_null_aa) {
      both_MH=(F_estimate[i,j]-g^2*F_estimate[source,source])/(1-g)^2 # prob coalesce given two lineages are type MH
      F_null[i,j]=F_null[j,i]=p_r^2*both_MH + (1-p_r)^2*F_estimate[i,j] + 2*p_r*(1-p_r)*(1-g)*both_MH
    }
    
    # btwn mh_null_aa and never possibly selected population
    for (k in neverSel) {
      F_null[i,k]=F_null[k,i]=p_r*F_estimate[i,k]/(1-g) + (1-p_r)*F_estimate[i,k]
    }
    
    # btwn mh_null_aa and source population
    F_null[i,source]=F_null[source,i]= (1-p_r)*F_estimate[i,source]
    
  }
  
  #### step 2: now modify the _Aa subpops of populations that have them
  for(i in pops[grepl("_Aa",pops)]) { # i = Aa subpop
    AA = gsub("_Aa","_AA",i) # selected subpop
    aa = gsub("_Aa","_aa",i) # neutral subpop
    
    if(AA %in% pops & aa %in% pops) {
      # estimates depend on whether each lineage (or single lineage) sampled in Aa is on same haplotype as selected or non-selected lineage
      F_null[i,i] = 0.25*F_null[aa,aa] + 0.25*F_null[AA,AA] + 0.5*F_null[AA,aa]
      
      for(m in setdiff(pops,i)) {
        F_null[i,m] = F_null[m,i] = 0.5*F_null[AA,m] + 0.5*F_null[aa,m]
      }
    }
  }
  
  return(F_null[real_pops,real_pops])
}

calc_F.N_all = function(F_est) {
  # Computes F^(N) for all genetic distances to the putative selected site
  #
  # Args:
  #	F_est: neutral F estimated from genome-wide data
  #
  # Returns:
  #	list of length numBins of matrices describing probabilities
  # of coalescing within/between pops under
  # scenarios of partitioning populations at a neutral site
  # according to Neanderthal ancestry (F^(N))
  #
  # Have list of length numBins because F^(N) differs
  # according to distance of the neutral site to the selected site
  
  FOmegas = lapply(mid_genetic_distances, function(i) calc_F.N(r=i, F_estimate=F_est))
  return(FOmegas)
}
