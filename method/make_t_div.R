# this script describes time since divergence of pairs of populations
# time is in units of real-time generations

# It makes a matrix describing divergence times (how long ago populations split)
# between pairs of populations. Rows and columns correspond to populations

# make divergence time matrix 

t_div=outer(intro_pops,intro_pops, Vectorize(function(pop1,pop2){
  if(gsub("_..$","",pop1)==gsub("_..$","",pop2)){
    return(0)
  } else if(nest_index[pop1]==nest_index[pop2]){
    div_times[nest_index[pop1]]
  } else if(nest_index[pop1] < nest_index[pop2]){
    return(div_times[nest_index[pop1]])
  } else { # nest_index[pop1] > nest_index[pop2]
    div_times[nest_index[pop2]]
  }
} ))
dimnames(t_div)=list(intro_pops,intro_pops)


