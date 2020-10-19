# Plot profile composite likelihood surface of t_b and s
# Also just print the values for the maximum composite likelihood estimates of t_b and s

SelModel_CLs=apply(results,c(1,2),as.numeric) # makes "null" and "invalid" entries NA, but selection always gets higher CL values anyway

# PROFILE LIKELIHOOD SURFACE FOR T_B

# maximum composite likelihood for each value of t_b (get maximum value of each column)
tb_profile = apply(SelModel_CLs,2,max,na.rm=T) 
tb_MCL_profile = data.frame(tb=times_btwn,MCL=tb_profile)

tb_MCL_profile %>% 
  ggplot(aes(x=tb,y=MCL))+
  geom_line()+
  labs(x="proposed waiting time until selection",
       y="log maximum composite likelihood") +
  scale_x_continuous(breaks=times_btwn) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1),
        panel.grid.minor = element_blank())

# PROFILE LIKELIHOOD SURFACE FOR S

# maximum composite likelihood for each value of s (get maximum value of each row)
sel_profile = apply(SelModel_CLs,1,max,na.rm=T) 
sel_MCL_profile = data.frame(sel=sels,MCL=sel_profile)

sel_MCL_profile %>% 
  ggplot(aes(x=sel,y=MCL))+
  geom_line()+
  labs(x="proposed selection coefficient",
       y="log maximum composite likelihood") +
  scale_x_continuous(breaks=sels) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1),
        panel.grid.minor = element_blank())

# ESTIMATES FOR T_B AND S
times_btwn[which.max(tb_profile)] # t_b estimate
sels[which.max(sel_profile)] # s esitmate

