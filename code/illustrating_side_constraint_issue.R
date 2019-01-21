library(mgcv)
library(dplyr)
library(tidyr)
library(ggplot2)


#Set the default theme for ggplot objects to theme_bw()
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

#creating a new CO2 data set with just 4 plant levels (to aid plotting)
CO2_new <- CO2%>%
  filter(Plant %in% c("Qn1", "Qc1","Mn1", "Mc1")) %>%
  mutate(Plant_uo=factor(Plant, ordered=FALSE, levels = unique(.$Plant)))
  

modG <- gam(log(uptake) ~ s(log(conc), k=5, bs="tp") +
                  s(Plant_uo, k=5, bs="re"),
                data=CO2_new, method="REML", family="gaussian")

modGS_fs <- gam(log(uptake) ~ s(log(conc), k=5, m=2) + 
                   s(log(conc), Plant_uo, k=5,  bs="fs", m=2),
                 data=CO2_new, method="REML")

#The tensor product model *with* a global smoother
modGS_te <- gam(log(uptake) ~ s(log(conc), k=5, m=2) + 
                   te(log(conc), Plant_uo, k=c(5,4),  bs=c("tp","re"), m=2),
                 data=CO2_new, method="REML")

#The tensor product model *without* a global smoother
modS_te <- gam(log(uptake) ~ te(log(conc), Plant_uo, k=c(5,4),  bs=c("tp","re"), m=2),
               data=CO2_new, method="REML")

modGI <- gam(log(uptake) ~ s(log(conc), k=5, m=2, bs="tp") +
                           s(log(conc), by=Plant_uo, k=5, m=2, bs="tp") +
                           s(Plant_uo, bs="re"),
                 data=CO2_new, method="REML")




plotting_data = crossing(conc = exp(seq(log(min(CO2$conc)), 
                                        log(max(CO2$conc)),
                                        length.out = 100)),
                         Plant_uo = unique(CO2_new$Plant_uo))



modGS_fs_indiv_basis = plotting_data %>%
  bind_cols(as_data_frame(predict(modGS_fs,newdata = ., type = "lpmatrix")))%>%
  select(-contains("(Intercept)"), 
         -contains("s(log(conc))"),
         -contains("s(Plant_ou)")
         )%>%
  mutate(model = "Shared 'fs' smoothers")%>%
  gather(key = basis, value = value, contains("s(log(conc),Plant_uo)"))%>%
  filter(!near(value,0))

modGS_te_indiv_basis = plotting_data %>%
  bind_cols(as_data_frame(predict(modGS_te,newdata = ., type = "lpmatrix")))%>%
  select(-contains("(Intercept)"), 
         -contains("s(log(conc))"),
         -contains("s(Plant_ou)")
  )%>%
  mutate(model = "Shared 'te' smoothers")%>%
  gather(key = basis, value = value, contains("te(log(conc),Plant_uo)"))%>%
  filter(!near(value,0))


modGI_indiv_basis = plotting_data %>%
  bind_cols(as_data_frame(predict(modGI,newdata = ., type = "lpmatrix")))%>%
  select(contains("s(log(conc)):"),Plant_uo, conc)%>%
  mutate(model = "Individual smoothers")%>%
  gather(key = basis, value = value, contains("s(log(conc)):"))%>%
  filter(!near(value,0))

modS_te_indiv_basis = plotting_data %>%
  bind_cols(as_data_frame(predict(modS_te,newdata = ., type = "lpmatrix")))%>%
  select(-contains("(Intercept)"), 
         -contains("s(log(conc))"),
         -contains("s(Plant_ou)")
  )%>%
  mutate(model = "No global term")%>%
  gather(key = basis, value = value, contains("te(log(conc),Plant_uo)"))%>%
  filter(!near(value,0))

basis_function_plot = bind_rows(modGS_fs_indiv_basis,
                                modGS_te_indiv_basis,
                                modGI_indiv_basis)%>%
  mutate(model = factor(model, 
                        levels = c("Shared 'te' smoothers",
                                   "Shared 'fs' smoothers",
                                   "Individual smoothers")))%>%
  ggplot(aes(x= conc, y= value, group = basis))+
  facet_grid(model~Plant_uo)+
  geom_line()+
  scale_x_log10()

print(basis_function_plot)

GS_S_basis_function_plot = bind_rows(modGS_te_indiv_basis,
                                     modS_te_indiv_basis)%>%
  mutate(model = ifelse(model =="Shared 'te' smoothers", 
                        "With global term present",
                        "No global term"))%>%
  ggplot(aes(x= conc, y= value, group = basis))+
  facet_grid(model~Plant_uo)+
  geom_line()+
  scale_x_log10()

print(GS_S_basis_function_plot)

#Plotting individual global and deviance curves ####
modGS_fs_terms = plotting_data %>%
  bind_cols(as_data_frame(predict(modGS_fs,newdata = ., type = "terms")))%>%
  mutate(model = "Shared 'fs' smoothers")%>%
  rename(global_smooth = `s(log(conc))`,
         Plant_deviation = `s(log(conc),Plant_uo)`)%>%
  mutate(total_term =  global_smooth + Plant_deviation)%>%
  gather(key =term, value = value,global_smooth,  Plant_deviation,total_term)

modGS_te_terms = plotting_data %>%
  bind_cols(as_data_frame(predict(modGS_te,newdata = ., type = "terms")))%>%
  mutate(model = "Shared 'te' smoothers")%>%
  rename(global_smooth = `s(log(conc))`,
         Plant_deviation = `te(log(conc),Plant_uo)`)%>%
  mutate(total_term =  global_smooth + Plant_deviation)%>%
  gather(key =term, value = value,global_smooth,  Plant_deviation,total_term)


modGI_terms = plotting_data %>%
  bind_cols(as_data_frame(predict(modGI,newdata = ., type = "terms")))%>%
  mutate(model = "Individual smoothers",
         Plant_deviation = rowSums(select(., `s(log(conc)):Plant_uoQn1`:`s(Plant_uo)`))
         )%>%
  rename(global_smooth = `s(log(conc))`)%>%
  mutate(total_term =  global_smooth + Plant_deviation)%>%
  gather(key =term, value = value,global_smooth,  Plant_deviation,total_term)
         

term_plot = bind_rows(modGS_fs_terms, modGS_te_terms, modGI_terms)%>%
  mutate(model = factor(model, 
                        levels = c("Shared 'te' smoothers",
                                   "Shared 'fs' smoothers",
                                   "Individual smoothers")),
         Plant_uo = ifelse(term=="global_smooth", 
                           "global smooth", 
                           as.character(Plant_uo)))%>%
  mutate(Plant_uo = factor(Plant_uo, 
                           levels = c("global smooth", 
                                      "Qn1",
                                      "Qc1", 
                                      "Mn1", 
                                      "Mc1")))%>%
  ggplot(aes(x= conc, y= value, group = Plant_uo,color = Plant_uo))+
  facet_grid(model~term)+
  geom_line()+
  scale_x_log10()+
  scale_color_brewer(palette = "Set1")

print(term_plot)
