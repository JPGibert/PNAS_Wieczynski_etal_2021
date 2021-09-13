## Main Analysis for Wieczynski et al. 2021, PNAS

### clear all data
rm(list=ls())
#
### set data directory / load packages / functions / color palettes -----

data_dir<-"~/PNAS_Wieczynski_etal_2021_Master/" # add master file to home directory

packages = c("tidyverse", "broom", "data.table", "RColorBrewer",  
             "colorRamps", "rstatix", "nls.multstart", "scales", "gtools", 
             "car", "yhat", "mixtools", "ggpubr", "lemon", 
             "Polychrome", "viridisLite", "lmodel2", "boot", "car",
             "ggridges") 

package.check <- lapply(packages, FUN = function(x) {
        if (!require(x, character.only = TRUE)) {
                install.packages(x, dependencies = TRUE)
                library(x, character.only = TRUE)
        }
})


p_val_func<-function(modelobject){
  f<-summary(modelobject)$fstatistic
  p<-pf(f[1],f[2],f[3],lower.tail=F)
  return(p)
}

Mean<-function(x){
  z<-na.omit(x)
  sum(z)/length(z)
}
SD<-function(x){
  z<-na.omit(x)
  n<-length(z)
  sqrt((sum((z-mean(z))^2))/n)
}
Variance<-function(x){
  z<-na.omit(x)
  n<-length(z)
  (sum((z-mean(z))^2))/n
}
Skewness<-function(x){
  z<-na.omit(x)
  n<-length(z)
  stddev<-SD(z)
  sum(((z-mean(z))/stddev)^3)/n
}
Kurtosis<-function(x){
  z<-na.omit(x)
  n<-length(z)
  stddev<-SD(z)
  (sum(((z-mean(z))/stddev)^4)/n)-3
}

cent_rescale<-function(val){(val-Mean(val))/SD(val)} # /SD(val)

T_color_palette<-brewer.pal(3, "Set1")[2:1]


spp_order<-c("Paramecium bursaria", "Paramecium multimicronucleatum", "Euplotes sp.","Blepharisma sp.", "Urocentrum turbo", 
             "Colpidium striatum", "Paramecium aurelia", "Paramecium caudatum", "Tillina magna", "Halteria grandinella", 
             "Cyclidium glaucoma", "Colpoda steinii", "Glaucoma sp.", "Tetrahymena pyriformis") # , "Strombidium", "Urocentrum turbo"
spp_order_abbr<-c("P. bursaria", "P. multimicronucleatum", "Euplotes sp.", "Blepharisma sp.", "U. turbo", 
                  "C. striatum", "P. aurelia", "P. caudatum", "T. magna", "H. grandinella", 
                  "C. glaucoma", "C. steinii", "Glaucoma sp.", "T. pyriformis") # , "Strombidium", "U. turbo"

spp_color_palette<-c("#A6CEE3", "#1F78B4", "#043f70", "#33A02C", "#46DBBD", "#FB9A99", "#E31A1C", "#8c1029", "#FFD500", "#C69354", "#FF7F00", "#c91797", "#CAB2D6", "#6A3D9A")

spp_name_df<-data.frame(Species=spp_order, Species_abbr=spp_order_abbr, color=spp_color_palette[1:length(spp_order)])
spp_labels<-spp_name_df$Species_abbr
names(spp_labels)<-spp_name_df$Species

#

###-----
###
### Demographic parameter analysis -----


## load data

demo_data<-read.csv(paste0(data_dir, "data/Demography.csv")) %>%
  filter(Species %in% spp_order, !N.initial %in% c(1, 31.3)) %>%
  group_by(Species, Temperature) %>%
  mutate(R=ifelse(R>=mean(R, na.rm=T)-3*sd(R, na.rm=T) & R<=mean(R, na.rm=T)+3*sd(R, na.rm=T), R, NA)) %>%
  ungroup %>%
  mutate(Temperature=as.factor(Temperature), Species=factor(Species, levels=spp_order)) 




## regressions

stats_lm_1<-demo_data %>%
  group_by(Species) %>% 
  do(main = lm(R ~ N.initial + Temperature, data = ., na.action=na.omit),
     interaction = lm(R ~ N.initial * Temperature, data = ., na.action=na.omit)) %>% 
  gather(model, fit, -Species) %>%
  rowwise %>%
  mutate(aic=AIC(fit), 
         Intercept=summary(fit)$coefficients[1,1], 
         Intercept_se=summary(fit)$coefficients[1,2], 
         N.init=summary(fit)$coefficients[2,1], 
         N.init_se=summary(fit)$coefficients[2,2], 
         Temp=summary(fit)$coefficients[3,1], 
         Temp_se=summary(fit)$coefficients[3,2], 
         N.init.Temp=ifelse(model=="main", NA, summary(fit)$coefficients[4,1]),
         N.init.Temp_se=ifelse(model=="main", NA, summary(fit)$coefficients[4,2]),
         
         Intercept_N.init_cov=vcov(fit)[2,1],
         Intercept_Temp_cov=vcov(fit)[3,1],
         N.init_N.init.Temp_cov=ifelse(model=="main", NA, vcov(fit)[4,2]),
         
         model_pval=p_val_func(fit),
         Intercept_pval=summary(fit)$coefficients[1,4], 
         N.init_pval=summary(fit)$coefficients[2,4], 
         Temp_pval=summary(fit)$coefficients[3,4], 
         N.init.Temp_pval=ifelse(model=="main", NA, summary(fit)$coefficients[4,4]),
         int_sig=ifelse(N.init.Temp_pval<=0.05, "*", "ns")) %>%
  group_by(Species) %>%
  mutate(min_aic=min(aic)) %>% 
  rowwise() %>%
  mutate(aic_diff=aic-min_aic, aic_support=ifelse(aic_diff<4, "*", "ns"), sig=ifelse(model_pval<=0.05, "*", "ns")) %>%
  group_by(Species) %>%
  mutate(n_supported=3-length(unique(aic_support))) %>%
  ungroup %>%
  arrange(Species) %>%
  mutate(support_sig=ifelse(aic_support=="*" & sig=="*", "*", NA), 
         best_model_overall=ifelse(support_sig=="*" & n_supported==2 & model=="interaction", NA, support_sig),
         best_model=ifelse(model=="main" & lead(N.init.Temp_pval)>0.05, "*", ifelse(model=="interaction" & N.init.Temp_pval<=0.05, "*", NA))
         )


## export regression summary
regression_export<-stats_lm_1 %>%
  filter(best_model=="*") %>%
  select(c(1, 2, 16, 5:12, 17:20)) %>%
  rename(Intercept_Estimate=Intercept, N.init_Estimate=N.init, Temp_Estimate=Temp, N.init.Temp_Estimate=N.init.Temp) %>%
  pivot_longer(!c(Species, model, model_pval), names_to=c("Variable", ".value"), names_pattern="(.+)_(.+$)") %>%
  arrange(Species, factor(Estimate, levels=c("Intercept", "N.init", "Temp", "N.init.Temp"))) %>%
  mutate(across(c(model_pval, Estimate, se, pval), ~sprintf("%.3f", round(.x, 3))),
         model_pval=ifelse(model_pval<0.001, "<0.001", model_pval), pval=ifelse(pval<0.001, "<0.001", pval)) %>%
  rename("Best model"=model, "Model p-value"=model_pval, "Std. Error"=se, "p-value"=pval)

fwrite(regression_export, paste0(data_dir, "data/demo_regression_summary.csv"))



## calculate "r" and "K"

axis_limits<-demo_data %>%
  group_by(Species) %>%
  summarize(min_y=round(min(R, na.rm=T), 1), max_y=round(max(R, na.rm=T), 1)) %>%
  mutate(min_x=ifelse(Species=="Cyclidium glaucoma", 125, 4), max_x=ifelse(Species=="Cyclidium glaucoma", 625, 20))

axis_range<-axis_limits %>%
  gather(y_lim, y, min_y, max_y) %>%
  gather(x_lim, x, min_x, max_x) %>%
  dplyr::select(Species, x, y) %>%
  distinct()


demo_params_full<-expand.grid(Species=unique(demo_data$Species), Temperature=c("22", "25")) %>%
  left_join(stats_lm_1) %>%
  mutate(b=ifelse(Temperature=="22", Intercept, Intercept+Temp),
         m=ifelse(Temperature=="25" & model=="interaction", N.init+N.init.Temp, N.init),
         r=b,
         K=-b/m,
         b_se=ifelse(Temperature=="22", Intercept_se, sqrt(Intercept_se^2 + Temp_se^2 + 2*Intercept_Temp_cov)),
         m_se=ifelse(Temperature=="25" & model=="interaction", sqrt(N.init_se^2 + N.init.Temp_se^2 + 2*N.init_N.init.Temp_cov), N.init_se),
         r_se=b_se,
         K_se=abs(K)*sqrt((b_se/b)^2 + (m_se/m)^2), 
         r_sig=ifelse(Intercept_pval<=0.05, "*", "ns"),
         K_sig=ifelse(Intercept_pval<=0.05 & N.init_pval<=0.05, "*", "ns"),
         K_sig=ifelse(Species=="Paramecium caudatum",  "*", K_sig)
         ) %>%
  left_join(axis_limits) %>%
  mutate(y1=m*min_x+b, y2=m*max_x+b) %>%
  mutate(Species=factor(Species, levels=spp_order)) %>%
  arrange(Species)



rK_sig<-demo_params_full %>%
  filter(best_model=="*") %>%
  select(Species, r_sig, K_sig) %>%
  distinct %>%
  gather(param, sig, r_sig, K_sig) %>%
  mutate(param=word(param, 1, sep="_"))


demo_params_reduced<-demo_params_full %>%
  filter(best_model=="*") %>% 
  dplyr::select(Species, Temperature, r, K) %>%
  gather(param, param_val, r, K) %>%
  mutate(log10_param_val=log10(param_val),
         Species=factor(Species, levels=spp_order),
         param=factor(param, levels=c("r", "K"))) %>%
  left_join(rK_sig)

T_diff_df<-demo_params_reduced %>%
  dplyr::select(-log10_param_val) %>%
  spread(Temperature, param_val) %>%
  rename(T25='25', T22='22') %>%
  mutate(T_diff=abs(T25-T22)) %>% 
  dplyr::select(-T22, -T25) %>%
  rename(param_val=T_diff) %>%
  mutate(param = dplyr::recode(param, "K" = "K_diff", "r" = "r_diff"), log10_param_val=log10(param_val)) 

demo_gdata<-demo_params_reduced %>%
  filter(param=="r") %>%
  mutate(Temperature=dplyr::recode(Temperature, "22"=22, "25"=25))





## calculate confidence intervals

conf_int_func<-function(model, N.initial, Temperature){predict(model, newdata=data.frame(N.initial=N.initial, Temperature=Temperature), interval="confidence", level = 0.95)}

conf_int_df<-expand_grid(Species=spp_order, N.initial=seq(4, 20, length.out=100), Temperature=c("22", "25")) %>%
  mutate(N.initial=ifelse(Species=="Cyclidium glaucoma", rescale(N.initial, to=c(125.2, 626), from=c(4,20)), N.initial)) %>%
  left_join(dplyr::select(stats_lm_1, Species, model, fit, best_model)) %>%
  filter(best_model=="*") %>%
  rowwise() %>%
  mutate(lwr=conf_int_func(fit, N.initial, Temperature)[2], upr=conf_int_func(fit, N.initial, Temperature)[3]) %>%
  dplyr::select(-fit) %>%
  mutate(Species=factor(Species, levels=spp_order))




## export demo parameters

demo_export<-demo_params_reduced %>%
  mutate(param=paste0(param, Temperature)) %>%
  dplyr::select(-log10_param_val, -Temperature, -sig)  %>%
  spread(param, param_val) %>%
  mutate(ID=row_number()) %>%
  dplyr::select(ID, Species, r22, r25, K22, K25)

fwrite(demo_export, paste0(data_dir, "data/demo_parameters.csv"))



# plot "r" vs. "N initial"  (individual panels)

for(i in 1:length(spp_order)){ 
  temp_sp<-spp_order[i]
  temp_conf_int_df<-filter(conf_int_df, Species==temp_sp)
  temp_demo_params_full<-filter(demo_params_full, Species==temp_sp) %>%
    filter(best_model=="*")
  temp_demo_data<-filter(demo_data, Species==temp_sp)
  x_breaks<-if(temp_sp=="Cyclidium glaucoma") {seq(125, 625, length.out=3)} else {seq(4, 20, length.out=3)}
  y_min_raw<-min(temp_demo_data$R, na.rm=T)
  y_max_raw<-max(temp_demo_data$R, na.rm=T)
  y_min<-ifelse(round(y_min_raw, 1)>y_min_raw, round(y_min_raw, 1)-0.1, round(y_min_raw, 1))
  y_max<-ifelse(round(y_max_raw, 1)<y_max_raw, round(y_max_raw, 1)+0.1, round(y_max_raw, 1))
  y_breaks<-seq(y_min, y_max, length.out=3)
  
  assign(paste("g", i, sep="_"),
         ggplot() +
           geom_ribbon(data=temp_conf_int_df, aes(x=N.initial, ymin=lwr, ymax=upr, group=Temperature), fill='gray', alpha=0.2) +
           geom_segment(data=temp_demo_params_full, aes(x=min_x, xend=max_x, y=y1, yend=y2, linetype=Temperature), color=spp_color_palette[i], size=.8) +
           geom_point(data=temp_demo_data, aes(N.initial, R, shape=Temperature), color=spp_color_palette[i], size=2)+
           scale_linetype_manual(values=c(1, 2)) +
           scale_shape_manual(values=c(19, 21)) +
           scale_x_continuous(breaks=x_breaks, expand = expansion(mult = c(.05, .05))) +
           scale_y_continuous(breaks=y_breaks, limits=c(y_min, y_max), expand = expansion(mult = c(.05, .05))) +
           coord_capped_cart(bottom='right', left='top', gap=0) +
           labs(x=Initial~density~(cells~mL^-1), y=Per~capita~growth~rate~(cells~cell^-1~d^-1)) + 
           theme(plot.background = element_blank(), panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank(),
                 axis.line=element_line(),
                 axis.text=element_text(size=13),
                 axis.title=element_blank(),
                 strip.text=element_blank(), strip.background=element_blank(), 
                 legend.position="none",
                 legend.key=element_blank()) 
  )
}

gg_blank<-ggplot() + theme_void()


svg(file=paste0(data_dir, "graphics/raw_graphics/rK_regressions.svg"), bg=NA, width=8, height=10) 
ggarrange(g_1, g_2, g_3, g_4, g_5, g_6, g_7, g_8, g_9, g_10, gg_blank, gg_blank, g_11, g_12, gg_blank, gg_blank, g_13, g_14,
          nrow=5, ncol=4, align="hv") 
dev.off()





## plot "r" & "K" 

gdata<-demo_params_reduced %>%
  mutate(param_label=ifelse(param=="r", "r", "K"),
         param_label=factor(param_label, levels=c("r", "K")))


svg(file=paste0(data_dir, "graphics/raw_graphics/rK_Temperatures.svg"), bg=NA, width=4, height=6.4) 
ggplot()+ 
  geom_line(data=gdata, aes(Temperature, param_val, color=Species, group=Species))+
  geom_point(data=gdata, aes(Temperature, param_val, color=Species, group=Species, shape=Temperature), size=2.5, stroke=1.5)+
  scale_shape_manual(values=c(19, 21)) +
  scale_x_discrete(expand=c(.2, .1)) +
  scale_y_log10(expand = expansion(mult = c(.05, .05))) +
  annotation_logticks(sides="l", short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm")) +
  scale_color_manual(values=spp_color_palette) +
  labs(x="Temperature (C)") +
  facet_wrap(~param_label, scales="free", strip.position="top", labeller=label_parsed) + 
  theme(plot.background = element_blank(), panel.background = element_blank(), panel.border = element_blank(),
        panel.grid.major.x=element_line(color="gray"),
        axis.line.y=element_line(),
        axis.text.x=element_text(size=16), axis.text.y=element_text(size=14),
        axis.title=element_text(size=16), 
        axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        strip.text=element_text(size=26, face="bold"), strip.placement = "outside", strip.background=element_blank(),
        legend.position="none",
        legend.text=element_text(size=14, face="italic"), legend.title=element_text(size=16), 
        legend.key=element_blank())
dev.off()









#
### TPC analysis -----


## load data
TPC_df_raw<-read.csv(paste0(data_dir, "data/TPCs.csv"))


# model fitting

fit_scale_factor<-10

TPC_df<-TPC_df_raw %>%
  filter(Species %in% spp_order,
         Nf>=1) %>%
  group_by(Species) %>%
  mutate(r_scale=10,
         log_r=log(r+r_scale)) %>%
  ungroup


nls_model<-TPC_df %>%
  group_by(Species) %>%
  do(TPC_fit = nls_multstart(log_r ~ a + (E_a/(8.6*10^-5))*(1/298.15-1/(Temperature+273.15)) - log(1+exp((E_d/(8.6*10^-5))*(1/Th-1/(Temperature+273.15)))),
                             data = .,
                             iter = 500,
                             start_lower = c(a=-10, E_a=0.1, E_d=0.5, Th=285),
                             start_upper = c(a=10, E_a=4, E_d=10, Th=330),
                             supp_errors = 'Y',
                             na.action = na.omit,
                             lower = c(a=-10, E_a=0, E_d=0, Th=0))) %>%
  rowwise() %>%
  mutate(a=coef(TPC_fit)[[1]], E_a=coef(TPC_fit)[[2]], E_d=coef(TPC_fit)[[3]], Th=coef(TPC_fit)[[4]], T_opt=E_d*Th/(E_d+8.6e-5*Th*log(E_d/E_a-1)))


## export TPC model parameters
TPC_model_export<-dplyr::select(nls_model, -TPC_fit)

fwrite(TPC_model_export, paste0(data_dir, "data/TPC_model_parameters.csv"))



range_r<-TPC_df %>%
  group_by(Species) %>%
  summarize(min_r=min(r, na.rm=T), max_r=max(r, na.rm=T), range_r=max_r-min_r)

TPC_eqn<-function(a, E_a, E_d, Th, Temperature, r_scale){exp(a + (E_a/(8.6*10^-5))*(1/298.15-1/(Temperature+273.15)) - log(1+exp((E_d/(8.6*10^-5))*(1/Th-1/(Temperature+273.15)))))-r_scale}

TPC_predicted<-expand.grid(Species=spp_order, Temperature=seq(0, 50, length.out=1000)) %>%
  left_join(dplyr::select(nls_model, -TPC_fit)) %>%
  left_join(distinct(dplyr::select(TPC_df, Species, r_scale))) %>%
  mutate(r=TPC_eqn(a, E_a, E_d, Th, Temperature, r_scale)) %>%
  filter(r>=-1)

TPC_predicted_trunc<-TPC_predicted %>%
  left_join(range_r) %>%
  filter(r>=min_r)

TPC_summary_spread<-TPC_predicted %>%
  group_by(Species) %>%
  mutate(CT_min=ifelse(lag(r)<0 & r>0, Temperature, NA),
         CT_max=ifelse(lag(r)>0 & r<0, Temperature, NA),
         r_peak=max(r)) %>%
  filter(!is.na(CT_min) | !is.na(CT_max)) %>%
  dplyr::select(Species, CT_min, CT_max, r_peak) %>%
  ungroup() %>%
  gather(param, param_val, -Species) %>%
  drop_na %>%
  distinct %>%
  spread(param, param_val) %>%
  left_join(dplyr::select(nls_model, Species, E_a, E_d, T_opt)) %>%
  mutate(T_opt=T_opt-273.15, T_range=CT_max-CT_min, TPC_asymmetry=abs((T_opt-CT_min)-(CT_max-T_opt))) %>%
  arrange(r_peak) 

TPC_params<-c("E_a", "E_d", "r_peak", "CT_max", "CT_min", "T_opt", "T_range", "TPC_asymmetry")

TPC_summary<-TPC_summary_spread %>%
  gather(param, param_val, -Species) %>%
  mutate(log10_param_val=log10(param_val),
         param=factor(param, levels=TPC_params)) %>%
  left_join(range_r)

TPC_summary_gdata<-TPC_summary %>%
  dplyr::select(Species, param, param_val) %>%
  spread(param, param_val) %>%
  left_join(spp_name_df) %>%
  arrange(T_opt) %>%
  mutate(Species=factor(Species, levels=rev(unique(.$Species)))) 



TPC_summary_spread$Species<-factor(TPC_summary_spread$Species, levels=spp_order)
TPC_summary$Species<-factor(TPC_summary$Species, levels=spp_order)
TPC_predicted$Species<-factor(TPC_predicted$Species, levels=spp_order)
TPC_predicted_trunc$Species<-factor(TPC_predicted_trunc$Species, levels=spp_order)
TPC_df$Species<-factor(TPC_df$Species, levels=spp_order)


## export TPC parameters
TPC_export<-TPC_summary_spread %>%
  mutate(across(-Species, ~sprintf("%.2f", round(.x, 2)))) %>%
  select(Species, r_peak, T_opt, CT_min, CT_max, T_range, TPC_asymmetry, E_a, E_d)

fwrite(TPC_export, paste0(data_dir, "data/TPC_parameters.csv"))




# plot all species' TPCs together
TPC_predicted$Species<-factor(TPC_predicted$Species, levels=rev(spp_order))

g_1<-
ggplot()+
  geom_hline(yintercept=0, color="gray", size=0.7, linetype=2)+
  geom_line(data=TPC_predicted, aes(Temperature, r, color=Species), size=1)+
  scale_color_manual(values=rev(spp_color_palette)) +
  scale_x_continuous(limits=c(0, 43), sec.axis = dup_axis()) +
  labs(x="Temperature (C)", y=Intrinsic~growth~rate~(r)~(cells~cell^-1~d^-1)) +
  theme(plot.background=element_blank(), panel.background=element_blank(),panel.grid=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        axis.text=element_text(size=14), axis.title=element_text(size=16),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.ticks.length.x.top=unit(0.15, "cm"),
        plot.margin=unit(c(0,0,0,0), "cm"),
        aspect.ratio=0.55,
        legend.key=element_blank(),
        legend.position="none",
        legend.text=element_text(size=12, face="italic"))


g_2<-
ggplot()+
  geom_segment(data=TPC_summary_gdata, aes(x=CT_min, xend=CT_max, y=Species, yend=Species, color=Species), lineend="round", size=1.2)+
  geom_point(data=TPC_summary_gdata, aes(T_opt, Species, color=Species), size=2.75) +
  scale_color_manual(values=rev(TPC_summary_gdata$color)) +
  xlim(0, 43) +
  theme(plot.background=element_blank(), panel.grid=element_blank(),
        panel.background=element_blank(), panel.border=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        legend.key=element_blank(),
        legend.position="none",
        plot.margin=unit(c(0,0,0,0), "cm"))


svg(file=paste0(data_dir, "graphics/raw_graphics/TPCs.svg"), width=7.5, height=6, bg=NA)
ggarrange(g_2, g_1,
          nrow=2, ncol=1, heights=c(.4, 1), widths=1, align="v") 
dev.off()







#
### compare results from rK and TPC experiments -----


TPC_predicted_22_25<-expand.grid(Species=spp_order, Temperature=c(22, 25)) %>%
  left_join(dplyr::select(nls_model, -TPC_fit)) %>%
  left_join(distinct(dplyr::select(TPC_df, Species, r_scale))) %>%
  mutate(r=TPC_eqn(a, E_a, E_d, Th, Temperature, r_scale)) %>%
  dplyr::select(Species, Temperature, r)

rK_TPC_comparison<-demo_gdata[1:4] %>%
  spread(param, param_val) %>%
  rename(r_rK=r) %>%
  left_join(TPC_predicted_22_25) %>%
  rename(r_TPC=r) %>%
  mutate(Temperature=as.factor(Temperature))
rK_TPC_comparison$Species<-factor(rK_TPC_comparison$Species, levels=spp_order)

fit<-lm( r_TPC ~ r_rK , data = rK_TPC_comparison, na.action=na.omit)
summary(fit)

stat_text<-list(paste('R^2 == ', signif(summary(fit)$adj.r.squared, 2)), 
                paste('p < ', '10^-11'), 
                paste('m == ', signif(coef(fit)[[2]], 3)))



svg(file=paste0(data_dir, "graphics/raw_graphics/rK_TPC_comparison.svg"), bg=NA, width=4.2, height=4.2)
ggplot() +
  geom_abline(intercept=0, slope=1, color="black", size=0.7) +
  geom_abline(intercept=coef(fit)[[1]], slope=coef(fit)[[2]], color="gray40", linetype=2, size=0.7) +
  geom_point(data=rK_TPC_comparison, aes(r_rK, r_TPC, shape=Temperature, color=Species), size=3, stroke=.9) +
  annotate("text", x = 3.0, y = c(1.4, 0.9, 0.4), label = stat_text, parse=T, hjust = 0, size=5) +
  scale_color_manual(values=spp_color_palette) +
  scale_shape_manual(values=c(19,1), guide=F) +
  lims(x=c(0, 4.2), y=c(0, 4.2)) +
  labs(x=italic(r)~" (density-dependence assays)", y=italic(r)~" (TPC assays)") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(panel.background = element_blank(), plot.background = element_blank(), 
        panel.border = element_rect(color = "black", size=1, fill=NA),
        panel.grid=element_blank(),
        aspect.ratio=1,
        axis.text=element_text(size=14), axis.title=element_text(size=16),
        legend.position="none",
        legend.key=element_blank(),
        legend.text=element_text(size=10, face="italic"), legend.title=element_text(size=12))
dev.off()






#
### Trait analysis -----


trait_list<-c("Biovolume_Spheroid_GEO", "Geodesic_Aspect_Ratio", "Sigma.Intensity")

## load and process raw trait data
phenotype_data<-
  read.csv(paste0(data_dir, "data/Traits.csv")) %>%
  distinct() %>%
  dplyr::select(-Original.Reference.ID, -Timestamp, -Date, -Filter.Score, -Feret.Angle.Min, -Feret.Angle.Max, -Fiber.Curl) %>% 
  rename(Species=Class) %>% 
  mutate(Magnification=word(Species, 2, sep="_"),
         Magnification=ifelse(is.na(Magnification), ifelse(Species %in% c("Euplotes sp.", "Tillina magna"), "x4", "x10"), Magnification),
         Species=word(Species, 1, sep="_"),
         
         Geodesic.Thickness=ifelse(Geodesic.Thickness==Geodesic.Length, NA, Geodesic.Thickness), # remove exact spheres
         Geodesic.Length=ifelse(Geodesic.Thickness==Geodesic.Length, NA, Geodesic.Length),       # remove exact spheres
         
         Symmetry=1.1-Symmetry,
         Convexity=1.1-Convexity,
         Circularity=1.1-Circularity,
         Eccentricity_GEO=1.1-sqrt(1-((Geodesic.Thickness/2)^2/(Geodesic.Length/2)^2)),
         
         Volume_ABD=(4/3)*pi*(Diameter..ABD./2)^3,
         Volume_ESD=(4/3)*pi*(Diameter..ESD./2)^3,
         Transparency=1-(Diameter..ABD./Diameter..ESD.),
         Biovolume_Cylinder_GEO=pi*(Geodesic.Thickness/2)^2*Geodesic.Length,
         Biovolume_Spheroid_GEO=(4/3)*pi*(Geodesic.Thickness/2)^2*(Geodesic.Length/2),
         Biovolume_Spheroid_FERET=(pi/6)*Width^2*Length,
         Biovolume_Sphere_ABD=(4/3)*pi*(Diameter..ABD./2)^3, 
         Geodesic_Aspect_Ratio=Geodesic.Thickness/Geodesic.Length,
         Aspect_Ratio_FERET=Width/Length 
         ) %>% 
  group_by(Species, Magnification) %>% 
  mutate(rep=row_number()) %>%
  ungroup() %>%
  gather(Metric, Metric_val, -c(Species, Magnification, rep)) %>%
  mutate(Metric_val=ifelse(Species %in% c("Euplotes sp.", "Tillina magna"),
                     ifelse(Magnification=="x10",
                            ifelse(str_detect(Metric, "Red|Green|Blue|Intensity"), Metric_val, NA),
                            ifelse(str_detect(Metric, "Red|Green|Blue|Intensity"), NA, Metric_val)),
                     Metric_val)) %>%
  spread(Metric, Metric_val) %>% 
  filter(Species %in% spp_order)


# log10 transformations and outlier removal
phenotype_df<-phenotype_data %>%
  gather(Trait, Trait_value, -c(Species, rep, Magnification)) %>%
  mutate(Trait_value=log10(Trait_value)) %>% # log-transform trait values (ALL traits)
  group_by(Species, Magnification, Trait) %>% 
  mutate(Trait_value=ifelse(Trait_value>=Mean(Trait_value)-3*SD(Trait_value) & Trait_value<=Mean(Trait_value)+3*SD(Trait_value), Trait_value, NA)) %>%
  ungroup() 


# calculate species-level trait moments (log10)
phenotype_summary<-phenotype_df %>%
  group_by(Species, Trait) %>% 
  summarize(Mean=Mean(Trait_value), Variance=Variance(Trait_value), Skewness=Skewness(Trait_value), Kurtosis=Kurtosis(Trait_value)) %>%
  gather(Moment, Moment_value, Mean:Kurtosis) %>%
  spread(Trait, Moment_value)


## export trait data (transformed back to original scale (i.e., not log10))
trait_min_max_data<-phenotype_df %>%
  filter(Trait %in% trait_list) %>%
  mutate(Trait_value=10^Trait_value,
         Species=factor(Species, levels=spp_order)) %>%
  group_by(Species, Trait) %>%
  summarize(Min=min(Trait_value, na.rm=T), Max=max(Trait_value, na.rm=T)) %>%
  ungroup()

trait_mean_data<-phenotype_summary %>%
  filter(Moment=="Mean") %>%
  select(Species, all_of(trait_list)) %>%
  mutate(across(all_of(trait_list), ~10^.x))

fwrite(trait_mean_data, paste0(data_dir, "data/traits_summary.csv"))



## plot trait distributions 

gg_func<-function(i){ 
  i<-i
  temp_t<-trait_list[[i]]
  temp_gdata<-phenotype_df %>%
    filter(Trait==temp_t) %>%
    mutate(Trait_value=10^Trait_value,
           Species=factor(Species, levels=rev(spp_order)))
  temp_summary_data<-phenotype_summary %>%
    filter(Moment=="Mean") %>%
    select(Species, temp_t) %>%
    rename(Trait_value=temp_t) %>%
    mutate(Trait_value=10^Trait_value)
  
  
  x_breaks<-if(str_detect(temp_t, "Biovolume")) {c(10^3, 10^4, 10^5, 10^6)} else {waiver()} 
  x_ticks<-if(str_detect(temp_t, "Biovolume")) {trans_format("log10", math_format(10^.x))} else {waiver()} 
  x_lab<-case_when(temp_t=="Biovolume_Spheroid_GEO" ~ "Volume~(mu*'m'^3)", 
                   temp_t=="Geodesic_Aspect_Ratio" ~ "Aspect~ratio", 
                   temp_t=="Sigma.Intensity" ~ "Contrast")
  max_density<-case_when(temp_t=="Biovolume_Spheroid_GEO" ~ 3.54,
                         temp_t=="Geodesic_Aspect_Ratio" ~ 7.68,
                         temp_t=="Sigma.Intensity" ~ 19.6)
  y_val<-max_density+max_density*0.075
  
  return(
    ggplot() +
      geom_density(data=temp_gdata, aes(Trait_value, ..density.., color=Species, fill=Species), size=0.8, alpha=0.2) +
      geom_segment(data=temp_summary_data, aes(x=Trait_value, xend=Trait_value, y=y_val, yend=y_val+max_density*0.08, color=Species), lineend="round", size=1.1) +
      scale_color_manual(values=rev(spp_color_palette)) +
      scale_fill_manual(values=rev(spp_color_palette)) +
      labs(x=parse(text=x_lab), y="Frequency") +
      scale_x_log10(limits=c(min(temp_gdata$Trait_value), max(temp_gdata$Trait_value)), labels=x_ticks, expand = c(0.0,0.0)) +
      coord_cartesian(clip = "off") +
      annotation_logticks(sides="b", short = unit(-0.1, "cm"), mid = unit(-0.1, "cm"), long = unit(-0.2, "cm")) +
      theme(axis.text.x=element_text(size=14), axis.title.x=element_text(size=16),
            axis.text.y=element_blank(), axis.ticks.y=element_blank(),
            axis.title.y=element_blank(),
            axis.line.x=element_line(),
            strip.text=element_blank(),
            strip.background=element_blank(),
            legend.position="none",
            plot.margin = unit(c(.5, 0, 0, 0), "cm"),
            plot.background=element_blank(), panel.grid=element_blank(), panel.background=element_blank(), panel.border=element_blank())
  )
}

for(i in 1:length(trait_list)){assign(paste("g", i, sep="_"), gg_func(i))}


svg(file=paste0(data_dir, "graphics/raw_graphics/trait_distributions.svg"), bg=NA, width=5, height=7)
ggarrange(g_1, g_2, g_3,
          nrow=3, ncol=1, align="hv") 
dev.off()





#

### Use traits to predict TPCs -----


traits_trueTPC<-phenotype_summary[c("Species", "Moment", unique(phenotype_df$Trait))] %>%
  dplyr::select(Species, Moment, all_of(trait_list)) %>%
  filter(Moment=="Mean") %>%
  left_join(TPC_model_export) %>%
  ungroup() %>%
  select(-Moment)



## simple regressions

ind_reg<-traits_trueTPC %>%
  pivot_longer(Biovolume_Spheroid_GEO:Sigma.Intensity, names_to="trait", values_to="trait_value") %>%
  pivot_longer(a:T_opt, names_to="tpc_param", values_to="tpc_param_value") %>%
  mutate(trait_tpc=paste(trait, tpc_param, sep="_"))

stats_lm_2<-ind_reg %>%
  group_by(trait, tpc_param) %>% 
  do(fit = lm(tpc_param_value ~ trait_value, data = ., na.action=na.omit)) %>% 
  rowwise() %>%
  mutate(intercept=summary(fit)$coefficients[1,1], 
         intercept_se=summary(fit)$coefficients[1,2], 
         intercept_p_val=summary(fit)$coefficients[1,4], 
         slope=summary(fit)$coefficients[2,1], 
         slope_se=summary(fit)$coefficients[2,2], 
         slope_p_val=summary(fit)$coefficients[2,4],
         adj_R2=summary(fit)$adj.r.squared, 
         model_p_val=p_val_func(fit),
         sig=ifelse(model_p_val<=0.05, 1, 0)) %>%
  ungroup() %>%
  dplyr::select(-fit)



## multiple regressions

true_param_list<-c("a", "E_a", "E_d", "Th")

out<-list()
for(i in true_param_list){
  out$true_param_list[[i]]<-summary(lm(eval(parse(text=i)) ~ Biovolume_Spheroid_GEO + Geodesic_Aspect_Ratio + Sigma.Intensity, data = traits_trueTPC, na.action=na.omit))
}
out


fit<-lm(K_22 ~ Biovolume_Spheroid_GEO, data=filter(demo_pheno_df, Moment=="Mean"), na.action=na.omit)
summary(fit)


summary(lm(K_22 ~ Biovolume_Spheroid_GEO + Geodesic_Aspect_Ratio + Sigma.Intensity, data=filter(demo_pheno_df, Moment=="Mean"), na.action=na.omit))


## trait-based TPC & K estimates

a_func<-function(size, shape, contents){1.99549 - 0.06188*size + 0.12502*shape + 0.5266*contents}
Ea_func<-function(contents){-0.6629233 + 0.48503528*contents}
Ed_func<-function(size){11.6972910 - 1.61148575*size}
Th_func<-function(size, shape, contents){355.996 + 3.881*size + 6.362*shape - 39.365*contents}
K22_func<-function(size, shape, contents){8.2847 - 0.4703*size + 0.15*shape - 2.6011*contents}

est_TPCs<-traits_trueTPC %>%
  mutate(est_a=a_func(Biovolume_Spheroid_GEO, Geodesic_Aspect_Ratio, Sigma.Intensity),
         est_Ea=Ea_func(Sigma.Intensity),
         est_Ed=Ed_func(Biovolume_Spheroid_GEO),
         est_Th=Th_func(Biovolume_Spheroid_GEO, Geodesic_Aspect_Ratio, Sigma.Intensity),
         est_K22=10^K22_func(Biovolume_Spheroid_GEO, Geodesic_Aspect_Ratio, Sigma.Intensity))

# export predicted TPC parameters
fwrite(est_TPCs, paste0(data_dir, "data/TPC_model_parameters_est.csv"))




# plot empirical and predicted TPCs (individual species)

TPC_predicted_2<-expand.grid(Species=spp_order, Temperature=seq(0, 42, length.out=1000)) %>%
  left_join(est_TPCs) %>%
  mutate(r=TPC_eqn(est_a, est_Ea, est_Ed, est_Th, Temperature, 10)) %>%
  left_join(range_r) %>%
  filter(r>=min_r)


svg(file=paste0(data_dir, "graphics/raw_graphics/TPCs_empirical_and_predicted.svg"), bg=NA, width=13, height=6.5)
ggplot() +
  geom_hline(yintercept=0, color="gray", size=0.7, linetype=2) +
  geom_line(data=TPC_predicted_trunc, aes(Temperature, r, color=Species)) +
  geom_line(data=TPC_predicted_2, aes(Temperature, r, color=Species), linetype=2) +
  geom_segment(data=TPC_summary_spread, aes(x=T_opt, y=0, xend=T_opt, yend=r_peak), color="black") +
  geom_segment(data=filter(TPC_summary, param %in% c("CT_min", "CT_max")), aes(x=param_val, y=0, xend=param_val, yend=min_r), color="black") +
  geom_point(data=TPC_df, aes(Temperature, r, color=Species)) +
  geom_point(data=filter(TPC_summary, param %in% c("CT_min", "CT_max", "T_opt")), aes(param_val, 0, shape=param), size=2) +
  geom_point(data=demo_gdata, aes(Temperature, param_val), shape=5, color="black", size=2, stroke=.8) +
  geom_text(data=left_join(TPC_summary_spread, range_r), aes(x=T_opt, y=-(range_r*0.15), label=round(T_opt, 1)), size=4)+
  geom_text(data=filter(TPC_summary, param %in% c("CT_min", "CT_max")), aes(x=param_val, y=min_r-(range_r*0.1), label=round(param_val, 1)), size=4) +
  scale_shape_manual(values=c(1, 1, 19)) +
  scale_color_manual(values=spp_color_palette) +
  scale_x_continuous(expand=expansion(mult = c(.1, .1))) +
  scale_y_continuous(expand=expansion(mult = c(.1, .1))) +
  labs(x="Temperature (C)", y=Intrinsic~growth~rate~(r)~(cells~cell^-1~d^-1)) +
  facet_wrap(~factor(Species), scales="free", nrow=3, labeller=labeller(Species=spp_labels)) +
  theme(panel.background=element_blank(), panel.border=element_rect(color = "black", fill=NA),
        plot.background=element_blank(), panel.grid=element_blank(),
        axis.text=element_text(size=14), axis.title=element_text(size=16),
        strip.text=element_text(size=13, face="italic"),
        strip.background=element_blank(), 
        legend.position="none",
        legend.text=element_text(size=12, face="italic"), legend.title=element_text(size=14))
dev.off()




#
### Demographic parameters vs. traits dataframe -----


## demo-pheno dataframe

demo_pheno_df_full<-demo_params_reduced %>%
  mutate(param=paste(param, Temperature, sep="_")) %>%
  bind_rows(T_diff_df) %>%
  filter(sig=="*") %>%
  mutate(param_val=ifelse(Species=="Colpidium striatum" & param %in% c("r_22", "r_25", "r_diff", "K_22", "K_25", "K_diff"), NA, param_val),
         log10_param_val=ifelse(Species=="Colpidium striatum" & param %in% c("r_22", "r_25", "r_diff", "K_22", "K_25", "K_diff"), NA, log10_param_val)) %>%
  bind_rows(TPC_summary) %>%
  mutate(param_val=log10_param_val) %>%  
  dplyr::select(Species, param, param_val) %>%
  spread(param, param_val) %>% 
  left_join(phenotype_summary[c("Species", "Moment", unique(phenotype_df$Trait))]) %>%
  ungroup()


trait_list<-c("Biovolume_Spheroid_GEO", "Geodesic_Aspect_Ratio", "Sigma.Intensity") 
param_list<-c("r_22", "r_25", "r_diff", "K_22", "K_25", "K_diff", "E_a", "E_d", "r_peak", "CT_max", "CT_min", "T_opt", "T_range", "TPC_asymmetry")
param_trait_list<-c(param_list, trait_list)

demo_pheno_df<-demo_pheno_df_full %>%
  dplyr::select(Species, Moment, all_of(param_list), all_of(trait_list))


#
### demo pheno correlation analysis -----


corr_data_1<-demo_pheno_df %>%
  dplyr::select(-Species) %>%
  group_by(Moment) %>%
  do(corr_df=replace_lower_triangle(cor_mat(dplyr::select(., -Moment)), by=NA),
     pval_df=replace_lower_triangle(cor_pmat(dplyr::select(., -Moment)), by=NA))

corr_df<-corr_data_1 %>%
  group_by(Moment) %>%
  do(data.frame(.$corr_df)) %>%
  rename(V1=rowname) %>%
  gather(V2, corr, -Moment, -V1)

pval_df<-corr_data_1 %>%
  group_by(Moment) %>%
  do(data.frame(.$pval_df)) %>%
  rename(V1=rowname) %>%
  gather(V2, p_val, -Moment, -V1)

corr_data<-left_join(corr_df, pval_df) %>% 
  drop_na(corr) %>%
  mutate(corr_sig=ifelse(p_val<=0.05, corr, NA),
         V1=factor(V1, levels=param_trait_list),
         V2=factor(V2, levels=rev(param_trait_list)),
         Moment=factor(Moment, levels=c("Mean", "Variance", "Skewness", "Kurtosis")))

corr_data_reduced<-corr_data %>%
  filter(V1 %in% param_list) %>%
  filter(V2 %in% trait_list) %>%
  mutate(V1=factor(V1, levels=param_list),
         V2=factor(V2, levels=rev(trait_list)))



# plot correlations summary table (demo params VS traits)

temp_moments<-c("Mean")

label_func<-function(x){
  case_when(x=="Biovolume_Spheroid_GEO" ~ "Volume", 
            x=="Geodesic_Aspect_Ratio" ~ "Aspect~ratio",
            x=="Sigma.Intensity" ~ "Contrast",
            T ~ paste0("'", word(x, 1, sep="_"), "'", "[", word(x, 2, sep="_"), "]"))
}

gdata<-corr_data %>% # corr_data OR corr_data_reduced   ## choose which correlation data to plot (corr_data = full ; corr_data_reduced = only trait-demo param correlations)
  filter(Moment %in% temp_moments) %>%
  mutate(corr_text=ifelse(corr==1, NA, sprintf("%.2f", round(corr, 2))),
         p_val_text=ifelse(corr==1, NA, ifelse(p_val<0.001, "(<0.001)", paste0("(", sprintf("%.3f", round(p_val, 3)), ")"))),
         g_text=paste0(corr_text, "\n", p_val_text), 
         x=label_func(V1),
         y=label_func(V2))




svg(file=paste0(data_dir, "graphics/raw_graphics/moment_corr_all.svg"), bg=NA, width=20, height=10)
ggplot() +
  geom_tile(data=gdata, aes(x=V1, y=V2, fill=corr), color="black", size=.25) +
  geom_text(data=filter(gdata, p_val<=0.05), aes(V1, V2, label=g_text), size=3, color="white") +
  scale_fill_gradient2(low="red", mid="white", high="blue", na.value = "white", midpoint=0, limits=c(-1,1), name="Correlation") + 
  geom_point(data=filter(gdata, p_val>0.05), aes(x=V1, y=V2), shape=4, color="black", size=3) +
  scale_x_discrete(labels=parse(text=unique(gdata$x))) +
  scale_y_discrete(labels=parse(text=rev(unique(gdata$y)))) +
  facet_rep_wrap(~Moment, ncol=2, repeat.tick.labels=T) + 
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=14), axis.text.y=element_text(size=14),
        strip.text=element_text(size=16),
        # strip.text=element_blank(),
        strip.background=element_blank(),
        axis.title=element_blank(),
        aspect.ratio=.7,
        legend.key.height=unit(1,"cm"), legend.key.width=unit(0.3,"cm"), legend.text=element_text(size=10), legend.title=element_text(size=12),
        panel.background=element_blank(), panel.border=element_blank(), plot.background=element_blank(), panel.grid=element_blank())
dev.off()




### param-trait linear regressions -----


# pairwise regressions
lm_demo_pheno_df<-demo_pheno_df %>%
  gather(V1, V1_val, -Species, -Moment) %>%
  left_join(demo_pheno_df) %>%
  gather(V2, V2_val, -Species, -Moment, -V1, -V1_val)


stats_lm_1<-lm_demo_pheno_df %>%
  group_by(Moment, V1, V2) %>% 
  do(fit = lm(V2_val ~ V1_val, data = ., na.action=na.omit)) %>% 
  rowwise() %>%
  mutate(intercept=summary(fit)$coefficients[1,1], 
         slope=summary(fit)$coefficients[2,1], 
         adj_R2=summary(fit)$adj.r.squared, 
         p_val=p_val_func(fit),
         slope_sig=ifelse(p_val<=0.05, slope, NA),
         sig=ifelse(p_val<=0.05, 1, 0)) %>%
  ungroup() %>%
  dplyr::select(-fit) %>%
  mutate(V1=factor(V1, levels=param_trait_list),
         V2=factor(V2, levels=rev(param_trait_list)),
         Moment=factor(Moment, levels=c("Mean", "Variance", "Skewness", "Kurtosis")))



#
### plot params VS traits -----


temp_moment<-"Mean"

stats_gdata<-stats_lm_1 %>%
  rename(param=V2, trait=V1) %>%
  filter(Moment==temp_moment, param %in% param_list, trait %in% trait_list) %>%
  mutate(param_trait=paste(param, trait, sep=" | ")) %>%
  mutate(g_text=paste0('R\U00B2 = ', round(adj_R2, 2),  ifelse(round(p_val, 2)==0, '\np < ', '\np = '), ifelse(round(p_val, 2)==0, "0.01", round(p_val, 2)), '\nm = ', round(slope, 2)),
         g_text_2=ifelse(word(param, 2, sep="_") %in% c("22", "25"), paste0("'m'", "[", word(param, 2, sep="_"),"]~", "' = '~", "'", sprintf("%.2f", round(slope, 2)), "'"), paste0("'m = '~", "'", sprintf("%.2f", round(slope, 2)), "'")))


gdata<-demo_pheno_df %>%
  gather(param, param_val, all_of(param_list)) %>%
  gather(trait, trait_val, all_of(trait_list)) %>%
  left_join(spp_name_df) %>%
  mutate(param_trait=paste(param, trait, sep=" | "),
         param=factor(param, levels=param_list),
         trait=factor(trait, levels=trait_list),
         Species_abbr=factor(Species_abbr, levels=spp_order_abbr)) %>%
  arrange(factor(trait, levels=trait_list), factor(param, levels=param_list)) %>%
  left_join(dplyr::select(stats_gdata, Moment, trait, param, sig)) %>%
  filter(Moment==temp_moment, sig==1) 



temp_pt_df<-data.frame(
  param=c("r_2", "r_2", "K_2", "K_diff", "r_diff",
          "r_peak", "r_peak", "r_peak", "CT_min", "CT_min", "E_a", "T_range"),
  trait=c("Biovolume_Spheroid_GEO", "Geodesic_Aspect_Ratio", "Biovolume_Spheroid_GEO", "Biovolume_Spheroid_GEO", "Sigma.Intensity",
          "Biovolume_Spheroid_GEO", "Geodesic_Aspect_Ratio", "Sigma.Intensity", "Biovolume_Spheroid_GEO", "Geodesic_Aspect_Ratio", "Sigma.Intensity", "Geodesic_Aspect_Ratio")
)


gg_func<-function(i){
  i<-i
  temp_param<-temp_pt_df[i,]$param
  temp_trait<-temp_pt_df[i,]$trait
  
  temp_gdata<-filter(gdata, str_detect(param, temp_param), trait==temp_trait) %>%
    mutate(param_val=10^param_val, trait_val=10^trait_val)
  if(temp_param %in% c("r_2", "K_2")) {
    temp_p1<-unique(temp_gdata$param)[[1]]
    temp_p2<-unique(temp_gdata$param)[[2]]
  }
  temp_stats_gdata<-filter(stats_gdata, str_detect(param, temp_param), trait==temp_trait)
  
  x_breaks<-if(str_detect(temp_trait, "Biovolume")) {c(10^4, 10^5)} else {waiver()} 
  x_ticks<-if(str_detect(temp_trait, "Biovolume")) {trans_format("log10", math_format(10^.x))} else {waiver()} 
  x_lab<-case_when(temp_trait=="Biovolume_Spheroid_GEO" ~ "Volume~(mu*'m'^3)", 
                   temp_trait=="Geodesic_Aspect_Ratio" ~ "Aspect~ratio", 
                   temp_trait=="Sigma.Intensity" ~ "Contrast")
  y_lab=case_when(temp_param %in% c("r_2", "K_2") ~ paste0("'", word(temp_param, 1, sep="_"), "'"),
                  temp_param %in% c("r_diff", "K_diff") ~ paste0("'|'*", "'", word(temp_param, 1, sep="_"), "'", "[", word(temp_param, 2, sep="_"), "]*", "'|'"),
                  T ~ paste0("'", word(temp_param, 1, sep="_"), "'", "[", word(temp_param, 2, sep="_"), "]"))

  stripe_range=10^seq(log10(min(temp_gdata$trait_val)), log10(max(temp_gdata$trait_val)), length.out=16)
  
  stat_lab_x<-c(min(temp_gdata$trait_val, na.rm=T), max(temp_gdata$trait_val, na.rm=T))
  stat_lab_y<-c((max(temp_gdata$param_val, na.rm=T)-min(temp_gdata$param_val, na.rm=T))*0.1)
  
  g_base<-
    ggplot() +
    scale_shape_manual(values=c(19,1), guide=F) +
    scale_linetype_manual(values=c(1,2), guide=F) +
    scale_color_manual(values=spp_color_palette, name="Species")+
    labs(x=parse(text=x_lab), y=parse(text=y_lab)) +
    scale_y_log10() +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    theme(panel.background = element_blank(), panel.border = element_rect(color="black", fill=NA),
          plot.background=element_blank(), panel.grid=element_blank(),
          axis.text.x=element_text(size=14, vjust=1), axis.text.y=element_text(size=14, hjust=1), 
          axis.title=element_text(size=16, margin = margin(0, 0, 0, 0)),
          strip.text=element_text(size=12),
          legend.text=element_text(size=8, face="italic"), legend.title=element_text(size=12, margin=margin(0, 0, 0, 0)),
          legend.key=element_blank(),
          legend.key.height=unit(4, "mm"), legend.key.width=unit(1, "mm"),
          legend.position="none",
          plot.margin = unit(c(0.0, 0.0, 0.0, 0.0), "cm"),
          aspect.ratio=1)
  
  
  g_temp<-if(temp_param %in% c("r_2", "K_2")){
    g_base +
      geom_smooth(data=filter(temp_gdata, param==temp_p2), aes(trait_val, param_val), color=NA, fill='gray', alpha=0.25, method="lm", formula="y~x", se=T) +
      geom_vline(xintercept = stripe_range, color="white", size=1) +
      geom_smooth(data=filter(temp_gdata, param==temp_p1), aes(trait_val, param_val), color=NA, fill='gray', alpha=0.25, method="lm", formula="y~x", se=T) +
      geom_text(data=temp_stats_gdata[2,], aes(x=ifelse(slope>=0, stat_lab_x[[1]], stat_lab_x[[2]]), y=Inf, label=g_text_2, hjust=ifelse(slope>=0, -0.08, 1.08)), vjust=3, parse=T, size=4) +
      scale_x_log10(breaks=x_breaks, labels=x_ticks) +
      annotation_logticks(sides="bl", short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))
  } else if(!temp_param %in% c("r_2", "K_2")){
    g_base +
      geom_smooth(data=temp_gdata, aes(trait_val, param_val), color=NA, fill='gray', alpha=0.25, method="lm", formula="y~x", se=T) +
      scale_x_log10(breaks=x_breaks, labels=x_ticks) +
      annotation_logticks(sides="bl", short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))
  }
  
  return(
         g_temp +
           geom_smooth(data=temp_gdata, aes(trait_val, param_val, linetype=param), color='black', size=0.7, method="lm", formula="y~x", se=F, fullrange=T) +
           geom_point(data=temp_gdata, aes(trait_val, param_val, color=Species_abbr, shape=param), size=2) +
           geom_text(data=temp_stats_gdata[1,], aes(x=ifelse(slope>=0, stat_lab_x[[1]], stat_lab_x[[2]]), y=Inf, label=g_text_2, hjust=ifelse(slope>=0, -0.08, 1.08)), vjust=1.5, parse=T, size=4)
  )
}

gg_assign<-function(i){assign(paste("g", i, sep="_"), gg_func(i))}

for(i in 1:nrow(temp_pt_df)){assign(paste("g", i, sep="_"), gg_func(i))}



gg_blank<-ggplot() + theme_void()

svg(file=paste0(data_dir, "graphics/raw_graphics/param_trait_regressions.svg"), bg=NA, width=17, height=7.5) 
ggarrange(g_1, g_2, gg_blank, g_6, g_7, g_8,
          g_3, gg_blank, gg_blank, g_9, g_10, g_11,
          g_4, g_5, gg_blank, gg_blank, g_12, gg_blank,
          nrow=3, ncol=6, align="hv")
dev.off()



#
### param-trait linear model w/ commonality analysis -----


lm_df<-demo_pheno_df %>%
  group_by(Moment) %>%
  mutate(across(c(trait_list), cent_rescale)) %>% 
  ungroup() %>%
  gather(param, param_val, all_of(param_list)) %>%
  filter(Moment=="Mean", param %in% c("r_peak", "CT_min", "r_22", "r_25")) %>%
  split(.$param) %>%
  map(~{
    temp_p<-unique(.x$param)
    temp_d<-.x
    lm_fit<-if(temp_p=="r_peak"){
      lm(param_val ~ Biovolume_Spheroid_GEO + Geodesic_Aspect_Ratio + Sigma.Intensity, data=temp_d, na.action=na.omit)
    } else {
      lm(param_val ~ Biovolume_Spheroid_GEO + Geodesic_Aspect_Ratio, data=temp_d, na.action=na.omit)
      }
    setNames(list(vif(lm_fit), regr(lm_fit)), c("VIF", "CA"))
    }) 


pval_df<-lm_df %>%
  map(~{
    temp_fstat<-.x$CA$LM_Output$fstatistic
    pf(temp_fstat[1], temp_fstat[2], temp_fstat[3], lower.tail=F)
  }) %>%
  map_dfr(~data.frame(model_pval=., row.names=NULL),
          .id="Param")

lm_summary<-lm_df %>%
  map_dfr(~data.frame(coef(.$CA$LM_Output),
                      R2=.$CA$LM_Output$r.squared),
          .id="Param") %>%
  rownames_to_column(var="Trait") %>%
  mutate(Trait=word(Trait, 1, sep=fixed(".."))) %>%
  filter(!Trait=="(Intercept)") %>%
  left_join(pval_df)

ca_summary<-lm_df %>%
  map_dfr(~data.frame(VIF=.$VIF,
                      .$CA$Commonality_Data$CCTotalbyVar),
          .id="Param") %>%
  rownames_to_column(var="Trait") %>%
  mutate(Trait=word(Trait, 1, sep=fixed("..")))



lm_summary_df<-left_join(lm_summary, ca_summary) %>%
  mutate('U%'=Unique/R2*100,
         Param=factor(Param, levels=param_list),
         Trait=factor(Trait, levels=trait_list),
         across(where(is.numeric), ~sprintf("%.3f", round(.x, 3)))) %>%
  select(Param, R2, model_pval, Trait, Estimate, Pr...t.., VIF, Unique, 'U%', Common, Total) %>%
  arrange(Param, Trait)

fwrite(lm_summary_df, paste0(data_dir, "data/param_trait_commonality_summary.csv"))




#
### param-trait bootstrapping analysis -----

temp_moment<-"Mean"

stats_gdata<-stats_lm_1 %>%
  rename(param=V2, trait=V1) %>%
  filter(Moment==temp_moment, param %in% param_list, trait %in% trait_list) %>%
  mutate(param_trait=paste(param, trait, sep=" | ")) %>%
  mutate(g_text=paste0('R\U00B2 = ', round(adj_R2, 2),  ifelse(round(p_val, 2)==0, '\np < ', '\np = '), ifelse(round(p_val, 2)==0, "0.01", round(p_val, 2)), '\nm = ', round(slope, 2)),
         g_text_2=ifelse(word(param, 2, sep="_") %in% c("22", "25"), paste0("'s'", "[", word(param, 2, sep="_"),"]~", "' = '~", "'", sprintf("%.2f", round(slope, 2)), "'"), paste0("'s = '~", "'", sprintf("%.2f", round(slope, 2)), "'")))


gdata<-demo_pheno_df %>%
  gather(param, param_val, all_of(param_list)) %>%
  gather(trait, trait_val, all_of(trait_list)) %>%
  left_join(spp_name_df) %>%
  mutate(param_trait=paste(param, trait, sep=" | "),
         param=factor(param, levels=param_list),
         trait=factor(trait, levels=trait_list),
         Species_abbr=factor(Species_abbr, levels=spp_order_abbr)) %>%
  arrange(factor(trait, levels=trait_list), factor(param, levels=param_list)) %>%
  left_join(dplyr::select(stats_gdata, Moment, trait, param, sig)) %>%
  filter(Moment==temp_moment, sig==1)  # , sig==1
head(gdata)



temp_pt_df<-data.frame(
  param=c("r_peak", "r_peak", "r_peak", "r_22", "r_25", "r_22", "r_25",
          "CT_min", "CT_min", "E_a", "K_22", "K_25", 
          "T_range", "K_diff", "r_diff"),
  trait=c("Biovolume_Spheroid_GEO", "Geodesic_Aspect_Ratio", "Sigma.Intensity", "Biovolume_Spheroid_GEO", "Biovolume_Spheroid_GEO", "Geodesic_Aspect_Ratio", "Geodesic_Aspect_Ratio",
          "Biovolume_Spheroid_GEO", "Geodesic_Aspect_Ratio", "Sigma.Intensity", "Biovolume_Spheroid_GEO", "Biovolume_Spheroid_GEO", 
          "Geodesic_Aspect_Ratio", "Biovolume_Spheroid_GEO", "Sigma.Intensity")
)

k_stat <- function(df, i) {
  resample=df[i,]
  fit<-lm(V2_val ~ V1_val, data = resample, na.action=na.omit)
  return(c(fit$coef[[1]], fit$coef[[2]]))
}

gg_func<-function(i){
  temp_param<-temp_pt_df[i,]$param
  temp_trait<-temp_pt_df[i,]$trait
  
  temp_bs_data<-lm_demo_pheno_df %>%
    filter(Moment=="Mean", V2==temp_param, V1==temp_trait)
  
  temp_b<-boot(temp_bs_data, k_stat, R=1000)
  
  temp_gdata_1<-data.frame(temp_b$t) %>% 
    rename(intercept=X1, slope=X2) %>%
    mutate(sign_test=ifelse(sign(slope)==sign(temp_b$t0[2]), 1, 0))
  
  pct_true_sign<-sum(temp_gdata_1$sign_test)/nrow(temp_gdata_1)*100
  
  temp_gdata_2<-filter(gdata, str_detect(param, temp_param), trait==temp_trait) %>%
    mutate(param_val=10^param_val, trait_val=10^trait_val)
  
  temp_gdata_3<-data.frame(slope=temp_b$t0[2],
                           pct_true=sum(temp_gdata_1$sign_test)/nrow(temp_gdata_1)*100) %>%
    mutate(g_text=paste0(pct_true, '%'))
  
  stat_lab_x<-c(min(temp_gdata_2$trait_val, na.rm=T), max(temp_gdata_2$trait_val, na.rm=T))
  
  x_breaks<-if(str_detect(temp_trait, "Biovolume")) {c(10^4, 10^5)} else {waiver()} 
  x_ticks<-if(str_detect(temp_trait, "Biovolume")) {trans_format("log10", math_format(10^.x))} else {waiver()} 
  x_lab<-case_when(temp_trait=="Biovolume_Spheroid_GEO" ~ "Volume~(mu*'m'^3)", 
                   temp_trait=="Geodesic_Aspect_Ratio" ~ "Aspect~ratio", 
                   temp_trait=="Sigma.Intensity" ~ "Contrast")
  y_lab=case_when(
    temp_param %in% c("r_diff", "K_diff") ~ paste0("'|'*", "'", word(temp_param, 1, sep="_"), "'", "[", word(temp_param, 2, sep="_"), "]*", "'|'"),
    T ~ paste0("'", word(temp_param, 1, sep="_"), "'", "[", word(temp_param, 2, sep="_"), "]"))
  
  return(
    ggplot() +
      geom_smooth(data=temp_gdata_2, aes(trait_val, param_val), color=NA, fill=NA, alpha=0.25, method="lm", formula="y~x", se=T) +
      geom_abline(data=temp_gdata_1, aes(intercept=intercept, slope=slope), color="gray", size=0.1) +
      geom_abline(intercept=temp_b$t0[1], slope=temp_b$t0[2], color='black', size=0.7) +
      geom_point(data=temp_gdata_2, aes(trait_val, param_val, color=Species_abbr, shape=param), size=2) +
      geom_text(data=temp_gdata_3, aes(x=ifelse(slope>=0, stat_lab_x[[1]], stat_lab_x[[2]]), y=Inf, label=g_text, hjust=ifelse(slope>=0, -0.08, 1.08)), vjust=3, parse=F, size=4) +
      
      scale_shape_manual(values=c(19,1), guide=F) +
      scale_linetype_manual(values=c(1,2), guide=F) +
      scale_color_manual(values=spp_color_palette, name="Species")+
      labs(x=parse(text=x_lab), y=parse(text=y_lab)) +
      scale_x_log10(breaks=x_breaks, labels=x_ticks) +
      scale_y_log10() +
      annotation_logticks(sides="bl", short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm")) +
      guides(color = guide_legend(override.aes = list(size = 3))) +
      theme(panel.background = element_blank(), panel.border = element_rect(color="black", fill=NA),
            plot.background=element_blank(), panel.grid=element_blank(),
            axis.text.x=element_text(size=14, vjust=1), axis.text.y=element_text(size=14, hjust=1), 
            axis.title=element_text(size=16, margin = margin(0, 0, 0, 0)),
            strip.text=element_text(size=12),
            legend.text=element_text(size=8, face="italic"), legend.title=element_text(size=12, margin=margin(0, 0, 0, 0)),
            legend.key=element_blank(),
            legend.key.height=unit(4, "mm"), legend.key.width=unit(1, "mm"),
            legend.position="none",
            plot.margin = unit(c(0.0, 0.0, 0.0, 0.0), "cm"),
            aspect.ratio=1)
  )
}
gg_assign<-function(i){assign(paste("g", i, sep="_"), gg_func(i))}
for(i in 1:nrow(temp_pt_df)){assign(paste("g", i, sep="_"), gg_func(i))}



gg_blank<-ggplot() + theme_void()

svg(file=paste0(data_dir, "graphics/raw_graphics/param_trait_regressions_boot.svg"), bg=NA, width=20, height=7.5)
ggarrange(g_1, g_2, g_3, g_4, g_5, g_6, g_7, 
          g_8, g_9, g_10, g_11, g_12, gg_blank, gg_blank,
          g_13, g_14, g_15, gg_blank, gg_blank, gg_blank, gg_blank,
          nrow=3, ncol=7, align="hv") 
dev.off()


#
### model results -----


pheno_data<-fread(paste0(data_dir, "data/traits_summary.csv")) %>%
  mutate(mass_g=Biovolume_Spheroid_GEO*10^-12,
         metabolic_rate=0.001520711*mass_g^0.97)



## plot model results for Figures 5/S4 A&E-G and Figure S5

model_rawdata<-fread(paste0(data_dir, "data/model_rawdata_small.csv"))
total_K<-3324.34 # For empirical TPCs use 3324.34 | For trait-based TPCs use 2054.25

density_gdata<-model_rawdata %>%
  gather(Species, density, -c(n, Temperature, alpha, delta)) %>%
  mutate(Species=factor(Species, levels=rev(spp_order)),
         density=density*total_K,
         alpha_label=paste0("alpha~'='~", alpha),
         delta_label=paste0("delta~'='~", delta),
         facet_label=paste0(alpha_label, "~'|'~", delta_label)) %>%
  arrange(delta, rev(alpha))
density_gdata$delta_label=factor(density_gdata$delta_label, levels=rev(unique(density_gdata$delta_label)))
density_gdata$facet_label=factor(density_gdata$facet_label, levels=rev(unique(density_gdata$facet_label)))

max_y<-max(c(10^3, density_gdata$density))

## plot model results (constant alphas)

# equilibrium species densities
svg(file=paste0(data_dir, "graphics/raw_graphics/model_densities.svg"), width=7, height=7, bg=NA) # Figures 5A/S4A: width=7, height=7    Figure S5: width=9, height=5
ggplot()+
  geom_vline(xintercept=c(10, 16, 22, 28, 34), color="gray", size=0.5) +  # Figures 5A & S4A only
  geom_line(data=density_gdata, aes(Temperature, density, color=Species), size=0.7) + 
  scale_color_manual(values=rev(spp_color_palette)) +
  scale_y_log10(labels=trans_format("log10", math_format(10^.x)), limits=c(NA, max_y), breaks=c(10^-1, 10^0, 10^1, 10^2, 10^3)) +
  coord_capped_cart(left='top', gap=0) +
  facet_wrap(~delta_label, ncol=1, strip.position="top", labeller=label_parsed) + # Figures 5A & S4A only
  # facet_wrap(~facet_label, ncol=3, strip.position="top", labeller=label_parsed) +  # Figure S5 only
  labs(x="Temperature (C)", y="Equilibrium density"~(ind~mL^-1)) +
  theme(panel.background = element_blank(), panel.border=element_blank(),
        panel.spacing=unit(2.5, "cm"),  # Figures 5A & S4A only
        legend.position="none",
        axis.line=element_line(),
        axis.text=element_text(size=11), axis.title=element_text(size=13),
        axis.text.y=element_text(size=11),
        strip.text.x=element_text(size=13, angle=0, hjust=0.5), 
        strip.background=element_blank(), plot.background=element_blank(), panel.grid=element_blank())
dev.off()




## plot model results for Figures 5/S4 B-D

model_rawdata<-fread(paste0(data_dir, "data/model_rawdata_large_empiricalTPCs.csv")) 
total_K<-3324.34 # For empirical TPCs use 3324.34 | For trait-based TPCs use 2054.25

model_spp_data<-model_rawdata %>%
  gather(Species, density, -c(Temperature, alpha, delta)) %>%
  left_join(pheno_data) %>%
  mutate(presence=sign(density),
         density=density*total_K,
         mass_present=ifelse(presence==0, NA, mass_g),
         biomass=density*mass_g,
         density_metabolic_rate=density*metabolic_rate, 
         biomass_metabolic_rate=biomass*metabolic_rate,
         Species=factor(Species, levels=rev(spp_order)))


model_comm_data<-model_spp_data %>%
  group_by(Temperature, alpha, delta) %>% 
  summarize(spp_richness=sum(presence), avg_mass=mean(mass_present, na.rm=T), comm_biomass=sum(biomass), 
            comm_metabolic_rate=sum(density_metabolic_rate), comm_metabolic_rate_2=sum(biomass_metabolic_rate)) %>% 
  ungroup


base_theme<-theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, color="black", size=1),
                  plot.background=element_blank(), panel.grid=element_blank(),
                  axis.text=element_text(size=11), 
                  legend.text=element_text(size=10),
                  legend.key.height=unit(0.75,"cm"), legend.key.width=unit(0.3,"cm"),
                  plot.margin = unit(c(0.1, 0.0, 0.7, 0.0), "cm"),
                  # legend.position="none",
                  plot.title=element_blank(),
                  axis.title=element_blank(),
                  legend.title=element_blank(),
                  aspect.ratio=1)



# plot equilibrium species richness (Figures 5/S4 B)
model_g1<-
  ggplot() +
  geom_raster(data=model_comm_data, aes(Temperature, delta, fill=spp_richness)) +
  scale_fill_viridis_c(option = "plasma", breaks=c(0, 7, 14)) +
  scale_x_continuous(expand=expansion(mult = c(0.0, 0.0))) +
  scale_y_continuous(expand=expansion(mult = c(0.0, 0.0)), limits=c(0, 1), breaks=c(0, 0.5, 1), labels=c("0.0", "0.5", "1.0")) +
  labs(x="Temperature (C)", y=Global~mortality~rate~(delta), title="Species richness") +
  base_theme


# plot average species mass (Figures 5/S4 C)
model_g2<-
  ggplot() +
  geom_raster(data=model_comm_data, aes(Temperature, delta, fill=avg_mass)) +
  scale_fill_viridis_c(option="plasma", trans="log10", na.value=plasma(1), limits=c(NA, 0.000000602), breaks=c(10^-8, 10^-7), labels=trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(expand=expansion(mult = c(0.0, 0.0))) +
  scale_y_continuous(expand=expansion(mult = c(0.0, 0.0)), limits=c(0, 1), breaks=c(0, 0.5, 1), labels=c("0.0", "0.5", "1.0")) +
  labs(x="Temperature (C)", y=Global~mortality~rate~(delta), title="Mean species mass (g)") +
  base_theme


# plot total community metabolic rate (Figures 5/S4 D)
model_g3<-
  ggplot() +
  geom_raster(data=model_comm_data, aes(Temperature, delta, fill=comm_metabolic_rate)) + 
  scale_fill_viridis_c(option="plasma", na.value=plasma(1), trans="log10", labels=trans_format("log10", math_format(10^.x)), limits=c(10^-9.5, max(model_comm_data$comm_metabolic_rate))) +
  scale_x_continuous(expand=expansion(mult = c(0.0, 0.0))) +
  scale_y_continuous(expand=expansion(mult = c(0.0, 0.0)), limits=c(0, 1), breaks=c(0, 0.5, 1), labels=c("0.0", "0.5", "1.0")) +
  labs(x="Temperature (C)", y=Global~mortality~rate~(delta), title="Respiration rate"~(W~mL^-1)) + 
  base_theme


svg(file=paste0(data_dir, "graphics/raw_graphics/model_comm_summary.svg"), bg=NA, width=4, height=7.9)
ggarrange(model_g1, model_g2, model_g3,
          nrow=3, ncol=1, align="hv")
dev.off()




## plot model results for Figure S6 (random alphas)

model_rawdata<-fread(paste0(data_dir, "data/model_rawdata_short.csv"))
total_K<-3324.34 # For empirical TPCs use 3324.34 | For trait-based TPCs use 2054.25

density_gdata<-model_rawdata %>%
  gather(Species, density, -c(n, Temperature, alpha, delta)) %>%
  mutate(Species=factor(Species, levels=rev(spp_order)),
         density=density*total_K,
         alpha_label=paste0("alpha~'='~", alpha),
         delta_label=paste0("delta~'='~", delta),
         facet_label=paste0(alpha_label, "~'|'~", delta_label)) %>%
  arrange(delta, rev(alpha))
density_gdata$facet_label=factor(density_gdata$facet_label, levels=rev(unique(density_gdata$facet_label)))

density_gdata_2<-density_gdata %>%
  group_by(Temperature, alpha, delta, Species, delta_label) %>%
  summarize(lower=quantile(density, 0.05), upper=quantile(density, 0.95)) %>%
  ungroup()



# equilibrium species densities
svg(file=paste0(data_dir, "graphics/raw_graphics/model_densities.svg"), width=7, height=3, bg=NA)
ggplot()+
  geom_ribbon(data=density_gdata_2, aes(Temperature, ymin=lower, ymax=upper, color=Species, fill=Species), size=0.4, alpha=0.8) + 
  scale_color_manual(values=rev(spp_color_palette)) +
  scale_fill_manual(values=rev(spp_color_palette)) +
  scale_y_log10(labels=trans_format("log10", math_format(10^.x))) +
  coord_capped_cart(left='top', gap=0) +
  facet_wrap(~rev(delta_label), ncol=1, strip.position="top", labeller=label_parsed) +
  labs(x="Temperature (C)", y="Equilibrium density"~(ind~mL^-1)) +
  theme(panel.background = element_blank(), panel.border=element_blank(),
        legend.position="none",
        axis.line=element_line(),
        axis.text=element_text(size=11), axis.title=element_text(size=13),
        axis.text.y=element_text(size=11),
        strip.text.x=element_blank(), 
        strip.background=element_blank(), plot.background=element_blank(), panel.grid=element_blank())
dev.off()



#
### plot community experiment results -----


community_rawdata<-fread(paste0(data_dir, "data/CommunityExperiment_data.csv")) %>%
  rename(Tillina='Tillina magna') %>%
  mutate(Tillina=ifelse(Tillina==0.5, 1, Tillina)) %>%
  rename('Tillina magna'=Tillina) %>%
  filter(!Temperature==10) %>%
  rename(Respiration_Rate="Respiration_Rate_(%O2.min-1)") %>%
  mutate(Respiration_Rate=ifelse(Respiration_Rate>=mean(Respiration_Rate, na.rm=T)-2*sd(Respiration_Rate, na.rm=T) & Respiration_Rate<=mean(Respiration_Rate, na.rm=T)+2*sd(Respiration_Rate, na.rm=T), Respiration_Rate, NA),
         Respiration_Rate=Respiration_Rate/100/60*20)

pheno_data<-fread(paste0(data_dir, "data/traits_summary.csv")) %>%
  mutate(mass_g=Biovolume_Spheroid_GEO*10^-12,
         metabolic_rate=0.001520711*mass_g^0.97)
total_K<-3324.34 # For empirical TPCs use 3324.34 | For trait-based TPCs use 2054.25

community_df<-community_rawdata %>%
  select(ID:"Halteria grandinella") %>%
  pivot_longer(-c(ID, Temperature), names_to="Species", values_to="presence") %>%
  left_join(pheno_data) %>%
  filter(presence==1) %>%
  group_by(ID, Temperature) %>%
  summarize(spp_richness=sum(presence, na.rm=T), avg_mass=mean(mass_g, na.rm=T)) %>%
  ungroup() %>%
  left_join(select(community_rawdata, ID, Temperature, Respiration_Rate))


model_expTest_rawdata<-fread(paste0(data_dir, "data/model_rawdata_small.csv")) 


model_expTest_spp_data<-model_expTest_rawdata %>%
  select(-n) %>%
  gather(Species, density, -c(Temperature, alpha, delta)) %>%
  left_join(pheno_data) %>%
  mutate(presence=sign(density),
         density=density*total_K,
         mass_present=ifelse(presence==0, NA, mass_g),
         biomass=density*mass_g,
         density_metabolic_rate=density*metabolic_rate,
         biomass_metabolic_rate=biomass*metabolic_rate,
         Species=factor(Species, levels=rev(spp_order)))


model_expTest_comm_data<-model_expTest_spp_data %>%
  group_by(Temperature, alpha, delta) %>%
  summarize(spp_richness=sum(presence), avg_mass=mean(mass_present, na.rm=T), comm_biomass=sum(biomass),
            comm_metabolic_rate=sum(density_metabolic_rate), comm_metabolic_rate_2=sum(biomass_metabolic_rate)) %>%
  ungroup() %>%
  filter(alpha==0.01)




exp_g1<-
  ggplot() +
  geom_line(data=model_expTest_comm_data, aes(Temperature, spp_richness, linetype=factor(delta)), lineend="round", color="gray50")  +
  geom_smooth(data=community_df, aes(Temperature, spp_richness), color=NA, fill="gray70", span=0.8) +
  geom_point(data=community_df, aes(Temperature, spp_richness), color="gray30") +
  geom_smooth(data=community_df, aes(Temperature, spp_richness), size=1.5, color="gold1", fill=NA, span=0.8) +
  scale_linetype_manual(values=c(1,5,2,3), name=bquote(delta)) +
  scale_x_continuous(expand=expansion(mult = c(0.0, 0.0))) +
  scale_y_continuous(limits=c(0,10), breaks=c(0, 5, 10), labels=c("0", "5", "10")) +
  labs(y="Species richness") +
  base_theme


exp_g2<-
  ggplot() +
  geom_line(data=model_expTest_comm_data, aes(Temperature, avg_mass, linetype=factor(delta)), lineend="round", color="gray50")  +
  geom_smooth(data=community_df, aes(Temperature, avg_mass), color=NA, fill="gray70", span=0.8) +
  geom_point(data=community_df, aes(Temperature, avg_mass), color="gray30") +
  geom_smooth(data=community_df, aes(Temperature, avg_mass), size=1.5, color="gold1", fill=NA, span=0.8) +
  scale_linetype_manual(values=c(1,5,2,3), name=bquote(delta)) +
  scale_x_continuous(expand=expansion(mult = c(0.0, 0.0))) +
  scale_y_continuous(breaks=c(0, 2*10^-7, 4*10^-7, 6*10^-7), labels=c("0", "2", "4", "6")) +
  labs(y="Mean species mass (g)"~('x'*10^-7)) +
  base_theme


exp_g3<-
  ggplot() +
  geom_line(data=model_expTest_comm_data, aes(Temperature, comm_metabolic_rate*20, linetype=factor(delta)), lineend="round", color="gray50") +
  geom_smooth(data=community_df, aes(Temperature, Respiration_Rate), color=NA, fill="gray70", span=0.8) +
  geom_point(data=community_df, aes(Temperature, Respiration_Rate), color="gray30") +
  geom_smooth(data=community_df, aes(Temperature, Respiration_Rate), size=1.5, color="gold1", se=F, span=0.8) +
  scale_linetype_manual(values=c(1,5,2,3), name=bquote(delta)) +
  scale_x_continuous(expand=expansion(mult = c(0.0, 0.0))) +
  scale_y_continuous(breaks=c(0, 1*10^-6, 2*10^-6, 3*10^-6), labels=c("0", "1", "2", "3"),
                     sec.axis = sec_axis(~ ./20, breaks=c(0, 0.5*10^-7, 1*10^-7, 1.5*10^-7), labels=c("0", "0.5", "1.0", "1.5"), name="Respiration rate (theory)"~(W~mL^-1)~('x'*10^-7))) +  
  labs(x="Temperature (C)", y="Respiration rate (exp.)"~(W~mL^-1)~('x'*10^-6)) +
  base_theme
  


svg(file=paste0(data_dir, "graphics/raw_graphics/exp_comm_summary.svg"), bg=NA, width=4, height=7.9) 
ggarrange(exp_g1, exp_g2, exp_g3,
          nrow=3, ncol=1, align="hv") 
dev.off()




#


### THE END -----

