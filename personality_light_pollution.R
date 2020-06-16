#This document is an executable version of a Rmarkdown file. As such, discussions and plots cannot be visualised here. For more information, consult "Experiment_2_IIV.Rmd" or "Experiment_2_IIV.html". The numbers of interactions for the MCMCglmm module have been modified to reduce processing times. 


#Install missing packages and load all R libraries used in this analysis

if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "ggpubr", "MCMCglmm", "gridExtra", "psycho", 
               "lme4", "ggsignif", "brms", "lattice", "lmerTest", "cowplot",
               "tidybayes", "latticeExtra", "htmlTable", "magrittr", "DataCombine", 
               "ggplot2", "kableExtra")
##All done! 

library(tidyverse)
library(ggpubr)
library(MCMCglmm)
library(gridExtra)
library(patchwork) ##library added to include title with plots
library(psycho)
library(lme4)
library(ggsignif)
library(brms)
library(lattice)
library(lmerTest)
library(cowplot)
library(tidybayes)
library(latticeExtra)
library(kableExtra)
library(htmlTable)
library(magrittr)
library(DataCombine)
library(ggplot2)


### First lets check take a look in raw data

lightdark=read_csv("data/Experiment_2_data.csv", col_types = cols()) %>% 
  dplyr::select( id="Animal", srtime="Time colected", treatment="Treatment",
                 period = "Treatment order", obs = "Occasion", wt = "weight", 
                 shellwt = "shell weight", temp = "temperature", sr = "startle response", 
                 haemoabs = "haemocyanin",haemoporc= "haemocyaninporc", 
                 mr="oxygen consumption")  %>% 
  mutate(id, factor(id), srtime= factor(srtime), treatment = factor(treatment))

#Add a second column merging time and treatment for posterior use
lightdark = data.frame(lightdark %>% 
       unite(time_treatment, c(srtime, treatment), 
             sep = "_", remove = FALSE)) 


D1 <- lightdark %>% 
  group_by(treatment) %>% 
  ggplot(aes(x = as.factor(id), y = sr, color =treatment)) +
  geom_point(alpha = 0.7) + labs(x = "Individual ID", 
               y = "Startle response duration") + 
  guides(fill=FALSE) + 
  theme_classic() + 
  theme(legend.position = "none") +
  scale_color_manual(values=c("#172c3d", "orange")) 

D2 <- lightdark %>% 
  group_by(treatment) %>% 
  ggplot(aes(x = as.factor(id), y = log10(sr+ 1), color =treatment)) +
  geom_point(alpha = 0.7) + labs(x = "Individual ID", 
               y = "Startle response duration (log 10)") + 
  guides(fill=FALSE)  + theme_classic() + 
  scale_color_manual(values=c("#172c3d", "orange")) 

D1 + D2 + 
  plot_annotation(tag_levels = c('A','B'), 
    title ='The effect of permanent light on statle response duration', 
    caption = 'Plot of startle reponse duration, raw data (A) and log10 transformed (B) per individual.') 

  
# 1. The effect of permanent light on the startle response duration

### Model selection: comparing priors


##adding differnt priors to test the best model fit
          
prior1.1<-list(R = list(V = diag(1), nu = 1.002), 
               G = list(G1 = list(V = diag(1), nu =1.002, alpha.V=100),
      G2 = list(V = diag(1), nu =1.002, alpha.V=100)))
          
prior1.2<-list(R = list(V = diag(1), nu = 0.002), 
               G = list(G1 = list(V = diag(1), nu =0.002, alpha.V=1000),
      G2 = list(V = diag(1), nu =0.002, alpha.V=1000)))
          
prior1.3<-list(R = list(V = diag(1), nu = 1.002), 
               G = list(G1 = list(V = diag(1), nu =1.002, alpha.V=1000),
      G2 = list(V = diag(1), nu =1.002, alpha.V=1000)))

prior1.4<-list(R = list(V = diag(1), nu = 0.002), 
               G = list(G1 = list(V = diag(1), nu =0.002, alpha.V=100),
      G2 = list(V = diag(1), nu =0.002, alpha.V=1000)))
          
M1.1 <- MCMCglmm(log(sr+1) ~ srtime + treatment +  obs + 
 wt + haemoporc + treatment*srtime
                 , random=~ idh(id):units + period
                 , prior= prior1.1, family="gaussian"
                 , pl=TRUE, pr= TRUE
                 , saveX = TRUE , saveZ = TRUE
                 , data=lightdark, DIC= T
                 , singular.ok = TRUE
                 , nitt=5000, thin=10, burnin=50,verbose = F)  
# nitt=5000000, thin=100, burnin=500000 (original)
          
M1.2 <- MCMCglmm(log(sr+1) ~ srtime + treatment +  obs +  
 wt + haemoporc + treatment*srtime
                 , random=~ idh(id):units + period
                 , prior= prior1.2, family="gaussian"
                 , pl=TRUE, pr= TRUE
                 , saveX = TRUE , saveZ = TRUE
                 , data=lightdark, DIC= T
                 , singular.ok = TRUE
                 , nitt=5000, thin=10, burnin=50,verbose = F)  
# nitt=5000000, thin=100, burnin=500000

M1.3 <- MCMCglmm(log(sr+1) ~ srtime + treatment +  obs + 
 wt + haemoporc + treatment*srtime
                 , random=~ idh(id):units + period
                 , prior= prior1.3, family="gaussian" 
                 , pl=TRUE, pr= TRUE
                 , saveX = TRUE , saveZ = TRUE
                 , data=lightdark, DIC= T
                 , singular.ok = TRUE
                 , nitt=5000, thin=10, burnin=50,verbose = F)  
# nitt=5000000, thin=100, burnin=500000
          
M1.4 <- MCMCglmm(log(sr+1) ~ srtime + treatment +  obs + 
 wt + haemoporc + treatment*srtime
                 , random=~ idh(id):units + period
                 , prior= prior1.4, family="gaussian"
                 , pl=TRUE, pr= TRUE
                 , saveX = TRUE , saveZ = TRUE
                 , data=lightdark, DIC= T
                 , singular.ok = TRUE
                 , nitt=5000, thin=10, burnin=50,verbose = F)  
# nitt=5000000, thin=100, burnin=500000

          
#### MCMCglmm prior comparison for fixed effects
          

MME <- c("Intercept", "Time","Treatment", "Occcasion", 
 "Weight", "Haemocyanin concentration", "Treatment x Time") 
RI = c("Hermit Crab ID", "Period", "Occcasion")
          
Prior_1_MME =as_tibble(cbind("MME" = MME, 
"Effective size" = effectiveSize(M1.1$Sol)[1:7], 
"Autocorrelation diagnostics" = (diag(autocorr(M1.1$Sol)[2,,]))[1:7]))
          
Prior_1_COV= as_tibble(cbind( "RI" = RI, 
"Effective size" = effectiveSize(M1.1$VCV), 
"Autocorrelation diagnostics" = diag(autocorr(M1.1$VCV)[2,,])))
          
Prior_2_MME =as_tibble(cbind("MME" = MME, 
"Effective size" = effectiveSize(M1.2$Sol)[1:7],
"Autocorrelation diagnostics" = (diag(autocorr(M1.2$Sol)[2,,]))[1:7]))
          
Prior_2_COV= as_tibble(cbind( "RI" = RI, 
"Effective size" = effectiveSize(M1.2$VCV), 
"Autocorrelation diagnostics" = diag(autocorr(M1.2$VCV)[2,,])))
          
Prior_3_MME =as_tibble(cbind("MME" = MME, 
"Effective size" = effectiveSize(M1.3$Sol)[1:7], 
"Autocorrelation diagnostics" = (diag(autocorr(M1.3$Sol)[2,,]))[1:7]))
          
Prior_3_COV= as_tibble(cbind( "RI" = RI, 
"Effective size" = effectiveSize(M1.3$VCV),
"Autocorrelation diagnostics" = diag(autocorr(M1.3$VCV)[2,,])))
          
Prior_4_MME =as_tibble(cbind("MME" = MME, 
"Effective size" = effectiveSize(M1.4$Sol)[1:7],
"Autocorrelation diagnostics" = (diag(autocorr(M1.4$Sol)[2,,]))[1:7]))
          
Prior_4_COV= as_tibble(cbind( "RI" = RI, 
"Effective size" = effectiveSize(M1.4$VCV), 
"Autocorrelation diagnostics" = diag(autocorr(M1.4$VCV)[2,,])))
          
Priors_MME<- as_tibble(rbind(Prior_1_MME,Prior_2_MME,Prior_3_MME,Prior_4_MME)) 
Priors_COV<- as_tibble(rbind(Prior_1_COV,Prior_2_COV,Prior_3_COV,Prior_4_COV)) 
          
Priors_MME<- Priors_MME %>% 
  mutate_at(c( "Effective size", "Autocorrelation diagnostics"), as.numeric) %>% 
            mutate_at(vars( "Effective size"), funs(round(., 0))) %>%
            mutate_at(vars("Autocorrelation diagnostics"), funs(round(., 4)))
          
Priors_COV<- Priors_COV %>%
  mutate_at(c( "Effective size", "Autocorrelation diagnostics"), as.numeric) %>% 
  mutate_at(vars( "Effective size"), funs(round(., 0))) %>% 
  mutate_at(vars("Autocorrelation diagnostics"), funs(round(., 4)))

htmlTable(as.data.frame(Priors_MME[,-1]), 
          rnames = paste(c(Priors_MME$MME)),
          rgroup = c("Prior 1","Prior 2", "Prior 3", "Prior 4"),
          n.rgroup = c(7, 7, 7, 7), 
          caption="Table S1: Effective size and autocorrelation diagnostic of fixed effects obtained from the MCMC model with four different priors",
          tfoot="To detect proper convergence, autocorrelation diagnostics should be < 0.1")
          
            
#### MCMCglmm prior comparison for random effects

htmlTable(as.data.frame(Priors_COV[,-1]), 
          rnames = paste(c(Priors_COV$RI)),
          rgroup = c("Prior 1","Prior 2", "Prior 3", "Prior 4"),
          n.rgroup = c(3,3,3,3), 
          caption="Table S2: Effective size and autocorrelation diagnostic of random effects obtained from the MCMC model with four different priors",
          tfoot="To detect proper convergence, autocorrelation diagnostics should be < 0.1")

### Model selection: comparing Deviance Information Criterion (DIC)
            

Value<- c(M1.1$DIC,  M1.2$DIC,  M1.3$DIC,  M1.4$DIC)
          
Prior<- c("Prior 1", "Prior 2", "Prior 3", "Prior 4")
          
Results = as.matrix(data.frame(Prior,Value)) 

htmlTable(Results, 
          ctable=c("solid", "double"),
          cspan.rgroup = 2,
          caption="Table S3: Comparion of DIC output between four different priors")

         
#### The effect of permanent light on the startle response duration of hermit crabs
          
preds<- c(predict(M1.4, marginal = NULL, type = "response")) ##Get the predicted startle response (for plotting)
          
lightdark$prediction = preds
          
g1<- lightdark %>% 
  group_by(treatment) %>% 
  ggplot(aes(x =time_treatment, y = prediction, fill =time_treatment)) +
  geom_boxplot(size = 0.75, width = 0.25) + 
  labs(x = "Treatment", y = "Startle response duration (predicted)") +
  scale_fill_manual(values=c("orange", "orange", "#172c3d", "orange")) + 
  guides(fill=FALSE) +
  scale_x_discrete(labels = c("LD during day", "PL during day", 
            "LD during night", "PL during night")) + 
  theme_classic() 
          
g2<- lightdark %>% 
  group_by(treatment) %>% 
  ggplot(aes(x =treatment, y = prediction, fill =treatment)) +
  geom_boxplot(size = 0.75, width = 0.25) +
  labs(x = "Treatment", y = "Startle response duration (predicted)") +
  scale_fill_manual(values=c("#172c3d", "orange")) + 
  guides(fill=FALSE)  +
  theme_classic() 
          

g2 + g1 + 
  plot_annotation(tag_levels = c('A', 'B'), 
    title ='The effect of permanent light on statle response duration', 
    caption = 'The effect of permanent light and the interaction (b) between 
    treatment and time (day or night) on the startle response duration') 



#### Summary of the results:
            

MME1 <- c("Intercept", "Time","Treatment","Occcasion", 
          "Weight", "Haemocyanin concentration", "Treatment x Time") 
          
M1.results<-as.data.frame(summary(M1.4)$solutions)
M1.results<- tibble::rownames_to_column(M1.results, "Parameter name")
          
M1.results<- M1.results %>% 
  dplyr::select("Parameter name", "Posterior mean"="post.mean", 
                "95% CI lower"="l-95% CI", "95% CI upper"="u-95% CI",
                "p"="pMCMC") %>% 
  mutate_at(c( "Posterior mean","95% CI lower","95% CI upper","p"),as.numeric) %>%
  mutate_at(vars("Posterior mean","95% CI lower","95% CI upper"), funs(round(., 4))) %>%
  mutate_at(vars("p"), funs(round(., 4)))
          
htmlTable(as.data.frame(M1.results[,-1]), 
          rnames = paste(c(MME1)),
          caption="Table S4: Posterior summary statistics for the mean effect of startle response, showing posterior mean, lower and upper 95% CIs and P-values. Results are given for fixed effects only")

Gcov<-summary(M1.4)$Gcovariances
Random_effects<- as_tibble(rbind(Gcov, "Rcov" = summary(M1.4)$Rcovariances))
RE<- c("Hermit Crab ID", "Period", "Occcasion")
          
Random_effects = Random_effects %>% 
  dplyr::select("Posterior mean"="post.mean", "95% CI lower"="l-95% CI", 
                "95% CI upper"="u-95% CI") %>% 
  mutate_at(c( "Posterior mean","95% CI lower","95% CI upper"),as.numeric) %>%
  mutate_at(vars("Posterior mean","95% CI lower","95% CI upper"), funs(round(., 4)))
          
htmlTable(as.data.frame(Random_effects), 
  rnames = paste(c(RE)),
  caption="Table S5: Posterior summary statistics for the mean effect of startle response, showing posterior mean, lower and upper 95% CIs. Results are given for random effects only.")

# 2. Comparing the repeatability and variance components of startle responses
            
## Model selection: comparing priors

prior2.1<- list(R=list(V=diag(4)/3.002, nu= 3.002),
                G=list(G1=list(V=diag(4)/3.002, nu= 3.002)))

prior2.2<- list(R=list(V=diag(4)/1.002, nu= 1.002), 
                G=list(G1=list(V=diag(4)/1.002, nu= 1.002)))

prior2.3<- list(R=list(V=diag(4)*(1e-6/3.002), nu= 3.002), 
                G=list(G1=list(V=diag(4)*(1e-6/3.002), nu= 3.002)))

prior2.4<- list(R=list(V=diag(4)/3, nu=3), 
                G=list(G1=list(V=diag(4)/3, nu=3)))

prior2.5<- list(R=list(V=diag(4), nu=1.002), 
                G=list(G1=list(V=diag(4), nu=1.002)))


IIV_matrix<-matrix(c("VID_Day_LD",0,0,0,
   0,"VID_Night_LD",0,0,
   0,0,"VID_Day_PL",0,
   0,0,0,"VID_Night_PL"),
 nrow=4,ncol=4,byrow=T)

###Model 1
options(warn = -1) #deactivate warning for the MCMCglmm

M2.1 <-MCMCglmm(log(sr+1)~  srtime + treatment + obs + 
                  wt + haemoporc + treatment*srtime + treatment*period
                , random=~ idh(time_treatment):id
                , rcov=~idh(time_treatment):units
                , prior= prior2.1
                , family="gaussian"
                , data=lightdark
                , singular.ok=TRUE
                , DIC= T
                , nitt=5000, thin=10, burnin=50,verbose = F)
#, nitt=5000000, thin=100, burnin=500000,verbose = F

###Model 2

M2.2 <-  MCMCglmm(log(sr+1)~  srtime + treatment + obs + 
                    wt + haemoporc + treatment*srtime + treatment*period
                  , random=~ idh(time_treatment):id
                  , rcov=~idh(time_treatment):units
                  , prior= prior2.2
                  , family="gaussian"
                  , data=lightdark
                  , singular.ok=TRUE
                  , DIC= T
                  , nitt=5000, thin=10, burnin=50,verbose = F)
#, nitt=5000000, thin=100, burnin=500000,verbose = F

### Model 3

M2.3 <- MCMCglmm(log(sr+1)~  srtime + treatment + obs + 
                   wt + haemoporc + treatment*srtime + treatment*period
                 , random=~ idh(time_treatment):id
                 , rcov=~idh(time_treatment):units
                 , prior= prior2.3
                 , family="gaussian"
                 , data=lightdark
                 , singular.ok=TRUE
                 , DIC= T
                 , nitt=5000, thin=10, burnin=50,verbose = F)
#, nitt=5000000, thin=100, burnin=500000,verbose = F

### Model 4

M2.4 <-MCMCglmm(log(sr+1)~  srtime + treatment + obs + 
                  wt + haemoporc + treatment*srtime + treatment*period
                , random=~ idh(time_treatment):id
                , rcov=~idh(time_treatment):units
                , prior= prior2.4
                , family="gaussian"
                , data=lightdark
                , singular.ok=TRUE
                , DIC= T
                , nitt=5000, thin=10, burnin=50,verbose = F)
#, nitt=5000000, thin=100, burnin=500000,verbose = F

### Model 5

M2.5 <- MCMCglmm(log(sr+1)~ srtime + treatment + obs +  
                   wt + haemoporc + treatment*srtime + treatment*period
                 , random=~ idh(time_treatment):id
                 , rcov=~idh(time_treatment):units
                 , prior= prior2.5
                 , family="gaussian"
                 , data=lightdark
                 , singular.ok=TRUE
                 , DIC= T
                 , nitt=5000, thin=10, burnin=50,verbose = F)
#, nitt=5000000, thin=100, burnin=500000,verbose = F



MME <- c("Intercept", "Time","Treatment", "Occcasion", 
         "Weight", "Haemocyanin concentration","Period", 
         "Treatment x Time", "Treatment x Period") 

Prior_1 =as_tibble(cbind("MME" = MME, "MME solutions (ES)" = effectiveSize(M2.1$Sol),
       "(co)variance matrices(ES)" = effectiveSize(M2.1$VCV), 
       "MME solutions (AC)" = diag(autocorr(M2.1$Sol)[2, , ]), 
       "(co)variance matrices (AC)" = diag(autocorr(M2.1$VCV)[2, , ])))

Prior_2=as_tibble(cbind("MME" = MME, "MME solutions (ES)" = effectiveSize(M2.2$Sol),
      "(co)variance matrices(ES)"= effectiveSize(M2.2$VCV), 
      "MME solutions (AC)" = diag(autocorr(M2.2$Sol)[2, , ]), 
      "(co)variance matrices (AC)" = diag(autocorr(M2.2$VCV)[2, , ])))

Prior_3=as_tibble(cbind("MME" = MME, "MME solutions (ES)" = effectiveSize(M2.3$Sol), 
      "(co)variance matrices(ES)" = effectiveSize(M2.3$VCV), 
      "MME solutions (AC)" = diag(autocorr(M2.3$Sol)[2, , ]), 
      "(co)variance matrices (AC)" = diag(autocorr(M2.3$VCV)[2, , ])))

Prior_4=as_tibble(cbind("MME" = MME, "MME solutions (ES)" = effectiveSize(M2.4$Sol), 
      "(co)variance matrices(ES)" = effectiveSize(M2.4$VCV),
      "MME solutions (AC)" = diag(autocorr(M2.4$Sol)[2, , ]), 
      "(co)variance matrices (AC)" = diag(autocorr(M2.4$VCV)[2, , ])))

Prior_5=as_tibble(cbind("MME" = MME, "MME solutions (ES)" = effectiveSize(M2.5$Sol), 
      "(co)variance matrices(ES)" = effectiveSize(M2.5$VCV), 
      "MME solutions (AC)" = diag(autocorr(M2.5$Sol)[2, , ]),
      "(co)variance matrices (AC)" = diag(autocorr(M2.5$VCV)[2, , ])))


Priors<- as_tibble(rbind(Prior_1,Prior_2,Prior_3,Prior_4,Prior_5)) 

Priors<-Priors %>% 
  mutate_at(c( "MME solutions (ES)", "(co)variance matrices(ES)", 
               "MME solutions (AC)","(co)variance matrices (AC)"),as.numeric) %>% 
  mutate_at(vars( "MME solutions (ES)", "(co)variance matrices(ES)"), funs(round(., 0))) %>%
  mutate_at(vars("MME solutions (AC)","(co)variance matrices (AC)"), funs(round(., 4)))

htmlTable(as.data.frame(Priors[,-1]), 
          header =   paste(c("MME solutions", "(co)variance matrices", 
           "MME solutions","(co)variance matrices")),
          rnames = paste(c(Priors$MME)),
          rgroup = c("Prior 1","Prior 2", "Prior 3", "Prior 4", "Prior 5"),
          n.rgroup = c(9, 9,9, 9, 9), 
          cgroup = c("Effective size","Autocorrelation diagnostics"),
          n.cgroup = c(2,2), 
          caption="Table S6: Effective size and autocorrelation diagnostic the MCMC model with five different priors",
          tfoot="To detect proper convergence, autocorrelation diagnostics should be < 0.1")

### Model selection: comparing Deviance Information Criterion (DIC)
  

Value<- c(M2.1$DIC,  M2.2$DIC,  M2.3$DIC,  M2.4$DIC,  M2.5$DIC)

Prior<- c("Prior 1", "Prior 2", "Prior 3", "Prior 4", "Prior 5")

Results = as.matrix(data.frame(Prior,Value)) 
htmlTable(Results,  
          ctable=c("solid", "double"),
          cspan.rgroup = 2,
          caption="Table S7: Comparion of DIC output between five different priors")



## Comparing repeatability across treatment groups and time on which startle response duration was measured
 
#Repeatability of day measures in the light and dark treatment
VID_Day_LD = M2.2$VCV[,1]; VR_Day_LD = M2.2$VCV[,5]
R_Day_LD = VID_Day_LD/(VID_Day_LD + VR_Day_LD)

#Repeatability of day measures in the permanent light treatment
VID_Day_PL = M2.2$VCV[,2]; VR_Day_PL = M2.2$VCV[,6] 
R_Night_LD = VID_Day_PL/(VID_Day_PL +  VR_Day_PL)

#Repeatability of  night measures in the light and dark treatment
VID_Night_LD = M2.2$VCV[,3]; VR_Night_LD = M2.2$VCV[,7]
R_Day_PL = VID_Night_LD/(VID_Night_LD + VR_Night_LD)

#Repeatability of night measures in the permanent light treatment

VID_Night_PL = M2.2$VCV[,4]; VR_Night_PL = M2.2$VCV[,8]
R_Night_PL = VID_Night_PL/(VID_Night_PL + VR_Night_PL)

#Difference in repeatability between day and night in the permanent light treatment
deltaR_PL = (R_Day_PL - R_Night_PL)

#Difference in repeatability between day and night in the light and dark treatment
deltaR_LD = (R_Day_LD - R_Night_LD)

#Difference in repeatability between permanent light and the light and dark treatment during day
deltaR_Day = (R_Day_PL - R_Day_LD)

#Difference in repeatability between day and night in the light and dark treatment
deltaR_Night = (R_Night_PL - R_Night_LD)


deltas<-  paste0("\u0394", c(" PL", " LD", " Day", " Night"))

Repeat<-data_frame(
  name = c("Day LD", "Night LD", "Day PL", "Night PL",deltas),
  mean =  c(posterior.mode(R_Day_LD), posterior.mode(R_Night_LD),
            posterior.mode(R_Day_PL), posterior.mode(R_Night_PL), 
            posterior.mode(deltaR_PL), posterior.mode(deltaR_LD),
            posterior.mode(deltaR_Day), posterior.mode(deltaR_Night)), 
  lower =  c(HPDinterval(R_Day_LD)[1], HPDinterval(R_Night_LD)[1],
             HPDinterval(R_Day_PL)[1], HPDinterval(R_Night_PL)[1], 
             HPDinterval(deltaR_PL)[1], HPDinterval(deltaR_LD)[1], 
             HPDinterval(deltaR_Day)[1], HPDinterval(deltaR_Night)[1]),
  upper =  c(HPDinterval(R_Day_LD)[2], HPDinterval(R_Night_LD)[2],
             HPDinterval(R_Day_PL)[2], HPDinterval(R_Night_PL)[2],
             HPDinterval(deltaR_PL)[2], HPDinterval(deltaR_LD)[2],
             HPDinterval(deltaR_Day)[2], HPDinterval(deltaR_Night)[2]),
  Type = c("Repeatability", "Repeatability", "Repeatability", "Repeatability", "\u0394", "\u0394", "\u0394", "\u0394"))



plt1 <- ggplot(Repeat, aes(x = name, y = mean, ymin = upper, ymax = lower)) +
  geom_linerange(aes(color =  name), size = 1, alpha = 0.5) +
  geom_point(aes(color = name), 
             position=position_dodge(width=c(0.6,0.4)), size = 3) +
  coord_flip() + labs(y = "Posterior mean", x = "Treatment") + 
  geom_hline(yintercept = 0)  + 
  labs(colour= "Repeatability estimates") +
  theme_bw()


plt1 +  scale_color_manual(values=c("orange", "orange","#172c3d", "orange",
"orange", "#a88f0f", "#a88f0f", "orange" )) +
  scale_fill_manual(values=c("orange", "orange","#172c3d", "orange", 
           "orange", "#a88f0f", "#a88f0f", "orange")) +
  theme(legend.position = "none") +
  labs(title = "Repeatability estimate of the startle response duration across treatment groups and time")


Repeat <- Repeat %>% 
  mutate_at(c( "mean","lower", "upper"),as.numeric) %>% 
  mutate_at(vars("mean","lower", "upper"), funs(round(., 3))) %>%
  select("name", "mean","lower", "upper")

repeat_L<- c("Repeatability estimate", paste0("\u0394", " Repeatability"))
htmlTable(as.data.frame(Repeat[,-1]), 
          header =   paste(c("Posterior mean","95% CI lower","95% CI upper")),
          rnames = paste(c(Repeat$name)),
          rgroup = c(repeat_L),
          n.rgroup = c(4, 4), 
          #          col.rgroup = c(rows=c(AD1),"#c95400"),
          #          col.rgroup = c(rows=c(AD2),"#c95400"),
          caption="Table S8: Repeatability estimates and differences in repeatability (∆R) in startle response duration between treatment groups within time",
          tfoot= "Differences between treatments and time (∆R) was estimated as PL - LD and was Day - Night")

  
## Comparing between-individual variation across treatment groups and time on which startle response duration was measured
  

#Between-individual variation during PL treatment morning versus night 
deltaVID_PL = VID_Day_PL - VID_Night_PL    
#Between-individual variation during LD treatment morning versus night 
deltaVID_LD = VID_Day_LD -VID_Night_LD   
#Between-individual variation morning PL treatment versus LD treatment  
deltaVID_Day = VID_Day_PL - VID_Day_LD       
#Between-individual variation night PL treatment versus LD treatment  
deltaVID_Night = VID_Night_PL - VID_Night_LD       

VID<-data_frame(
  name = c("Day LD", "Night LD", "Day PL", "Night PL",deltas),
  mean =  c(posterior.mode(VID_Day_LD), posterior.mode(VID_Night_LD),
            posterior.mode(VID_Day_PL), posterior.mode(VID_Night_PL),
            posterior.mode(deltaVID_PL), posterior.mode(deltaVID_LD),
            posterior.mode(deltaVID_Day), posterior.mode(deltaVID_Night)), 
  lower =  c(HPDinterval(VID_Day_LD)[1], HPDinterval(VID_Night_LD)[1],
             HPDinterval(VID_Day_PL)[1], HPDinterval(VID_Night_PL)[1],
             HPDinterval(deltaVID_PL)[1], HPDinterval(deltaVID_LD)[1],
             HPDinterval(deltaVID_Day)[1], HPDinterval(deltaVID_Night)[1]),
  upper =  c(HPDinterval(VID_Day_LD)[2], HPDinterval(VID_Night_LD)[2],
             HPDinterval(VID_Day_PL)[2], HPDinterval(VID_Night_PL)[2],
             HPDinterval(deltaVID_PL)[2], HPDinterval(deltaVID_LD)[2],
             HPDinterval(deltaVID_Day)[2], HPDinterval(deltaVID_Night)[2]),
  Type = c("Repeatability", "Repeatability", "Repeatability", 
           "Repeatability", "\u0394", "\u0394", "\u0394", "\u0394"))

plt2 <- ggplot(VID, aes(x = name, y = mean, ymin = upper, ymax = lower)) +
  geom_linerange(aes(color =  name), size = 1, alpha = 0.5) +
  geom_point(aes(color = name), 
             position=position_dodge(width=c(0.6,0.4)), size = 3) +
  coord_flip() + 
  labs(y = "Posterior mean", x = "Treatment") +  
  
  labs(colour= "Between-individual variation") +
  geom_hline(yintercept = 0)  + 
  theme_bw()


plt2 +  scale_color_manual(values=c("orange", "orange","#172c3d", "orange",
"orange", "#a88f0f", "#a88f0f", "orange" )) +
  scale_fill_manual(values=c("orange", "orange","#172c3d", "orange",
           "orange", "#a88f0f", "#a88f0f", "orange")) +
  theme(legend.position = "none") +
  labs(title = "Between-individual variation of the startle response duration across treatment groups and time")


VID <- VID %>% 
  mutate_at(c( "mean","lower", "upper"),as.numeric) %>% 
  mutate_at(vars("mean","lower", "upper"), funs(round(., 3))) %>% 
  select("name", "mean","lower", "upper")

repeat_VID<- c("Repeatability estimate", paste0("\u0394", " Variance"))

htmlTable(as.data.frame(VID[,-1]), 
          header =   paste(c("Posterior mean","95% CI lower","95% CI upper")),
          rnames = paste(c(Repeat$name)),
          rgroup = c(repeat_VID),
          n.rgroup = c(4, 4), 
          caption="Table S9:Between-individual variation (VBI) and differences in VBI (∆V) in startle response duration between treatment groups within time",
          tfoot= "Differences between treatments and time (∆V) was estimated as PL - LD and was Day - Night")

  
## Comparing within-individual variation across treatment groups and time on which startle response duration was measured
  

# Within-individual variation during PL treatment morning versus night  
deltaVR_PL = VR_Day_PL - VR_Night_PL 
# Within-individual variation during LD treatment morning versus night   
deltaVR_LD = VR_Day_LD - VR_Night_LD     
# Within-individual variation morning PL treatment versus LD treatment   
deltaVR_Day = VID_Day_PL - VID_Day_LD        
# Within-individual variation night PL treatment versus LD treatment   
deltaVR_Night = VR_Night_PL - VR_Night_LD       

WID<-data_frame(
  name = c("Day LD", "Night LD", "Day PL", "Night PL",deltas),
  mean =  c(posterior.mode(VID_Day_LD), posterior.mode(VID_Night_LD),
            posterior.mode(VID_Day_PL), posterior.mode(VID_Night_PL),
            posterior.mode(deltaVR_PL), posterior.mode(deltaVR_LD), 
            posterior.mode(deltaVR_Day), posterior.mode(deltaVR_Night)), 
  lower =  c(HPDinterval(VID_Day_LD)[1], HPDinterval(VID_Night_LD)[1],
             HPDinterval(VID_Day_PL)[1], HPDinterval(VID_Night_PL)[1],
             HPDinterval(deltaVR_PL)[1], HPDinterval(deltaVR_LD)[1], 
             HPDinterval(deltaVR_Day)[1], HPDinterval(deltaVR_Night)[1]),
  upper =  c(HPDinterval(VID_Day_LD)[2], HPDinterval(VID_Night_LD)[2],
             HPDinterval(VID_Day_PL)[2], HPDinterval(VID_Night_PL)[2],
             HPDinterval(deltaVR_PL)[2], HPDinterval(deltaVR_LD)[2],
             HPDinterval(deltaVR_Day)[2], HPDinterval(deltaVR_Night)[2]),
  Type = c("Repeatability", "Repeatability", "Repeatability", 
           "Repeatability", "\u0394", "\u0394", "\u0394", "\u0394"), 
  Color = c("#a37d0b", "#52514e","#f2b90f", "#d19f0a", "#eda909", 
            "#7d786e", "#eda909", "#7d786e" ))



plt3 <- ggplot(WID, aes(x = name, y = mean, ymin = upper, ymax = lower)) +
  geom_linerange(aes(color =  name), size = 1, alpha = 0.5) +
  geom_point(aes(color = name), position=position_dodge(width=c(0.6,0.4)), size = 3) +
  coord_flip() + labs(y = "Posterior mean", x = "Treatment") +  
  labs(colour= "Within-individual variation") +
  geom_hline(yintercept = 0)  + 
  theme_bw()

plt3 +  scale_color_manual(values=c("orange", "orange","#172c3d", "orange", 
"orange", "#a88f0f", "#a88f0f", "orange" )) +
  scale_fill_manual(values=c("orange", "orange","#172c3d", "orange",
           "orange", "#a88f0f", "#a88f0f", "orange")) +
  theme(legend.position = "none") +
  labs(title = "Within-individual variation of the startle response duration across treatment groups and time")


WID <- WID %>% 
  mutate_at(c( "mean","lower", "upper"),as.numeric) %>% 
  mutate_at(vars("mean","lower", "upper"), funs(round(., 3))) %>% 
  select("name", "mean","lower", "upper")

repeat_VID<- c("Repeatability estimate", paste0("\u0394", " Variance"))

htmlTable(as.data.frame(WID[,-1]), 
          header =   paste(c("Posterior mean","95% CI lower","95% CI upper")),
          rnames = paste(c(Repeat$name)),
          rgroup = c(repeat_VID),
          n.rgroup = c(4, 4), 
          caption="Table S10: Within-individual variation (VWI) and differences in VWI (∆V) in startle response duration between treatment groups within time",
          tfoot= "Differences between treatments and time (∆V) was estimated as PL - LD and was Day - Night")


  
# 3. The effect of permanent light on metabolic rate
  
MR=read.csv("data/MR1.csv", header=TRUE)
id<-as.factor(MR$ID)
period<-as.factor(MR$Period)
Treat<-as.factor(MR$Treatment)
time<-as.factor(MR$Time)
metR<-as.numeric(MR$metabolic.rate)
Period_AOV <- lmer(log10(metR+1) ~ period + (1|id), data=MR)
Pr_result<- as.data.frame(anova(Period_AOV))

Period_L <- c("Sum Sq  ","Mean Sq  ", "NumDF  " , "DenDF  "  , "F value  ", "p  ") 

Pr_result = Pr_result %>% mutate_at(c( "Sum Sq","Mean Sq", "NumDF" ,
    "DenDF", "F value", "Pr(>F)" ),as.numeric) %>%
  mutate_at(vars( "Sum Sq","Mean Sq", "NumDF" ,"DenDF", "F value", "Pr(>F)"), funs(round(., 4)))

htmlTable(as.data.frame(Pr_result), 
          header =   paste(c(Period_L)),
          caption="Table S11: Results of repeated measures ANOVA, testing whether metabolic rate varied according with the period on which it was collected (A or B)",
          rnames = paste("Perid"))

MR_AOV <- lmer(log10(metR+1) ~ time*Treat + (1|id), data=MR)

MR_result<- as.data.frame(anova(MR_AOV))

MR__L <- c("Sum Sq  ","Mean Sq  ", "NumDF  " , "DenDF  "  , "F value  ", "p  ") 

MR_result = MR_result %>% mutate_at(c( "Sum Sq","Mean Sq", "NumDF" ,"DenDF", 
                                       "F value", "Pr(>F)" ),as.numeric) %>% 
  mutate_at(vars( "Sum Sq","Mean Sq", "NumDF" ,"DenDF", "F value", "Pr(>F)"), funs(round(., 4)))

htmlTable(as.data.frame(MR_result), 
          header =   paste(c(Period_L)),
          caption="Table S12: Results of repeated measures ANOVA, testing the effect of permanent light on metabolic rate",
          rnames = paste(c("Time", "Treatment", "Treatment x Time")))

MR %>% 
  group_by(Treatment) %>% 
  ggplot(aes(y =log10(metabolic.rate +1), x = Treatment, fill =Treatment)) +
  geom_boxplot(size = 0.75, width = 0.25) + labs(x = "Treatment", y = "Metabolic rate") +
  theme_classic() + coord_flip() +
  scale_fill_manual(values=c("#172c3d", "orange")) + guides(fill=FALSE) +
  scale_color_manual(values=c("#a37d0b", "#d1a615")) +
  scale_x_discrete(labels = c("LD treatment", "PL Treatment")) +
  labs(title = "The effect of permanent light on the metabolic rate (MO2) ")


