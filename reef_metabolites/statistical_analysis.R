# statistical analysis
# bucket table (glmer) & curated feature data (glmer and indicator analysis including heatmap Fig. 5 & 6 manuscript)

#___________________________________________________________
# install packages
library(tidyverse)
library(lme4)
library(multcomp)
library(dplyr)
library(vegan)
library(performance)
library(rstatix)
library(car)
library(indicspecies)
library(ggpubr)
library(grid)
#____________________________________________________________
# M2NT
all <- read_csv("data_prep/data_metabolites_nodes_prep.csv")
all <- all %>% select(-c(...1))

#turn long for analysis
all_long <- all %>% 
  pivot_longer(!Fragment_ID : Treatment, names_to = "Features", values_to = "intensity")
#turn intensity binary
all_long$intensity <- ifelse(all_long$intensity >= 1, 1, 0)

#________________overall__________________
# treatment
#model_treat_lmer <- lmer((log_intensity+1) ~ Treatment + (1|Day), data = all_long)
#shapiro.test(residuals(model_treat_lmer))
#hist(residuals(model_treat_lmer), breaks = 20, col = "lightblue", main = "Histogram of Residuals")
#residuals not normally distributed -> glmer

model_treat <- glmer((intensity) ~ Treatment + (1|Day), 
                     family = "binomial", data = all_long)
cftest(model_treat)
summary(glht(model_treat, linfct = mcp(Treatment = "Tukey")),
        test = adjusted("holm"))
rm(model_treat)

# species
#model_spec_lmer <- lmer((log_intensity+1) ~ Species + (1|Day), data = all_long)
#shapiro.test(residuals(model_spec_lmer))
#hist(residuals(model_spec_lmer), breaks = 20, col = "lightblue", main = "Histogram of Residuals")
# residuals not normally distributed -> glmer

model_spec <- glmer((intensity) ~ Species + (1|Day), 
                    family = "binomial", data = all_long)
cftest(model_spec)
summary(glht(model_spec, linfct = mcp(Species = "Tukey")),
        test = adjusted("holm"))
rm(model_spec)

# organism group
#model_org_lmer <- lmer((log_intensity+1) ~ Organism_group + (1|Day), data = all_long)
#hist(residuals(model_org_lmer), breaks = 20, col = "lightblue", main = "Histogram of Residuals")
# residuals not normally distributed

model_org <- glmer((intensity) ~ Organism_group + (1|Day), 
                   family = "binomial", data = all_long)
cftest(model_org)
summary(glht(model_org, linfct = mcp(Organism_group = "Tukey")),
        test = adjusted("holm"))
rm(model_org)

#____________within treatment_____________
mono <- all_long %>% 
  subset(Treatment != "controlled_Polyculture") %>%
  subset(Treatment != "seminatural_Polyculture")

Pc <- all_long %>% 
  subset(Treatment !="Monoculture") %>%
  subset(Treatment != "seminatural_Polyculture")

Ps <- all_long %>% 
  subset(Treatment !="Monoculture")%>%
  subset(Treatment != "controlled_Polyculture")

# Monoculture
#model_mono_lmer <- lmer((log_intensity+1) ~ Species + (1|Day), data = mono)
#shapiro.test(residuals(model_mono_lmer))
#hist(residuals(model_mono_lmer), breaks = 20, col = "lightblue", main = "Histogram of Residuals")
# residuals not normally distributed

model_mono <- glmer((intensity) ~ Species + (1|Day), 
                    family = "binomial", data = mono)
cftest(model_mono)
summary(glht(model_mono, linfct = mcp(Species = "Tukey")),
        test = adjusted("holm"))
rm(model_mono)

# Polyculture controlled
#model_poly_lmer <- lmer((log_intensity+1) ~ Species + (1|Day), data = Pc)
#hist(residuals(model_poly_lmer), breaks = 20, col = "lightblue", main = "Histogram of Residuals")
# residuals not normally distributed

model_poly <- glmer((intensity) ~ Species + (1|Day), 
                    family = "binomial", data = Pc)
cftest(model_poly)
summary(glht(model_poly, linfct = mcp(Species = "Tukey")),
        test = adjusted("holm"))
rm(model_poly)

# Polyculture seminatural
#model_semi_lmer <- lmer((log_intensity+1) ~ Species + (1|Day), data = Ps)
#hist(residuals(model_semi_lmer), breaks = 20, col = "lightblue", main = "Histogram of Residuals")
# residuals not normally distributed

model_semi <- glmer((intensity) ~ Species + (1|Day), 
                    family = "binomial", data = Ps)
cftest(model_semi)
summary(glht(model_semi, linfct = mcp(Species = "Tukey")),
        test = adjusted("holm"))
rm(model_semi)
rm(mono, Ps, Pc)

#________within species___________________
Cau <- all_long %>% 
  subset(Species != "Peyssonnelia sp.") %>%
  subset(Species != "Sinularia sp.") %>%
  subset(Species != "Xenia sp.") %>%
  subset(Species != "Montipora digitata") %>%
  subset(Species != "Pocillopora verrucosa") %>%
  subset(Species != "Haliclona cnidata")

Pey <- all_long %>% 
  subset(Species != "Caulerpa sp.") %>%
  subset(Species != "Sinularia sp.") %>%
  subset(Species != "Xenia sp.") %>%
  subset(Species != "Montipora digitata") %>%
  subset(Species != "Pocillopora verrucosa") %>%
  subset(Species != "Haliclona cnidata")

Sin <- all_long %>% 
  subset(Species != "Peyssonnelia sp.") %>%
  subset(Species != "Caulerpa sp.") %>%
  subset(Species != "Xenia sp.") %>%
  subset(Species != "Montipora digitata") %>%
  subset(Species != "Pocillopora verrucosa") %>%
  subset(Species != "Haliclona cnidata")

Xen <- all_long %>% 
  subset(Species != "Peyssonnelia sp.") %>%
  subset(Species != "Sinularia sp.") %>%
  subset(Species != "Caulerpa sp.") %>%
  subset(Species != "Montipora digitata") %>%
  subset(Species != "Pocillopora verrucosa") %>%
  subset(Species != "Haliclona cnidata")

Mdi <- all_long %>% 
  subset(Species != "Peyssonnelia sp.") %>%
  subset(Species != "Sinularia sp.") %>%
  subset(Species != "Xenia sp.") %>%
  subset(Species != "Caulerpa sp.") %>%
  subset(Species != "Pocillopora verrucosa") %>%
  subset(Species != "Haliclona cnidata")

Pve <- all_long %>% 
  subset(Species != "Peyssonnelia sp.") %>%
  subset(Species != "Sinularia sp.") %>%
  subset(Species != "Xenia sp.") %>%
  subset(Species != "Montipora digitata") %>%
  subset(Species != "Caulerpa sp.") %>%
  subset(Species != "Haliclona cnidata")

Hcn <- all_long %>% 
  subset(Species != "Peyssonnelia sp.") %>%
  subset(Species != "Sinularia sp.") %>%
  subset(Species != "Xenia sp.") %>%
  subset(Species != "Montipora digitata") %>%
  subset(Species != "Pocillopora verrucosa") %>%
  subset(Species != "Caulerpa sp.")

# Caulerpa
model_cau <- glmer((intensity) ~ Treatment + (1|Day), 
                   family = "binomial", data = Cau)
cftest(model_cau)
summary(glht(model_cau, linfct = mcp(Treatment = "Tukey")),
        test = adjusted("holm"))
rm(model_cau, Cau)


#______________________________________
# Peyssonnelia sp.

model_pey <- glmer((intensity) ~ Treatment + (1|Day), 
                   family = "binomial", data = Pey)
cftest(model_pey)
summary(glht(model_pey, linfct = mcp(Treatment = "Tukey")),
        test = adjusted("holm"))
rm(Pey, model_pey)

#________________________________________________
# Sinularia
model_sin <- glmer((intensity) ~ Treatment + (1|Day), 
                   family = "binomial", data = Sin)
cftest(model_sin)
summary(glht(model_sin, linfct = mcp(Treatment = "Tukey")),
        test = adjusted("holm"))
rm(model_sin, Sin)

#________________________________________________
# Xenia

model_xen <- glmer((intensity) ~ Treatment + (1|Day), 
                   family = "binomial", data = Xen)
cftest(model_xen)
summary(glht(model_xen, linfct = mcp(Treatment = "Tukey")),
        test = adjusted("holm"))
rm(model_xen, Xen)

#________________________________________________
# Mdi

model_mdi <- glmer((intensity) ~ Treatment + (1|Day), 
                   family = "binomial", data = Mdi)
cftest(model_mdi)
summary(glht(model_mdi, linfct = mcp(Treatment = "Tukey")),
        test = adjusted("holm"))
rm(model_mdi, Mdi)

#________________________________________________
# Pve

model_pve <- glmer((intensity) ~ Treatment + (1|Day), 
                   family = "binomial", data = Pve)
cftest(model_pve)
summary(glht(model_pve, linfct = mcp(Treatment = "Tukey")),
        test = adjusted("holm"))
rm(model_pve, Pve)

#________________________________________________
# Hcn

model_hcn <- glmer((intensity) ~ Treatment + (1|Day), 
                   family = "binomial", data = Hcn)
cftest(model_hcn)
summary(glht(model_hcn, linfct = mcp(Treatment = "Tukey")),
        test = adjusted("holm"))
rm(model_hcn, Hcn)

rm(all_long, all)

# ___________________________________________________________
# CFT
all <- read_csv("data_prep/all_decontaminated.csv")
all_long <- all %>% 
  pivot_longer(!Fragment_ID : Treatment, names_to = "Features", values_to = "intensity")
all_long$log_intensity <- log(all_long$intensity+1)


#________________overall__________________
# treatment
#model_treat_lmer <- lmer((log_intensity+1) ~ Treatment + (1|Day), data = all_long)
#shapiro.test(residuals(model_treat_lmer))
# residuals not normally distributed

model_treat <- glmer((log_intensity + 1) ~ Treatment + (1|Day), 
                     family = "Gamma", data = all_long)
cftest(model_treat)
summary(glht(model_treat, linfct = mcp(Treatment = "Tukey")),
        test = adjusted("holm"))
rm(model_treat)

# species
#model_spec_lmer <- lmer((log_intensity+1) ~ Species + (1|Day), data = all_long)
#shapiro.test(residuals(model_spec_lmer))
#hist(residuals(model_spec_lmer), breaks = 20, col = "lightblue", main = "Histogram of Residuals")
# residuals not normally distributed

model_spec <- glmer((log_intensity + 1) ~ Species + (1|Day), 
                    family = "Gamma", data = all_long)
cftest(model_spec)
summary(glht(model_spec, linfct = mcp(Species = "Tukey")),
        test = adjusted("holm"))
rm(model_spec)

# organism group
#model_org_lmer <- lmer((log_intensity+1) ~ Organism_group + (1|Day), data = all_long)
#hist(residuals(model_org_lmer), breaks = 20, col = "lightblue", main = "Histogram of Residuals")
# residuals not normally distributed

model_org <- glmer((log_intensity + 1) ~ Organism_group + (1|Day), 
                   family = "Gamma", data = all_long)
cftest(model_org)
summary(glht(model_org, linfct = mcp(Organism_group = "Tukey")),
        test = adjusted("holm"))
rm(model_org)

#____________within treatment_____________
mono <- subset(all_long, Treatment == "Monoculture")
poly <- subset(all_long, Treatment == "controlled_Polyculture")
semi <- subset(all_long, Treatment == "seminatural_Polyculture")

# Monoculture
#model_mono_lmer <- lmer((log_intensity+1) ~ Species + (1|Day), data = mono)
#shapiro.test(residuals(model_mono_lmer))
#hist(residuals(model_mono_lmer), breaks = 20, col = "lightblue", main = "Histogram of Residuals")
# residuals not normally distributed

model_mono <- glmer((log_intensity + 1) ~ Species + (1|Day), 
                    family = "Gamma", data = mono)
cftest(model_mono)
summary(glht(model_mono, linfct = mcp(Species = "Tukey")),
        test = adjusted("holm"))
rm(model_mono)


# Polyculture controlled
#model_poly_lmer <- lmer((log_intensity+1) ~ Species + (1|Day), data = poly)
#hist(residuals(model_poly_lmer), breaks = 20, col = "lightblue", main = "Histogram of Residuals")
# residuals not normally distributed

model_poly <- glmer((log_intensity + 1) ~ Species + (1|Day), 
                    family = "Gamma", data = poly)
cftest(model_poly)
summary(glht(model_poly, linfct = mcp(Species = "Tukey")),
        test = adjusted("holm"))
rm(model_poly)

# Polyculture seminatural
#model_semi_lmer <- lmer((log_intensity+1) ~ Species + (1|Day), data = semi)
#hist(residuals(model_semi_lmer), breaks = 20, col = "lightblue", main = "Histogram of Residuals")
# residuals not normally distributed

model_semi <- glmer((log_intensity + 1) ~ Species + (1|Day), 
                    family = "Gamma", data = semi)
cftest(model_semi)
summary(glht(model_semi, linfct = mcp(Species = "Tukey")),
        test = adjusted("holm"))
rm(model_semi)
rm(mono, poly, semi)

#________within species___________________
Cau <- subset(all_long, Species == "Caulerpa sp.")
Pey <- subset(all_long, Species == "Peyssonnelia sp.")
Xen <- subset(all_long, Species == "Xenia sp.")
Sin <- subset(all_long, Species == "Sinularia sp.")
Mdi <- subset(all_long, Species == "Montipora digitata")
Pve <- subset(all_long, Species == "Pocillopora verrucosa")
Hcn <- subset(all_long, Species == "Haliclona cnidata")

# Caulerpa
Cau$log_intensity <- log(Cau$intensity+1)
#model_cau_lmer <- lmer(intensity ~ Treatment + (1|Features), data = Cau)
#shapiro.test(residuals(model_cau_lmer))
# residuals not normally distributed

model_cau <- glmer((log_intensity+1) ~ Treatment + (1|Day), 
                   family = "Gamma", data = Cau)
cftest(model_cau)
summary(glht(model_cau, linfct = mcp(Treatment = "Tukey")),
        test = adjusted("holm"))
rm(model_cau, Cau)

#______________________________________
# Peyssonnelia sp.
Pey$log_intensity <- log(Pey$intensity+1) 
#model_pey_lmer <- lmer(intensity ~ Treatment + (1|Features), data = Pey)
#model_pey_log <- lmer((log_intensity + 100) ~ Treatment + (1|Features), data = Pey)
#model_pey_sqrt <- lmer((sqrt_intensity + 100) ~ Treatment + (1|Features), data = Pey)
#shapiro.test(residuals(model_pey_lmer))
#shapiro.test(residuals(model_pey_log))
#shapiro.test(residuals(model_pey_sqrt))
# data tranformation does not help with the residuals in the shapiro test

model_pey <- glmer((log_intensity + 1) ~ Treatment + (1|Day), 
                   family = "Gamma", data = Pey)
cftest(model_pey)
summary(glht(model_pey, linfct = mcp(Treatment = "Tukey")),
        test = adjusted("holm"))
rm(Pey, model_pey)

#________________________________________________
# Sinularia
Sin$log_intensity <- log(Sin$intensity+1)
#model_sin_lmer <- lmer(intensity ~ Treatment + (1|Features), data = Sin)
#shapiro.test(residuals(model_sin_lmer))

model_sin <- glmer((log_intensity + 1) ~ Treatment + (1|Day), 
                   family = "Gamma", data = Sin)
cftest(model_sin)
summary(glht(model_sin, linfct = mcp(Treatment = "Tukey")),
        test = adjusted("holm"))
rm(model_sin, Sin)

#________________________________________________
# Xenia
Xen$log_intensity <- log(Xen$intensity+1)
#model_xen_lmer <- lmer(intensity ~ Treatment + (1|Features), data = Xen)
#shapiro.test(residuals(model_xen_lmer))
model_xen <- glmer((log_intensity + 1) ~ Treatment + (1|Day) + (1|Features), 
                   family = "Gamma", data = Xen)
cftest(model_xen)
summary(glht(model_xen, linfct = mcp(Treatment = "Tukey")),
        test = adjusted("holm"))
rm(model_xen, Xen)

#________________________________________________
# Mdi
Mdi$log_intensity <- log(Mdi$intensity+1)
#model_mdi_lmer <- lmer(intensity ~ Treatment + (1|Features), data = Mdi)
#shapiro.test(residuals(model_mdi_lmer))
model_mdi <- glmer((log_intensity + 1) ~ Treatment + (1|Day), 
                   family = "Gamma", data = Mdi)
cftest(model_mdi)
summary(glht(model_mdi, linfct = mcp(Treatment = "Tukey")),
        test = adjusted("holm"))
rm(model_mdi, Mdi)

#________________________________________________
# Pve
Pve$log_intensity <- log(Pve$intensity+1)
#model_pve_lmer <- lmer(intensity ~ Treatment + (1|Features), data = Pve)
#shapiro.test(residuals(model_pve_lmer))
model_pve <- glmer((log_intensity + 1) ~ Treatment + (1|Day), 
                   family = "Gamma", data = Pve)
cftest(model_pve)
summary(glht(model_pve, linfct = mcp(Treatment = "Tukey")),
        test = adjusted("holm"))
rm(Pve, model_pve)

#________________________________________________
# Hcn
Hcn$log_intensity <- log(Hcn$intensity+1)
#model_hcn_lmer <- lmer(intensity ~ Treatment + (1|Features), data = Hcn)
#shapiro.test(residuals(model_hcn_lmer))
model_hcn <- glmer((log_intensity + 1) ~ Treatment + (1|Day), 
                   family = "Gamma", data = Hcn)
cftest(model_hcn)
summary(glht(model_hcn, linfct = mcp(Treatment = "Tukey")),
        test = adjusted("holm"))

rm(model_hcn, Hcn)
rm(all, all_long)


#__________________________________________________
# indicator analysis
all <- read_csv("data_prep/all_decontaminated.csv")

all <- all %>% select(-c(...1))
all <- all %>% 
  mutate(Treatment = gsub("_d1", "", Treatment)) %>%
  mutate(Treatment = gsub("_d2", "", Treatment)) %>%
  mutate(Treatment = gsub("_d3", "", Treatment)) %>%
  mutate(Treatment = gsub("_d", "", Treatment)) 
all_nc <- all %>% subset(Treatment != "control")

# this is necessary to compare the peaks
names(all_nc) = gsub(pattern = " : ", replacement = "", x = names(all_nc))
names(all_nc) = gsub(pattern = " ", replacement = "", x = names(all_nc))
names(all_nc) = gsub(pattern = "/", replacement = "", x = names(all_nc))
rm(all)
# ______________________________________________
# indicator metabolites/peaks for each treatment
set.seed(1102007)
all_nc_multi <-data.frame(all_nc[,c(6:343)]) 

cluster<-factor(all_nc$Treatment)
Multi<-multipatt(all_nc_multi, cluster)
#capture.output(Multi, file = "Peak_treat.txt")
summary(Multi, alpha = 0.01) 
Multi2 <- data.frame(Multi$sign)
Multi2 <- subset(Multi2,p.value<=0.01)
Multi2$Peak <- rownames(Multi2)
write_csv2(Multi2, "indicator/Peak_ind_treatment_deco.csv")

# create heatmap
sig_peaks  <- read_csv2("indicator/Peak_ind_treatment_deco.csv")
sig_peaks <- sig_peaks %>% mutate(Peak = gsub("X", "", Peak))
all_nc_long <- all_nc %>% pivot_longer(!Fragment_ID : Treatment, names_to = "Peak", values_to = "intensity")
data_filtered <- subset(all_nc_long, Peak %in% sig_peaks$Peak)
names(data_filtered)[7]<-paste("peak_area")

data_filtered <- data_filtered %>% 
  mutate(Fragment_ID = gsub("HcnMo4", "HcnMo_4", Fragment_ID))

data_filtered$Fragment_ID<- factor(data_filtered$Fragment_ID, levels = c("MdiMo_1", "MdiMo_2", "MdiMo_4", 
                                                                         "PveMo_1", "PveMo_2", "PveMo_3", "PveMo_4",
                                                                         "CauMo_2", "CauMo_3", "CauMo_4", 
                                                                         "PeyMo_1", "PeyMo_2", "PeyMo_3", "PeyMo_4",
                                                                         "SinMo_1", "SinMo_2", "SinMo_3", "SinMo_4", 
                                                                         "XenMo_1", "XenMo_2", "XenMo_3", "XenMo_4",
                                                                         "HcnMo_1", "HcnMo_2", "HcnMo_3", "HcnMo_4",
                                                                         "MdiPc_1", "MdiPc_2", "MdiPc_3", "MdiPc_4", 
                                                                         "PvePc_2", "PvePc_3", "PvePc_4", "PvePc_5",
                                                                         "CauPc_1", "CauPc_2", "CauPc_3", "CauPc_4", 
                                                                         "PeyPc_1", "PeyPc_2", "PeyPc_3", "PeyPc_4",
                                                                         "SinPc_1", "SinPc_2", "SinPc_3", "SinPc_4", 
                                                                         "XenPc_1", "XenPc_2", "XenPc_3", "XenPc_4",
                                                                         "HcnPc_1", "HcnPc_2", "HcnPc_3", "HcnPc_4",
                                                                         "MdiPs_1", "MdiPs_2", "MdiPs_3", "MdiPs_4", 
                                                                         "PvePs_1", "PvePs_2", "PvePs_3", "PvePs_4",
                                                                         "CauPs_1", "CauPs_3", 
                                                                         "PeyPs_1", "PeyPs_2", "PeyPs_3", "PeyPs_4",
                                                                         "SinPs_1", "SinPs_2", "SinPs_3", "SinPs_4", 
                                                                         "XenPs_1", "XenPs_2", "XenPs_4",
                                                                         "HcnPs_1", "HcnPs_2", "HcnPs_3", "HcnPs_4"))


data_filtered$Peak<- factor(data_filtered$Peak, levels = c("961.1s315.252mz", "923.1s487.360mz", "923.6s438.378mz",
                                                           "923.5s443.334mz", "923.1s421.352mz", "922.3s531.386mz", "888.3s438.378mz",
                                                           "872.0s438.378mz", "148.5s137.107mz", "689.5s332.222mz", "1005.1s109.101mz",
                                                           "45.1s132.102mz", "233.2s155.107mz", "773.8s890.384mz", "909.3s211.205mz",
                                                           "912.1s747.464mz", "912.5s675.445mz", "936.8s283.263mz"))

text_mo <- textGrob("Monoculture", gp=gpar(fontsize=10, fontface="plain"))
text_pc <- textGrob("Controlled Polyculture", gp=gpar(fontsize=10, fontface="plain"))
text_ps <- textGrob("Seminatural Polyculture", gp=gpar(fontsize=10, fontface="plain"))

text_mdi <-textGrob("M.digitata", gp=gpar(fontsize=8.5, fontface="italic"), rot = 45, just = 1)
text_pve <-textGrob("P. verrucosa", gp=gpar(fontsize=8.5, fontface="italic"), rot = 45, just = 1)
text_cau <-textGrob("Caulerpa sp.", gp=gpar(fontsize=8.5, fontface="italic"), rot = 45, just = 1)
text_pey <-textGrob("Peyssonnelia sp.", gp=gpar(fontsize=8.5, fontface="italic"), rot = 45, just = 1)  
text_sin <-textGrob("Sinularia sp.", gp=gpar(fontsize=8.5, fontface="italic"), rot = 45, just = 1) 
text_xen <-textGrob("Xenia sp.", gp=gpar(fontsize=8.5, fontface="italic"), rot = 45, just = 1) 
text_hcn <-textGrob("H. cnidata", gp=gpar(fontsize=8.5, fontface="italic"), rot = 45, just = 1)

line2 <-textGrob("___", gp=gpar(fontsize=8, fontface="plain"))
line3 <-textGrob("____", gp=gpar(fontsize=8, fontface="plain"))
line4 <-textGrob("_____", gp=gpar(fontsize=8, fontface="plain"))

group1 <-textGrob("Controlled 
Polyculture", gp=gpar(fontsize=8.5, fontface="plain", color = "lightgrey"), just = 0)
group2 <-textGrob("Monoculture & 
Controlled 
Polyculture", gp=gpar(fontsize=8.5, fontface="plain", color = "lightgrey"), just = 0)
group3 <-textGrob("Monoculture & 
Seminatural 
Polyculture", gp=gpar(fontsize=8.5, fontface="plain", color = "lightgrey"), just = 0)


# heatmap
treat_all <- ggplot(data_filtered, aes(x = Fragment_ID, y = Peak, fill = peak_area)) +
  geom_tile(colour="lightgrey", linewidth=0.25) +
  scale_fill_binned(breaks = c(0, 10, 1000, 50000, 200000, 500000, 1000000, 1200000)) +
  binned_scale(aesthetics = "fill", scale_name = "custom", 
               palette = ggplot2:::binned_pal(scales::manual_pal(values = c("#F4F4F4" ,"#FFCCCC", "#FF9999", "#FF6666", 
                                                                            "#FF3333", "#FF0000", "#CC0000"))),
               guide = "bins",
               breaks = c(0, 10, 1000, 50000, 200000, 500000, 1000000, 1200000), limits = c(0, 1200000), show.limits = T) +
  ylab("Feature") +
  geom_vline(xintercept = 26.5, linetype = 2, color = "black") +
  geom_vline(xintercept = 54.5, linetype = 2, color = "black") +
  geom_hline(yintercept = 11.5, linetype = 1, linewidth = 4, color = "white") +
  geom_hline(yintercept = 8.5, linetype = 1, linewidth = 4,color = "white") +
  geom_hline(yintercept = 11.5, linetype = 2, linewidth = 1, color = "lightgrey") +
  geom_hline(yintercept = 8.5, linetype = 2, linewidth = 1,color = "lightgrey") +
  theme(legend.title = element_text(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "right") +
  theme(plot.margin = unit(c(1.5,1,5,5), "lines")) +
  annotation_custom(text_mo,xmin=3.8,xmax=3.8,ymin=19,ymax=19) + 
  annotation_custom(text_pc,xmin=32.5,xmax=32.5,ymin=19,ymax=19) +
  annotation_custom(text_ps,xmin=61,xmax=61,ymin=19,ymax=19) +
  annotation_custom(line3,xmin = 2, xmax = 2, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_mdi,xmin = 2, xmax = 2, ymin = -0.2, ymax = -0.2) +
  annotation_custom(line4,xmin = 5.5, xmax = 5.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_pve,xmin = 5.5, xmax = 5.5, ymin = -0.2, ymax = -0.2) +
  annotation_custom(line3,xmin = 9, xmax = 9, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_cau,xmin = 9, xmax = 9, ymin = -0.2, ymax = -0.2) +
  annotation_custom(line4,xmin = 12.5, xmax = 12.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_pey,xmin = 12.5, xmax = 12.5, ymin = -0.2, ymax = -0.2) +
  annotation_custom(line4,xmin = 16.5, xmax = 16.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_sin,xmin = 16.5, xmax = 16.5, ymin = -0.2, ymax = -0.2) +
  annotation_custom(line4,xmin = 20.5, xmax = 20.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_xen,xmin = 20.5, xmax = 20.5, ymin = -0.2, ymax = -0.2) +
  annotation_custom(line4,xmin = 24.5, xmax = 24.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_hcn,xmin = 24.5, xmax = 24.5, ymin = -0.2, ymax = -0.2) +
  annotation_custom(line4,xmin = 28.5, xmax = 28.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_mdi,xmin = 28.5, xmax = 28.5, ymin = -0.2, ymax = -0.2) + 
  annotation_custom(line4,xmin = 32.5, xmax = 32.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_pve,xmin = 32.5, xmax = 32.5, ymin = -0.2, ymax = -0.2) +
  annotation_custom(line4,xmin = 36.5, xmax = 36.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_cau,xmin = 36.5, xmax = 36.5, ymin = -0.2, ymax = -0.2) +
  annotation_custom(line4,xmin = 40.5, xmax = 40.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_pey,xmin = 40.5, xmax = 40.5, ymin = -0.2, ymax = -0.2) + 
  annotation_custom(line4,xmin = 44.5, xmax = 44.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_sin,xmin = 44.5, xmax = 44.5, ymin = -0.2, ymax = -0.2) +
  annotation_custom(line4,xmin = 48.5, xmax = 48.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_xen,xmin = 48.5, xmax = 48.5, ymin = -0.2, ymax = -0.2) +
  annotation_custom(line4,xmin = 52.5, xmax = 52.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_hcn,xmin = 52.5, xmax = 52.5, ymin = -0.2, ymax = -0.2) +
  annotation_custom(line4,xmin = 56.5, xmax = 56.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_mdi,xmin = 56.5, xmax = 56.5, ymin = -0.2, ymax = -0.2) + 
  annotation_custom(line4,xmin = 60.5, xmax = 60.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_pve,xmin = 60.5, xmax = 60.5, ymin = -0.2, ymax = -0.2) +
  annotation_custom(line2,xmin = 63.5, xmax = 63.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_cau,xmin = 63.5, xmax = 63.5, ymin = -0.2, ymax = -0.2) +
  annotation_custom(line4,xmin = 66.5, xmax = 66.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_pey,xmin = 66.5, xmax = 66.5, ymin = -0.2, ymax = -0.2) + 
  annotation_custom(line4,xmin = 70.5, xmax = 70.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_sin,xmin = 70.5, xmax = 70.5, ymin = -0.2, ymax = -0.2) +
  annotation_custom(line3,xmin = 74, xmax = 74, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_xen,xmin = 74, xmax = 74, ymin = -0.2, ymax = -0.2) +
  annotation_custom(line4,xmin = 77.5, xmax = 77.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_hcn,xmin = 77.5, xmax = 77.5, ymin = -0.2, ymax = -0.2) +
  annotation_custom(group1,xmin = -7, xmax = -7, ymin = 17.5, ymax = 17.5) +
  annotation_custom(group2,xmin = -7, xmax = -7, ymin = 10.5, ymax = 10.5) +
  annotation_custom(group3,xmin = -7, xmax = -7, ymin = 7.5, ymax = 7.5) +
  coord_cartesian(clip = "off") 

treat_all

rm(all_nc_multi, Multi, Multi2, sig_peaks, cluster, all_nc_long, treat_all, data_filtered, text_mo, text_pc, text_ps,
   line2, line3, line4, text_cau, text_pey, text_pve, text_mdi, text_sin, text_xen, text_hcn, group1, group2, group3)


# indicator features for species
set.seed(1102007)
all_nc_multi <-data.frame(all_nc[,c(6:343)]) 

cluster<-factor(all_nc$Species)
Multi<-multipatt(all_nc_multi, cluster)
#capture.output(Multi, file = "Peak_species.txt")
summary(Multi, alpha = 0.01) 
Multi2 <- data.frame(Multi$sign)
Multi2 <- subset(Multi2,p.value<=0.01)
Multi2$Peak <- rownames(Multi2)
write_csv2(Multi2, "indicator/Peak_ind_species_deco.csv")

# create heatmap
sig_peaks  <- read_csv2("indicator/Peak_ind_species_deco.csv")
sig_peaks <- sig_peaks %>% mutate(Peak = gsub("X", "", Peak))
#bring peaks in one column and intensities in another using pivot_longer
all_nc_long <- all_nc %>% pivot_longer(!Fragment_ID : Treatment, names_to = "Peak", values_to = "intensity")
# delete all Peaks in all_nc_long that do not apear in sig_peaks
data_filtered <- subset(all_nc_long, Peak %in% sig_peaks$Peak)
names(data_filtered)[7]<-paste("peak_area")

data_filtered <- data_filtered %>% 
  mutate(Fragment_ID = gsub("HcnMo4", "HcnMo_4", Fragment_ID))

data_filtered$Peak<- factor(data_filtered$Peak, levels = c("264.3s245.138mz", "294.9s259.154mz", "434.3s387.201mz", "470.5s445.243mz",
                                                           "538.9s373.185mz", "539.3s390.212mz", "567.4s387.201mz", "567.4s404.227mz",
                                                           "567.4s409.183mz", "587.1s401.216mz", "587.0s418.243mz", "586.9s423.198mz",
                                                           "688.2s129.055mz", "688.2s147.065mz", "688.4s281.172mz", "731.4s299.221mz",
                                                           "795.9s445.279mz", "796.0s467.261mz", "811.7s459.295mz",
                                                           "408.0s373.185mz", "688.2s259.190mz",
                                                           "1024.8s310.310mz", "539.8s395.167mz", "872.9s482.404mz",
                                                           "408.3s329.174mz", "512.9s329.174mz", "549.4s329.174mz",
                                                           "211.0s388.254mz", "211.1s371.228mz", "229.6s415.254mz",
                                                           "229.6s432.280mz", "246.0s459.280mz", "464.2s149.023mz",
                                                           "689.4s331.188mz", "689.5s165.055mz"))

text_cau <- textGrob("Caulerpa sp.", gp=gpar(fontsize=8.5, fontface="italic"))
text_hcn <- textGrob("H. cnidata", gp=gpar(fontsize=8.5, fontface="italic"))
text_mdi <- textGrob("M. digitata", gp=gpar(fontsize=8.5, fontface="italic"))
text_pey <- textGrob("Peyssonnelia sp.", gp=gpar(fontsize=8.5, fontface="italic"))
text_pve <- textGrob("P. verrucosa", gp=gpar(fontsize=8.5, fontface="italic"))
text_sin <- textGrob("Sinularia sp.", gp=gpar(fontsize=8.5, fontface="italic"))
text_xen <- textGrob("Xenia sp.", gp=gpar(fontsize=8.5, fontface="italic"))

line2 <-textGrob("___", gp=gpar(fontsize=8, fontface="plain"))
line3 <-textGrob("____", gp=gpar(fontsize=8, fontface="plain"))
line4 <-textGrob("_____", gp=gpar(fontsize=8, fontface="plain"))

text_pc <- textGrob("Controlled Polyculture", gp=gpar(fontsize=7.5, fontface="plain"), rot = 45, just = 1)
text_mo <- textGrob("Monoculture", gp=gpar(fontsize=7.5, fontface="plain"), rot = 45, just = 1)
text_ps <- textGrob("Seminatural Polyculture", gp=gpar(fontsize=7.5, fontface="plain"), rot = 45, just = 1)

group1 <-textGrob("Cau", gp=gpar(fontsize=7.5, fontface="plain", color = "lightgrey"), just = 1)
group2 <-textGrob("Xen", gp=gpar(fontsize=7.5, fontface="plain", color = "lightgrey"), just = 1)
group3 <-textGrob("Cau, Mdi, Sin", gp=gpar(fontsize=7.5, fontface="plain", color = "lightgrey"), just = 1)
group4 <-textGrob("Hcn, Pey, Pve, Sin", gp=gpar(fontsize=7.5, fontface="plain", color = "lightgrey"), just = 1)
group5 <-textGrob("Cau, Hcn, Pey, Pve, Xen", gp=gpar(fontsize=7.5, fontface="plain", color = "lightgrey"), just = 1)
group6 <-textGrob("Hcn, Mdi, Pey, Pve, Sin", gp=gpar(fontsize=7.5, fontface="plain", color = "lightgrey"), just = 1)
group7 <-textGrob("Hcn, Mdi, Pey, Pve,
                  Sin, Xen", gp=gpar(fontsize=7.5, fontface="plain", color = "lightgrey"), just = 1)

spec_all <- ggplot(data_filtered, aes(x = Fragment_ID, y = Peak)) +
  geom_tile(colour="lightgrey", linewidth=0.25, aes(fill = peak_area)) +
  ylab("Feature") +
  scale_fill_binned(breaks = c(0, 10, 1000, 50000, 200000, 500000, 1000000, 1500000, 2000000)) +
  binned_scale(aesthetics = "fill", scale_name = "custom", 
               palette = ggplot2:::binned_pal(scales::manual_pal(values = c("#F4F4F4", "#FFCCCC", "#FF9999", "#FF6666", 
                                                                            "#FF3333", "#FF0000", "#CC0000", "#990000", "#660000"))),
               guide = "bins",
               breaks = c(0, 10, 1000, 50000, 200000, 500000, 1000000, 1500000, 2000000), limits = c(0, 2000000), show.limits = T) +
  geom_vline(xintercept = 9.5, linetype = 2, color = "black") +
  geom_vline(xintercept = 21.5, linetype = 2, color = "black") +
  geom_vline(xintercept = 32.5, linetype = 2, color = "black") +
  geom_vline(xintercept = 44.5, linetype = 2, color = "black") +
  geom_vline(xintercept = 56.5, linetype = 2, color = "black") +
  geom_vline(xintercept = 68.5, linetype = 2, color = "black") +
  geom_hline(yintercept = 27.5, linetype = 1, linewidth = 2, color = "white") +
  geom_hline(yintercept = 24.5, linetype = 1, linewidth = 2, color = "white") +
  geom_hline(yintercept = 23.5, linetype = 1, linewidth = 2, color = "white") +
  geom_hline(yintercept = 22.5, linetype = 1, linewidth = 2, color = "white") +
  geom_hline(yintercept = 21.5, linetype = 1, linewidth = 2, color = "white") +
  geom_hline(yintercept = 19.5, linetype = 1, linewidth = 2, color = "white") +
  geom_hline(yintercept = 27.5, linetype = 2, linewidth = 0.5, color = "lightgrey") +
  geom_hline(yintercept = 24.5, linetype = 2, linewidth = 0.5, color = "lightgrey") +
  geom_hline(yintercept = 23.5, linetype = 2, linewidth = 0.5, color = "lightgrey") +
  geom_hline(yintercept = 22.5, linetype = 2, linewidth = 0.5, color = "lightgrey") +
  geom_hline(yintercept = 21.5, linetype = 2, linewidth = 0.5, color = "lightgrey") +
  geom_hline(yintercept = 19.5, linetype = 2, linewidth = 0.5, color = "lightgrey") +
  theme(legend.title = element_text(),
      legend.position = "right",
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank()) +
  theme(plot.margin = unit(c(1.5,0.5,5,5.5), "lines")) +
  annotation_custom(text_cau,xmin=4,xmax=4,ymin=36,ymax=36) + 
  annotation_custom(text_hcn,xmin=12,xmax=12,ymin=36,ymax=36) +
  annotation_custom(text_mdi,xmin=24,xmax=24,ymin=36,ymax=36) +
  annotation_custom(text_pey,xmin=37,xmax=37,ymin=36,ymax=36) +
  annotation_custom(text_pve,xmin=48,xmax=48,ymin=36,ymax=36) +
  annotation_custom(text_sin,xmin=60,xmax=60,ymin=36,ymax=36) +
  annotation_custom(text_xen,xmin=71,xmax=71,ymin=36,ymax=36) +
  annotation_custom(line3,xmin = 2, xmax = 2, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_mo,xmin = 2, xmax = 2, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line4,xmin = 5.5, xmax = 5.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_pc,xmin = 5.5, xmax = 5.5, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line2,xmin = 8.5, xmax = 8.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_ps,xmin = 8.5, xmax = 8.5, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line4,xmin = 11.5, xmax = 11.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_mo,xmin = 11.5, xmax = 11.5, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line4,xmin = 15.5, xmax = 15.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_pc,xmin = 15.5, xmax = 15.5, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line4,xmin = 19.5, xmax = 19.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_ps,xmin = 19.5, xmax = 19.5, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line3,xmin = 23, xmax = 23, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_mo,xmin = 23, xmax = 23, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line4,xmin = 26.5, xmax = 26.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_pc,xmin = 26.5, xmax = 26.5, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line4,xmin = 30.5, xmax = 30.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_ps,xmin = 30.5, xmax = 30.5, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line4,xmin = 34.5, xmax = 34.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_mo,xmin = 34.5, xmax = 34.5, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line4,xmin = 38.5, xmax = 38.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_pc,xmin = 38.5, xmax = 38.5, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line4,xmin = 42.5, xmax = 42.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_ps,xmin = 42.5, xmax = 42.5, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line4,xmin = 46.5, xmax = 46.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_mo,xmin = 46.5, xmax = 46.5, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line4,xmin = 50.5, xmax = 50.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_pc,xmin = 50.5, xmax = 50.5, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line4,xmin = 54.5, xmax = 54.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_ps,xmin = 54.5, xmax = 54.5, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line4,xmin = 58.5, xmax = 58.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_mo,xmin = 58.5, xmax = 58.5, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line4,xmin = 62.5, xmax = 62.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_pc,xmin = 62.5, xmax = 62.5, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line4,xmin = 66.5, xmax = 66.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_ps,xmin = 66.5, xmax = 66.5, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line4,xmin = 70.5, xmax = 70.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_mo,xmin = 70.5, xmax = 70.5, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line4,xmin = 74.5, xmax = 74.5, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_pc,xmin = 74.5, xmax = 74.5, ymin = -0.5, ymax = -0.5) +
  annotation_custom(line3,xmin = 78, xmax = 78, ymin = 0.2, ymax = 0.2) +
  annotation_custom(text_ps,xmin = 78, xmax = 78, ymin = -0.5, ymax = -0.5) +
  annotation_custom(group1,xmin = -0, xmax = -0, ymin = 35, ymax = 35) +
  annotation_custom(group2,xmin = -0, xmax = -0, ymin = 27, ymax = 27) +
  annotation_custom(group3,xmin = -0, xmax = -0, ymin = 24, ymax = 24) +
  annotation_custom(group4,xmin = -0, xmax = -0, ymin = 23, ymax = 23) +
  annotation_custom(group5,xmin = -0, xmax = -0, ymin = 22, ymax = 22) +
  annotation_custom(group6,xmin = -0, xmax = -0, ymin = 21, ymax = 21) +
  annotation_custom(group7,xmin = -1.5, xmax = -1.5, ymin = 18.5, ymax = 18.5) +
  coord_cartesian(clip = "off") 

spec_all

rm(all_nc_multi, Multi, Multi2, sig_peaks, cluster, all_nc_long, text_cau, text_hcn, text_mdi, text_pey,
   text_pve, text_sin, text_xen, spec_all, data_filtered, all_nc, group1, group2, group3, group4, group5, group6, group7,
   line2, line3, line4, text_mo, text_pc, text_ps)
