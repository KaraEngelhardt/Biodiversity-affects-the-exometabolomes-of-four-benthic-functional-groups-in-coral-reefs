# data preparation
# bucket table and curated feature data
# __________________________________________________________
# install packages
library(tidyverse)
library(decontam)
#___________________________________________________________
# M2NT
all <- read_csv2("raw_data/Node_table.csv")

# transpose tibble
transpose_df <- function(df) {
  t_df <- data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  t_df <- t_df %>%
    tibble::rownames_to_column(.data = .) %>%
    tibble::as_tibble(.)
  return(t_df)}

all <- transpose_df(all)

# include header 
names(all) <- all[1,]
all <- all[-1,]

# add metadata
# Species
all$Species = all$OTU_ID

all <- all %>% 
  mutate(Species = gsub("PcPcMo_1", "", Species)) %>%
  mutate(Species = gsub("PcPcMo_2", "", Species)) %>%
  mutate(Species = gsub("PcPcMo_3", "", Species)) %>%
  mutate(Species = gsub("PcPcMo_4", "", Species)) %>%
  mutate(Species = gsub("PcPcPc_1", "", Species)) %>%
  mutate(Species = gsub("PcPcPc_2", "", Species)) %>%
  mutate(Species = gsub("PcPcPc_3", "", Species)) %>%
  mutate(Species = gsub("PcPcPc_4", "", Species)) %>%
  mutate(Species = gsub("PcPcPs_1", "", Species)) %>%
  mutate(Species = gsub("PcPcPs_2", "", Species)) %>%
  mutate(Species = gsub("PcPcPs_3", "", Species)) %>%
  mutate(Species = gsub("PcPcPs_4", "", Species)) %>%
  mutate(Species = gsub("PcPcMo4", "", Species)) %>%
  mutate(Species = gsub("PcPcPc_5", "", Species)) %>%
  mutate(Species = gsub("_d1_1", "", Species)) %>%
  mutate(Species = gsub("_d1_2", "", Species)) %>%
  mutate(Species = gsub("_d1_3", "", Species)) %>%
  mutate(Species = gsub("_d2_1", "", Species)) %>%
  mutate(Species = gsub("_d2_2", "", Species)) %>%
  mutate(Species = gsub("_d2_3", "", Species)) %>%
  mutate(Species = gsub("_d3_1", "", Species)) %>%
  mutate(Species = gsub("_d3_2", "", Species)) %>%
  mutate(Species = gsub("_d3_3", "", Species)) %>%
  mutate(Species = gsub("_d4_1", "", Species)) %>%
  mutate(Species = gsub("_d4_2", "", Species)) %>%
  mutate(Species = gsub("_d4_3", "", Species))

all <- all %>%
  mutate(Species = gsub("Cau", "Caulerpa sp.",Species)) %>%
  mutate(Species = gsub("Pey", "Peyssonnelia sp.", Species)) %>%
  mutate(Species = gsub("Mdi", "Montipora digitata", Species)) %>%
  mutate(Species = gsub("Pve", "Pocillopora verrucosa", Species)) %>%
  mutate(Species = gsub("Xen", "Xenia sp.", Species)) %>%
  mutate(Species = gsub("Sin", "Sinularia sp.", Species)) %>%
  mutate(Species = gsub("Hcn", "Haliclona cnidata", Species))

# treatment
all$Treatment = all$OTU_ID

all <- all%>%
  mutate(Treatment = gsub("CauPcPc", "", Treatment)) %>%
  mutate(Treatment = gsub("PeyPcPc", "", Treatment)) %>%
  mutate(Treatment = gsub("SinPcPc", "", Treatment)) %>%
  mutate(Treatment = gsub("XenPcPc", "", Treatment)) %>%
  mutate(Treatment = gsub("MdiPcPc", "", Treatment)) %>%
  mutate(Treatment = gsub("PvePcPc", "", Treatment)) %>%
  mutate(Treatment = gsub("HcnPcPc", "", Treatment)) %>%
  mutate(Treatment = gsub("_d1_1", "", Treatment)) %>%
  mutate(Treatment = gsub("_d1_2", "", Treatment)) %>%
  mutate(Treatment = gsub("_d1_3", "", Treatment)) %>%
  mutate(Treatment = gsub("_d2_1", "", Treatment)) %>%
  mutate(Treatment = gsub("_d2_2", "", Treatment)) %>%
  mutate(Treatment = gsub("_d2_3", "", Treatment)) %>%
  mutate(Treatment = gsub("_d3_1", "", Treatment)) %>%
  mutate(Treatment = gsub("_d3_2", "", Treatment)) %>%
  mutate(Treatment = gsub("_d3_3", "", Treatment)) %>%
  mutate(Treatment = gsub("_d4_1", "", Treatment)) %>%
  mutate(Treatment = gsub("_d4_2", "", Treatment)) %>%
  mutate(Treatment = gsub("_d4_3", "", Treatment))

all <- all%>%
  mutate(Treatment = gsub("_1", "", Treatment)) %>%
  mutate(Treatment = gsub("_2", "", Treatment)) %>%
  mutate(Treatment = gsub("_3", "", Treatment)) %>%
  mutate(Treatment = gsub("_4", "", Treatment)) %>%
  mutate(Treatment = gsub("_5", "", Treatment)) 

all <- all%>%
  mutate(Treatment = gsub("Pc", "controlled_Polyculture", Treatment)) %>%
  mutate(Treatment = gsub("Mo", "Monoculture", Treatment)) %>%
  mutate(Treatment = gsub("Ps", "seminatural_Polyculture", Treatment))

# organism group
all$Organism_group = all$Species

all <- all%>%
  mutate(Organism_group = gsub("Xenia sp.", "Soft coral", Organism_group)) %>%
  mutate(Organism_group = gsub("Sinularia sp.", "Soft coral", Organism_group)) %>%
  mutate(Organism_group = gsub("Peyssonnelia sp.", "Macroalgae", Organism_group)) %>%
  mutate(Organism_group = gsub("Caulerpa sp.", "Macroalgae", Organism_group)) %>%
  mutate(Organism_group = gsub("Montipora digitata", "Stony coral", Organism_group)) %>%
  mutate(Organism_group = gsub("Pocillopora verrucosa", "Stony coral", Organism_group)) %>%
  mutate(Organism_group = gsub("Haliclona cnidata", "Sponge", Organism_group))

#add day infomration
day <- read_csv2("raw_data/Day_information.csv")
day <- transpose_df(day)
names(day) <- day[1,]
day <- day[-1,]
all <- dplyr::rename(all, c("Fragment_ID" = "OTU_ID"))
all <- merge(x = all, y = day, by = "Fragment_ID")

#change order of columns
all <- all %>% 
  relocate(Species, .after = Fragment_ID)%>%
  relocate(Treatment, .after = Species) %>%
  relocate(Organism_group, .after = Species) %>%
  relocate(Day, .after = Fragment_ID)

#save prep data as csv.
write.csv(all, file = "data_prep/data_metabolites_nodes_prep.csv")
rm(all, day)
#___________________________________________________________
#___________________________________________________________
# CFT
#load in table
bolit <- read_csv2("raw_data/data_metabolites_curated.csv")

bolit <- transpose_df(bolit) # function from above
bolit <- bolit %>% select(-c(rowname))
names(bolit) <- bolit[1,]
bolit <- bolit[-1,]
str(bolit) 
bolit[, c(2:339)] <- sapply(bolit[, c(2:339)], as.numeric)
bolit <- bolit %>% mutate(Fragment_ID = gsub("PcPc", "", Fragment_ID))

# outlier: MKZE001_B-06, MKZE001_A-06, MKZE001_A-01 -> this corresponds to:CauPcPcMo_1, MdiPcPcMo_3, XenPcPcPs_3
# remove methanol controls and outlier (detected by Christoph Hartwig)
bolit <- bolit %>% 
  filter(Fragment_ID != "Methanol_control") %>%
  filter(Fragment_ID != "CauMo_1") %>%
  filter(Fragment_ID != "XenPs_3") %>%
  filter(Fragment_ID != "MdiMo_3")

#add metadata
#Species
bolit$Species = bolit$Fragment_ID

bolit <- bolit %>% 
  mutate(Species = gsub("Mo_1", "", Species)) %>%
  mutate(Species = gsub("Mo_2", "", Species)) %>%
  mutate(Species = gsub("Mo_3", "", Species)) %>%
  mutate(Species = gsub("Mo_4", "", Species)) %>%
  mutate(Species = gsub("Pc_1", "", Species)) %>%
  mutate(Species = gsub("Pc_2", "", Species)) %>%
  mutate(Species = gsub("Pc_3", "", Species)) %>%
  mutate(Species = gsub("Pc_4", "", Species)) %>%
  mutate(Species = gsub("Ps_1", "", Species)) %>%
  mutate(Species = gsub("Ps_2", "", Species)) %>%
  mutate(Species = gsub("Ps_3", "", Species)) %>%
  mutate(Species = gsub("Ps_4", "", Species)) %>%
  mutate(Species = gsub("Mo4", "", Species)) %>%
  mutate(Species = gsub("Pc_5", "", Species))

bolit <- bolit %>%
  mutate(Species = gsub("Cau", "Caulerpa sp.", Species)) %>%
  mutate(Species = gsub("Pey", "Peyssonnelia sp.", Species)) %>%
  mutate(Species = gsub("Mdi", "Montipora digitata", Species)) %>%
  mutate(Species = gsub("Pve", "Pocillopora verrucosa", Species)) %>%
  mutate(Species = gsub("Xen", "Xenia sp.", Species)) %>%
  mutate(Species = gsub("Sin", "Sinularia sp.", Species)) %>%
  mutate(Species = gsub("Hcn", "Haliclona cnidata", Species))

#organism_group
bolit$Organism_group = bolit$Species

bolit <- bolit %>%
  mutate(Organism_group = gsub("Caulerpa sp.", "Macroalgae", Organism_group)) %>%
  mutate(Organism_group = gsub("Peyssonnelia sp.", "Macroalgae", Organism_group)) %>%
  mutate(Organism_group = gsub("Montipora digitata", "Stony coral", Organism_group)) %>%
  mutate(Organism_group = gsub("Pocillopora verrucosa", "Stony coral", Organism_group)) %>%
  mutate(Organism_group = gsub("Xenia sp.", "Soft coral", Organism_group)) %>%
  mutate(Organism_group = gsub("Sinularia sp.", "Soft coral", Organism_group)) %>%
  mutate(Organism_group = gsub("Haliclona cnidata", "Sponge", Organism_group))

#treatment
bolit$Treatment = bolit$Fragment_ID

bolit <- bolit%>%
  mutate(Treatment = gsub("Cau", "", Treatment)) %>%
  mutate(Treatment = gsub("Pey", "", Treatment)) %>%
  mutate(Treatment = gsub("Sin", "", Treatment)) %>%
  mutate(Treatment = gsub("Xen", "", Treatment)) %>%
  mutate(Treatment = gsub("Mdi", "", Treatment)) %>%
  mutate(Treatment = gsub("Pve", "", Treatment)) %>%
  mutate(Treatment = gsub("Hcn", "", Treatment))

bolit <- bolit%>%
  mutate(Treatment = gsub("_1", "", Treatment)) %>%
  mutate(Treatment = gsub("_2", "", Treatment)) %>%
  mutate(Treatment = gsub("_3", "", Treatment)) %>%
  mutate(Treatment = gsub("_4", "", Treatment)) %>%
  mutate(Treatment = gsub("_5", "", Treatment)) %>%
  mutate(Treatment = gsub("4", "", Treatment))

bolit <- bolit%>%
  mutate(Treatment = gsub("Pc", "controlled_Polyculture", Treatment)) %>%
  mutate(Treatment = gsub("Mo", "Monoculture", Treatment)) %>%
  mutate(Treatment = gsub("Ps", "seminatural_Polyculture", Treatment))

#change order of columns
bolit <- bolit %>% 
  relocate(Species, .after = Fragment_ID)%>%
  relocate(Organism_group, .after = Species)%>%
  relocate(Treatment, .after = Organism_group)

# apply decontam
data_all_days <- read_csv2("raw_data/metadata_features.csv")
day1 <- data_all_days %>% subset(Day == "1")
day2 <- data_all_days %>% subset(Day == "2")
day3 <- data_all_days %>% subset(Day == "3")
day4 <- data_all_days %>% subset(Day == "4")
rm(data_all_days)

# decontam day 1
Control <- day1$is.neg <- day1$Sample_or_Control == "TRUE Sample"
day1 <- day1[, -1]
day1 <- day1[, -1]
day1 <- day1[, -1]
day1 <- day1[,-339]
str(day1)
day1 <- as.matrix(day1)
cont_prev_day1 <- isContaminant(day1, method="prevalence", neg=Control)
table(cont_prev_day1$contaminant)
# 13 contaminants
rm(Control, day1)

#_____________________________________
# decontam day 2
Control <- day2$is.neg <- day2$Sample_or_Control == "TRUE Sample"
day2 <- day2[, -1]
day2 <- day2[, -1]
day2 <- day2[, -1]
day2 <- day2[,-339]
str(day2)
day2 <- as.matrix(day2)
cont_prev_day2 <- isContaminant(day2, method="prevalence", neg=Control)
table(cont_prev_day2$contaminant)
# 16 contaminants
rm(Control, day2)

#_____________________________________
# decontam day 3
Control <- day3$is.neg <- day3$Sample_or_Control == "TRUE Sample"
day3 <- day3[, -1]
day3 <- day3[, -1]
day3 <- day3[, -1]
day3 <- day3[,-339]
str(day3)
day3 <- as.matrix(day3)
cont_prev_day3 <- isContaminant(day3, method="prevalence", neg=Control)
table(cont_prev_day3$contaminant)
# 10 contaminants
rm(Control, day3)

#_____________________________________
# decontam day 4 
Control <- day4$is.neg <- day4$Sample_or_Control == "TRUE Sample"
day4 <- day4[, -1]
day4 <- day4[, -1]
day4 <- day4[, -1]
day4 <- day4[,-339]
str(day4)
day4 <- as.matrix(day4)
cont_prev_day4 <- isContaminant(day4, method="prevalence", neg=Control)
table(cont_prev_day4$contaminant)
# 24 contaminants
rm(Control, day4)

#______________________________________
# normalise contaminants from bolit
day <- read_csv2("raw_data/Day_information.csv")
day <- transpose_df(day)
names(day) <- day[1,]
day <- day[-1,]
day <- day %>% mutate(Fragment_ID = gsub("PcPc", "", Fragment_ID))
all <- merge(x = bolit, y = day, by = "Fragment_ID")
all <-all %>% relocate(Day, .after = Fragment_ID)
#split into days
day1 <- all %>% subset(Day == "1")
day2 <- all %>% subset(Day == "2")
day3 <- all %>% subset(Day == "3")
day4 <- all %>% subset(Day == "4")
rm(day, bolit, transpose_df, all)

#__________________________________
# remove contaminants day 1
cont_prev_day1 <- rownames_to_column(cont_prev_day1, var = "Features")
day1_long <- day1 %>% 
  pivot_longer(!Fragment_ID : Treatment, names_to = "Features", values_to = "intensity")

merged_day1 <- left_join(day1_long, cont_prev_day1, by = c("Features" = "Features"))

merged_day1$intensity <- ifelse(merged_day1$contaminant, 0, merged_day1$intensity)

# remove controls (no longer needed since their are now normalized)
merged_day1 <- merged_day1 %>%
  filter(Treatment !="control_d1")

merged_day1 <- merged_day1[,c(1,2,3,4,5,6,7)]

day1_norm <- merged_day1 %>%
  pivot_wider(names_from = Features, values_from = intensity)

rm(day1, day1_long, merged_day1, cont_prev_day1)

#__________________________________
# remove contaminants day 2
cont_prev_day2 <- rownames_to_column(cont_prev_day2, var = "Features")

day2_long <- day2 %>% 
  pivot_longer(!Fragment_ID : Treatment, names_to = "Features", values_to = "intensity")

merged_day2 <- left_join(day2_long, cont_prev_day2, by = c("Features" = "Features"))

merged_day2$intensity <- ifelse(merged_day2$contaminant, 0, merged_day2$intensity)

# remove controls (no longer needed since their are now normalized)
merged_day2 <- merged_day2 %>%
  filter(Treatment !="control_d2")

merged_day2 <- merged_day2[,c(1,2,3,4,5,6,7)]

day2_norm <- merged_day2 %>%
  pivot_wider(names_from = Features, values_from = intensity)

rm(day2, day2_long, merged_day2, cont_prev_day2)

#__________________________________
# remove contaminants day 3
cont_prev_day3 <- rownames_to_column(cont_prev_day3, var = "Features")

day3_long <- day3 %>% 
  pivot_longer(!Fragment_ID : Treatment, names_to = "Features", values_to = "intensity")

merged_day3 <- left_join(day3_long, cont_prev_day3, by = c("Features" = "Features"))

merged_day3$intensity <- ifelse(merged_day3$contaminant, 0, merged_day3$intensity)

# remove controls (no longer needed since their are now normalized)
merged_day3 <- merged_day3 %>%
  filter(Treatment !="control_d3")

merged_day3 <- merged_day3[,c(1,2,3,4,5,6,7)]

day3_norm <- merged_day3 %>%
  pivot_wider(names_from = Features, values_from = intensity)

rm(day3, day3_long, merged_day3, cont_prev_day3)

#__________________________________
# remove contaminants day 4
cont_prev_day4 <- rownames_to_column(cont_prev_day4, var = "Features")

day4_long <- day4 %>% 
  pivot_longer(!Fragment_ID : Treatment, names_to = "Features", values_to = "intensity")

merged_day4 <- left_join(day4_long, cont_prev_day4, by = c("Features" = "Features"))

merged_day4$intensity <- ifelse(merged_day4$contaminant, 0, merged_day4$intensity)

merged_day4 <- merged_day4 %>%
  filter(Treatment !="control_d")

merged_day4 <- merged_day4[,c(1,2,3,4,5,6,7)]

day4_norm <- merged_day4 %>%
  pivot_wider(names_from = Features, values_from = intensity)

rm(day4, day4_long, merged_day4, cont_prev_day4)

# rebind to have a complete normalized table
all_deco <- rbind(day1_norm, day2_norm, day3_norm, day4_norm)
rm(day1_norm, day2_norm, day3_norm, day4_norm)
write.csv(all_deco, file = "data_prep/all_decontaminated.csv")
rm(all_deco)
