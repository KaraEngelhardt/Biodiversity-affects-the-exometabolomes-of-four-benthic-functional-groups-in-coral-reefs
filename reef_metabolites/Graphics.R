# Graphs for glmer visualisation
# code for Fig. 2, 3 & 4 (in manuscript) as well as Fig. 28, 29 & 30

#____________________________________________________________________
# install packages
library(vegan)
library(ggplot2)
library(tidyverse)
library(ggpubr)
#____________________________________________________________________
# M2NT
# Manuscript Fig. 2
all <- read_csv("data_prep/data_metabolites_nodes_prep.csv")
all <- all %>% select(-c(...1))
all_matrix <- select(all, -c(Day, Species, Treatment, Organism_group))
all_matrix <- all_matrix%>%as.data.frame
rownames(all_matrix) <- all_matrix$Fragment_ID
all_matrix <- all_matrix[,-1]
all_matrix <- as.matrix(all_matrix)

dist_nodes <- vegdist(all_matrix, method = "jaccard", binary = TRUE)
dist_nodes
nmds_treat <- metaMDS(dist_nodes)
nmds_treat
sppscores(nmds_treat) <- dist_nodes

data_nmds_treat <- as.data.frame(scores(nmds_treat)) 
data_nmds_treat
data_nmds_treat <- rownames_to_column(data_nmds_treat, var = "fragment_id")
data_nmds_treat <- as_tibble(data_nmds_treat)
data_nmds_treat$Treatment = all$Treatment
data_nmds_treat$Species = all$Species

# prepare for plot
data_nmds_treat$Treatment<- factor(data_nmds_treat$Treatment, levels = c("Monoculture",
                                                                         "controlled_Polyculture", 
                                                                         "seminatural_Polyculture",
                                                                         "control"))

centroid_treat <- data_nmds_treat%>%
  group_by(Treatment) %>%
  summarize(sites.NMDS1=mean(sites.NMDS1), sites.NMDS2=mean(sites.NMDS2))

treat <- ggplot(data_nmds_treat, aes(x=sites.NMDS1, y=sites.NMDS2, color = Treatment)) +
  geom_point(aes(colour=Treatment, fill=Treatment), size = 3, shape = 21, color = "black",
             position=position_jitter(0.01)) +
  geom_point(data = centroid_treat, size = 3, shape = 4, stroke = 2, show_guide = FALSE) +
  scale_fill_manual(values=c("#99FFFF", "#3333FF", "#000033", "#999999")) +
  scale_color_manual(values=c("#99FFFF", "#3333FF", "#000033", "#999999")) +
  stat_ellipse(aes(color = Treatment)) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 8),
        legend.position = "bottom") +
  annotate(geom = "text", x = -0.8, y = 0.35, size = 3,
           label = "Treatment
           2D Stress = 0.19",
           hjust = 0) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE))

treat
rm(all_matrix, all_treat, centroid_treat, nmds_treat, dist_nodes, groups_treat)

#______________________________________________________
# species analysis
all_matrix <- select(all, -c(Day, Species, Treatment, Organism_group))

all_matrix <- all_matrix%>%as.data.frame

rownames(all_matrix) <- all_matrix$Fragment_ID
all_matrix <- all_matrix[,-1]
all_matrix <- as.matrix(all_matrix)

dist_nodes <- vegdist(all_matrix, method = "jaccard", binary = TRUE)
dist_nodes
nmds_spec <- metaMDS(dist_nodes)
nmds_spec
sppscores(nmds_spec) <- dist_nodes

data_nmds_spec <- as.data.frame(scores(nmds_spec)) 
data_nmds_spec 
data_nmds_spec <- rownames_to_column(data_nmds_spec, var = "fragment_id")
data_nmds_spec <- as_tibble(data_nmds_spec)
data_nmds_spec$Treatment = all$Treatment
data_nmds_spec$Species = all$Species

# prepare for plot
data_nmds_spec$Species<- factor(data_nmds_spec$Species, levels = c("Montipora digitata", "Pocillopora verrucosa",
                                                                   "Xenia sp.", "Sinularia sp.",
                                                                   "Caulerpa sp.","Peyssonnelia sp.",
                                                                   "Haliclona cnidata", "control"))

centroid_spec <- data_nmds_spec%>%
  group_by(Species) %>%
  summarize(sites.NMDS1=mean(sites.NMDS1), sites.NMDS2=mean(sites.NMDS2))

species <- ggplot(data_nmds_spec, aes(x=sites.NMDS1, y=sites.NMDS2, color = Species)) +
  geom_point(aes(colour=Species, fill=Species), size = 3, shape = 21, color = "black",
             position=position_jitter(0.01)) +
  geom_point(data = centroid_spec, size = 3, shape = 4, stroke = 2, show_guide = FALSE) +
  scale_fill_manual(values=c("#FFCC00", "#FFFF33", "#FF3399", "#FF99CC", "#009933", "#00FF00", "#3399FF", "#999999")) +
  scale_colour_manual(values=c( "#FFCC00", "#FFFF33","#FF3399", "#FF99CC",  "#009933", "#00FF00", "#3399FF", "#999999")) +
  stat_ellipse(aes(color = Species)) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 8),
        legend.position = "bottom") +
  annotate(geom = "text", x = -0.8, y = 0.45, size = 3,
           label = "Species
           2D Stress = 0.19",
           hjust = 0) +
  guides(fill=guide_legend(nrow=4, byrow=TRUE))

species
rm(all_matrix, all_spec, centroid_spec, nmds_spec, dist_nodes, groups_species)

#______________________________________________________
# organism group analysis

all_matrix <- select(all, -c(Day, Species, Treatment, Organism_group))
all_matrix <- all_matrix%>%as.data.frame
rownames(all_matrix) <- all_matrix$Fragment_ID
all_matrix <- all_matrix[,-1]
all_matrix <- as.matrix(all_matrix)
dist_nodes <- vegdist(all_matrix, method = "jaccard", binary = TRUE)
dist_nodes
nmds_org <- metaMDS(dist_nodes)
nmds_org
sppscores(nmds_org) <- dist_nodes

data_nmds_org <- as.data.frame(scores(nmds_org)) 
data_nmds_org 
data_nmds_org <- rownames_to_column(data_nmds_org, var = "fragment_id")
data_nmds_org <- as_tibble(data_nmds_org)
data_nmds_org$Treatment = all$Treatment
data_nmds_org$Species = all$Species
data_nmds_org$Organism_group = all$Organism_group

#prepare for plot
data_nmds_org$Organism_group<- factor(data_nmds_org$Organism_group, levels = c("Stony coral", "Soft coral", 
                                                                               "Macroalgae", "Sponge", "control"))

centroid_org <- data_nmds_org %>%
  group_by(Organism_group) %>%
  summarize(sites.NMDS1=mean(sites.NMDS1), sites.NMDS2=mean(sites.NMDS2))

Org <- ggplot(data_nmds_org, aes(x=sites.NMDS1, y=sites.NMDS2, colour = Organism_group)) +
  geom_point(aes(colour=Organism_group, fill=Organism_group), size = 3, shape = 21, color = "black",
             position=position_jitter(0.01)) +
  geom_point(data = centroid_org, size = 3, shape = 4, stroke = 2, show_guide = FALSE) +
  scale_fill_manual(values=c("#FF9900", "#FF3399", "#009933", "#3399FF", "#999999")) +
  scale_colour_manual(values=c("#FF9900", "#FF3399", "#009933", "#3399FF", "#999999")) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  stat_ellipse() +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 8),
        legend.position = "bottom") +
  guides(color = guide_legend(nrow = 3)) +
  annotate(geom = "text", x = -0.8, y = 0.42, size = 3,
           label = "Organism group
           2D Stress = 0.19",
           hjust = 0) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE))

Org
rm(all_matrix, centroid_org, nmds_org, dist_nodes)

#_________
# merge
arranged_overview_nodes <- ggarrange(treat, Org, species, 
                                     labels = c("a)", "b)", "c)"),
                                     align = "hv",
                                     label.x = 0,
                                     label.y = 0.95,
                                     common.legend = FALSE,
                                     # legend = c("bottom"),
                                     ncol = 3, nrow = 1,
                                     hjust = -0.5,
                                     vjust = 0) 


arranged_overview_nodes
rm(arranged_overview_nodes, species, treat, Org, data_nmds_org, data_nmds_spec, data_nmds_treat)

#________________________________________________________________________________
# within treatments
# Manuscript Fig. 3
mono <- all %>% 
  subset(Treatment != "controlled_Polyculture") %>%
  subset(Treatment != "seminatural_Polyculture")

Pc <- all %>% 
  subset(Treatment !="Monoculture") %>%
  subset(Treatment != "seminatural_Polyculture")

Ps <- all %>% 
  subset(Treatment !="Monoculture")%>%
  subset(Treatment != "controlled_Polyculture")

mono <- mono %>% as.data.frame()
rownames(mono) <- mono$Fragment_ID
mono <- mono %>% select(-c(Fragment_ID, Species, Organism_group, Treatment))

mono_matrix <- as.matrix(mono)
rm(mono)

set.seed(12345)
mono_dist <- vegdist(mono_matrix, method="jaccard", binary = TRUE)
mono_dist
set.seed(12345)
nmds_mono <- metaMDS(mono_dist)
sppscores(nmds_mono) <- mono_dist
nmds_mono

data_nmds_mono <- as.data.frame(scores(nmds_mono))
data_nmds_mono <- rownames_to_column(data_nmds_mono, var = "Fragment_id")

mono <- all %>% 
  subset(Treatment != "controlled_Polyculture") %>%
  subset(Treatment != "seminatural_Polyculture")

data_nmds_mono$Species <- mono$Species
data_nmds_mono$Treatment <- mono$Treatment
data_nmds_mono$Organism_group <- mono$Organism_group

#plot
data_nmds_mono$Organism_group<- factor(data_nmds_mono$Organism_group, levels = c("Stony coral", "Soft coral", 
                                                                                 "Macroalgae", "Sponge", "control"))

data_nmds_mono$Species<- factor(data_nmds_mono$Species, levels = c("Montipora digitata", "Pocillopora verrucosa", 
                                                                   "Xenia sp.", "Sinularia sp.",
                                                                   "Caulerpa sp.", "Peyssonnelia sp.",
                                                                   "Haliclona cnidata", "control"))

mono_species <- ggplot(data_nmds_mono, aes(x=sites.NMDS1, y=sites.NMDS2, col = Organism_group, shape = Species)) +
  geom_point(size = 3.5, stroke = 1.5) +
  scale_shape_manual(name = "Species", label=c("Montipora digitata" = expression(italic("Montipora digitata")),
                                               "Pocillopora verrucosa" = expression(italic("Pocillopora verrucosa")),
                                               "Xenia sp." = expression(paste(italic("Xenia")," sp.")),
                                               "Sinularia sp." = expression(paste(italic("Sinularia")," sp.")),
                                               "Caulerpa sp." = expression(paste(italic("Caulerpa")," sp.")),
                                               "Peyssonnelia sp." = expression(paste(italic("Peyssonnelia")," sp.")),
                                               "Haliclona cnidata" = expression(italic("Haliclona cnidata")),
                                               "control" = expression("control")),
                     values=c(3, 4, 8, 15, 16, 17, 18, 1)) +
  scale_colour_manual(name = "Organism group", values=c("#FF9900", "#FF3399", "#009933", "#3399FF", "#999999")) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  annotate(geom = "text", x = -0.5, y = 0.35, size = 3,
           label = "Monoculture
           2D Stress = 0.15",
           hjust = 0) +
  theme(plot.title = element_text(size = 10, face = "bold")) +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(face = "bold", size = 8))

mono_species  
rm(mono, mono_matrix, nmds_mono, mono_dist)

#______________________________________
# test for differences within polyculture controlled
poly <- Pc %>% as.data.frame()
rownames(poly) <- poly$Fragment_ID
poly <- poly %>% select(-c(Day, Fragment_ID, Species, Organism_group, Treatment))

poly_matrix <- as.matrix(poly)
rm(poly)

set.seed(12345)
poly_dist <- vegdist(poly_matrix, method="jaccard", binary = TRUE)
poly_dist
set.seed(12345)
nmds_poly <- metaMDS(poly_dist)
sppscores(nmds_poly) <- poly_dist
nmds_poly

data_nmds_poly <- as.data.frame(scores(nmds_poly))
data_nmds_poly <- rownames_to_column(data_nmds_poly, var = "Fragment_id")

data_nmds_poly$Species <- Pc$Species
data_nmds_poly$Treatment <- Pc$Treatment
data_nmds_poly$Organism_group <- Pc$Organism_group

#plot
data_nmds_poly$Organism_group<- factor(data_nmds_poly$Organism_group, levels = c("Stony coral", "Soft coral", 
                                                                                 "Macroalgae", "Sponge", "control"))

data_nmds_poly$Species<- factor(data_nmds_poly$Species, levels = c("Montipora digitata", "Pocillopora verrucosa", 
                                                                   "Xenia sp.", "Sinularia sp.",
                                                                   "Caulerpa sp.", "Peyssonnelia sp.",
                                                                   "Haliclona cnidata", "control"))

poly_species <- ggplot(data_nmds_poly, aes(x=sites.NMDS1, y=sites.NMDS2, col = Organism_group, shape = Species)) +
  geom_point(size = 3.5, stroke = 1.5) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_colour_manual(values=c("#FF9900", "#FF3399", "#009933", "#3399FF", "#999999")) +
  scale_shape_manual(name = "Species", label=c("Montipora digitata" = expression(italic("Montipora digitata")),
                                               "Pocillopora verrucosa" = expression(italic("Pocillopora verrucosa")),
                                               "Xenia sp." = expression(paste(italic("Xenia")," sp.")),
                                               "Sinularia sp." = expression(paste(italic("Sinularia")," sp.")),
                                               "Caulerpa sp." = expression(paste(italic("Caulerpa")," sp.")),
                                               "Peyssonnelia sp." = expression(paste(italic("Peyssonnelia")," sp.")),
                                               "Haliclona cnidata" = expression(italic("Haliclona cnidata")),
                                               "control" = expression("control")),
                     values=c(3, 4, 8, 15, 16, 17, 18, 1)) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  annotate(geom = "text", x = -0.5, y = 0.34, size = 3,
           label = "Controlled Polyculture
           2D Stress = 0.17 ",
           hjust = 0) +
  theme(legend.text = element_text( size = 8),
        legend.title = element_text(face = "bold", size = 8))

poly_species
rm(Pc, poly_matrix, nmds_poly, poly_dist)

#______________________________________
# test within polyculture seminatural
semi <- Ps %>% as.data.frame()
rownames(semi) <- semi$Fragment_ID
semi <- semi %>% select(-c(Day, Fragment_ID, Species, Organism_group, Treatment))

semi_matrix <- as.matrix(semi)
rm(semi)
set.seed(12345)
semi_dist <- vegdist(semi_matrix, method="jaccard", binary = TRUE)
semi_dist
set.seed(12345)
nmds_semi <- metaMDS(semi_dist)
sppscores(nmds_semi) <- semi_dist
nmds_semi

data_nmds_semi <- as.data.frame(scores(nmds_semi))

data_nmds_semi <- rownames_to_column(data_nmds_semi, var = "Fragment_id")

data_nmds_semi$Species <- Ps$Species
data_nmds_semi$Treatment <- Ps$Treatment
data_nmds_semi$Organism_group <- Ps$Organism_group

#plot
data_nmds_semi$Organism_group<- factor(data_nmds_semi$Organism_group, levels = c("Stony coral", "Soft coral", 
                                                                                 "Macroalgae", "Sponge", "control"))

data_nmds_semi$Species<- factor(data_nmds_semi$Species, levels = c("Montipora digitata", "Pocillopora verrucosa", 
                                                                   "Xenia sp.", "Sinularia sp.",
                                                                   "Caulerpa sp.", "Peyssonnelia sp.",
                                                                   "Haliclona cnidata", "control"))

semi_species <- ggplot(data_nmds_semi, aes(x=sites.NMDS1, y=sites.NMDS2, col = Organism_group, shape = Species)) +
  geom_point(size = 3.5, stroke = 1.5) +
  scale_shape_manual(name = "Species", label=c("Montipora digitata" = expression(italic("Montipora digitata")),
                                               "Pocillopora verrucosa" = expression(italic("Pocillopora verrucosa")),
                                               "Xenia sp." = expression(paste(italic("Xenia")," sp.")),
                                               "Sinularia sp." = expression(paste(italic("Sinularia")," sp.")),
                                               "Caulerpa sp." = expression(paste(italic("Caulerpa")," sp.")),
                                               "Peyssonnelia sp." = expression(paste(italic("Peyssonnelia")," sp.")),
                                               "Haliclona cnidata" = expression(italic("Haliclona cnidata")),
                                               "control" = expression("control")),
                     values=c(3, 4, 8, 15, 16, 17, 18, 1)) +
  scale_colour_manual(values=c("#FF9900", "#FF3399","#009933", "#3399FF", "#999999")) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  annotate(geom = "text", x = 0.75, y = 0.2, size = 3,
           label = "Seminatural Polyculture
           2D Stress = 0.17",
           hjust = 1) +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(face = "bold", size = 8))

semi_species
rm(Ps, semi_matrix, nmds_semi, semi_dist, groups_species_ps)

#______________________
#create plot with the three treatments

arranged_treatment <- ggarrange(mono_species, poly_species, semi_species,
                                labels = c("a)", "b)", "c)"),
                                label.x = 0,
                                label.y = 0.95,
                                common.legend = TRUE,
                                legend = "bottom",
                                ncol = 3, nrow = 1,
                                hjust = -0.5,
                                vjust = 0)

arranged_treatment
rm(arranged_treatment, mono_species, poly_species, semi_species, data_nmds_mono, data_nmds_poly, data_nmds_semi)

#___________________________________________________________________
#within species
# Manuscript Fig. 4
# Caulerpa sp.
cau <- all %>% 
  subset(Species != "Peyssonnelia sp.") %>%
  subset(Species != "Sinularia sp.") %>%
  subset(Species != "Xenia sp.") %>%
  subset(Species != "Montipora digitata") %>%
  subset(Species != "Pocillopora verrucosa") %>%
  subset(Species != "Haliclona cnidata")

cau_data <- cau %>% as.data.frame()
rownames(cau_data) <- cau_data$Fragment_ID
cau_data <- cau_data %>% select(-c(Day, Fragment_ID, Species, Treatment, Organism_group))

cau_matrix <- as.matrix(cau_data)
rm(cau_data)
set.seed(12345)
cau_dist <- vegdist(cau_matrix, method="jaccard", binary = TRUE) 
cau_dist
set.seed(12345)
nmds_cau <- metaMDS(cau_dist)
sppscores(nmds_cau) <- cau_dist
nmds_cau

data_nmds_cau <- as.data.frame(scores(nmds_cau))
data_nmds_cau <- rownames_to_column(data_nmds_cau, var = "Fragment_id")
data_nmds_cau$Species <- cau$Species
data_nmds_cau$Treatment <- cau$Treatment

#plot
data_nmds_cau$Treatment<- factor(data_nmds_cau$Treatment, levels = c("Monoculture", "controlled_Polyculture", 
                                                                     "seminatural_Polyculture", "control"))

cau_treat <- ggplot(data_nmds_cau, aes(x=sites.NMDS1, y=sites.NMDS2, color = Treatment, shape = Treatment)) +
  geom_point(aes(colour=Treatment, fill=Treatment), size = 3, color = "black") +
  scale_fill_manual(values=c("#99FFFF", "#3333FF", "#000033", "#999999")) +
  scale_shape_manual(values=c(21, 21, 21, 4)) +
  theme_classic() +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  annotate(geom = "text", x = -0.7, y = 0.4, size = 3,
           label = "Caulerpa sp.
           2D Stress = 0.14",
           hjust = 0)

cau_treat
rm(cau, cau_matrix, cau_dist, groups_cau, nmds_cau)

#___________________________________________________________________
#Peyssonnelia sp.
pey <- all %>% 
  subset(Species != "Caulerpa sp.") %>%
  subset(Species != "Sinularia sp.") %>%
  subset(Species != "Xenia sp.") %>%
  subset(Species != "Montipora digitata") %>%
  subset(Species != "Pocillopora verrucosa") %>%
  subset(Species != "Haliclona cnidata")

pey <- pey %>% arrange(Treatment)
pey_data <- pey %>% as.data.frame()
rownames(pey_data) <- pey_data$Fragment_ID
pey_data <- pey_data %>% select(-c(Fragment_ID, Species, Treatment, Organism_group))
pey_matrix <- as.matrix(pey_data)
rm(pey_data)

set.seed(12345)
pey_dist <- vegdist(pey_matrix, method="jaccard", binary = TRUE) 
pey_dist
set.seed(12345)
nmds_pey <- metaMDS(pey_dist)
sppscores(nmds_pey) <- pey_dist
nmds_pey

data_nmds_pey <- as.data.frame(scores(nmds_pey))
data_nmds_pey <- rownames_to_column(data_nmds_pey, var = "Fragment_id")
data_nmds_pey$Species <- pey$Species
data_nmds_pey$Treatment <- pey$Treatment

#plot 
data_nmds_pey$Treatment<- factor(data_nmds_pey$Treatment, levels = c("Monoculture", "controlled_Polyculture", 
                                                                     "seminatural_Polyculture", "control"))

pey_treat <- ggplot(data_nmds_pey, aes(x=sites.NMDS1, y=sites.NMDS2, color = Treatment, shape = Treatment)) +
  geom_point(aes(colour=Treatment, fill=Treatment), size = 3, color = "black") +
  scale_fill_manual(values=c("#99FFFF", "#3333FF", "#000033", "#999999")) +
  scale_shape_manual(values=c(21, 21, 21, 4)) +
  theme_classic() +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  annotate(geom = "text", x = -0.8, y = 0.2, size = 3,
           label = "Peyssonnelia sp.
           2D Stress = 0.11",
           hjust = 0)

pey_treat
rm(pey, pey_matrix, pey_dist, groups_pey, nmds_pey)
#___________________________________________________________________
# Sinularia sp.
sin <- all %>% 
  subset(Species != "Peyssonnelia sp.") %>%
  subset(Species != "Caulerpa sp.") %>%
  subset(Species != "Xenia sp.") %>%
  subset(Species != "Montipora digitata") %>%
  subset(Species != "Pocillopora verrucosa") %>%
  subset(Species != "Haliclona cnidata")

sin <- sin %>% arrange(Treatment)
sin_data <- sin %>% as.data.frame()
rownames(sin_data) <- sin_data$Fragment_ID
sin_data <- sin_data %>% select(-c(Day, Fragment_ID, Species, Treatment, Organism_group))

sin_matrix <- as.matrix(sin_data)
rm(sin_data)
set.seed(12345)
sin_dist <- vegdist(sin_matrix, method="jaccard", binary = TRUE) 
sin_dist
set.seed(12345)
nmds_sin <- metaMDS(sin_dist)
sppscores(nmds_sin) <- sin_dist
nmds_sin

data_nmds_sin <- as.data.frame(scores(nmds_sin))
data_nmds_sin <- rownames_to_column(data_nmds_sin, var = "Fragment_id")
data_nmds_sin$Species <- sin$Species
data_nmds_sin$Treatment <- sin$Treatment

#plot
data_nmds_sin$Treatment<- factor(data_nmds_sin$Treatment, levels = c("Monoculture", "controlled_Polyculture", 
                                                                     "seminatural_Polyculture", "control"))

sin_treat <- ggplot(data_nmds_sin, aes(x=sites.NMDS1, y=sites.NMDS2, color = Treatment, shape = Treatment)) +
  geom_point(aes(colour=Treatment, fill=Treatment), size = 3, color = "black") +
  scale_fill_manual(values=c("#99FFFF", "#3333FF", "#000033", "#999999")) +
  scale_shape_manual(values=c(21, 21, 21, 4)) +
  theme_classic() +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  annotate(geom = "text", x = -0.7, y = 0.22, size = 3,
           label = "Sinularia sp.
           2D Stress = 0.10",
           hjust = 0)

sin_treat
rm(sin, sin_matrix, sin_dist, groups_sin, nmds_sin)

#___________________________________________________________________
# Xenia sp.
xen <- all %>% 
  subset(Species != "Peyssonnelia sp.") %>%
  subset(Species != "Sinularia sp.") %>%
  subset(Species != "Caulerpa sp.") %>%
  subset(Species != "Montipora digitata") %>%
  subset(Species != "Pocillopora verrucosa") %>%
  subset(Species != "Haliclona cnidata")

xen <- xen %>% arrange(Treatment)
xen_data <- xen %>% as.data.frame()
rownames(xen_data) <- xen_data$Fragment_ID
xen_data <- xen_data %>% select(-c(Day, Fragment_ID, Species, Treatment, Organism_group))
xen_matrix <- as.matrix(xen_data)
rm(xen_data)
set.seed(12345)
xen_dist <- vegdist(xen_matrix, method="jaccard", binary = TRUE) 
xen_dist
set.seed(12345)
nmds_xen <- metaMDS(xen_dist)
sppscores(nmds_xen) <- xen_dist
nmds_xen

data_nmds_xen <- as.data.frame(scores(nmds_xen))
data_nmds_xen <- rownames_to_column(data_nmds_xen, var = "Fragment_id")
data_nmds_xen$Species <- xen$Species
data_nmds_xen$Treatment <- xen$Treatment

#plot
data_nmds_xen$Treatment<- factor(data_nmds_xen$Treatment, levels = c("Monoculture", "controlled_Polyculture", 
                                                                     "seminatural_Polyculture", "control"))

xen_treat <- ggplot(data_nmds_xen, aes(x=sites.NMDS1, y=sites.NMDS2, color = Treatment, shape = Treatment)) +
  geom_point(aes(colour=Treatment, fill=Treatment), size = 3, color = "black") +
  scale_fill_manual(values=c("#99FFFF", "#3333FF", "#000033", "#999999")) +
  scale_shape_manual(values=c(21, 21, 21, 4)) +
  theme_classic() +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  annotate(geom = "text", x = -0.48, y = 0.35, size = 3,
           label = "Xenia sp.
           2D Stress = 0.12",
           hjust = 0)

xen_treat
rm(xen, xen_matrix, xen_dist, groups_xen, nmds_xen)

#___________________________________________________________________
# M. digitata
mdi <- all %>% 
  subset(Species != "Peyssonnelia sp.") %>%
  subset(Species != "Sinularia sp.") %>%
  subset(Species != "Xenia sp.") %>%
  subset(Species != "Caulerpa sp.") %>%
  subset(Species != "Pocillopora verrucosa") %>%
  subset(Species != "Haliclona cnidata")

mdi <- mdi %>% arrange(Treatment)
mdi_data <- mdi %>% as.data.frame()
rownames(mdi_data) <- mdi_data$Fragment_ID
mdi_data <- mdi_data %>% select(-c(Day, Fragment_ID, Species, Treatment, Organism_group))
mdi_matrix <- as.matrix(mdi_data)
rm(mdi_data)
set.seed(12345)
mdi_dist <- vegdist(mdi_matrix, method="jaccard", binary = TRUE) 
mdi_dist
set.seed(12345)
nmds_mdi <- metaMDS(mdi_dist)
sppscores(nmds_mdi) <- mdi_dist
nmds_mdi

data_nmds_mdi <- as.data.frame(scores(nmds_mdi))
data_nmds_mdi <- rownames_to_column(data_nmds_mdi, var = "Fragment_id")
data_nmds_mdi$Species <- mdi$Species
data_nmds_mdi$Treatment <- mdi$Treatment

#plot
data_nmds_mdi$Treatment<- factor(data_nmds_mdi$Treatment, levels = c("Monoculture", "controlled_Polyculture", 
                                                                     "seminatural_Polyculture", "control"))

mdi_treat <- ggplot(data_nmds_mdi, aes(x=sites.NMDS1, y=sites.NMDS2, color = Treatment, shape = Treatment)) +
  geom_point(aes(colour=Treatment, fill=Treatment), size = 3, color = "black") +
  scale_fill_manual(values=c("#99FFFF", "#3333FF", "#000033", "#999999")) +
  scale_shape_manual(values=c(21, 21, 21, 4)) +
  theme_classic() +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  annotate(geom = "text", x = -0.7, y = 0.15, size = 3,
           label = "M. digitata
           2D Stress = 0.12",
           hjust = 0)

mdi_treat
rm(mdi, mdi_matrix, mdi_dist, groups_mdi, nmds_mdi)

#___________________________________________________________________
# P. verrucosa
pve <- all %>% 
  subset(Species != "Peyssonnelia sp.") %>%
  subset(Species != "Sinularia sp.") %>%
  subset(Species != "Xenia sp.") %>%
  subset(Species != "Montipora digitata") %>%
  subset(Species != "Caulerpa sp.") %>%
  subset(Species != "Haliclona cnidata")

pve <- pve %>% arrange(Treatment)
pve_data <- pve %>% as.data.frame()
rownames(pve_data) <- pve_data$Fragment_ID
pve_data <- pve_data %>% select(-c(Day, Fragment_ID, Species, Treatment, Organism_group))
pve_matrix <- as.matrix(pve_data)
rm(pve_data)

set.seed(12345)
pve_dist <- vegdist(pve_matrix, method="jaccard", binary = TRUE) 
pve_dist
set.seed(12345)
nmds_pve <- metaMDS(pve_dist)
sppscores(nmds_pve) <- pve_dist
nmds_pve

data_nmds_pve <- as.data.frame(scores(nmds_pve))
data_nmds_pve <- rownames_to_column(data_nmds_pve, var = "Fragment_id") 
data_nmds_pve$Species <- pve$Species
data_nmds_pve$Treatment <- pve$Treatment

#plot
data_nmds_pve$Treatment<- factor(data_nmds_pve$Treatment, levels = c("Monoculture", "controlled_Polyculture", 
                                                                     "seminatural_Polyculture", "control"))

pve_treat <- ggplot(data_nmds_pve, aes(x=sites.NMDS1, y=sites.NMDS2, color = Treatment, shape = Treatment)) +
  geom_point(aes(colour=Treatment, fill=Treatment), size = 3, color = "black") +
  scale_fill_manual(values=c("#99FFFF", "#3333FF", "#000033", "#999999")) +
  scale_shape_manual(values=c(21, 21, 21, 4)) +
  theme_classic() +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  annotate(geom = "text", x = -0.48, y = 0.2, size = 3,
           label = "P. verrucosa
           2D Stress = 0.14",
           hjust = 0)

pve_treat
rm(pve, pve_matrix, pve_dist, groups_pve, nmds_pve)

#___________________________________________________________________
# H. cnidata
hcn <- all %>% 
  subset(Species != "Peyssonnelia sp.") %>%
  subset(Species != "Sinularia sp.") %>%
  subset(Species != "Xenia sp.") %>%
  subset(Species != "Montipora digitata") %>%
  subset(Species != "Pocillopora verrucosa") %>%
  subset(Species != "Caulerpa sp.")

hcn <- hcn %>% arrange(Treatment)
hcn_data <- hcn %>% as.data.frame()
rownames(hcn_data) <- hcn_data$Fragment_ID
hcn_data <- hcn_data %>% select(-c(Day, Fragment_ID, Species, Treatment, Organism_group))
hcn_matrix <- as.matrix(hcn_data)
rm(hcn_data)

set.seed(12345)
hcn_dist <- vegdist(hcn_matrix, method="jaccard", binary = TRUE) 
hcn_dist
set.seed(12345)
nmds_hcn <- metaMDS(hcn_dist)
sppscores(nmds_hcn) <- hcn_dist
nmds_hcn

data_nmds_hcn <- as.data.frame(scores(nmds_hcn))
data_nmds_hcn <- rownames_to_column(data_nmds_hcn, var = "Fragment_id")
data_nmds_hcn$Species <- hcn$Species
data_nmds_hcn$Treatment <- hcn$Treatment

#plot
data_nmds_hcn$Treatment<- factor(data_nmds_hcn$Treatment, levels = c("Monoculture", "controlled_Polyculture", 
                                                                     "seminatural_Polyculture", "control"))

hcn_treat <- ggplot(data_nmds_hcn, aes(x=sites.NMDS1, y=sites.NMDS2, color = Treatment, shape = Treatment)) +
  geom_point(aes(colour=Treatment, fill=Treatment), size = 3, color = "black") +
  scale_fill_manual(values=c("#99FFFF", "#3333FF", "#000033", "#999999")) +
  scale_shape_manual(values=c(21, 21, 21, 4)) +
  theme_classic() +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  annotate(geom = "text", x = -0.48, y = 0.28, size = 3,
           label = "H. cnidata
           2D Stress = 0.15",
           hjust = 0)

hcn_treat
rm(hcn, hcn_matrix, hcn_dist, groups_hcn, nmds_hcn)

#____________________________________________________
#merge
arranged_single_species_nodes <- ggarrange(cau_treat, pey_treat, sin_treat, xen_treat,
                                           mdi_treat, pve_treat, hcn_treat,
                                           labels = c("a)", "b)", "c)", "d)", "e)", "f)", "g)"),
                                           label.x = 0,
                                           label.y = 0.95,
                                           align = "hv",
                                           common.legend = TRUE,
                                           legend = "bottom",
                                           ncol = 2, nrow = 4,
                                           hjust = -0.5,
                                           vjust = 0) 

arranged_single_species_nodes
rm(list = ls())

#______________________________________________________________
#______________________________________________________________
# CFT (supplementary information)
# Supplementary information Fig. 28

all <- read_csv("data_prep/all_decontaminated.csv")
all <- all %>% select(-c(...1))

#change order for organism group and treatment
all <- all %>%
  arrange(Species)
all_org <- all %>%
  arrange(Organism_group)
all_treat <- all %>%
  arrange(Treatment)

all_nc_matrix <- select(all, -c(Species, Organism_group, Treatment, Day))
all_nc_matrix <- all_nc_matrix%>%as.data.frame
rownames(all_nc_matrix) <- all_nc_matrix$Fragment_ID
all_nc_matrix <- all_nc_matrix[, -1]
all_nc_matrix <- as.matrix(all_nc_matrix)

dist_all_nc <- vegdist(all_nc_matrix, method = "bray")
dist_all_nc

nmds_spec_nc <- metaMDS(dist_all_nc, k=2, trymax = 500)
nmds_spec_nc
sppscores(nmds_spec_nc) <- dist_all_nc

data_nmds_spec_nc <- as.data.frame(scores(nmds_spec_nc)) 
data_nmds_spec_nc
data_nmds_spec_nc <- rownames_to_column(data_nmds_spec_nc, var = "fragment_id")
data_nmds_spec_nc <- as_tibble(data_nmds_spec_nc)
data_nmds_spec_nc$Organism_group = all$Organism_group
data_nmds_spec_nc$Treatment = all$Treatment
data_nmds_spec_nc$Species = all$Species

# create plot
data_nmds_spec_nc$Species<- factor(data_nmds_spec_nc$Species, levels = c("Montipora digitata", "Pocillopora verrucosa",
                                                                         "Xenia sp.", "Sinularia sp.",
                                                                         "Caulerpa sp.","Peyssonnelia sp.",
                                                                         "Haliclona cnidata"))
data_nmds_spec_nc <- data_nmds_spec_nc %>% 
  relocate(Species, .before = fragment_id)

centroid_spec <- data_nmds_spec_nc %>%
  group_by(Species) %>%
  summarize(sites.NMDS1=mean(sites.NMDS1), sites.NMDS2=mean(sites.NMDS2))

Species <- ggplot(data_nmds_spec_nc, aes(x=sites.NMDS1, y=sites.NMDS2, color = Species)) +
  geom_point(aes(color=Species, fill=Species), size = 3, shape = 21, color = "black", 
             position=position_jitter(0.01)) +
  geom_point(data = centroid_spec, size = 3, shape = 4, stroke = 2, show_guide = FALSE) +
  scale_fill_manual(values=c("#FFCC00", "#FFFF33", "#FF3399", "#FF99CC", "#009933", "#00FF00", "#3399FF")) +
  scale_colour_manual(values=c( "#FFCC00", "#FFFF33","#FF3399", "#FF99CC",  "#009933", "#00FF00", "#3399FF")) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  stat_ellipse() +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 8, face = "italic"),
        legend.position = "bottom") +
  annotate(geom = "text", x = -0.8, y = 0.4, size = 2.5, 
           label = "Species 
           2D Stress = 0.22
           Cau ~ Hcn p = 0.006**
           Cau ~ Pey p = 0.014 *
           Xen ~ Hcn p = 0.016 *
           Xen ~ Pey p = 0.083 *",
           hjust = 0) +
  guides(fill=guide_legend(nrow=4, byrow=TRUE)) 

Species
#remove all species specific data
rm(all, dist_all_nc, all_nc_matrix, centroid_spec, nmds_spec_nc)

##____________________________________________
# treatment
treat_nc_matrix <- select(all_treat, -c(Species, Organism_group, Treatment, Day))
treat_nc_matrix <- treat_nc_matrix%>%as.data.frame
rownames(treat_nc_matrix) <- treat_nc_matrix$Fragment_ID
treat_nc_matrix <- treat_nc_matrix[, -1]
treat_nc_matrix <- as.matrix(treat_nc_matrix)
dist_treat_nc <- vegdist(treat_nc_matrix)
dist_treat_nc
nmds_treat_nc <- metaMDS(dist_treat_nc, k=2, trymax = 500)
nmds_treat_nc
sppscores(nmds_treat_nc) <- dist_treat_nc

data_nmds_treat_nc <- as.data.frame(scores(nmds_treat_nc)) 
data_nmds_treat_nc
data_nmds_treat_nc <- rownames_to_column(data_nmds_treat_nc, var = "fragment_id")
data_nmds_treat_nc <- as_tibble(data_nmds_treat_nc)
data_nmds_treat_nc$Organism_group = all_treat$Organism_group
data_nmds_treat_nc$Treatment = all_treat$Treatment
data_nmds_treat_nc$Species = all_treat$Species

data_nmds_treat_nc$Treatment<- factor(data_nmds_treat_nc$Treatment, levels = c("Monoculture",
                                                                               "controlled_Polyculture", 
                                                                               "seminatural_Polyculture"))
centroid_treat <- data_nmds_treat_nc%>%
  group_by(Treatment) %>%
  summarize(sites.NMDS1=mean(sites.NMDS1), sites.NMDS2=mean(sites.NMDS2))

treat <- ggplot(data_nmds_treat_nc, aes(x=sites.NMDS1, y=sites.NMDS2, color = Treatment)) +
  geom_point(aes(colour=Treatment, fill=Treatment), size = 3, shape = 21, color = "black",
             position=position_jitter(0.01)) +
  geom_point(data = centroid_treat, size = 3, shape = 4, stroke = 2, show_guide = FALSE) +
  stat_ellipse(aes(color = Treatment)) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_fill_manual(values=c("#99FFFF", "#3333FF", "#000033")) +
  scale_color_manual(values=c("#99FFFF", "#3333FF", "#000033")) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 8),
        legend.position = "bottom") +
  annotate(geom = "text", x = -0.6, y = 0.47, size = 2.5,
           label = "Treatment
           2D Stress = 0.22
           p > 0.05",
           hjust = 0) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE))

treat
#remove
rm(all_treat, treat_nc_matrix, dist_treat_nc, centroid_treat, nmds_treat_nc)

#____________________________________________
# organism group
org_nc_matrix <- select(all_org, -c(Species, Organism_group, Treatment, Day))
org_nc_matrix <- org_nc_matrix%>%as.data.frame
rownames(org_nc_matrix) <- org_nc_matrix$Fragment_ID
org_nc_matrix <- org_nc_matrix[, -1]
org_nc_matrix <- as.matrix(org_nc_matrix)
dist_org_nc <- vegdist(org_nc_matrix)
dist_org_nc

nmds_org_nc <- metaMDS(dist_org_nc, k=2, trymax = 500)
nmds_org_nc
sppscores(nmds_org_nc) <- dist_org_nc

data_nmds_org_nc <- as.data.frame(scores(nmds_org_nc)) 
data_nmds_org_nc
data_nmds_org_nc <- rownames_to_column(data_nmds_org_nc, var = "fragment_id")
data_nmds_org_nc <- as_tibble(data_nmds_org_nc)
data_nmds_org_nc$Organism_group = all_org$Organism_group
data_nmds_org_nc$Treatment = all_org$Treatment
data_nmds_org_nc$Species = all_org$Species

data_nmds_org_nc$Organism_group<- factor(data_nmds_org_nc$Organism_group, levels = c("Stony coral", "Soft coral", 
                                                                                     "Macroalgae", "Sponge"))

centroid_semi <- data_nmds_org_nc %>%
  group_by(Organism_group) %>%
  summarize(sites.NMDS1=mean(sites.NMDS1), sites.NMDS2=mean(sites.NMDS2))

Org_nc <- ggplot(data_nmds_org_nc, aes(x=sites.NMDS1, y=sites.NMDS2, colour = Organism_group)) +
  geom_point(aes(colour=Organism_group, fill=Organism_group), size = 3, shape = 21, color = "black",
             position=position_jitter(0.01)) +
  geom_point(data = centroid_semi, size = 3, shape = 4, stroke = 2, show_guide = FALSE) +
  scale_fill_manual(values=c("#FF9900", "#FF3399", "#009933", "#3399FF")) +
  scale_colour_manual(values=c("#FF9900", "#FF3399", "#009933", "#3399FF")) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  stat_ellipse() +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 8),
        legend.position = "bottom") +
  guides(color = guide_legend(nrow = 3)) +
  annotate(geom = "text", x = -0.6, y = 0.47, size = 2.5,
           label = "Organism group
           2D Stress = 0.22
           Sponge ~ Soft coral p = 0.042*",
           hjust = 0) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE))

Org_nc

#merge plots
arranged_overview <- ggarrange(treat, Org_nc, Species, 
                               labels = c("a)", "b)", "c)"),
                               align = "hv",
                               label.x = 0,
                               label.y = 0.95,
                               common.legend = FALSE,
                               # legend = c("bottom"),
                               ncol = 3, nrow = 1,
                               hjust = -0.5,
                               vjust = 0) 

arranged_overview
rm(all_org, org_nc_matrix, dist_org_nc, data_nmds_org_nc, nmds_org_nc, arranged_overview, centroid_semi, data_nmds_spec_nc,
   data_nmds_treat_nc, Org_nc, Species, treat)

#________________________________________________
# within treatment
# Supplementary information Fig. 29
all <- read_csv("data_prep/all_decontaminated.csv")
all <- all %>% select(-c(...1))

all <- all %>%
  arrange(Species)
mono <- all %>% 
  subset(Treatment =="Monoculture")
Pc <- all %>% 
  subset(Treatment =="controlled_Polyculture")
Ps <- all %>% 
  subset(Treatment =="seminatural_Polyculture")

mono <- mono %>% as.data.frame()
rownames(mono) <- mono$Fragment_ID
mono <- mono %>% select(-c(Fragment_ID, Species, Organism_group, Treatment, Day))
mono_matrix <- as.matrix(mono)
rm(mono)

set.seed(12345)
mono_dist <- vegdist(mono_matrix, method="bray")
mono_dist
set.seed(12345)
nmds_mono <- metaMDS(mono_dist)
sppscores(nmds_mono) <- mono_dist
nmds_mono

data_nmds_mono <- as.data.frame(scores(nmds_mono))
data_nmds_mono <- rownames_to_column(data_nmds_mono, var = "Fragment_id")
mono <- all %>% 
  subset(Treatment =="Monoculture")

data_nmds_mono$Species <- mono$Species
data_nmds_mono$Treatment <- mono$Treatment
data_nmds_mono$Organism_group <- mono$Organism_group

#plot
data_nmds_mono$Organism_group<- factor(data_nmds_mono$Organism_group, levels = c("Stony coral", "Soft coral", 
                                                                                 "Macroalgae", "Sponge"))

data_nmds_mono$Species<- factor(data_nmds_mono$Species, levels = c("Montipora digitata", "Pocillopora verrucosa", 
                                                                   "Xenia sp.", "Sinularia sp.",
                                                                   "Caulerpa sp.", "Peyssonnelia sp.",
                                                                   "Haliclona cnidata"))

mono_species <- ggplot(data_nmds_mono, aes(x=sites.NMDS1, y=sites.NMDS2, col = Organism_group, shape = Species)) +
  geom_point(size = 3.5, stroke = 1.5) +
  scale_shape_manual(name = "Species", label=c("Montipora digitata" = expression(italic("Montipora digitata")),
                                               "Pocillopora verrucosa" = expression(italic("Pocillopora verrucosa")),
                                               "Xenia sp." = expression(paste(italic("Xenia")," sp.")),
                                               "Sinularia sp." = expression(paste(italic("Sinularia")," sp.")),
                                               "Caulerpa sp." = expression(paste(italic("Caulerpa")," sp.")),
                                               "Peyssonnelia sp." = expression(paste(italic("Peyssonnelia")," sp.")),
                                               "Haliclona cnidata" = expression(italic("Haliclona cnidata"))),
                     values=c(3, 4, 8, 15, 16, 17, 18)) +
  scale_colour_manual(name = "Organism group", values=c("#FF9900", "#FF3399", "#009933", "#3399FF")) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  annotate(geom = "text", x = -0.4, y = 0.27, size = 3,
           label = "Monoculture
           2D Stress = 0.18",
           hjust = 0) +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(face = "bold", size = 8))

mono_species  
rm(mono, data_nmds_mono, mono_matrix, nmds_mono, mono_dist, groups_species_mo)

#________________________________________________
# polyculture controlled
poly <- Pc %>% as.data.frame()
rownames(poly) <- poly$Fragment_ID
poly <- poly %>% select(-c(Fragment_ID, Species, Organism_group, Treatment, Day))
poly_matrix <- as.matrix(poly)
rm(poly)
set.seed(12345)
poly_dist <- vegdist(poly_matrix, method="bray")
poly_dist
set.seed(12345)
nmds_poly <- metaMDS(poly_dist)
sppscores(nmds_poly) <- poly_dist
nmds_poly

data_nmds_poly <- as.data.frame(scores(nmds_poly))
data_nmds_poly <- rownames_to_column(data_nmds_poly, var = "Fragment_id")
data_nmds_poly$Species <- Pc$Species
data_nmds_poly$Treatment <- Pc$Treatment
data_nmds_poly$Organism_group <- Pc$Organism_group

#plot
data_nmds_poly$Organism_group<- factor(data_nmds_poly$Organism_group, levels = c("Stony coral", "Soft coral", 
                                                                                 "Macroalgae", "Sponge"))

data_nmds_poly$Species<- factor(data_nmds_poly$Species, levels = c("Montipora digitata", "Pocillopora verrucosa", 
                                                                   "Xenia sp.", "Sinularia sp.",
                                                                   "Caulerpa sp.", "Peyssonnelia sp.",
                                                                   "Haliclona cnidata"))

poly_species <- ggplot(data_nmds_poly, aes(x=sites.NMDS1, y=sites.NMDS2, col = Organism_group, shape = Species)) +
  geom_point(size = 3.5, stroke = 1.5) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_colour_manual(values=c("#FF9900", "#FF3399", "#009933", "#3399FF")) +
  scale_shape_manual(name = "Species", label=c("Montipora digitata" = expression(italic("Montipora digitata")),
                                               "Pocillopora verrucosa" = expression(italic("Pocillopora verrucosa")),
                                               "Xenia sp." = expression(paste(italic("Xenia")," sp.")),
                                               "Sinularia sp." = expression(paste(italic("Sinularia")," sp.")),
                                               "Caulerpa sp." = expression(paste(italic("Caulerpa")," sp.")),
                                               "Peyssonnelia sp." = expression(paste(italic("Peyssonnelia")," sp.")),
                                               "Haliclona cnidata" = expression(italic("Haliclona cnidata"))),
                     values=c(3, 4, 8, 15, 16, 17, 18)) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  annotate(geom = "text", x = -0.6, y = 0.38, size = 3,
           label = "controlled Polyculture
           2D Stress = 0.19",
           hjust = 0) +
  theme(legend.text = element_text( size = 8),
        legend.title = element_text(face = "bold", size = 8))

poly_species
rm(Pc, data_nmds_poly, poly_matrix, nmds_poly, poly_dist)

#______________________________________
# seminatural polyculture
semi <- Ps %>% as.data.frame()
rownames(semi) <- semi$Fragment_ID
semi <- semi %>% select(-c(Fragment_ID, Species, Organism_group, Treatment, Day))
semi_matrix <- as.matrix(semi)
rm(semi)
set.seed(12345)
semi_dist <- vegdist(semi_matrix, method="bray")
semi_dist
set.seed(12345)
nmds_semi <- metaMDS(semi_dist)
sppscores(nmds_semi) <- semi_dist
nmds_semi
data_nmds_semi <- as.data.frame(scores(nmds_semi))
data_nmds_semi <- rownames_to_column(data_nmds_semi, var = "Fragment_id")
data_nmds_semi$Species <- Ps$Species
data_nmds_semi$Treatment <- Ps$Treatment
data_nmds_semi$Organism_group <- Ps$Organism_group

#plot
data_nmds_semi$Organism_group<- factor(data_nmds_semi$Organism_group, levels = c("Stony coral", "Soft coral", 
                                                                                 "Macroalgae", "Sponge"))

data_nmds_semi$Species<- factor(data_nmds_semi$Species, levels = c("Montipora digitata", "Pocillopora verrucosa", 
                                                                   "Xenia sp.", "Sinularia sp.",
                                                                   "Caulerpa sp.", "Peyssonnelia sp.",
                                                                   "Haliclona cnidata"))

semi_species <- ggplot(data_nmds_semi, aes(x=sites.NMDS1, y=sites.NMDS2, col = Organism_group, shape = Species)) +
  geom_point(size = 3.5, stroke = 1.5) +
  scale_shape_manual(name = "Species", label=c("Montipora digitata" = expression(italic("Montipora digitata")),
                                               "Pocillopora verrucosa" = expression(italic("Pocillopora verrucosa")),
                                               "Xenia sp." = expression(paste(italic("Xenia")," sp.")),
                                               "Sinularia sp." = expression(paste(italic("Sinularia")," sp.")),
                                               "Caulerpa sp." = expression(paste(italic("Caulerpa")," sp.")),
                                               "Peyssonnelia sp." = expression(paste(italic("Peyssonnelia")," sp.")),
                                               "Haliclona cnidata" = expression(italic("Haliclona cnidata"))),
                     values=c(3, 4, 8, 15, 16, 17, 18)) +
  scale_colour_manual(values=c("#FF9900", "#FF3399","#009933", "#3399FF")) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  annotate(geom = "text", x = -0.4, y = 0.2, size = 3,
           label = "seminatural Polyculture
           2D Stress = 0.19",
           hjust = 0) +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(face = "bold", size = 8))

semi_species

# merge plots (Supplementary information Fig. 29)
arranged_treatment <- ggarrange(mono_species, poly_species, semi_species,
                                labels = c("a)", "b)", "c)"),
                                label.x = 0,
                                label.y = 0.95,
                                common.legend = TRUE,
                                legend = "bottom",
                                ncol = 3, nrow = 1,
                                hjust = -0.5,
                                vjust = 0) 

arranged_treatment
rm(mono_species, poly_species, semi_species, all, arranged_treatment, data_nmds_semi, nmds_semi, Ps, semi_matrix, semi_dist)

#_____________________________________________
# within species
# Supplementary information Fig. 30
all <- read_csv("data_prep/all_decontaminated.csv")
all <- all %>% select(-c(...1))

#___________________________________________________________________
#Caulerpa sp. 
cau <- all %>% subset(Species=="Caulerpa sp.")
cau_data <- cau %>% as.data.frame()
rownames(cau_data) <- cau_data$Fragment_ID
cau_data <- cau_data %>% select(-c(Fragment_ID, Species, Organism_group, Treatment, Day))
cau_matrix <- as.matrix(cau_data)
rm(cau_data)

set.seed(12345)
cau_dist <- vegdist(cau_matrix, method="bray")
cau_dist
set.seed(12345)
nmds_cau <- metaMDS(cau_dist)
sppscores(nmds_cau) <- cau_dist
nmds_cau

data_nmds_cau <- as.data.frame(scores(nmds_cau))
data_nmds_cau <- rownames_to_column(data_nmds_cau, var = "Fragment_id")
data_nmds_cau$Species <- cau$Species
data_nmds_cau$Treatment <- cau$Treatment
data_nmds_cau$Organism_group <- cau$Organism_group

#plot
cau_treat <- ggplot(data_nmds_cau, aes(x=sites.NMDS1, y=sites.NMDS2, color = Treatment, shape = Treatment)) +
  geom_point(aes(colour=Treatment, fill=Treatment), size = 3, color = "black") +
  scale_fill_manual(values=c("#99FFFF", "#3333FF", "#000033")) +
  scale_shape_manual(values=c(21, 21, 21)) +
  theme_classic() +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  annotate(geom = "text", x = -0.48, y = 0.16, size = 3,
           label = "Caulerpa sp.
           2D Stress = 0.04",
           hjust = 0)

cau_treat
#remove now unnecessary materials
rm(cau, cau_dist, cau_matrix, nmds_cau, data_nmds_cau)

#______________________________
#####Peyssonnelia
pey <- all %>% subset(Species=="Peyssonnelia sp.")
pey_data <- pey %>% as.data.frame()
rownames(pey_data) <- pey_data$Fragment_ID
pey_data <- pey_data %>% select(-c(Fragment_ID, Species, Organism_group, Treatment, Day))
pey_matrix <- as.matrix(pey_data)
rm(pey_data)
set.seed(12345)
pey_dist <- vegdist(pey_matrix, method="bray") 
pey_dist
set.seed(12345)
nmds_pey <- metaMDS(pey_dist)
sppscores(nmds_pey) <- pey_dist
nmds_pey

data_nmds_pey <- as.data.frame(scores(nmds_pey))
data_nmds_pey <- rownames_to_column(data_nmds_pey, var = "Fragment_id")
data_nmds_pey$Species <- pey$Species
data_nmds_pey$Treatment <- pey$Treatment
data_nmds_pey$Organism_group <- pey$Organism_group

#plot
pey_treat <- ggplot(data_nmds_pey, aes(x=sites.NMDS1, y=sites.NMDS2, colour = Treatment, shape = Treatment)) +
  geom_point(aes(colour=Treatment, fill=Treatment), size = 3, color = "black") +
  scale_fill_manual(values=c("#99FFFF", "#3333FF", "#000033")) +
  scale_shape_manual(values=c(21, 21, 21)) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_classic() +
  #stat_ellipse() +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  annotate(geom = "text", x = -0.3, y = 0.28, size = 3,
           label = "Peyssonnelia sp.,
           2D Stress = 0.14",
           hjust = 0)

pey_treat
#remove 
rm(pey, pey_dist, pey_matrix, nmds_pey, data_nmds_pey)

#______________________________________
# Xenia
xen <- all %>% subset(Species=="Xenia sp.")
xen_data <- xen %>% as.data.frame()
rownames(xen_data) <- xen_data$Fragment_ID
xen_data <- xen_data %>% select(-c(Fragment_ID, Species, Organism_group, Treatment, Day))
xen_matrix <- as.matrix(xen_data)
rm(xen_data)

set.seed(12345)
xen_dist <- vegdist(xen_matrix, method="bray")  
xen_dist

set.seed(12345)
nmds_xen <- metaMDS(xen_dist)
sppscores(nmds_xen) <- xen_dist
nmds_xen

data_nmds_xen <- as.data.frame(scores(nmds_xen))
data_nmds_xen <- rownames_to_column(data_nmds_xen, var = "Fragment_id")
data_nmds_xen$Species <- xen$Species
data_nmds_xen$Treatment <- xen$Treatment
data_nmds_xen$Organism_group <- xen$Organism_group

#plot
xen_treat <- ggplot(data_nmds_xen, aes(x=sites.NMDS1, y=sites.NMDS2, color = Treatment, shape = Treatment)) +
  geom_point(aes(colour=Treatment, fill=Treatment), size = 3, color = "black") +
  scale_fill_manual(values=c("#99FFFF", "#3333FF", "#000033")) +
  scale_shape_manual(values=c(21, 21, 21)) + 
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  annotate(geom = "text", x = -0.6, y = 0.2, size = 3,
           label = "Xenia sp.
           2D Stress = 0.08",
           hjust = 0)

xen_treat
#remove now unnecessary materials
rm(xen, xen_dist, xen_matrix, nmds_xen, data_nmds_xen)

#_________________________________________________________________________________
#Sinularia sp.
# and once without control (Treatments)
sin <- all %>% subset(Species=="Sinularia sp.") 
sin_data <- sin %>% as.data.frame()
rownames(sin_data) <- sin_data$Fragment_ID
sin_data <- sin_data %>% select(-c(Fragment_ID, Species, Organism_group, Treatment, Day))
sin_matrix <- as.matrix(sin_data)
rm(sin_data)

set.seed(12345)
sin_dist <- vegdist(sin_matrix, method="bray") 
sin_dist

set.seed(12345)
nmds_sin <- metaMDS(sin_dist)
sppscores(nmds_sin) <- sin_dist
nmds_sin

data_nmds_sin <- as.data.frame(scores(nmds_sin))
data_nmds_sin <- rownames_to_column(data_nmds_sin, var = "Fragment_id") 
data_nmds_sin$Species <- sin$Species
data_nmds_sin$Treatment <- sin$Treatment
data_nmds_sin$Organism_group <- sin$Organism_group

#plot
sin_treat <- ggplot(data_nmds_sin, aes(x=sites.NMDS1, y=sites.NMDS2, color = Treatment, shape = Treatment)) +
  geom_point(aes(colour=Treatment, fill=Treatment), size = 3, color = "black") +
  scale_fill_manual(values=c("#99FFFF", "#3333FF", "#000033")) +
  scale_shape_manual(values=c(21, 21, 21)) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  annotate(geom = "text", x = -0.4, y = 0.28, size = 3,
           label = "Sinularia sp.
           2D Stress = 0.09",
           hjust = 0)

sin_treat
#remove 
rm(sin, sin_dist, sin_matrix, nmds_sin, data_nmds_sin)


#_____________________________________________________________________________________________________
# M. digitata
mdi <- all %>% subset(Species=="Montipora digitata") 
mdi_data <- mdi %>% as.data.frame()
rownames(mdi_data) <- mdi_data$Fragment_ID
mdi_data <- mdi_data %>% select(-c(Fragment_ID, Species, Organism_group, Treatment, Day))
mdi_matrix <- as.matrix(mdi_data)
rm(mdi_data)
set.seed(12345)
mdi_dist <- vegdist(mdi_matrix, method="bray") 
mdi_dist

set.seed(12345)
nmds_mdi <- metaMDS(mdi_dist)
sppscores(nmds_mdi) <- mdi_dist
nmds_mdi

data_nmds_mdi <- as.data.frame(scores(nmds_mdi))
data_nmds_mdi <- rownames_to_column(data_nmds_mdi, var = "Fragment_id")
data_nmds_mdi$Species <- mdi$Species
data_nmds_mdi$Treatment <- mdi$Treatment
data_nmds_mdi$Organism_group <- mdi$Organism_group

#plot
mdi_treat <- ggplot(data_nmds_mdi, aes(x=sites.NMDS1, y=sites.NMDS2, color = Treatment, shape = Treatment)) +
  geom_point(aes(colour=Treatment, fill=Treatment), size = 3, color = "black") +
  scale_fill_manual(values=c("#99FFFF", "#3333FF", "#000033")) +
  scale_shape_manual(values=c(21, 21, 21)) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  annotate(geom = "text", x = -0.3, y = 0.3, size = 3,
           label = "M. digitata
           2D Stress = 0.09",
           hjust = 0)

mdi_treat
#remove
rm(mdi, mdi_dist, mdi_matrix, nmds_mdi, data_nmds_mdi)

#________________________________________ 
# P. verrucosa
pve <- all %>% subset(Species=="Pocillopora verrucosa")
pve_data <- pve %>% as.data.frame()
rownames(pve_data) <- pve_data$Fragment_ID
pve_data <- pve_data %>% select(-c(Fragment_ID, Species, Organism_group, Treatment, Day))
pve_matrix <- as.matrix(pve_data)
rm(pve_data)
set.seed(12345)
pve_dist <- vegdist(pve_matrix, method="bray") 
pve_dist
set.seed(12345)
nmds_pve <- metaMDS(pve_dist)
sppscores(nmds_pve) <- pve_dist
nmds_pve

data_nmds_pve <- as.data.frame(scores(nmds_pve))
data_nmds_pve <- rownames_to_column(data_nmds_pve, var = "Fragment_id")
data_nmds_pve$Species <- pve$Species
data_nmds_pve$Treatment <- pve$Treatment
data_nmds_pve$Organism_group <- pve$Organism_group

#plot 
pve_treat <- ggplot(data_nmds_pve, aes(x=sites.NMDS1, y=sites.NMDS2, color = Treatment, shape = Treatment)) +
  geom_point(aes(colour=Treatment, fill=Treatment), size = 3, color = "black") +
  scale_fill_manual(values=c("#99FFFF", "#3333FF", "#000033")) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_shape_manual(values=c(21, 21, 21)) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  annotate(geom = "text", x = -0.35, y = 0.3, size = 3,
           label = "P. verrucosa
           2D Stress = 0.08",
           hjust = 0)

pve_treat
#remove
rm(pve, pve_dist, pve_matrix, nmds_pve, data_nmds_pve)

#__________________________________________ 
# H. cnidata
hcn <- all %>% subset(Species=="Haliclona cnidata")
hcn_data <- hcn %>% as.data.frame()
rownames(hcn_data) <- hcn_data$Fragment_ID
hcn_data <- hcn_data %>% select(-c(Fragment_ID, Species, Organism_group, Treatment, Day))
hcn_matrix <- as.matrix(hcn_data)
rm(hcn_data)

set.seed(12345)
hcn_dist <- vegdist(hcn_matrix, method="bray") 
hcn_dist
set.seed(12345)
nmds_hcn <- metaMDS(hcn_dist)
sppscores(nmds_hcn) <- hcn_dist
nmds_hcn

data_nmds_hcn <- as.data.frame(scores(nmds_hcn))
data_nmds_hcn <- rownames_to_column(data_nmds_hcn, var = "Fragment_id")
data_nmds_hcn$Species <- hcn$Species
data_nmds_hcn$Treatment <- hcn$Treatment
data_nmds_hcn$Organism_group <- hcn$Organism_group

#plot
hcn_treat <- ggplot(data_nmds_hcn, aes(x=sites.NMDS1, y=sites.NMDS2, color = Treatment, shape = Treatment)) +
  geom_point(aes(colour=Treatment, fill=Treatment), size = 3, color = "black") +
  scale_fill_manual(values=c("#99FFFF", "#3333FF", "#000033")) +
  scale_shape_manual(values=c(21, 21, 21)) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  annotate(geom = "text", x = -0.4, y = 0.34, size = 3,
           label = "H. cnidata
           2D Stress = 0.14",
           hjust = 0)

hcn_treat
#remove
rm(hcn, hcn_dist, hcn_matrix, nmds_hcn, data_nmds_hcn)

#arrange all 7 plots to one overview plot
arranged_single_species <- ggarrange(cau_treat, pey_treat, sin_treat, xen_treat,
                                     mdi_treat, pve_treat, hcn_treat,
                                     labels = c("a)", "b)", "c)", "d)", "e)", "f)", "g)"),
                                     font.label = list(size = 12, color = "black", face = "italic", family = NULL),
                                     label.x = 0,
                                     label.y = 0.04,
                                     common.legend = TRUE,
                                     legend = "bottom",
                                     ncol = 2, nrow = 4,
                                     hjust = -0.5,
                                     vjust = -1) +
  theme(legend.text = element_text(size = 22))

arranged_single_species
rm(all, arranged_single_species, cau_treat, hcn_treat, mdi_treat, pey_treat, pve_treat, sin_treat, xen_treat)
