# Script created by I'ah Donovan-Banfield to parse DiversiTools amino acid outputs for 
# data visualisation for Figs. 2, 3 and Supplementary Figs. 2-4 in the manuscript 
# "Characterisation of SARS-CoV-2 genomic variations in response to molnupiravir treatment in the AGILE Phase IIa clinical trial."

library(tidyverse)
library(ggplot2)
library(cowplot)
library(scales)


# set working directory to directory with your diversitools outputs from samples that passed the genome criteria
setwd("~/data/agile-nimagen/easyseq_bams/easyseq_forR/easyseq_DT_outputs")

##### read in for DiversiTools outputs function - amino acids "sample-name_AA.txt" ####
diversitools_readin_aa <- function(filepath ="./", file_pattern, file_list) {
  
  temp <- list.files(filepath, pattern=file_pattern) # creates list with filenames as in the directory
  
  file_list <- list() # creates empty list to add to in the for loop
  
  for (i in 1:length(temp)) { # read in all data files as DFs, add Proportion(nonsyn) column to all 
    
    file_list[[i]] <- read.delim(temp[i]) 
    file_list[[i]]$AAcoverage <- as.numeric(file_list[[i]]$AAcoverage)
    file_list[[i]]$TopAAcnt <- as.numeric(file_list[[i]]$TopAAcnt) # changes class to numeric
    file_list[[i]]$SndAAcnt <- as.numeric(file_list[[i]]$SndAAcnt)
    file_list[[i]]$TrdAAcnt <- as.numeric(file_list[[i]]$TrdAAcnt)
    file_list[[i]]$ProportionTopAA <-  file_list[[i]]$TopAAcnt / file_list[[i]]$AAcoverage
    file_list[[i]]$ProportionSndAA <-  file_list[[i]]$SndAAcnt / file_list[[i]]$AAcoverage
    file_list[[i]]$ProportionTrdAA <-  file_list[[i]]$TrdAAcnt / file_list[[i]]$AAcoverage
    file_list[[i]]$AAPosition <- as.numeric(file_list[[i]]$AAPosition)
    file_list[[i]][,7:14] <- lapply(file_list[[i]][,7:14], as.character)
    file_list[[i]]$TopCodoncnt <- as.numeric(file_list[[i]]$TopCodoncnt) # changes class to numeric
    file_list[[i]]$SndCodoncnt <- as.numeric(file_list[[i]]$SndCodoncnt)
    file_list[[i]]$TrdCodoncnt <- as.numeric(file_list[[i]]$TrdCodoncnt)
    file_list[[i]]$AAPsn <- 1:nrow(file_list[[i]]) #adding another column to set the positions as the number of rows so that they are not filtered out later on
    file_list[[i]]$Protein <- factor(file_list[[i]]$Protein, levels=c("nsp1", "nsp2", "nsp3", "nsp4", "nsp5", "nsp6", "nsp7", "nsp8", "nsp9",
                                                                      "nsp10", "nsp12", "nsp12_2", "nsp13", "nsp14", "nsp15", "nsp16", "S", 
                                                                      "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10"))
    #adding levels to the Protein column so that they come up in the right order along the output graph and in the key
    file_list[[i]] <- file_list[[i]]  %>% filter(Sample != "Sample")
    }
  
  myfilenames <- gsub("_AA.txt","", temp) # give names to each file, removing unnecessary strings
  names(file_list) <- myfilenames
  return(file_list)
}

####AGILE data ####
# create list of your dataframes 
aa_list <- diversitools_readin_aa(filepath = "./", file_pattern = "_AA.txt", file_list = aa_list)
# convert to large df 
all_aa <- bind_rows(aa_list)
# change column label to match metadata spreadsheet
all_aa <- rename(all_aa, kit_id=Sample)
# if your filename contained the barcodes - this removes them 
all_aa$kit_id <- gsub("_[^_]*", "", all_aa$kit_id)
# read in metadata
metadata <- read.csv("~/projects/agile_mpv_seq_project/sample_info/AGILE-nimagen-full-180_final-metadata.csv")
# merge sample dataframe with metadata 
metadata <- read.csv("~/projects/agile_mpv_seq_project/sample_info/AGILE-nimagen-full-180_final-metadata.csv")
# add vax status 
vax <- read.csv("~/projects/agile_mpv_seq_project/sample_info/AGILE_CST2_group_info.csv")
vax <- rename(vax, tx_group=Group)
metadata <- merge(metadata, vax, by=c("Subject", "tx_group"))
# merge sample dataframe with metadata 
all_aa_meta <- merge(all_aa, metadata, by="kit_id")
# select only the useful columns
all_aa_mv <- select(all_aa_meta, 3:6, 14, 16, 18, 26:33, 35:37, 43:44)
# create one column for AA ranks
all_aa_mv_m <- pivot_longer(all_aa_mv, cols=c("TopAA", "SndAA", "TrdAA"), names_to = "AA_rank", values_to = "AA")
# creat column for the proportion values
all_aa_mv_m <- pivot_longer(all_aa_mv_m, cols=c("ProportionTopAA", "ProportionSndAA", "ProportionTrdAA"), names_to = "Prop_AA_name", values_to = "Prop_AA")
# these steps filter out the duplicate rows so you only have three rows per AA position (TOp, 2nd and 3rd)
all_aa_mv_m$Prop_AA_name <- gsub("Proportion", "", all_aa_mv_m$Prop_AA_name)
all_aa_mv_m_final <- all_aa_mv_m %>% filter(AA_rank == Prop_AA_name)
all_aa_mv_m_final$Subject <- factor(all_aa_mv_m_final$Subject)
all_aa_mv_m_final$Vaccinated <- factor(all_aa_mv_m_final$Vaccinated, levels= c("No", "Yes"))
# make tx_group and AA_rank columns factors for plotting
all_aa_mv_m_final$tx_group <- factor(all_aa_mv_m_final$tx_group, levels=c("Group B", "Group A"))
# rename Snd and Trd and set class to factor
all_aa_mv_m_final$AA_rank <- gsub("SndAA", "2ndAA", all_aa_mv_m_final$AA_rank)
all_aa_mv_m_final$AA_rank <- gsub("TrdAA", "3rdAA", all_aa_mv_m_final$AA_rank)
all_aa_mv_m_final$AA_rank <- factor(all_aa_mv_m_final$AA_rank, levels=c("TopAA", "2ndAA", "3rdAA"))

### Plotting ####
# all plots filter out positions with <200 Coverage

# labels 
group.labs <- c( "placebo", "molnupiravir") 
names(group.labs) <- c("Group B", "Group A")

day.labs <- c("Day 1", "Day 3", "Day 5")
names(day.labs) <- c("1", "3", "5")

## subset melted dataframe into lineage groups for speed of plotting ####
alpha_final <- all_aa_mv_m_final %>% filter(AA!="") %>% filter(AA!="<NA>") %>% 
  filter(lineage_group == "alpha") %>% filter(AAcoverage>200)
eu1_final <- all_aa_mv_m_final %>% filter(AA!="") %>% filter(AA!="<NA>") %>% 
  filter(lineage_group == "B.1.177/EU1") %>% filter(AAcoverage>200)
delta_final <- all_aa_mv_m_final %>% filter(AA!="") %>% filter(AA!="<NA>") %>% 
  filter(lineage_group == "delta") %>% filter(AAcoverage>200)
ba1_final <- all_aa_mv_m_final %>% filter(AA!="") %>% filter(AA!="<NA>") %>% 
  filter(lineage_group == "BA.1") %>% filter(AAcoverage>200)


# plot whole genome minor variations ####
plot_wg_mv_alpha <- ggplot(alpha_final, aes(x=RefSite, y=Prop_AA, colour=AA_rank)) + 
  geom_point(aes(size=AAcoverage), alpha=0.5) +
  facet_grid(tx_group~visit_day,
             labeller = labeller(tx_group = group.labs, visit_day=day.labs)) +
  ylim(-0.05,1.05) +
  ylab("Proportion of amino acid") + xlab("Codon position (1000s)") +
  ggtitle("B.1.1.7 / Alpha: whole genome \n") +
  theme_bw() +
  scale_colour_viridis_d(begin = 0, end = 0.8, direction=-1) +
  scale_x_continuous(breaks = seq(0,30000,5000), 
                     labels = label_number(scale = 1e-3)) + 
  labs(colour="AA rank", size="Read depth") +
  guides(colour=guide_legend(order=1), size = guide_legend(order=2)) +
  geom_hline(yintercept=c(0.10, 0.5), linetype= "dashed", colour="grey") +
  scale_size_area(breaks=c(1000, 7500, 15000, 22500, 30000)) +
  theme(axis.text = element_text(colour="black", size=10),
        axis.title.y = element_text(vjust=2),
        axis.title = element_text(size=14),
        plot.title = element_text(face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(hjust=0.5, size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.box="vertical", legend.margin=margin(), legend.key.size =unit(1, "cm"),
        legend.spacing.y = unit(-0.15, "cm"), legend.key = element_rect(size = 5, color = 'NA'),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), units = "pt"))

plot_wg_mv_EU1 <- ggplot(eu1_final, aes(x=RefSite, y=Prop_AA, colour=AA_rank)) + 
  geom_point(aes(size=AAcoverage), alpha=0.5) +
  facet_grid(tx_group~visit_day,
             labeller = labeller(tx_group = group.labs, visit_day=day.labs)) +
  ylim(-0.05,1.05) + 
  ylab("Proportion of amino acid") + xlab("Codon position (1000s)") +
  ggtitle("B.1.177 / EU1: whole genome \n") +
  theme_bw() +
  scale_colour_viridis_d(begin = 0, end = 0.8, direction=-1) +
  scale_x_continuous(breaks = seq(0,30000,5000), 
                     labels = label_number(scale = 1e-3)) + 
  labs(colour="AA rank", size="Read depth") +
  guides(colour=guide_legend(order=1), size = guide_legend(order=2)) +
  geom_hline(yintercept=c(0.10, 0.5), linetype= "dashed", colour="grey") +
  scale_size_area(breaks=c(1000, 7500, 15000, 22500, 30000)) +
  theme(axis.text = element_text(colour="black", size=10),
        axis.title.y = element_text(vjust=2),
        axis.title = element_text(size=14),
        plot.title = element_text(face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(hjust=0.5, size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.box="vertical", legend.margin=margin(), legend.key.size =unit(1, "cm"),
        legend.spacing.y = unit(-0.15, "cm"), legend.key = element_rect(size = 5, color = 'NA'),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), units = "pt"))

plot_wg_mv_delta <- ggplot(delta_final, aes(x=RefSite, y=Prop_AA, colour=AA_rank)) + 
  geom_point(aes(size=AAcoverage), alpha=0.5) +
  facet_grid(tx_group~visit_day,
             labeller = labeller(tx_group = group.labs, visit_day=day.labs)) +
  ylim(-0.05,1.05) + 
  ylab("Proportion of amino acid") + xlab("Codon position (1000s)") +
  ggtitle("B.1.617.2 / Delta: whole genome \n") +
  theme_bw() +
  scale_colour_viridis_d(begin = 0, end = 0.8, direction=-1) +
  scale_x_continuous(breaks = seq(0,30000,5000), 
                     labels = label_number(scale = 1e-3)) + 
  labs(colour="AA rank", size="Read depth") +
  guides(colour=guide_legend(order=1), size = guide_legend(order=2)) +
  geom_hline(yintercept=c(0.10, 0.5), linetype= "dashed", colour="grey") +
  scale_size_area(breaks=c(1000, 7500, 15000, 22500, 30000)) +
  theme(axis.text = element_text(colour="black", size=10),
        axis.title.y = element_text(vjust=2),
        axis.title = element_text(size=14),
        plot.title = element_text(face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(hjust=0.5, size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.box="vertical", legend.margin=margin(), legend.key.size =unit(1, "cm"),
        legend.spacing.y = unit(-0.15, "cm"), legend.key = element_rect(size = 5, color = 'NA'),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), units = "pt"))

plot_wg_mv_BA1 <- ggplot(ba1_final, aes(x=RefSite, y=Prop_AA, colour=AA_rank)) + 
  geom_point(aes(size=AAcoverage), alpha=0.5) +
  facet_grid(tx_group~visit_day,
             labeller = labeller(tx_group = group.labs, visit_day=day.labs)) +
  ylim(-0.05,1.05) + 
  ylab("Proportion of amino acid") + xlab("Codon position (1000s)") +
  ggtitle("BA.1 / Omicron: whole genome \n") +
  theme_bw() +
  scale_colour_viridis_d(begin = 0, end = 0.8, direction=-1) +
  scale_x_continuous(breaks = seq(0,30000,5000), 
                     labels = label_number(scale = 1e-3)) + 
  labs(colour="AA rank", size="Read depth") +
  guides(colour=guide_legend(order=1), size = guide_legend(order=2)) +
  geom_hline(yintercept=c(0.10, 0.5), linetype= "dashed", colour="grey") +
  scale_size_area(breaks=c(1000, 7500, 15000, 22500, 30000)) +
  theme(axis.text = element_text(colour="black", size=10),
        axis.title.y = element_text(vjust=2),
        axis.title = element_text(size=14),
        plot.title = element_text(face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(hjust=0.5, size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.box="vertical", legend.margin=margin(), legend.key.size =unit(1, "cm"),
        legend.spacing.y = unit(-0.15, "cm"), legend.key = element_rect(size = 5, color = 'NA'),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), units = "pt"))

#write out 
ggsave(plot_wg_mv_alpha, filename = "~/projects/agile_mpv_seq_project/Reports/whole-genome_mv_alpha.png", device="png", dpi=300, width=12, height=6)
ggsave(plot_wg_mv_EU1, filename = "~/projects/agile_mpv_seq_project/Reports/whole-genome_mv_EU1.png", device="png", dpi=300, width=12, height=6)
ggsave(plot_wg_mv_delta, filename = "~/projects/agile_mpv_seq_project/Reports/whole-genome_mv_delta.png", device="png", dpi=300, width=12, height=6)
ggsave(plot_wg_mv_BA1, filename = "~/projects/agile_mpv_seq_project/Reports/whole-genome_mv_BA1.png", device="png", dpi=300, width=12, height=6)

# NSP14 data #### 
#pulling out nsp14 data + 200 cov cut off
nsp14_sum <- all_aa_mv_m_final %>% filter(Protein=="nsp14") %>% filter(AAcoverage>200)

#plots for all variants' nsp14 
plot_nsp14_mv_alpha <- ggplot(nsp14_sum %>% filter(AA!="") %>% filter(AA!="<NA>") %>% filter(lineage_group == "alpha"), 
                              aes(x=RefSite, y=Prop_AA, colour=AA_rank)) + 
  geom_point(aes(size=AAcoverage), alpha=0.5) +
  facet_grid(tx_group~visit_day,
             labeller = labeller(tx_group = group.labs, visit_day=day.labs)) +
  ylim(-0.05,1.05) + 
  ylab("Proportion of amino acid") + xlab("Codon position (1000s)") +
  ggtitle("Nsp14") +
  theme_bw() +
  scale_colour_viridis_d(begin = 0, end = 0.8, direction=-1) +
  scale_x_continuous(labels = label_number(scale = 1e-3)) + 
  labs(colour="AA rank", size="Read depth") +
  guides(colour=guide_legend(order=1), size = guide_legend(order=2)) +
  geom_hline(yintercept=c(0.10, 0.5), linetype= "dashed", colour="grey") +
  scale_size_area(breaks=c(1000, 7500, 15000, 22500, 30000)) +
  theme(axis.text = element_text(colour="black", size=10),
        axis.title.y = element_text(vjust=2),
        axis.title = element_text(size=14),
        plot.title = element_text(face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(hjust=0.5, size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.box="vertical", legend.margin=margin(), legend.key.size =unit(1, "cm"),
        legend.spacing.y = unit(-0.15, "cm"), legend.key = element_rect(size = 5, color = 'NA'),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), units = "pt"))

plot_nsp14_mv_delta <- ggplot(nsp14_sum %>% filter(AA!="") %>% filter(AA!="<NA>") %>% filter(lineage_group == "delta"),
                              aes(x=RefSite, y=Prop_AA, colour=AA_rank)) + 
  geom_point(aes(size=AAcoverage), alpha=0.5) +
  facet_grid(tx_group~visit_day,
             labeller = labeller(tx_group = group.labs, visit_day=day.labs)) +
  ylim(-0.05,1.05) + 
  ylab("Proportion of amino acid") + xlab("Codon position (1000s)") +
  ggtitle("Nsp14") +
  theme_bw() +
  scale_colour_viridis_d(begin = 0, end = 0.8, direction=-1) +
  scale_x_continuous(labels = label_number(scale = 1e-3)) + 
  labs(colour="AA rank", size="Read depth") +
  guides(colour=guide_legend(order=1), size = guide_legend(order=2)) +
  geom_hline(yintercept=c(0.10, 0.5), linetype= "dashed", colour="grey") +
  scale_size_area(breaks=c(1000, 7500, 15000, 22500, 30000)) +
  theme(axis.text = element_text(colour="black", size=10),
        axis.title.y = element_text(vjust=2),
        axis.title = element_text(size=14),
        plot.title = element_text(face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(hjust=0.5, size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.box="vertical", legend.margin=margin(), legend.key.size =unit(1, "cm"),
        legend.spacing.y = unit(-0.15, "cm"), legend.key = element_rect(size = 5, color = 'NA'),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), units = "pt"))

plot_nsp14_mv_EU1 <- ggplot(nsp14_sum %>% filter(AA!="") %>% filter(AA!="<NA>") %>% filter(lineage_group == "B.1.177/EU1"),
                            aes(x=RefSite, y=Prop_AA, colour=AA_rank)) + 
  geom_point(aes(size=AAcoverage), alpha=0.5) +
  facet_grid(tx_group~visit_day,
             labeller = labeller(tx_group = group.labs, visit_day=day.labs)) +
  ylim(-0.05,1.05) + 
  ylab("Proportion of amino acid") + xlab("Codon position (1000s)") +
  ggtitle("Nsp14") +
  theme_bw() +
  scale_colour_viridis_d(begin = 0, end = 0.8, direction=-1) +
  scale_x_continuous(labels = label_number(scale = 1e-3)) + 
  labs(colour="AA rank", size="Read depth") +
  guides(colour=guide_legend(order=1), size = guide_legend(order=2)) +
  geom_hline(yintercept=c(0.10, 0.5), linetype= "dashed", colour="grey") +
  scale_size_area(breaks=c(1000, 7500, 15000, 22500, 30000)) +
  theme(axis.text = element_text(colour="black", size=10),
        axis.title.y = element_text(vjust=2),
        axis.title = element_text(size=14),
        plot.title = element_text(face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(hjust=0.5, size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.box="vertical", legend.margin=margin(), legend.key.size =unit(1, "cm"),
        legend.spacing.y = unit(-0.15, "cm"), legend.key = element_rect(size = 5, color = 'NA'),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), units = "pt"))

plot_nsp14_mv_BA1 <- ggplot(nsp14_sum %>% filter(AA!="") %>% filter(AA!="<NA>") %>% filter(lineage_group == "BA.1"), 
                            aes(x=RefSite, y=Prop_AA, colour=AA_rank)) + 
  geom_point(aes(size=AAcoverage), alpha=0.5) +
  facet_grid(tx_group~visit_day,
             labeller = labeller(tx_group = group.labs, visit_day=day.labs)) +
  ylim(-0.05,1.05) + 
  ylab("Proportion of amino acid") + xlab("Codon position (1000s)") +
  ggtitle("Nsp14") +
  theme_bw() +
  scale_colour_viridis_d(begin = 0, end = 0.8, direction=-1) +
  scale_x_continuous(labels = label_number(scale = 1e-3)) + 
  labs(colour="AA rank", size="Read depth") +
  guides(colour=guide_legend(order=1), size = guide_legend(order=2)) +
  geom_hline(yintercept=c(0.10, 0.5), linetype= "dashed", colour="grey") +
  scale_size_area(breaks=c(1000, 7500, 15000, 22500, 30000)) +
  theme(axis.text = element_text(colour="black", size=10),
        axis.title.y = element_text(vjust=2),
        axis.title = element_text(size=14),
        plot.title = element_text(face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(hjust=0.5, size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.box="vertical", legend.margin=margin(), legend.key.size =unit(1, "cm"),
        legend.spacing.y = unit(-0.15, "cm"), legend.key = element_rect(size = 5, color = 'NA'),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), units = "pt"))

## write out NSP14 plots 
# for git 
ggsave(plot_nsp14_mv_alpha, filename = "~/projects/agile_mpv_seq_project/Reports/nsp14_mv_alpha_agile.png", device="png", dpi=300, width=8, height=12)
ggsave(plot_nsp14_mv_delta, filename = "~/projects/agile_mpv_seq_project/Reports/nsp14_mv_delta_agile.png", device="png", dpi=300, width=8, height=12)
ggsave(plot_nsp14_mv_EU1, filename = "~/projects/agile_mpv_seq_project/Reports/nsp14_mv_EU1_agile.png", device="png", dpi=300, width=8, height=12)
ggsave(plot_nsp14_mv_BA1, filename = "~/projects/agile_mpv_seq_project/Reports/nsp14_mv_BA1_agile.png", device="png", dpi=300, width=8, height=12)

# nsp12 data ####
# make nsp12 one thing 
all_aa_mv_m_final[all_aa_mv_m_final$Protein %in% c("nsp12", "nsp12_2"), ]$Protein <- "nsp12"

nsp12_sum <- all_aa_mv_m_final %>% filter(Protein=="nsp12") %>% filter(AAcoverage>200)

#plot for each lineage 
plot_nsp12_mv_alpha <- ggplot(nsp12_sum %>% filter(AA!="") %>% filter(AA!="<NA>") %>% filter(lineage_group == "alpha"),
                              aes(x=RefSite, y=Prop_AA, colour=AA_rank)) + 
  geom_point(aes(size=AAcoverage), alpha=0.5) +
  facet_grid(tx_group~visit_day,
             labeller = labeller(tx_group = group.labs, visit_day=day.labs)) +
  ylim(-0.05,1.05) + 
  ylab("Proportion of amino acid") + xlab("Codon position (1000s)") +
  ggtitle("Nsp12") +
  theme_bw() +
  scale_colour_viridis_d(begin = 0, end = 0.8, direction=-1) +
  scale_x_continuous(labels = label_number(scale = 1e-3)) + 
  labs(colour="AA rank", size="Read depth") +
  guides(colour=guide_legend(order=1), size = guide_legend(order=2)) +
  geom_hline(yintercept=c(0.10, 0.5), linetype= "dashed", colour="grey") +
  scale_size_area(breaks=c(1000, 7500, 15000, 22500, 30000)) +
  theme(axis.text = element_text(colour="black", size=10),
        axis.title.y = element_text(vjust=2),
        axis.title = element_text(size=14),
        plot.title = element_text(face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(hjust=0.5, size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.box="vertical", legend.margin=margin(), legend.key.size =unit(1, "cm"),
        legend.spacing.y = unit(-0.15, "cm"), legend.key = element_rect(size = 5, color = 'NA'),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), units = "pt"))

plot_nsp12_mv_delta <- ggplot(nsp12_sum %>% filter(AA!="") %>% filter(AA!="<NA>") %>% filter(lineage_group == "delta"), 
                              aes(x=RefSite, y=Prop_AA, colour=AA_rank)) + 
  geom_point(aes(size=AAcoverage), alpha=0.5) +
  facet_grid(tx_group~visit_day,
             labeller = labeller(tx_group = group.labs, visit_day=day.labs)) +
  ylim(-0.05,1.05) + 
  ylab("Proportion of amino acid") + xlab("Codon position (1000s)") +
  ggtitle("Nsp12") +
  theme_bw() +
  scale_colour_viridis_d(begin = 0, end = 0.8, direction=-1) +
  scale_x_continuous(labels = label_number(scale = 1e-3)) + 
  labs(colour="AA rank", size="Read depth") +
  guides(colour=guide_legend(order=1), size = guide_legend(order=2)) +
  geom_hline(yintercept=c(0.10, 0.5), linetype= "dashed", colour="grey") +
  scale_size_area(breaks=c(1000, 7500, 15000, 22500, 30000)) +
  theme(axis.text = element_text(colour="black", size=10),
        axis.title.y = element_text(vjust=2),
        axis.title = element_text(size=14),
        plot.title = element_text(face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(hjust=0.5, size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.box="vertical", legend.margin=margin(), legend.key.size =unit(1, "cm"),
        legend.spacing.y = unit(-0.15, "cm"), legend.key = element_rect(size = 5, color = 'NA'),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), units = "pt"))

plot_nsp12_mv_EU1 <- ggplot(nsp12_sum %>% filter(AA!="") %>% filter(AA!="<NA>") %>% filter(lineage_group == "B.1.177/EU1"),
                            aes(x=RefSite, y=Prop_AA, colour=AA_rank)) + 
  geom_point(aes(size=AAcoverage), alpha=0.5) +
  facet_grid(tx_group~visit_day,
             labeller = labeller(tx_group = group.labs, visit_day=day.labs)) +
  ylim(-0.05,1.05) + 
  ylab("Proportion of amino acid") + xlab("Codon position (1000s)") +
  ggtitle("Nsp12") +
  theme_bw() +
  scale_colour_viridis_d(begin = 0, end = 0.8, direction=-1) +
  scale_x_continuous(labels = label_number(scale = 1e-3)) + 
  labs(colour="AA rank", size="Read depth") +
  guides(colour=guide_legend(order=1), size = guide_legend(order=2)) +
  geom_hline(yintercept=c(0.10, 0.5), linetype= "dashed", colour="grey") +
  scale_size_area(breaks=c(1000, 7500, 15000, 22500, 30000)) +
  theme(axis.text = element_text(colour="black", size=10),
        axis.title.y = element_text(vjust=2),
        axis.title = element_text(size=14),
        plot.title = element_text(face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(hjust=0.5, size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.box="vertical", legend.margin=margin(), legend.key.size =unit(1, "cm"),
        legend.spacing.y = unit(-0.15, "cm"), legend.key = element_rect(size = 5, color = 'NA'),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), units = "pt"))

plot_nsp12_mv_BA1 <- ggplot(nsp12_sum %>% filter(AA!="") %>% filter(AA!="<NA>") %>% filter(lineage_group == "BA.1"),
                            aes(x=RefSite, y=Prop_AA, colour=AA_rank)) + 
  geom_point(aes(size=AAcoverage), alpha=0.5) +
  facet_grid(tx_group~visit_day,
             labeller = labeller(tx_group = group.labs, visit_day=day.labs)) +
  ylim(-0.05,1.05) + 
  ylab("Proportion of amino acid") + xlab("Codon position (1000s)") +
  ggtitle("Nsp12") +
  theme_bw() +
  scale_colour_viridis_d(begin = 0, end = 0.8, direction=-1) +
  scale_x_continuous(labels = label_number(scale = 1e-3)) + 
  labs(colour="AA rank", size="Read depth") +
  guides(colour=guide_legend(order=1), size = guide_legend(order=2)) +
  geom_hline(yintercept=c(0.10, 0.5), linetype= "dashed", colour="grey") +
  scale_size_area(breaks=c(1000, 7500, 15000, 22500, 30000)) +
  theme(axis.text = element_text(colour="black", size=10),
        axis.title.y = element_text(vjust=2),
        axis.title = element_text(size=14),
        plot.title = element_text(face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(hjust=0.5, size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.box="vertical", legend.margin=margin(), legend.key.size =unit(1, "cm"),
        legend.spacing.y = unit(-0.15, "cm"), legend.key = element_rect(size = 5, color = 'NA'),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), units = "pt"))

#for git 
ggsave(plot_nsp12_mv_alpha, filename = "~/projects/agile_mpv_seq_project/Reports/nsp12_mv_alpha_agile.png", device="png", dpi=300, width=8, height=12)
ggsave(plot_nsp12_mv_delta, filename = "~/projects/agile_mpv_seq_project/Reports/nsp12_mv_delta_agile.png", device="png", dpi=300, width=8, height=12)
ggsave(plot_nsp12_mv_EU1, filename = "~/projects/agile_mpv_seq_project/Reports/nsp12_mv_EU1_agile.png", device="png", dpi=300, width=8, height=12)
ggsave(plot_nsp12_mv_BA1, filename = "~/projects/agile_mpv_seq_project/Reports/nsp12_mv_BA1_agile.png", device="png", dpi=300, width=8, height=12)

# Spike data ####
spike_sum <- all_aa_mv_m_final %>% filter(Protein=="S") %>% filter(AAcoverage>200)

#plot for each lineage
plot_spike_mv_alpha <- ggplot(spike_sum %>% filter(AA!="") %>% filter(AA!="<NA>") %>% filter(lineage_group == "alpha"),
                              aes(x=RefSite, y=Prop_AA, colour=AA_rank)) + 
  geom_point(aes(size=AAcoverage), alpha=0.5) +
  facet_grid(tx_group~visit_day,
             labeller = labeller(tx_group = group.labs, visit_day=day.labs)) +
  ylim(-0.05,1.05) + 
  ylab("Proportion of amino acid") + xlab("Codon position (1000s)") +
  ggtitle("Spike") +
  theme_bw() +
  scale_colour_viridis_d(begin = 0, end = 0.8, direction=-1) +
  scale_x_continuous(labels = label_number(scale = 1e-3)) + 
  labs(colour="AA rank", size="Read depth") +
  guides(colour=guide_legend(order=1), size = guide_legend(order=2)) +
  geom_hline(yintercept=c(0.10, 0.5), linetype= "dashed", colour="grey") +
  scale_size_area(breaks=c(1000, 7500, 15000, 22500, 30000)) +
  theme(axis.text = element_text(colour="black", size=10),
        axis.title.y = element_text(vjust=2),
        axis.title = element_text(size=14),
        plot.title = element_text(face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(hjust=0.5, size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.box="vertical", legend.margin=margin(), legend.key.size =unit(1, "cm"),
        legend.spacing.y = unit(-0.15, "cm"), legend.key = element_rect(size = 5, color = 'NA'),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), units = "pt"))

plot_spike_mv_delta <- ggplot(spike_sum %>% filter(AA!="") %>% filter(AA!="<NA>") %>% filter(lineage_group == "delta"),
                              aes(x=RefSite, y=Prop_AA, colour=AA_rank)) + 
  geom_point(aes(size=AAcoverage), alpha=0.5) +
  facet_grid(tx_group~visit_day,
             labeller = labeller(tx_group = group.labs, visit_day=day.labs)) +
  ylim(-0.05,1.05) + 
  ylab("Proportion of amino acid") + xlab("Codon position (1000s)") +
  ggtitle("Spike") +
  theme_bw() +
  scale_colour_viridis_d(begin = 0, end = 0.8, direction=-1) +
  scale_x_continuous(labels = label_number(scale = 1e-3)) + 
  labs(colour="AA rank", size="Read depth") +
  guides(colour=guide_legend(order=1), size = guide_legend(order=2)) +
  geom_hline(yintercept=c(0.10, 0.5), linetype= "dashed", colour="grey") +
  scale_size_area(breaks=c(1000, 7500, 15000, 22500, 30000)) +
  theme(axis.text = element_text(colour="black", size=10),
        axis.title.y = element_text(vjust=2),
        axis.title = element_text(size=14),
        plot.title = element_text(face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(hjust=0.5, size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.box="vertical", legend.margin=margin(), legend.key.size =unit(1, "cm"),
        legend.spacing.y = unit(-0.15, "cm"), legend.key = element_rect(size = 5, color = 'NA'),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), units = "pt"))

plot_spike_mv_EU1 <- ggplot(spike_sum %>% filter(AA!="") %>% filter(AA!="<NA>") %>% filter(lineage_group == "B.1.177/EU1"),
                            aes(x=RefSite, y=Prop_AA, colour=AA_rank)) + 
  geom_point(aes(size=AAcoverage), alpha=0.5) +
  facet_grid(tx_group~visit_day,
             labeller = labeller(tx_group = group.labs, visit_day=day.labs)) +
  ylim(-0.05,1.05) + 
  ylab("Proportion of amino acid") + xlab("Codon position (1000s)") +
  ggtitle("Spike") +
  scale_colour_viridis_d(begin = 0, end = 0.8, direction=-1) +
  theme_bw() +
  labs(colour="AA rank", size="Read depth") +
  guides(colour=guide_legend(order=1), size = guide_legend(order=2)) +
  geom_hline(yintercept=c(0.10, 0.5), linetype= "dashed", colour="grey") +
  scale_x_continuous(labels = label_number(scale = 1e-3)) + 
  theme(axis.text = element_text(colour="black", size=12, vjust=0.2),
        axis.title.y = element_text( vjust=2),
        axis.title = element_text(size=14),
        plot.title = element_text(face="bold"),
        strip.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(hjust=0.5, size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.box="vertical", legend.margin=margin(), legend.key.size =unit(1, "cm"),
        legend.spacing.y = unit(-0.15, "cm"), legend.key = element_rect(size = 5, color = 'NA'))

plot_spike_mv_BA1 <- ggplot(spike_sum %>% filter(AA!="") %>% filter(AA!="<NA>") %>% filter(lineage_group == "BA.1"),
                            aes(x=RefSite, y=Prop_AA, colour=AA_rank)) + 
  geom_point(aes(size=AAcoverage), alpha=0.5) +
  facet_grid(tx_group~visit_day,
             labeller = labeller(tx_group = group.labs, visit_day=day.labs)) +
  ylim(-0.05,1.05) + 
  ylab("Proportion of amino acid") + xlab("Codon position (1000s)") +
  ggtitle("Spike") +
  theme_bw() +
  scale_colour_viridis_d(begin = 0, end = 0.8, direction=-1) +
  scale_x_continuous(labels = label_number(scale = 1e-3)) + 
  labs(colour="AA rank", size="Read depth") +
  guides(colour=guide_legend(order=1), size = guide_legend(order=2)) +
  geom_hline(yintercept=c(0.10, 0.5), linetype= "dashed", colour="grey") +
  scale_size_area(breaks=c(1000, 7500, 15000, 22500, 30000)) +
  theme(axis.text = element_text(colour="black", size=10),
        axis.title.y = element_text(vjust=2),
        axis.title = element_text(size=14),
        plot.title = element_text(face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(hjust=0.5, size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.box="vertical", legend.margin=margin(), legend.key.size =unit(1, "cm"),
        legend.spacing.y = unit(-0.15, "cm"), legend.key = element_rect(size = 5, color = 'NA'),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), units = "pt"))

## write out 
#for git 
ggsave(plot_spike_mv_alpha, filename = "~/projects/agile_mpv_seq_project/Reports/spike_mv_alpha_agile.png", device="png", dpi=300, width=8, height=12)
ggsave(plot_spike_mv_delta, filename = "~/projects/agile_mpv_seq_project/Reports/spike_mv_delta_agile.png", device="png", dpi=300, width=8, height=12)
ggsave(plot_spike_mv_EU1, filename = "~/projects/agile_mpv_seq_project/Reports/spike_mv_EU1_agile.png", device="png", dpi=300, width=8, height=12)
ggsave(plot_spike_mv_BA1, filename = "~/projects/agile_mpv_seq_project/Reports/spike_mv_BA1_agile.png", device="png", dpi=300, width=8, height=12)


### figure making ####

# legend for all aa diversity plots - same scale used for all plots
legend_b <- get_legend(
  plot_spike_mv_delta + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", 
          legend.title = element_text(size=14),
          legend.text = element_text(size=14), legend.margin=margin(), legend.key.size =unit(1.5, "cm"))
)


fig2_abcd <- plot_grid(
    plot_wg_mv_delta + theme(legend.position = "none",
                             axis.title.x = element_blank(),
                             strip.text = element_text(size=10),
                             axis.title.y = element_text(size=12)),
    plot_nsp12_mv_delta + theme(axis.title.y = element_text(size=12),
                              strip.text = element_text(size=10),
                              legend.position = "none", 
                              axis.title.x = element_blank(),
                              plot.title=element_text(vjust=-3, face="plain")),
    plot_nsp14_mv_delta + theme(axis.title.x = element_blank(),
                                axis.title.y = element_text(size=12),
                                strip.text = element_text(size=10),
                                legend.position = "none",
                                plot.title=element_text(vjust=-3, face="plain")),
    plot_spike_mv_delta + theme(legend.position = "none",
                                axis.title.y = element_text(size=12),
                                strip.text = element_text(size=10),
                                plot.title=element_text(vjust=-3, face="plain")),
  labels = "auto", label_size = 14,
  rel_widths = c(1,1,1,1),rel_heights = c(1,1,1,1), scale=1,
  hjust = -0.5,
  nrow = 4, align = "vh", axis = "l")

fig2_wleg <- plot_grid(fig2_abcd, legend_b, ncol = 1, rel_heights = c(1,0.1), rel_widths = c(1, 0.5),
                       align = "vh", axis = "l")

# write out
save_plot(fig2_wleg, filename = ("~/projects/agile_mpv_seq_project/Reports/Agile_paper_figs/fig2_final.png"), base_width=10, base_height = 12, dpi=300, bg="white")

# Fig3 - all whole genome for rest of approx. equally matched lineages (alpha, EU1, BA.1)
fig3_wg_all <- plot_grid(
  plot_wg_mv_alpha + theme(legend.position = "none",
                           axis.title.x = element_blank(),
                           strip.text = element_text(size=10),
                           axis.title.y = element_text(size=12),
                           plot.title=element_text(vjust=-3, face="plain")),
  plot_wg_mv_EU1 + theme(axis.title.y = element_text(size=12),
                              strip.text = element_text(size=10),
                              legend.position = "none", 
                              axis.title.x = element_blank(),
                              plot.title=element_text(vjust=-3, face="plain")),
  plot_wg_mv_BA1 + theme(axis.title.y = element_text(size=12),
                              strip.text = element_text(size=10),
                              legend.position = "none",
                              plot.title=element_text(vjust=-3, face="plain")),
  labels = "auto", label_size = 14,
  rel_widths = c(1,1,1),rel_heights = c(1,1,1), scale=1,
  hjust = -0.5,
  nrow = 3, align = "vh", axis = "l")

fig3_wleg <- plot_grid(fig3_wg_all, legend_b, ncol = 1, rel_heights = c(1,0.1), rel_widths = c(1, 0.5),
                       align = "vh", axis = "l")

save_plot(fig3_wleg, filename = ("~/projects/agile_mpv_seq_project/Reports/Agile_paper_figs/fig3_final.png"), base_width=10, base_height = 12, dpi=300, bg="white")


### supplemental figures
# S2 - alpha nsp12, nsp14 and S
S2_alpha_nspS <- plot_grid(
  plot_nsp12_mv_alpha + theme(legend.position = "none",
                          axis.title.x = element_blank(),
                           strip.text = element_text(size=12),
                           axis.title.y = element_text(size=12),
                           plot.title=element_text(vjust=-1, face="plain")),
  plot_nsp14_mv_alpha + theme(axis.title.y = element_text(size=12),
                         strip.text = element_text(size=12),
                         legend.position = "none", 
                         axis.title.x = element_blank(),
                         plot.title=element_text(vjust=-1, face="plain")),
  plot_spike_mv_alpha + theme(axis.title.y = element_text(size=12),
                         strip.text = element_text(size=12),
                         legend.position = "none",
                         plot.title=element_text(vjust=-1, face="plain")),
  labels = "auto", label_size = 14,
  rel_widths = c(1,1,1),rel_heights = c(1,1,1), scale=1,
  hjust = -0.5,
  nrow = 3, align = "vh", axis = "l")

S2_wleg <- plot_grid(S2_alpha_nspS, legend_b, ncol = 1, rel_heights = c(1,0.1), rel_widths = c(1, 0.5),
                     align = "vh", axis = "l")

save_plot(S2_wleg, filename = ("~/projects/agile_mpv_seq_project/Reports/Agile_paper_figs/S2_alpha_aa-diversity_final.png"), base_width=10, base_height = 12, dpi=300, bg="white")

# S3 - EU1 nsp12, nsp14 and S
S3_EU1_nspS <- plot_grid(
  plot_nsp12_mv_EU1 + theme(legend.position = "none",
                              axis.title.x = element_blank(),
                              strip.text = element_text(size=12),
                              axis.title.y = element_text(size=12),
                              plot.title=element_text(vjust=-1, face="plain")),
  plot_nsp14_mv_EU1 + theme(axis.title.y = element_text(size=12),
                              strip.text = element_text(size=12),
                              legend.position = "none", 
                              axis.title.x = element_blank(),
                              plot.title=element_text(vjust=-1, face="plain")),
  plot_spike_mv_EU1 + theme(axis.title.y = element_text(size=12),
                              strip.text = element_text(size=12),
                              legend.position = "none",
                              plot.title=element_text(vjust=-1, face="plain")),
  labels = "auto", label_size = 14,
  rel_widths = c(1,1,1),rel_heights = c(1,1,1), scale=1,
  hjust = -0.5,
  nrow = 3, align = "vh", axis = "l")

S3_wleg <- plot_grid(S3_EU1_nspS, legend_b, ncol = 1, rel_heights = c(1,0.1), rel_widths = c(1, 0.5),
                     align = "vh", axis = "l")

save_plot(S3_wleg, filename = ("~/projects/agile_mpv_seq_project/Reports/Agile_paper_figs/S3_EU1_aa-diversity_final.png"), base_width=10, base_height = 12, dpi=300, bg="white")

# S4 - BA1 nsp12, nsp14 and S
S4_BA1_nspS <- plot_grid(
  plot_nsp12_mv_BA1 + theme(legend.position = "none",
                            axis.title.x = element_blank(),
                            strip.text = element_text(size=12),
                            axis.title.y = element_text(size=12),
                            plot.title=element_text(vjust=-1, face="plain")),
  plot_nsp14_mv_BA1 + theme(axis.title.y = element_text(size=12),
                            strip.text = element_text(size=12),
                            legend.position = "none", 
                            axis.title.x = element_blank(),
                            plot.title=element_text(vjust=-1, face="plain")),
  plot_spike_mv_BA1 + theme(axis.title.y = element_text(size=12),
                            strip.text = element_text(size=12),
                            legend.position = "none",
                            plot.title=element_text(vjust=-1, face="plain")),
  labels = "auto", label_size = 14,
  rel_widths = c(1,1,1),rel_heights = c(1,1,1), scale=1,
  hjust = -0.5,
  nrow = 3, align = "vh", axis = "l")

S4_wleg <- plot_grid(S4_BA1_nspS, legend_b, ncol = 1, rel_heights = c(1,0.1), rel_widths = c(1, 0.5),
                     align = "vh", axis = "l")

save_plot(S4_wleg, filename = ("~/projects/agile_mpv_seq_project/Reports/Agile_paper_figs/S4_BA1_aa-diversity_final.png"), base_width=10, base_height = 12, dpi=300, bg="white")
