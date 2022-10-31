# Script created by I'ah Donovan-Banfield to parse DiversiTools nucleotide outputs for 
# data visualisation for figures 1c and 1d in the manuscript "Characterisation of SARS-CoV-2 genomic variations in response to molnupiravir treatment in the AGILE Phase IIa clinical trial."

#set wd as directory where your Diversitools outputs are. If there are any outputs that are empty due to the sample not running well, remove them, as empty files cause readin function to fail.
setwd("~/data/agile-nimagen/easyseq_bams/easyseq_forR/easyseq_DT_outputs/")

library(scales)
library(ggplot2)
library(tidyverse)
library(gggenes)
library(patchwork)
library(ggrepel)
library(reshape2)
library(ggpubr)
library(rstatix)
library(ggbreak)
library(ggsignif)
library(cowplot)
library(magick) 



##### Read in nucleotide Diversitools output files "sample-name_entropy.txt" ####
diversitools_readin_tvts <- function(filepath, file_pattern, file_list) {
  
  temp <- list.files(filepath, pattern=file_pattern) # creates list variable  with file names as they are in the directory
  file_list <- list() # creates empty list to add to in the for loop
  
  for (i in 1:length(temp)) { # read in all data files as DFs, add Tv and Ts ratio column to all 
    
    file_list[[i]] <- read.delim(temp[i]) 
    file_list[[i]]$Coverage <- as.numeric(file_list[[i]]$Coverage)
    file_list[[i]]$CntTs <- as.numeric(file_list[[i]]$CntTs)
    file_list[[i]]$CntTv <- as.numeric(file_list[[i]]$CntTv)
    file_list[[i]]$TsRatio <- file_list[[i]]$CntTs / file_list[[i]]$Coverage
    file_list[[i]]$TvRatio <- file_list[[i]]$CntTv / file_list[[i]]$Coverage
    file_list[[i]]$ins_cnt <- as.numeric(file_list[[i]]$ins_cnt)
    file_list[[i]]$del_cnt <- as.numeric(file_list[[i]]$del_cnt)
    file_list[[i]]$insRatio <- file_list[[i]]$ins_cnt / file_list[[i]]$Coverage
    file_list[[i]]$delRatio <- file_list[[i]]$del_cnt / file_list[[i]]$Coverage
    file_list[[i]][,6:16] <- sapply(file_list[[i]][,6:16], as.character)
    file_list[[i]]$Ncnt <- as.character(file_list[[i]]$Ncnt)
    
  }
  
  #these two lines will be specific to your DiversiTools output file naming structure
  myfilenames <- gsub("_entropy.txt", "", temp) #removes common strings
  names(file_list) <- myfilenames
  return(file_list)
}


##### read in files #### 
setwd("~/data/agile-nimagen/easyseq_bams/easyseq_forR/easyseq_DT_outputs")
mut_list <- diversitools_readin_tvts(filepath = "./", file_pattern = "entropy.txt", file_list = mut_list)

# adding metadata 
# You need to make a metadata csv file that shares a column for unique identifiers ie sample number or patient randomisation number 
# and all the other metadata you might want to separate your data by i.e. treatment group, time point, ct value, lineage
metadata <- read.csv("~/projects/agile_mpv_seq_project/sample_info/AGILE-nimagen-full-180_final-metadata.csv")

# minor changes to lineage names for better labelling later on
metadata$lineage_group <- gsub("alpha", "B.1.1.7 / Alpha", metadata$lineage_group)
metadata$lineage_group <- gsub("B.1.177/EU1", "B.1.177 / EU1", metadata$lineage_group)
metadata$lineage_group <- gsub("delta", "B.1.617.2 / Delta", metadata$lineage_group)
metadata$lineage_group <- gsub("BA.1", "BA.1 / Omicron", metadata$lineage_group)
metadata$lineage_group <- gsub("BA.2", "BA.2 / Omicron", metadata$lineage_group)

# changing DiversiTools Sample column name and removing unnecessary strings from descriptor
for (i in 1:length(mut_list)) {
  mut_list[[i]] <- rename(mut_list[[i]], kit_id=Sample)
  # if your file name structure is not <sampleID>_<nimagenIndex_entropy.txt this second line is not needed
  mut_list[[i]]$kit_id <- gsub("_[^_]*", "", mut_list[[i]]$kit_id)
}

#bind all data together, add metadata and convert to long form
all_mut <- bind_rows(mut_list, .id = "data.frame")
all_mut <- merge(metadata, all_mut, by = "kit_id")
melt_mut <- reshape2::melt(all_mut, na.rm=F, id.vars=c("Subject", "visit_day", "tx_group", "lineage_group", "kit_id",
                                                       "Position", "Coverage", "ct_group"), 
                           measure.vars=c("TsRatio", "TvRatio"))
# set class for columns 
melt_mut$variable <- as.factor(melt_mut$variable)
melt_mut$visit_day <- factor(melt_mut$visit_day, levels= c(1, 3, 5))
melt_mut$tx_group <- factor(melt_mut$tx_group, levels=c("Group B", "Group A"))
melt_mut$Subject <- factor(melt_mut$Subject)
melt_mut$lineage_group <- factor(melt_mut$lineage_group)
melt_mut$ct_group <- factor(melt_mut$ct_group, levels=c("> 28", "25-28", "20-24", "< 20"))
melt_mut$value <- as.numeric(melt_mut$value)

#for changing facet wrap labels
group.labs <- c( "placebo", "molnupiravir") 
names(group.labs) <- c("Group B", "Group A")

dodge <- position_dodge(width=1)# create variable for spacing of boxplots so they don't overlap

### get Ts/Tv ratio per sample ####
meanTs <- all_mut %>% filter(Coverage >200) %>% group_by(tx_group, lineage_group, visit_day, ct_group, kit_id) %>% 
  summarize(Transition=mean(TsRatio, na.rm=T ), stdev=sd(TsRatio, na.rm=T)) 

meanTv <- all_mut %>% filter(Coverage >200) %>% group_by(tx_group, lineage_group, visit_day, ct_group, kit_id) %>% 
  summarize(Transversion=mean(TvRatio, na.rm=T ), stdev=sd(TvRatio, na.rm=T)) 
#merge together
meanTsTv <- merge(meanTs, meanTv, by =c("kit_id","tx_group", "lineage_group", "visit_day","ct_group"))
meanTsTv$TsTvRatio <- meanTsTv$Transition / meanTsTv$Transversion
melt_meanTsTv <- reshape2::melt(meanTsTv, na.rm=F, id.vars=c("kit_id","tx_group", "lineage_group", "visit_day","ct_group"), measure.vars=c("TsTvRatio", "Transition", "Transversion"))
melt_meanTsTv$visit_day <- factor(melt_meanTsTv$visit_day, levels=c("1", "3", "5"))
melt_meanTsTv$variable <- factor(melt_meanTsTv$variable, levels=c("TsTvRatio", "Transversion", "Transition"))
melt_meanTsTv$tx_group <- factor(melt_meanTsTv$tx_group, levels=c("Group B", "Group A"))


# subset out just the TsTv ratios
melt_meanTsTv_ratio <- subset(melt_meanTsTv, variable=="TsTvRatio")
melt_meanTsTv_ratio <- rename(melt_meanTsTv_ratio, ratio=variable)

#generate stats [rstatix package]
stats.test.tvts <- melt_meanTsTv_ratio %>% group_by(visit_day, ratio) %>%
  wilcox_test(value~tx_group, comparisons = list(c("1", "3"), c("1", "5"), c("3", "5") )) %>% 
  adjust_pvalue(method="bonferroni") %>% add_significance("p.adj") %>% add_xy_position(x = "visit_day", dodge = 0.8)

plot_tvtsRatio_stats <- ggplot(melt_meanTsTv_ratio, 
                               aes(x=visit_day, y=value, colour=tx_group)) + 
  geom_boxplot(outlier.colour=NA, fill="white") +
  geom_point(position=position_jitterdodge(), size=1.2, alpha=0.5) +
  ylim(0.2, 1.24) +
  ylab("Mean Ts/Tv Ratio") +
  xlab("Day of swab sample") +
  labs(colour="Treatment allocation") +
  theme_bw() +
  # scale_colour_manual(values=c("#27AD81FF", "#404788FF"), labels=group.labs)+
  scale_colour_manual(values=c( "#7ad151ff","#440154FF"), labels=group.labs)+
  theme(axis.text = element_text(colour="black", size=12, hjust=0.5),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(7.5, 5.5, 0, 5.5), units = "pt")) +
        stat_pvalue_manual(stats.test.tvts,  label = "p.adj.signif", y.position = 1.24, tip.length = 0)
#write out 
#png for gitlab
ggsave(plot_tvtsRatio_stats, filename = "~/projects/agile_mpv_seq_project/Reports/Agile_paper_figs/mean-tvts_stats_colour-update.png", device="png", width = 8, height=6, dpi=300)
#tiff for pubs
ggsave(plot_tvtsRatio_stats, filename = "~/data/agile-nimagen/Reports/mean-tvts_stats_colour-update.tiff", device="tiff", width = 8, height=6, dpi=300)


#### individual base changes ####

#genome nucleotide position of proteins
nsps <- read.csv("/home/hannahg/projects/dstl_project/data/below_15%_N/protein_position.csv") 
names(nsps)[1]<-'Position'
sapply(nsps,class) #checking class of positions
nsps$Position<-as.numeric(nsps$Position) # change pos to numeric

all_mut<-merge(nsps,all_mut)
all_mut$Acnt <- as.numeric(all_mut$Acnt)
all_mut$Ccnt <- as.numeric(all_mut$Ccnt)
all_mut$Gcnt <- as.numeric(all_mut$Gcnt)
all_mut$Tcnt <- as.numeric(all_mut$Tcnt)

# mutations -  makes individual tibble for each nucleotide 
# G nucleotide 
all_mut_G <- subset(all_mut, RefBase == "G")
all_mut_G < as_tibble(all_mut_G)
# A nucleotide
all_mut_A <- subset(all_mut, RefBase == "A")
all_mut_A < as_tibble(all_mut_A)
# T nucleotide (changing to U for RNA)
all_mut_U <-subset(all_mut, RefBase =="T")
all_mut_U <- as_tibble(all_mut_U)
# C nucleotide
all_mut_C <-subset(all_mut, RefBase =="C")
all_mut_C <- as_tibble(all_mut_C)

## add columns for each base change 
all_mut_G <- all_mut_G %>% mutate(
  G_to_A = (Acnt / Coverage),
  G_to_C = (Ccnt / Coverage),
  G_to_U = (Tcnt / Coverage)
)

all_mut_A <- all_mut_A %>% mutate(
  A_to_G = (Gcnt / Coverage),
  A_to_C = (Ccnt / Coverage),
  A_to_U = (Tcnt / Coverage)
)

all_mut_U <- all_mut_U %>% mutate(
  U_to_G = (Gcnt / Coverage),
  U_to_A = (Acnt / Coverage),
  U_to_C = (Ccnt / Coverage)
)

all_mut_C <- all_mut_C %>%mutate(
  C_to_G = (Gcnt / Coverage),
  C_to_A = (Acnt / Coverage),
  C_to_U = (Tcnt / Coverage)
)

# melts mutation rows to be variables in "basechange" with their value in the next column
all_mut_A <-all_mut_A %>% gather(basechange, value, 45:47)
all_mut_C <-all_mut_C %>% gather(basechange, value, 45:47)
all_mut_U <-all_mut_U %>% gather(basechange, value, 45:47)
all_mut_G <-all_mut_G %>% gather(basechange, value, 45:47)

# bind together 
basechange_all <- bind_rows(all_mut_A, all_mut_C, all_mut_G, all_mut_U)

# filter only positions that have min 200 coverage
basechange_all_200 <- subset(basechange_all, Coverage > 200)

# class to factor for plotting
basechange_all_200$visit_day <- factor(basechange_all_200$visit_day, levels=c("1", "3", "5"))
basechange_all_200$tx_group <- factor(basechange_all_200$tx_group, levels=c("Group B", "Group A"))


# getting mean basechange value per sample for neater plotting
meanGU <- basechange_all_200 %>% filter(basechange=="G_to_U") %>% group_by(tx_group, lineage_group, visit_day, ct_group, kit_id) %>% 
  summarize(G_to_U=mean(value, na.rm=T ), stdev=sd(value, na.rm=T))
meanGC <- basechange_all_200 %>% filter(basechange=="G_to_C") %>% group_by(tx_group, lineage_group, visit_day, ct_group, kit_id) %>% 
  summarize(G_to_C=mean(value, na.rm=T ), stdev=sd(value, na.rm=T))
meanGA <- basechange_all_200 %>% filter(basechange=="G_to_A") %>% group_by(tx_group, lineage_group, visit_day, ct_group, kit_id) %>% 
  summarize(G_to_A=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 

meanCU <- basechange_all_200 %>% filter(basechange=="C_to_U") %>% group_by(tx_group, lineage_group, visit_day, ct_group, kit_id) %>% 
  summarize(C_to_U=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 
meanCA <- basechange_all_200 %>% filter(basechange=="C_to_A") %>% group_by(tx_group, lineage_group, visit_day, ct_group, kit_id) %>% 
  summarize(C_to_A=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 
meanCG <- basechange_all_200 %>% filter(basechange=="C_to_G") %>% group_by(tx_group, lineage_group, visit_day, ct_group, kit_id) %>% 
  summarize(C_to_G=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 

meanAU <- basechange_all_200 %>% filter(basechange=="A_to_U") %>% group_by(tx_group, lineage_group, visit_day, ct_group, kit_id) %>% 
  summarize(A_to_U=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 
meanAC <- basechange_all_200 %>% filter(basechange=="A_to_C") %>% group_by(tx_group, lineage_group, visit_day, ct_group, kit_id) %>% 
  summarize(A_to_C=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 
meanAG <- basechange_all_200 %>% filter(basechange=="A_to_G") %>% group_by(tx_group, lineage_group, visit_day, ct_group, kit_id) %>% 
  summarize(A_to_G=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 

meanUC <- basechange_all_200 %>% filter(basechange=="U_to_C") %>% group_by(tx_group, lineage_group, visit_day, ct_group, kit_id) %>% 
  summarize(U_to_C=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 
meanUG <- basechange_all_200 %>% filter(basechange=="U_to_G") %>% group_by(tx_group, lineage_group, visit_day, ct_group, kit_id) %>% 
  summarize(U_to_G=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 
meanUA <- basechange_all_200 %>% filter(basechange=="U_to_A") %>% group_by(tx_group, lineage_group, visit_day, ct_group, kit_id) %>% 
  summarize(U_to_A=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 

#put back together for plotting 
mean_GCGU <- merge(meanGC, meanGU,
                             by=c("kit_id","tx_group", "lineage_group", "visit_day","ct_group"))
mean_CACG <- merge(meanCA, meanCG,
                   by=c("kit_id","tx_group", "lineage_group", "visit_day","ct_group"))
mean_UCAU <- merge(meanUC, meanAU,
                   by=c("kit_id","tx_group", "lineage_group", "visit_day","ct_group"))
mean_AGAC <- merge(meanAG, meanAC,
                   by=c("kit_id","tx_group", "lineage_group", "visit_day","ct_group"))
mean_UAUG <- merge(meanUA, meanUG,
                   by=c("kit_id","tx_group", "lineage_group", "visit_day","ct_group"))
meanGACU <- merge(meanGA, meanCU, 
                  by =c("kit_id","tx_group", "lineage_group", "visit_day","ct_group"))


mean_basechange <- merge(merge(merge(merge(merge(
  mean_GCGU, 
  mean_CACG, by=c("kit_id","tx_group", "lineage_group", "visit_day","ct_group")),
  mean_UCAU, by=c("kit_id","tx_group", "lineage_group", "visit_day","ct_group")),
  mean_AGAC, by=c("kit_id","tx_group", "lineage_group", "visit_day","ct_group")),
  mean_UAUG, by=c("kit_id","tx_group", "lineage_group", "visit_day","ct_group")),
  meanGACU, by=c("kit_id","tx_group", "lineage_group", "visit_day","ct_group"))

melt_allbc <- reshape2::melt(mean_basechange, na.rm=F, id.vars=c("kit_id","tx_group", "lineage_group", "visit_day","ct_group"), 
                             measure.vars=c("G_to_A", "G_to_U", "G_to_C", "C_to_G", "C_to_A", "C_to_U", 
                                            "A_to_G", "A_to_U", "A_to_C", "U_to_A", "U_to_G", "U_to_C"))
melt_allbc$visit_day <- factor(melt_allbc$visit_day, levels=c("1", "3", "5"))
melt_allbc$variable <- factor(melt_allbc$variable)
melt_allbc$tx_group <- factor(melt_allbc$tx_group, levels=c("Group B", "Group A"))
melt_allbc <- rename(melt_allbc, basechange=variable)
melt_allbc$basechange <- factor(melt_allbc$basechange, levels=c("G_to_A", "G_to_U", "G_to_C", "C_to_U", "C_to_A", "C_to_G", 
  "A_to_G", "A_to_U", "A_to_C", "U_to_A", "U_to_G", "U_to_C"))
melt_allbc$basechange <- gsub("_to_", " > ", melt_allbc$basechange)

### plotting ###
# generate stats
stat.test.base.all <- melt_allbc %>%
  group_by(visit_day, basechange, comparisons=list(c("Group B", "Group A"))) %>% #have to also group by the basechange otherwise you do not get individual tests
  wilcox_test(value ~ tx_group) %>% adjust_pvalue(method = "bonferroni") %>% 
  add_significance("p.adj") %>% add_xy_position(x = "visit_day", dodge = 0.8)


base.change.all.plot <- ggplot(melt_allbc, aes(x=visit_day, y=value, colour= tx_group)) +
  geom_boxplot(outlier.colour=NA) +
  geom_point(position=position_jitterdodge(), size=0.5, alpha=0.5) +
  ylab("Mean proportion of base changes") + xlab("Day of swab sample") +
  theme(axis.text.x = element_text( vjust = 0.5, hjust=1))+
  facet_wrap(~basechange, ncol=4, dir = "v") + 
  scale_y_log10() + ylim(0,0.0065) +
  theme_bw() +
  scale_colour_manual(values=c( "#7ad151ff","#440154FF"), labels=group.labs)+
  theme(axis.text = element_text(colour="black", size=12, hjust=0.5),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        strip.background = element_rect(size=1, linetype="solid", fill="white"),
        strip.text.x = element_text(size = 12, color = "black", face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # panel.border = element_rect(size=1),
        # axis.line = element_line(colour = "black"),
        legend.title = element_text(hjust=0.5),
        legend.position="bottom",
        plot.margin = unit(c(7.5, 5.5, 0, 5.5), units = "pt")) +
  labs(colour="Treatment allocation") +
  stat_pvalue_manual(
    stat.test.base.all,  label = "p.adj.signif", tip.length = 0.0, y.position = 0.0062)

ggsave(base.change.all.plot, filename = "~/projects/agile_mpv_seq_project/Reports/Agile_paper_figs/S1_allbase_jit.pdf", device="pdf", dpi=300, width=10, height=8)

# isolating GA, AG, CU and UC (molnupiravir changes) for fig 1 
meanGACU <- merge(meanGA, meanCU, by =c("kit_id","tx_group", "lineage_group", "visit_day","ct_group"))
meanAGAC <- merge(meanAG, meanUC,
                   by=c("kit_id","tx_group", "lineage_group", "visit_day","ct_group"))
mean_mpv_changes <- merge(meanGACU, meanAGAC,
                          by=c("kit_id","tx_group", "lineage_group", "visit_day","ct_group"))

melt_mean.mpv <- reshape2::melt(mean_mpv_changes , na.rm=F, id.vars=c("kit_id","tx_group", "lineage_group", "visit_day","ct_group"), measure.vars=c("G_to_A", "A_to_G","C_to_U", "U_to_C"))
melt_mean.mpv$visit_day <- factor(melt_mean.mpv$visit_day, levels=c("1", "3", "5"))
melt_mean.mpv$variable <- gsub("_to_", " > ", melt_mean.mpv$variable )
melt_mean.mpv$variable <- factor(melt_mean.mpv$variable, levels=c("G > A", "A > G", "C > U","U > C"))
melt_mean.mpv$tx_group <- factor(melt_mean.mpv$tx_group, levels=c("Group B", "Group A"))
melt_mean.mpv <- rename(melt_mean.mpv, basechange=variable)

# get stats for comparing treatment groups 
stat.test.base <- melt_mean.mpv %>%
  group_by(visit_day, basechange, comparisons=list(c("Group B", "Group A"))) %>% #have to also group by the basechange otherwise you do not get individual tests
  wilcox_test(value ~ tx_group) %>% adjust_pvalue(method = "bonferroni") %>% #may have to check this appropriate stats test
  add_significance("p.adj") %>% add_xy_position(x = "visit_day", dodge = 0.8)

# plotting for figure 1d
mpv.base.change.plot <- ggplot(melt_mean.mpv, aes(x=visit_day, y=value, colour= tx_group)) +
  geom_boxplot(outlier.colour=NA) +
  geom_point(position=position_jitterdodge(), size=0.5, alpha=0.5) +
  facet_wrap(~basechange, dir= "h") + #, labeller = labeller(basechange = base.labs)
  ylab("Mean proportion of base changes") + xlab("Day of swab sample") + ylim(0, 0.0065) +
  theme(axis.text.x = element_text( vjust = 0.5, hjust=1))+
  theme_bw() +
  scale_colour_manual(values=c( "#7ad151ff","#440154FF"), labels=group.labs)+
  theme(axis.text = element_text(colour="black", size=12, hjust=0.5),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        strip.background = element_rect(size=1, linetype="solid", fill="white"),
        strip.text.x = element_text(size = 12, color = "black", face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # panel.border = element_rect(size=1),
        # axis.line = element_line(colour = "black"),
        legend.title = element_text(hjust=0.5),
        legend.position="bottom",
        plot.margin = unit(c(7.5, 5.5, 0, 5.5), units = "pt")) +
  labs(colour="Treatment allocation") +
  stat_pvalue_manual(
    stat.test.base,  label = "p.adj.signif", tip.length = 0.0, y.position = 0.0062 )

#write out
#for git 
ggsave(mpv.base.change.plot, filename = ("~/projects/agile_mpv_seq_project/Reports/mpv-basechange-boxplot_mean.png"), device="png", dpi=300, width=8,height=6)


#### Figure creation ####

fig1_a <- cowplot::ggdraw() + 
  cowplot::draw_image("~/projects/agile_mpv_seq_project/Reports/nimagen-pipeline_revised.png", scale = 1) 

fig1_b <- cowplot::ggdraw() +
  cowplot::draw_image("~/projects/agile_mpv_seq_project/Reports/molnu_template_graphic.png", scale = 0.95)

fig1_ab <- plot_grid(fig1_a, fig1_b, 
                     labels=c("a", "b"), label_size = 14, scale=1,
                     hjust = -1, vjust = 1.5,
                     nrow = 1,
                     align = "vh", axis = "b",
                     rel_widths = c(1, 0.4), rel_heights = c(1,1))

fig1_cd <-  plot_grid(
  plot_tvtsRatio_stats + theme(legend.position="none"),
  mpv.base.change.plot + theme(legend.position="none"),
  labels = c("c", "d"), label_size = 14, scale=0.95,
  hjust = -2, vjust = 1,
  nrow = 1,
  align = "vh", axis = "b", rel_widths = c(0.5, 1)
)

legend_cd <- get_legend(
  mpv.base.change.plot + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", 
          legend.title = element_text(size=12),
          legend.text = element_text(size=12), legend.margin=margin(), legend.key.size =unit(1, "cm"))
)

fig1_cd_leg <- plot_grid(fig1_cd, legend_cd, ncol = 1, rel_heights = c(1, .1))



fig1_all <- plot_grid(fig1_ab, fig1_cd_leg,
                     hjust = -1, vjust = 1,
                     nrow = 2, 
                     align = "vh", axis = "b", 
                     rel_widths = c(1, 1.4), rel_heights = c(1, 1.4)
)

save_plot(fig1_all, filename = "~/projects/agile_mpv_seq_project/Reports/Agile_paper_figs/fig1_panel-jit_revised.pdf", dpi=300, base_width=8, base_height = 10, device="pdf", bg="white")
