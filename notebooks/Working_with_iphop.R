library(ggplot2)
library(dplyr)
library(janitor)
library(tidyverse)
library(utils)
library(ggvenn)
library(ggpubr)
library(stringi)
library(DataEditR)

###VIRSORTER2###
virsorter1 <-read.delim("../Data/run12special.txt", header = T, check.names = F)
virsorter2 <- virsorter1 %>% filter(!grepl("seqname", `139536_S27_vs2.tsv seqname`))
virsorter2$`139536_S27_vs2.tsv seqname` <- str_remove(virsorter2$`139536_S27_vs2.tsv seqname`, "_vs2.tsv")
virsorter3 <- virsorter2 %>% separate(`139536_S27_vs2.tsv seqname`, into = c("sampleID","contig"), sep = " ")
virsorter3$seq_name <- paste0(virsorter3$sampleID, "_" , virsorter3$contig)
virsorter3[3] <- sapply(virsorter3[3], as.numeric)
colnames(virsorter3)[3] <- "length"
sum(duplicated(virsorter3$seq_name))

#GENOMAD_PR
genomad_PR <- read.delim("../Data/PR/PR_genomad.tsv", header = T)
genomad_PR2 <- genomad_PR %>% filter(!grepl("seq_name", seq_name))
genomad_PR3 <- genomad_PR2 %>% separate(taxonomy.PR24_final_contigs_virus_summary.tsv, into = c("contig", "sampleID"), sep = " ")
genomad_PR3[2] <- sapply(genomad_PR3[2], as.numeric)
genomad_PR3$sampleID <- str_remove(genomad_PR3$sampleID, "_final_contigs_virus_summary.tsv")
genomad_PR3$seqname <- paste0(genomad_PR3$sampleID, "_", genomad_PR3$seq_name)
sum(duplicated(genomad_PR3$seqname))
genomad_PR4 <- genomad_PR3 %>% separate(seqname, into = c("seqname", "provirus"), sep = ";")

#write.table(genomad_PR4, "../Data/genomad_PR4.tsv", quote = F, row.names = F, sep = "\t")

###GENOMAD###
genomad <- read.delim("../Data/genomad_MG_SD.tsv", header = T)
genomad2 <- genomad %>% filter(!grepl("seq_name", seq_name))
genomad3 <- genomad2 %>% separate(taxonomy.139536_S27_genomad.tsv, into = c("contig", "sampleID"), sep = " ")
genomad3[2] <- sapply(genomad3[2], as.numeric)
genomad3$sampleID <- str_remove(genomad3$sampleID, "_genomad.tsv")
genomad3$seqname <- paste0(genomad3$sampleID, "_", genomad3$seq_name)
sum(duplicated(genomad3$seqname))

###Dereplication
derep_list <- read.delim("../Data/genomad_clusters_f1.tsv", header = F)

###Combine rows

genomad_bind <- bind_rows(genomad_PR4[,-14], genomad3)
genomad_bind_10k <- genomad_bind %>% filter(genomad_bind$length >= 10000)

###Filter 

genomad_filter <- genomad_bind_10k %>% filter(genomad_bind_10k$seqname %in% derep_list$V1) %>% distinct(seqname, .keep_all = T)
sum(duplicated(genomad_filter$seqname))


##Filter
colnames(genomad_filter)[2]

genomad_filter <- genomad_filter %>% mutate(classify = 
                                              ifelse(length < 50000 ,
                                                     "<50k",
                                                     ifelse(length >= 50000 & length <= 100000,
                                                            "100-200k",
                                                            ">1000k")))


###CheckV

checkv <- read.delim("../Data/quality_summary.tsv", header = T, check.names = F)
checkv_PR <- read.delim("../Data/PR/quality_summary.tsv", header = T, check.names = F) 

checkv_bind <- bind_rows(checkv, checkv_PR)
checkv_bind_sep <- checkv_bind %>% separate(contig_id, into = c("seqname", "provirus"), sep = ";")

derep_checkv_bind <- checkv_bind_sep %>% filter(checkv_bind_sep$seqname %in% derep_list$V1)
###Join columns
genomad_checkv <- left_join(genomad_filter, checkv_bind_sep, by = "seqname")


##########Filtering
genomad_checkv <- genomad_checkv %>% mutate(classify = 
                                              ifelse(length < 50000 ,
                                                     "<50k",
                                                     ifelse(length >= 50000 & length <= 100000,
                                                            "100-200k",
                                                            ">1000k")))

genomad_checkv <- genomad_checkv %>% mutate(checkv = 
                                              ifelse(completeness < 50 ,
                                                     "Low quality (<50%)",
                                                     ifelse(completeness >= 50 & completeness <= 90,
                                                            "Medium quality (50-90%)",
                                                            ifelse(completeness >= 100 & checkv_quality %in% c("Complete"),
                                                                   "Complete (100%)",
                                                                   "High quality (>90%)"))))
genomad_checkv[, 28][is.na(genomad_checkv[, 28])] <- "Not-determined"

genomad_checkv <- genomad_checkv %>% mutate(score = 
                                              ifelse(virus_score < 0.7 ,
                                                     "Low",
                                                     ifelse(virus_score >= 0.7 & virus_score<=0.9,
                                                            "Good",
                                                            "High")))
#Plot!!!!!
options(scipen = 100, digits = 4)
give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}
genomad_checkv$checkv_quality <- factor(genomad_checkv$checkv_quality, levels = c( "Not-determined", "Low-quality", "Medium-quality", "High-quality", "Complete"))
# auto positioning function
give.n <- function(x){
  return(c(y = median(x)*1.175, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

gvs_l <- ggplot(genomad_checkv, aes(x = checkv_quality ,y = length)) + 
  geom_boxplot(
    
    # custom boxes
    color="darkblue",
    fill="darkblue",
    alpha=0.2,
    
    # Notch?
    notch=T,
    notchwidth = 0.1,
    
    # custom outliers
    outlier.colour="black",
    outlier.fill="black",
    outlier.size=0.1 ) +
  theme_classic() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(),
        axis.text.y = element_text(size = 12)) +
  coord_flip()+
  ylim(0, 100000)+
  stat_summary(fun.data = give.n, geom = "text", fun = median, color = "black")

gvs_l

gvs_m <- ggplot(genomad_checkv, aes(x = topology ,y = length)) + 
  geom_boxplot(
    
    # custom boxes
    color="darkblue",
    fill="darkblue",
    alpha=0.2,
    
    # Notch?
    notch=FALSE,
    notchwidth = 0.5,
    
    # custom outliers
    outlier.colour="black",
    outlier.fill="black",
    outlier.size=1 ) +
  theme_classic() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect()) +
  coord_flip() +
  ylim(0, 100000)
#stat_summary(fun.data = give.n, geom = "text", fun = count, color = "black")

gvs_m

#ggsave("../Figures/fig1b.png", plot = gvs_m, dpi = 600, units = c("in"), width = 7, height = 3)



blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"),
    legend.position = "bottom",
    legend.title = element_blank()
  )

topology <- genomad_checkv %>% group_by(topology) %>% tally()

genomad_iphop2 <- genomad_iphop2 %>% separate(contig, into = c("Type", "Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")

family <- genomad_iphop2 %>% group_by(Class) %>% tally()

#####IPHOP
iphop_r214 <- read.csv("../Data/r214/Host_prediction_to_genome_m90.csv", header = T)
iphop_PR <- read.csv("../Data/PR/Host_prediction_to_genome_m90.csv", header = T)

iphop <- rbind(iphop_r214, iphop_PR)

iphop_derep <- iphop  %>% filter(iphop$Virus %in% derep_list$V1) 

iphop_split <- separate(iphop_derep ,Host.taxonomy,into = c("Domain", "Phylum","Class", "Order", "Family", "Genus", "Species"),sep = ";",remove = FALSE,extra = "merge")
iphop_split$Order.Genus <- factor(paste0("(", iphop_split$Order, ")", iphop_split$Genus))
iphop_split$Order.Genus <- gsub("o__", "", iphop_split$Order.Genus)
iphop_split$Order.Genus <- gsub("g__", " ", iphop_split$Order.Genus)

iphop_split %>% filter(!grepl("GB_GCA_|RS_GCF_", Host.genome)) %>% group_by(Domain) %>% tally()
iphop_split %>% group_by(Virus) %>% tally()

#write.table(iphop_split, "../Data/iphop_split.tsv", quote = F, row.names = F, sep = "\t")
######
hostphylum <-  iphop_split %>% filter(Domain == "d__Archaea") %>% group_by(Phylum) %>% tally()
H10 <- data.frame("Others", 7077)
colnames(H10)[1] <- "Phylum"
colnames(H10)[2] <- "n"

hp <- rbind(hostphylum[1:10,], H10[1,])

######
iphop_checkv <- left_join(iphop_split, derep_checkv_bind ,
                          by = c("Virus" = "seqname"))

genomad_iphop <- left_join(iphop_checkv, genomad_checkv[, 1:14] ,
                           by = c("Virus" = "seqname"))

genomad_iphop2 <- genomad_iphop %>% filter(!str_detect(contig,"Nucleocytoviricota"))

#write.table(genomad_iphop2, "../Data/genomad_iphop.tsv", quote = F, row.names = F, sep = "\t")

topology <- genomad_iphop2 %>% #filter(checkv_quality == "Not-determined") %>%  
  group_by(topology) %>% tally()
domain <- genomad_iphop2 %>% #filter(checkv_quality == "Not-determined") %>%  
  group_by(Domain) %>% tally()
host <- genomad_iphop2 %>% #filter(checkv_quality == "Not-determined") %>%  
  group_by(Host.genome) %>% tally()
provirus <- genomad_iphop2 %>% filter(Domain == "d__Archaea") %>%  
  filter(!proviral_length == "NA")
Thermoplasmatota <- genomad_iphop2 %>% filter(Phylum == "p__Thermoplasmatota") %>%  
  filter(checkv_quality == c("High-quality", "Complete"))
complete <- genomad_iphop2 %>%  filter(Domain == "d__Archaea") %>%
  filter(completeness >= 50) 

#write.table(complete, "../Data/TableS8.tsv", quote = F, row.names = F, sep = "\t")

#write.table(provirus, "../Data/TableS9.tsv", quote = F, row.names = F, sep = "\t")

######
pie <- ggplot(topology, aes(x="", y=log(n), fill=topology))+
  geom_bar(alpha = 0.8, width = 1, stat = 'identity') +
  coord_polar("y", start=0) +
  scale_fill_manual(values = c("No terminal repeats" = "#c75a8c",
                               "Low-quality" = "#d39132",
                               "DTR" = "#a361c7",
                               "ITR" = "#6587cd",
                               "Provirus" = "#5db754"
  )) + 
  theme_void() + 
  theme(axis.text.x=element_blank())
pie



################


pie <- ggplot(hp, aes(x="", y=log(n), fill=Phylum))+
  geom_bar(alpha = 0.8, width = 1, stat = 'identity') +
  coord_polar("y", start=0) +
  scale_fill_manual(values = c("#51afdf",
                               "#377cb1",
                               "#5197e8",
                               "#345ea0",
                               "#889fdd",
                               "#6088ec",
                               "#6777b3",
                               "#4067be",
                               "#2f65d0",
                               "#4a5589",
                               "#396ced")) + 
  blank_theme + 
  theme(axis.text.x=element_blank(),
        legend.position = "right")
pie

pie <- ggplot(hostphylum, aes(x="", y=log(n), fill=Phylum))+
  geom_bar(alpha = 0.8, width = 1, stat = 'identity') +
  coord_polar("y", start=0) +
  scale_fill_manual(values = c("#a1bb2c",
                               "#789432",
                               "#54701c",
                               "#9abb68",
                               "#7fc14f",
                               "#50ca36",
                               "#419e32",
                               "#387b34",
                               "#54cc7d",
                               "#3fa066",
                               "#3dc5ac")) + 
  blank_theme + 
  theme(axis.text.x=element_blank(),
        legend.position = "right")
pie

##############Archaea


Archaea_checkv <- genomad_iphop2 %>% filter(Domain == "d__Archaea")
Archaea_HQ <- Archaea_checkv %>% filter(contig_length >= 10000 & completeness >= 50)
Archaea_HQ$Phylum <- gsub("p__", "", Archaea_HQ$Phylum)

Archaea_checkv %>% group_by(checkv_quality == "High-quality") %>% tally()
###Count viruses
Archaea_count <- Archaea_HQ %>% group_by(Phylum) %>% tally()

Archaea_count$Phylum<- factor(Archaea_count$Phylum, levels = c( "Nanoarchaeota",  
                                                 "Asgardarchaeota", 
                                                 "Thermoproteota",
                                                 "Halobacteriota", 
                                                 "Thermoplasmatota", 
                                                 "Methanobacteriota"))
Archaea_genomad_checkv <- left_join(Archaea_checkv[,c(1,2,1,12,14,16,17,19,20,21,23,25)], genomad_checkv[,c(2,11:19,21)],
                                    by = c("Virus" = "seqname"))
##########Barplot
Ac0 <- ggplot(Archaea_HQ, aes(x = Order.Genus,  fill = Phylum, alpha = 0.5)) + geom_bar(colour="black")+
  theme_bw() + 
  facet_grid(Phylum ~ ., scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_text(size = 14),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_blank(),
        strip.background = element_blank( ),
        panel.background = element_blank(),
        # Hide panel borders and remove grid lines
        #panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  ylab("Number of viral contigs") +
  xlab("Host Taxonomy by genus") +
  scale_fill_manual(values = c("#b4a7d6","#6598c6",  "#ffe599",  "#45818e",  "#d5a6bd", "#b6d7a8"
                                        ))+
                                          scale_y_log10() +
  coord_flip() +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50,5))

Ac0

##############################Boxplot
Archaea_HQ$Phylum<- factor(Archaea_HQ$Phylum, levels = c( "Nanoarchaeota",  
                                                 "Asgardarchaeota", 
                                                 "Thermoproteota",
                                                 "Halobacteriota", 
                                                 "Thermoplasmatota", 
                                                 "Methanobacteriota"))


Ac1 <- ggplot(Archaea_HQ, aes(x = Order.Genus,y = contig_length/1000,  fill = Phylum, alpha = 0.5)) + geom_boxplot() + geom_jitter()+
  theme_bw() + 
  facet_grid(Phylum ~ ., scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_blank(),
        strip.background = element_blank( ),
        panel.background = element_blank(),
        # Hide panel borders and remove grid lines
        #panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  ylab("Viral genome length (kb)") +
  xlab("Host Taxonomy by genus") +
  #stat_summary(fun.data = give.n, geom = "text", fun = count, color = "darkblue", position = "fill") +
  #ggtitle("genomad-virsorter2 consensus") + 
  scale_fill_manual(values = c("#b4a7d6","#6598c6",  "#ffe599",  "#45818e",  "#d5a6bd", "#b6d7a8"
  ))+
                                          scale_y_log10() +
  coord_flip()

Ac1

############
pcm <- melt(Archaea_HQ[,c(5,14,19,20)])

pcm$Phylum<- factor(pcm$Phylum, levels = c( "Nanoarchaeota",  
                                                          "Asgardarchaeota", 
                                                          "Thermoproteota",
                                                          "Halobacteriota", 
                                                          "Thermoplasmatota", 
                                                          "Methanobacteriota"))
Ac2 <- ggplot(pcm, aes(x = Order.Genus,y = value,  fill = variable, alpha = 0.5)) + 
  geom_col(position = "fill") + 
  theme_bw() + 
  facet_grid(Phylum~ ., scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_blank(),
        strip.background = element_blank( ),
        panel.background = element_blank(),
        # Hide panel borders and remove grid lines
        #panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  ylab("Host:Viral genes") +
  xlab("Number of genes") +
  #stat_summary(fun.data = give.n, geom = "text", fun = count, color = "darkblue", position = "dodge") +
  #ggtitle("genomad-virsorter2 consensus") + 
  scale_fill_manual(values = c("darkgreen", "maroon"))+
                                          scale_y_log10() +
  coord_flip()

Ac2
