library(ggplot2)
library(dplyr)
library(janitor)
library(tidyverse)
library(utils)
library(ggvenn)
library(ggpubr)
library(stringi)
library(DataEditR)
library(ggbreak)
library(reshape2)

#Working wth genomad results from hybrid metagenomes
genomad_PR <- read.delim("../Data/PR/PR_genomad.tsv", header = T)
genomad_PR2 <- genomad_PR %>% filter(!grepl("seq_name", seq_name))
genomad_PR3 <- genomad_PR2 %>% separate(taxonomy.PR24_final_contigs_virus_summary.tsv, into = c("contig", "sampleID"), sep = " ")
genomad_PR3[2] <- sapply(genomad_PR3[2], as.numeric)
genomad_PR3$sampleID <- str_remove(genomad_PR3$sampleID, "_final_contigs_virus_summary.tsv")
genomad_PR3$seqname <- paste0(genomad_PR3$sampleID, "_", genomad_PR3$seq_name)
sum(duplicated(genomad_PR3$seqname))
genomad_PR4 <- genomad_PR3 %>% separate(seqname, into = c("seqname", "provirus"), sep = ";")

#Working with genomad results from short read metagenomes
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

gvs_l <- ggplot(genomad_checkv, aes(x = checkv_quality ,y = length)) + 
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
  coord_flip()+
  ylim(0, 100000)+
   stat_summary(fun.data = give.n, geom = "text", fun = count, color = "black")

gvs_l

#ggsave("../Figures/fig1c.png", plot = gvs_l, dpi = 600, units = c("in"), width = 7, height = 3)

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

genomad_checkv <- genomad_checkv %>% separate(contig, into = c("Type", "Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")

family <- genomad_checkv %>% group_by(Class) %>% tally()

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


######
hostphylum <-  iphop_split %>% filter(Domain == "d__Bacteria") %>% group_by(Phylum) %>% tally()
H10 <- data.frame("Others", 7077)
colnames(H10)[1] <- "Phylum"
colnames(H10)[2] <- "n"

hp <- rbind(hostphylum[1:10,], H10[1,])

######

topology <- genomad_checkv %>% #filter(checkv_quality == "Not-determined") %>%  
  group_by(topology) %>% tally()

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
                                          blank_theme + 
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

#ggsave("../Figures/bacteriaphylum.png", plot = pie, dpi = 600, units = c("in"), width = 5, height = 5)

library(ggplot2)
cc <- ggplot(genomad_checkv, aes(x = viral_genes, y = host_genes)) +
  geom_point(aes(color = factor(classify)), position = "jitter",size=2, alpha = 0.7) +  
  scale_colour_manual(values=c("<50k"  = "darkblue" ,"100-200k"="#a353d7", ">1000k" = "darkgreen")) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank()) +
  xlab("Viral_genes") +
  ylab("Host_genes") +
  labs(color= "Genome quality") + 
  xlim(0,100) + ylim(0,30)+
  geom_smooth(color="red", formula = y ~ x) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  guides(colour = guide_legend(override.aes = list(size=5)))
cc



###############Archaea
iphop_checkv <- left_join(iphop_split, derep_checkv_bind ,
                          by = c("Virus" = "seqname"))

Archaea_checkv <- iphop_checkv %>% filter(Domain == "d__Archaea")
Archaea_HQ <- Archaea_checkv %>% filter(contig_length >= 10000 & completeness >= 50)
Archaea_HQ$Phylum <- gsub("p__", "", Archaea_HQ$Phylum)


###Count viruses
Archaea_count <- Archaea_HQ %>% group_by(Phylum) %>% tally()

Archaea_count$Phylum<- factor(Archaea_count$Phylum, levels = c( "Nanoarchaeota",  
                                                 "Asgardarchaeota", 
                                                 "Thermoproteota",
                                                 "Halobacteriota", 
                                                 "Thermoplasmatota", 
                                                 "Methanobacteriota"))
##########Barplot
Ac0 <- ggplot(Archaea_count, aes(x = Order.Genus,  fill = Phylum, alpha = 0.5)) + geom_bar(colour="black")+
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
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30,5))

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

###Read virus and host abundance

virus_tpm <- read.delim("../../../Viruses/Feb2024/filtered/viruses/virus_tpm_table.tsv", header = T, row.names = 1)
host_tpm <- read.delim("../../../Viruses/Feb2024/filtered/host/host_tpm_table.tsv", header = T)

mean_virus_tpm <- data.frame(rowMeans(virus_tpm))
mean_virus_tpm$Virus = rownames(mean_virus_tpm)
#mean_host_tpm <- data.frame(rowMeans(host_tpm))


join_Avir <- mean_virus_tpm %>% filter( mean_virus_tpm$Virus %in% Archaea_HQ$Virus)
Avir <- left_join(Archaea_HQ, join_Avir, by = c("Virus"))
Archaea_HQ$Host.genome <- gsub("GB_", "", Archaea_HQ$Host.genome)
Archaea_HQ$Host.genome <- gsub("RS_", "", Archaea_HQ$Host.genome)

br_headers <- read.delim("../../../Viruses/Feb2024/filtered/host/BR_headers_filename.tzt", header = F)
br_headers$V2 <- br_headers$V1
br_headers <- separate(br_headers, V2, into = c("V2", "contig"), sep = ":" )
br_headers$V2 <- gsub(".fa", "", br_headers$V2)
br_headers$V1 <- gsub(".fa:>", "_", br_headers$V1)

gtdb_headers <- read.delim("../../../Viruses/Feb2024/filtered/host/GTDB_headers.txt", header = F)
gtdb_headers$V2 <- gtdb_headers$V1
gtdb_headers$V2 <- substr(gtdb_headers$V2, 1, 16)
gtdb_headers$V2 <- gsub(">", "", gtdb_headers$V2)
gtdb_headers$V1 <- gsub(">", "", gtdb_headers$V1)
gtdb_headers$V1 <- sub(" .*", "", gtdb_headers$V1)

headers <- rbind(br_headers[, c(1:2)], gtdb_headers)

############Read host

host_tpm2 <- left_join(host_tpm, headers, by = c("Contig" = "V1"))
host_tpm3 <- host_tpm2[, -1] %>% group_by(V2) %>% summarise(across(everything(), list(sum)))
mean_host_tpm <- data.frame(rowMeans(host_tpm3[,-1]))
mean_host_tpm$Genome <- host_tpm3$V2


join_Ahost <- mean_host_tpm %>% filter(mean_host_tpm$Genome %in% Archaea_HQ$Host.genome)
Avir_host <- left_join(Avir, join_Ahost, by = c("Host.genome" = "Genome"))
################


Ac4 <- ggarrange(Ac0, Ac1, Ac2, nrow = 1)
Ac4
ggsave("../Figures/fig2b.png", plot = Ac4, dpi = 600, units = c("in"), width = 15, height = 7)



######################
Avir_host$ratio <- Avir_host$rowMeans.host_tpm3....1../Avir_host$rowMeans.virus_tpm.
pcm_3 <- melt(Avir_host[,c(5,14,28,29)])

Ac5 <- ggplot(pcm_3, aes(x = Order.Genus,y = value,  fill = Phylum, pattern = variable)) + 
  geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", pattern_fill = "white", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  theme_bw() + 
  facet_grid( Phylum ~ ., scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_text(angle = 0, size = 13),
        strip.background = element_rect(fill = "white" ),
        panel.background = element_rect(fill = "white"),
        # Hide panel borders and remove grid lines
        #panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  ylab("Viral abundance") +
  xlab("Number of genes") +
  #stat_summary(fun.data = give.n, geom = "text", fun = count, color = "darkblue", position = "dodge") +
  #ggtitle("genomad-virsorter2 consensus") + 
  scale_fill_manual(values = c("#6598c6",
                                        "#45818e",
                                        "#b6d7a8",
                                        "#b4a7d6",
                                        "#d5a6bd",
                                        "#ffe599"))+
                                          scale_y_log10() +
  scale_pattern_manual(values = c(rowMeans.virus_tpm. = "none", 	
                                  rowMeans.host_tpm3....1.. = "stripe"))+
  coord_flip()

Ac5
###########################

Ac4 <- Ac0 + Ac1 + Ac2 #ggarrange(Ac0, Ac1, Ac2, nrow = 1)

ggsave("../Figures/fig2b.png", plot = Ac4, dpi = 600, units = c("in"), width = 15, height = 7)


#############
hostphylum <-  Archaea_HQ %>% group_by(Phylum) %>% tally()
H10 <- data.frame("Others", 7077)
colnames(H10)[1] <- "Phylum"
colnames(H10)[2] <- "n"

hp <- rbind(hostphylum[1:10,], H10[1,])


Avirus <-  Archaea_HQ %>% group_by(Virus) %>% tally()
Avirus_genomad_checkv <- genomad_checkv %>% filter(genomad_checkv$seqname %in% Avirus$Virus)
###############Adding abundance


write.table(iphop_derep, file = "../Data/iphop_derep_output.tsv", sep = "\t", quote = F, row.names = F)

write.table(Avirus, file = "../Data/Avirus_output_90x.tsv", sep = "\t", quote = F, row.names = F)

viral_abundance <- read.delim("../Data/count_table.tsv", header = T )


###################Readmapping
virus = iphop_derep %>% group_by(Virus) %>% tally()
host = iphop_derep %>% group_by(Host.genome) %>% tally()

genomad50 = genomad_checkv %>% filter(completeness >= 50) %>% group_by(seqname) %>% tally()
genomad_host50 = iphop_derep %>% filter(Virus %in% genomad50$seqname) %>% group_by(Host.genome) %>% tally()


write.table(genomad50 , file = "../Data/genomad50.tsv", sep = "\t", quote = F, row.names = F)
write.table(genomad_host50 , file = "../Data/genomad_host50.tsv", sep = "\t", quote = F, row.names = F)

##################################
colnames(Avir_host)[28] <- "Virus abundance (TPM)"
colnames(Avir_host)[29] <- "Host abundance (TPM)"

pcm_3 <- melt(Avir_host[,c(5,14,28,29)])
pcm_3$Phylum <- factor(pcm_3$Phylum, levels = c( "Nanoarchaeota",  
                                                 "Asgardarchaeota", 
                                                 "Thermoproteota",
                                                 "Halobacteriota", 
                                                 "Thermoplasmatota", 
                                                 "Methanobacteriota"))


Ac3 <- ggplot(pcm_3, aes(x = Order.Genus,y = value,  fill = variable, alpha = 0.5)) + 
  geom_boxplot()+ 
  theme_bw() + 
  facet_grid(Phylum  ~ ., scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_text(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_text(angle = 0, size = 13),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        strip.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent', color=NA),
       # Hide panel borders and remove grid lines
        #panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  ylab("Host abundance") +
  xlab("Number of genes") +
  #stat_summary(fun.data = give.n, geom = "text", fun = count, color = "darkblue", position = "dodge") +
  #ggtitle("genomad-virsorter2 consensus") + 
  scale_fill_manual(values = c("darkgreen","maroon"))+
  scale_y_log10() +
  coord_flip()

Ac3

ggsave("../Figures/fig2c.png", bg = "transparent", plot = Ac3, dpi = 600, units = c("in"), width = 8, height = 8)

################
min = Archaea_HQ %>% select(Order.Genus, completeness) %>% group_by(Order.Genus)  %>% summarise(across(where(is.numeric), ~ min(.x, na.rm = TRUE)))
max = Archaea_HQ %>% select(Order.Genus, completeness) %>% group_by(Order.Genus)  %>% summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)))

prophages = Archaea_HQ %>% select(Order.Genus, completeness) %>% group_by(Order.Genus)  %>% summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)))

##########For archaea
write.table(Archaea_HQ , file = "../Data/Archaea_HQcheckv_50.tsv", sep = "\t", quote = F, row.names = F)

write.table(Archaea_checkv, file = "../Data/Archaea_checkv.tsv", sep = "\t", quote = F, row.names = F)

#############Archaea_crispr
archaea_crispr <- read.delim("../Data/crispr_archea.tsv.txt", header = T)

archaea_crispr_hq <- Archaea_HQ %>% filter(Archaea_HQ$Virus %in% archaea_crispr$Virus)

archaea_checkv_hq <- Archaea_checkv %>% filter(Archaea_checkv$Virus %in% archaea_crispr$Virus)

archaea_crispr_genomad <- genomad_checkv %>% filter(genomad_checkv$seqname %in% archaea_crispr$Virus) %>% filter(completeness >= 100) %>% 
  filter(!grepl("Phycodnaviridae", Family))

archaea_crispr_genomad_iphop <- archaea_crispr %>% filter(archaea_crispr$Virus %in%   archaea_crispr_genomad$seqname)

annotations_genomad <- read.delim("../Data/annotations/virus_genes_out.txt", header = T)

archaeaviruses <- Archaea_crispr %>% group_by(Virus) %>% tally()

capsid <- annotations_genomad %>% filter(grepl("Major capsid protein" , annotation_description)) %>% filter(taxname == "Caudoviricetes") 

capsid$Gene <- capsid$X139536_S27.gene

capsid$Gene <- gsub("_17", "", capsid$Gene)        

crisprnoderep  <- read.delim("../Data/crispr-noderep.tsv", header = F)

checkv_crispr_noderep <- checkv %>% filter(checkv$contig_id %in% crisprnoderep$V1) %>% filter(completeness >= 100)

checkv_crispr_noderep_iphop <- iphop_checkv %>% filter(iphop_checkv$Virus %in% checkv_crispr_noderep$contig_id )

###################################Virus_hsot abundance 
virustpm <- read.delim("../Data/New folder/virus_tpm_table.tsv", header = T)
filtered <- c("139572_S173_NODE_115_length_108255_cov_48.245018", "139544_COM_NODE_136_length_101154_cov_77.364376", 
                "139581_S182_NODE_221_length_79954_cov_16.268514", "PR24_contig_10065_pilon")
virustpm_filtered <- virustpm %>% filter(virustpm$Contig %in% filtered)
colnames(virustpm_filtered) <- gsub("viruses_3097x_sr.mmi.", "", colnames(virustpm_filtered))
colnames(virustpm_filtered) <- gsub("_R1_001.fastq.gz_filtered.TPM", "", colnames(virustpm_filtered))
virustpm_f <- virustpm_filtered[, -c(1,65:70)]

rownames(virustpm_f) <- virustpm_filtered$Contig
metadata <- read.delim("../Data/New folder/metadata.txt", header = T)

colnames(virustpm_f) <- metadata$Sample_Season
virustpm_f$contig <- rownames(virustpm_f)

pcm = melt(virustpm_f)

plot_virus <- ggplot(pcm, aes(x = variable, y = value, fill = contig)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c("orange",
                                        "maroon",
                                        "darkgreen",
                                        "darkblue"
                                        )) + 
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 90)
  )

plot_virus

######################################################################
###########Co-occurrence###########################################
#####################################################################
library(igraph)
library(tidyverse)

sp_ratio<- read.delim("../Data/virus_ref_coverage/virus_ref_tpm_table.tsv", header = T)
colnames(sp_ratio) <- gsub("allgenomes_nodup_sl_sr.mmi.", "", colnames(sp_ratio))
colnames(sp_ratio) <- gsub("_R1_001.fastq.gz_filtered.TPM", "", colnames(sp_ratio))

taxa_table <- read.delim("../Data/virus_ref_coverage/taxa_table.txt", header = T, row.names = 1)
metadata <- read.delim("../Data/virus_ref_coverage/env_data.txt", header = T, row.names = 1)

sp_ratio <- sp_ratio %>% 
  as_tibble() %>% 
column_to_rownames(var = "Species")
min.prevalence=1
incidence=sp_ratio
incidence[incidence>0]=1
sp_ratio_filtered <- sp_ratio[which(rowSums(incidence)>=min.prevalence),] ### end of prevalence filtering

sp_correl <- sp_ratio_filtered %>% 
  t() %>% 
  cor(method = "spearman")  ### correlation calculations
sp_correl[abs(sp_correl)<0.65]=0 
net_work <- graph_from_adjacency_matrix(sp_correl,mode="lower",weighted=TRUE, diag=FALSE)  ## create network file
net_work <- delete.vertices(net_work,degree(net_work)==0) #remove nodes without edges

plot(net_work, vertex.label = NA, edge.width = 5, vertex.size = 10) ## first plot 

taxa_tbl <- taxa_table %>% 
  as_tibble()  #tidy taxonomic infos

net_work_used <- V(net_work)$name %>% 
  as_tibble() %>% 
  mutate(species = value) %>% 
  select(species) #extract species represented in the network

v_attr <- sp_ratio %>% ### we create a table of attributes for nodes (vertex)
  rownames_to_column( var = "species") %>% 
  as_tibble() %>% 
  pivot_longer(-species, names_to = "sample_id", values_to = "ratio" )

library(dplyr)

v_attr2 <- v_attr %>%  group_by(species) %>%
  summarise(rel_abundance = sum(ratio)) 


%>% 
  inner_join(net_work_used, by = "species") %>% 
  inner_join(taxa_tbl, by = "species") %>%  
  mutate(rel_abundance = abs(exp (rel_abundance))) # we join taxonomic infos, relative abundance and species that were only represented in the network to create the attribute table
 
