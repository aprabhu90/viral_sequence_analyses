library(igraph)
library(tidyverse)

sp_ratio<- read.delim("../Data/virus_ref_coverage/virus_ref_tpm_table.tsv", header = T)
colnames(sp_ratio) <- gsub("allgenomes_nodup_sl_sr.mmi.", "", colnames(sp_ratio))
colnames(sp_ratio) <- gsub("_R1_001.fastq.gz_filtered.TPM", "", colnames(sp_ratio))


taxa_table <- read.delim("../Data/virus_ref_coverage/taxa_table3.txt", header = T)
metadata <- read.delim("../Data/virus_ref_coverage/env_data.txt", header = T, row.names = 1)


###################First network
sp_ratio <- sp_ratio %>% 
  as_tibble() %>%  
  filter(sp_ratio$Species %in% taxa_table$species) 


colnames(host_Pos_v1)[1] <- "Species"
sp_ratio2 = bind_rows(host_Pos_v1, sp_ratio)

sp_ratio = sp_ratio2 %>% column_to_rownames(var = "Species")

colnames(sp_ratio) = metadata$Site

sp_ratio <- sp_ratio[,c(7,8,9,16,17,18,25,26,27,34,35,36,43,44,45,52,53,54,61,62,63,67,68,69)]

min.prevalence=1
incidence=sp_ratio
incidence[incidence>0.1]=1
sp_ratio_filtered <- sp_ratio[which(rowSums(incidence)>=min.prevalence),] ### end of prevalence filtering

sp_correl <- sp_ratio_filtered %>% 
  t() %>% 
  cor(method = "spearman")  ### correlation calculations
sp_correl[abs(sp_correl)<0.61]=0 
net_work <- graph_from_adjacency_matrix(sp_correl,mode="lower",weighted=TRUE, diag=FALSE)  ## create network file
net_work <- delete.vertices(net_work,degree(net_work)==0) #remove nodes without edges

plot(net_work, vertex.label = NA, edge.width = 5, vertex.size = 10) ## first plot 

taxa_tbl <- taxa_table 
host_tax_Pos2 <- host_tax_Pos[,c(2,122,120)]
host_tax_Pos2$Genus_214[is.na(host_tax_Pos2$Genus_214)] <- "Poseidonia"
host_tax_Pos2$Order_214[is.na(host_tax_Pos2$Order_214)] <- "Poseioniales"
colnames(host_tax_Pos2) = colnames(taxa_table)
taxa_table$Group <- gsub("g__", "" ,taxa_table$Group)

taxa_tbl = bind_rows(taxa_table, host_tax_Pos2)
taxa_tbl$Group <- gsub("g__", "" ,taxa_tbl$Group)
net_work_used <- V(net_work)$name %>% 
  as_tibble() %>% 
  mutate(species = value) %>% 
  select(species) #extract species represented in the network

v_attr <- sp_ratio %>% ### we create a table of attributes for nodes (vertex)
  rownames_to_column( var = "species") %>% 
  as_tibble() %>% 
  pivot_longer(-species, names_to = "sample_id", values_to = "ratio" ) %>%  
  group_by(species) %>%
  summarise(rel_abundance = log10(sum(ratio))) %>% 
  inner_join(net_work_used, by = "species") %>% 
  inner_join(taxa_tbl, by = "species") %>%  
  mutate(rel_abundance = abs(exp(rel_abundance))) # we join taxonomic infos, relative abundance and species that were only represented in the network to create the attribute table

network_table <- igraph::as_data_frame(net_work, 'both') ##we convert the network in data frames

network_table$vertices <- network_table$vertices %>%
  as_tibble() %>% 
  inner_join(v_attr, by = c("name"="species")) 

net_work1 <- graph_from_data_frame(network_table$edges,
                                   directed = F,
                                   vertices = network_table$vertices) # we convert back data frames to a network


mom_data <- V(net_work1)$Group %>% 
  as_tibble_col(column_name = "species") #formating the class variable as factor (is needed for coloring edges)
mom_data$species <- as_factor(mom_data$species)
color_easy <-  c("#c95bb3",
                          "#50aa64",
                          "#8561d0",
                          "#8a9f3d",
                          "#7d7fc5",
                          "#c18b40",
                          "#45b1c4",
                          "#cc5643",
                          "#c45a7f",
                          "darkblue")[mom_data$species] #creating a color palette to represent each class levels

V(net_work1)$color <-  color_easy ## we have now the color attributes based on class

E(net_work1)$sign <- E(net_work1)$weight ##we create another attribute from the weight

E(net_work1)$weight <- abs(E(net_work1)$weight) ## we then use absolute value of weight because of the specific layout we will use

#png("../Figures/netnames_corr2.png", units="in", width=15, height=10, res=1000)

plot(net_work1, vertex.size= 6, #size of nodes
     edge.width=abs(E(net_work1)$weight)*1, #width of edges (edge is a segment between nodes)
     vertex.label = NA, #V(net_work1)$species, #remove microbes names
     vertex.label = NA, #remove microbes names
     edge.color=ifelse(E(net_work1)$sign > 0.7, 
                              "blue", 
                              "black"), # color edges based on weight sign
     layout=layout_with_kk(net_work1)) #specific layout

title_legend1 <- V(net_work1)$Group %>% 
  as_factor() %>% 
  levels() #extracting class levels in an object for the legend

legend(x = 1, y = 1, title_legend1,
       pch = 21, pt.bg = c("#c95bb3",
                                    "#50aa64",
                                    "#8561d0",
                                    "#8a9f3d",
                                    "#7d7fc5",
                                    "#c18b40",
                                    "#45b1c4",
                                    "#cc5643",
                                    "#c45a7f",
                                    "darkblue"), 
       pt.cex = 1.5, bty = "n", ncol = 1) #adding legend


dev.off()
