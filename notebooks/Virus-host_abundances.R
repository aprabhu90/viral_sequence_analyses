###Read virus and host abundance

virus_tpm <- read.delim("../viruses/virus_count_table.tsv", header = T, row.names = 1)
host_tpm <- read.delim("../host/host_tpm_table.tsv", header = T)

mean_virus_tpm <- data.frame(rowMeans(virus_tpm))
mean_virus_tpm$Virus = rownames(mean_virus_tpm)
#mean_host_tpm <- data.frame(rowMeans(host_tpm))


join_Avir <- mean_virus_tpm %>% filter( mean_virus_tpm$Virus %in% Archaea_HQ$Virus)
Avir <- left_join(Archaea_HQ, join_Avir, by = c("Virus"))
Archaea_HQ$Host.genome <- gsub("GB_", "", Archaea_HQ$Host.genome)
Archaea_HQ$Host.genome <- gsub("RS_", "", Archaea_HQ$Host.genome)

br_headers <- read.delim("../host/BR_headers_filename.tzt", header = F)
br_headers$V2 <- br_headers$V1
br_headers <- separate(br_headers, V2, into = c("V2", "contig"), sep = ":" )
br_headers$V2 <- gsub(".fa", "", br_headers$V2)
br_headers$V1 <- gsub(".fa:>", "_", br_headers$V1)

gtdb_headers <- read.delim("../host/GTDB_headers.txt", header = F)
gtdb_headers$V2 <- gtdb_headers$V1
gtdb_headers$V2 <- substr(gtdb_headers$V2, 1, 16)
gtdb_headers$V2 <- gsub(">", "", gtdb_headers$V2)
gtdb_headers$V1 <- gsub(">", "", gtdb_headers$V1)
gtdb_headers$V1 <- sub(" .*", "", gtdb_headers$V1)

headers <- rbind(br_headers[, c(1:2)], gtdb_headers)

############Read host

host_tpm2 <- left_join(host_tpm, headers, by = c("Contig" = "V1"))

colnames(host_tpm2)

hosttpm_filtered <- host_tpm2 #[rowSums(host_tpm2[,c(-1,-71)])>0,]
colnames(hosttpm_filtered) <- gsub("BR_GTDB_1027x_sr.mmi.", "", colnames(hosttpm_filtered))
colnames(hosttpm_filtered) <- gsub("_R1_001.fastq.gz_filtered.TPM", "", colnames(hosttpm_filtered))
host_tpm3 <-hosttpm_filtered[, -1] %>% group_by(V2) %>% summarise_all(sum)

gtdbtk_tax <-  read.delim("../Data/ar53_taxonomy_r214.tsv", header = F)
gtdbtk_tax$V1 <- gsub("RS_", "", gtdbtk_tax$V1)
gtdbtk_tax$V1 <- gsub("GB_", "", gtdbtk_tax$V1)
br_tax <-  read.delim("../Data/brgenomes.txt", header = T)
br_tax2 <- br_tax[,c(1,7)]
colnames(gtdbtk_tax) = colnames(br_tax2)
gtbtk_br_tax <- bind_rows(br_tax, gtdbtk_tax)


host_tax <- left_join(host_tpm3, gtbtk_br_tax, by = c("V2" = "Bin.Id")) %>% relocate(GTDB_214_taxonomy)

host_tax_Pos <- host_tax %>% filter(str_detect(host_tax$GTDB_214_taxonomy, "Poseidoniales"))

host_Pos_v1 <- host_tax_Pos[,c(2:71)]
dim(host_Pos_v1)
################


Ac4 <- ggarrange(Ac0, Ac1, Ac2, nrow = 1)
Ac4
ggsave("../Figures/fig2b_v2.png", plot = Ac4, dpi = 600, units = c("in"), width = 20, height = 7)

