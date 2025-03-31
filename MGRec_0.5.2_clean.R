#MGRec0.5.1_clean
#Besides the recruitment plotting functions, this script includes functions to calculate breadth of coverage for each species,
#which adds another option to filter interesting species and facilitates deciding what species are likely present.
#The input files in this script have been filtered using the "random best" approach, where a random hit among the best ones is chosen, 
#if there are several "equally best" hits.

setwd("/home/urahara/Documents/Projekti/RLadies/MGRec/")
library(scales)
library(tidyverse)
library(grid)
library(gridExtra)

###### FUNCTIONS - DON'T TOUCH ######
union.all <- function(x,y){
  union(x,y)
}

cumsum0 <- function(lengths){
  cumlengths <- cumsum(as.numeric(lengths))
  cumsum0 <- cumlengths - lengths
  return(cumsum0)
}

#The following function takes the list of contigs (accession, length, genus, species) used in the blast database
#and adds numbers used to help create the X coordinates plotting the species
List_sp_Offset <- function(file){
  A <- read_tsv(file, col_names = c("Accession", "Length", "Genus", "Species")) %>%
    arrange(Genus, Species) %>% 
    #This is very necessary, because tapply is not intrinsically connected to the tibble 
    #and apparently R follows a different alphabetic order than Python (Which is first? Bacillus soli or Bacillus solimangrovi?).
    #Due to the ordering difference, the Offset_sp vector is crap unless you order the tibble to R order.
    mutate(Offset_sp = unlist(tapply(.$Length, .$Species, cumsum0), use.names = F))
  return(A)
}

#The following function takes the output of List_sp_Offset and is used to calculate offsets for genus plots
List_ge_Offset <- function(df){
  B <- df %>%
    arrange(Genus, Species) %>%
    group_by(Species, Genus) %>% 
    summarise(Length = sum(Length)) %>% 
    ungroup() %>%
    mutate(Offset_ge_basis = unlist(tapply(.$Length, .$Genus, cumsum0), use.names=F),
           Cumsum = unlist(tapply(.$Length, .$Genus, function(x) cumsum(as.numeric(x))), use.names=F))
  return(B)
}

#The following function reads the filtered blast output, uses the contig/species lists with genus/species offsets,
#and creates the X coordinates for genus-level and species-level figures
processBLASTbest <- function(file, List_sp_Offset, List_ge_Offset, qcovlimit){
  BL <- read_tsv(file, col_names = c("query", "Accession", "pid", "len", "mismatch",
                                     "gapopen", "qstart", "qend", "sstart", "send",
                                     "evalue", "bitscore", "qlen", "slen")) %>%
    left_join(List_sp_Offset, by = "Accession") %>%
    left_join(List_ge_Offset, by = "Species") %>%
    mutate(X_sp = (sstart + send)/2 + Offset_sp,
           X_ge = (sstart + send)/2 + Offset_ge_basis + Offset_sp,
           qcovmod = (qend - qstart)/qlen*100) %>%
    select(Query = query, Accession = Accession, Percentage_id = pid, Length = len, Mismatch = mismatch, GapOpen = gapopen,
           Query_start = qstart, Query_end = qend, Subject_start = sstart, Subject_end = send, Evalue = evalue, Bitscore = bitscore, Query_length = qlen, Subject_length = slen, Query_coverage = qcovmod,
           Genus = Genus.x, Species = Species, Offset_species = Offset_sp, Offset_genus_basis = Offset_ge_basis, X_species = X_sp, X_genus = X_ge) %>%
    filter(Query_coverage > qcovlimit)
  return(BL)
}

#The same function, but retaining fewer columns to save memory
processBLASTbest_memsave <- function(file, List_sp_Offset, List_ge_Offset, qcovlimit){
  BL <- read_tsv(file, col_names = c("query", "Accession", "pid", "len", "mismatch",
                                     "gapopen", "qstart", "qend", "sstart", "send",
                                     "evalue", "bitscore", "qlen", "slen")) %>%
    left_join(List_sp_Offset, by = "Accession") %>%
    left_join(List_ge_Offset, by = "Species") %>%
    mutate(X_sp = (sstart + send)/2 + Offset_sp,
           X_ge = (sstart + send)/2 + Offset_ge_basis + Offset_sp,
           qcovmod = (qend - qstart)/qlen*100) %>%
    select(Query = query, Accession = Accession, Percentage_id = pid, len, Mismatch = mismatch, GapOpen = gapopen,
           Query_start = qstart, Query_end = qend, Subject_start = sstart, Subject_end = send, Evalue = evalue, Bitscore = bitscore, 
           Query_length = qlen, Subject_length = slen, Query_coverage = qcovmod,
           Genus = Genus.x, Species = Species, Offset_species = Offset_sp, Offset_genus_basis = Offset_ge_basis, X_species = X_sp, X_genus = X_ge, Genome_length = Length.y) %>%
    select(-len, -Mismatch, -GapOpen, -Query_start, -Query_end, -Offset_species, -Offset_genus_basis, -Evalue, -Bitscore) %>%
    filter(Query_coverage > qcovlimit)
  return(BL)
}

#The following two functions take the processed BLAST tibble, a genus/species name, and makes a genus- or species-level plot

RecruitmentPlot_Genus <- function(df, sample, genus, Sp_gen_Genome,
                                  binwidth.x = tail(Sp_gen_Genome$Cumsum[Sp_gen_Genome$Genus == genus], 1)/10000, binwidth.y = 0.25,
                                  limits = c(0, tail(Sp_gen_Genome$Cumsum[Sp_gen_Genome$Genus == genus], 1)),
                                  low_col = "lightblue", mid_col = "lightblue", high_col = "black", full_y_axis = F, species_number = T, guide = T, title = T) {
  df <- filter(df, Genus == genus)
  Sp_gen_Genome <- filter(Sp_gen_Genome, Genus == genus)
  theme_base <- theme(axis.ticks.x = element_blank(),
                      axis.title.x = element_blank(),
                      axis.text.x = if(species_number == T) element_text(colour = "black") else element_blank(),
                      axis.line.y = element_line(colour = "black"),
                      axis.ticks.y = if(full_y_axis == T) element_line(colour = "black") else element_blank(),
                      axis.title.y = if(full_y_axis == T) element_text(colour = "black") else element_blank(),
                      axis.text.y = if(full_y_axis == T) element_text(colour = "black") else element_blank(),
                      panel.background = element_blank(),
                      panel.grid.major.x = element_blank(),
                      panel.grid.minor.x = element_line(colour = "gray50"),
                      panel.grid.major.y = element_blank(),
                      panel.grid.minor.y = element_blank(),
                      plot.title = if(title == T) element_text(hjust = 0.5) else element_blank())
  
  PLOT <- ggplot(df, aes(x = X_genus, y = Percentage_id)) +
    geom_bin2d(binwidth = c(binwidth.x, binwidth.y), aes(fill = ..count../max(..count..))) +
    ggtitle(substitute(a~", "~italic(b), list(a = sample, b = genus))) +
    theme_base +
    scale_y_continuous(labels = comma, limits = c(60, 100.5), expand = c(0, 0)) +
    scale_x_continuous(name = NULL, breaks = c((Sp_gen_Genome$Offset_ge_basis + Sp_gen_Genome$Cumsum)/2),
                       labels = as.character(row.names(Sp_gen_Genome)),
                       minor_breaks = c(Sp_gen_Genome$Offset_ge_basis, tail(Sp_gen_Genome$Cumsum)),
                       limits = limits,
                       expand = c(0,0))
  
  if (guide == T) {
    PLOT <- PLOT +
      scale_fill_gradientn(name = "Normalised count", colors = c(low_col, mid_col, high_col))
  } else {
    PLOT <- PLOT +
      scale_fill_gradientn(colors = c(low_col, mid_col, high_col), guide = F)
  }
  
  if (full_y_axis == T) {
    PLOT <- PLOT +
      scale_y_continuous("Identity [%]", labels = comma, limits = c(60, 100.5), expand = c(0,0))
  }
  
  return(PLOT)
}

RecruitmentPlot_Species <- function(df, sample, species, title, Sp_gen_Genome,
                                    binwidth.x = c(Sp_gen_Genome$Length[Sp_gen_Genome$Species == species]/1000), binwidth.y = 0.25,
                                    limits = c(0, Sp_gen_Genome$Length[Sp_gen_Genome$Species == species]*1.001),
                                    low_col = "lightblue", mid_col = "lightblue", high_col = "black", full_y_axis = F, guide = T) {
  df <- filter(df, Species == species)
  theme_base <- theme(axis.ticks.x = element_blank(),
                      axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                      axis.line.y = element_line(colour = "black"),
                      axis.ticks.y = if(full_y_axis == T) element_line(colour = "black") else element_blank(),
                      axis.title.y = if(full_y_axis == T) element_text(colour = "black") else element_blank(),
                      axis.text.y = if(full_y_axis == T) element_text(colour = "black") else element_blank(),
                      panel.background = element_blank(),
                      panel.grid.major.x = element_blank(),
                      panel.grid.minor.x = element_line(colour = "gray50"),
                      panel.grid.major.y = element_blank(),
                      panel.grid.minor.y = element_blank(),
                      plot.title = element_text(hjust = 0.5))
  
  PLOT<- ggplot(df, aes(x = X_species, y = Percentage_id))+
    geom_bin2d(binwidth = c(binwidth.x, binwidth.y), aes(fill = ..count../max(..count..))) +
    ggtitle(substitute(italic(b), list(b=title))) +
    theme_base +
    scale_y_continuous(labels = comma, limits = c(60, 100.5), expand = c(0, 0)) +
    scale_x_continuous(name = NULL, 
                       breaks = NULL,
                       labels = NULL,
                       minor_breaks = NULL,
                       limits = limits,
                       expand = c(0, 0))
  
  if (guide == T) {
    PLOT <- PLOT +
      scale_fill_gradientn(name = "Normalised count", colors = c(low_col, mid_col, high_col))
  }
  else {
    PLOT <- PLOT +
      scale_fill_gradientn(colors = c(low_col, mid_col, high_col), guide = F)
  }
  
  if (full_y_axis == T) {
    PLOT <- PLOT +
      scale_y_continuous("Identity [%]", labels = comma, limits = c(60, 100.5), expand = c(0,0))
  }
  
  return(PLOT)
}

RecruitmentPlot_Contig <- function(df, contigAcc, title, Acc_len_gen_sp,
                                    binwidth.x = 500, binwidth.y = 0.25,
                                    limits = c(0, Acc_len_gen_sp %>% filter(Accession == contigAcc) %>% pull(Length)),
                                    low_col = "lightblue", mid_col = "lightblue", high_col = "black", full_y_axis = F, guide = T) {
  df <- filter(df, Accession == contigAcc)
  theme_base <- theme(axis.ticks.x = element_blank(),
                      axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                      axis.line.y = element_line(colour="black"),
                      axis.ticks.y = if(full_y_axis==T) element_line(colour="black") else element_blank(),
                      axis.title.y = if(full_y_axis==T) element_text(colour="black") else element_blank(),
                      axis.text.y = if(full_y_axis==T) element_text(colour="black") else element_blank(),
                      panel.background = element_blank(),
                      panel.grid.major.x = element_blank(),
                      panel.grid.minor.x = element_line(colour="gray50"),
                      panel.grid.major.y = element_blank(),
                      panel.grid.minor.y = element_blank(),
                      plot.title = element_text(hjust = 0.5))
  
  PLOT<- ggplot(df, aes(x = (Subject_start + Subject_end)/2, y = Percentage_id)) +
    geom_bin2d(binwidth = c(binwidth.x, binwidth.y), aes(fill = ..count../max(..count..))) +
    ggtitle(substitute(italic(b), list(b = title))) +
    theme_base +
    scale_y_continuous(labels = comma, limits = c(60, 100.5), expand = c(0, 0)) +
    scale_x_continuous(name = NULL, 
                       breaks = NULL,
                       labels = NULL,
                       minor_breaks = NULL,
                       limits = limits,
                       expand = c(0, 0))
  
  if (guide == T) {
    PLOT <- PLOT +
      scale_fill_gradientn(name = "Normalised count", colors = c(low_col, mid_col, high_col))
  }
  else {
    PLOT <- PLOT +
      scale_fill_gradientn(colors = c(low_col, mid_col, high_col), guide = F)
  }
  
  if (full_y_axis == T) {
    PLOT <- PLOT +
      scale_y_continuous("Identity [%]", labels = comma, limits = c(60, 100.5), expand = c(0, 0))
  }
  
  return(PLOT)
}


sp_histogram <- function(df, sample, species, minReads){
  d <- df[df$Species == species & df$Percentage_id >= 90,]
  if (nrow(d) > minReads) {
    ggplot(d, aes(Percentage_id))+
      geom_histogram(binwidth = 0.25)+
      ggtitle(substitute(italic(b), list(b = species)))+
      scale_x_reverse(limits = c(100.25, 89.75), name = "Identity [%]", expand = c(0, 0.05))+
      scale_y_continuous(name = "Count", expand = c(0.01, 0))+
      theme(
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5))
  }}

sp_gen_table <- function(df, taxlevel, readnumber){
  df <- as.data.frame(table(df[[taxlevel]])/readnumber*100)
  colnames(df) <- c(taxlevel, "Percentage")
  df <- df[order(df$Percentage, decreasing = T),]
  return(df)
}

sp_gen_table_filtered2 <- function(df, taxlevel, readnumber){
  taxlevel <- enquo(taxlevel)
  df <- df %>% 
    group_by_at(vars(!!taxlevel)) %>% 
    summarise(Percentage = length(unique(Query))/readnumber*100) %>%
    ungroup() %>%
    arrange(desc(Percentage))
  return(df)
}

merge.all.species <- function(x,y){
  z <- merge(x, y, by = "Species", all = T)
  return(z)
}
merge.all.genus <- function(x, y){
  z <- merge(x, y, by = "Genus", all = T)
  return(z)
}

X_sp_static_bins_RSD <- function(vector_of_X_sp, num_bins = 10, max_length){
  reads_per_bin <- hist(vector_of_X_sp, breaks = seq(0, max_length, max_length/num_bins), plot = F)$counts
  return(sd(reads_per_bin)/mean(reads_per_bin))
}

Species_RSD_static_bins_ANI2 <- function(BLAST_filtered, num_bins = 10){ #If the BLAST_filtered includes the genome lengths of the reference genomes
  Sp_RSD_ANI <- BLAST_filtered %>% 
    group_by(Species) %>% 
    summarise(RSD = X_sp_static_bins_RSD(X_species, num_bins, unique(Length)),
              ANI = mean(Percentage_id),
              Reads = n()) %>% 
    ungroup() %>% 
    select(Species, RSD, ANI, Reads)
}

get_coverage_for_species_MGRec <- function(df, species, acc_len_gen_sp, pidmin, qcovmin){
  df <- df %>%
    filter(Species == species, Percentage_id >= pidmin, Query_coverage > qcovmin) %>%
    left_join(acc_len_gen_sp, by="Accession") %>%
    select(-Length, -Genus.y, -Species.y, Genus = Genus.x, Species = Species.x) %>%
    rowwise()%>%
    mutate(osstart = min(Subject_start, Subject_end) + Offset_sp, osend = max(Subject_start, Subject_end) + Offset_sp)
  hist_data_table <- df %>%
    rowwise() %>% transmute(base = list(seq(osstart, osend, 1))) %>% ungroup() %>%
    unnest(base) %>%
    group_by(base) %>% summarize(N = n()) %>% ungroup()
  return(hist_data_table)
}

#The barebone version of the function is useful to get coverage regardless of Percentage_id and Query coverage and 
#is intended to only be of use in the create_breadth_table_complete
get_coverage_for_species_MGRec_essential <- function(df, species, acc_len_gen_sp){
  df <- df %>%
    filter(Species == species) %>%
    left_join(acc_len_gen_sp,by="Accession") %>%
    select(-Length,-Genus.y,-Species.y,Genus=Genus.x,Species=Species.x) %>%
    rowwise()%>%
    mutate(osstart = min(Subject_start, Subject_end) + Offset_sp, osend = max(Subject_start, Subject_end) + Offset_sp)
  hist_data_table <- df %>%
    rowwise() %>% transmute(base = list(seq(osstart,osend,1))) %>% ungroup() %>%
    unnest(base) %>%
    group_by(base) %>% summarize(N=n()) %>% ungroup()
  return(hist_data_table)
}

get_breadth_for_species <- function(coverage_data_species_MGRec, genome_length, min_coverage){
  df <- coverage_data_species_MGRec %>%
    filter(N >= min_coverage)
  breadth <- nrow(df) / genome_length * 100
  return(breadth)
}

get_breadth_for_species_from_vector <- function(coverage_data_species_MGRec, genome_length, min_coverage){
  df <- coverage_data_species_MGRec %>%
    filter(N >= min_coverage)
  breadth <- nrow(df) / genome_length * 100
  return(breadth)
}

create_breadth_table_complete <- function(df, Acc_len_gen_sp, Sp_gen_genome, pidmin, qcovmin){ 
  df_filtered <- filter(df, Percentage_id >= pidmin, Query_coverage >= qcovmin)
  Accession <- c()
  Average_coverage <- numeric()
  Breadth <- numeric()
  for (i in unique(df_filtered$Species)){
    Genome_length <- filter(Sp_gen_Genome, Species == i) %>% pull(Length)
    coverage_table <- get_coverage_for_species_MGRec_essential(df_filtered, i, Acc_len_gen_sp)
    Accession <- c(Accession, i)
    Average_coverage <- c(Average_coverage, sum(coverage_table$N)/Genome_length)
    Breadth <- c(Breadth, get_breadth_for_species(coverage_table, Genome_length, 1))
  }
  df_breadth <- tibble(Species = Accession, Average_coverage = Average_coverage, Breadth = Breadth) %>%
    left_join(Sp_gen_Genome, by = "Species") %>%
    select(-Length, -Offset_ge_basis, -Cumsum) %>%
    arrange(desc(Breadth))
}

#The following function works if the BLAST_filtered includes the genome lengths of the reference genomes.
#On a sample of 100000 first reads of CRF1D7_MGRec, it produces identical results to the previous version, but at > 100x speed
create_breadth_table_complete2 <- function(df, Acc_len_gen_sp, pidmin = 95, qcovmin = 95, minreads = 2){ #If the BLAST_filtered includes the genome lengths of the reference genomes
  
  breadths <- filter(df, Percentage_id >= pidmin, Query_coverage >= qcovmin) %>% 
    left_join((Acc_len_gen_sp %>% select(Accession, Offset_sp)), by = "Accession") %>% 
    group_by(Query) %>% 
    mutate(osstart = min(Subject_start, Subject_end) + Offset_sp, osend = max(Subject_start, Subject_end) + Offset_sp) %>% 
    ungroup() %>% 
    group_by(Species) %>% 
    filter(n() >= minreads) %>% #Species with only one read would cause a problem
    arrange(osend, osstart) %>% #This line and the next one are there to avoid a small bug that occurs when two equal osend values that are also in the cumsum set are separated by a larger value
    filter(!duplicated(osend)) %>% 
    arrange(osstart) %>% 
    filter(osend %in% unique(cummax(osend))) %>%  #This row prevents premature ending of the interval due to short reads that are contained within other reads, even if there are more reads like that in a row
    mutate(Eshifted = c(0, osend[1:length(osend)-1]),
           max_osend = max(osend)) %>%
    filter(((Eshifted < osstart) & (Eshifted < osend))) %>% #The first part of the fix for the problem of the last interval
    mutate(Eint = c(tail(Eshifted, n()-1), unique(max_osend)),
           Interval_length = Eint - osstart + 1) %>%
    summarise(Covered_bases = sum(Interval_length),
              Breadth = sum(Interval_length)/unique(Genome_length) * 100) %>% 
    ungroup()
  
  coverages <- filter(df, Percentage_id >= pidmin, Query_coverage >= qcovmin) %>%  
    group_by(Species) %>% 
    mutate(QuasiAlignmentLength = abs(Subject_end - Subject_start) + 1) %>% 
    summarise(TotalBases = sum(QuasiAlignmentLength),
              Average_coverage = TotalBases/unique(Genome_length)) %>% 
    ungroup() %>% 
    select(Species, TotalBases, Average_coverage)
  
  breadths <- breadths %>% 
    left_join(coverages, by = "Species") %>% 
    group_by(Species) %>% 
    mutate(Average_coverage_on_covered_area = TotalBases/Covered_bases) %>% 
    ungroup()
  
  return(breadths)
}


### UTILITY TABLES #######

#Genera used: 
#The tabular output (with qlen and slen) was filtered
#Only hits retained have:
#minimum query length: 50
#minimum percentage identity: 60

#List of Contig names, lengths, genus, species:
Lactobacillus_reclassification <- read_tsv("Lactobacillus_reclassification_corr3.txt", 
                                           col_names = c("Species", "Species_new", "Genus_new"))

Acc_len_gen_sp <- List_sp_Offset("MGRec_CR3_contig_len_gen_sp.txt") %>%
  left_join(Lactobacillus_reclassification, by = "Species") %>% 
  mutate(Species_new = coalesce(Species_new, Species),
         Genus_new = coalesce(Genus_new, Genus)) %>% 
  select(Accession, Length, Genus = Genus_new, Species = Species_new, Offset_sp)

Sp_gen_Genome <- List_ge_Offset(Acc_len_gen_sp)
Sp_gen_Genome_order <- Sp_gen_Genome %>%
  group_by(Genus) %>%
  mutate(Sp_order = seq(1, n())) %>%
  ungroup()


#### INPUT DATA PROCESSING ####
#Example:
CRF2D20_MGRec <- processBLASTbest_memsave("BLAST_CR3_MGRec5/CRF2D20_MGRec_CR3_filtered_50_60_randombest.blastn", Acc_len_gen_sp, Sp_gen_Genome, 60)

CRF2D20_MGRec <- CRF2D20_MGRec %>% filter(Genus != "Theobroma")

##### RSD_ANI #####
#Example:
CRF2D20_RSD_ANI<- Species_RSD_static_bins_ANI2(CRF2D20_MGRec %>% left_join(Sp_gen_Genome %>% select(Species, Length), by = "Species")) %>% mutate(Sample="CRF2D20")

CRALL_RSD_ANI <- read_tsv("CRALL_MGRec5_RSD_ANI_filtered_Apr2020.txt",col_names = T) %>%
  left_join(Sp_gen_Genome_order, by = "Species") %>%
  select(-Length, -Offset_ge_basis, -Cumsum)

CRALL_Genus_Reads <- CRALL_RSD_ANI %>%
  group_by(Genus, Sample) %>%
  summarise(Reads = sum(Reads)) %>%
  ungroup()

CRALL_Species_Reads <- CRALL_RSD_ANI %>%
  group_by(Species, Sample) %>%
  summarise(Genus = unique(Genus), Reads = sum(Reads)) %>%
  ungroup()

CRALL_RSD08_1000 <- CRALL_RSD_ANI %>% filter(RSD <= 0.8, Reads >= 1000)

MaxReads_relevant_species_noANI <- CRALL_RSD_ANI %>% 
  group_by(Species) %>% 
  filter(RSD <= 0.5, Reads >= 1000) %>% 
  filter(Reads == max(Reads)) %>% 
  ungroup() %>% 
  arrange(Species)

relReadsOver90_MaxReads_in_CRF2D20 <- CRF2D20_MGRec %>%
  filter(Species %in% (MaxReads_relevant_species_noANI %>% filter(Sample == "CRF2D20") %>% pull(Species))) %>% 
  group_by(Species) %>% 
  summarise(relReads_over_pid90 = sum(Percentage_id > 90)/n()) %>% 
  ungroup()
  
### CALCULATE BREADTH AND COVERAGE ####
# Example:
# CRF2D20_Breadths_pidmin95_qcovmin99 <- create_breadth_table_complete2(CRF2D20_MGRec, Acc_len_gen_sp, 95, 99)
# write_tsv(CRF2D20_Breadths_pidmin95_qcovmin99,"CRF2D20_Breadths_pidmin95_qcovmin99_Apr2020.txt", col_names = T)

CRF2D20_Breadths_pidmin95_qcovmin99 <- read_tsv("CRF2D20_Breadths_pidmin95_qcovmin99_Apr2020.txt") %>% mutate(Sample = "CRF2D20")

CR_Breadths_pidmin95_qcovmin99 <- read_tsv("CR_Breadths_pidmin95_qcovmin99_Apr2020.txt")

#I'll try to plot them and see the distribution of breadths
ggplot(CR_Breadths_pidmin95_qcovmin99, aes(x = Sample, y = Breadth)) +
  geom_jitter()

#Histogram per sample
ggplot(CR_Breadths_pidmin95_qcovmin99 %>% filter(Breadth > 1), aes(x = Breadth)) +
  geom_histogram(binwidth = 5) +
  facet_wrap(vars(Sample))

#Relation between RSD and Breadth
CR_Breadth_RSD <- left_join(CR_Breadths_pidmin95_qcovmin99, CRALL_RSD_ANI, by = c("Species", "Sample")) %>% 
  select(Species, Genus = Genus.x, Sample, Reads, Average_coverage, Average_coverage_on_covered_area, Breadth, RSD, ANI, Sp_order)
ggplot(CR_Breadth_RSD, aes(y = Breadth, x = log(RSD, base = 10))) +
  geom_point(size = 0.1)

### WHICH SPECIES COULD BE PRESENT? ###

RSD_Species_pass <- function(RSD_df, genus, maxRSD, minANI, minReads){
  df <- RSD_df %>% 
    filter(Genus == genus, 
           RSD <= maxRSD, 
           ANI >= minANI, 
           Reads >= minReads) %>% 
    pull(Species) %>% 
    sort() %>% 
    unique()
  return(df)
}

Breadth_Species_pass <- function(Breadth_df, genus, minBreadth){
  df <- Breadth_df %>%
    filter(Genus == genus, 
           Breadth >= minBreadth) %>%
    pull(Species) %>% 
    sort() %>% 
    unique()
  return(df)
}

RSD_AND_Breadth_Species_pass <- function(RSD_df, Breadth_df, genus, maxRSD, minANI, minReads, minBreadth){
  RSD_df <- RSD_Species_pass(RSD_df, genus, maxRSD, minANI, minReads)
  Breadth_df <- Breadth_Species_pass(Breadth_df, genus, minBreadth)
  Both_pass <- intersect(RSD_df, Breadth_df)
  return(Both_pass)
}

theme_single_species <- theme(axis.title.y = element_text(size = 6),
                              axis.text.y = element_text(size = 6),
                              title = element_text(size = 10),
                              legend.key.size = unit(c(6), units = "points"),
                              legend.text = element_text(size = 6),
                              legend.title = element_text(size = 6),
                              legend.margin = margin(c(0,0,0,-10)))

RecruitmentPlot_Species_single_figure <- function(filename, width, height, units, res, df, species, title, Sp_gen_Genome){
  tiff(filename, width = width, height = height, units = units, res = res)
  print(RecruitmentPlot_Species(df, "", species, title, Sp_gen_Genome, full_y_axis = T) +
    theme_single_species)
  dev.off()
}


##### MAJOR GENERA - more than 1.0 % of all reads ####
#Example

#Limosilactobacillus (max. 27.2 % of all reads)
Visual_Limosilactobacillus <- Sp_gen_Genome_order %>% filter(Genus == "Limosilactobacillus", Sp_order %in% c(4, 15))
RecruitmentPlot_Genus(CRF2D20_MGRec, "CRFD20", "Limosilactobacillus", Sp_gen_Genome = Sp_gen_Genome)
#Visually: Limosilactobacillus fermentum, Limosilactobacillus vaginalis
RSD_AND_Breadth_Species_pass(CRALL_RSD_ANI, CR_Breadths_pidmin95_qcovmin99, "Limosilactobacillus", 0.4, 90, 1000, 66) %>% paste(collapse = ", ")
#Min. 66 %: Limosilactobacillus fermentum
RSD_AND_Breadth_Species_pass(CRALL_RSD_ANI, CR_Breadths_pidmin95_qcovmin99, "Limosilactobacillus", 0.4, 90, 1000, 10) %>% paste(collapse = ", ")
#Min 10 % Breadth: Limosilactobacillus fermentum
RSD_AND_Breadth_Species_pass(CRALL_RSD_ANI, CR_Breadths_pidmin95_qcovmin99, "Limosilactobacillus", 0.4, 90, 1000, 5) %>% paste(collapse = ", ")
#Min 5 % Breadth: Limosilactobacillus fermentum
CRALL_RSD_ANI %>% filter(Genus == "Limosilactobacillus") %>% filter(RSD == min(RSD) | ANI == max(ANI) | Reads == max(Reads))

##### Limosilactobacillus single species plots #####
#Example:
RecruitmentPlot_Species_single_figure("Limosilactobacillus_fermentum_single_sp.tiff", width = 8, height = 5, units = "cm", res = 300, CRF2D20_MGRec, "Limosilactobacillus fermentum", "Limosilactobacillus fermentum", Sp_gen_Genome)
