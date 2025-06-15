library(tidyverse)
library(conflicted)

# Import NCBI fasta
# as string set
ncbi <- Biostrings::readDNAStringSet(here::here("sequences/bee_seqs_NCBI_nucleotide_20250612.fa"), format="fasta")

ncbi_tibble <- ncbi %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  # split accession (first value before space) into one column and description into another column
  tidyr::separate(rowname, into = c("accession"), sep = "\\s", extra = "drop") %>%
  dplyr::rename("sequence" = "x")

# grab accessions
accession <- ncbi_tibble$accession
taxaId <- taxonomizr::accessionToTaxa(accession,"/Users/ksjones/taxonomizr_database/accessionTaxa.sql")

# get taxonomy
ncbi_clean <- taxonomizr::getTaxonomy(taxaId,"/Users/ksjones/taxonomizr_database/accessionTaxa.sql",
                                desiredTaxa = c("phylum", "class", "order", "family", "genus",
                                                "species", "subspecies")) %>%
  as.data.frame() %>%
  bind_cols(as.data.frame(accession)) %>%
  bind_cols(as.data.frame(taxaId)) %>%
  left_join(ncbi_tibble, by = join_by(accession)) %>%
  dplyr::filter(order == "Hymenoptera") %>%
  select(species:accession, sequence) %>%
  distinct()
  
# Read in species list
species_list <- read.csv(here::here("species_lists/corrected_species_list_ITIS.csv"), header = TRUE)

# Join with NCBI data
ncbi_final <- species_list %>%
  left_join(ncbi_clean, by = join_by(species_synonym == species), relationship = "many-to-many") %>%
  # remove entries with no sequence
  dplyr::filter(
    !is.na(sequence), 
    sequence != "",
    stringr::str_length(sequence) > 100
  )  %>%
  distinct(beeid, sequence, .keep_all = TRUE) 

# Export BOLD data
write_tsv(ncbi_final, here::here("sequences/NCBI_data_for_cropping.tsv"))

# Create fasta file
ncbi_fasta <- ncbi_final %>%
  # Create fasta header
  tidyr::unite(col = header, c(accession, family, species_final, subspecies, beeid), sep = ";", remove = TRUE, na.rm = FALSE) %>%
  # Replace underscores with spaces
  dplyr::mutate(header = str_replace_all(header, " ", "_"))

# convert to DNAStringSet
seqs <- Biostrings::DNAStringSet(ncbi_fasta$sequence)
names(seqs) <- ncbi_fasta$header

Biostrings::writeXStringSet(seqs, here::here("sequences/NCBI_data_for_cropping.fa"), append=FALSE, compress=FALSE, format="fasta")

