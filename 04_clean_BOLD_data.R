library(tidyverse)
library(conflicted)

bold <- read.table(here::here("BOLD_data/bold_results_2025-06-10.tsv"), sep = "\t", header = TRUE, fill = TRUE, quote = "")

bold_clean <- bold %>%
  dplyr::select("accession" = processid, species, subspecies, "sequence" = nuc) %>%  
    # remove gaps from sequences and add length of sequence
    dplyr::mutate(
      sequence = stringr::str_to_upper(stringr::str_remove_all(sequence, "-")),
      sequence = stringr::str_remove(sequence, "^N+")
      ) %>%
    # remove entries with no sequence or species info
    dplyr::filter(
      !is.na(sequence), 
      sequence != "",
      !is.na(species),
      species != "",
      stringr::str_length(sequence) > 100
      )   
                
# Read in species list
species_list <- read.csv(here::here("species_lists/corrected_species_list_ITIS.csv"), header = TRUE)

# Join with BOLD data
bold_final <- species_list %>%
  left_join(bold_clean, by = join_by(species_synonym == species), relationship = "many-to-many")  %>%
  # remove entries with no sequence
  dplyr::filter(
    !is.na(sequence), 
    sequence != "",
    stringr::str_length(sequence) > 100
  )  %>%
  distinct(beeid, sequence, .keep_all = TRUE) 

# Export BOLD data
write_tsv(bold_final, here::here("sequences/BOLD_data_for_cropping.tsv"))

# Create fasta file
bold_fasta <- bold_final %>%
  # Create fasta header
  tidyr::unite(col = header, c(accession, family, species_final, subspecies, beeid), sep = ";", remove = TRUE, na.rm = FALSE) %>%
  # Replace underscores with spaces
  dplyr::mutate(header = str_replace_all(header, " ", "_"))

# convert to DNAStringSet
seqs <- Biostrings::DNAStringSet(bold_fasta$sequence)
names(seqs) <- bold_fasta$header

Biostrings::writeXStringSet(seqs, here::here("sequences/BOLD_data_for_cropping.fa"), append=FALSE, compress=FALSE, format="fasta")


