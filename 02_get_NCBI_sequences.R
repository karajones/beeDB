library(tidyverse)
library(conflicted)

# Read in species list
species_list <- read.csv(here::here("species_lists/corrected_species_list_ITIS.csv"), header = TRUE)


# Get NCBI accessions for species ----------------------------------------------

# Create list of all names
species_names <- c(species_list$species_final, species_list$species_synonym) %>% unique()

# Get taxonomic IDs
NCBI_taxids <- taxonomizr::getId(species_names, sqlFile = "/Users/ksjones/taxonomizr_database/accessionTaxa.sql", onlyScientific = TRUE)

# Get subspecies information
descendents <- taxonomizr::getDescendants(na.omit(NCBI_taxids), sqlFile = "/Users/ksjones/taxonomizr_database/accessionTaxa.sql", desiredTaxa = "subspecies")
descendents_taxids <- taxonomizr::getId(descendents, sqlFile = "/Users/ksjones/taxonomizr_database/accessionTaxa.sql", onlyScientific = TRUE)

# Combine subspecies names, taxids, and parent species names
NCBI_descendents <- tibble::tibble(subspecies = descendents, NCBI_taxids = descendents_taxids) %>%
  dplyr::mutate(parent = paste(stringr::str_split_i(descendents, " ", 1), stringr::str_split_i(descendents, " ", 2), sep = " "))

accessions <- taxonomizr::getAccessions(na.omit(c(NCBI_taxids, descendents_taxids)), sqlFile = "/Users/ksjones/taxonomizr_database/accessionTaxa.sql", version = "version")
write_tsv(accessions, here::here("sequences/all_accessions_NCBI.tsv"))

# Filter accessions (see: https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/)
accessions_filtered <- accessions %>% 
  # keep accessions coded for nucleotides
  dplyr::filter(str_detect(accession, "^[A-Za-z]{1}\\d{4,}") | str_detect(accession, "^[A-Za-z]{2}\\d{4,}") | str_detect(accession, "^[A-Za-z]{2}_\\d{6,}")) %>%
  # remove accessions coded for WGS/scaffold
  dplyr::filter(!str_detect(accession, "^AN|BA|CH|CM|DF|DG|DS|EM|EN|EP|EQ|GG|GL|JH|KB|KD|KE|KI|KK|KL|KN|KQ|KV|KZ|LD|ML|MU|PS|FA")) %>%
  # remove third party annotated (TPA) and transcriptome (TSA) sequences
  dplyr::filter(!str_detect(accession, "^BK|BN|BR|BL|GJ|HT|HU|GK|EZ|FX|HP|JI|JL|JO|JP|JR|JT|JU|JV|JW|KA|LA|LE|LH|LI|LJ]")) %>%
  # remove refseq sequences (they're duplicates on sequences already in the database)
  dplyr::filter(!str_detect(accession, "_"))
  
write.table(accessions_filtered$accession, here::here("sequences/filtered_accessions_NCBI.txt"), row.names = FALSE, quote = FALSE, col.names = FALSE)

# command line:
#efetch -input filtered_accessions_NCBI.txt -db nucleotide -format fasta >> bee_seqs_NCBI_nucleotide_20250612.fa


# not used ---------------------------------------------------------------------

data_summary <- rentrez::entrez_summary(db = "nucleotide", id = accessions$accession, retmode = 10)

# combine species with taxids
NCBI_species <- tibble::tibble(species_final = species_names, NCBI_taxids) %>%
  dplyr::right_join(species_list, by = join_by(species_final)) %>%
  dplyr::left_join(NCBI_descendents, by = join_by(species_final), relationship = "many-to-many") 
  dplyr::select()
  
  rename("Species_ITIS" = all_names_deduped) %>%
  full_join(species_list) %>%
  # remove any species where there with no hits in NCBI AND ITIS
  dplyr::filter(!(is.na(taxid_NCBI) & is.na(taxid_ITIS)))

# Check names against NCBI records
NCBI <- taxize::gna_verifier(species_names, data_sources = 4) %>%
  dplyr::filter(matchType == "Exact")

GBIF <- taxize::gna_verifier(species_list_clean$species_final, data_sources = 11)

