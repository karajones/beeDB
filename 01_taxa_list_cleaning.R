library(tidyverse)

# Initial checklist cleaning ---------------------------------------------------

# Import raw species list
species_list <- read.table(here::here("species_lists/raw_species_lists.tsv"), 
                           sep = "\t", header = TRUE) %>%
  # clean up white space issues
  dplyr::mutate(species_final = str_trim(str_squish(species_final)),
                species_synonym = str_trim(str_squish(species_synonym)))

# De-duplicate species list
species_list_clean <- species_list %>%
  dplyr::distinct(species_final)


# Check names against ITIS -----------------------------------------------------
# check each species name against Global Names Verifier database
# data source: 3 = ITIS; 4 = NCBI; 11 = GBIF
ITIS <- taxize::gna_verifier(species_list_clean$species_final, data_sources = 3) %>%
  dplyr::filter(matchType == "Exact")  %>%
    dplyr::select("species_final" = submittedName, 
                  "taxid_TSN" = currentRecordId, 
                  "valid_authority_ITIS" = currentName, 
                  "valid_species_ITIS" = currentCanonicalSimple,
                  taxonomicStatus)

# Merge ITIS data with species list
species_list_ITIS <- dplyr::left_join(species_list, ITIS, by = join_by(species_final)) %>%
  # replace invalid names with valid ITIS names
  dplyr::mutate(species_final = case_when(taxonomicStatus == "Synonym" ~ valid_species_ITIS,
                                          .default = species_final),
                # create a column indicating whether name is accepted or a synonym
                taxonomic_status = case_when(species_final == species_synonym ~ "accepted",
                                             .default = "synonym")) %>%
  dplyr::select(-taxonomicStatus)

# Export correct list with valid names
write.csv(species_list_ITIS, here::here("species_lists/corrected_species_list_ITIS.csv"), row.names = FALSE)


# Summarize species list -------------------------------------------------------

# Re-import corrected list with valid names
species_list_ITIS <- read.csv(here::here("species_lists/corrected_species_list_ITIS.csv"), header = TRUE) %>%
  dplyr::group_by(species_final) %>%
  # add bee id to each species
  mutate(beeid = paste0("B",str_pad(cur_group_id(), 5, pad = "0")))

# Save corrected list with valid names and beeIDs
write.csv(species_list_ITIS, here::here("species_lists/corrected_species_list_ITIS.csv"), row.names = FALSE)

# Summarize bees by species
species_list_summary <- species_list_ITIS %>%
  dplyr::mutate(taxonomic_status = case_when(species_final == species_synonym ~ "accepted",
                                             .default = "synonym")) %>%
  dplyr::group_by(beeid) %>%
  dplyr::mutate(n_synonyms = n_distinct(species_synonym)-1) %>%
    dplyr::summarise(across(everything(), ~ toString(na.omit(unique(.x)))), .groups = "drop") %>%
  dplyr::select(beeid, family, species_final, species_synonym, n_synonyms, source, country)
  
# Save summary species list
write.csv(species_list_summary, here::here("species_lists/species_list_summary.csv"), row.names = FALSE)

species_summary <- species_list_summary %>%
  dplyr::mutate(country = case_when(str_detect(country, ",") ~ "both", 
                .default = country)) %>%
  #dplyr::group_by(country) %>%
  dplyr::group_by(family) %>%
  dplyr::summarise(n_distinct(species_final))

