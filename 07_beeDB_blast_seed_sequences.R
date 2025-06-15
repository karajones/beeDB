library(tidyverse)
library(progress) # add progress bar to for loops
library(conflicted) # provides error when conflicted namespaces are encountered
library(here)
library(Biostrings) # DNA sequence manipulation tools
library(rBLAST) # BLAST searches without leaving R
library(conflicted)

# Path to sequences ------------------------------------------------------------

# import original sequences
path <- here::here("processed_sequences/cropped_seqs")
blast_db_path <- here::here("unprocessed_sequences/combined_data_for_cropping.fa")
primers <- read.table(here::here("metadata/primers.tsv"), sep = "\t", header = TRUE) 

#rBLAST::makeblastdb(file = blast_db_path, dbtype = "nucl")

# Import BLAST database

# as string set
seqs_str <- Biostrings::readDNAStringSet(blast_db_path, format="fasta")
seqs_str_names <- names(seqs_str) %>% 
  as.data.frame() %>%
  dplyr::rename("accession" = ".")

# as data frame
seqs_str_df <- as.data.frame(seqs_str) %>%
  dplyr::bind_cols(seqs_str_names) %>%
  dplyr::rename("sequence" = "x") %>%
  distinct(accession, .keep_all = TRUE)


# as_dna
seqs <- bioseq::as_dna(seqs_str_df$sequence)
names(seqs) <- seqs_str_df$accession

# fasta files for seed sequences
files <- sort(list.files(
  here::here("processed_sequences/seed_seqs"),
  pattern = ".fa",
  full.names = TRUE
))


# BEGIN MAIN LOOP ------------------------------------------------

for (i in seq_along(files)) {
  
  # get primer name from file name
  primer_name <- stringr::str_split_i(basename(files[i]), pattern = "_", 4) %>%
    str_remove(".fa")
  
  primer_metadata <- dplyr::filter(primers, primer == primer_name)
  gene <-  pull(primer_metadata, gene)
  expected_length <- pull(primer_metadata, length)
  
  # calculate min/max length (±5% for COX, ±15% for others)
  if (gene == "COX1") {
  min_length <- expected_length-(expected_length*0.05)
  max_length <- expected_length+(expected_length*0.05)
  } else {
    min_length <- expected_length-(expected_length*0.15)
    max_length <- expected_length+(expected_length*0.15)
  }
  
  # import sequences from file
  seed_seqs <- Biostrings::readDNAStringSet(files[i], format = "fasta")
  
  ## BLAST seeds against database -----------------------------------------------
  
  # prep for blastn
  blastn <- rBLAST::blast(db = blast_db_path, type = "blastn")
  
  # run BLAST
  blast_results <- predict(
    blastn,
    seed_seqs,
    custom_format = "qseqid sseqid sstart send length qcovs bitscore"
  )
  
  # remove duplicates
  blast_results_filtered <- blast_results %>%
    # remove qseqid
    dplyr::select(-qseqid) %>%
    # sort by longest length
    # there may be duplicate entries and this will ensure only the longest are kept
    dplyr::arrange(desc(bitscore)) %>%
    # remove duplicate entries
    dplyr::distinct(sseqid, .keep_all = TRUE)
  
  cli::cli_alert_info(
    "BLAST identified {nrow(blast_results_filtered)} candidate sequences for {primer_name}."
  )
  
  ## Add taxonomy to sequences and save -----------------------------------------
  
  # Begin loop: Check if there are results -------------------------------------------------
  
  if (nrow(blast_results_filtered > 0)) {

  taxa_accessions <- blast_results_filtered %>%
    dplyr::left_join(seqs_str_df, by = join_by(sseqid == accession)) %>%
    tidyr::separate(sseqid, into = c("accession","family","species","subspecies","beeid"), sep = ";", extra = "drop") %>%
    dplyr::mutate(species = str_replace_all(species, "_", " "),
                  subspecies = str_replace_all(subspecies, "_", " "),
                  genus = str_split_i(species, " ", 1),
                  subspecies = case_when(subspecies == "" ~ NA_character_,
                                         subspecies == "NA" ~ NA_character_,
                                         .default = subspecies)) %>%
    # check for reversed sequences
    dplyr::mutate(start = dplyr::case_when(sstart > send ~ send, .default = sstart),
                  end = dplyr::case_when(sstart > send ~ sstart, .default = send),
                  reversed = dplyr::case_when(sstart > send ~ TRUE, .default = FALSE))
    
  # write CSV of candidate sequences
  write.csv(
    taxa_accessions,
    paste0(path, "/02_", primer_name, "_blast_results_unfiltered.csv"),
    row.names = FALSE)
  
    # Crop sequences by position ---------------------------------------------------
    
    # initiate tibble for cropped sequences
    cropped_sequences <- tibble::tibble(this_accession = character(), sequence = character())
    
    # cropping progress bar
    cli::cli_progress_bar("Cropping sequences",
                          type = "iterator",
                          total = nrow(taxa_accessions))
    
    # extract sequence based on BLAST position
    for (i in 1:nrow(taxa_accessions)) {
      cli::cli_progress_update()
      
      this_accession <- taxa_accessions[i, ]$accession
      this_sstart <- taxa_accessions[i, ]$start
      this_ssend <- taxa_accessions[i, ]$end
      this_seq <- taxa_accessions[i, ]$sequence
      
      sequence <- stringr::str_sub(this_seq, start = this_sstart, end = this_ssend)
      new_seq <- cbind(this_accession, sequence)
      cropped_sequences <- rbind(cropped_sequences, new_seq)
    }
    
    cli::cli_alert_success("Completed cropping candidate sequences.")
    
    # merge cropped sequences with the rest of the data
    sequence_database <- taxa_accessions %>%
      # remove uncropped sequence
      dplyr::select(-c(sequence, sstart, send, bitscore, qcovs)) %>%
      # add cropped sequence
      dplyr::left_join(cropped_sequences, by = join_by(accession == this_accession), relationship = "many-to-many") %>%
      dplyr::distinct() %>%
      dplyr::filter(length > min_length,
                    length < max_length) 

    # pull out sequences that need to be reversed
    sequence_database_reversed <- sequence_database %>%
      dplyr::filter(reversed == TRUE)
    
    # reverse complement sequences
    reverse_seq <- bioseq::as_dna(sequence_database_reversed$sequence) %>%
      bioseq::seq_complement() %>%
      bioseq::seq_reverse() %>%
      tibble::as_tibble()
    
    sequence_database_revcomp <- cbind(sequence_database_reversed, reverse_seq) %>%
      # remove old sequence
      dplyr::select(-sequence) %>%
      # rename new sequence column
      dplyr::rename("sequence" = value)
    
    # merge sequences back with data
    sequence_database_fixed <- sequence_database %>%
      dplyr::filter(reversed == FALSE) %>%
      dplyr::bind_rows(sequence_database_revcomp) %>%
      dplyr::mutate(primer_name = primer_name, gene = gene)
    
    write.csv(
      sequence_database_fixed,
      paste0(path, "/03_", primer_name, "_blast_results_filtered.csv"),
      row.names = FALSE
    )
    
    cli::cli_alert_info(
      "There were {nrow(sequence_database_fixed)} candidate sequences for {primer_name} after filtering."
    )
    
    # Remove duplicates ------------------------------------------------------------
    
    # identify duplicate sequences with same taxonomy
    duplicate_sequences <- sequence_database_fixed %>%
      # identify all duplicates
      janitor::get_dupes(sequence, species)
    
    duplicates_to_keep <- duplicate_sequences %>%
      # group identical sequences
      dplyr::group_by(sequence) %>%
      # select one sequence to represent each group of duplicates
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup() %>%
      # get vector of accessions
      dplyr::pull(accession)
    
    duplicates_to_remove <- duplicate_sequences %>%
      # remove accessions to keep
      dplyr::filter(!accession %in% duplicates_to_keep) %>%
      dplyr::pull(accession)
    
    # remove duplicates from database
    sequence_database_dedup <- sequence_database_fixed %>%
      dplyr::filter(!accession %in% duplicates_to_remove)
    
    
  # Add sequence conflicts -----------------------------------------------------
    
  # identify duplicate sequences with different taxonomy
  sequence_conflicts <- sequence_database_dedup %>%
    janitor::get_dupes(sequence) %>%
    dplyr::group_by(sequence) %>%
    # create column with count of number of conflicts
    dplyr::mutate(species_count = dplyr::n_distinct(species, na.rm = TRUE),
                  genus_count = dplyr::n_distinct(genus, na.rm = TRUE),
                  family_count = dplyr::n_distinct(family, na.rm = TRUE)) %>%
    # identify highest taxonomic level of conflict
    dplyr::mutate(conflict_level = case_when(family_count > 1 ~ "family",
                                             genus_count > 1 ~ "genus",
                                             species_count > 1 ~ "species",
                                             .default = NA_character_)) %>%
    # combine names of conflicts into a single column
    # stri_remove_empty_na to remove empty strings and NAs
    # stri_unique to only list each conflict once
    dplyr::mutate(conflict_taxa = case_when(conflict_level == "species" ~ paste0(stringi::stri_remove_empty_na(stringi::stri_unique(species)), collapse = ", "),
                                            conflict_level == "genus" ~ paste0(stringi::stri_remove_empty_na(stringi::stri_unique(genus)), collapse = ", "),
                                            conflict_level == "family" ~ paste0(stringi::stri_remove_empty_na(stringi::stri_unique(family)), collapse = ", "),
                                            .default = NA_character_)) %>%
    dplyr::select(sequence, accession:gene, conflict_count = dupe_count, conflict_level, conflict_taxa)

  # add conflict flags to database
  sequence_database_with_conflicts <- sequence_database_dedup %>%
    dplyr::left_join(sequence_conflicts, by = join_by(accession, family, species, subspecies, beeid, length, genus, start,
                                                       end, reversed, sequence, primer_name, gene)) 
  
  write.csv(sequence_database_with_conflicts, 
              paste0(path, "/04_", primer_name, "_sequence_database.csv"),
              row.names = FALSE)
  }
  else {
  cli::cli_alert_warning("No sequences found.")
  } # End loop to check if there are results -----------------------------------
  
  
} # END MAIN LOOP -----------------------------------------------------------------




