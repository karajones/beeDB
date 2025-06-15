library(tidyverse)
library(progress) # add progress bar to for loops
library(conflicted) # provides error when conflicted namespaces are encountered
library(here)
library(Biostrings) # DNA sequence manipulation tools

# Path to sequences ------------------------------------------------------------

# NCBI
path <- here::here("processed_sequences")
fasta_path <- here::here("sequences/combined_data_for_cropping.fa")
suffix <- "2025-06-12"

# import primer metadata
primers <- read.table(here::here("metadata/primers.tsv"), sep = "\t", header = TRUE) 

for (i in 1:nrow(primers)) {
  
  # Look for primers in sequences ----------------------------------------------
  ## define primer sequences ----
  primer_name <- primers[i,]$primer
  primer_gene <- primers[i,]$gene
  forward_primer <- primers[i,]$forward_sequence
  reverse_primer <- primers[i,]$reverse_sequence
  
  # get primer lengths
  forward_primer_length <- str_count(forward_primer)
  reverse_primer_length <- str_count(reverse_primer)
  expected_length <- primers[i,]$length
  
  # get reverse complement of reverse primer (3' --> 5')
  reverse_primer_reversecomp <- bioseq::as_dna(reverse_primer) %>%
    bioseq::seq_complement() %>%
    bioseq::seq_reverse()
  
  # create dictionary of disambiguated sequences for forward and reverse primers
  forward_primer_dictionary <- bioseq::seq_disambiguate_IUPAC(bioseq::as_dna(forward_primer)) %>%
    unlist() %>%
    as.data.frame() %>%
    dplyr::rename("forward_primer_disambiguated" = ".") %>%
    mutate(forward_primer = forward_primer)
  
  reverse_primer_dictionary <- bioseq::seq_disambiguate_IUPAC(bioseq::as_dna(reverse_primer_reversecomp)) %>%
    unlist() %>%
    as.data.frame() %>%
    dplyr::rename("reverse_primer_disambiguated" = ".") %>%
    mutate(reverse_primer = reverse_primer_reversecomp)
  
  # calculate min/max length (±20%)
  min_length <- expected_length-(expected_length*0.2)
  max_length <- expected_length+(expected_length*0.2)
  
  ## Import fasta ----
  # as string set
  seqs_str <- Biostrings::readDNAStringSet(fasta_path, format="fasta")
  seqs_str_names <- names(seqs_str) %>% 
    as.data.frame() %>%
    dplyr::rename("header" = ".")
  
  # as data frame
  seqs_str_df <- as.data.frame(seqs_str) %>%
    dplyr::bind_cols(seqs_str_names) %>%
    # split accession (first value before space) into one column and description into another column
    tidyr::separate(header, into = c("accession","family","species","subspecies","beeid"), sep = ";", extra = "drop") %>%
    dplyr::rename("sequence" = "x") %>%
    distinct(sequence, beeid, .keep_all = TRUE)
  
  # as_dna
  seqs <- bioseq::as_dna(seqs_str_df$sequence)
  names(seqs) <- seqs_str_df$accession
  
  
  # Extract amplicons using primer sequences ------------------------------------- 
  # allows error in primers
  # return without primer portion of sequence
  
  cli::cli_alert_info("Looking for sequences for {primer_name}.")
  
  # extract sequences using primers; allow 20% error in primers
  seed_seqs <- 
    bioseq::seq_crop_pattern(
      seqs,
      forward_primer,
      reverse_primer_reversecomp,
      max_error_in = 0.2,
      max_error_out = 0.2,
      include_patterns = TRUE
      # remove sequences that don't have matches 
    ) %>% na.omit() %>%
    # convert to tibble
    bioseq::as_tibble.bioseq_dna() 
  
  # if there were sequences extracted, then continue...
  if (nrow(seed_seqs) > 0) {
  seed_seqs_with_primers <- seed_seqs %>%
    dplyr::rename(accession = label) %>%
    # extract primer sequences
    mutate(sequence_length = str_count(sequence),
           forward_primer_match = str_sub(sequence, 1, forward_primer_length),
           reverse_primer_match = str_sub(sequence, start = sequence_length-reverse_primer_length, end = sequence_length),
           length_without_primers = sequence_length-(forward_primer_length+reverse_primer_length)) %>%
    # join by dictionary and calculate number of errors in forward primer sequence
    fuzzyjoin::stringdist_full_join(forward_primer_dictionary, by = c("forward_primer_match" = "forward_primer_disambiguated"), 
                                    max_dist = forward_primer_length,
                                    ignore_case = TRUE,
                                    distance_col = "forward_primer_errors") %>%
    # keep only the matches with the lowest number of errors
    dplyr::arrange(forward_primer_errors) %>%
    dplyr::distinct(accession, .keep_all = TRUE) %>%
    # join by dictionary and calculate number of errors in forward reverse sequence
    # the stringdist call will throw an error if there are no sequences that pass the length filter
    fuzzyjoin::stringdist_full_join(reverse_primer_dictionary, by = c("reverse_primer_match" = "reverse_primer_disambiguated"), 
                                    max_dist = reverse_primer_length,
                                    ignore_case = TRUE,
                                    distance_col = "reverse_primer_errors") %>%
    dplyr::filter(length_without_primers > min_length,
                  length_without_primers < max_length) %>%
    dplyr::arrange(reverse_primer_errors) %>%
    dplyr::distinct(accession, .keep_all = TRUE) %>%
    dplyr::mutate(sequence_without_primers = str_sub(sequence, 1+forward_primer_length, sequence_length-reverse_primer_length)) %>%
    dplyr::select(accession, sequence = sequence_without_primers, length = length_without_primers, forward_primer_match, forward_primer_errors, reverse_primer_match, reverse_primer_errors)
  
  # write CSV of seed sequences
  write.csv(seed_seqs_with_primers, paste0(path, "/seed_seqs/01_", primer_name, "_seed_sequences_", suffix, ".csv"), row.names = FALSE)
  
  # print number of sequences
  cli::cli_alert_info("There were {nrow(seed_seqs_with_primers)} initial seed sequences for {primer_name}.")
  
  
  } else 
    
    # ...else print there were no sequences
    cli::cli_alert_info("{primer_name} had no sequences.")

}

  