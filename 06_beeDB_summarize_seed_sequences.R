library(tidyverse)
library(conflicted) # provides error when conflicted namespaces are encountered
library(here)


# Compile seed seq info --------------------------------------------------------

primers <- read.table(here::here("metadata/primers.tsv"), sep = "\t", header = TRUE) %>%
  dplyr::rename(expected_length = length)

primer_gene <- primers %>%
  dplyr::select(primer, gene, expected_length)

# seed seq files
files <- sort(list.files(here::here("processed_sequences/seed_seqs"), pattern="01_", full.names = TRUE))

# import original sequences
fasta_path <- here::here("sequences/combined_data_for_cropping.fa")

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
  dplyr::rename("original_sequence" = "x") %>%
  dplyr::mutate(species = str_replace_all(species, "_", " "),
                subspecies = str_replace_all(subspecies, "_", " "),
                genus = str_split_i(species, " ", 1)) %>%
  distinct(original_sequence, beeid, .keep_all = TRUE)

# read in files and bind into one tibble
seed_seqs <- readr::read_csv(files, id = "primer_name", show_col_types = FALSE) %>%
  # get name of primer from filename
  dplyr::mutate(primer_name = str_split_i(primer_name, "_", 5),
                source = case_when(str_detect(accession, "-") ~ "BOLD",
                                   .default = "NCBI")) %>%
  dplyr::left_join(primers, by = join_by(primer_name == primer))  %>%
  dplyr::left_join(seqs_str_df, by = join_by(accession), relationship = "many-to-many") %>%
  dplyr::distinct(primer_name, accession, sequence, .keep_all = TRUE)


# Primer statistics ------------------------------------------------------------

summary_basic <- seed_seqs %>%
  group_by(primer_name) %>%
  mutate(n_accession = n_distinct(accession),
         n_seqs = n_distinct(sequence),
         n_families = n_distinct(family),
         n_genera = n_distinct(genus),
         n_species = n_distinct(species),
    min_length = min(length),
         max_length = max(length),
    length_diff = max_length-min_length,
         mean_length = round(mean(length), 1),
         median_length = median(length),
         is_median = case_when(length == median_length ~ 1, .default = 0),
         n_median = sum(is_median),
          percent_median = round((n_median/n_seqs)*100, 1),
         median_diff = expected_length - median_length,
         n_accession = n_distinct(accession),
         n_sequences = n_distinct(sequence),
         min_forward_errors = min(forward_primer_errors),
         max_forward_errors = max(forward_primer_errors),
         mean_forward_errors = round(mean(forward_primer_errors), 1),
         median_forward_errors = round(median(forward_primer_errors), 0),
         min_reverse_errors = min(reverse_primer_errors),
         max_reverse_errors = max(reverse_primer_errors),
         mean_reverse_errors = round(mean(reverse_primer_errors), 1),
         median_reverse_errors = round(median(reverse_primer_errors), 0),
    mean_primer_errors = round(mean(forward_primer_errors+reverse_primer_errors), 1)
  ) %>%
  distinct(primer_name, .keep_all = TRUE) %>%
  select(primer_name, gene, expected_length, n_accession:mean_primer_errors, -is_median)

write_tsv(summary_basic, here::here("results/summary_seed_seqs_data.tsv"))

# deduplicate seed sequences
seed_seqs_dedup <- seed_seqs %>%
  distinct(primer_name, sequence, .keep_all = TRUE) %>%
  dplyr::group_by(primer_name) %>%
  dplyr::group_split()

# write seed seqs to separate files
for(i in seq_along(seed_seqs_dedup)) {
  
  # pull out tibble
  seeds <- seed_seqs_dedup[[i,]]
  
  # get primer name
  primer <- seeds %>%
    distinct(primer_name) %>%
    pull()
  
  # convert to DNAStringSet
  seqs <- Biostrings::DNAStringSet(seeds$sequence)
  names(seqs) <- seeds$accession
  
  # save sequences to fasta file
  Biostrings::writeXStringSet(seqs, paste0(here::here("processed_sequences/seed_seqs/"), "01_seed_seqs_", primer,".fa"), append=FALSE, compress=FALSE, format="fasta")
  
}

# COMMAND LINE -----------------------------------------------------------------

# Extract gene sequences
#seqkit grep -nri -p "gene=COX1" refseq_bee_mitogenome_genes.fa > refseq_bee_mitogenome_genes.COX1.fa
#seqkit grep -nri -p "gene=COX1" refseq_bee_mitogenome_genes.fa > refseq_bee_mitogenome_genes.COX1.fa  
# Align seed sequence to mitogenomes using MAFFT
# for f in *.fa; do mafft --thread 8 --reorder --keeplength --addfragments $f --auto ../refseq/refseq_bee_mitogenome_genes.COX1.mafft.fa  > ${f%.*}.mafft.fa ; done
# for f in *.fa; do mafft --thread 8 --reorder --keeplength --addfragments $f --auto ../refseq/refseq_bee_mitogenomes_2025-04-01.mafft.fa  > ${f%.*}.mito.mafft.fa ; done

# Align seed sequences 
# for f in *.fa; do mafft --globalpair --thread 8 --leavegappyregion $f > ${f%.*}.mafft.fa; done




