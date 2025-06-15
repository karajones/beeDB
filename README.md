## R Scripts for Nearctic Bee Database (beeDB)

> [!CAUTION]
> This repository is currently under development. This repository is not currently intended for use by anyone other than me.

The purpose of these scripts is to create a comprehensive metabarcoding reference sequence database for all known bee species in the Nearctic (i.e., most of the US and Canada) for all primer pairs that could work on bees. The end goals are to identify how many bees have reference sequences and which primer pairs are most likely to be effective in amplifying bees so that researchers can make more informed choices for metabarcoding.

Currently uploaded scripts (need to finish documentation):
1) Compile a checklist of all known bee species in the Nearctic and their possible synonyms
2) Identify sequence data available for those bee species from NCBI's nucleotide database
3) Clean up and add taxonomic info to NCBI sequence data
4) Prepare data downloaded from BOLD and merge with NCBI data
5) Use simplified in silico PCR to extract homologous seed sequences with primers
6) BLAST seed sequences against database to identify homologous sequences that lack primer regions
7) Clean up resulting database to remove duplicates and identify potential taxonomic conflicts

Scripts for the next steps are still being validated:
- Validate sequence homology using trees
- Calculate barcode gap
- Summarize results
- Shiny interface for visualizing results and accessing beeDB
