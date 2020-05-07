README for ncov analysis presented in:
"SARS-CoV-2 genome evolution exposes early human adaptations"

# AUTHOR
Erik S. Wright <eswright@pitt.edu>

# DEPENDENCIES
R >= 4.0
DECIPHER >= 2.16.1

# TO RUN
Install R then DECIPHER
Open AnalyzeSequences_v6.R
Set working directory to path to ncov
Source the code in R

# FILES
AnalyzeSequences_v6.R = MAIN SCRIPT to perform all analyses
coordinates_v1.R = Positions of features in reference genome
gisaid_cov2020_sequences-May2.fasta.gz = FASTA file with all genomes
map_v1.R = Function for mapping substitutions on the phylogenetic tree
movavg_v1.R = Function for performing center-point exponential moving averaging
NC_045512.2.fas = FASTA file with reference genome
results_v5.csv = Results of the analysis
metadata_v1.tsv = Tab delimited matrix of metadata and acknowledgements

# CHANGING THE DATASET
To change the dataset, simply change the file name of:
gisaid_cov2020_sequences-May2.fasta.gz
and:
metadata_v1.tsv
in:
AnalyzeSequences_v6.R
All analyses should (hopefully) still work, but figures will need some adjustment.
