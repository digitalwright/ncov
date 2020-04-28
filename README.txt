README for ncov analysis presented in:
"SARS-CoV-2 genome evolution exposes early human adaptations"

# DEPENDENCIES
R >= 3.15
DECIPHER >= 2.14

# TO RUN
Open AnalyzeSequences_v5.R
Set working directory to path to ncov
Source the code in R

# FILES
AnalyzeSequences_v5.R = MAIN SCRIPT to perform all analyses
coordinates_v1.R = Positions of features in reference genome
gisaid_cov2020_sequences-Apr14.fasta.gz = FASTA file with all genomes except reference
map_v1.R = Function for mapping substitutions on the phylogenetic tree
NC_045512.2.fas = FASTA file with reference genome
results_v3.csv = Results of the analysis

# CHANGING THE DATASET
To change the dataset, simply change the file name of:
gisaid_cov2020_sequences-Apr14.fasta.gz
in:
AnalyzeSequences_v5.R
All analyses should (hopefully) work, but figures might need some adjustment.
