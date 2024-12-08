# Project Background 
This is my repository for the bioinformatic pipeline for the raw reads produced following the RNA-Seq analysis. These data were sent by AZENTA—GeneWiz and will be analyzed using a combination of an HPC with BASH and a local desktop in R Studio.

# About these Data
The coral species Montipora capitata is dominant among Hawaiian reefs and is present throughout the Indo-Pacific. This species belongs to the dominant coral family, Arocporidae, which is dominant globally. Corals grow by synthesizing calcium carbonate (CaCO3) through various physiological pathways that extract calcium and dissolved inorganic carbon from the seawater. Thus, corals rely on stable seawater chemistry to maintain this calcification growth process. However, ocean acidification chemistry changes following the increase of atmospheric carbon dioxide (CO2), has been shown to reduce calcification in a large portion of coral species. Montipora capitata is not a species we consider vulnerable to ocean acidification due to measured calcification responses to increased seawater pCO2, and could have physiological mechanisms to maintain this growth in acidified seawater. Therefore, we conducted an RNA sampling procedure and differential gene expression analysis on M. capitata exposed to ambient and acidifed seawater conditions to identify potential physiological pathways involved in calcification in ocean acidification.

# Meta-data
The raw reads are in the "raw_reads" file in this repo. The sample labeling scheme is as follows: 1-3 are replicates, B and D are treatments (B = acidified and D = ambient), and S is the summer season.

# Accessing the HPC at TAMUCC
Connect to the Admin wifi and Cisco Anyconnect VPN, then on your local device ssh into the HPC
```bash
ssh your_email@crest-login.tamucc.edu
```
# Copying repo
Copy this repo into your HPC environment
```bash
git@github.com:darmstrong4islander/Montipora_capitata_ocean_acidifcation_RNA_seq.git
```
# Upload your data from local machine
On your local device scp raw reads from a local directory to the HPC directory where you copied this repo, you may have to use passcodes specific to you local ssh or HPC ssh keys
```bash
scp -r "./Mcap.genome_assembly.fa" your_user@crest-login.tamucc.edu:~/Montipora_capitata_ocean_acidifcation_RNA_seq/reference_genome
scp -r "./Mcap.GFFannotation.putnam.gff" your_user@crest-login.tamucc.edu:~/Montipora_capitata_ocean_acidifcation_RNA_seq/reference_genome
scp -r "./00_fastq" your_user@crest-login.tamucc.edu:~/Montipora_capitata_ocean_acidifcation_RNA_seq/00_fastq
```
# Load modules
```bash
module load fastp                
module load hisat2               
module load samtools             
```
# Verify Checksums 
Run this script, this should populate an output checksum report, and an error file where any errors that occured when running the scrip will show
```bash
sbatch checksum_file.sh
```

# RNA data QC and Alignment
Run this script which integrates tools on the HPC and in this repo, this should generate preliminary files of count matrices in which we will then merge in the following script
```bash
sbatch mcap_rna_cnt.sh
```
