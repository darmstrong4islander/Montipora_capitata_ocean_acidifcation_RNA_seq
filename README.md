# Project Background and Preliminary Steps
This is my repository for the bioinformatic pipeline for the raw reads produced following the RNA-Seq analysis. These data were sent by AZENTA—GeneWiz and will be analyzed using a combination of an HPC with BASH and a local desktop in R Studio.

# About these Data
The coral species Montipora capitata is dominant among Hawaiian reefs and is present throughout the Indo-Pacific. This species belongs to the dominant coral family, Arocporidae, which is dominant globally. Corals grow by synthesizing calcium carbonate (CaCO3) through various physiological pathways that extract calcium and dissolved inorganic carbon from the seawater. Thus, corals rely on stable seawater chemistry to maintain this calcification growth process. However, ocean acidification chemistry changes following the increase of atmospheric carbon dioxide (CO2), has been shown to reduce calcification in a large portion of coral species. Montipora capitata is not a species we consider vulnerable to ocean acidification due to measured calcification responses to increased seawater pCO2, and could have physiological mechanisms to maintain this growth in acidified seawater. Therefore, we conducted an RNA sampling procedure and differential gene expression analysis on M. capitata exposed to ambient and acidifed seawater conditions to identify potential physiological pathways involved in calcification in ocean acidification.

# Meta-data
The raw reads are in the "raw_reads" file in this repo. The sample labeling scheme is as follows: 1-3 are replicates, B and D are treatments (B = acidified and D = ambient), and S is the summer season.

# Step 1 - Accessing the HPC at TAMUCC
Connect to the Admin wifi and Cisco Anyconnect VPN, then on your local device ssh into the HPC
```bash
ssh your_email@crest-login.tamucc.edu
```
# Step 2 - Copying repo
Copy this repo into your HPC environment
```bash
git@github.com:darmstrong4islander/mcap_rna_seq_2024.git
```
