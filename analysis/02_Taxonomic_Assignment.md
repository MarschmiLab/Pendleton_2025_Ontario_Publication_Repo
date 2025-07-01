---
title: "Taxonomic Assignment" 
author: "Augustus Pendleton"
date: "26 June, 2025"
output:
  html_document:
    code_folding: show
    highlight: default
    keep_md: yes
    theme: journal
    toc: yes
    toc_float:
      collapsed: no
      smooth_scroll: yes
      toc_depth: 3
editor_options: 
  chunk_output_type: console
---
<style>
pre code, pre, code {
  white-space: pre !important;
  overflow-x: scroll !important;
  word-break: keep-all !important;
  word-wrap: initial !important;
}
</style>





# Purpose

The purpose of this document is to take the feature count (otu) table we generated with dada2 and assign taxonomy to our ASVs. 


# Load packages, source files, and R objects from DADA2_Processing


```r
pacman::p_load(tidyverse, phyloseq, Biostrings, dada2, install = FALSE)

knitr::write_bib(file = "data/02_taxonomy_exports/packages.bib")

# Load in the colors and shapes
source("code/R/plotting_aesthetics.R")

# Load in our feature count table
load("data/01_dada2_exports/seqtab_nc_len.RData")

# And load our meta_data
load("data/01_dada2_exports/meta_track.RData")
```

# Assign Taxonomy with TaxAss


```bash

mkdir -p data/TaxAss_analysis

# First, clone Taxass into Taxass folder

cd data/TaxAss_analysis # Move into Taxass data folder

git clone git@github.com:McMahonLab/TaxAss.git # Clone Taxass Github Repo

cd TaxAss # Move into git repo

cd FreshTrain-files # Move to database zip folders

unzip FreshTrain15Jun2020silva138.zip # Unzip newest version

cd FreshTrain15Jun2020silva138 # Move into it

# Move files into TaxAss_analysis directory
mv -f silva_nr_v138_taxass.fasta ../../../

mv -f silva_nr_v138_taxass.taxonomy ../../../

mv -f FreshTrain15Jun2020silva138.fasta ../../../

mv -f FreshTrain15Jun2020silva138.taxonomy ../../../

# Move back to root of TaxAss git repo
cd ../../

# Move into tax-scripts
cd tax-scripts

# Move everyting into TaxAss_analysis directory
mv -f * ../../

# Move to TaxAss_analysis directory
cd ../../

# Remove TaxAss git repo
rm -r -f TaxAss


# Put dependencies in path
export PATH=/programs/mothur:$PATH
export PATH=/programs/ncbi-blast-2.13.0+/bin:$PATH

cat $PATH
# Run taxass with 50 threads
./RunSteps_quickie.sh ASVs FreshTrain15Jun2020silva138 silva_nr_v138_taxass 98 80 80 50

# Make directory for taxonomy exports if it doesn't exist yet
mkdir -p ../02_taxonomy_exports

# Move our final result file there
mv -f ASVs.98.80.80.taxonomy ../02_taxonomy_exports
mv -f FreshTrain15Jun2020silva138.taxonomy ../02_taxonomy_exports


# Remove everything else besides our ASVs.fasta file (so repo looks exactly like it did when we started)
ls | grep -xv "ASVs.fasta" | xargs rm

# Move back to project root
cd ../../../
```

```
## Cloning into 'TaxAss'...
## Archive:  FreshTrain15Jun2020silva138.zip
##    creating: FreshTrain15Jun2020silva138/
##   inflating: __MACOSX/._FreshTrain15Jun2020silva138  
##   inflating: FreshTrain15Jun2020silva138/Compare_FreshTrain_and_Silva138_Names.xlsx  
##   inflating: __MACOSX/FreshTrain15Jun2020silva138/._Compare_FreshTrain_and_Silva138_Names.xlsx  
##   inflating: FreshTrain15Jun2020silva138/.DS_Store  
##   inflating: __MACOSX/FreshTrain15Jun2020silva138/._.DS_Store  
##   inflating: FreshTrain15Jun2020silva138/README-138.html  
##   inflating: __MACOSX/FreshTrain15Jun2020silva138/._README-138.html  
##   inflating: FreshTrain15Jun2020silva138/silva_nr_v138_taxass.taxonomy  
##   inflating: __MACOSX/FreshTrain15Jun2020silva138/._silva_nr_v138_taxass.taxonomy  
##   inflating: FreshTrain15Jun2020silva138/README-138_cyano_edits.html  
##   inflating: __MACOSX/FreshTrain15Jun2020silva138/._README-138_cyano_edits.html  
##   inflating: FreshTrain15Jun2020silva138/FreshTrain15Jun2020silva138.taxonomy  
##   inflating: FreshTrain15Jun2020silva138/Compare_Edited_Cyano_and_Silva138_Names.xlsx  
##   inflating: __MACOSX/FreshTrain15Jun2020silva138/._Compare_Edited_Cyano_and_Silva138_Names.xlsx  
##   inflating: FreshTrain15Jun2020silva138/silva_nr_v138_taxass.fasta  
##   inflating: FreshTrain15Jun2020silva138/silva_nr_v138_taxass_cyano_edits.taxonomy  
##   inflating: __MACOSX/FreshTrain15Jun2020silva138/._silva_nr_v138_taxass_cyano_edits.taxonomy  
##   inflating: FreshTrain15Jun2020silva138/FreshTrain15Jun2020silva138.fasta  
## cat: '/programs/ncbi-blast-2.13.0+/bin:/programs/mothur:/workdir/arp277/2022_Muskegon_GrowthRate/miniconda3/condabin:/home/arp277/.local/bin:/home/arp277/bin:/usr/share/Modules/bin:/programs/docker/bin:/programs/geos-3.11.0/bin:/programs/gdal-3.5.2/bin:/programs/geos-3.11.0/bin:/programs/gdal-3.5.2/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/usr/lib/rstudio-server/bin/quarto/bin:/usr/lib/rstudio-server/bin/postback:/programs/bin:/programs/bin/mummer:/programs/bin/util:/programs/bin/bowtie:/programs/bin/bwa:/programs/bin/cufflinks:/programs/bin/tophat:/programs/bin/fastx:/programs/bin/blast:/programs/bin/blat:/programs/bin/perlscripts:/programs/bin/labutils:/programs/bin/gnuplot:/programs/bin/seqclean:/programs/bin/blast+:/programs/bin/sra:/programs/bin/bedtools/bin:/programs/bin/plink:/programs/bin/fastqc:/programs/bin/bowtie2:/programs/bin/clustalw:/programs/bin/rsem:/programs/bin/vcftools:/programs/RepeatMasker:/programs/bin/exonerate/bin:/programs/augustus/bin:/programs/bin/structure:/programs/bin/irods:/programs/bin/bedops:/programs/iv/x86_64/bin:/usr/lib64/openmpi/bin:/programs/texlive/bin/x86_64-linux:/programs/R-4.0.5-r9/bin:/programs/samtools-1.20/bin:/programs/bcftools-1.20/bin:/programs/htslib-1.20/bin:/usr/Arcconf': No such file or directory
## Running TaxAss as quickly as possible- just generate the final taxonomy table!
## 
## otu filename: ASVs.fasta
## custom database filenames: FreshTrain15Jun2020silva138.fasta and FreshTrain15Jun2020silva138.taxonomy
## general database filenames: silva_nr_v138_taxass.fasta and silva_nr_v138_taxass.taxonomy
## 
## Creating the final taxonomy table with pident: 98
## 
## Classification bootstrap confidence cutoffs are 80 % for the custom classification and 80 % for the general classification.
## Using 50 processors.
## 
## 
## Building a new DB, current time: 06/26/2025 13:05:10
## New DB name:   /local/workdir/arp277/Pendleton_2025_Ontario_Publication_Repo/data/TaxAss_analysis/FreshTrain15Jun2020silva138.db
## New DB title:  FreshTrain15Jun2020silva138.fasta
## Sequence type: Nucleotide
## Keep MBits: T
## Maximum file size: 3000000000B
## Adding sequences from FASTA; added 1275 sequences in 0.119254 seconds.
## 
## 
## WARNING: ignoring environment value of R_HOME
## 
## From March 1979
##    by Tomas Tranströmer
##    translated to english by Robin Fulton
## 
## Weary of all who come with words, words but no language
## I make my way to the snow-covered island.
## The untamed has no words.
## The unwritten pages spread out on every side!
## I come upon the tracks of deer in the snow.
## Language but no words.
## 
## WARNING: ignoring environment value of R_HOME
## 
## Added 353 sequences matching the custom database with >=  98 % sequence identity to ids.above.98 
## Assign taxonomy to these sequences using the Small Custom database
## 
## WARNING: ignoring environment value of R_HOME
## 
## Added 5540 sequences matching the custom database with <  98 % sequence identity to ids.below.98 
## Assign taxonomy to these sequences using the Large General database
## 
## Linux version
## 
## Using Boost,HDF5,GSL
## mothur v.1.48.0
## Last updated: 5/20/22
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## 
## University of Michigan
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## For questions and analysis support, please visit our forum at https://forum.mothur.org
## 
## Type 'quit()' to exit program
## 
## [NOTE]: Setting random seed to 19760620.
## 
## Script Mode
## 
## 
## 
## mothur > classify.seqs(fasta=ASVs.above.98.fasta, template=FreshTrain15Jun2020silva138.fasta,  taxonomy=FreshTrain15Jun2020silva138.taxonomy, method=wang, probs=T, processors=50, cutoff=0)
## 
## Using 50 processors.
## Generating search database...    DONE.
## It took 1 seconds generate search database.
## 
## Reading in the FreshTrain15Jun2020silva138.taxonomy taxonomy...	DONE.
## Calculating template taxonomy tree...     DONE.
## Calculating template probabilities...     DONE.
## It took 1 seconds get probabilities.
## Classifying sequences from ASVs.above.98.fasta ...
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 8
## 7
## 7
## 8
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 7
## 8
## 7
## 7
## 7
## 7
## 7
## 
## It took 0 secs to classify 353 sequences.
## 
## 
## It took 0 secs to create the summary file for 353 sequences.
## 
## 
## Output File Names: 
## ASVs.above.98.FreshTrain15Jun2020silva138.wang.taxonomy
## ASVs.above.98.FreshTrain15Jun2020silva138.wang.tax.summary
## 
## 
## mothur > quit()
## 
## 
## It took 1 seconds to run 2 commands from your script.
## 
## Logfile : mothur.1750957522.logfile
## 
## Linux version
## 
## Using Boost,HDF5,GSL
## mothur v.1.48.0
## Last updated: 5/20/22
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## 
## University of Michigan
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## For questions and analysis support, please visit our forum at https://forum.mothur.org
## 
## Type 'quit()' to exit program
## 
## [NOTE]: Setting random seed to 19760620.
## 
## Script Mode
## 
## 
## 
## mothur > classify.seqs(fasta=ASVs.below.98.fasta, template=silva_nr_v138_taxass.fasta, taxonomy=silva_nr_v138_taxass.taxonomy, method=wang, probs=T, processors=50, cutoff=0)
## 
## Using 50 processors.
## Generating search database...    DONE.
## It took 47 seconds generate search database.
## 
## Reading in the silva_nr_v138_taxass.taxonomy taxonomy...	DONE.
## Calculating template taxonomy tree...     DONE.
## Calculating template probabilities...     DONE.
## It took 123 seconds get probabilities.
## Classifying sequences from ASVs.below.98.fasta ...
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 100
## 168
## 168
## 168
## 168
## 168
## 168
## 168
## 168
## 168
## 169
## 168
## 168
## 168
## 169
## 168
## 169
## 168
## 168
## 168
## 169
## 169
## 169
## 169
## 169
## 168
## 168
## 169
## 169
## 168
## 169
## 168
## 169
## 169
## 169
## 169
## 168
## 169
## 168
## 168
## 170
## 168
## 168
## 168
## 168
## 168
## 169
## 168
## 168
## 168
## 168
## 
## [WARNING]: mothur reversed some your sequences for a better classification.  If you would like to take a closer look, please check ASVs.below.98.silva_nr_v138_taxass.wang.flip.accnos for the list of the sequences.
## 
## It took 4 secs to classify 8419 sequences.
## 
## 
## It took 0 secs to create the summary file for 8419 sequences.
## 
## 
## Output File Names: 
## ASVs.below.98.silva_nr_v138_taxass.wang.taxonomy
## ASVs.below.98.silva_nr_v138_taxass.wang.tax.summary
## ASVs.below.98.silva_nr_v138_taxass.wang.flip.accnos
## 
## 
## mothur > quit()
## 
## 
## It took 127 seconds to run 2 commands from your script.
## 
## Logfile : mothur.1750957523.logfile
## 
## 
## 
## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Detected 1 [WARNING] messages, please review.
## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## WARNING: ignoring environment value of R_HOME
## 
## And the Days Are Not Full Enough
## by Ezra Pound
## 
## And the days are not full enough
## And the nights are not full enough
## And life slips by like a field mouse
## 	Not shaking the grass.
## 
## Made file:  ASVs.98.80.80.taxonomy 
## All done. If this finished without errors, RunStep_16.sh to delete intermediate files and organize the rest into convenient folders. 
##  
```

# Assign taxonomy using Dada/Silva

Here, we will use the silva database version 138. You'll need to download these databases into the 02_taxonomy_exports folder. 


```r
# The next line took 2 mins to run
taxa <- assignTaxonomy(seqtab_nc_len, "data/02_taxonomy_exports/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

# the next line took 3 minutes 
taxa <- addSpecies(taxa, "data/02_taxonomy_exports/silva_species_assignment_v138.1.fa.gz")

# Inspect the taxonomy 
taxa_print <- taxa # Removing sequence rownames for display only
rownames(taxa_print) <- NULL
head(taxa_print, 20)
```

```
##       Kingdom    Phylum              Class                 Order                 Family                 Genus                       Species
##  [1,] "Bacteria" "Proteobacteria"    "Alphaproteobacteria" "SAR11 clade"         "Clade III"            NA                          NA     
##  [2,] "Bacteria" "Actinobacteriota"  "Actinobacteria"      "Frankiales"          "Sporichthyaceae"      "hgcI clade"                NA     
##  [3,] "Bacteria" "Chloroflexi"       "Anaerolineae"        "Anaerolineales"      "Anaerolineaceae"      NA                          NA     
##  [4,] "Bacteria" "Cyanobacteria"     "Cyanobacteriia"      "Chloroplast"         NA                     NA                          NA     
##  [5,] "Bacteria" "Proteobacteria"    "Gammaproteobacteria" "Burkholderiales"     "Methylophilaceae"     "Candidatus Methylopumilus" NA     
##  [6,] "Bacteria" "Actinobacteriota"  "Acidimicrobiia"      "Microtrichales"      "Ilumatobacteraceae"   "CL500-29 marine group"     NA     
##  [7,] "Bacteria" "Cyanobacteria"     "Cyanobacteriia"      "Synechococcales"     "Cyanobiaceae"         "Cyanobium PCC-6307"        NA     
##  [8,] "Bacteria" "Verrucomicrobiota" "Verrucomicrobiae"    "Methylacidiphilales" "Methylacidiphilaceae" NA                          NA     
##  [9,] "Bacteria" "Verrucomicrobiota" "Verrucomicrobiae"    NA                    NA                     NA                          NA     
## [10,] "Bacteria" "Actinobacteriota"  "Actinobacteria"      "Frankiales"          "Sporichthyaceae"      "hgcI clade"                NA     
## [11,] "Bacteria" "Actinobacteriota"  "Actinobacteria"      "Frankiales"          "Sporichthyaceae"      "hgcI clade"                NA     
## [12,] "Bacteria" "Actinobacteriota"  "Actinobacteria"      "Frankiales"          "Sporichthyaceae"      "Candidatus Planktophila"   NA     
## [13,] "Bacteria" "Chloroflexi"       "SL56 marine group"   NA                    NA                     NA                          NA     
## [14,] "Bacteria" "Verrucomicrobiota" "Verrucomicrobiae"    NA                    NA                     NA                          NA     
## [15,] "Bacteria" "Actinobacteriota"  "Acidimicrobiia"      "Microtrichales"      "Ilumatobacteraceae"   "CL500-29 marine group"     NA     
## [16,] "Bacteria" "Cyanobacteria"     "Cyanobacteriia"      "Chloroplast"         NA                     NA                          NA     
## [17,] "Bacteria" "Verrucomicrobiota" "Verrucomicrobiae"    "Pedosphaerales"      "Pedosphaeraceae"      "SH3-11"                    NA     
## [18,] "Bacteria" "Cyanobacteria"     "Cyanobacteriia"      "Chloroplast"         NA                     NA                          NA     
## [19,] "Bacteria" "Bacteroidota"      "Bacteroidia"         "Flavobacteriales"    "Crocinitomicaceae"    "Fluviicola"                NA     
## [20,] "Bacteria" "Verrucomicrobiota" "Verrucomicrobiae"    NA                    NA                     NA                          NA
```


## 2. Taxonomy Table 


```r
##### Prepare tax table 
# Add the ASV sequences from the rownames to a column 
new_tax_tab <- taxa %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASVseqs") 
head(new_tax_tab)
```

```
##                                                                                                                                                                                                                                                        ASVseqs
## 1 ACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCGTAAAGCGCGCGTAGGCCGTTGAGTTAGTTAATTGTGAAATCCCAAAGCTTAACTTTGGAACTGCAATTAAAACTGCTCGACTAGAGTTTGATAGAGGAAAGCGGAATACATAGTGTAGAGGTGAAATTCGTAGATATTATGTAGAACACCAGTTGCGAAGGCGGCTTTCTGGATCAACACTGACGCTGAGGCGCGAAAGTATGGGTAGCAAAGAGG
## 2 ACATAGGGTGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCTCGTAGGTCGTTTGTTACGTCGGATGTGAAAACCTGAGGCTCAACCTCAGGCCTGCATTCGATACGGGCAAACTAGAGTTTGGTAGGGGAGACTGGAATTCCTGGTGTAGCGGTGGAATGCGCAGATATCAGGAGGAACACCAATGGCGAAGGCAGGTCTCTGGGCCAATACTGACACTGAGGAGCGAAAGTCTGGGGAGCGAACAGG
## 3 ACGTAGGATCCAAGCGTTATCCGGAATTACTGGGCGTAAAGGGCGTGTAGGAGGTTGGGCAAGTCGGCCATGAAAGCTCCCGGCTCAACTGGGAGAGGCTGGTCGATACTGCCTGGCTAGAGGGCAAGAGAGGGAGGTGGAATTCCCGGTGTAGTGGTGAAATGCGTAGATATCGGGAGGAACACCAGTGGCGAAGGCGGCCTCCTGGCTTGTACCTGACTCTGAAACGCGAAAGCATGGGTAGCAAACAGG
## 4 ACAGAGGATGCAAGCGTTATCCGGAATCACTGGGCATAAAGCGTCTGTAGGTGGTTTGGTAAGTCTGCTGTTAAAGACTGGGGCTCAACCCCAGAAAAGCAGTGGAAACTGCCAGACTTGAGTGTGGTAGAGGTAGAGGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAAGAACACCAATGGCGAAGGCACTCTACTGGACCATAACTGACACTGAGAGACGACAGCTAGGGGAGCAAATGGG
## 5 ACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTTTGTAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCGTTCGAAACTGCAAGGCTAGAGTGTGTCAGAGGGGGGTAGAATTCCACGTGTAGCAGTGAAATGCGTAGAGATGTGGAGGAATACCAATGGCGAAGGCAGCCCCCTGGGATAACACTGACGCTCATGCACGAAAGCGTGGGGAGCAAACAGG
## 6 ACATAGGCTTCAAGCGTTGTCCGGATTTATTGGGCGTAAAGAGTTCGTAGGCGGTCGAGTAAGTCGGGTGTGAAAATTCTGGGCTCAACCCAGAAACGCCACCCGATACTGCTTAACTTGAGTTCGATAGGGGAGTGGGGAATTCCTAGTGTAGCGGTGAAATGCGCAGATATTAGGAGGAACACCGGTGGCGAAGGCGCCACTCTGGATCGATACTGACGCTGAGGAACGAAAGCATGGGTAGCAAACAGG
##    Kingdom           Phylum               Class           Order             Family                     Genus Species
## 1 Bacteria   Proteobacteria Alphaproteobacteria     SAR11 clade          Clade III                      <NA>    <NA>
## 2 Bacteria Actinobacteriota      Actinobacteria      Frankiales    Sporichthyaceae                hgcI clade    <NA>
## 3 Bacteria      Chloroflexi        Anaerolineae  Anaerolineales    Anaerolineaceae                      <NA>    <NA>
## 4 Bacteria    Cyanobacteria      Cyanobacteriia     Chloroplast               <NA>                      <NA>    <NA>
## 5 Bacteria   Proteobacteria Gammaproteobacteria Burkholderiales   Methylophilaceae Candidatus Methylopumilus    <NA>
## 6 Bacteria Actinobacteriota      Acidimicrobiia  Microtrichales Ilumatobacteraceae     CL500-29 marine group    <NA>
```

```r
# intution check 
stopifnot(new_tax_tab$ASVseqs == colnames(seqtab_nc_len))

asv_headers <- vector(dim(seqtab_nc_len)[2], mode = "character")

# loop through vector and fill it in with ASV names 

for (i in 1:dim(seqtab_nc_len)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}

# intitution check
head(asv_headers, 20)
```

```
##  [1] ">ASV_1"  ">ASV_2"  ">ASV_3"  ">ASV_4"  ">ASV_5"  ">ASV_6"  ">ASV_7"  ">ASV_8"  ">ASV_9"  ">ASV_10" ">ASV_11" ">ASV_12" ">ASV_13" ">ASV_14" ">ASV_15" ">ASV_16" ">ASV_17" ">ASV_18" ">ASV_19" ">ASV_20"
```

```r
##### Rename ASVs in table then write out our ASV fasta file! 
#View(seqtab_nochim)
asv_tab <- t(seqtab_nc_len)
#View(asv_tab)

## Rename our asvs! 
row.names(asv_tab) <- sub(">", "", asv_headers)

# Now let's add the ASV names 
rownames(new_tax_tab) <- rownames(asv_tab)
#View(new_tax_tab)

### Final prep of tax table. Add new column with ASV names 
asv_tax <- 
  new_tax_tab %>%
  # add rownames from count table for phyloseq handoff
  mutate(ASV = rownames(asv_tab)) %>%
  # Resort the columns with select
  dplyr::select(Kingdom, Phylum, Class, Order, Family, Genus, Species, ASV, ASVseqs)

# Intution check
stopifnot(asv_tax$ASV == rownames(asv_tax), rownames(asv_tax) == rownames(asv_tab))
```



```r
taxass_taxa <- read_csv("data/02_taxonomy_exports/ASVs.98.80.80.taxonomy") %>%
  mutate(across(everything(), ~ifelse(str_detect(.x, "unclassified."), NA, .x)), # Convert unclassified to NAs
         across(everything(), ~str_remove_all(.x, "\\(\\d+\\)"))) %>% # Remove the boostrap values after names
  rename_with(.cols = everything(), ~paste0(.x, "_taxass")) # Change column names

cat("Using newest silva, we failed to classify",sum(is.na(asv_tax)), "slots")
```

```
## Using newest silva, we failed to classify 20588 slots
```

```r
cat("Using taxass, we failed to classify",sum(is.na(taxass_taxa)), "slots")
```

```
## Using taxass, we failed to classify 27325 slots
```

```r
# Here, we figure out which ASVs were assigned from within Taxass
taxass_lineages <- read_tsv("data/02_taxonomy_exports/FreshTrain15Jun2020silva138.taxonomy", col_names = FALSE) %>% # Read in the tax file from the TaxAss custom database
  separate(X2, sep = ";", into = c("kingdom_taxass","phylum_taxass","class_taxass","order_taxass","lineage_taxass","clade_taxass","tribe_taxass","x")) %>% #Do some name cleaning
  mutate(across(everything(), ~ifelse(str_detect(.x, "unnamed."), NA, .x))) %>% # Instead of unclassified, convert to NA (to match our later output)
  inner_join(taxass_taxa) %>% # Innerjoin with our taxonomy file - only ASVs that were classified with the Taxass custom database (assumed to be well curated) stay
  pull(seqID_taxass) # Find those ASV ids.
taxa_joined <- inner_join(asv_tax, taxass_taxa, by = c("ASV" = "seqID_taxass")) %>%
  mutate(from_taxass = ASV %in% taxass_lineages)

# Next, I want to see if when we assign Taxa from TaxAss, do we get better lineage resolution?
taxa_joined %>%
  mutate(Order_Winner = case_when(is.na(Order)&!is.na(order_taxass)~"Taxass",
                                  !is.na(Order)&is.na(order_taxass)~"Silva"),
         Family_Winner = case_when(is.na(Family)&!is.na(lineage_taxass)~"Taxass",
                                  !is.na(Family)&is.na(lineage_taxass)~"Silva"),
         Genus_Winner = case_when( is.na(Genus)&!is.na(clade_taxass)~"Taxass",
                                  !is.na(Genus)&is.na(clade_taxass)~"Silva"),
         Species_Winner = case_when(is.na(Species)&!is.na(tribe_taxass)~"Taxass",
                                   !is.na(Species)&is.na(tribe_taxass)~"Silva")) %>% # Detect whether TaxAss or Silva assigned an ASV to a level that the other couldn't
  group_by(from_taxass) %>%
  summarize(Order_Win_Perc =   sum(str_detect(Order_Winner, "Taxass"), na.rm = T) / sum(!is.na(Order_Winner)),
            Family_Win_Perc = sum(str_detect(Family_Winner, "Taxass"), na.rm = T) / sum(!is.na(Family_Winner)),
            Genus_Win_Perc =  sum(str_detect(Genus_Winner, "Taxass"), na.rm = T) / sum(!is.na(Genus_Winner)),
            Species_Win_Perc = sum(str_detect(Species_Winner, "Taxass"), na.rm = T) / sum(!is.na(Species_Winner))) # Calculate what percentage of time TaxAss "Won" when is used its custom database
```

```
## # A tibble: 2 × 5
##   from_taxass Order_Win_Perc Family_Win_Perc Genus_Win_Perc Species_Win_Perc
##   <lgl>                <dbl>           <dbl>          <dbl>            <dbl>
## 1 FALSE              0.00179         0.00306         0.0429            0    
## 2 TRUE             NaN               1               0.25              0.914
```

Generally, if taxa were assigned out of the Freshtrain custom database, we have better assignment to the Order, Family, and Species level. It is slightly worse than Silva at the genus level. Some things, however, we just can't control.


```r
taxa_joined %>% 
  pivot_longer(cols = !c(ASV, ASVseqs, from_taxass), names_to = "Level", values_to = "Classification") %>%
  group_by(Level) %>%
  summarize(nas = sum(is.na(Classification))) %>%
  mutate(Classifier = ifelse(str_detect(Level, "taxass"), "TaxAss", "Silva"),
         Level = str_remove(Level, "_taxass"),
         Level = str_to_title(Level),
         Level = case_match(Level,
                            "Lineage" ~ "Family",
                          "Clade"~ "Genus",
                          "Tribe" ~ "Species",
                          .default = Level),
         Level = factor(Level, levels = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))) %>%
  ggplot(aes(x = Level, y = nas, color = Classifier, group = Classifier)) + 
  geom_point() + 
  geom_line() +
  theme_classic() + 
  labs(x = "Phylogenetic Level", y = "Number of Unclassified")
```

<img src="../figures/02_Taxonomic_Assignment/taxass vs silva-1.png" style="display: block; margin: auto;" />

```r
taxa_decided <- taxa_joined %>%
  transmute(ASV = ASV,
            Kingdom = ifelse(from_taxass, kingdom_taxass, Kingdom),
            Phylum = ifelse(from_taxass, phylum_taxass, Phylum),
            Class = ifelse(from_taxass, class_taxass, Class),
            Order = ifelse(from_taxass, order_taxass, Order),
            Family = ifelse(from_taxass, lineage_taxass, Family),
            Genus = ifelse(from_taxass, clade_taxass, Genus),
            Species = ifelse(from_taxass, tribe_taxass, Species),
            ASVseqs = ASVseqs)

cat("Prioritizing TaxAss and then the newest silva, we failed to classify",sum(is.na(taxa_decided)), "slots")
```

```
## Prioritizing TaxAss and then the newest silva, we failed to classify 20511 slots
```

```r
rownames(taxa_decided) <- taxa_decided$ASV
```




```r
# Write the table 
write.table(taxa_joined, "data/02_taxonomy_exports/ASV_both_taxonomy.tsv", sep = "\t", quote = FALSE, col.names = NA)
write.table(taxa_decided, "data/02_taxonomy_exports/ASV_final_taxonomy.tsv", sep = "\t", quote = FALSE, col.names = NA)
```


## Handoff to phyloseq 

```r
# Add names to rownames for phyloseq happiness 
rownames(meta_track) <-meta_track$DNA_ID

raw_physeq <- phyloseq(otu_table(asv_tab, taxa_are_rows = TRUE),
                         sample_data(meta_track),
                         tax_table(as.matrix(taxa_decided)))

raw_physeq
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 8772 taxa and 151 samples ]
## sample_data() Sample Data:       [ 151 samples by 47 sample variables ]
## tax_table()   Taxonomy Table:    [ 8772 taxa by 9 taxonomic ranks ]
```

```r
save(raw_physeq, file = "data/02_taxonomy_exports/Raw_Physeq.RData")
```

# Session Information 

```r
# Reproducibility
devtools::session_info()
```

```
## ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
##  setting  value
##  version  R version 4.3.3 (2024-02-29)
##  os       Rocky Linux 9.5 (Blue Onyx)
##  system   x86_64, linux-gnu
##  ui       X11
##  language (EN)
##  collate  en_US.UTF-8
##  ctype    en_US.UTF-8
##  tz       America/New_York
##  date     2025-06-26
##  pandoc   3.1.1 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/ (via rmarkdown)
## 
## ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
##  ! package              * version   date (UTC) lib source
##  P abind                  1.4-5     2016-07-21 [?] CRAN (R 4.3.2)
##  P ade4                   1.7-22    2023-02-06 [?] CRAN (R 4.3.2)
##  P ape                    5.7-1     2023-03-13 [?] CRAN (R 4.3.2)
##  P Biobase                2.62.0    2023-10-24 [?] Bioconductor
##  P BiocGenerics         * 0.48.1    2023-11-01 [?] Bioconductor
##  P BiocManager            1.30.22   2023-08-08 [?] CRAN (R 4.3.2)
##  P BiocParallel           1.36.0    2023-10-24 [?] Bioconductor
##  P biomformat             1.30.0    2023-10-24 [?] Bioconductor
##  P Biostrings           * 2.70.1    2023-10-25 [?] Bioconductor
##  P bit                    4.0.5     2022-11-15 [?] CRAN (R 4.3.2)
##  P bit64                  4.0.5     2020-08-30 [?] CRAN (R 4.3.2)
##  P bitops                 1.0-7     2021-04-24 [?] CRAN (R 4.3.2)
##  P bslib                  0.5.1     2023-08-11 [?] CRAN (R 4.3.2)
##  P cachem                 1.0.8     2023-05-01 [?] CRAN (R 4.3.2)
##  P callr                  3.7.3     2022-11-02 [?] CRAN (R 4.3.2)
##  P cli                    3.6.1     2023-03-23 [?] CRAN (R 4.3.2)
##  P cluster                2.1.4     2022-08-22 [?] CRAN (R 4.3.2)
##  P codetools              0.2-19    2023-02-01 [?] CRAN (R 4.3.3)
##  P colorspace             2.1-0     2023-01-23 [?] CRAN (R 4.3.2)
##  P crayon                 1.5.2     2022-09-29 [?] CRAN (R 4.3.2)
##  P dada2                * 1.30.0    2023-10-24 [?] Bioconductor
##  P data.table             1.15.2    2024-02-29 [?] CRAN (R 4.3.2)
##  P DelayedArray           0.28.0    2023-10-24 [?] Bioconductor
##  P deldir                 2.0-4     2024-02-28 [?] CRAN (R 4.3.2)
##  P devtools               2.4.4     2022-07-20 [?] CRAN (R 4.2.1)
##  P digest                 0.6.33    2023-07-07 [?] CRAN (R 4.3.2)
##  P dplyr                * 1.1.3     2023-09-03 [?] CRAN (R 4.3.2)
##  P ellipsis               0.3.2     2021-04-29 [?] CRAN (R 4.3.2)
##  P evaluate               0.23      2023-11-01 [?] CRAN (R 4.3.2)
##  P fansi                  1.0.5     2023-10-08 [?] CRAN (R 4.3.2)
##  P farver                 2.1.1     2022-07-06 [?] CRAN (R 4.3.2)
##  P fastmap                1.1.1     2023-02-24 [?] CRAN (R 4.3.2)
##  P forcats              * 1.0.0     2023-01-29 [?] CRAN (R 4.3.2)
##  P foreach                1.5.2     2022-02-02 [?] CRAN (R 4.3.2)
##  P fs                     1.6.3     2023-07-20 [?] CRAN (R 4.3.2)
##  P generics               0.1.3     2022-07-05 [?] CRAN (R 4.3.2)
##  P GenomeInfoDb         * 1.38.0    2023-10-24 [?] Bioconductor
##  P GenomeInfoDbData       1.2.11    2023-11-07 [?] Bioconductor
##  P GenomicAlignments      1.38.2    2024-01-16 [?] Bioconduc~
##  P GenomicRanges          1.54.1    2023-10-29 [?] Bioconductor
##  P ggplot2              * 3.5.0     2024-02-23 [?] CRAN (R 4.3.2)
##  P glue                   1.6.2     2022-02-24 [?] CRAN (R 4.3.2)
##  P gtable                 0.3.4     2023-08-21 [?] CRAN (R 4.3.2)
##  P highr                  0.10      2022-12-22 [?] CRAN (R 4.3.2)
##  P hms                    1.1.3     2023-03-21 [?] CRAN (R 4.3.2)
##  P htmltools              0.5.7     2023-11-03 [?] CRAN (R 4.3.2)
##  P htmlwidgets            1.6.2     2023-03-17 [?] CRAN (R 4.3.2)
##  P httpuv                 1.6.12    2023-10-23 [?] CRAN (R 4.3.2)
##  P hwriter                1.3.2.1   2022-04-08 [?] CRAN (R 4.3.2)
##  P igraph                 1.5.1     2023-08-10 [?] CRAN (R 4.3.2)
##  P interp                 1.1-6     2024-01-26 [?] CRAN (R 4.3.2)
##  P IRanges              * 2.36.0    2023-10-24 [?] Bioconductor
##  P iterators              1.0.14    2022-02-05 [?] CRAN (R 4.3.2)
##  P jpeg                   0.1-10    2022-11-29 [?] CRAN (R 4.3.2)
##  P jquerylib              0.1.4     2021-04-26 [?] CRAN (R 4.3.2)
##  P jsonlite               1.8.7     2023-06-29 [?] CRAN (R 4.3.2)
##  P knitr                  1.45      2023-10-30 [?] CRAN (R 4.3.2)
##  P labeling               0.4.3     2023-08-29 [?] CRAN (R 4.3.2)
##  P later                  1.3.1     2023-05-02 [?] CRAN (R 4.3.2)
##  P lattice                0.21-9    2023-10-01 [?] CRAN (R 4.3.2)
##  P latticeExtra           0.6-30    2022-07-04 [?] CRAN (R 4.3.2)
##  P lifecycle              1.0.3     2022-10-07 [?] CRAN (R 4.3.2)
##  P lubridate            * 1.9.3     2023-09-27 [?] CRAN (R 4.3.2)
##  P magrittr               2.0.3     2022-03-30 [?] CRAN (R 4.3.2)
##  P MASS                   7.3-60    2023-05-04 [?] CRAN (R 4.3.2)
##  P Matrix                 1.6-1.1   2023-09-18 [?] CRAN (R 4.3.2)
##  P MatrixGenerics         1.14.0    2023-10-24 [?] Bioconductor
##  P matrixStats            1.2.0     2023-12-11 [?] CRAN (R 4.3.2)
##  P memoise                2.0.1     2021-11-26 [?] CRAN (R 4.3.2)
##  P mgcv                   1.9-0     2023-07-11 [?] CRAN (R 4.3.2)
##  P mime                   0.12      2021-09-28 [?] CRAN (R 4.3.2)
##  P miniUI                 0.1.1.1   2018-05-18 [?] CRAN (R 4.3.2)
##  P multtest               2.58.0    2023-10-24 [?] Bioconductor
##  P munsell                0.5.0     2018-06-12 [?] CRAN (R 4.3.2)
##  P NatParksPalettes     * 0.2.0     2022-10-09 [?] CRAN (R 4.3.2)
##  P nlme                   3.1-163   2023-08-09 [?] CRAN (R 4.3.2)
##  P pacman                 0.5.1     2019-03-11 [?] CRAN (R 4.3.2)
##  P permute                0.9-7     2022-01-27 [?] CRAN (R 4.3.2)
##  P phyloseq             * 1.46.0    2023-10-24 [?] Bioconductor
##  P pillar                 1.9.0     2023-03-22 [?] CRAN (R 4.3.2)
##  P pkgbuild               1.4.2     2023-06-26 [?] CRAN (R 4.3.2)
##  P pkgconfig              2.0.3     2019-09-22 [?] CRAN (R 4.3.2)
##  P pkgload                1.3.3     2023-09-22 [?] CRAN (R 4.3.2)
##  P plyr                   1.8.9     2023-10-02 [?] CRAN (R 4.3.2)
##  P png                    0.1-8     2022-11-29 [?] CRAN (R 4.3.2)
##  P prettyunits            1.2.0     2023-09-24 [?] CRAN (R 4.3.2)
##  P processx               3.8.2     2023-06-30 [?] CRAN (R 4.3.2)
##  P profvis                0.3.8     2023-05-02 [?] CRAN (R 4.3.2)
##  P promises               1.2.1     2023-08-10 [?] CRAN (R 4.3.2)
##  P ps                     1.7.5     2023-04-18 [?] CRAN (R 4.3.2)
##  P purrr                * 1.0.2     2023-08-10 [?] CRAN (R 4.3.2)
##  P R6                     2.5.1     2021-08-19 [?] CRAN (R 4.3.2)
##  P RColorBrewer           1.1-3     2022-04-03 [?] CRAN (R 4.3.2)
##  P Rcpp                 * 1.0.11    2023-07-06 [?] CRAN (R 4.3.2)
##  P RcppParallel           5.1.7     2023-02-27 [?] CRAN (R 4.3.2)
##  P RCurl                  1.98-1.13 2023-11-02 [?] CRAN (R 4.3.2)
##  P readr                * 2.1.5     2024-01-10 [?] CRAN (R 4.3.2)
##  P remotes                2.4.2.1   2023-07-18 [?] CRAN (R 4.3.2)
##    renv                   1.0.5     2024-02-29 [1] CRAN (R 4.3.2)
##  P reshape2               1.4.4     2020-04-09 [?] CRAN (R 4.3.2)
##  P rhdf5                  2.46.1    2023-11-29 [?] Bioconduc~
##  P rhdf5filters           1.14.1    2023-11-06 [?] Bioconductor
##  P Rhdf5lib               1.24.2    2024-02-07 [?] Bioconduc~
##  P rlang                  1.1.2     2023-11-04 [?] CRAN (R 4.3.2)
##  P rmarkdown              2.25      2023-09-18 [?] CRAN (R 4.3.2)
##  P Rsamtools              2.18.0    2023-10-24 [?] Bioconductor
##  P rstudioapi             0.15.0    2023-07-07 [?] CRAN (R 4.3.2)
##  P S4Arrays               1.2.0     2023-10-24 [?] Bioconductor
##  P S4Vectors            * 0.40.1    2023-10-26 [?] Bioconductor
##  P sass                   0.4.7     2023-07-15 [?] CRAN (R 4.3.2)
##  P scales                 1.3.0     2023-11-28 [?] CRAN (R 4.3.2)
##  P sessioninfo            1.2.2     2021-12-06 [?] CRAN (R 4.3.2)
##  P shiny                  1.7.5.1   2023-10-14 [?] CRAN (R 4.3.2)
##  P ShortRead              1.60.0    2023-10-24 [?] Bioconductor
##  P SparseArray            1.2.4     2024-02-11 [?] Bioconduc~
##  P stringi                1.7.12    2023-01-11 [?] CRAN (R 4.3.2)
##  P stringr              * 1.5.0     2022-12-02 [?] CRAN (R 4.3.2)
##  P SummarizedExperiment   1.32.0    2023-10-24 [?] Bioconductor
##  P survival               3.5-8     2024-02-14 [?] CRAN (R 4.3.3)
##  P tibble               * 3.2.1     2023-03-20 [?] CRAN (R 4.3.2)
##  P tidyr                * 1.3.1     2024-01-24 [?] CRAN (R 4.3.2)
##  P tidyselect             1.2.0     2022-10-10 [?] CRAN (R 4.3.2)
##  P tidyverse            * 2.0.0     2023-02-22 [?] CRAN (R 4.3.2)
##  P timechange             0.3.0     2024-01-18 [?] CRAN (R 4.3.2)
##  P tzdb                   0.4.0     2023-05-12 [?] CRAN (R 4.3.2)
##  P urlchecker             1.0.1     2021-11-30 [?] CRAN (R 4.3.2)
##  P usethis                2.2.2     2023-07-06 [?] CRAN (R 4.3.2)
##  P utf8                   1.2.4     2023-10-22 [?] CRAN (R 4.3.2)
##  P vctrs                  0.6.4     2023-10-12 [?] CRAN (R 4.3.2)
##  P vegan                  2.6-4     2022-10-11 [?] CRAN (R 4.3.2)
##  P vroom                  1.6.5     2023-12-05 [?] CRAN (R 4.3.2)
##  P withr                  2.5.2     2023-10-30 [?] CRAN (R 4.3.2)
##  P xfun                   0.52      2025-04-02 [?] CRAN (R 4.3.3)
##  P xtable                 1.8-4     2019-04-21 [?] CRAN (R 4.3.2)
##  P XVector              * 0.42.0    2023-10-24 [?] Bioconductor
##  P yaml                   2.3.7     2023-01-23 [?] CRAN (R 4.3.2)
##  P zlibbioc               1.48.0    2023-10-24 [?] Bioconductor
## 
##  [1] /local/workdir/arp277/Pendleton_2025_Ontario_Publication_Repo/renv/library/R-4.3/x86_64-pc-linux-gnu
##  [2] /home/arp277/.cache/R/renv/sandbox/R-4.3/x86_64-pc-linux-gnu/fd835031
## 
##  P ── Loaded and on-disk path mismatch.
## 
## ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
```
