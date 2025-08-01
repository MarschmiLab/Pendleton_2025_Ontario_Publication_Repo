---
title: "DADA2 Processing" 
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




# Loading packages 

```r
# Efficient loading of the packages 
pacman::p_load(dada2, tidyverse, patchwork, phyloseq, Biostrings, fun.gus, install = FALSE)

knitr::write_bib(file = "data/01_dada2_exports/packages.bib")

# Load in the colors and shapes
source("code/R/plotting_aesthetics.R")
```

# Goals of this document

In this document, I read in raw reads from run "231220_M05470_0458_000000000-L9NR8", Order #10471351, where I sequencing the V4 region of the 16S gene for whole-fraction samples from Lake Ontario. I clean-up metadata, examine read quality, filter, trim, remove chimeras, and assess how many reads per sample we retain. I also analyze accuracy in regards to the Mock Community.

To note, this entire process is run within the AAH_ONT_2023 home directory, which in the future will contain multiple projects that resulted from this cruise.

# Loading Metadata

In this section, I load data sheets from the lab data log, join them together, and calculate some stats including how long it took for us to filter and how much water was extracted for each sample.


```r
full_metadata <- read_csv("data/01_dada2_exports/ontario_metadata.csv")
```

# Set the path to the seq files 

This assumes that you've downloaded the raw fastq.gz files from PRJNA1212049 into the directory data/00_raw_reads


```r
# Set path to the gzipped files 
path <- "data/00_raw_reads"
path
```

```
## [1] "data/00_raw_reads"
```

```r
# What files do we have?
list.files(path)
```

```
##   [1] "Plate1_A01_AP_D1_10471351_L9NR8_R1.fastq.gz"          "Plate1_A01_AP_D1_10471351_L9NR8_R2.fastq.gz"          "Plate1_A02_AP_D9_10471351_L9NR8_R1.fastq.gz"          "Plate1_A02_AP_D9_10471351_L9NR8_R2.fastq.gz"         
##   [5] "Plate1_A03_AP_D17_10471351_L9NR8_R1.fastq.gz"         "Plate1_A03_AP_D17_10471351_L9NR8_R2.fastq.gz"         "Plate1_A04_AP_D25_10471351_L9NR8_R1.fastq.gz"         "Plate1_A04_AP_D25_10471351_L9NR8_R2.fastq.gz"        
##   [9] "Plate1_A05_AP_D33_10471351_L9NR8_R1.fastq.gz"         "Plate1_A05_AP_D33_10471351_L9NR8_R2.fastq.gz"         "Plate1_A06_AP_D41_10471351_L9NR8_R1.fastq.gz"         "Plate1_A06_AP_D41_10471351_L9NR8_R2.fastq.gz"        
##  [13] "Plate1_A07_AP_D49_10471351_L9NR8_R1.fastq.gz"         "Plate1_A07_AP_D49_10471351_L9NR8_R2.fastq.gz"         "Plate1_A08_AP_D58_10471351_L9NR8_R1.fastq.gz"         "Plate1_A08_AP_D58_10471351_L9NR8_R2.fastq.gz"        
##  [17] "Plate1_A09_AP_D72_10471351_L9NR8_R1.fastq.gz"         "Plate1_A09_AP_D72_10471351_L9NR8_R2.fastq.gz"         "Plate1_A10_AP_D80_10471351_L9NR8_R1.fastq.gz"         "Plate1_A10_AP_D80_10471351_L9NR8_R2.fastq.gz"        
##  [21] "Plate1_A11_AP_D88_10471351_L9NR8_R1.fastq.gz"         "Plate1_A11_AP_D88_10471351_L9NR8_R2.fastq.gz"         "Plate1_A12_AP_D96_10471351_L9NR8_R1.fastq.gz"         "Plate1_A12_AP_D96_10471351_L9NR8_R2.fastq.gz"        
##  [25] "Plate1_B01_AP_D2_10471351_L9NR8_R1.fastq.gz"          "Plate1_B01_AP_D2_10471351_L9NR8_R2.fastq.gz"          "Plate1_B02_AP_D10_10471351_L9NR8_R1.fastq.gz"         "Plate1_B02_AP_D10_10471351_L9NR8_R2.fastq.gz"        
##  [29] "Plate1_B03_AP_D18_10471351_L9NR8_R1.fastq.gz"         "Plate1_B03_AP_D18_10471351_L9NR8_R2.fastq.gz"         "Plate1_B04_AP_D26_10471351_L9NR8_R1.fastq.gz"         "Plate1_B04_AP_D26_10471351_L9NR8_R2.fastq.gz"        
##  [33] "Plate1_B05_AP_D34_10471351_L9NR8_R1.fastq.gz"         "Plate1_B05_AP_D34_10471351_L9NR8_R2.fastq.gz"         "Plate1_B06_AP_D42_10471351_L9NR8_R1.fastq.gz"         "Plate1_B06_AP_D42_10471351_L9NR8_R2.fastq.gz"        
##  [37] "Plate1_B07_AP_D50_10471351_L9NR8_R1.fastq.gz"         "Plate1_B07_AP_D50_10471351_L9NR8_R2.fastq.gz"         "Plate1_B08_AP_D59_10471351_L9NR8_R1.fastq.gz"         "Plate1_B08_AP_D59_10471351_L9NR8_R2.fastq.gz"        
##  [41] "Plate1_B09_AP_D73_10471351_L9NR8_R1.fastq.gz"         "Plate1_B09_AP_D73_10471351_L9NR8_R2.fastq.gz"         "Plate1_B10_AP_D81_10471351_L9NR8_R1.fastq.gz"         "Plate1_B10_AP_D81_10471351_L9NR8_R2.fastq.gz"        
##  [45] "Plate1_B11_AP_D89_10471351_L9NR8_R1.fastq.gz"         "Plate1_B11_AP_D89_10471351_L9NR8_R2.fastq.gz"         "Plate1_B12_AP_D97_10471351_L9NR8_R1.fastq.gz"         "Plate1_B12_AP_D97_10471351_L9NR8_R2.fastq.gz"        
##  [49] "Plate1_C01_AP_D3_10471351_L9NR8_R1.fastq.gz"          "Plate1_C01_AP_D3_10471351_L9NR8_R2.fastq.gz"          "Plate1_C02_AP_D11_10471351_L9NR8_R1.fastq.gz"         "Plate1_C02_AP_D11_10471351_L9NR8_R2.fastq.gz"        
##  [53] "Plate1_C03_AP_D19_10471351_L9NR8_R1.fastq.gz"         "Plate1_C03_AP_D19_10471351_L9NR8_R2.fastq.gz"         "Plate1_C04_AP_D27_10471351_L9NR8_R1.fastq.gz"         "Plate1_C04_AP_D27_10471351_L9NR8_R2.fastq.gz"        
##  [57] "Plate1_C05_AP_D35_10471351_L9NR8_R1.fastq.gz"         "Plate1_C05_AP_D35_10471351_L9NR8_R2.fastq.gz"         "Plate1_C06_AP_D43_10471351_L9NR8_R1.fastq.gz"         "Plate1_C06_AP_D43_10471351_L9NR8_R2.fastq.gz"        
##  [61] "Plate1_C07_AP_D51_10471351_L9NR8_R1.fastq.gz"         "Plate1_C07_AP_D51_10471351_L9NR8_R2.fastq.gz"         "Plate1_C08_AP_D60_10471351_L9NR8_R1.fastq.gz"         "Plate1_C08_AP_D60_10471351_L9NR8_R2.fastq.gz"        
##  [65] "Plate1_C09_AP_D74_10471351_L9NR8_R1.fastq.gz"         "Plate1_C09_AP_D74_10471351_L9NR8_R2.fastq.gz"         "Plate1_C10_AP_D82_10471351_L9NR8_R1.fastq.gz"         "Plate1_C10_AP_D82_10471351_L9NR8_R2.fastq.gz"        
##  [69] "Plate1_C11_AP_D90_10471351_L9NR8_R1.fastq.gz"         "Plate1_C11_AP_D90_10471351_L9NR8_R2.fastq.gz"         "Plate1_C12_AP_D98_10471351_L9NR8_R1.fastq.gz"         "Plate1_C12_AP_D98_10471351_L9NR8_R2.fastq.gz"        
##  [73] "Plate1_D01_AP_D4_10471351_L9NR8_R1.fastq.gz"          "Plate1_D01_AP_D4_10471351_L9NR8_R2.fastq.gz"          "Plate1_D02_AP_D12_10471351_L9NR8_R1.fastq.gz"         "Plate1_D02_AP_D12_10471351_L9NR8_R2.fastq.gz"        
##  [77] "Plate1_D03_AP_D20_10471351_L9NR8_R1.fastq.gz"         "Plate1_D03_AP_D20_10471351_L9NR8_R2.fastq.gz"         "Plate1_D04_AP_D28_10471351_L9NR8_R1.fastq.gz"         "Plate1_D04_AP_D28_10471351_L9NR8_R2.fastq.gz"        
##  [81] "Plate1_D05_AP_D36_10471351_L9NR8_R1.fastq.gz"         "Plate1_D05_AP_D36_10471351_L9NR8_R2.fastq.gz"         "Plate1_D06_AP_D44_10471351_L9NR8_R1.fastq.gz"         "Plate1_D06_AP_D44_10471351_L9NR8_R2.fastq.gz"        
##  [85] "Plate1_D07_AP_D52_10471351_L9NR8_R1.fastq.gz"         "Plate1_D07_AP_D52_10471351_L9NR8_R2.fastq.gz"         "Plate1_D08_AP_D61_10471351_L9NR8_R1.fastq.gz"         "Plate1_D08_AP_D61_10471351_L9NR8_R2.fastq.gz"        
##  [89] "Plate1_D09_AP_D75_10471351_L9NR8_R1.fastq.gz"         "Plate1_D09_AP_D75_10471351_L9NR8_R2.fastq.gz"         "Plate1_D10_AP_D83_10471351_L9NR8_R1.fastq.gz"         "Plate1_D10_AP_D83_10471351_L9NR8_R2.fastq.gz"        
##  [93] "Plate1_D11_AP_D91_10471351_L9NR8_R1.fastq.gz"         "Plate1_D11_AP_D91_10471351_L9NR8_R2.fastq.gz"         "Plate1_D12_AP_D99_10471351_L9NR8_R1.fastq.gz"         "Plate1_D12_AP_D99_10471351_L9NR8_R2.fastq.gz"        
##  [97] "Plate1_E01_AP_D5_10471351_L9NR8_R1.fastq.gz"          "Plate1_E01_AP_D5_10471351_L9NR8_R2.fastq.gz"          "Plate1_E02_AP_D13_10471351_L9NR8_R1.fastq.gz"         "Plate1_E02_AP_D13_10471351_L9NR8_R2.fastq.gz"        
## [101] "Plate1_E03_AP_D21_10471351_L9NR8_R1.fastq.gz"         "Plate1_E03_AP_D21_10471351_L9NR8_R2.fastq.gz"         "Plate1_E04_AP_D29_10471351_L9NR8_R1.fastq.gz"         "Plate1_E04_AP_D29_10471351_L9NR8_R2.fastq.gz"        
## [105] "Plate1_E05_AP_D37_10471351_L9NR8_R1.fastq.gz"         "Plate1_E05_AP_D37_10471351_L9NR8_R2.fastq.gz"         "Plate1_E06_AP_D45_10471351_L9NR8_R1.fastq.gz"         "Plate1_E06_AP_D45_10471351_L9NR8_R2.fastq.gz"        
## [109] "Plate1_E07_AP_D53_10471351_L9NR8_R1.fastq.gz"         "Plate1_E07_AP_D53_10471351_L9NR8_R2.fastq.gz"         "Plate1_E08_AP_D62_10471351_L9NR8_R1.fastq.gz"         "Plate1_E08_AP_D62_10471351_L9NR8_R2.fastq.gz"        
## [113] "Plate1_E09_AP_D76_10471351_L9NR8_R1.fastq.gz"         "Plate1_E09_AP_D76_10471351_L9NR8_R2.fastq.gz"         "Plate1_E10_AP_D84_10471351_L9NR8_R1.fastq.gz"         "Plate1_E10_AP_D84_10471351_L9NR8_R2.fastq.gz"        
## [117] "Plate1_E11_AP_D92_10471351_L9NR8_R1.fastq.gz"         "Plate1_E11_AP_D92_10471351_L9NR8_R2.fastq.gz"         "Plate1_E12_AP_D100_10471351_L9NR8_R1.fastq.gz"        "Plate1_E12_AP_D100_10471351_L9NR8_R2.fastq.gz"       
## [121] "Plate1_F01_AP_D6_10471351_L9NR8_R1.fastq.gz"          "Plate1_F01_AP_D6_10471351_L9NR8_R2.fastq.gz"          "Plate1_F02_AP_D14_10471351_L9NR8_R1.fastq.gz"         "Plate1_F02_AP_D14_10471351_L9NR8_R2.fastq.gz"        
## [125] "Plate1_F03_AP_D22_10471351_L9NR8_R1.fastq.gz"         "Plate1_F03_AP_D22_10471351_L9NR8_R2.fastq.gz"         "Plate1_F04_AP_D30_10471351_L9NR8_R1.fastq.gz"         "Plate1_F04_AP_D30_10471351_L9NR8_R2.fastq.gz"        
## [129] "Plate1_F05_AP_D38_10471351_L9NR8_R1.fastq.gz"         "Plate1_F05_AP_D38_10471351_L9NR8_R2.fastq.gz"         "Plate1_F06_AP_D46_10471351_L9NR8_R1.fastq.gz"         "Plate1_F06_AP_D46_10471351_L9NR8_R2.fastq.gz"        
## [133] "Plate1_F07_AP_D54_10471351_L9NR8_R1.fastq.gz"         "Plate1_F07_AP_D54_10471351_L9NR8_R2.fastq.gz"         "Plate1_F08_AP_D69_10471351_L9NR8_R1.fastq.gz"         "Plate1_F08_AP_D69_10471351_L9NR8_R2.fastq.gz"        
## [137] "Plate1_F09_AP_D77_10471351_L9NR8_R1.fastq.gz"         "Plate1_F09_AP_D77_10471351_L9NR8_R2.fastq.gz"         "Plate1_F10_AP_D85_10471351_L9NR8_R1.fastq.gz"         "Plate1_F10_AP_D85_10471351_L9NR8_R2.fastq.gz"        
## [141] "Plate1_F11_AP_D93_10471351_L9NR8_R1.fastq.gz"         "Plate1_F11_AP_D93_10471351_L9NR8_R2.fastq.gz"         "Plate1_F12_AP_D101_10471351_L9NR8_R1.fastq.gz"        "Plate1_F12_AP_D101_10471351_L9NR8_R2.fastq.gz"       
## [145] "Plate1_G01_AP_D7_10471351_L9NR8_R1.fastq.gz"          "Plate1_G01_AP_D7_10471351_L9NR8_R2.fastq.gz"          "Plate1_G02_AP_D15_10471351_L9NR8_R1.fastq.gz"         "Plate1_G02_AP_D15_10471351_L9NR8_R2.fastq.gz"        
## [149] "Plate1_G03_AP_D23_10471351_L9NR8_R1.fastq.gz"         "Plate1_G03_AP_D23_10471351_L9NR8_R2.fastq.gz"         "Plate1_G04_AP_D31_10471351_L9NR8_R1.fastq.gz"         "Plate1_G04_AP_D31_10471351_L9NR8_R2.fastq.gz"        
## [153] "Plate1_G05_AP_D39_10471351_L9NR8_R1.fastq.gz"         "Plate1_G05_AP_D39_10471351_L9NR8_R2.fastq.gz"         "Plate1_G06_AP_D47_10471351_L9NR8_R1.fastq.gz"         "Plate1_G06_AP_D47_10471351_L9NR8_R2.fastq.gz"        
## [157] "Plate1_G07_AP_D55_10471351_L9NR8_R1.fastq.gz"         "Plate1_G07_AP_D55_10471351_L9NR8_R2.fastq.gz"         "Plate1_G08_AP_D70_10471351_L9NR8_R1.fastq.gz"         "Plate1_G08_AP_D70_10471351_L9NR8_R2.fastq.gz"        
## [161] "Plate1_G09_AP_D78_10471351_L9NR8_R1.fastq.gz"         "Plate1_G09_AP_D78_10471351_L9NR8_R2.fastq.gz"         "Plate1_G10_AP_D86_10471351_L9NR8_R1.fastq.gz"         "Plate1_G10_AP_D86_10471351_L9NR8_R2.fastq.gz"        
## [165] "Plate1_G11_AP_D94_10471351_L9NR8_R1.fastq.gz"         "Plate1_G11_AP_D94_10471351_L9NR8_R2.fastq.gz"         "Plate1_G12_AP_D102_10471351_L9NR8_R1.fastq.gz"        "Plate1_G12_AP_D102_10471351_L9NR8_R2.fastq.gz"       
## [169] "Plate1_H01_AP_D8_10471351_L9NR8_R1.fastq.gz"          "Plate1_H01_AP_D8_10471351_L9NR8_R2.fastq.gz"          "Plate1_H02_AP_D16_10471351_L9NR8_R1.fastq.gz"         "Plate1_H02_AP_D16_10471351_L9NR8_R2.fastq.gz"        
## [173] "Plate1_H03_AP_D24_10471351_L9NR8_R1.fastq.gz"         "Plate1_H03_AP_D24_10471351_L9NR8_R2.fastq.gz"         "Plate1_H04_AP_D32_10471351_L9NR8_R1.fastq.gz"         "Plate1_H04_AP_D32_10471351_L9NR8_R2.fastq.gz"        
## [177] "Plate1_H05_AP_D40_10471351_L9NR8_R1.fastq.gz"         "Plate1_H05_AP_D40_10471351_L9NR8_R2.fastq.gz"         "Plate1_H06_AP_D48_10471351_L9NR8_R1.fastq.gz"         "Plate1_H06_AP_D48_10471351_L9NR8_R2.fastq.gz"        
## [181] "Plate1_H07_AP_D57_10471351_L9NR8_R1.fastq.gz"         "Plate1_H07_AP_D57_10471351_L9NR8_R2.fastq.gz"         "Plate1_H08_AP_D71_10471351_L9NR8_R1.fastq.gz"         "Plate1_H08_AP_D71_10471351_L9NR8_R2.fastq.gz"        
## [185] "Plate1_H09_AP_D79_10471351_L9NR8_R1.fastq.gz"         "Plate1_H09_AP_D79_10471351_L9NR8_R2.fastq.gz"         "Plate1_H10_AP_D87_10471351_L9NR8_R1.fastq.gz"         "Plate1_H10_AP_D87_10471351_L9NR8_R2.fastq.gz"        
## [189] "Plate1_H11_AP_D95_10471351_L9NR8_R1.fastq.gz"         "Plate1_H11_AP_D95_10471351_L9NR8_R2.fastq.gz"         "Plate1_H12_PCR_Blanks_10471351_L9NR8_R1.fastq.gz"     "Plate1_H12_PCR_Blanks_10471351_L9NR8_R2.fastq.gz"    
## [193] "Plate2_1351_A01_AP_D103_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_A01_AP_D103_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_A02_AP_D111_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_A02_AP_D111_10471351_L9NR8_R2.fastq.gz"  
## [197] "Plate2_1351_A03_AP_D119_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_A03_AP_D119_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_A04_AP_D127_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_A04_AP_D127_10471351_L9NR8_R2.fastq.gz"  
## [201] "Plate2_1351_A05_AP_D135_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_A05_AP_D135_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_A06_AP_D143_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_A06_AP_D143_10471351_L9NR8_R2.fastq.gz"  
## [205] "Plate2_1351_A07_AP_D151_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_A07_AP_D151_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_B01_AP_D104_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_B01_AP_D104_10471351_L9NR8_R2.fastq.gz"  
## [209] "Plate2_1351_B02_AP_D112_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_B02_AP_D112_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_B03_AP_D120_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_B03_AP_D120_10471351_L9NR8_R2.fastq.gz"  
## [213] "Plate2_1351_B04_AP_D128_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_B04_AP_D128_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_B05_AP_D136_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_B05_AP_D136_10471351_L9NR8_R2.fastq.gz"  
## [217] "Plate2_1351_B06_AP_D144_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_B06_AP_D144_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_B07_AP_D152_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_B07_AP_D152_10471351_L9NR8_R2.fastq.gz"  
## [221] "Plate2_1351_C01_AP_D105_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_C01_AP_D105_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_C02_AP_D113_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_C02_AP_D113_10471351_L9NR8_R2.fastq.gz"  
## [225] "Plate2_1351_C03_AP_D121_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_C03_AP_D121_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_C04_AP_D129_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_C04_AP_D129_10471351_L9NR8_R2.fastq.gz"  
## [229] "Plate2_1351_C05_AP_D137_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_C05_AP_D137_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_C06_AP_D145_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_C06_AP_D145_10471351_L9NR8_R2.fastq.gz"  
## [233] "Plate2_1351_C07_AP_D153_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_C07_AP_D153_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_D01_AP_D106_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_D01_AP_D106_10471351_L9NR8_R2.fastq.gz"  
## [237] "Plate2_1351_D02_AP_D114_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_D02_AP_D114_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_D03_AP_D122_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_D03_AP_D122_10471351_L9NR8_R2.fastq.gz"  
## [241] "Plate2_1351_D04_AP_D130_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_D04_AP_D130_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_D05_AP_D138_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_D05_AP_D138_10471351_L9NR8_R2.fastq.gz"  
## [245] "Plate2_1351_D06_AP_D146_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_D06_AP_D146_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_D07_AP_D154_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_D07_AP_D154_10471351_L9NR8_R2.fastq.gz"  
## [249] "Plate2_1351_E01_AP_D107_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_E01_AP_D107_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_E02_AP_D115_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_E02_AP_D115_10471351_L9NR8_R2.fastq.gz"  
## [253] "Plate2_1351_E03_AP_D123_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_E03_AP_D123_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_E04_AP_D131_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_E04_AP_D131_10471351_L9NR8_R2.fastq.gz"  
## [257] "Plate2_1351_E05_AP_D139_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_E05_AP_D139_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_E06_AP_D147_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_E06_AP_D147_10471351_L9NR8_R2.fastq.gz"  
## [261] "Plate2_1351_E07_Zymo_Mock_10471351_L9NR8_R1.fastq.gz" "Plate2_1351_E07_Zymo_Mock_10471351_L9NR8_R2.fastq.gz" "Plate2_1351_F01_AP_D108_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_F01_AP_D108_10471351_L9NR8_R2.fastq.gz"  
## [265] "Plate2_1351_F02_AP_D116_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_F02_AP_D116_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_F03_AP_D124_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_F03_AP_D124_10471351_L9NR8_R2.fastq.gz"  
## [269] "Plate2_1351_F04_AP_D132_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_F04_AP_D132_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_F05_AP_D140_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_F05_AP_D140_10471351_L9NR8_R2.fastq.gz"  
## [273] "Plate2_1351_F06_AP_D148_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_F06_AP_D148_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_F07_Blank1_10471351_L9NR8_R1.fastq.gz"    "Plate2_1351_F07_Blank1_10471351_L9NR8_R2.fastq.gz"   
## [277] "Plate2_1351_G01_AP_D109_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_G01_AP_D109_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_G02_AP_D117_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_G02_AP_D117_10471351_L9NR8_R2.fastq.gz"  
## [281] "Plate2_1351_G03_AP_D125_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_G03_AP_D125_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_G04_AP_D133_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_G04_AP_D133_10471351_L9NR8_R2.fastq.gz"  
## [285] "Plate2_1351_G05_AP_D141_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_G05_AP_D141_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_G06_AP_D149_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_G06_AP_D149_10471351_L9NR8_R2.fastq.gz"  
## [289] "Plate2_1351_G07_Blank2_10471351_L9NR8_R1.fastq.gz"    "Plate2_1351_G07_Blank2_10471351_L9NR8_R2.fastq.gz"    "Plate2_1351_H01_AP_D110_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_H01_AP_D110_10471351_L9NR8_R2.fastq.gz"  
## [293] "Plate2_1351_H02_AP_D118_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_H02_AP_D118_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_H03_AP_D126_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_H03_AP_D126_10471351_L9NR8_R2.fastq.gz"  
## [297] "Plate2_1351_H04_AP_D134_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_H04_AP_D134_10471351_L9NR8_R2.fastq.gz"   "Plate2_1351_H05_AP_D142_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_H05_AP_D142_10471351_L9NR8_R2.fastq.gz"  
## [301] "Plate2_1351_H06_AP_D150_10471351_L9NR8_R1.fastq.gz"   "Plate2_1351_H06_AP_D150_10471351_L9NR8_R2.fastq.gz"
```

We have 151 samples, with forward and reverse reads for each. These include smaples included in our DNA_Log, as well as PCR_Blanks (from first-round PCR), "Blank1" and "Blank2" which are blanks from the indexing step, and "Zymo_Mock" which is our Zymo Microbiomics Mock Community.

# Load in Forward and Reverse reads and assess the quality


```r
# Create variable for the forward and the reverse reads

# 1. Forward read variable 
forward_reads <- sort(list.files(path, pattern = "R1.fastq.gz", 
                      full.names = TRUE))
forward_reads
```

```
##   [1] "data/00_raw_reads/Plate1_A01_AP_D1_10471351_L9NR8_R1.fastq.gz"          "data/00_raw_reads/Plate1_A02_AP_D9_10471351_L9NR8_R1.fastq.gz"          "data/00_raw_reads/Plate1_A03_AP_D17_10471351_L9NR8_R1.fastq.gz"        
##   [4] "data/00_raw_reads/Plate1_A04_AP_D25_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_A05_AP_D33_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_A06_AP_D41_10471351_L9NR8_R1.fastq.gz"        
##   [7] "data/00_raw_reads/Plate1_A07_AP_D49_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_A08_AP_D58_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_A09_AP_D72_10471351_L9NR8_R1.fastq.gz"        
##  [10] "data/00_raw_reads/Plate1_A10_AP_D80_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_A11_AP_D88_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_A12_AP_D96_10471351_L9NR8_R1.fastq.gz"        
##  [13] "data/00_raw_reads/Plate1_B01_AP_D2_10471351_L9NR8_R1.fastq.gz"          "data/00_raw_reads/Plate1_B02_AP_D10_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_B03_AP_D18_10471351_L9NR8_R1.fastq.gz"        
##  [16] "data/00_raw_reads/Plate1_B04_AP_D26_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_B05_AP_D34_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_B06_AP_D42_10471351_L9NR8_R1.fastq.gz"        
##  [19] "data/00_raw_reads/Plate1_B07_AP_D50_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_B08_AP_D59_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_B09_AP_D73_10471351_L9NR8_R1.fastq.gz"        
##  [22] "data/00_raw_reads/Plate1_B10_AP_D81_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_B11_AP_D89_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_B12_AP_D97_10471351_L9NR8_R1.fastq.gz"        
##  [25] "data/00_raw_reads/Plate1_C01_AP_D3_10471351_L9NR8_R1.fastq.gz"          "data/00_raw_reads/Plate1_C02_AP_D11_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_C03_AP_D19_10471351_L9NR8_R1.fastq.gz"        
##  [28] "data/00_raw_reads/Plate1_C04_AP_D27_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_C05_AP_D35_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_C06_AP_D43_10471351_L9NR8_R1.fastq.gz"        
##  [31] "data/00_raw_reads/Plate1_C07_AP_D51_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_C08_AP_D60_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_C09_AP_D74_10471351_L9NR8_R1.fastq.gz"        
##  [34] "data/00_raw_reads/Plate1_C10_AP_D82_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_C11_AP_D90_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_C12_AP_D98_10471351_L9NR8_R1.fastq.gz"        
##  [37] "data/00_raw_reads/Plate1_D01_AP_D4_10471351_L9NR8_R1.fastq.gz"          "data/00_raw_reads/Plate1_D02_AP_D12_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_D03_AP_D20_10471351_L9NR8_R1.fastq.gz"        
##  [40] "data/00_raw_reads/Plate1_D04_AP_D28_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_D05_AP_D36_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_D06_AP_D44_10471351_L9NR8_R1.fastq.gz"        
##  [43] "data/00_raw_reads/Plate1_D07_AP_D52_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_D08_AP_D61_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_D09_AP_D75_10471351_L9NR8_R1.fastq.gz"        
##  [46] "data/00_raw_reads/Plate1_D10_AP_D83_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_D11_AP_D91_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_D12_AP_D99_10471351_L9NR8_R1.fastq.gz"        
##  [49] "data/00_raw_reads/Plate1_E01_AP_D5_10471351_L9NR8_R1.fastq.gz"          "data/00_raw_reads/Plate1_E02_AP_D13_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_E03_AP_D21_10471351_L9NR8_R1.fastq.gz"        
##  [52] "data/00_raw_reads/Plate1_E04_AP_D29_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_E05_AP_D37_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_E06_AP_D45_10471351_L9NR8_R1.fastq.gz"        
##  [55] "data/00_raw_reads/Plate1_E07_AP_D53_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_E08_AP_D62_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_E09_AP_D76_10471351_L9NR8_R1.fastq.gz"        
##  [58] "data/00_raw_reads/Plate1_E10_AP_D84_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_E11_AP_D92_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_E12_AP_D100_10471351_L9NR8_R1.fastq.gz"       
##  [61] "data/00_raw_reads/Plate1_F01_AP_D6_10471351_L9NR8_R1.fastq.gz"          "data/00_raw_reads/Plate1_F02_AP_D14_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_F03_AP_D22_10471351_L9NR8_R1.fastq.gz"        
##  [64] "data/00_raw_reads/Plate1_F04_AP_D30_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_F05_AP_D38_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_F06_AP_D46_10471351_L9NR8_R1.fastq.gz"        
##  [67] "data/00_raw_reads/Plate1_F07_AP_D54_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_F08_AP_D69_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_F09_AP_D77_10471351_L9NR8_R1.fastq.gz"        
##  [70] "data/00_raw_reads/Plate1_F10_AP_D85_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_F11_AP_D93_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_F12_AP_D101_10471351_L9NR8_R1.fastq.gz"       
##  [73] "data/00_raw_reads/Plate1_G01_AP_D7_10471351_L9NR8_R1.fastq.gz"          "data/00_raw_reads/Plate1_G02_AP_D15_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_G03_AP_D23_10471351_L9NR8_R1.fastq.gz"        
##  [76] "data/00_raw_reads/Plate1_G04_AP_D31_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_G05_AP_D39_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_G06_AP_D47_10471351_L9NR8_R1.fastq.gz"        
##  [79] "data/00_raw_reads/Plate1_G07_AP_D55_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_G08_AP_D70_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_G09_AP_D78_10471351_L9NR8_R1.fastq.gz"        
##  [82] "data/00_raw_reads/Plate1_G10_AP_D86_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_G11_AP_D94_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_G12_AP_D102_10471351_L9NR8_R1.fastq.gz"       
##  [85] "data/00_raw_reads/Plate1_H01_AP_D8_10471351_L9NR8_R1.fastq.gz"          "data/00_raw_reads/Plate1_H02_AP_D16_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_H03_AP_D24_10471351_L9NR8_R1.fastq.gz"        
##  [88] "data/00_raw_reads/Plate1_H04_AP_D32_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_H05_AP_D40_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_H06_AP_D48_10471351_L9NR8_R1.fastq.gz"        
##  [91] "data/00_raw_reads/Plate1_H07_AP_D57_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_H08_AP_D71_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_H09_AP_D79_10471351_L9NR8_R1.fastq.gz"        
##  [94] "data/00_raw_reads/Plate1_H10_AP_D87_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_H11_AP_D95_10471351_L9NR8_R1.fastq.gz"         "data/00_raw_reads/Plate1_H12_PCR_Blanks_10471351_L9NR8_R1.fastq.gz"    
##  [97] "data/00_raw_reads/Plate2_1351_A01_AP_D103_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_A02_AP_D111_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_A03_AP_D119_10471351_L9NR8_R1.fastq.gz"  
## [100] "data/00_raw_reads/Plate2_1351_A04_AP_D127_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_A05_AP_D135_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_A06_AP_D143_10471351_L9NR8_R1.fastq.gz"  
## [103] "data/00_raw_reads/Plate2_1351_A07_AP_D151_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_B01_AP_D104_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_B02_AP_D112_10471351_L9NR8_R1.fastq.gz"  
## [106] "data/00_raw_reads/Plate2_1351_B03_AP_D120_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_B04_AP_D128_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_B05_AP_D136_10471351_L9NR8_R1.fastq.gz"  
## [109] "data/00_raw_reads/Plate2_1351_B06_AP_D144_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_B07_AP_D152_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_C01_AP_D105_10471351_L9NR8_R1.fastq.gz"  
## [112] "data/00_raw_reads/Plate2_1351_C02_AP_D113_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_C03_AP_D121_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_C04_AP_D129_10471351_L9NR8_R1.fastq.gz"  
## [115] "data/00_raw_reads/Plate2_1351_C05_AP_D137_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_C06_AP_D145_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_C07_AP_D153_10471351_L9NR8_R1.fastq.gz"  
## [118] "data/00_raw_reads/Plate2_1351_D01_AP_D106_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_D02_AP_D114_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_D03_AP_D122_10471351_L9NR8_R1.fastq.gz"  
## [121] "data/00_raw_reads/Plate2_1351_D04_AP_D130_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_D05_AP_D138_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_D06_AP_D146_10471351_L9NR8_R1.fastq.gz"  
## [124] "data/00_raw_reads/Plate2_1351_D07_AP_D154_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_E01_AP_D107_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_E02_AP_D115_10471351_L9NR8_R1.fastq.gz"  
## [127] "data/00_raw_reads/Plate2_1351_E03_AP_D123_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_E04_AP_D131_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_E05_AP_D139_10471351_L9NR8_R1.fastq.gz"  
## [130] "data/00_raw_reads/Plate2_1351_E06_AP_D147_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_E07_Zymo_Mock_10471351_L9NR8_R1.fastq.gz" "data/00_raw_reads/Plate2_1351_F01_AP_D108_10471351_L9NR8_R1.fastq.gz"  
## [133] "data/00_raw_reads/Plate2_1351_F02_AP_D116_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_F03_AP_D124_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_F04_AP_D132_10471351_L9NR8_R1.fastq.gz"  
## [136] "data/00_raw_reads/Plate2_1351_F05_AP_D140_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_F06_AP_D148_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_F07_Blank1_10471351_L9NR8_R1.fastq.gz"   
## [139] "data/00_raw_reads/Plate2_1351_G01_AP_D109_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_G02_AP_D117_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_G03_AP_D125_10471351_L9NR8_R1.fastq.gz"  
## [142] "data/00_raw_reads/Plate2_1351_G04_AP_D133_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_G05_AP_D141_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_G06_AP_D149_10471351_L9NR8_R1.fastq.gz"  
## [145] "data/00_raw_reads/Plate2_1351_G07_Blank2_10471351_L9NR8_R1.fastq.gz"    "data/00_raw_reads/Plate2_1351_H01_AP_D110_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_H02_AP_D118_10471351_L9NR8_R1.fastq.gz"  
## [148] "data/00_raw_reads/Plate2_1351_H03_AP_D126_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_H04_AP_D134_10471351_L9NR8_R1.fastq.gz"   "data/00_raw_reads/Plate2_1351_H05_AP_D142_10471351_L9NR8_R1.fastq.gz"  
## [151] "data/00_raw_reads/Plate2_1351_H06_AP_D150_10471351_L9NR8_R1.fastq.gz"
```

```r
# 2. Reverse read variable 
reverse_reads <- sort(list.files(path, pattern = "R2.fastq.gz", 
                      full.names = TRUE))
reverse_reads
```

```
##   [1] "data/00_raw_reads/Plate1_A01_AP_D1_10471351_L9NR8_R2.fastq.gz"          "data/00_raw_reads/Plate1_A02_AP_D9_10471351_L9NR8_R2.fastq.gz"          "data/00_raw_reads/Plate1_A03_AP_D17_10471351_L9NR8_R2.fastq.gz"        
##   [4] "data/00_raw_reads/Plate1_A04_AP_D25_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_A05_AP_D33_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_A06_AP_D41_10471351_L9NR8_R2.fastq.gz"        
##   [7] "data/00_raw_reads/Plate1_A07_AP_D49_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_A08_AP_D58_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_A09_AP_D72_10471351_L9NR8_R2.fastq.gz"        
##  [10] "data/00_raw_reads/Plate1_A10_AP_D80_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_A11_AP_D88_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_A12_AP_D96_10471351_L9NR8_R2.fastq.gz"        
##  [13] "data/00_raw_reads/Plate1_B01_AP_D2_10471351_L9NR8_R2.fastq.gz"          "data/00_raw_reads/Plate1_B02_AP_D10_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_B03_AP_D18_10471351_L9NR8_R2.fastq.gz"        
##  [16] "data/00_raw_reads/Plate1_B04_AP_D26_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_B05_AP_D34_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_B06_AP_D42_10471351_L9NR8_R2.fastq.gz"        
##  [19] "data/00_raw_reads/Plate1_B07_AP_D50_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_B08_AP_D59_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_B09_AP_D73_10471351_L9NR8_R2.fastq.gz"        
##  [22] "data/00_raw_reads/Plate1_B10_AP_D81_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_B11_AP_D89_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_B12_AP_D97_10471351_L9NR8_R2.fastq.gz"        
##  [25] "data/00_raw_reads/Plate1_C01_AP_D3_10471351_L9NR8_R2.fastq.gz"          "data/00_raw_reads/Plate1_C02_AP_D11_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_C03_AP_D19_10471351_L9NR8_R2.fastq.gz"        
##  [28] "data/00_raw_reads/Plate1_C04_AP_D27_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_C05_AP_D35_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_C06_AP_D43_10471351_L9NR8_R2.fastq.gz"        
##  [31] "data/00_raw_reads/Plate1_C07_AP_D51_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_C08_AP_D60_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_C09_AP_D74_10471351_L9NR8_R2.fastq.gz"        
##  [34] "data/00_raw_reads/Plate1_C10_AP_D82_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_C11_AP_D90_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_C12_AP_D98_10471351_L9NR8_R2.fastq.gz"        
##  [37] "data/00_raw_reads/Plate1_D01_AP_D4_10471351_L9NR8_R2.fastq.gz"          "data/00_raw_reads/Plate1_D02_AP_D12_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_D03_AP_D20_10471351_L9NR8_R2.fastq.gz"        
##  [40] "data/00_raw_reads/Plate1_D04_AP_D28_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_D05_AP_D36_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_D06_AP_D44_10471351_L9NR8_R2.fastq.gz"        
##  [43] "data/00_raw_reads/Plate1_D07_AP_D52_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_D08_AP_D61_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_D09_AP_D75_10471351_L9NR8_R2.fastq.gz"        
##  [46] "data/00_raw_reads/Plate1_D10_AP_D83_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_D11_AP_D91_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_D12_AP_D99_10471351_L9NR8_R2.fastq.gz"        
##  [49] "data/00_raw_reads/Plate1_E01_AP_D5_10471351_L9NR8_R2.fastq.gz"          "data/00_raw_reads/Plate1_E02_AP_D13_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_E03_AP_D21_10471351_L9NR8_R2.fastq.gz"        
##  [52] "data/00_raw_reads/Plate1_E04_AP_D29_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_E05_AP_D37_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_E06_AP_D45_10471351_L9NR8_R2.fastq.gz"        
##  [55] "data/00_raw_reads/Plate1_E07_AP_D53_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_E08_AP_D62_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_E09_AP_D76_10471351_L9NR8_R2.fastq.gz"        
##  [58] "data/00_raw_reads/Plate1_E10_AP_D84_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_E11_AP_D92_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_E12_AP_D100_10471351_L9NR8_R2.fastq.gz"       
##  [61] "data/00_raw_reads/Plate1_F01_AP_D6_10471351_L9NR8_R2.fastq.gz"          "data/00_raw_reads/Plate1_F02_AP_D14_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_F03_AP_D22_10471351_L9NR8_R2.fastq.gz"        
##  [64] "data/00_raw_reads/Plate1_F04_AP_D30_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_F05_AP_D38_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_F06_AP_D46_10471351_L9NR8_R2.fastq.gz"        
##  [67] "data/00_raw_reads/Plate1_F07_AP_D54_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_F08_AP_D69_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_F09_AP_D77_10471351_L9NR8_R2.fastq.gz"        
##  [70] "data/00_raw_reads/Plate1_F10_AP_D85_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_F11_AP_D93_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_F12_AP_D101_10471351_L9NR8_R2.fastq.gz"       
##  [73] "data/00_raw_reads/Plate1_G01_AP_D7_10471351_L9NR8_R2.fastq.gz"          "data/00_raw_reads/Plate1_G02_AP_D15_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_G03_AP_D23_10471351_L9NR8_R2.fastq.gz"        
##  [76] "data/00_raw_reads/Plate1_G04_AP_D31_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_G05_AP_D39_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_G06_AP_D47_10471351_L9NR8_R2.fastq.gz"        
##  [79] "data/00_raw_reads/Plate1_G07_AP_D55_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_G08_AP_D70_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_G09_AP_D78_10471351_L9NR8_R2.fastq.gz"        
##  [82] "data/00_raw_reads/Plate1_G10_AP_D86_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_G11_AP_D94_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_G12_AP_D102_10471351_L9NR8_R2.fastq.gz"       
##  [85] "data/00_raw_reads/Plate1_H01_AP_D8_10471351_L9NR8_R2.fastq.gz"          "data/00_raw_reads/Plate1_H02_AP_D16_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_H03_AP_D24_10471351_L9NR8_R2.fastq.gz"        
##  [88] "data/00_raw_reads/Plate1_H04_AP_D32_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_H05_AP_D40_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_H06_AP_D48_10471351_L9NR8_R2.fastq.gz"        
##  [91] "data/00_raw_reads/Plate1_H07_AP_D57_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_H08_AP_D71_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_H09_AP_D79_10471351_L9NR8_R2.fastq.gz"        
##  [94] "data/00_raw_reads/Plate1_H10_AP_D87_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_H11_AP_D95_10471351_L9NR8_R2.fastq.gz"         "data/00_raw_reads/Plate1_H12_PCR_Blanks_10471351_L9NR8_R2.fastq.gz"    
##  [97] "data/00_raw_reads/Plate2_1351_A01_AP_D103_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_A02_AP_D111_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_A03_AP_D119_10471351_L9NR8_R2.fastq.gz"  
## [100] "data/00_raw_reads/Plate2_1351_A04_AP_D127_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_A05_AP_D135_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_A06_AP_D143_10471351_L9NR8_R2.fastq.gz"  
## [103] "data/00_raw_reads/Plate2_1351_A07_AP_D151_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_B01_AP_D104_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_B02_AP_D112_10471351_L9NR8_R2.fastq.gz"  
## [106] "data/00_raw_reads/Plate2_1351_B03_AP_D120_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_B04_AP_D128_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_B05_AP_D136_10471351_L9NR8_R2.fastq.gz"  
## [109] "data/00_raw_reads/Plate2_1351_B06_AP_D144_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_B07_AP_D152_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_C01_AP_D105_10471351_L9NR8_R2.fastq.gz"  
## [112] "data/00_raw_reads/Plate2_1351_C02_AP_D113_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_C03_AP_D121_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_C04_AP_D129_10471351_L9NR8_R2.fastq.gz"  
## [115] "data/00_raw_reads/Plate2_1351_C05_AP_D137_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_C06_AP_D145_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_C07_AP_D153_10471351_L9NR8_R2.fastq.gz"  
## [118] "data/00_raw_reads/Plate2_1351_D01_AP_D106_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_D02_AP_D114_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_D03_AP_D122_10471351_L9NR8_R2.fastq.gz"  
## [121] "data/00_raw_reads/Plate2_1351_D04_AP_D130_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_D05_AP_D138_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_D06_AP_D146_10471351_L9NR8_R2.fastq.gz"  
## [124] "data/00_raw_reads/Plate2_1351_D07_AP_D154_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_E01_AP_D107_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_E02_AP_D115_10471351_L9NR8_R2.fastq.gz"  
## [127] "data/00_raw_reads/Plate2_1351_E03_AP_D123_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_E04_AP_D131_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_E05_AP_D139_10471351_L9NR8_R2.fastq.gz"  
## [130] "data/00_raw_reads/Plate2_1351_E06_AP_D147_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_E07_Zymo_Mock_10471351_L9NR8_R2.fastq.gz" "data/00_raw_reads/Plate2_1351_F01_AP_D108_10471351_L9NR8_R2.fastq.gz"  
## [133] "data/00_raw_reads/Plate2_1351_F02_AP_D116_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_F03_AP_D124_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_F04_AP_D132_10471351_L9NR8_R2.fastq.gz"  
## [136] "data/00_raw_reads/Plate2_1351_F05_AP_D140_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_F06_AP_D148_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_F07_Blank1_10471351_L9NR8_R2.fastq.gz"   
## [139] "data/00_raw_reads/Plate2_1351_G01_AP_D109_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_G02_AP_D117_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_G03_AP_D125_10471351_L9NR8_R2.fastq.gz"  
## [142] "data/00_raw_reads/Plate2_1351_G04_AP_D133_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_G05_AP_D141_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_G06_AP_D149_10471351_L9NR8_R2.fastq.gz"  
## [145] "data/00_raw_reads/Plate2_1351_G07_Blank2_10471351_L9NR8_R2.fastq.gz"    "data/00_raw_reads/Plate2_1351_H01_AP_D110_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_H02_AP_D118_10471351_L9NR8_R2.fastq.gz"  
## [148] "data/00_raw_reads/Plate2_1351_H03_AP_D126_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_H04_AP_D134_10471351_L9NR8_R2.fastq.gz"   "data/00_raw_reads/Plate2_1351_H05_AP_D142_10471351_L9NR8_R2.fastq.gz"  
## [151] "data/00_raw_reads/Plate2_1351_H06_AP_D150_10471351_L9NR8_R2.fastq.gz"
```

```r
# Get clean names
sample_names <- data.frame(bname = basename(forward_reads),
            rname = basename(reverse_reads)) %>% # Extract the basename (file name)
  separate_wider_delim(bname, delim = "_", names_sep = "_", too_few = "align_start") %>% # Split up by the underscores
  transmute(sample_name = ifelse(bname_1 == "Plate1",
                              paste(bname_3, bname_4, sep = "_"),
                              paste(bname_4, bname_5, sep = "_")))%>% # Plate 1 and Plate 2 had a diff. number of underscores until the DNA idea - hence the ifelse
  pull()

# 3. Prep filepaths for filtered files.
# Create a variable holding file names for the Forward and Reverse filtered reads 

filtered_forward_reads <- file.path("data", "filtered", paste0(sample_names, "_R1_filtered.fastq.gz"))
filtered_reverse_reads <- file.path("data", "filtered", paste0(sample_names, "_R2_filtered.fastq.gz"))

# An small view into the quality of some samples forward and reverse
set.seed(031491)
random_plots <- sample(1:151, size = 20)

plotQualityProfile(forward_reads[random_plots])
```

<img src="../figures/01_DADA2_Processing/read-and-quality-1.png" style="display: block; margin: auto;" />

```r
plotQualityProfile(reverse_reads[random_plots])
```

<img src="../figures/01_DADA2_Processing/read-and-quality-2.png" style="display: block; margin: auto;" />

forward reads: 
- Pretty clean, in general
- Slight lower quality at beginning - maybe 15 bases in it gets pretty clean?
- The only samples that have a significant drop in quality halfway in are the Blanks/Negative Controls and AP_D33.
- The Mock community was also lower quality than all of our samples (annoying and a little strange)
- Read counts look good

reverse reads: 
- Worse quality overall than forward (expected)
- Most have some reads that start to drop ~220

aggregate plot:
- Overall pretty happy. 
- Forward and reverse get up to high quality within 20 bp
- Tbh it doesn't look like forward need to be right-trimmed, but will do a bit
- Reverse gonna right-trim to 220bp

# Filter and Trim

```r
set.seed(031491)
# Major filtering step!
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads,
                              truncLen = c(240,220), trimLeft = c(20,20),
                              maxN = 0, maxEE = c(1,1), truncQ = 2, 
                              rm.phix = TRUE, compress = TRUE, 
                              multithread = 10,
                              orient.fwd = "GTG")

# Plot the quality of trimmed reads! 
plotQualityProfile(filtered_forward_reads[random_plots])
```

<img src="../figures/01_DADA2_Processing/filter-trim-1.png" style="display: block; margin: auto;" />

```r
plotQualityProfile(filtered_reverse_reads[random_plots])
```

<img src="../figures/01_DADA2_Processing/filter-trim-2.png" style="display: block; margin: auto;" />

```r
#lines stay above 30 on all plots, median reads.in around 60k
filter_df <- as.data.frame(filtered_out)
cat("Median ",median(filter_df$reads.in), " reads in")
```

```
## Median  70037  reads in
```

```r
cat("Median ", median(filter_df$reads.out), "reads out")
```

```
## Median  49405 reads out
```

```r
cat("Median", median(filter_df$reads.in) - median(filter_df$reads.out), " reads removed")
```

```
## Median 20632  reads removed
```

Overall, dropped about 20,000 reads from each sample.

# Generate an error model 

```r
# Learn errors
err_forward_reads <- learnErrors(filtered_forward_reads, multithread = TRUE) 
```

```
## 100769460 total bases in 458043 reads from 11 samples will be used for learning the error rates.
```

```r
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread = TRUE)
```

```
## 109653400 total bases in 548267 reads from 12 samples will be used for learning the error rates.
```

```r
# Plot the errors
# the estimated error rates (black line) are a good fit to the observed rates (points), and the error rates drop with increased quality as expected
plotErrors(err_forward_reads, nominalQ = TRUE) 
```

<img src="../figures/01_DADA2_Processing/learn-errors-1.png" style="display: block; margin: auto;" />

```r
plotErrors(err_reverse_reads, nominalQ = TRUE)
```

<img src="../figures/01_DADA2_Processing/learn-errors-2.png" style="display: block; margin: auto;" />

Predicted error hews fairly close to observed error, which is good.

# Inferring ASVs on the forward and reverse sequences 


```r
# run dada2 on the forward sequences
dada_forward <- dada(filtered_forward_reads, err = err_forward_reads, multithread = TRUE)
```

```
## Sample 1 - 48217 reads in 8321 unique sequences.
## Sample 2 - 28469 reads in 4750 unique sequences.
## Sample 3 - 31754 reads in 5395 unique sequences.
## Sample 4 - 26624 reads in 4502 unique sequences.
## Sample 5 - 11515 reads in 2662 unique sequences.
## Sample 6 - 36082 reads in 5596 unique sequences.
## Sample 7 - 39115 reads in 6272 unique sequences.
## Sample 8 - 33406 reads in 4975 unique sequences.
## Sample 9 - 63256 reads in 8696 unique sequences.
## Sample 10 - 64171 reads in 9420 unique sequences.
## Sample 11 - 75434 reads in 9547 unique sequences.
## Sample 12 - 90224 reads in 12488 unique sequences.
## Sample 13 - 45932 reads in 7967 unique sequences.
## Sample 14 - 31211 reads in 5490 unique sequences.
## Sample 15 - 40948 reads in 6738 unique sequences.
## Sample 16 - 35308 reads in 5643 unique sequences.
## Sample 17 - 21382 reads in 4235 unique sequences.
## Sample 18 - 35753 reads in 5434 unique sequences.
## Sample 19 - 46557 reads in 7411 unique sequences.
## Sample 20 - 38005 reads in 6194 unique sequences.
## Sample 21 - 52096 reads in 6875 unique sequences.
## Sample 22 - 53240 reads in 7652 unique sequences.
## Sample 23 - 48606 reads in 6294 unique sequences.
## Sample 24 - 63832 reads in 10380 unique sequences.
## Sample 25 - 37113 reads in 5755 unique sequences.
## Sample 26 - 31422 reads in 4867 unique sequences.
## Sample 27 - 35794 reads in 5365 unique sequences.
## Sample 28 - 17992 reads in 3620 unique sequences.
## Sample 29 - 31745 reads in 4600 unique sequences.
## Sample 30 - 36426 reads in 5716 unique sequences.
## Sample 31 - 44744 reads in 6206 unique sequences.
## Sample 32 - 41653 reads in 6660 unique sequences.
## Sample 33 - 34166 reads in 4919 unique sequences.
## Sample 34 - 48288 reads in 6766 unique sequences.
## Sample 35 - 50110 reads in 6805 unique sequences.
## Sample 36 - 52011 reads in 7606 unique sequences.
## Sample 37 - 42465 reads in 6200 unique sequences.
## Sample 38 - 32523 reads in 4867 unique sequences.
## Sample 39 - 29796 reads in 4656 unique sequences.
## Sample 40 - 31946 reads in 5274 unique sequences.
## Sample 41 - 30514 reads in 4768 unique sequences.
## Sample 42 - 32653 reads in 5444 unique sequences.
## Sample 43 - 44551 reads in 6459 unique sequences.
## Sample 44 - 71 reads in 42 unique sequences.
## Sample 45 - 40333 reads in 6574 unique sequences.
## Sample 46 - 46683 reads in 6280 unique sequences.
## Sample 47 - 24968 reads in 5314 unique sequences.
## Sample 48 - 60020 reads in 9386 unique sequences.
## Sample 49 - 64062 reads in 8785 unique sequences.
## Sample 50 - 36920 reads in 5762 unique sequences.
## Sample 51 - 51463 reads in 6818 unique sequences.
## Sample 52 - 36502 reads in 5669 unique sequences.
## Sample 53 - 66001 reads in 8646 unique sequences.
## Sample 54 - 57261 reads in 8692 unique sequences.
## Sample 55 - 59245 reads in 7282 unique sequences.
## Sample 56 - 481 reads in 207 unique sequences.
## Sample 57 - 53535 reads in 8502 unique sequences.
## Sample 58 - 59175 reads in 6887 unique sequences.
## Sample 59 - 53969 reads in 8235 unique sequences.
## Sample 60 - 62115 reads in 8807 unique sequences.
## Sample 61 - 49405 reads in 7263 unique sequences.
## Sample 62 - 30748 reads in 5153 unique sequences.
## Sample 63 - 50942 reads in 8450 unique sequences.
## Sample 64 - 37547 reads in 6294 unique sequences.
## Sample 65 - 40970 reads in 6534 unique sequences.
## Sample 66 - 47023 reads in 7781 unique sequences.
## Sample 67 - 39663 reads in 6068 unique sequences.
## Sample 68 - 39517 reads in 6453 unique sequences.
## Sample 69 - 39819 reads in 6081 unique sequences.
## Sample 70 - 49004 reads in 7289 unique sequences.
## Sample 71 - 31812 reads in 5045 unique sequences.
## Sample 72 - 64139 reads in 13632 unique sequences.
## Sample 73 - 59274 reads in 8140 unique sequences.
## Sample 74 - 25833 reads in 4770 unique sequences.
## Sample 75 - 40865 reads in 6505 unique sequences.
## Sample 76 - 38541 reads in 6102 unique sequences.
## Sample 77 - 63708 reads in 8425 unique sequences.
## Sample 78 - 47856 reads in 7330 unique sequences.
## Sample 79 - 47549 reads in 6958 unique sequences.
## Sample 80 - 48030 reads in 7269 unique sequences.
## Sample 81 - 55035 reads in 7347 unique sequences.
## Sample 82 - 64747 reads in 8337 unique sequences.
## Sample 83 - 38017 reads in 5519 unique sequences.
## Sample 84 - 85082 reads in 15586 unique sequences.
## Sample 85 - 52439 reads in 7487 unique sequences.
## Sample 86 - 14375 reads in 2947 unique sequences.
## Sample 87 - 39359 reads in 6497 unique sequences.
## Sample 88 - 31545 reads in 5742 unique sequences.
## Sample 89 - 56397 reads in 7989 unique sequences.
## Sample 90 - 59990 reads in 8933 unique sequences.
## Sample 91 - 61815 reads in 8802 unique sequences.
## Sample 92 - 46650 reads in 6843 unique sequences.
## Sample 93 - 64826 reads in 9459 unique sequences.
## Sample 94 - 66077 reads in 8875 unique sequences.
## Sample 95 - 50118 reads in 8111 unique sequences.
## Sample 96 - 50 reads in 30 unique sequences.
## Sample 97 - 67981 reads in 10710 unique sequences.
## Sample 98 - 62672 reads in 8438 unique sequences.
## Sample 99 - 63098 reads in 7474 unique sequences.
## Sample 100 - 67712 reads in 8914 unique sequences.
## Sample 101 - 53514 reads in 7090 unique sequences.
## Sample 102 - 52786 reads in 8054 unique sequences.
## Sample 103 - 48591 reads in 6547 unique sequences.
## Sample 104 - 77235 reads in 12803 unique sequences.
## Sample 105 - 77233 reads in 10940 unique sequences.
## Sample 106 - 63232 reads in 7540 unique sequences.
## Sample 107 - 76281 reads in 10205 unique sequences.
## Sample 108 - 67302 reads in 8252 unique sequences.
## Sample 109 - 55348 reads in 9705 unique sequences.
## Sample 110 - 54247 reads in 7217 unique sequences.
## Sample 111 - 54678 reads in 7058 unique sequences.
## Sample 112 - 57973 reads in 6902 unique sequences.
## Sample 113 - 54594 reads in 6765 unique sequences.
## Sample 114 - 40049 reads in 7000 unique sequences.
## Sample 115 - 65994 reads in 7835 unique sequences.
## Sample 116 - 53527 reads in 8292 unique sequences.
## Sample 117 - 4 reads in 4 unique sequences.
## Sample 118 - 74505 reads in 8460 unique sequences.
## Sample 119 - 56489 reads in 7201 unique sequences.
## Sample 120 - 56688 reads in 6817 unique sequences.
## Sample 121 - 62259 reads in 9287 unique sequences.
## Sample 122 - 68860 reads in 8625 unique sequences.
## Sample 123 - 49170 reads in 7675 unique sequences.
## Sample 124 - 21 reads in 17 unique sequences.
## Sample 125 - 69851 reads in 10184 unique sequences.
## Sample 126 - 84626 reads in 9684 unique sequences.
## Sample 127 - 47399 reads in 6947 unique sequences.
## Sample 128 - 65358 reads in 8726 unique sequences.
## Sample 129 - 52672 reads in 6873 unique sequences.
## Sample 130 - 51587 reads in 6893 unique sequences.
## Sample 131 - 53153 reads in 6091 unique sequences.
## Sample 132 - 68799 reads in 10934 unique sequences.
## Sample 133 - 66528 reads in 7879 unique sequences.
## Sample 134 - 50314 reads in 7145 unique sequences.
## Sample 135 - 50221 reads in 7229 unique sequences.
## Sample 136 - 56045 reads in 7377 unique sequences.
## Sample 137 - 36732 reads in 5489 unique sequences.
## Sample 138 - 11 reads in 9 unique sequences.
## Sample 139 - 51243 reads in 8037 unique sequences.
## Sample 140 - 55336 reads in 7478 unique sequences.
## Sample 141 - 49304 reads in 7424 unique sequences.
## Sample 142 - 42559 reads in 5732 unique sequences.
## Sample 143 - 55527 reads in 8050 unique sequences.
## Sample 144 - 42904 reads in 6810 unique sequences.
## Sample 145 - 5 reads in 4 unique sequences.
## Sample 146 - 64277 reads in 9223 unique sequences.
## Sample 147 - 60605 reads in 7942 unique sequences.
## Sample 148 - 51170 reads in 7245 unique sequences.
## Sample 149 - 57969 reads in 7504 unique sequences.
## Sample 150 - 91587 reads in 11588 unique sequences.
## Sample 151 - 48527 reads in 8276 unique sequences.
```

```r
# run dada2 on the reverse sequences 
dada_reverse <- dada(filtered_reverse_reads, err = err_reverse_reads, multithread = TRUE)
```

```
## Sample 1 - 48217 reads in 8710 unique sequences.
## Sample 2 - 28469 reads in 5085 unique sequences.
## Sample 3 - 31754 reads in 5573 unique sequences.
## Sample 4 - 26624 reads in 4582 unique sequences.
## Sample 5 - 11515 reads in 2844 unique sequences.
## Sample 6 - 36082 reads in 6418 unique sequences.
## Sample 7 - 39115 reads in 6748 unique sequences.
## Sample 8 - 33406 reads in 5776 unique sequences.
## Sample 9 - 63256 reads in 9260 unique sequences.
## Sample 10 - 64171 reads in 9874 unique sequences.
## Sample 11 - 75434 reads in 10594 unique sequences.
## Sample 12 - 90224 reads in 13955 unique sequences.
## Sample 13 - 45932 reads in 7981 unique sequences.
## Sample 14 - 31211 reads in 5451 unique sequences.
## Sample 15 - 40948 reads in 6269 unique sequences.
## Sample 16 - 35308 reads in 5374 unique sequences.
## Sample 17 - 21382 reads in 4377 unique sequences.
## Sample 18 - 35753 reads in 6158 unique sequences.
## Sample 19 - 46557 reads in 7632 unique sequences.
## Sample 20 - 38005 reads in 6763 unique sequences.
## Sample 21 - 52096 reads in 7216 unique sequences.
## Sample 22 - 53240 reads in 7820 unique sequences.
## Sample 23 - 48606 reads in 6943 unique sequences.
## Sample 24 - 63832 reads in 11062 unique sequences.
## Sample 25 - 37113 reads in 6136 unique sequences.
## Sample 26 - 31422 reads in 5153 unique sequences.
## Sample 27 - 35794 reads in 5864 unique sequences.
## Sample 28 - 17992 reads in 3789 unique sequences.
## Sample 29 - 31745 reads in 4994 unique sequences.
## Sample 30 - 36426 reads in 6886 unique sequences.
## Sample 31 - 44744 reads in 6680 unique sequences.
## Sample 32 - 41653 reads in 7185 unique sequences.
## Sample 33 - 34166 reads in 5092 unique sequences.
## Sample 34 - 48288 reads in 7518 unique sequences.
## Sample 35 - 50110 reads in 7574 unique sequences.
## Sample 36 - 52011 reads in 10243 unique sequences.
## Sample 37 - 42465 reads in 6485 unique sequences.
## Sample 38 - 32523 reads in 5019 unique sequences.
## Sample 39 - 29796 reads in 4344 unique sequences.
## Sample 40 - 31946 reads in 5521 unique sequences.
## Sample 41 - 30514 reads in 5316 unique sequences.
## Sample 42 - 32653 reads in 6235 unique sequences.
## Sample 43 - 44551 reads in 6390 unique sequences.
## Sample 44 - 71 reads in 45 unique sequences.
## Sample 45 - 40333 reads in 6842 unique sequences.
## Sample 46 - 46683 reads in 6401 unique sequences.
## Sample 47 - 24968 reads in 4395 unique sequences.
## Sample 48 - 60020 reads in 10229 unique sequences.
## Sample 49 - 64062 reads in 9281 unique sequences.
## Sample 50 - 36920 reads in 6462 unique sequences.
## Sample 51 - 51463 reads in 8428 unique sequences.
## Sample 52 - 36502 reads in 6136 unique sequences.
## Sample 53 - 66001 reads in 10471 unique sequences.
## Sample 54 - 57261 reads in 10732 unique sequences.
## Sample 55 - 59245 reads in 7752 unique sequences.
## Sample 56 - 481 reads in 235 unique sequences.
## Sample 57 - 53535 reads in 9389 unique sequences.
## Sample 58 - 59175 reads in 7353 unique sequences.
## Sample 59 - 53969 reads in 9401 unique sequences.
## Sample 60 - 62115 reads in 12233 unique sequences.
## Sample 61 - 49405 reads in 6895 unique sequences.
## Sample 62 - 30748 reads in 4953 unique sequences.
## Sample 63 - 50942 reads in 7601 unique sequences.
## Sample 64 - 37547 reads in 5582 unique sequences.
## Sample 65 - 40970 reads in 5868 unique sequences.
## Sample 66 - 47023 reads in 7773 unique sequences.
## Sample 67 - 39663 reads in 5463 unique sequences.
## Sample 68 - 39517 reads in 6481 unique sequences.
## Sample 69 - 39819 reads in 5656 unique sequences.
## Sample 70 - 49004 reads in 6710 unique sequences.
## Sample 71 - 31812 reads in 4730 unique sequences.
## Sample 72 - 64139 reads in 13384 unique sequences.
## Sample 73 - 59274 reads in 8402 unique sequences.
## Sample 74 - 25833 reads in 4972 unique sequences.
## Sample 75 - 40865 reads in 7398 unique sequences.
## Sample 76 - 38541 reads in 6346 unique sequences.
## Sample 77 - 63708 reads in 8603 unique sequences.
## Sample 78 - 47856 reads in 8420 unique sequences.
## Sample 79 - 47549 reads in 7380 unique sequences.
## Sample 80 - 48030 reads in 8130 unique sequences.
## Sample 81 - 55035 reads in 7571 unique sequences.
## Sample 82 - 64747 reads in 8872 unique sequences.
## Sample 83 - 38017 reads in 5698 unique sequences.
## Sample 84 - 85082 reads in 16279 unique sequences.
## Sample 85 - 52439 reads in 8021 unique sequences.
## Sample 86 - 14375 reads in 4302 unique sequences.
## Sample 87 - 39359 reads in 6856 unique sequences.
## Sample 88 - 31545 reads in 5287 unique sequences.
## Sample 89 - 56397 reads in 8285 unique sequences.
## Sample 90 - 59990 reads in 10353 unique sequences.
## Sample 91 - 61815 reads in 9534 unique sequences.
## Sample 92 - 46650 reads in 7815 unique sequences.
## Sample 93 - 64826 reads in 10329 unique sequences.
## Sample 94 - 66077 reads in 9206 unique sequences.
## Sample 95 - 50118 reads in 9183 unique sequences.
## Sample 96 - 50 reads in 34 unique sequences.
## Sample 97 - 67981 reads in 14951 unique sequences.
## Sample 98 - 62672 reads in 9271 unique sequences.
## Sample 99 - 63098 reads in 10163 unique sequences.
## Sample 100 - 67712 reads in 9430 unique sequences.
## Sample 101 - 53514 reads in 7666 unique sequences.
## Sample 102 - 52786 reads in 11717 unique sequences.
## Sample 103 - 48591 reads in 6886 unique sequences.
## Sample 104 - 77235 reads in 12341 unique sequences.
## Sample 105 - 77233 reads in 11517 unique sequences.
## Sample 106 - 63232 reads in 9581 unique sequences.
## Sample 107 - 76281 reads in 10239 unique sequences.
## Sample 108 - 67302 reads in 8750 unique sequences.
## Sample 109 - 55348 reads in 9270 unique sequences.
## Sample 110 - 54247 reads in 6898 unique sequences.
## Sample 111 - 54678 reads in 8031 unique sequences.
## Sample 112 - 57973 reads in 8436 unique sequences.
## Sample 113 - 54594 reads in 9330 unique sequences.
## Sample 114 - 40049 reads in 6557 unique sequences.
## Sample 115 - 65994 reads in 10459 unique sequences.
## Sample 116 - 53527 reads in 9272 unique sequences.
## Sample 117 - 4 reads in 4 unique sequences.
## Sample 118 - 74505 reads in 10802 unique sequences.
## Sample 119 - 56489 reads in 9036 unique sequences.
## Sample 120 - 56688 reads in 9953 unique sequences.
## Sample 121 - 62259 reads in 10949 unique sequences.
## Sample 122 - 68860 reads in 10715 unique sequences.
## Sample 123 - 49170 reads in 9240 unique sequences.
## Sample 124 - 21 reads in 17 unique sequences.
## Sample 125 - 69851 reads in 11744 unique sequences.
## Sample 126 - 84626 reads in 11862 unique sequences.
## Sample 127 - 47399 reads in 8876 unique sequences.
## Sample 128 - 65358 reads in 9466 unique sequences.
## Sample 129 - 52672 reads in 7850 unique sequences.
## Sample 130 - 51587 reads in 7765 unique sequences.
## Sample 131 - 53153 reads in 4863 unique sequences.
## Sample 132 - 68799 reads in 12191 unique sequences.
## Sample 133 - 66528 reads in 9705 unique sequences.
## Sample 134 - 50314 reads in 9341 unique sequences.
## Sample 135 - 50221 reads in 7672 unique sequences.
## Sample 136 - 56045 reads in 8152 unique sequences.
## Sample 137 - 36732 reads in 6009 unique sequences.
## Sample 138 - 11 reads in 9 unique sequences.
## Sample 139 - 51243 reads in 9160 unique sequences.
## Sample 140 - 55336 reads in 8703 unique sequences.
## Sample 141 - 49304 reads in 8811 unique sequences.
## Sample 142 - 42559 reads in 5883 unique sequences.
## Sample 143 - 55527 reads in 8139 unique sequences.
## Sample 144 - 42904 reads in 7334 unique sequences.
## Sample 145 - 5 reads in 4 unique sequences.
## Sample 146 - 64277 reads in 9788 unique sequences.
## Sample 147 - 60605 reads in 9016 unique sequences.
## Sample 148 - 51170 reads in 9088 unique sequences.
## Sample 149 - 57969 reads in 7857 unique sequences.
## Sample 150 - 91587 reads in 12493 unique sequences.
## Sample 151 - 48527 reads in 8885 unique sequences.
```


# Merge forward and reverse ASVs 

```r
# Merge the forward ASVs and the reverse ASVs
merged_amplicons <- mergePairs(dada_forward, filtered_forward_reads, 
                               dada_reverse, filtered_reverse_reads,
                               verbose = TRUE)
# Evaluate the output 
length(merged_amplicons) #151
```

```
## [1] 151
```

```r
names(merged_amplicons)
```

```
##   [1] "AP_D1_R1_filtered.fastq.gz"           "AP_D9_R1_filtered.fastq.gz"           "AP_D17_R1_filtered.fastq.gz"          "AP_D25_R1_filtered.fastq.gz"          "AP_D33_R1_filtered.fastq.gz"          "AP_D41_R1_filtered.fastq.gz"         
##   [7] "AP_D49_R1_filtered.fastq.gz"          "AP_D58_R1_filtered.fastq.gz"          "AP_D72_R1_filtered.fastq.gz"          "AP_D80_R1_filtered.fastq.gz"          "AP_D88_R1_filtered.fastq.gz"          "AP_D96_R1_filtered.fastq.gz"         
##  [13] "AP_D2_R1_filtered.fastq.gz"           "AP_D10_R1_filtered.fastq.gz"          "AP_D18_R1_filtered.fastq.gz"          "AP_D26_R1_filtered.fastq.gz"          "AP_D34_R1_filtered.fastq.gz"          "AP_D42_R1_filtered.fastq.gz"         
##  [19] "AP_D50_R1_filtered.fastq.gz"          "AP_D59_R1_filtered.fastq.gz"          "AP_D73_R1_filtered.fastq.gz"          "AP_D81_R1_filtered.fastq.gz"          "AP_D89_R1_filtered.fastq.gz"          "AP_D97_R1_filtered.fastq.gz"         
##  [25] "AP_D3_R1_filtered.fastq.gz"           "AP_D11_R1_filtered.fastq.gz"          "AP_D19_R1_filtered.fastq.gz"          "AP_D27_R1_filtered.fastq.gz"          "AP_D35_R1_filtered.fastq.gz"          "AP_D43_R1_filtered.fastq.gz"         
##  [31] "AP_D51_R1_filtered.fastq.gz"          "AP_D60_R1_filtered.fastq.gz"          "AP_D74_R1_filtered.fastq.gz"          "AP_D82_R1_filtered.fastq.gz"          "AP_D90_R1_filtered.fastq.gz"          "AP_D98_R1_filtered.fastq.gz"         
##  [37] "AP_D4_R1_filtered.fastq.gz"           "AP_D12_R1_filtered.fastq.gz"          "AP_D20_R1_filtered.fastq.gz"          "AP_D28_R1_filtered.fastq.gz"          "AP_D36_R1_filtered.fastq.gz"          "AP_D44_R1_filtered.fastq.gz"         
##  [43] "AP_D52_R1_filtered.fastq.gz"          "AP_D61_R1_filtered.fastq.gz"          "AP_D75_R1_filtered.fastq.gz"          "AP_D83_R1_filtered.fastq.gz"          "AP_D91_R1_filtered.fastq.gz"          "AP_D99_R1_filtered.fastq.gz"         
##  [49] "AP_D5_R1_filtered.fastq.gz"           "AP_D13_R1_filtered.fastq.gz"          "AP_D21_R1_filtered.fastq.gz"          "AP_D29_R1_filtered.fastq.gz"          "AP_D37_R1_filtered.fastq.gz"          "AP_D45_R1_filtered.fastq.gz"         
##  [55] "AP_D53_R1_filtered.fastq.gz"          "AP_D62_R1_filtered.fastq.gz"          "AP_D76_R1_filtered.fastq.gz"          "AP_D84_R1_filtered.fastq.gz"          "AP_D92_R1_filtered.fastq.gz"          "AP_D100_R1_filtered.fastq.gz"        
##  [61] "AP_D6_R1_filtered.fastq.gz"           "AP_D14_R1_filtered.fastq.gz"          "AP_D22_R1_filtered.fastq.gz"          "AP_D30_R1_filtered.fastq.gz"          "AP_D38_R1_filtered.fastq.gz"          "AP_D46_R1_filtered.fastq.gz"         
##  [67] "AP_D54_R1_filtered.fastq.gz"          "AP_D69_R1_filtered.fastq.gz"          "AP_D77_R1_filtered.fastq.gz"          "AP_D85_R1_filtered.fastq.gz"          "AP_D93_R1_filtered.fastq.gz"          "AP_D101_R1_filtered.fastq.gz"        
##  [73] "AP_D7_R1_filtered.fastq.gz"           "AP_D15_R1_filtered.fastq.gz"          "AP_D23_R1_filtered.fastq.gz"          "AP_D31_R1_filtered.fastq.gz"          "AP_D39_R1_filtered.fastq.gz"          "AP_D47_R1_filtered.fastq.gz"         
##  [79] "AP_D55_R1_filtered.fastq.gz"          "AP_D70_R1_filtered.fastq.gz"          "AP_D78_R1_filtered.fastq.gz"          "AP_D86_R1_filtered.fastq.gz"          "AP_D94_R1_filtered.fastq.gz"          "AP_D102_R1_filtered.fastq.gz"        
##  [85] "AP_D8_R1_filtered.fastq.gz"           "AP_D16_R1_filtered.fastq.gz"          "AP_D24_R1_filtered.fastq.gz"          "AP_D32_R1_filtered.fastq.gz"          "AP_D40_R1_filtered.fastq.gz"          "AP_D48_R1_filtered.fastq.gz"         
##  [91] "AP_D57_R1_filtered.fastq.gz"          "AP_D71_R1_filtered.fastq.gz"          "AP_D79_R1_filtered.fastq.gz"          "AP_D87_R1_filtered.fastq.gz"          "AP_D95_R1_filtered.fastq.gz"          "PCR_Blanks_R1_filtered.fastq.gz"     
##  [97] "AP_D103_R1_filtered.fastq.gz"         "AP_D111_R1_filtered.fastq.gz"         "AP_D119_R1_filtered.fastq.gz"         "AP_D127_R1_filtered.fastq.gz"         "AP_D135_R1_filtered.fastq.gz"         "AP_D143_R1_filtered.fastq.gz"        
## [103] "AP_D151_R1_filtered.fastq.gz"         "AP_D104_R1_filtered.fastq.gz"         "AP_D112_R1_filtered.fastq.gz"         "AP_D120_R1_filtered.fastq.gz"         "AP_D128_R1_filtered.fastq.gz"         "AP_D136_R1_filtered.fastq.gz"        
## [109] "AP_D144_R1_filtered.fastq.gz"         "AP_D152_R1_filtered.fastq.gz"         "AP_D105_R1_filtered.fastq.gz"         "AP_D113_R1_filtered.fastq.gz"         "AP_D121_R1_filtered.fastq.gz"         "AP_D129_R1_filtered.fastq.gz"        
## [115] "AP_D137_R1_filtered.fastq.gz"         "AP_D145_R1_filtered.fastq.gz"         "AP_D153_R1_filtered.fastq.gz"         "AP_D106_R1_filtered.fastq.gz"         "AP_D114_R1_filtered.fastq.gz"         "AP_D122_R1_filtered.fastq.gz"        
## [121] "AP_D130_R1_filtered.fastq.gz"         "AP_D138_R1_filtered.fastq.gz"         "AP_D146_R1_filtered.fastq.gz"         "AP_D154_R1_filtered.fastq.gz"         "AP_D107_R1_filtered.fastq.gz"         "AP_D115_R1_filtered.fastq.gz"        
## [127] "AP_D123_R1_filtered.fastq.gz"         "AP_D131_R1_filtered.fastq.gz"         "AP_D139_R1_filtered.fastq.gz"         "AP_D147_R1_filtered.fastq.gz"         "Zymo_Mock_R1_filtered.fastq.gz"       "AP_D108_R1_filtered.fastq.gz"        
## [133] "AP_D116_R1_filtered.fastq.gz"         "AP_D124_R1_filtered.fastq.gz"         "AP_D132_R1_filtered.fastq.gz"         "AP_D140_R1_filtered.fastq.gz"         "AP_D148_R1_filtered.fastq.gz"         "Blank1_10471351_R1_filtered.fastq.gz"
## [139] "AP_D109_R1_filtered.fastq.gz"         "AP_D117_R1_filtered.fastq.gz"         "AP_D125_R1_filtered.fastq.gz"         "AP_D133_R1_filtered.fastq.gz"         "AP_D141_R1_filtered.fastq.gz"         "AP_D149_R1_filtered.fastq.gz"        
## [145] "Blank2_10471351_R1_filtered.fastq.gz" "AP_D110_R1_filtered.fastq.gz"         "AP_D118_R1_filtered.fastq.gz"         "AP_D126_R1_filtered.fastq.gz"         "AP_D134_R1_filtered.fastq.gz"         "AP_D142_R1_filtered.fastq.gz"        
## [151] "AP_D150_R1_filtered.fastq.gz"
```


# Generate a count table! 


```r
seqtab <- makeSequenceTable(merged_amplicons)

dim(seqtab) # 149 20329
```

```
## [1]   151 19879
```

```r
# Inspect the distribution of sequence lengths of all ASVs in dataset 
# all sequences fall between 220 and 409
data.frame(Seq_Length = nchar(getSequences(seqtab))) %>%
  ggplot(aes(x = Seq_Length )) + 
  geom_histogram()
```

<img src="../figures/01_DADA2_Processing/gen-countTable-seqTab-1.png" style="display: block; margin: auto;" />

```r
# Inspect the distribution of sequence lengths of all ASVs in dataset 
table(nchar(getSequences(seqtab))) %>%sort()
```

```
## 
##  355  398  388  358  389  391  392  368  323  343  325  369  390  399  341  353  381  385  387  397  333  349  361  404  405  340  374  380  382  306  352  370  375  319  386  304  328  335  337  342  346  364  365  366  384  400  315  329  334  330 
##    8   10   13   16   16   16   16   17   18   18   19   19   19   19   20   20   20   20   20   20   21   21   21   21   21   22   22   22   22   23   23   23   23   24   24   25   25   25   25   25   25   25   25   25   25   25   26   26   26   27 
##  336  393  295  296  301  309  310  313  344  347  402  406  408  292  322  339  373  383  403  305  326  281  324  331  289  302  307  312  317  288  363  367  396  297  376  407  287  321  351  354  279  300  318  283  294  298  311  332  348  350 
##   27   27   28   28   28   28   28   28   28   28   28   28   28   29   29   29   29   29   29   30   30   31   31   31   32   32   32   32   32   33   33   33   33   34   34   34   35   35   35   35   36   36   36   37   37   37   37   37   37   37 
##  308  314  316  356  372  377  378  395  401  274  345  265  299  303  224  327  359  379  275  236  239  262  338  271  394  371  244  278  284  242  286  293  360  234  245  280  320  268  261  276  277  223  260  221  238  290  282  273  235  357 
##   38   38   38   38   38   38   38   38   38   39   39   40   40   40   41   41   41   41   42   43   43   43   43   44   44   45   46   46   46   48   49   49   49   50   50   50   50   51   52   52   52   53   54   55   55   55   56   58   59   59 
##  228  243  285  266  237  257  272  230  233  270  362  229  263  267  246  247  232  240  259  264  248  225  227  231  291  269  226  222  258  241  256  249  255  250  254  251  220  253  252 
##   61   61   62   63   66   67   67   68   68   68   68   70   70   74   75   76   78   78   79   81   85   86   86   86   86   89   92   95   99  102  108  135  164  190  312  520  794 1224 9184
```

Okay so it looks like a majority of ASVs are length 252 and 253


# Check & Remove for Chimeras (Bimeras)

Check size after bimera removal - trimming step
Add table() command

Also look at tweaking trimming to see if that cleans up the Mock


```r
# Identify and remove chimeras 
seqtab_nochim <- removeBimeraDenovo(seqtab, verbose = TRUE)
# Identified 1663 bimeras out of 21502 input sequences.


data.frame(Seq_Length_NoChim = nchar(getSequences(seqtab_nochim))) %>%
  ggplot(aes(x = Seq_Length_NoChim )) + 
  geom_histogram()
```

<img src="../figures/01_DADA2_Processing/check-chimeras-1.png" style="display: block; margin: auto;" />

```r
asvs <- dim(seqtab_nochim)[2]

# What proportion of counts were removed? 
chim_check <- sum(seqtab_nochim)/sum(seqtab) #  0.9835776
frac_removed <- (1-chim_check)*100
frac_removed # 1.642241
```

```
## [1] 1.653096
```

Chimeras represented 1.653096 percent of the data, with parameters trimLeft(20,20).

Dataset includes 18200 ASVs.

# Size Selection

I also want to get rid of ASVs that too big.


```r
asv_keeps <- nchar(getSequences(seqtab_nochim)) %in% c(252,253) # Find which sequences have length 252 or 253

seqtab_nc_len <- seqtab_nochim[,asv_keeps] # Remove other sequences

data.frame(Seq_Length_NoChim_Len = nchar(getSequences(seqtab_nc_len))) %>%
  ggplot(aes(x = Seq_Length_NoChim_Len )) + 
  geom_histogram()
```

<img src="../figures/01_DADA2_Processing/size selection-1.png" style="display: block; margin: auto;" />

We removed all ASVs that weren't size 252 or 253.

# Track the sequences through the pipeline 


```r
# create a little function to identify number seqs 
getN <- function(x) sum(getUniques(x))

# Make the table to track the seqs 
track <- cbind(filtered_out, 
               sapply(dada_forward, getN),
               sapply(dada_reverse, getN),
               sapply(merged_amplicons, getN),
               rowSums(seqtab_nc_len))

head(track)
```

```
##                                              reads.in reads.out                        
## Plate1_A01_AP_D1_10471351_L9NR8_R1.fastq.gz     69027     48217 46611 46651 45363 44190
## Plate1_A02_AP_D9_10471351_L9NR8_R1.fastq.gz     41411     28469 27641 27785 26654 25778
## Plate1_A03_AP_D17_10471351_L9NR8_R1.fastq.gz    44873     31754 30810 30791 29370 28214
## Plate1_A04_AP_D25_10471351_L9NR8_R1.fastq.gz    39900     26624 25751 25768 24746 24078
## Plate1_A05_AP_D33_10471351_L9NR8_R1.fastq.gz    22430     11515 10839 10784 10106  9653
## Plate1_A06_AP_D41_10471351_L9NR8_R1.fastq.gz    53787     36082 35262 35145 34301 33014
```

```r
# Change column names 
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochim")
rownames(track) <- sample_names

# Generate a plot to track the reads through our DADA2 pipeline
track %>%
  # make it a dataframe
  as.data.frame() %>%
  rownames_to_column(var = "sample_name") %>%
  pivot_longer(input:nochim, names_to = "read_type", values_to = "num_reads") %>%
  #left_join(metadata, by = "sample_name") %>% 
  mutate(read_type = fct_relevel(read_type, 
                                 "input", "filtered", "denoisedF", "denoisedR", "merged", "nochim"),
         is_blank = sample_name%in%c("PCR_Blanks","AP_D61","AP_D62","AP_D153","AP_D154")) %>%
  ggplot(aes(x = read_type, y = num_reads, fill = read_type, color = is_blank)) + 
  #facet_grid(~strata) + 
  geom_line(aes(group = sample_name)) + 
  geom_point(shape = 21, size = 3, alpha = 0.8, color = "grey") + 
  scale_fill_brewer(palette = "Spectral") + 
  theme_bw() + 
  labs(x = "Filtering Step", y = "Number of Sequences") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
```

<img src="../figures/01_DADA2_Processing/track-seqs-1.png" style="display: block; margin: auto;" />

```r
# add metadata
# make track dataframe and edit names to match metadata
track_df <- as.data.frame(track) %>% 
  rownames_to_column(var = "DNA_ID") %>%
  mutate(perc_reads_retained = 100 * nochim / input)

#merge
meta_track <- left_join(track_df, full_metadata, by = "DNA_ID")

#DNA Concentrations graph
dna_nochim <- meta_track %>% 
  pivot_longer(input:nochim, names_to = "read_type", values_to = "num_reads") %>%
  mutate(read_type = fct_relevel(read_type, 
                                 "input", "filtered", "denoisedF", "denoisedR", "merged", "nochim")) %>%
  filter(read_type == "nochim")

dna_nochim %>% 
  ggplot(aes(x = Concentration_ng_ul, y = num_reads)) + 
  geom_point() + 
  facet_wrap(~Depth_Class)+
  theme_classic() + 
  labs(x = "DNA Concentration (ng/uL)", y = "Number of reads") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
```

<img src="../figures/01_DADA2_Processing/track-seqs-2.png" style="display: block; margin: auto;" />

```r
dna_nochim %>% 
  ggplot(aes(x = Concentration_ng_ul, y = perc_reads_retained)) + 
  geom_point() + 
  facet_wrap(~Depth_Class)+
  theme_classic() + 
  labs(x = "DNA Concentration (ng/uL)", y = "Percent reads retained") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
```

<img src="../figures/01_DADA2_Processing/track-seqs-3.png" style="display: block; margin: auto;" />

Read retention looks good, and more importantly, consistent between depths and across original DNA concentrations. There is some drop in read counts as DNA concentrations increase. My hypothesis is that those are very high biomass (brown) samples whose DNA extractions may have been a little dirtier, and so in the end produced slighlty worse read quality. Or, higher sample diversity within them could lead to worse read quality/ability to distinguish ASVs.


# Prepare the data for export! 
## 1. ASV Table 


```r
# Prep the asv table! 

samples_out <- rownames(seqtab_nc_len)

# Pull out sample names from the fastq file name 

sample_names_reformatted <- str_remove(samples_out, "_R1_filtered.fastq.gz")

# Replace the names in our seqtable 
rownames(seqtab_nc_len) <- sample_names_reformatted

### intuition check 
stopifnot(rownames(seqtab_nc_len) == sample_names_reformatted)

############## Modify the ASV names and then save a fasta file! 
# Give headers more manageable names
# First pull the ASV sequences
asv_seqs <- colnames(seqtab_nc_len)

# make headers for our ASV seq fasta file, which will be our asv names
asv_headers <- vector(dim(seqtab_nc_len)[2], mode = "character")

# loop through vector and fill it in with ASV names 

for (i in 1:dim(seqtab_nc_len)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}

# intitution check
asv_headers
```

```
##    [1] ">ASV_1"    ">ASV_2"    ">ASV_3"    ">ASV_4"    ">ASV_5"    ">ASV_6"    ">ASV_7"    ">ASV_8"    ">ASV_9"    ">ASV_10"   ">ASV_11"   ">ASV_12"   ">ASV_13"   ">ASV_14"   ">ASV_15"   ">ASV_16"   ">ASV_17"   ">ASV_18"   ">ASV_19"   ">ASV_20"  
##   [21] ">ASV_21"   ">ASV_22"   ">ASV_23"   ">ASV_24"   ">ASV_25"   ">ASV_26"   ">ASV_27"   ">ASV_28"   ">ASV_29"   ">ASV_30"   ">ASV_31"   ">ASV_32"   ">ASV_33"   ">ASV_34"   ">ASV_35"   ">ASV_36"   ">ASV_37"   ">ASV_38"   ">ASV_39"   ">ASV_40"  
##   [41] ">ASV_41"   ">ASV_42"   ">ASV_43"   ">ASV_44"   ">ASV_45"   ">ASV_46"   ">ASV_47"   ">ASV_48"   ">ASV_49"   ">ASV_50"   ">ASV_51"   ">ASV_52"   ">ASV_53"   ">ASV_54"   ">ASV_55"   ">ASV_56"   ">ASV_57"   ">ASV_58"   ">ASV_59"   ">ASV_60"  
##   [61] ">ASV_61"   ">ASV_62"   ">ASV_63"   ">ASV_64"   ">ASV_65"   ">ASV_66"   ">ASV_67"   ">ASV_68"   ">ASV_69"   ">ASV_70"   ">ASV_71"   ">ASV_72"   ">ASV_73"   ">ASV_74"   ">ASV_75"   ">ASV_76"   ">ASV_77"   ">ASV_78"   ">ASV_79"   ">ASV_80"  
##   [81] ">ASV_81"   ">ASV_82"   ">ASV_83"   ">ASV_84"   ">ASV_85"   ">ASV_86"   ">ASV_87"   ">ASV_88"   ">ASV_89"   ">ASV_90"   ">ASV_91"   ">ASV_92"   ">ASV_93"   ">ASV_94"   ">ASV_95"   ">ASV_96"   ">ASV_97"   ">ASV_98"   ">ASV_99"   ">ASV_100" 
##  [101] ">ASV_101"  ">ASV_102"  ">ASV_103"  ">ASV_104"  ">ASV_105"  ">ASV_106"  ">ASV_107"  ">ASV_108"  ">ASV_109"  ">ASV_110"  ">ASV_111"  ">ASV_112"  ">ASV_113"  ">ASV_114"  ">ASV_115"  ">ASV_116"  ">ASV_117"  ">ASV_118"  ">ASV_119"  ">ASV_120" 
##  [121] ">ASV_121"  ">ASV_122"  ">ASV_123"  ">ASV_124"  ">ASV_125"  ">ASV_126"  ">ASV_127"  ">ASV_128"  ">ASV_129"  ">ASV_130"  ">ASV_131"  ">ASV_132"  ">ASV_133"  ">ASV_134"  ">ASV_135"  ">ASV_136"  ">ASV_137"  ">ASV_138"  ">ASV_139"  ">ASV_140" 
##  [141] ">ASV_141"  ">ASV_142"  ">ASV_143"  ">ASV_144"  ">ASV_145"  ">ASV_146"  ">ASV_147"  ">ASV_148"  ">ASV_149"  ">ASV_150"  ">ASV_151"  ">ASV_152"  ">ASV_153"  ">ASV_154"  ">ASV_155"  ">ASV_156"  ">ASV_157"  ">ASV_158"  ">ASV_159"  ">ASV_160" 
##  [161] ">ASV_161"  ">ASV_162"  ">ASV_163"  ">ASV_164"  ">ASV_165"  ">ASV_166"  ">ASV_167"  ">ASV_168"  ">ASV_169"  ">ASV_170"  ">ASV_171"  ">ASV_172"  ">ASV_173"  ">ASV_174"  ">ASV_175"  ">ASV_176"  ">ASV_177"  ">ASV_178"  ">ASV_179"  ">ASV_180" 
##  [181] ">ASV_181"  ">ASV_182"  ">ASV_183"  ">ASV_184"  ">ASV_185"  ">ASV_186"  ">ASV_187"  ">ASV_188"  ">ASV_189"  ">ASV_190"  ">ASV_191"  ">ASV_192"  ">ASV_193"  ">ASV_194"  ">ASV_195"  ">ASV_196"  ">ASV_197"  ">ASV_198"  ">ASV_199"  ">ASV_200" 
##  [201] ">ASV_201"  ">ASV_202"  ">ASV_203"  ">ASV_204"  ">ASV_205"  ">ASV_206"  ">ASV_207"  ">ASV_208"  ">ASV_209"  ">ASV_210"  ">ASV_211"  ">ASV_212"  ">ASV_213"  ">ASV_214"  ">ASV_215"  ">ASV_216"  ">ASV_217"  ">ASV_218"  ">ASV_219"  ">ASV_220" 
##  [221] ">ASV_221"  ">ASV_222"  ">ASV_223"  ">ASV_224"  ">ASV_225"  ">ASV_226"  ">ASV_227"  ">ASV_228"  ">ASV_229"  ">ASV_230"  ">ASV_231"  ">ASV_232"  ">ASV_233"  ">ASV_234"  ">ASV_235"  ">ASV_236"  ">ASV_237"  ">ASV_238"  ">ASV_239"  ">ASV_240" 
##  [241] ">ASV_241"  ">ASV_242"  ">ASV_243"  ">ASV_244"  ">ASV_245"  ">ASV_246"  ">ASV_247"  ">ASV_248"  ">ASV_249"  ">ASV_250"  ">ASV_251"  ">ASV_252"  ">ASV_253"  ">ASV_254"  ">ASV_255"  ">ASV_256"  ">ASV_257"  ">ASV_258"  ">ASV_259"  ">ASV_260" 
##  [261] ">ASV_261"  ">ASV_262"  ">ASV_263"  ">ASV_264"  ">ASV_265"  ">ASV_266"  ">ASV_267"  ">ASV_268"  ">ASV_269"  ">ASV_270"  ">ASV_271"  ">ASV_272"  ">ASV_273"  ">ASV_274"  ">ASV_275"  ">ASV_276"  ">ASV_277"  ">ASV_278"  ">ASV_279"  ">ASV_280" 
##  [281] ">ASV_281"  ">ASV_282"  ">ASV_283"  ">ASV_284"  ">ASV_285"  ">ASV_286"  ">ASV_287"  ">ASV_288"  ">ASV_289"  ">ASV_290"  ">ASV_291"  ">ASV_292"  ">ASV_293"  ">ASV_294"  ">ASV_295"  ">ASV_296"  ">ASV_297"  ">ASV_298"  ">ASV_299"  ">ASV_300" 
##  [301] ">ASV_301"  ">ASV_302"  ">ASV_303"  ">ASV_304"  ">ASV_305"  ">ASV_306"  ">ASV_307"  ">ASV_308"  ">ASV_309"  ">ASV_310"  ">ASV_311"  ">ASV_312"  ">ASV_313"  ">ASV_314"  ">ASV_315"  ">ASV_316"  ">ASV_317"  ">ASV_318"  ">ASV_319"  ">ASV_320" 
##  [321] ">ASV_321"  ">ASV_322"  ">ASV_323"  ">ASV_324"  ">ASV_325"  ">ASV_326"  ">ASV_327"  ">ASV_328"  ">ASV_329"  ">ASV_330"  ">ASV_331"  ">ASV_332"  ">ASV_333"  ">ASV_334"  ">ASV_335"  ">ASV_336"  ">ASV_337"  ">ASV_338"  ">ASV_339"  ">ASV_340" 
##  [341] ">ASV_341"  ">ASV_342"  ">ASV_343"  ">ASV_344"  ">ASV_345"  ">ASV_346"  ">ASV_347"  ">ASV_348"  ">ASV_349"  ">ASV_350"  ">ASV_351"  ">ASV_352"  ">ASV_353"  ">ASV_354"  ">ASV_355"  ">ASV_356"  ">ASV_357"  ">ASV_358"  ">ASV_359"  ">ASV_360" 
##  [361] ">ASV_361"  ">ASV_362"  ">ASV_363"  ">ASV_364"  ">ASV_365"  ">ASV_366"  ">ASV_367"  ">ASV_368"  ">ASV_369"  ">ASV_370"  ">ASV_371"  ">ASV_372"  ">ASV_373"  ">ASV_374"  ">ASV_375"  ">ASV_376"  ">ASV_377"  ">ASV_378"  ">ASV_379"  ">ASV_380" 
##  [381] ">ASV_381"  ">ASV_382"  ">ASV_383"  ">ASV_384"  ">ASV_385"  ">ASV_386"  ">ASV_387"  ">ASV_388"  ">ASV_389"  ">ASV_390"  ">ASV_391"  ">ASV_392"  ">ASV_393"  ">ASV_394"  ">ASV_395"  ">ASV_396"  ">ASV_397"  ">ASV_398"  ">ASV_399"  ">ASV_400" 
##  [401] ">ASV_401"  ">ASV_402"  ">ASV_403"  ">ASV_404"  ">ASV_405"  ">ASV_406"  ">ASV_407"  ">ASV_408"  ">ASV_409"  ">ASV_410"  ">ASV_411"  ">ASV_412"  ">ASV_413"  ">ASV_414"  ">ASV_415"  ">ASV_416"  ">ASV_417"  ">ASV_418"  ">ASV_419"  ">ASV_420" 
##  [421] ">ASV_421"  ">ASV_422"  ">ASV_423"  ">ASV_424"  ">ASV_425"  ">ASV_426"  ">ASV_427"  ">ASV_428"  ">ASV_429"  ">ASV_430"  ">ASV_431"  ">ASV_432"  ">ASV_433"  ">ASV_434"  ">ASV_435"  ">ASV_436"  ">ASV_437"  ">ASV_438"  ">ASV_439"  ">ASV_440" 
##  [441] ">ASV_441"  ">ASV_442"  ">ASV_443"  ">ASV_444"  ">ASV_445"  ">ASV_446"  ">ASV_447"  ">ASV_448"  ">ASV_449"  ">ASV_450"  ">ASV_451"  ">ASV_452"  ">ASV_453"  ">ASV_454"  ">ASV_455"  ">ASV_456"  ">ASV_457"  ">ASV_458"  ">ASV_459"  ">ASV_460" 
##  [461] ">ASV_461"  ">ASV_462"  ">ASV_463"  ">ASV_464"  ">ASV_465"  ">ASV_466"  ">ASV_467"  ">ASV_468"  ">ASV_469"  ">ASV_470"  ">ASV_471"  ">ASV_472"  ">ASV_473"  ">ASV_474"  ">ASV_475"  ">ASV_476"  ">ASV_477"  ">ASV_478"  ">ASV_479"  ">ASV_480" 
##  [481] ">ASV_481"  ">ASV_482"  ">ASV_483"  ">ASV_484"  ">ASV_485"  ">ASV_486"  ">ASV_487"  ">ASV_488"  ">ASV_489"  ">ASV_490"  ">ASV_491"  ">ASV_492"  ">ASV_493"  ">ASV_494"  ">ASV_495"  ">ASV_496"  ">ASV_497"  ">ASV_498"  ">ASV_499"  ">ASV_500" 
##  [501] ">ASV_501"  ">ASV_502"  ">ASV_503"  ">ASV_504"  ">ASV_505"  ">ASV_506"  ">ASV_507"  ">ASV_508"  ">ASV_509"  ">ASV_510"  ">ASV_511"  ">ASV_512"  ">ASV_513"  ">ASV_514"  ">ASV_515"  ">ASV_516"  ">ASV_517"  ">ASV_518"  ">ASV_519"  ">ASV_520" 
##  [521] ">ASV_521"  ">ASV_522"  ">ASV_523"  ">ASV_524"  ">ASV_525"  ">ASV_526"  ">ASV_527"  ">ASV_528"  ">ASV_529"  ">ASV_530"  ">ASV_531"  ">ASV_532"  ">ASV_533"  ">ASV_534"  ">ASV_535"  ">ASV_536"  ">ASV_537"  ">ASV_538"  ">ASV_539"  ">ASV_540" 
##  [541] ">ASV_541"  ">ASV_542"  ">ASV_543"  ">ASV_544"  ">ASV_545"  ">ASV_546"  ">ASV_547"  ">ASV_548"  ">ASV_549"  ">ASV_550"  ">ASV_551"  ">ASV_552"  ">ASV_553"  ">ASV_554"  ">ASV_555"  ">ASV_556"  ">ASV_557"  ">ASV_558"  ">ASV_559"  ">ASV_560" 
##  [561] ">ASV_561"  ">ASV_562"  ">ASV_563"  ">ASV_564"  ">ASV_565"  ">ASV_566"  ">ASV_567"  ">ASV_568"  ">ASV_569"  ">ASV_570"  ">ASV_571"  ">ASV_572"  ">ASV_573"  ">ASV_574"  ">ASV_575"  ">ASV_576"  ">ASV_577"  ">ASV_578"  ">ASV_579"  ">ASV_580" 
##  [581] ">ASV_581"  ">ASV_582"  ">ASV_583"  ">ASV_584"  ">ASV_585"  ">ASV_586"  ">ASV_587"  ">ASV_588"  ">ASV_589"  ">ASV_590"  ">ASV_591"  ">ASV_592"  ">ASV_593"  ">ASV_594"  ">ASV_595"  ">ASV_596"  ">ASV_597"  ">ASV_598"  ">ASV_599"  ">ASV_600" 
##  [601] ">ASV_601"  ">ASV_602"  ">ASV_603"  ">ASV_604"  ">ASV_605"  ">ASV_606"  ">ASV_607"  ">ASV_608"  ">ASV_609"  ">ASV_610"  ">ASV_611"  ">ASV_612"  ">ASV_613"  ">ASV_614"  ">ASV_615"  ">ASV_616"  ">ASV_617"  ">ASV_618"  ">ASV_619"  ">ASV_620" 
##  [621] ">ASV_621"  ">ASV_622"  ">ASV_623"  ">ASV_624"  ">ASV_625"  ">ASV_626"  ">ASV_627"  ">ASV_628"  ">ASV_629"  ">ASV_630"  ">ASV_631"  ">ASV_632"  ">ASV_633"  ">ASV_634"  ">ASV_635"  ">ASV_636"  ">ASV_637"  ">ASV_638"  ">ASV_639"  ">ASV_640" 
##  [641] ">ASV_641"  ">ASV_642"  ">ASV_643"  ">ASV_644"  ">ASV_645"  ">ASV_646"  ">ASV_647"  ">ASV_648"  ">ASV_649"  ">ASV_650"  ">ASV_651"  ">ASV_652"  ">ASV_653"  ">ASV_654"  ">ASV_655"  ">ASV_656"  ">ASV_657"  ">ASV_658"  ">ASV_659"  ">ASV_660" 
##  [661] ">ASV_661"  ">ASV_662"  ">ASV_663"  ">ASV_664"  ">ASV_665"  ">ASV_666"  ">ASV_667"  ">ASV_668"  ">ASV_669"  ">ASV_670"  ">ASV_671"  ">ASV_672"  ">ASV_673"  ">ASV_674"  ">ASV_675"  ">ASV_676"  ">ASV_677"  ">ASV_678"  ">ASV_679"  ">ASV_680" 
##  [681] ">ASV_681"  ">ASV_682"  ">ASV_683"  ">ASV_684"  ">ASV_685"  ">ASV_686"  ">ASV_687"  ">ASV_688"  ">ASV_689"  ">ASV_690"  ">ASV_691"  ">ASV_692"  ">ASV_693"  ">ASV_694"  ">ASV_695"  ">ASV_696"  ">ASV_697"  ">ASV_698"  ">ASV_699"  ">ASV_700" 
##  [701] ">ASV_701"  ">ASV_702"  ">ASV_703"  ">ASV_704"  ">ASV_705"  ">ASV_706"  ">ASV_707"  ">ASV_708"  ">ASV_709"  ">ASV_710"  ">ASV_711"  ">ASV_712"  ">ASV_713"  ">ASV_714"  ">ASV_715"  ">ASV_716"  ">ASV_717"  ">ASV_718"  ">ASV_719"  ">ASV_720" 
##  [721] ">ASV_721"  ">ASV_722"  ">ASV_723"  ">ASV_724"  ">ASV_725"  ">ASV_726"  ">ASV_727"  ">ASV_728"  ">ASV_729"  ">ASV_730"  ">ASV_731"  ">ASV_732"  ">ASV_733"  ">ASV_734"  ">ASV_735"  ">ASV_736"  ">ASV_737"  ">ASV_738"  ">ASV_739"  ">ASV_740" 
##  [741] ">ASV_741"  ">ASV_742"  ">ASV_743"  ">ASV_744"  ">ASV_745"  ">ASV_746"  ">ASV_747"  ">ASV_748"  ">ASV_749"  ">ASV_750"  ">ASV_751"  ">ASV_752"  ">ASV_753"  ">ASV_754"  ">ASV_755"  ">ASV_756"  ">ASV_757"  ">ASV_758"  ">ASV_759"  ">ASV_760" 
##  [761] ">ASV_761"  ">ASV_762"  ">ASV_763"  ">ASV_764"  ">ASV_765"  ">ASV_766"  ">ASV_767"  ">ASV_768"  ">ASV_769"  ">ASV_770"  ">ASV_771"  ">ASV_772"  ">ASV_773"  ">ASV_774"  ">ASV_775"  ">ASV_776"  ">ASV_777"  ">ASV_778"  ">ASV_779"  ">ASV_780" 
##  [781] ">ASV_781"  ">ASV_782"  ">ASV_783"  ">ASV_784"  ">ASV_785"  ">ASV_786"  ">ASV_787"  ">ASV_788"  ">ASV_789"  ">ASV_790"  ">ASV_791"  ">ASV_792"  ">ASV_793"  ">ASV_794"  ">ASV_795"  ">ASV_796"  ">ASV_797"  ">ASV_798"  ">ASV_799"  ">ASV_800" 
##  [801] ">ASV_801"  ">ASV_802"  ">ASV_803"  ">ASV_804"  ">ASV_805"  ">ASV_806"  ">ASV_807"  ">ASV_808"  ">ASV_809"  ">ASV_810"  ">ASV_811"  ">ASV_812"  ">ASV_813"  ">ASV_814"  ">ASV_815"  ">ASV_816"  ">ASV_817"  ">ASV_818"  ">ASV_819"  ">ASV_820" 
##  [821] ">ASV_821"  ">ASV_822"  ">ASV_823"  ">ASV_824"  ">ASV_825"  ">ASV_826"  ">ASV_827"  ">ASV_828"  ">ASV_829"  ">ASV_830"  ">ASV_831"  ">ASV_832"  ">ASV_833"  ">ASV_834"  ">ASV_835"  ">ASV_836"  ">ASV_837"  ">ASV_838"  ">ASV_839"  ">ASV_840" 
##  [841] ">ASV_841"  ">ASV_842"  ">ASV_843"  ">ASV_844"  ">ASV_845"  ">ASV_846"  ">ASV_847"  ">ASV_848"  ">ASV_849"  ">ASV_850"  ">ASV_851"  ">ASV_852"  ">ASV_853"  ">ASV_854"  ">ASV_855"  ">ASV_856"  ">ASV_857"  ">ASV_858"  ">ASV_859"  ">ASV_860" 
##  [861] ">ASV_861"  ">ASV_862"  ">ASV_863"  ">ASV_864"  ">ASV_865"  ">ASV_866"  ">ASV_867"  ">ASV_868"  ">ASV_869"  ">ASV_870"  ">ASV_871"  ">ASV_872"  ">ASV_873"  ">ASV_874"  ">ASV_875"  ">ASV_876"  ">ASV_877"  ">ASV_878"  ">ASV_879"  ">ASV_880" 
##  [881] ">ASV_881"  ">ASV_882"  ">ASV_883"  ">ASV_884"  ">ASV_885"  ">ASV_886"  ">ASV_887"  ">ASV_888"  ">ASV_889"  ">ASV_890"  ">ASV_891"  ">ASV_892"  ">ASV_893"  ">ASV_894"  ">ASV_895"  ">ASV_896"  ">ASV_897"  ">ASV_898"  ">ASV_899"  ">ASV_900" 
##  [901] ">ASV_901"  ">ASV_902"  ">ASV_903"  ">ASV_904"  ">ASV_905"  ">ASV_906"  ">ASV_907"  ">ASV_908"  ">ASV_909"  ">ASV_910"  ">ASV_911"  ">ASV_912"  ">ASV_913"  ">ASV_914"  ">ASV_915"  ">ASV_916"  ">ASV_917"  ">ASV_918"  ">ASV_919"  ">ASV_920" 
##  [921] ">ASV_921"  ">ASV_922"  ">ASV_923"  ">ASV_924"  ">ASV_925"  ">ASV_926"  ">ASV_927"  ">ASV_928"  ">ASV_929"  ">ASV_930"  ">ASV_931"  ">ASV_932"  ">ASV_933"  ">ASV_934"  ">ASV_935"  ">ASV_936"  ">ASV_937"  ">ASV_938"  ">ASV_939"  ">ASV_940" 
##  [941] ">ASV_941"  ">ASV_942"  ">ASV_943"  ">ASV_944"  ">ASV_945"  ">ASV_946"  ">ASV_947"  ">ASV_948"  ">ASV_949"  ">ASV_950"  ">ASV_951"  ">ASV_952"  ">ASV_953"  ">ASV_954"  ">ASV_955"  ">ASV_956"  ">ASV_957"  ">ASV_958"  ">ASV_959"  ">ASV_960" 
##  [961] ">ASV_961"  ">ASV_962"  ">ASV_963"  ">ASV_964"  ">ASV_965"  ">ASV_966"  ">ASV_967"  ">ASV_968"  ">ASV_969"  ">ASV_970"  ">ASV_971"  ">ASV_972"  ">ASV_973"  ">ASV_974"  ">ASV_975"  ">ASV_976"  ">ASV_977"  ">ASV_978"  ">ASV_979"  ">ASV_980" 
##  [981] ">ASV_981"  ">ASV_982"  ">ASV_983"  ">ASV_984"  ">ASV_985"  ">ASV_986"  ">ASV_987"  ">ASV_988"  ">ASV_989"  ">ASV_990"  ">ASV_991"  ">ASV_992"  ">ASV_993"  ">ASV_994"  ">ASV_995"  ">ASV_996"  ">ASV_997"  ">ASV_998"  ">ASV_999"  ">ASV_1000"
## [1001] ">ASV_1001" ">ASV_1002" ">ASV_1003" ">ASV_1004" ">ASV_1005" ">ASV_1006" ">ASV_1007" ">ASV_1008" ">ASV_1009" ">ASV_1010" ">ASV_1011" ">ASV_1012" ">ASV_1013" ">ASV_1014" ">ASV_1015" ">ASV_1016" ">ASV_1017" ">ASV_1018" ">ASV_1019" ">ASV_1020"
## [1021] ">ASV_1021" ">ASV_1022" ">ASV_1023" ">ASV_1024" ">ASV_1025" ">ASV_1026" ">ASV_1027" ">ASV_1028" ">ASV_1029" ">ASV_1030" ">ASV_1031" ">ASV_1032" ">ASV_1033" ">ASV_1034" ">ASV_1035" ">ASV_1036" ">ASV_1037" ">ASV_1038" ">ASV_1039" ">ASV_1040"
## [1041] ">ASV_1041" ">ASV_1042" ">ASV_1043" ">ASV_1044" ">ASV_1045" ">ASV_1046" ">ASV_1047" ">ASV_1048" ">ASV_1049" ">ASV_1050" ">ASV_1051" ">ASV_1052" ">ASV_1053" ">ASV_1054" ">ASV_1055" ">ASV_1056" ">ASV_1057" ">ASV_1058" ">ASV_1059" ">ASV_1060"
## [1061] ">ASV_1061" ">ASV_1062" ">ASV_1063" ">ASV_1064" ">ASV_1065" ">ASV_1066" ">ASV_1067" ">ASV_1068" ">ASV_1069" ">ASV_1070" ">ASV_1071" ">ASV_1072" ">ASV_1073" ">ASV_1074" ">ASV_1075" ">ASV_1076" ">ASV_1077" ">ASV_1078" ">ASV_1079" ">ASV_1080"
## [1081] ">ASV_1081" ">ASV_1082" ">ASV_1083" ">ASV_1084" ">ASV_1085" ">ASV_1086" ">ASV_1087" ">ASV_1088" ">ASV_1089" ">ASV_1090" ">ASV_1091" ">ASV_1092" ">ASV_1093" ">ASV_1094" ">ASV_1095" ">ASV_1096" ">ASV_1097" ">ASV_1098" ">ASV_1099" ">ASV_1100"
## [1101] ">ASV_1101" ">ASV_1102" ">ASV_1103" ">ASV_1104" ">ASV_1105" ">ASV_1106" ">ASV_1107" ">ASV_1108" ">ASV_1109" ">ASV_1110" ">ASV_1111" ">ASV_1112" ">ASV_1113" ">ASV_1114" ">ASV_1115" ">ASV_1116" ">ASV_1117" ">ASV_1118" ">ASV_1119" ">ASV_1120"
## [1121] ">ASV_1121" ">ASV_1122" ">ASV_1123" ">ASV_1124" ">ASV_1125" ">ASV_1126" ">ASV_1127" ">ASV_1128" ">ASV_1129" ">ASV_1130" ">ASV_1131" ">ASV_1132" ">ASV_1133" ">ASV_1134" ">ASV_1135" ">ASV_1136" ">ASV_1137" ">ASV_1138" ">ASV_1139" ">ASV_1140"
## [1141] ">ASV_1141" ">ASV_1142" ">ASV_1143" ">ASV_1144" ">ASV_1145" ">ASV_1146" ">ASV_1147" ">ASV_1148" ">ASV_1149" ">ASV_1150" ">ASV_1151" ">ASV_1152" ">ASV_1153" ">ASV_1154" ">ASV_1155" ">ASV_1156" ">ASV_1157" ">ASV_1158" ">ASV_1159" ">ASV_1160"
## [1161] ">ASV_1161" ">ASV_1162" ">ASV_1163" ">ASV_1164" ">ASV_1165" ">ASV_1166" ">ASV_1167" ">ASV_1168" ">ASV_1169" ">ASV_1170" ">ASV_1171" ">ASV_1172" ">ASV_1173" ">ASV_1174" ">ASV_1175" ">ASV_1176" ">ASV_1177" ">ASV_1178" ">ASV_1179" ">ASV_1180"
## [1181] ">ASV_1181" ">ASV_1182" ">ASV_1183" ">ASV_1184" ">ASV_1185" ">ASV_1186" ">ASV_1187" ">ASV_1188" ">ASV_1189" ">ASV_1190" ">ASV_1191" ">ASV_1192" ">ASV_1193" ">ASV_1194" ">ASV_1195" ">ASV_1196" ">ASV_1197" ">ASV_1198" ">ASV_1199" ">ASV_1200"
## [1201] ">ASV_1201" ">ASV_1202" ">ASV_1203" ">ASV_1204" ">ASV_1205" ">ASV_1206" ">ASV_1207" ">ASV_1208" ">ASV_1209" ">ASV_1210" ">ASV_1211" ">ASV_1212" ">ASV_1213" ">ASV_1214" ">ASV_1215" ">ASV_1216" ">ASV_1217" ">ASV_1218" ">ASV_1219" ">ASV_1220"
## [1221] ">ASV_1221" ">ASV_1222" ">ASV_1223" ">ASV_1224" ">ASV_1225" ">ASV_1226" ">ASV_1227" ">ASV_1228" ">ASV_1229" ">ASV_1230" ">ASV_1231" ">ASV_1232" ">ASV_1233" ">ASV_1234" ">ASV_1235" ">ASV_1236" ">ASV_1237" ">ASV_1238" ">ASV_1239" ">ASV_1240"
## [1241] ">ASV_1241" ">ASV_1242" ">ASV_1243" ">ASV_1244" ">ASV_1245" ">ASV_1246" ">ASV_1247" ">ASV_1248" ">ASV_1249" ">ASV_1250" ">ASV_1251" ">ASV_1252" ">ASV_1253" ">ASV_1254" ">ASV_1255" ">ASV_1256" ">ASV_1257" ">ASV_1258" ">ASV_1259" ">ASV_1260"
## [1261] ">ASV_1261" ">ASV_1262" ">ASV_1263" ">ASV_1264" ">ASV_1265" ">ASV_1266" ">ASV_1267" ">ASV_1268" ">ASV_1269" ">ASV_1270" ">ASV_1271" ">ASV_1272" ">ASV_1273" ">ASV_1274" ">ASV_1275" ">ASV_1276" ">ASV_1277" ">ASV_1278" ">ASV_1279" ">ASV_1280"
## [1281] ">ASV_1281" ">ASV_1282" ">ASV_1283" ">ASV_1284" ">ASV_1285" ">ASV_1286" ">ASV_1287" ">ASV_1288" ">ASV_1289" ">ASV_1290" ">ASV_1291" ">ASV_1292" ">ASV_1293" ">ASV_1294" ">ASV_1295" ">ASV_1296" ">ASV_1297" ">ASV_1298" ">ASV_1299" ">ASV_1300"
## [1301] ">ASV_1301" ">ASV_1302" ">ASV_1303" ">ASV_1304" ">ASV_1305" ">ASV_1306" ">ASV_1307" ">ASV_1308" ">ASV_1309" ">ASV_1310" ">ASV_1311" ">ASV_1312" ">ASV_1313" ">ASV_1314" ">ASV_1315" ">ASV_1316" ">ASV_1317" ">ASV_1318" ">ASV_1319" ">ASV_1320"
## [1321] ">ASV_1321" ">ASV_1322" ">ASV_1323" ">ASV_1324" ">ASV_1325" ">ASV_1326" ">ASV_1327" ">ASV_1328" ">ASV_1329" ">ASV_1330" ">ASV_1331" ">ASV_1332" ">ASV_1333" ">ASV_1334" ">ASV_1335" ">ASV_1336" ">ASV_1337" ">ASV_1338" ">ASV_1339" ">ASV_1340"
## [1341] ">ASV_1341" ">ASV_1342" ">ASV_1343" ">ASV_1344" ">ASV_1345" ">ASV_1346" ">ASV_1347" ">ASV_1348" ">ASV_1349" ">ASV_1350" ">ASV_1351" ">ASV_1352" ">ASV_1353" ">ASV_1354" ">ASV_1355" ">ASV_1356" ">ASV_1357" ">ASV_1358" ">ASV_1359" ">ASV_1360"
## [1361] ">ASV_1361" ">ASV_1362" ">ASV_1363" ">ASV_1364" ">ASV_1365" ">ASV_1366" ">ASV_1367" ">ASV_1368" ">ASV_1369" ">ASV_1370" ">ASV_1371" ">ASV_1372" ">ASV_1373" ">ASV_1374" ">ASV_1375" ">ASV_1376" ">ASV_1377" ">ASV_1378" ">ASV_1379" ">ASV_1380"
## [1381] ">ASV_1381" ">ASV_1382" ">ASV_1383" ">ASV_1384" ">ASV_1385" ">ASV_1386" ">ASV_1387" ">ASV_1388" ">ASV_1389" ">ASV_1390" ">ASV_1391" ">ASV_1392" ">ASV_1393" ">ASV_1394" ">ASV_1395" ">ASV_1396" ">ASV_1397" ">ASV_1398" ">ASV_1399" ">ASV_1400"
## [1401] ">ASV_1401" ">ASV_1402" ">ASV_1403" ">ASV_1404" ">ASV_1405" ">ASV_1406" ">ASV_1407" ">ASV_1408" ">ASV_1409" ">ASV_1410" ">ASV_1411" ">ASV_1412" ">ASV_1413" ">ASV_1414" ">ASV_1415" ">ASV_1416" ">ASV_1417" ">ASV_1418" ">ASV_1419" ">ASV_1420"
## [1421] ">ASV_1421" ">ASV_1422" ">ASV_1423" ">ASV_1424" ">ASV_1425" ">ASV_1426" ">ASV_1427" ">ASV_1428" ">ASV_1429" ">ASV_1430" ">ASV_1431" ">ASV_1432" ">ASV_1433" ">ASV_1434" ">ASV_1435" ">ASV_1436" ">ASV_1437" ">ASV_1438" ">ASV_1439" ">ASV_1440"
## [1441] ">ASV_1441" ">ASV_1442" ">ASV_1443" ">ASV_1444" ">ASV_1445" ">ASV_1446" ">ASV_1447" ">ASV_1448" ">ASV_1449" ">ASV_1450" ">ASV_1451" ">ASV_1452" ">ASV_1453" ">ASV_1454" ">ASV_1455" ">ASV_1456" ">ASV_1457" ">ASV_1458" ">ASV_1459" ">ASV_1460"
## [1461] ">ASV_1461" ">ASV_1462" ">ASV_1463" ">ASV_1464" ">ASV_1465" ">ASV_1466" ">ASV_1467" ">ASV_1468" ">ASV_1469" ">ASV_1470" ">ASV_1471" ">ASV_1472" ">ASV_1473" ">ASV_1474" ">ASV_1475" ">ASV_1476" ">ASV_1477" ">ASV_1478" ">ASV_1479" ">ASV_1480"
## [1481] ">ASV_1481" ">ASV_1482" ">ASV_1483" ">ASV_1484" ">ASV_1485" ">ASV_1486" ">ASV_1487" ">ASV_1488" ">ASV_1489" ">ASV_1490" ">ASV_1491" ">ASV_1492" ">ASV_1493" ">ASV_1494" ">ASV_1495" ">ASV_1496" ">ASV_1497" ">ASV_1498" ">ASV_1499" ">ASV_1500"
## [1501] ">ASV_1501" ">ASV_1502" ">ASV_1503" ">ASV_1504" ">ASV_1505" ">ASV_1506" ">ASV_1507" ">ASV_1508" ">ASV_1509" ">ASV_1510" ">ASV_1511" ">ASV_1512" ">ASV_1513" ">ASV_1514" ">ASV_1515" ">ASV_1516" ">ASV_1517" ">ASV_1518" ">ASV_1519" ">ASV_1520"
## [1521] ">ASV_1521" ">ASV_1522" ">ASV_1523" ">ASV_1524" ">ASV_1525" ">ASV_1526" ">ASV_1527" ">ASV_1528" ">ASV_1529" ">ASV_1530" ">ASV_1531" ">ASV_1532" ">ASV_1533" ">ASV_1534" ">ASV_1535" ">ASV_1536" ">ASV_1537" ">ASV_1538" ">ASV_1539" ">ASV_1540"
## [1541] ">ASV_1541" ">ASV_1542" ">ASV_1543" ">ASV_1544" ">ASV_1545" ">ASV_1546" ">ASV_1547" ">ASV_1548" ">ASV_1549" ">ASV_1550" ">ASV_1551" ">ASV_1552" ">ASV_1553" ">ASV_1554" ">ASV_1555" ">ASV_1556" ">ASV_1557" ">ASV_1558" ">ASV_1559" ">ASV_1560"
## [1561] ">ASV_1561" ">ASV_1562" ">ASV_1563" ">ASV_1564" ">ASV_1565" ">ASV_1566" ">ASV_1567" ">ASV_1568" ">ASV_1569" ">ASV_1570" ">ASV_1571" ">ASV_1572" ">ASV_1573" ">ASV_1574" ">ASV_1575" ">ASV_1576" ">ASV_1577" ">ASV_1578" ">ASV_1579" ">ASV_1580"
## [1581] ">ASV_1581" ">ASV_1582" ">ASV_1583" ">ASV_1584" ">ASV_1585" ">ASV_1586" ">ASV_1587" ">ASV_1588" ">ASV_1589" ">ASV_1590" ">ASV_1591" ">ASV_1592" ">ASV_1593" ">ASV_1594" ">ASV_1595" ">ASV_1596" ">ASV_1597" ">ASV_1598" ">ASV_1599" ">ASV_1600"
## [1601] ">ASV_1601" ">ASV_1602" ">ASV_1603" ">ASV_1604" ">ASV_1605" ">ASV_1606" ">ASV_1607" ">ASV_1608" ">ASV_1609" ">ASV_1610" ">ASV_1611" ">ASV_1612" ">ASV_1613" ">ASV_1614" ">ASV_1615" ">ASV_1616" ">ASV_1617" ">ASV_1618" ">ASV_1619" ">ASV_1620"
## [1621] ">ASV_1621" ">ASV_1622" ">ASV_1623" ">ASV_1624" ">ASV_1625" ">ASV_1626" ">ASV_1627" ">ASV_1628" ">ASV_1629" ">ASV_1630" ">ASV_1631" ">ASV_1632" ">ASV_1633" ">ASV_1634" ">ASV_1635" ">ASV_1636" ">ASV_1637" ">ASV_1638" ">ASV_1639" ">ASV_1640"
## [1641] ">ASV_1641" ">ASV_1642" ">ASV_1643" ">ASV_1644" ">ASV_1645" ">ASV_1646" ">ASV_1647" ">ASV_1648" ">ASV_1649" ">ASV_1650" ">ASV_1651" ">ASV_1652" ">ASV_1653" ">ASV_1654" ">ASV_1655" ">ASV_1656" ">ASV_1657" ">ASV_1658" ">ASV_1659" ">ASV_1660"
## [1661] ">ASV_1661" ">ASV_1662" ">ASV_1663" ">ASV_1664" ">ASV_1665" ">ASV_1666" ">ASV_1667" ">ASV_1668" ">ASV_1669" ">ASV_1670" ">ASV_1671" ">ASV_1672" ">ASV_1673" ">ASV_1674" ">ASV_1675" ">ASV_1676" ">ASV_1677" ">ASV_1678" ">ASV_1679" ">ASV_1680"
## [1681] ">ASV_1681" ">ASV_1682" ">ASV_1683" ">ASV_1684" ">ASV_1685" ">ASV_1686" ">ASV_1687" ">ASV_1688" ">ASV_1689" ">ASV_1690" ">ASV_1691" ">ASV_1692" ">ASV_1693" ">ASV_1694" ">ASV_1695" ">ASV_1696" ">ASV_1697" ">ASV_1698" ">ASV_1699" ">ASV_1700"
## [1701] ">ASV_1701" ">ASV_1702" ">ASV_1703" ">ASV_1704" ">ASV_1705" ">ASV_1706" ">ASV_1707" ">ASV_1708" ">ASV_1709" ">ASV_1710" ">ASV_1711" ">ASV_1712" ">ASV_1713" ">ASV_1714" ">ASV_1715" ">ASV_1716" ">ASV_1717" ">ASV_1718" ">ASV_1719" ">ASV_1720"
## [1721] ">ASV_1721" ">ASV_1722" ">ASV_1723" ">ASV_1724" ">ASV_1725" ">ASV_1726" ">ASV_1727" ">ASV_1728" ">ASV_1729" ">ASV_1730" ">ASV_1731" ">ASV_1732" ">ASV_1733" ">ASV_1734" ">ASV_1735" ">ASV_1736" ">ASV_1737" ">ASV_1738" ">ASV_1739" ">ASV_1740"
## [1741] ">ASV_1741" ">ASV_1742" ">ASV_1743" ">ASV_1744" ">ASV_1745" ">ASV_1746" ">ASV_1747" ">ASV_1748" ">ASV_1749" ">ASV_1750" ">ASV_1751" ">ASV_1752" ">ASV_1753" ">ASV_1754" ">ASV_1755" ">ASV_1756" ">ASV_1757" ">ASV_1758" ">ASV_1759" ">ASV_1760"
## [1761] ">ASV_1761" ">ASV_1762" ">ASV_1763" ">ASV_1764" ">ASV_1765" ">ASV_1766" ">ASV_1767" ">ASV_1768" ">ASV_1769" ">ASV_1770" ">ASV_1771" ">ASV_1772" ">ASV_1773" ">ASV_1774" ">ASV_1775" ">ASV_1776" ">ASV_1777" ">ASV_1778" ">ASV_1779" ">ASV_1780"
## [1781] ">ASV_1781" ">ASV_1782" ">ASV_1783" ">ASV_1784" ">ASV_1785" ">ASV_1786" ">ASV_1787" ">ASV_1788" ">ASV_1789" ">ASV_1790" ">ASV_1791" ">ASV_1792" ">ASV_1793" ">ASV_1794" ">ASV_1795" ">ASV_1796" ">ASV_1797" ">ASV_1798" ">ASV_1799" ">ASV_1800"
## [1801] ">ASV_1801" ">ASV_1802" ">ASV_1803" ">ASV_1804" ">ASV_1805" ">ASV_1806" ">ASV_1807" ">ASV_1808" ">ASV_1809" ">ASV_1810" ">ASV_1811" ">ASV_1812" ">ASV_1813" ">ASV_1814" ">ASV_1815" ">ASV_1816" ">ASV_1817" ">ASV_1818" ">ASV_1819" ">ASV_1820"
## [1821] ">ASV_1821" ">ASV_1822" ">ASV_1823" ">ASV_1824" ">ASV_1825" ">ASV_1826" ">ASV_1827" ">ASV_1828" ">ASV_1829" ">ASV_1830" ">ASV_1831" ">ASV_1832" ">ASV_1833" ">ASV_1834" ">ASV_1835" ">ASV_1836" ">ASV_1837" ">ASV_1838" ">ASV_1839" ">ASV_1840"
## [1841] ">ASV_1841" ">ASV_1842" ">ASV_1843" ">ASV_1844" ">ASV_1845" ">ASV_1846" ">ASV_1847" ">ASV_1848" ">ASV_1849" ">ASV_1850" ">ASV_1851" ">ASV_1852" ">ASV_1853" ">ASV_1854" ">ASV_1855" ">ASV_1856" ">ASV_1857" ">ASV_1858" ">ASV_1859" ">ASV_1860"
## [1861] ">ASV_1861" ">ASV_1862" ">ASV_1863" ">ASV_1864" ">ASV_1865" ">ASV_1866" ">ASV_1867" ">ASV_1868" ">ASV_1869" ">ASV_1870" ">ASV_1871" ">ASV_1872" ">ASV_1873" ">ASV_1874" ">ASV_1875" ">ASV_1876" ">ASV_1877" ">ASV_1878" ">ASV_1879" ">ASV_1880"
## [1881] ">ASV_1881" ">ASV_1882" ">ASV_1883" ">ASV_1884" ">ASV_1885" ">ASV_1886" ">ASV_1887" ">ASV_1888" ">ASV_1889" ">ASV_1890" ">ASV_1891" ">ASV_1892" ">ASV_1893" ">ASV_1894" ">ASV_1895" ">ASV_1896" ">ASV_1897" ">ASV_1898" ">ASV_1899" ">ASV_1900"
## [1901] ">ASV_1901" ">ASV_1902" ">ASV_1903" ">ASV_1904" ">ASV_1905" ">ASV_1906" ">ASV_1907" ">ASV_1908" ">ASV_1909" ">ASV_1910" ">ASV_1911" ">ASV_1912" ">ASV_1913" ">ASV_1914" ">ASV_1915" ">ASV_1916" ">ASV_1917" ">ASV_1918" ">ASV_1919" ">ASV_1920"
## [1921] ">ASV_1921" ">ASV_1922" ">ASV_1923" ">ASV_1924" ">ASV_1925" ">ASV_1926" ">ASV_1927" ">ASV_1928" ">ASV_1929" ">ASV_1930" ">ASV_1931" ">ASV_1932" ">ASV_1933" ">ASV_1934" ">ASV_1935" ">ASV_1936" ">ASV_1937" ">ASV_1938" ">ASV_1939" ">ASV_1940"
## [1941] ">ASV_1941" ">ASV_1942" ">ASV_1943" ">ASV_1944" ">ASV_1945" ">ASV_1946" ">ASV_1947" ">ASV_1948" ">ASV_1949" ">ASV_1950" ">ASV_1951" ">ASV_1952" ">ASV_1953" ">ASV_1954" ">ASV_1955" ">ASV_1956" ">ASV_1957" ">ASV_1958" ">ASV_1959" ">ASV_1960"
## [1961] ">ASV_1961" ">ASV_1962" ">ASV_1963" ">ASV_1964" ">ASV_1965" ">ASV_1966" ">ASV_1967" ">ASV_1968" ">ASV_1969" ">ASV_1970" ">ASV_1971" ">ASV_1972" ">ASV_1973" ">ASV_1974" ">ASV_1975" ">ASV_1976" ">ASV_1977" ">ASV_1978" ">ASV_1979" ">ASV_1980"
## [1981] ">ASV_1981" ">ASV_1982" ">ASV_1983" ">ASV_1984" ">ASV_1985" ">ASV_1986" ">ASV_1987" ">ASV_1988" ">ASV_1989" ">ASV_1990" ">ASV_1991" ">ASV_1992" ">ASV_1993" ">ASV_1994" ">ASV_1995" ">ASV_1996" ">ASV_1997" ">ASV_1998" ">ASV_1999" ">ASV_2000"
## [2001] ">ASV_2001" ">ASV_2002" ">ASV_2003" ">ASV_2004" ">ASV_2005" ">ASV_2006" ">ASV_2007" ">ASV_2008" ">ASV_2009" ">ASV_2010" ">ASV_2011" ">ASV_2012" ">ASV_2013" ">ASV_2014" ">ASV_2015" ">ASV_2016" ">ASV_2017" ">ASV_2018" ">ASV_2019" ">ASV_2020"
## [2021] ">ASV_2021" ">ASV_2022" ">ASV_2023" ">ASV_2024" ">ASV_2025" ">ASV_2026" ">ASV_2027" ">ASV_2028" ">ASV_2029" ">ASV_2030" ">ASV_2031" ">ASV_2032" ">ASV_2033" ">ASV_2034" ">ASV_2035" ">ASV_2036" ">ASV_2037" ">ASV_2038" ">ASV_2039" ">ASV_2040"
## [2041] ">ASV_2041" ">ASV_2042" ">ASV_2043" ">ASV_2044" ">ASV_2045" ">ASV_2046" ">ASV_2047" ">ASV_2048" ">ASV_2049" ">ASV_2050" ">ASV_2051" ">ASV_2052" ">ASV_2053" ">ASV_2054" ">ASV_2055" ">ASV_2056" ">ASV_2057" ">ASV_2058" ">ASV_2059" ">ASV_2060"
## [2061] ">ASV_2061" ">ASV_2062" ">ASV_2063" ">ASV_2064" ">ASV_2065" ">ASV_2066" ">ASV_2067" ">ASV_2068" ">ASV_2069" ">ASV_2070" ">ASV_2071" ">ASV_2072" ">ASV_2073" ">ASV_2074" ">ASV_2075" ">ASV_2076" ">ASV_2077" ">ASV_2078" ">ASV_2079" ">ASV_2080"
## [2081] ">ASV_2081" ">ASV_2082" ">ASV_2083" ">ASV_2084" ">ASV_2085" ">ASV_2086" ">ASV_2087" ">ASV_2088" ">ASV_2089" ">ASV_2090" ">ASV_2091" ">ASV_2092" ">ASV_2093" ">ASV_2094" ">ASV_2095" ">ASV_2096" ">ASV_2097" ">ASV_2098" ">ASV_2099" ">ASV_2100"
## [2101] ">ASV_2101" ">ASV_2102" ">ASV_2103" ">ASV_2104" ">ASV_2105" ">ASV_2106" ">ASV_2107" ">ASV_2108" ">ASV_2109" ">ASV_2110" ">ASV_2111" ">ASV_2112" ">ASV_2113" ">ASV_2114" ">ASV_2115" ">ASV_2116" ">ASV_2117" ">ASV_2118" ">ASV_2119" ">ASV_2120"
## [2121] ">ASV_2121" ">ASV_2122" ">ASV_2123" ">ASV_2124" ">ASV_2125" ">ASV_2126" ">ASV_2127" ">ASV_2128" ">ASV_2129" ">ASV_2130" ">ASV_2131" ">ASV_2132" ">ASV_2133" ">ASV_2134" ">ASV_2135" ">ASV_2136" ">ASV_2137" ">ASV_2138" ">ASV_2139" ">ASV_2140"
## [2141] ">ASV_2141" ">ASV_2142" ">ASV_2143" ">ASV_2144" ">ASV_2145" ">ASV_2146" ">ASV_2147" ">ASV_2148" ">ASV_2149" ">ASV_2150" ">ASV_2151" ">ASV_2152" ">ASV_2153" ">ASV_2154" ">ASV_2155" ">ASV_2156" ">ASV_2157" ">ASV_2158" ">ASV_2159" ">ASV_2160"
## [2161] ">ASV_2161" ">ASV_2162" ">ASV_2163" ">ASV_2164" ">ASV_2165" ">ASV_2166" ">ASV_2167" ">ASV_2168" ">ASV_2169" ">ASV_2170" ">ASV_2171" ">ASV_2172" ">ASV_2173" ">ASV_2174" ">ASV_2175" ">ASV_2176" ">ASV_2177" ">ASV_2178" ">ASV_2179" ">ASV_2180"
## [2181] ">ASV_2181" ">ASV_2182" ">ASV_2183" ">ASV_2184" ">ASV_2185" ">ASV_2186" ">ASV_2187" ">ASV_2188" ">ASV_2189" ">ASV_2190" ">ASV_2191" ">ASV_2192" ">ASV_2193" ">ASV_2194" ">ASV_2195" ">ASV_2196" ">ASV_2197" ">ASV_2198" ">ASV_2199" ">ASV_2200"
## [2201] ">ASV_2201" ">ASV_2202" ">ASV_2203" ">ASV_2204" ">ASV_2205" ">ASV_2206" ">ASV_2207" ">ASV_2208" ">ASV_2209" ">ASV_2210" ">ASV_2211" ">ASV_2212" ">ASV_2213" ">ASV_2214" ">ASV_2215" ">ASV_2216" ">ASV_2217" ">ASV_2218" ">ASV_2219" ">ASV_2220"
## [2221] ">ASV_2221" ">ASV_2222" ">ASV_2223" ">ASV_2224" ">ASV_2225" ">ASV_2226" ">ASV_2227" ">ASV_2228" ">ASV_2229" ">ASV_2230" ">ASV_2231" ">ASV_2232" ">ASV_2233" ">ASV_2234" ">ASV_2235" ">ASV_2236" ">ASV_2237" ">ASV_2238" ">ASV_2239" ">ASV_2240"
## [2241] ">ASV_2241" ">ASV_2242" ">ASV_2243" ">ASV_2244" ">ASV_2245" ">ASV_2246" ">ASV_2247" ">ASV_2248" ">ASV_2249" ">ASV_2250" ">ASV_2251" ">ASV_2252" ">ASV_2253" ">ASV_2254" ">ASV_2255" ">ASV_2256" ">ASV_2257" ">ASV_2258" ">ASV_2259" ">ASV_2260"
## [2261] ">ASV_2261" ">ASV_2262" ">ASV_2263" ">ASV_2264" ">ASV_2265" ">ASV_2266" ">ASV_2267" ">ASV_2268" ">ASV_2269" ">ASV_2270" ">ASV_2271" ">ASV_2272" ">ASV_2273" ">ASV_2274" ">ASV_2275" ">ASV_2276" ">ASV_2277" ">ASV_2278" ">ASV_2279" ">ASV_2280"
## [2281] ">ASV_2281" ">ASV_2282" ">ASV_2283" ">ASV_2284" ">ASV_2285" ">ASV_2286" ">ASV_2287" ">ASV_2288" ">ASV_2289" ">ASV_2290" ">ASV_2291" ">ASV_2292" ">ASV_2293" ">ASV_2294" ">ASV_2295" ">ASV_2296" ">ASV_2297" ">ASV_2298" ">ASV_2299" ">ASV_2300"
## [2301] ">ASV_2301" ">ASV_2302" ">ASV_2303" ">ASV_2304" ">ASV_2305" ">ASV_2306" ">ASV_2307" ">ASV_2308" ">ASV_2309" ">ASV_2310" ">ASV_2311" ">ASV_2312" ">ASV_2313" ">ASV_2314" ">ASV_2315" ">ASV_2316" ">ASV_2317" ">ASV_2318" ">ASV_2319" ">ASV_2320"
## [2321] ">ASV_2321" ">ASV_2322" ">ASV_2323" ">ASV_2324" ">ASV_2325" ">ASV_2326" ">ASV_2327" ">ASV_2328" ">ASV_2329" ">ASV_2330" ">ASV_2331" ">ASV_2332" ">ASV_2333" ">ASV_2334" ">ASV_2335" ">ASV_2336" ">ASV_2337" ">ASV_2338" ">ASV_2339" ">ASV_2340"
## [2341] ">ASV_2341" ">ASV_2342" ">ASV_2343" ">ASV_2344" ">ASV_2345" ">ASV_2346" ">ASV_2347" ">ASV_2348" ">ASV_2349" ">ASV_2350" ">ASV_2351" ">ASV_2352" ">ASV_2353" ">ASV_2354" ">ASV_2355" ">ASV_2356" ">ASV_2357" ">ASV_2358" ">ASV_2359" ">ASV_2360"
## [2361] ">ASV_2361" ">ASV_2362" ">ASV_2363" ">ASV_2364" ">ASV_2365" ">ASV_2366" ">ASV_2367" ">ASV_2368" ">ASV_2369" ">ASV_2370" ">ASV_2371" ">ASV_2372" ">ASV_2373" ">ASV_2374" ">ASV_2375" ">ASV_2376" ">ASV_2377" ">ASV_2378" ">ASV_2379" ">ASV_2380"
## [2381] ">ASV_2381" ">ASV_2382" ">ASV_2383" ">ASV_2384" ">ASV_2385" ">ASV_2386" ">ASV_2387" ">ASV_2388" ">ASV_2389" ">ASV_2390" ">ASV_2391" ">ASV_2392" ">ASV_2393" ">ASV_2394" ">ASV_2395" ">ASV_2396" ">ASV_2397" ">ASV_2398" ">ASV_2399" ">ASV_2400"
## [2401] ">ASV_2401" ">ASV_2402" ">ASV_2403" ">ASV_2404" ">ASV_2405" ">ASV_2406" ">ASV_2407" ">ASV_2408" ">ASV_2409" ">ASV_2410" ">ASV_2411" ">ASV_2412" ">ASV_2413" ">ASV_2414" ">ASV_2415" ">ASV_2416" ">ASV_2417" ">ASV_2418" ">ASV_2419" ">ASV_2420"
## [2421] ">ASV_2421" ">ASV_2422" ">ASV_2423" ">ASV_2424" ">ASV_2425" ">ASV_2426" ">ASV_2427" ">ASV_2428" ">ASV_2429" ">ASV_2430" ">ASV_2431" ">ASV_2432" ">ASV_2433" ">ASV_2434" ">ASV_2435" ">ASV_2436" ">ASV_2437" ">ASV_2438" ">ASV_2439" ">ASV_2440"
## [2441] ">ASV_2441" ">ASV_2442" ">ASV_2443" ">ASV_2444" ">ASV_2445" ">ASV_2446" ">ASV_2447" ">ASV_2448" ">ASV_2449" ">ASV_2450" ">ASV_2451" ">ASV_2452" ">ASV_2453" ">ASV_2454" ">ASV_2455" ">ASV_2456" ">ASV_2457" ">ASV_2458" ">ASV_2459" ">ASV_2460"
## [2461] ">ASV_2461" ">ASV_2462" ">ASV_2463" ">ASV_2464" ">ASV_2465" ">ASV_2466" ">ASV_2467" ">ASV_2468" ">ASV_2469" ">ASV_2470" ">ASV_2471" ">ASV_2472" ">ASV_2473" ">ASV_2474" ">ASV_2475" ">ASV_2476" ">ASV_2477" ">ASV_2478" ">ASV_2479" ">ASV_2480"
## [2481] ">ASV_2481" ">ASV_2482" ">ASV_2483" ">ASV_2484" ">ASV_2485" ">ASV_2486" ">ASV_2487" ">ASV_2488" ">ASV_2489" ">ASV_2490" ">ASV_2491" ">ASV_2492" ">ASV_2493" ">ASV_2494" ">ASV_2495" ">ASV_2496" ">ASV_2497" ">ASV_2498" ">ASV_2499" ">ASV_2500"
## [2501] ">ASV_2501" ">ASV_2502" ">ASV_2503" ">ASV_2504" ">ASV_2505" ">ASV_2506" ">ASV_2507" ">ASV_2508" ">ASV_2509" ">ASV_2510" ">ASV_2511" ">ASV_2512" ">ASV_2513" ">ASV_2514" ">ASV_2515" ">ASV_2516" ">ASV_2517" ">ASV_2518" ">ASV_2519" ">ASV_2520"
## [2521] ">ASV_2521" ">ASV_2522" ">ASV_2523" ">ASV_2524" ">ASV_2525" ">ASV_2526" ">ASV_2527" ">ASV_2528" ">ASV_2529" ">ASV_2530" ">ASV_2531" ">ASV_2532" ">ASV_2533" ">ASV_2534" ">ASV_2535" ">ASV_2536" ">ASV_2537" ">ASV_2538" ">ASV_2539" ">ASV_2540"
## [2541] ">ASV_2541" ">ASV_2542" ">ASV_2543" ">ASV_2544" ">ASV_2545" ">ASV_2546" ">ASV_2547" ">ASV_2548" ">ASV_2549" ">ASV_2550" ">ASV_2551" ">ASV_2552" ">ASV_2553" ">ASV_2554" ">ASV_2555" ">ASV_2556" ">ASV_2557" ">ASV_2558" ">ASV_2559" ">ASV_2560"
## [2561] ">ASV_2561" ">ASV_2562" ">ASV_2563" ">ASV_2564" ">ASV_2565" ">ASV_2566" ">ASV_2567" ">ASV_2568" ">ASV_2569" ">ASV_2570" ">ASV_2571" ">ASV_2572" ">ASV_2573" ">ASV_2574" ">ASV_2575" ">ASV_2576" ">ASV_2577" ">ASV_2578" ">ASV_2579" ">ASV_2580"
## [2581] ">ASV_2581" ">ASV_2582" ">ASV_2583" ">ASV_2584" ">ASV_2585" ">ASV_2586" ">ASV_2587" ">ASV_2588" ">ASV_2589" ">ASV_2590" ">ASV_2591" ">ASV_2592" ">ASV_2593" ">ASV_2594" ">ASV_2595" ">ASV_2596" ">ASV_2597" ">ASV_2598" ">ASV_2599" ">ASV_2600"
## [2601] ">ASV_2601" ">ASV_2602" ">ASV_2603" ">ASV_2604" ">ASV_2605" ">ASV_2606" ">ASV_2607" ">ASV_2608" ">ASV_2609" ">ASV_2610" ">ASV_2611" ">ASV_2612" ">ASV_2613" ">ASV_2614" ">ASV_2615" ">ASV_2616" ">ASV_2617" ">ASV_2618" ">ASV_2619" ">ASV_2620"
## [2621] ">ASV_2621" ">ASV_2622" ">ASV_2623" ">ASV_2624" ">ASV_2625" ">ASV_2626" ">ASV_2627" ">ASV_2628" ">ASV_2629" ">ASV_2630" ">ASV_2631" ">ASV_2632" ">ASV_2633" ">ASV_2634" ">ASV_2635" ">ASV_2636" ">ASV_2637" ">ASV_2638" ">ASV_2639" ">ASV_2640"
## [2641] ">ASV_2641" ">ASV_2642" ">ASV_2643" ">ASV_2644" ">ASV_2645" ">ASV_2646" ">ASV_2647" ">ASV_2648" ">ASV_2649" ">ASV_2650" ">ASV_2651" ">ASV_2652" ">ASV_2653" ">ASV_2654" ">ASV_2655" ">ASV_2656" ">ASV_2657" ">ASV_2658" ">ASV_2659" ">ASV_2660"
## [2661] ">ASV_2661" ">ASV_2662" ">ASV_2663" ">ASV_2664" ">ASV_2665" ">ASV_2666" ">ASV_2667" ">ASV_2668" ">ASV_2669" ">ASV_2670" ">ASV_2671" ">ASV_2672" ">ASV_2673" ">ASV_2674" ">ASV_2675" ">ASV_2676" ">ASV_2677" ">ASV_2678" ">ASV_2679" ">ASV_2680"
## [2681] ">ASV_2681" ">ASV_2682" ">ASV_2683" ">ASV_2684" ">ASV_2685" ">ASV_2686" ">ASV_2687" ">ASV_2688" ">ASV_2689" ">ASV_2690" ">ASV_2691" ">ASV_2692" ">ASV_2693" ">ASV_2694" ">ASV_2695" ">ASV_2696" ">ASV_2697" ">ASV_2698" ">ASV_2699" ">ASV_2700"
## [2701] ">ASV_2701" ">ASV_2702" ">ASV_2703" ">ASV_2704" ">ASV_2705" ">ASV_2706" ">ASV_2707" ">ASV_2708" ">ASV_2709" ">ASV_2710" ">ASV_2711" ">ASV_2712" ">ASV_2713" ">ASV_2714" ">ASV_2715" ">ASV_2716" ">ASV_2717" ">ASV_2718" ">ASV_2719" ">ASV_2720"
## [2721] ">ASV_2721" ">ASV_2722" ">ASV_2723" ">ASV_2724" ">ASV_2725" ">ASV_2726" ">ASV_2727" ">ASV_2728" ">ASV_2729" ">ASV_2730" ">ASV_2731" ">ASV_2732" ">ASV_2733" ">ASV_2734" ">ASV_2735" ">ASV_2736" ">ASV_2737" ">ASV_2738" ">ASV_2739" ">ASV_2740"
## [2741] ">ASV_2741" ">ASV_2742" ">ASV_2743" ">ASV_2744" ">ASV_2745" ">ASV_2746" ">ASV_2747" ">ASV_2748" ">ASV_2749" ">ASV_2750" ">ASV_2751" ">ASV_2752" ">ASV_2753" ">ASV_2754" ">ASV_2755" ">ASV_2756" ">ASV_2757" ">ASV_2758" ">ASV_2759" ">ASV_2760"
## [2761] ">ASV_2761" ">ASV_2762" ">ASV_2763" ">ASV_2764" ">ASV_2765" ">ASV_2766" ">ASV_2767" ">ASV_2768" ">ASV_2769" ">ASV_2770" ">ASV_2771" ">ASV_2772" ">ASV_2773" ">ASV_2774" ">ASV_2775" ">ASV_2776" ">ASV_2777" ">ASV_2778" ">ASV_2779" ">ASV_2780"
## [2781] ">ASV_2781" ">ASV_2782" ">ASV_2783" ">ASV_2784" ">ASV_2785" ">ASV_2786" ">ASV_2787" ">ASV_2788" ">ASV_2789" ">ASV_2790" ">ASV_2791" ">ASV_2792" ">ASV_2793" ">ASV_2794" ">ASV_2795" ">ASV_2796" ">ASV_2797" ">ASV_2798" ">ASV_2799" ">ASV_2800"
## [2801] ">ASV_2801" ">ASV_2802" ">ASV_2803" ">ASV_2804" ">ASV_2805" ">ASV_2806" ">ASV_2807" ">ASV_2808" ">ASV_2809" ">ASV_2810" ">ASV_2811" ">ASV_2812" ">ASV_2813" ">ASV_2814" ">ASV_2815" ">ASV_2816" ">ASV_2817" ">ASV_2818" ">ASV_2819" ">ASV_2820"
## [2821] ">ASV_2821" ">ASV_2822" ">ASV_2823" ">ASV_2824" ">ASV_2825" ">ASV_2826" ">ASV_2827" ">ASV_2828" ">ASV_2829" ">ASV_2830" ">ASV_2831" ">ASV_2832" ">ASV_2833" ">ASV_2834" ">ASV_2835" ">ASV_2836" ">ASV_2837" ">ASV_2838" ">ASV_2839" ">ASV_2840"
## [2841] ">ASV_2841" ">ASV_2842" ">ASV_2843" ">ASV_2844" ">ASV_2845" ">ASV_2846" ">ASV_2847" ">ASV_2848" ">ASV_2849" ">ASV_2850" ">ASV_2851" ">ASV_2852" ">ASV_2853" ">ASV_2854" ">ASV_2855" ">ASV_2856" ">ASV_2857" ">ASV_2858" ">ASV_2859" ">ASV_2860"
## [2861] ">ASV_2861" ">ASV_2862" ">ASV_2863" ">ASV_2864" ">ASV_2865" ">ASV_2866" ">ASV_2867" ">ASV_2868" ">ASV_2869" ">ASV_2870" ">ASV_2871" ">ASV_2872" ">ASV_2873" ">ASV_2874" ">ASV_2875" ">ASV_2876" ">ASV_2877" ">ASV_2878" ">ASV_2879" ">ASV_2880"
## [2881] ">ASV_2881" ">ASV_2882" ">ASV_2883" ">ASV_2884" ">ASV_2885" ">ASV_2886" ">ASV_2887" ">ASV_2888" ">ASV_2889" ">ASV_2890" ">ASV_2891" ">ASV_2892" ">ASV_2893" ">ASV_2894" ">ASV_2895" ">ASV_2896" ">ASV_2897" ">ASV_2898" ">ASV_2899" ">ASV_2900"
## [2901] ">ASV_2901" ">ASV_2902" ">ASV_2903" ">ASV_2904" ">ASV_2905" ">ASV_2906" ">ASV_2907" ">ASV_2908" ">ASV_2909" ">ASV_2910" ">ASV_2911" ">ASV_2912" ">ASV_2913" ">ASV_2914" ">ASV_2915" ">ASV_2916" ">ASV_2917" ">ASV_2918" ">ASV_2919" ">ASV_2920"
## [2921] ">ASV_2921" ">ASV_2922" ">ASV_2923" ">ASV_2924" ">ASV_2925" ">ASV_2926" ">ASV_2927" ">ASV_2928" ">ASV_2929" ">ASV_2930" ">ASV_2931" ">ASV_2932" ">ASV_2933" ">ASV_2934" ">ASV_2935" ">ASV_2936" ">ASV_2937" ">ASV_2938" ">ASV_2939" ">ASV_2940"
## [2941] ">ASV_2941" ">ASV_2942" ">ASV_2943" ">ASV_2944" ">ASV_2945" ">ASV_2946" ">ASV_2947" ">ASV_2948" ">ASV_2949" ">ASV_2950" ">ASV_2951" ">ASV_2952" ">ASV_2953" ">ASV_2954" ">ASV_2955" ">ASV_2956" ">ASV_2957" ">ASV_2958" ">ASV_2959" ">ASV_2960"
## [2961] ">ASV_2961" ">ASV_2962" ">ASV_2963" ">ASV_2964" ">ASV_2965" ">ASV_2966" ">ASV_2967" ">ASV_2968" ">ASV_2969" ">ASV_2970" ">ASV_2971" ">ASV_2972" ">ASV_2973" ">ASV_2974" ">ASV_2975" ">ASV_2976" ">ASV_2977" ">ASV_2978" ">ASV_2979" ">ASV_2980"
## [2981] ">ASV_2981" ">ASV_2982" ">ASV_2983" ">ASV_2984" ">ASV_2985" ">ASV_2986" ">ASV_2987" ">ASV_2988" ">ASV_2989" ">ASV_2990" ">ASV_2991" ">ASV_2992" ">ASV_2993" ">ASV_2994" ">ASV_2995" ">ASV_2996" ">ASV_2997" ">ASV_2998" ">ASV_2999" ">ASV_3000"
## [3001] ">ASV_3001" ">ASV_3002" ">ASV_3003" ">ASV_3004" ">ASV_3005" ">ASV_3006" ">ASV_3007" ">ASV_3008" ">ASV_3009" ">ASV_3010" ">ASV_3011" ">ASV_3012" ">ASV_3013" ">ASV_3014" ">ASV_3015" ">ASV_3016" ">ASV_3017" ">ASV_3018" ">ASV_3019" ">ASV_3020"
## [3021] ">ASV_3021" ">ASV_3022" ">ASV_3023" ">ASV_3024" ">ASV_3025" ">ASV_3026" ">ASV_3027" ">ASV_3028" ">ASV_3029" ">ASV_3030" ">ASV_3031" ">ASV_3032" ">ASV_3033" ">ASV_3034" ">ASV_3035" ">ASV_3036" ">ASV_3037" ">ASV_3038" ">ASV_3039" ">ASV_3040"
## [3041] ">ASV_3041" ">ASV_3042" ">ASV_3043" ">ASV_3044" ">ASV_3045" ">ASV_3046" ">ASV_3047" ">ASV_3048" ">ASV_3049" ">ASV_3050" ">ASV_3051" ">ASV_3052" ">ASV_3053" ">ASV_3054" ">ASV_3055" ">ASV_3056" ">ASV_3057" ">ASV_3058" ">ASV_3059" ">ASV_3060"
## [3061] ">ASV_3061" ">ASV_3062" ">ASV_3063" ">ASV_3064" ">ASV_3065" ">ASV_3066" ">ASV_3067" ">ASV_3068" ">ASV_3069" ">ASV_3070" ">ASV_3071" ">ASV_3072" ">ASV_3073" ">ASV_3074" ">ASV_3075" ">ASV_3076" ">ASV_3077" ">ASV_3078" ">ASV_3079" ">ASV_3080"
## [3081] ">ASV_3081" ">ASV_3082" ">ASV_3083" ">ASV_3084" ">ASV_3085" ">ASV_3086" ">ASV_3087" ">ASV_3088" ">ASV_3089" ">ASV_3090" ">ASV_3091" ">ASV_3092" ">ASV_3093" ">ASV_3094" ">ASV_3095" ">ASV_3096" ">ASV_3097" ">ASV_3098" ">ASV_3099" ">ASV_3100"
## [3101] ">ASV_3101" ">ASV_3102" ">ASV_3103" ">ASV_3104" ">ASV_3105" ">ASV_3106" ">ASV_3107" ">ASV_3108" ">ASV_3109" ">ASV_3110" ">ASV_3111" ">ASV_3112" ">ASV_3113" ">ASV_3114" ">ASV_3115" ">ASV_3116" ">ASV_3117" ">ASV_3118" ">ASV_3119" ">ASV_3120"
## [3121] ">ASV_3121" ">ASV_3122" ">ASV_3123" ">ASV_3124" ">ASV_3125" ">ASV_3126" ">ASV_3127" ">ASV_3128" ">ASV_3129" ">ASV_3130" ">ASV_3131" ">ASV_3132" ">ASV_3133" ">ASV_3134" ">ASV_3135" ">ASV_3136" ">ASV_3137" ">ASV_3138" ">ASV_3139" ">ASV_3140"
## [3141] ">ASV_3141" ">ASV_3142" ">ASV_3143" ">ASV_3144" ">ASV_3145" ">ASV_3146" ">ASV_3147" ">ASV_3148" ">ASV_3149" ">ASV_3150" ">ASV_3151" ">ASV_3152" ">ASV_3153" ">ASV_3154" ">ASV_3155" ">ASV_3156" ">ASV_3157" ">ASV_3158" ">ASV_3159" ">ASV_3160"
## [3161] ">ASV_3161" ">ASV_3162" ">ASV_3163" ">ASV_3164" ">ASV_3165" ">ASV_3166" ">ASV_3167" ">ASV_3168" ">ASV_3169" ">ASV_3170" ">ASV_3171" ">ASV_3172" ">ASV_3173" ">ASV_3174" ">ASV_3175" ">ASV_3176" ">ASV_3177" ">ASV_3178" ">ASV_3179" ">ASV_3180"
## [3181] ">ASV_3181" ">ASV_3182" ">ASV_3183" ">ASV_3184" ">ASV_3185" ">ASV_3186" ">ASV_3187" ">ASV_3188" ">ASV_3189" ">ASV_3190" ">ASV_3191" ">ASV_3192" ">ASV_3193" ">ASV_3194" ">ASV_3195" ">ASV_3196" ">ASV_3197" ">ASV_3198" ">ASV_3199" ">ASV_3200"
## [3201] ">ASV_3201" ">ASV_3202" ">ASV_3203" ">ASV_3204" ">ASV_3205" ">ASV_3206" ">ASV_3207" ">ASV_3208" ">ASV_3209" ">ASV_3210" ">ASV_3211" ">ASV_3212" ">ASV_3213" ">ASV_3214" ">ASV_3215" ">ASV_3216" ">ASV_3217" ">ASV_3218" ">ASV_3219" ">ASV_3220"
## [3221] ">ASV_3221" ">ASV_3222" ">ASV_3223" ">ASV_3224" ">ASV_3225" ">ASV_3226" ">ASV_3227" ">ASV_3228" ">ASV_3229" ">ASV_3230" ">ASV_3231" ">ASV_3232" ">ASV_3233" ">ASV_3234" ">ASV_3235" ">ASV_3236" ">ASV_3237" ">ASV_3238" ">ASV_3239" ">ASV_3240"
## [3241] ">ASV_3241" ">ASV_3242" ">ASV_3243" ">ASV_3244" ">ASV_3245" ">ASV_3246" ">ASV_3247" ">ASV_3248" ">ASV_3249" ">ASV_3250" ">ASV_3251" ">ASV_3252" ">ASV_3253" ">ASV_3254" ">ASV_3255" ">ASV_3256" ">ASV_3257" ">ASV_3258" ">ASV_3259" ">ASV_3260"
## [3261] ">ASV_3261" ">ASV_3262" ">ASV_3263" ">ASV_3264" ">ASV_3265" ">ASV_3266" ">ASV_3267" ">ASV_3268" ">ASV_3269" ">ASV_3270" ">ASV_3271" ">ASV_3272" ">ASV_3273" ">ASV_3274" ">ASV_3275" ">ASV_3276" ">ASV_3277" ">ASV_3278" ">ASV_3279" ">ASV_3280"
## [3281] ">ASV_3281" ">ASV_3282" ">ASV_3283" ">ASV_3284" ">ASV_3285" ">ASV_3286" ">ASV_3287" ">ASV_3288" ">ASV_3289" ">ASV_3290" ">ASV_3291" ">ASV_3292" ">ASV_3293" ">ASV_3294" ">ASV_3295" ">ASV_3296" ">ASV_3297" ">ASV_3298" ">ASV_3299" ">ASV_3300"
## [3301] ">ASV_3301" ">ASV_3302" ">ASV_3303" ">ASV_3304" ">ASV_3305" ">ASV_3306" ">ASV_3307" ">ASV_3308" ">ASV_3309" ">ASV_3310" ">ASV_3311" ">ASV_3312" ">ASV_3313" ">ASV_3314" ">ASV_3315" ">ASV_3316" ">ASV_3317" ">ASV_3318" ">ASV_3319" ">ASV_3320"
## [3321] ">ASV_3321" ">ASV_3322" ">ASV_3323" ">ASV_3324" ">ASV_3325" ">ASV_3326" ">ASV_3327" ">ASV_3328" ">ASV_3329" ">ASV_3330" ">ASV_3331" ">ASV_3332" ">ASV_3333" ">ASV_3334" ">ASV_3335" ">ASV_3336" ">ASV_3337" ">ASV_3338" ">ASV_3339" ">ASV_3340"
## [3341] ">ASV_3341" ">ASV_3342" ">ASV_3343" ">ASV_3344" ">ASV_3345" ">ASV_3346" ">ASV_3347" ">ASV_3348" ">ASV_3349" ">ASV_3350" ">ASV_3351" ">ASV_3352" ">ASV_3353" ">ASV_3354" ">ASV_3355" ">ASV_3356" ">ASV_3357" ">ASV_3358" ">ASV_3359" ">ASV_3360"
## [3361] ">ASV_3361" ">ASV_3362" ">ASV_3363" ">ASV_3364" ">ASV_3365" ">ASV_3366" ">ASV_3367" ">ASV_3368" ">ASV_3369" ">ASV_3370" ">ASV_3371" ">ASV_3372" ">ASV_3373" ">ASV_3374" ">ASV_3375" ">ASV_3376" ">ASV_3377" ">ASV_3378" ">ASV_3379" ">ASV_3380"
## [3381] ">ASV_3381" ">ASV_3382" ">ASV_3383" ">ASV_3384" ">ASV_3385" ">ASV_3386" ">ASV_3387" ">ASV_3388" ">ASV_3389" ">ASV_3390" ">ASV_3391" ">ASV_3392" ">ASV_3393" ">ASV_3394" ">ASV_3395" ">ASV_3396" ">ASV_3397" ">ASV_3398" ">ASV_3399" ">ASV_3400"
## [3401] ">ASV_3401" ">ASV_3402" ">ASV_3403" ">ASV_3404" ">ASV_3405" ">ASV_3406" ">ASV_3407" ">ASV_3408" ">ASV_3409" ">ASV_3410" ">ASV_3411" ">ASV_3412" ">ASV_3413" ">ASV_3414" ">ASV_3415" ">ASV_3416" ">ASV_3417" ">ASV_3418" ">ASV_3419" ">ASV_3420"
## [3421] ">ASV_3421" ">ASV_3422" ">ASV_3423" ">ASV_3424" ">ASV_3425" ">ASV_3426" ">ASV_3427" ">ASV_3428" ">ASV_3429" ">ASV_3430" ">ASV_3431" ">ASV_3432" ">ASV_3433" ">ASV_3434" ">ASV_3435" ">ASV_3436" ">ASV_3437" ">ASV_3438" ">ASV_3439" ">ASV_3440"
## [3441] ">ASV_3441" ">ASV_3442" ">ASV_3443" ">ASV_3444" ">ASV_3445" ">ASV_3446" ">ASV_3447" ">ASV_3448" ">ASV_3449" ">ASV_3450" ">ASV_3451" ">ASV_3452" ">ASV_3453" ">ASV_3454" ">ASV_3455" ">ASV_3456" ">ASV_3457" ">ASV_3458" ">ASV_3459" ">ASV_3460"
## [3461] ">ASV_3461" ">ASV_3462" ">ASV_3463" ">ASV_3464" ">ASV_3465" ">ASV_3466" ">ASV_3467" ">ASV_3468" ">ASV_3469" ">ASV_3470" ">ASV_3471" ">ASV_3472" ">ASV_3473" ">ASV_3474" ">ASV_3475" ">ASV_3476" ">ASV_3477" ">ASV_3478" ">ASV_3479" ">ASV_3480"
## [3481] ">ASV_3481" ">ASV_3482" ">ASV_3483" ">ASV_3484" ">ASV_3485" ">ASV_3486" ">ASV_3487" ">ASV_3488" ">ASV_3489" ">ASV_3490" ">ASV_3491" ">ASV_3492" ">ASV_3493" ">ASV_3494" ">ASV_3495" ">ASV_3496" ">ASV_3497" ">ASV_3498" ">ASV_3499" ">ASV_3500"
## [3501] ">ASV_3501" ">ASV_3502" ">ASV_3503" ">ASV_3504" ">ASV_3505" ">ASV_3506" ">ASV_3507" ">ASV_3508" ">ASV_3509" ">ASV_3510" ">ASV_3511" ">ASV_3512" ">ASV_3513" ">ASV_3514" ">ASV_3515" ">ASV_3516" ">ASV_3517" ">ASV_3518" ">ASV_3519" ">ASV_3520"
## [3521] ">ASV_3521" ">ASV_3522" ">ASV_3523" ">ASV_3524" ">ASV_3525" ">ASV_3526" ">ASV_3527" ">ASV_3528" ">ASV_3529" ">ASV_3530" ">ASV_3531" ">ASV_3532" ">ASV_3533" ">ASV_3534" ">ASV_3535" ">ASV_3536" ">ASV_3537" ">ASV_3538" ">ASV_3539" ">ASV_3540"
## [3541] ">ASV_3541" ">ASV_3542" ">ASV_3543" ">ASV_3544" ">ASV_3545" ">ASV_3546" ">ASV_3547" ">ASV_3548" ">ASV_3549" ">ASV_3550" ">ASV_3551" ">ASV_3552" ">ASV_3553" ">ASV_3554" ">ASV_3555" ">ASV_3556" ">ASV_3557" ">ASV_3558" ">ASV_3559" ">ASV_3560"
## [3561] ">ASV_3561" ">ASV_3562" ">ASV_3563" ">ASV_3564" ">ASV_3565" ">ASV_3566" ">ASV_3567" ">ASV_3568" ">ASV_3569" ">ASV_3570" ">ASV_3571" ">ASV_3572" ">ASV_3573" ">ASV_3574" ">ASV_3575" ">ASV_3576" ">ASV_3577" ">ASV_3578" ">ASV_3579" ">ASV_3580"
## [3581] ">ASV_3581" ">ASV_3582" ">ASV_3583" ">ASV_3584" ">ASV_3585" ">ASV_3586" ">ASV_3587" ">ASV_3588" ">ASV_3589" ">ASV_3590" ">ASV_3591" ">ASV_3592" ">ASV_3593" ">ASV_3594" ">ASV_3595" ">ASV_3596" ">ASV_3597" ">ASV_3598" ">ASV_3599" ">ASV_3600"
## [3601] ">ASV_3601" ">ASV_3602" ">ASV_3603" ">ASV_3604" ">ASV_3605" ">ASV_3606" ">ASV_3607" ">ASV_3608" ">ASV_3609" ">ASV_3610" ">ASV_3611" ">ASV_3612" ">ASV_3613" ">ASV_3614" ">ASV_3615" ">ASV_3616" ">ASV_3617" ">ASV_3618" ">ASV_3619" ">ASV_3620"
## [3621] ">ASV_3621" ">ASV_3622" ">ASV_3623" ">ASV_3624" ">ASV_3625" ">ASV_3626" ">ASV_3627" ">ASV_3628" ">ASV_3629" ">ASV_3630" ">ASV_3631" ">ASV_3632" ">ASV_3633" ">ASV_3634" ">ASV_3635" ">ASV_3636" ">ASV_3637" ">ASV_3638" ">ASV_3639" ">ASV_3640"
## [3641] ">ASV_3641" ">ASV_3642" ">ASV_3643" ">ASV_3644" ">ASV_3645" ">ASV_3646" ">ASV_3647" ">ASV_3648" ">ASV_3649" ">ASV_3650" ">ASV_3651" ">ASV_3652" ">ASV_3653" ">ASV_3654" ">ASV_3655" ">ASV_3656" ">ASV_3657" ">ASV_3658" ">ASV_3659" ">ASV_3660"
## [3661] ">ASV_3661" ">ASV_3662" ">ASV_3663" ">ASV_3664" ">ASV_3665" ">ASV_3666" ">ASV_3667" ">ASV_3668" ">ASV_3669" ">ASV_3670" ">ASV_3671" ">ASV_3672" ">ASV_3673" ">ASV_3674" ">ASV_3675" ">ASV_3676" ">ASV_3677" ">ASV_3678" ">ASV_3679" ">ASV_3680"
## [3681] ">ASV_3681" ">ASV_3682" ">ASV_3683" ">ASV_3684" ">ASV_3685" ">ASV_3686" ">ASV_3687" ">ASV_3688" ">ASV_3689" ">ASV_3690" ">ASV_3691" ">ASV_3692" ">ASV_3693" ">ASV_3694" ">ASV_3695" ">ASV_3696" ">ASV_3697" ">ASV_3698" ">ASV_3699" ">ASV_3700"
## [3701] ">ASV_3701" ">ASV_3702" ">ASV_3703" ">ASV_3704" ">ASV_3705" ">ASV_3706" ">ASV_3707" ">ASV_3708" ">ASV_3709" ">ASV_3710" ">ASV_3711" ">ASV_3712" ">ASV_3713" ">ASV_3714" ">ASV_3715" ">ASV_3716" ">ASV_3717" ">ASV_3718" ">ASV_3719" ">ASV_3720"
## [3721] ">ASV_3721" ">ASV_3722" ">ASV_3723" ">ASV_3724" ">ASV_3725" ">ASV_3726" ">ASV_3727" ">ASV_3728" ">ASV_3729" ">ASV_3730" ">ASV_3731" ">ASV_3732" ">ASV_3733" ">ASV_3734" ">ASV_3735" ">ASV_3736" ">ASV_3737" ">ASV_3738" ">ASV_3739" ">ASV_3740"
## [3741] ">ASV_3741" ">ASV_3742" ">ASV_3743" ">ASV_3744" ">ASV_3745" ">ASV_3746" ">ASV_3747" ">ASV_3748" ">ASV_3749" ">ASV_3750" ">ASV_3751" ">ASV_3752" ">ASV_3753" ">ASV_3754" ">ASV_3755" ">ASV_3756" ">ASV_3757" ">ASV_3758" ">ASV_3759" ">ASV_3760"
## [3761] ">ASV_3761" ">ASV_3762" ">ASV_3763" ">ASV_3764" ">ASV_3765" ">ASV_3766" ">ASV_3767" ">ASV_3768" ">ASV_3769" ">ASV_3770" ">ASV_3771" ">ASV_3772" ">ASV_3773" ">ASV_3774" ">ASV_3775" ">ASV_3776" ">ASV_3777" ">ASV_3778" ">ASV_3779" ">ASV_3780"
## [3781] ">ASV_3781" ">ASV_3782" ">ASV_3783" ">ASV_3784" ">ASV_3785" ">ASV_3786" ">ASV_3787" ">ASV_3788" ">ASV_3789" ">ASV_3790" ">ASV_3791" ">ASV_3792" ">ASV_3793" ">ASV_3794" ">ASV_3795" ">ASV_3796" ">ASV_3797" ">ASV_3798" ">ASV_3799" ">ASV_3800"
## [3801] ">ASV_3801" ">ASV_3802" ">ASV_3803" ">ASV_3804" ">ASV_3805" ">ASV_3806" ">ASV_3807" ">ASV_3808" ">ASV_3809" ">ASV_3810" ">ASV_3811" ">ASV_3812" ">ASV_3813" ">ASV_3814" ">ASV_3815" ">ASV_3816" ">ASV_3817" ">ASV_3818" ">ASV_3819" ">ASV_3820"
## [3821] ">ASV_3821" ">ASV_3822" ">ASV_3823" ">ASV_3824" ">ASV_3825" ">ASV_3826" ">ASV_3827" ">ASV_3828" ">ASV_3829" ">ASV_3830" ">ASV_3831" ">ASV_3832" ">ASV_3833" ">ASV_3834" ">ASV_3835" ">ASV_3836" ">ASV_3837" ">ASV_3838" ">ASV_3839" ">ASV_3840"
## [3841] ">ASV_3841" ">ASV_3842" ">ASV_3843" ">ASV_3844" ">ASV_3845" ">ASV_3846" ">ASV_3847" ">ASV_3848" ">ASV_3849" ">ASV_3850" ">ASV_3851" ">ASV_3852" ">ASV_3853" ">ASV_3854" ">ASV_3855" ">ASV_3856" ">ASV_3857" ">ASV_3858" ">ASV_3859" ">ASV_3860"
## [3861] ">ASV_3861" ">ASV_3862" ">ASV_3863" ">ASV_3864" ">ASV_3865" ">ASV_3866" ">ASV_3867" ">ASV_3868" ">ASV_3869" ">ASV_3870" ">ASV_3871" ">ASV_3872" ">ASV_3873" ">ASV_3874" ">ASV_3875" ">ASV_3876" ">ASV_3877" ">ASV_3878" ">ASV_3879" ">ASV_3880"
## [3881] ">ASV_3881" ">ASV_3882" ">ASV_3883" ">ASV_3884" ">ASV_3885" ">ASV_3886" ">ASV_3887" ">ASV_3888" ">ASV_3889" ">ASV_3890" ">ASV_3891" ">ASV_3892" ">ASV_3893" ">ASV_3894" ">ASV_3895" ">ASV_3896" ">ASV_3897" ">ASV_3898" ">ASV_3899" ">ASV_3900"
## [3901] ">ASV_3901" ">ASV_3902" ">ASV_3903" ">ASV_3904" ">ASV_3905" ">ASV_3906" ">ASV_3907" ">ASV_3908" ">ASV_3909" ">ASV_3910" ">ASV_3911" ">ASV_3912" ">ASV_3913" ">ASV_3914" ">ASV_3915" ">ASV_3916" ">ASV_3917" ">ASV_3918" ">ASV_3919" ">ASV_3920"
## [3921] ">ASV_3921" ">ASV_3922" ">ASV_3923" ">ASV_3924" ">ASV_3925" ">ASV_3926" ">ASV_3927" ">ASV_3928" ">ASV_3929" ">ASV_3930" ">ASV_3931" ">ASV_3932" ">ASV_3933" ">ASV_3934" ">ASV_3935" ">ASV_3936" ">ASV_3937" ">ASV_3938" ">ASV_3939" ">ASV_3940"
## [3941] ">ASV_3941" ">ASV_3942" ">ASV_3943" ">ASV_3944" ">ASV_3945" ">ASV_3946" ">ASV_3947" ">ASV_3948" ">ASV_3949" ">ASV_3950" ">ASV_3951" ">ASV_3952" ">ASV_3953" ">ASV_3954" ">ASV_3955" ">ASV_3956" ">ASV_3957" ">ASV_3958" ">ASV_3959" ">ASV_3960"
## [3961] ">ASV_3961" ">ASV_3962" ">ASV_3963" ">ASV_3964" ">ASV_3965" ">ASV_3966" ">ASV_3967" ">ASV_3968" ">ASV_3969" ">ASV_3970" ">ASV_3971" ">ASV_3972" ">ASV_3973" ">ASV_3974" ">ASV_3975" ">ASV_3976" ">ASV_3977" ">ASV_3978" ">ASV_3979" ">ASV_3980"
## [3981] ">ASV_3981" ">ASV_3982" ">ASV_3983" ">ASV_3984" ">ASV_3985" ">ASV_3986" ">ASV_3987" ">ASV_3988" ">ASV_3989" ">ASV_3990" ">ASV_3991" ">ASV_3992" ">ASV_3993" ">ASV_3994" ">ASV_3995" ">ASV_3996" ">ASV_3997" ">ASV_3998" ">ASV_3999" ">ASV_4000"
## [4001] ">ASV_4001" ">ASV_4002" ">ASV_4003" ">ASV_4004" ">ASV_4005" ">ASV_4006" ">ASV_4007" ">ASV_4008" ">ASV_4009" ">ASV_4010" ">ASV_4011" ">ASV_4012" ">ASV_4013" ">ASV_4014" ">ASV_4015" ">ASV_4016" ">ASV_4017" ">ASV_4018" ">ASV_4019" ">ASV_4020"
## [4021] ">ASV_4021" ">ASV_4022" ">ASV_4023" ">ASV_4024" ">ASV_4025" ">ASV_4026" ">ASV_4027" ">ASV_4028" ">ASV_4029" ">ASV_4030" ">ASV_4031" ">ASV_4032" ">ASV_4033" ">ASV_4034" ">ASV_4035" ">ASV_4036" ">ASV_4037" ">ASV_4038" ">ASV_4039" ">ASV_4040"
## [4041] ">ASV_4041" ">ASV_4042" ">ASV_4043" ">ASV_4044" ">ASV_4045" ">ASV_4046" ">ASV_4047" ">ASV_4048" ">ASV_4049" ">ASV_4050" ">ASV_4051" ">ASV_4052" ">ASV_4053" ">ASV_4054" ">ASV_4055" ">ASV_4056" ">ASV_4057" ">ASV_4058" ">ASV_4059" ">ASV_4060"
## [4061] ">ASV_4061" ">ASV_4062" ">ASV_4063" ">ASV_4064" ">ASV_4065" ">ASV_4066" ">ASV_4067" ">ASV_4068" ">ASV_4069" ">ASV_4070" ">ASV_4071" ">ASV_4072" ">ASV_4073" ">ASV_4074" ">ASV_4075" ">ASV_4076" ">ASV_4077" ">ASV_4078" ">ASV_4079" ">ASV_4080"
## [4081] ">ASV_4081" ">ASV_4082" ">ASV_4083" ">ASV_4084" ">ASV_4085" ">ASV_4086" ">ASV_4087" ">ASV_4088" ">ASV_4089" ">ASV_4090" ">ASV_4091" ">ASV_4092" ">ASV_4093" ">ASV_4094" ">ASV_4095" ">ASV_4096" ">ASV_4097" ">ASV_4098" ">ASV_4099" ">ASV_4100"
## [4101] ">ASV_4101" ">ASV_4102" ">ASV_4103" ">ASV_4104" ">ASV_4105" ">ASV_4106" ">ASV_4107" ">ASV_4108" ">ASV_4109" ">ASV_4110" ">ASV_4111" ">ASV_4112" ">ASV_4113" ">ASV_4114" ">ASV_4115" ">ASV_4116" ">ASV_4117" ">ASV_4118" ">ASV_4119" ">ASV_4120"
## [4121] ">ASV_4121" ">ASV_4122" ">ASV_4123" ">ASV_4124" ">ASV_4125" ">ASV_4126" ">ASV_4127" ">ASV_4128" ">ASV_4129" ">ASV_4130" ">ASV_4131" ">ASV_4132" ">ASV_4133" ">ASV_4134" ">ASV_4135" ">ASV_4136" ">ASV_4137" ">ASV_4138" ">ASV_4139" ">ASV_4140"
## [4141] ">ASV_4141" ">ASV_4142" ">ASV_4143" ">ASV_4144" ">ASV_4145" ">ASV_4146" ">ASV_4147" ">ASV_4148" ">ASV_4149" ">ASV_4150" ">ASV_4151" ">ASV_4152" ">ASV_4153" ">ASV_4154" ">ASV_4155" ">ASV_4156" ">ASV_4157" ">ASV_4158" ">ASV_4159" ">ASV_4160"
## [4161] ">ASV_4161" ">ASV_4162" ">ASV_4163" ">ASV_4164" ">ASV_4165" ">ASV_4166" ">ASV_4167" ">ASV_4168" ">ASV_4169" ">ASV_4170" ">ASV_4171" ">ASV_4172" ">ASV_4173" ">ASV_4174" ">ASV_4175" ">ASV_4176" ">ASV_4177" ">ASV_4178" ">ASV_4179" ">ASV_4180"
## [4181] ">ASV_4181" ">ASV_4182" ">ASV_4183" ">ASV_4184" ">ASV_4185" ">ASV_4186" ">ASV_4187" ">ASV_4188" ">ASV_4189" ">ASV_4190" ">ASV_4191" ">ASV_4192" ">ASV_4193" ">ASV_4194" ">ASV_4195" ">ASV_4196" ">ASV_4197" ">ASV_4198" ">ASV_4199" ">ASV_4200"
## [4201] ">ASV_4201" ">ASV_4202" ">ASV_4203" ">ASV_4204" ">ASV_4205" ">ASV_4206" ">ASV_4207" ">ASV_4208" ">ASV_4209" ">ASV_4210" ">ASV_4211" ">ASV_4212" ">ASV_4213" ">ASV_4214" ">ASV_4215" ">ASV_4216" ">ASV_4217" ">ASV_4218" ">ASV_4219" ">ASV_4220"
## [4221] ">ASV_4221" ">ASV_4222" ">ASV_4223" ">ASV_4224" ">ASV_4225" ">ASV_4226" ">ASV_4227" ">ASV_4228" ">ASV_4229" ">ASV_4230" ">ASV_4231" ">ASV_4232" ">ASV_4233" ">ASV_4234" ">ASV_4235" ">ASV_4236" ">ASV_4237" ">ASV_4238" ">ASV_4239" ">ASV_4240"
## [4241] ">ASV_4241" ">ASV_4242" ">ASV_4243" ">ASV_4244" ">ASV_4245" ">ASV_4246" ">ASV_4247" ">ASV_4248" ">ASV_4249" ">ASV_4250" ">ASV_4251" ">ASV_4252" ">ASV_4253" ">ASV_4254" ">ASV_4255" ">ASV_4256" ">ASV_4257" ">ASV_4258" ">ASV_4259" ">ASV_4260"
## [4261] ">ASV_4261" ">ASV_4262" ">ASV_4263" ">ASV_4264" ">ASV_4265" ">ASV_4266" ">ASV_4267" ">ASV_4268" ">ASV_4269" ">ASV_4270" ">ASV_4271" ">ASV_4272" ">ASV_4273" ">ASV_4274" ">ASV_4275" ">ASV_4276" ">ASV_4277" ">ASV_4278" ">ASV_4279" ">ASV_4280"
## [4281] ">ASV_4281" ">ASV_4282" ">ASV_4283" ">ASV_4284" ">ASV_4285" ">ASV_4286" ">ASV_4287" ">ASV_4288" ">ASV_4289" ">ASV_4290" ">ASV_4291" ">ASV_4292" ">ASV_4293" ">ASV_4294" ">ASV_4295" ">ASV_4296" ">ASV_4297" ">ASV_4298" ">ASV_4299" ">ASV_4300"
## [4301] ">ASV_4301" ">ASV_4302" ">ASV_4303" ">ASV_4304" ">ASV_4305" ">ASV_4306" ">ASV_4307" ">ASV_4308" ">ASV_4309" ">ASV_4310" ">ASV_4311" ">ASV_4312" ">ASV_4313" ">ASV_4314" ">ASV_4315" ">ASV_4316" ">ASV_4317" ">ASV_4318" ">ASV_4319" ">ASV_4320"
## [4321] ">ASV_4321" ">ASV_4322" ">ASV_4323" ">ASV_4324" ">ASV_4325" ">ASV_4326" ">ASV_4327" ">ASV_4328" ">ASV_4329" ">ASV_4330" ">ASV_4331" ">ASV_4332" ">ASV_4333" ">ASV_4334" ">ASV_4335" ">ASV_4336" ">ASV_4337" ">ASV_4338" ">ASV_4339" ">ASV_4340"
## [4341] ">ASV_4341" ">ASV_4342" ">ASV_4343" ">ASV_4344" ">ASV_4345" ">ASV_4346" ">ASV_4347" ">ASV_4348" ">ASV_4349" ">ASV_4350" ">ASV_4351" ">ASV_4352" ">ASV_4353" ">ASV_4354" ">ASV_4355" ">ASV_4356" ">ASV_4357" ">ASV_4358" ">ASV_4359" ">ASV_4360"
## [4361] ">ASV_4361" ">ASV_4362" ">ASV_4363" ">ASV_4364" ">ASV_4365" ">ASV_4366" ">ASV_4367" ">ASV_4368" ">ASV_4369" ">ASV_4370" ">ASV_4371" ">ASV_4372" ">ASV_4373" ">ASV_4374" ">ASV_4375" ">ASV_4376" ">ASV_4377" ">ASV_4378" ">ASV_4379" ">ASV_4380"
## [4381] ">ASV_4381" ">ASV_4382" ">ASV_4383" ">ASV_4384" ">ASV_4385" ">ASV_4386" ">ASV_4387" ">ASV_4388" ">ASV_4389" ">ASV_4390" ">ASV_4391" ">ASV_4392" ">ASV_4393" ">ASV_4394" ">ASV_4395" ">ASV_4396" ">ASV_4397" ">ASV_4398" ">ASV_4399" ">ASV_4400"
## [4401] ">ASV_4401" ">ASV_4402" ">ASV_4403" ">ASV_4404" ">ASV_4405" ">ASV_4406" ">ASV_4407" ">ASV_4408" ">ASV_4409" ">ASV_4410" ">ASV_4411" ">ASV_4412" ">ASV_4413" ">ASV_4414" ">ASV_4415" ">ASV_4416" ">ASV_4417" ">ASV_4418" ">ASV_4419" ">ASV_4420"
## [4421] ">ASV_4421" ">ASV_4422" ">ASV_4423" ">ASV_4424" ">ASV_4425" ">ASV_4426" ">ASV_4427" ">ASV_4428" ">ASV_4429" ">ASV_4430" ">ASV_4431" ">ASV_4432" ">ASV_4433" ">ASV_4434" ">ASV_4435" ">ASV_4436" ">ASV_4437" ">ASV_4438" ">ASV_4439" ">ASV_4440"
## [4441] ">ASV_4441" ">ASV_4442" ">ASV_4443" ">ASV_4444" ">ASV_4445" ">ASV_4446" ">ASV_4447" ">ASV_4448" ">ASV_4449" ">ASV_4450" ">ASV_4451" ">ASV_4452" ">ASV_4453" ">ASV_4454" ">ASV_4455" ">ASV_4456" ">ASV_4457" ">ASV_4458" ">ASV_4459" ">ASV_4460"
## [4461] ">ASV_4461" ">ASV_4462" ">ASV_4463" ">ASV_4464" ">ASV_4465" ">ASV_4466" ">ASV_4467" ">ASV_4468" ">ASV_4469" ">ASV_4470" ">ASV_4471" ">ASV_4472" ">ASV_4473" ">ASV_4474" ">ASV_4475" ">ASV_4476" ">ASV_4477" ">ASV_4478" ">ASV_4479" ">ASV_4480"
## [4481] ">ASV_4481" ">ASV_4482" ">ASV_4483" ">ASV_4484" ">ASV_4485" ">ASV_4486" ">ASV_4487" ">ASV_4488" ">ASV_4489" ">ASV_4490" ">ASV_4491" ">ASV_4492" ">ASV_4493" ">ASV_4494" ">ASV_4495" ">ASV_4496" ">ASV_4497" ">ASV_4498" ">ASV_4499" ">ASV_4500"
## [4501] ">ASV_4501" ">ASV_4502" ">ASV_4503" ">ASV_4504" ">ASV_4505" ">ASV_4506" ">ASV_4507" ">ASV_4508" ">ASV_4509" ">ASV_4510" ">ASV_4511" ">ASV_4512" ">ASV_4513" ">ASV_4514" ">ASV_4515" ">ASV_4516" ">ASV_4517" ">ASV_4518" ">ASV_4519" ">ASV_4520"
## [4521] ">ASV_4521" ">ASV_4522" ">ASV_4523" ">ASV_4524" ">ASV_4525" ">ASV_4526" ">ASV_4527" ">ASV_4528" ">ASV_4529" ">ASV_4530" ">ASV_4531" ">ASV_4532" ">ASV_4533" ">ASV_4534" ">ASV_4535" ">ASV_4536" ">ASV_4537" ">ASV_4538" ">ASV_4539" ">ASV_4540"
## [4541] ">ASV_4541" ">ASV_4542" ">ASV_4543" ">ASV_4544" ">ASV_4545" ">ASV_4546" ">ASV_4547" ">ASV_4548" ">ASV_4549" ">ASV_4550" ">ASV_4551" ">ASV_4552" ">ASV_4553" ">ASV_4554" ">ASV_4555" ">ASV_4556" ">ASV_4557" ">ASV_4558" ">ASV_4559" ">ASV_4560"
## [4561] ">ASV_4561" ">ASV_4562" ">ASV_4563" ">ASV_4564" ">ASV_4565" ">ASV_4566" ">ASV_4567" ">ASV_4568" ">ASV_4569" ">ASV_4570" ">ASV_4571" ">ASV_4572" ">ASV_4573" ">ASV_4574" ">ASV_4575" ">ASV_4576" ">ASV_4577" ">ASV_4578" ">ASV_4579" ">ASV_4580"
## [4581] ">ASV_4581" ">ASV_4582" ">ASV_4583" ">ASV_4584" ">ASV_4585" ">ASV_4586" ">ASV_4587" ">ASV_4588" ">ASV_4589" ">ASV_4590" ">ASV_4591" ">ASV_4592" ">ASV_4593" ">ASV_4594" ">ASV_4595" ">ASV_4596" ">ASV_4597" ">ASV_4598" ">ASV_4599" ">ASV_4600"
## [4601] ">ASV_4601" ">ASV_4602" ">ASV_4603" ">ASV_4604" ">ASV_4605" ">ASV_4606" ">ASV_4607" ">ASV_4608" ">ASV_4609" ">ASV_4610" ">ASV_4611" ">ASV_4612" ">ASV_4613" ">ASV_4614" ">ASV_4615" ">ASV_4616" ">ASV_4617" ">ASV_4618" ">ASV_4619" ">ASV_4620"
## [4621] ">ASV_4621" ">ASV_4622" ">ASV_4623" ">ASV_4624" ">ASV_4625" ">ASV_4626" ">ASV_4627" ">ASV_4628" ">ASV_4629" ">ASV_4630" ">ASV_4631" ">ASV_4632" ">ASV_4633" ">ASV_4634" ">ASV_4635" ">ASV_4636" ">ASV_4637" ">ASV_4638" ">ASV_4639" ">ASV_4640"
## [4641] ">ASV_4641" ">ASV_4642" ">ASV_4643" ">ASV_4644" ">ASV_4645" ">ASV_4646" ">ASV_4647" ">ASV_4648" ">ASV_4649" ">ASV_4650" ">ASV_4651" ">ASV_4652" ">ASV_4653" ">ASV_4654" ">ASV_4655" ">ASV_4656" ">ASV_4657" ">ASV_4658" ">ASV_4659" ">ASV_4660"
## [4661] ">ASV_4661" ">ASV_4662" ">ASV_4663" ">ASV_4664" ">ASV_4665" ">ASV_4666" ">ASV_4667" ">ASV_4668" ">ASV_4669" ">ASV_4670" ">ASV_4671" ">ASV_4672" ">ASV_4673" ">ASV_4674" ">ASV_4675" ">ASV_4676" ">ASV_4677" ">ASV_4678" ">ASV_4679" ">ASV_4680"
## [4681] ">ASV_4681" ">ASV_4682" ">ASV_4683" ">ASV_4684" ">ASV_4685" ">ASV_4686" ">ASV_4687" ">ASV_4688" ">ASV_4689" ">ASV_4690" ">ASV_4691" ">ASV_4692" ">ASV_4693" ">ASV_4694" ">ASV_4695" ">ASV_4696" ">ASV_4697" ">ASV_4698" ">ASV_4699" ">ASV_4700"
## [4701] ">ASV_4701" ">ASV_4702" ">ASV_4703" ">ASV_4704" ">ASV_4705" ">ASV_4706" ">ASV_4707" ">ASV_4708" ">ASV_4709" ">ASV_4710" ">ASV_4711" ">ASV_4712" ">ASV_4713" ">ASV_4714" ">ASV_4715" ">ASV_4716" ">ASV_4717" ">ASV_4718" ">ASV_4719" ">ASV_4720"
## [4721] ">ASV_4721" ">ASV_4722" ">ASV_4723" ">ASV_4724" ">ASV_4725" ">ASV_4726" ">ASV_4727" ">ASV_4728" ">ASV_4729" ">ASV_4730" ">ASV_4731" ">ASV_4732" ">ASV_4733" ">ASV_4734" ">ASV_4735" ">ASV_4736" ">ASV_4737" ">ASV_4738" ">ASV_4739" ">ASV_4740"
## [4741] ">ASV_4741" ">ASV_4742" ">ASV_4743" ">ASV_4744" ">ASV_4745" ">ASV_4746" ">ASV_4747" ">ASV_4748" ">ASV_4749" ">ASV_4750" ">ASV_4751" ">ASV_4752" ">ASV_4753" ">ASV_4754" ">ASV_4755" ">ASV_4756" ">ASV_4757" ">ASV_4758" ">ASV_4759" ">ASV_4760"
## [4761] ">ASV_4761" ">ASV_4762" ">ASV_4763" ">ASV_4764" ">ASV_4765" ">ASV_4766" ">ASV_4767" ">ASV_4768" ">ASV_4769" ">ASV_4770" ">ASV_4771" ">ASV_4772" ">ASV_4773" ">ASV_4774" ">ASV_4775" ">ASV_4776" ">ASV_4777" ">ASV_4778" ">ASV_4779" ">ASV_4780"
## [4781] ">ASV_4781" ">ASV_4782" ">ASV_4783" ">ASV_4784" ">ASV_4785" ">ASV_4786" ">ASV_4787" ">ASV_4788" ">ASV_4789" ">ASV_4790" ">ASV_4791" ">ASV_4792" ">ASV_4793" ">ASV_4794" ">ASV_4795" ">ASV_4796" ">ASV_4797" ">ASV_4798" ">ASV_4799" ">ASV_4800"
## [4801] ">ASV_4801" ">ASV_4802" ">ASV_4803" ">ASV_4804" ">ASV_4805" ">ASV_4806" ">ASV_4807" ">ASV_4808" ">ASV_4809" ">ASV_4810" ">ASV_4811" ">ASV_4812" ">ASV_4813" ">ASV_4814" ">ASV_4815" ">ASV_4816" ">ASV_4817" ">ASV_4818" ">ASV_4819" ">ASV_4820"
## [4821] ">ASV_4821" ">ASV_4822" ">ASV_4823" ">ASV_4824" ">ASV_4825" ">ASV_4826" ">ASV_4827" ">ASV_4828" ">ASV_4829" ">ASV_4830" ">ASV_4831" ">ASV_4832" ">ASV_4833" ">ASV_4834" ">ASV_4835" ">ASV_4836" ">ASV_4837" ">ASV_4838" ">ASV_4839" ">ASV_4840"
## [4841] ">ASV_4841" ">ASV_4842" ">ASV_4843" ">ASV_4844" ">ASV_4845" ">ASV_4846" ">ASV_4847" ">ASV_4848" ">ASV_4849" ">ASV_4850" ">ASV_4851" ">ASV_4852" ">ASV_4853" ">ASV_4854" ">ASV_4855" ">ASV_4856" ">ASV_4857" ">ASV_4858" ">ASV_4859" ">ASV_4860"
## [4861] ">ASV_4861" ">ASV_4862" ">ASV_4863" ">ASV_4864" ">ASV_4865" ">ASV_4866" ">ASV_4867" ">ASV_4868" ">ASV_4869" ">ASV_4870" ">ASV_4871" ">ASV_4872" ">ASV_4873" ">ASV_4874" ">ASV_4875" ">ASV_4876" ">ASV_4877" ">ASV_4878" ">ASV_4879" ">ASV_4880"
## [4881] ">ASV_4881" ">ASV_4882" ">ASV_4883" ">ASV_4884" ">ASV_4885" ">ASV_4886" ">ASV_4887" ">ASV_4888" ">ASV_4889" ">ASV_4890" ">ASV_4891" ">ASV_4892" ">ASV_4893" ">ASV_4894" ">ASV_4895" ">ASV_4896" ">ASV_4897" ">ASV_4898" ">ASV_4899" ">ASV_4900"
## [4901] ">ASV_4901" ">ASV_4902" ">ASV_4903" ">ASV_4904" ">ASV_4905" ">ASV_4906" ">ASV_4907" ">ASV_4908" ">ASV_4909" ">ASV_4910" ">ASV_4911" ">ASV_4912" ">ASV_4913" ">ASV_4914" ">ASV_4915" ">ASV_4916" ">ASV_4917" ">ASV_4918" ">ASV_4919" ">ASV_4920"
## [4921] ">ASV_4921" ">ASV_4922" ">ASV_4923" ">ASV_4924" ">ASV_4925" ">ASV_4926" ">ASV_4927" ">ASV_4928" ">ASV_4929" ">ASV_4930" ">ASV_4931" ">ASV_4932" ">ASV_4933" ">ASV_4934" ">ASV_4935" ">ASV_4936" ">ASV_4937" ">ASV_4938" ">ASV_4939" ">ASV_4940"
## [4941] ">ASV_4941" ">ASV_4942" ">ASV_4943" ">ASV_4944" ">ASV_4945" ">ASV_4946" ">ASV_4947" ">ASV_4948" ">ASV_4949" ">ASV_4950" ">ASV_4951" ">ASV_4952" ">ASV_4953" ">ASV_4954" ">ASV_4955" ">ASV_4956" ">ASV_4957" ">ASV_4958" ">ASV_4959" ">ASV_4960"
## [4961] ">ASV_4961" ">ASV_4962" ">ASV_4963" ">ASV_4964" ">ASV_4965" ">ASV_4966" ">ASV_4967" ">ASV_4968" ">ASV_4969" ">ASV_4970" ">ASV_4971" ">ASV_4972" ">ASV_4973" ">ASV_4974" ">ASV_4975" ">ASV_4976" ">ASV_4977" ">ASV_4978" ">ASV_4979" ">ASV_4980"
## [4981] ">ASV_4981" ">ASV_4982" ">ASV_4983" ">ASV_4984" ">ASV_4985" ">ASV_4986" ">ASV_4987" ">ASV_4988" ">ASV_4989" ">ASV_4990" ">ASV_4991" ">ASV_4992" ">ASV_4993" ">ASV_4994" ">ASV_4995" ">ASV_4996" ">ASV_4997" ">ASV_4998" ">ASV_4999" ">ASV_5000"
## [5001] ">ASV_5001" ">ASV_5002" ">ASV_5003" ">ASV_5004" ">ASV_5005" ">ASV_5006" ">ASV_5007" ">ASV_5008" ">ASV_5009" ">ASV_5010" ">ASV_5011" ">ASV_5012" ">ASV_5013" ">ASV_5014" ">ASV_5015" ">ASV_5016" ">ASV_5017" ">ASV_5018" ">ASV_5019" ">ASV_5020"
## [5021] ">ASV_5021" ">ASV_5022" ">ASV_5023" ">ASV_5024" ">ASV_5025" ">ASV_5026" ">ASV_5027" ">ASV_5028" ">ASV_5029" ">ASV_5030" ">ASV_5031" ">ASV_5032" ">ASV_5033" ">ASV_5034" ">ASV_5035" ">ASV_5036" ">ASV_5037" ">ASV_5038" ">ASV_5039" ">ASV_5040"
## [5041] ">ASV_5041" ">ASV_5042" ">ASV_5043" ">ASV_5044" ">ASV_5045" ">ASV_5046" ">ASV_5047" ">ASV_5048" ">ASV_5049" ">ASV_5050" ">ASV_5051" ">ASV_5052" ">ASV_5053" ">ASV_5054" ">ASV_5055" ">ASV_5056" ">ASV_5057" ">ASV_5058" ">ASV_5059" ">ASV_5060"
## [5061] ">ASV_5061" ">ASV_5062" ">ASV_5063" ">ASV_5064" ">ASV_5065" ">ASV_5066" ">ASV_5067" ">ASV_5068" ">ASV_5069" ">ASV_5070" ">ASV_5071" ">ASV_5072" ">ASV_5073" ">ASV_5074" ">ASV_5075" ">ASV_5076" ">ASV_5077" ">ASV_5078" ">ASV_5079" ">ASV_5080"
## [5081] ">ASV_5081" ">ASV_5082" ">ASV_5083" ">ASV_5084" ">ASV_5085" ">ASV_5086" ">ASV_5087" ">ASV_5088" ">ASV_5089" ">ASV_5090" ">ASV_5091" ">ASV_5092" ">ASV_5093" ">ASV_5094" ">ASV_5095" ">ASV_5096" ">ASV_5097" ">ASV_5098" ">ASV_5099" ">ASV_5100"
## [5101] ">ASV_5101" ">ASV_5102" ">ASV_5103" ">ASV_5104" ">ASV_5105" ">ASV_5106" ">ASV_5107" ">ASV_5108" ">ASV_5109" ">ASV_5110" ">ASV_5111" ">ASV_5112" ">ASV_5113" ">ASV_5114" ">ASV_5115" ">ASV_5116" ">ASV_5117" ">ASV_5118" ">ASV_5119" ">ASV_5120"
## [5121] ">ASV_5121" ">ASV_5122" ">ASV_5123" ">ASV_5124" ">ASV_5125" ">ASV_5126" ">ASV_5127" ">ASV_5128" ">ASV_5129" ">ASV_5130" ">ASV_5131" ">ASV_5132" ">ASV_5133" ">ASV_5134" ">ASV_5135" ">ASV_5136" ">ASV_5137" ">ASV_5138" ">ASV_5139" ">ASV_5140"
## [5141] ">ASV_5141" ">ASV_5142" ">ASV_5143" ">ASV_5144" ">ASV_5145" ">ASV_5146" ">ASV_5147" ">ASV_5148" ">ASV_5149" ">ASV_5150" ">ASV_5151" ">ASV_5152" ">ASV_5153" ">ASV_5154" ">ASV_5155" ">ASV_5156" ">ASV_5157" ">ASV_5158" ">ASV_5159" ">ASV_5160"
## [5161] ">ASV_5161" ">ASV_5162" ">ASV_5163" ">ASV_5164" ">ASV_5165" ">ASV_5166" ">ASV_5167" ">ASV_5168" ">ASV_5169" ">ASV_5170" ">ASV_5171" ">ASV_5172" ">ASV_5173" ">ASV_5174" ">ASV_5175" ">ASV_5176" ">ASV_5177" ">ASV_5178" ">ASV_5179" ">ASV_5180"
## [5181] ">ASV_5181" ">ASV_5182" ">ASV_5183" ">ASV_5184" ">ASV_5185" ">ASV_5186" ">ASV_5187" ">ASV_5188" ">ASV_5189" ">ASV_5190" ">ASV_5191" ">ASV_5192" ">ASV_5193" ">ASV_5194" ">ASV_5195" ">ASV_5196" ">ASV_5197" ">ASV_5198" ">ASV_5199" ">ASV_5200"
## [5201] ">ASV_5201" ">ASV_5202" ">ASV_5203" ">ASV_5204" ">ASV_5205" ">ASV_5206" ">ASV_5207" ">ASV_5208" ">ASV_5209" ">ASV_5210" ">ASV_5211" ">ASV_5212" ">ASV_5213" ">ASV_5214" ">ASV_5215" ">ASV_5216" ">ASV_5217" ">ASV_5218" ">ASV_5219" ">ASV_5220"
## [5221] ">ASV_5221" ">ASV_5222" ">ASV_5223" ">ASV_5224" ">ASV_5225" ">ASV_5226" ">ASV_5227" ">ASV_5228" ">ASV_5229" ">ASV_5230" ">ASV_5231" ">ASV_5232" ">ASV_5233" ">ASV_5234" ">ASV_5235" ">ASV_5236" ">ASV_5237" ">ASV_5238" ">ASV_5239" ">ASV_5240"
## [5241] ">ASV_5241" ">ASV_5242" ">ASV_5243" ">ASV_5244" ">ASV_5245" ">ASV_5246" ">ASV_5247" ">ASV_5248" ">ASV_5249" ">ASV_5250" ">ASV_5251" ">ASV_5252" ">ASV_5253" ">ASV_5254" ">ASV_5255" ">ASV_5256" ">ASV_5257" ">ASV_5258" ">ASV_5259" ">ASV_5260"
## [5261] ">ASV_5261" ">ASV_5262" ">ASV_5263" ">ASV_5264" ">ASV_5265" ">ASV_5266" ">ASV_5267" ">ASV_5268" ">ASV_5269" ">ASV_5270" ">ASV_5271" ">ASV_5272" ">ASV_5273" ">ASV_5274" ">ASV_5275" ">ASV_5276" ">ASV_5277" ">ASV_5278" ">ASV_5279" ">ASV_5280"
## [5281] ">ASV_5281" ">ASV_5282" ">ASV_5283" ">ASV_5284" ">ASV_5285" ">ASV_5286" ">ASV_5287" ">ASV_5288" ">ASV_5289" ">ASV_5290" ">ASV_5291" ">ASV_5292" ">ASV_5293" ">ASV_5294" ">ASV_5295" ">ASV_5296" ">ASV_5297" ">ASV_5298" ">ASV_5299" ">ASV_5300"
## [5301] ">ASV_5301" ">ASV_5302" ">ASV_5303" ">ASV_5304" ">ASV_5305" ">ASV_5306" ">ASV_5307" ">ASV_5308" ">ASV_5309" ">ASV_5310" ">ASV_5311" ">ASV_5312" ">ASV_5313" ">ASV_5314" ">ASV_5315" ">ASV_5316" ">ASV_5317" ">ASV_5318" ">ASV_5319" ">ASV_5320"
## [5321] ">ASV_5321" ">ASV_5322" ">ASV_5323" ">ASV_5324" ">ASV_5325" ">ASV_5326" ">ASV_5327" ">ASV_5328" ">ASV_5329" ">ASV_5330" ">ASV_5331" ">ASV_5332" ">ASV_5333" ">ASV_5334" ">ASV_5335" ">ASV_5336" ">ASV_5337" ">ASV_5338" ">ASV_5339" ">ASV_5340"
## [5341] ">ASV_5341" ">ASV_5342" ">ASV_5343" ">ASV_5344" ">ASV_5345" ">ASV_5346" ">ASV_5347" ">ASV_5348" ">ASV_5349" ">ASV_5350" ">ASV_5351" ">ASV_5352" ">ASV_5353" ">ASV_5354" ">ASV_5355" ">ASV_5356" ">ASV_5357" ">ASV_5358" ">ASV_5359" ">ASV_5360"
## [5361] ">ASV_5361" ">ASV_5362" ">ASV_5363" ">ASV_5364" ">ASV_5365" ">ASV_5366" ">ASV_5367" ">ASV_5368" ">ASV_5369" ">ASV_5370" ">ASV_5371" ">ASV_5372" ">ASV_5373" ">ASV_5374" ">ASV_5375" ">ASV_5376" ">ASV_5377" ">ASV_5378" ">ASV_5379" ">ASV_5380"
## [5381] ">ASV_5381" ">ASV_5382" ">ASV_5383" ">ASV_5384" ">ASV_5385" ">ASV_5386" ">ASV_5387" ">ASV_5388" ">ASV_5389" ">ASV_5390" ">ASV_5391" ">ASV_5392" ">ASV_5393" ">ASV_5394" ">ASV_5395" ">ASV_5396" ">ASV_5397" ">ASV_5398" ">ASV_5399" ">ASV_5400"
## [5401] ">ASV_5401" ">ASV_5402" ">ASV_5403" ">ASV_5404" ">ASV_5405" ">ASV_5406" ">ASV_5407" ">ASV_5408" ">ASV_5409" ">ASV_5410" ">ASV_5411" ">ASV_5412" ">ASV_5413" ">ASV_5414" ">ASV_5415" ">ASV_5416" ">ASV_5417" ">ASV_5418" ">ASV_5419" ">ASV_5420"
## [5421] ">ASV_5421" ">ASV_5422" ">ASV_5423" ">ASV_5424" ">ASV_5425" ">ASV_5426" ">ASV_5427" ">ASV_5428" ">ASV_5429" ">ASV_5430" ">ASV_5431" ">ASV_5432" ">ASV_5433" ">ASV_5434" ">ASV_5435" ">ASV_5436" ">ASV_5437" ">ASV_5438" ">ASV_5439" ">ASV_5440"
## [5441] ">ASV_5441" ">ASV_5442" ">ASV_5443" ">ASV_5444" ">ASV_5445" ">ASV_5446" ">ASV_5447" ">ASV_5448" ">ASV_5449" ">ASV_5450" ">ASV_5451" ">ASV_5452" ">ASV_5453" ">ASV_5454" ">ASV_5455" ">ASV_5456" ">ASV_5457" ">ASV_5458" ">ASV_5459" ">ASV_5460"
## [5461] ">ASV_5461" ">ASV_5462" ">ASV_5463" ">ASV_5464" ">ASV_5465" ">ASV_5466" ">ASV_5467" ">ASV_5468" ">ASV_5469" ">ASV_5470" ">ASV_5471" ">ASV_5472" ">ASV_5473" ">ASV_5474" ">ASV_5475" ">ASV_5476" ">ASV_5477" ">ASV_5478" ">ASV_5479" ">ASV_5480"
## [5481] ">ASV_5481" ">ASV_5482" ">ASV_5483" ">ASV_5484" ">ASV_5485" ">ASV_5486" ">ASV_5487" ">ASV_5488" ">ASV_5489" ">ASV_5490" ">ASV_5491" ">ASV_5492" ">ASV_5493" ">ASV_5494" ">ASV_5495" ">ASV_5496" ">ASV_5497" ">ASV_5498" ">ASV_5499" ">ASV_5500"
## [5501] ">ASV_5501" ">ASV_5502" ">ASV_5503" ">ASV_5504" ">ASV_5505" ">ASV_5506" ">ASV_5507" ">ASV_5508" ">ASV_5509" ">ASV_5510" ">ASV_5511" ">ASV_5512" ">ASV_5513" ">ASV_5514" ">ASV_5515" ">ASV_5516" ">ASV_5517" ">ASV_5518" ">ASV_5519" ">ASV_5520"
## [5521] ">ASV_5521" ">ASV_5522" ">ASV_5523" ">ASV_5524" ">ASV_5525" ">ASV_5526" ">ASV_5527" ">ASV_5528" ">ASV_5529" ">ASV_5530" ">ASV_5531" ">ASV_5532" ">ASV_5533" ">ASV_5534" ">ASV_5535" ">ASV_5536" ">ASV_5537" ">ASV_5538" ">ASV_5539" ">ASV_5540"
## [5541] ">ASV_5541" ">ASV_5542" ">ASV_5543" ">ASV_5544" ">ASV_5545" ">ASV_5546" ">ASV_5547" ">ASV_5548" ">ASV_5549" ">ASV_5550" ">ASV_5551" ">ASV_5552" ">ASV_5553" ">ASV_5554" ">ASV_5555" ">ASV_5556" ">ASV_5557" ">ASV_5558" ">ASV_5559" ">ASV_5560"
## [5561] ">ASV_5561" ">ASV_5562" ">ASV_5563" ">ASV_5564" ">ASV_5565" ">ASV_5566" ">ASV_5567" ">ASV_5568" ">ASV_5569" ">ASV_5570" ">ASV_5571" ">ASV_5572" ">ASV_5573" ">ASV_5574" ">ASV_5575" ">ASV_5576" ">ASV_5577" ">ASV_5578" ">ASV_5579" ">ASV_5580"
## [5581] ">ASV_5581" ">ASV_5582" ">ASV_5583" ">ASV_5584" ">ASV_5585" ">ASV_5586" ">ASV_5587" ">ASV_5588" ">ASV_5589" ">ASV_5590" ">ASV_5591" ">ASV_5592" ">ASV_5593" ">ASV_5594" ">ASV_5595" ">ASV_5596" ">ASV_5597" ">ASV_5598" ">ASV_5599" ">ASV_5600"
## [5601] ">ASV_5601" ">ASV_5602" ">ASV_5603" ">ASV_5604" ">ASV_5605" ">ASV_5606" ">ASV_5607" ">ASV_5608" ">ASV_5609" ">ASV_5610" ">ASV_5611" ">ASV_5612" ">ASV_5613" ">ASV_5614" ">ASV_5615" ">ASV_5616" ">ASV_5617" ">ASV_5618" ">ASV_5619" ">ASV_5620"
## [5621] ">ASV_5621" ">ASV_5622" ">ASV_5623" ">ASV_5624" ">ASV_5625" ">ASV_5626" ">ASV_5627" ">ASV_5628" ">ASV_5629" ">ASV_5630" ">ASV_5631" ">ASV_5632" ">ASV_5633" ">ASV_5634" ">ASV_5635" ">ASV_5636" ">ASV_5637" ">ASV_5638" ">ASV_5639" ">ASV_5640"
## [5641] ">ASV_5641" ">ASV_5642" ">ASV_5643" ">ASV_5644" ">ASV_5645" ">ASV_5646" ">ASV_5647" ">ASV_5648" ">ASV_5649" ">ASV_5650" ">ASV_5651" ">ASV_5652" ">ASV_5653" ">ASV_5654" ">ASV_5655" ">ASV_5656" ">ASV_5657" ">ASV_5658" ">ASV_5659" ">ASV_5660"
## [5661] ">ASV_5661" ">ASV_5662" ">ASV_5663" ">ASV_5664" ">ASV_5665" ">ASV_5666" ">ASV_5667" ">ASV_5668" ">ASV_5669" ">ASV_5670" ">ASV_5671" ">ASV_5672" ">ASV_5673" ">ASV_5674" ">ASV_5675" ">ASV_5676" ">ASV_5677" ">ASV_5678" ">ASV_5679" ">ASV_5680"
## [5681] ">ASV_5681" ">ASV_5682" ">ASV_5683" ">ASV_5684" ">ASV_5685" ">ASV_5686" ">ASV_5687" ">ASV_5688" ">ASV_5689" ">ASV_5690" ">ASV_5691" ">ASV_5692" ">ASV_5693" ">ASV_5694" ">ASV_5695" ">ASV_5696" ">ASV_5697" ">ASV_5698" ">ASV_5699" ">ASV_5700"
## [5701] ">ASV_5701" ">ASV_5702" ">ASV_5703" ">ASV_5704" ">ASV_5705" ">ASV_5706" ">ASV_5707" ">ASV_5708" ">ASV_5709" ">ASV_5710" ">ASV_5711" ">ASV_5712" ">ASV_5713" ">ASV_5714" ">ASV_5715" ">ASV_5716" ">ASV_5717" ">ASV_5718" ">ASV_5719" ">ASV_5720"
## [5721] ">ASV_5721" ">ASV_5722" ">ASV_5723" ">ASV_5724" ">ASV_5725" ">ASV_5726" ">ASV_5727" ">ASV_5728" ">ASV_5729" ">ASV_5730" ">ASV_5731" ">ASV_5732" ">ASV_5733" ">ASV_5734" ">ASV_5735" ">ASV_5736" ">ASV_5737" ">ASV_5738" ">ASV_5739" ">ASV_5740"
## [5741] ">ASV_5741" ">ASV_5742" ">ASV_5743" ">ASV_5744" ">ASV_5745" ">ASV_5746" ">ASV_5747" ">ASV_5748" ">ASV_5749" ">ASV_5750" ">ASV_5751" ">ASV_5752" ">ASV_5753" ">ASV_5754" ">ASV_5755" ">ASV_5756" ">ASV_5757" ">ASV_5758" ">ASV_5759" ">ASV_5760"
## [5761] ">ASV_5761" ">ASV_5762" ">ASV_5763" ">ASV_5764" ">ASV_5765" ">ASV_5766" ">ASV_5767" ">ASV_5768" ">ASV_5769" ">ASV_5770" ">ASV_5771" ">ASV_5772" ">ASV_5773" ">ASV_5774" ">ASV_5775" ">ASV_5776" ">ASV_5777" ">ASV_5778" ">ASV_5779" ">ASV_5780"
## [5781] ">ASV_5781" ">ASV_5782" ">ASV_5783" ">ASV_5784" ">ASV_5785" ">ASV_5786" ">ASV_5787" ">ASV_5788" ">ASV_5789" ">ASV_5790" ">ASV_5791" ">ASV_5792" ">ASV_5793" ">ASV_5794" ">ASV_5795" ">ASV_5796" ">ASV_5797" ">ASV_5798" ">ASV_5799" ">ASV_5800"
## [5801] ">ASV_5801" ">ASV_5802" ">ASV_5803" ">ASV_5804" ">ASV_5805" ">ASV_5806" ">ASV_5807" ">ASV_5808" ">ASV_5809" ">ASV_5810" ">ASV_5811" ">ASV_5812" ">ASV_5813" ">ASV_5814" ">ASV_5815" ">ASV_5816" ">ASV_5817" ">ASV_5818" ">ASV_5819" ">ASV_5820"
## [5821] ">ASV_5821" ">ASV_5822" ">ASV_5823" ">ASV_5824" ">ASV_5825" ">ASV_5826" ">ASV_5827" ">ASV_5828" ">ASV_5829" ">ASV_5830" ">ASV_5831" ">ASV_5832" ">ASV_5833" ">ASV_5834" ">ASV_5835" ">ASV_5836" ">ASV_5837" ">ASV_5838" ">ASV_5839" ">ASV_5840"
## [5841] ">ASV_5841" ">ASV_5842" ">ASV_5843" ">ASV_5844" ">ASV_5845" ">ASV_5846" ">ASV_5847" ">ASV_5848" ">ASV_5849" ">ASV_5850" ">ASV_5851" ">ASV_5852" ">ASV_5853" ">ASV_5854" ">ASV_5855" ">ASV_5856" ">ASV_5857" ">ASV_5858" ">ASV_5859" ">ASV_5860"
## [5861] ">ASV_5861" ">ASV_5862" ">ASV_5863" ">ASV_5864" ">ASV_5865" ">ASV_5866" ">ASV_5867" ">ASV_5868" ">ASV_5869" ">ASV_5870" ">ASV_5871" ">ASV_5872" ">ASV_5873" ">ASV_5874" ">ASV_5875" ">ASV_5876" ">ASV_5877" ">ASV_5878" ">ASV_5879" ">ASV_5880"
## [5881] ">ASV_5881" ">ASV_5882" ">ASV_5883" ">ASV_5884" ">ASV_5885" ">ASV_5886" ">ASV_5887" ">ASV_5888" ">ASV_5889" ">ASV_5890" ">ASV_5891" ">ASV_5892" ">ASV_5893" ">ASV_5894" ">ASV_5895" ">ASV_5896" ">ASV_5897" ">ASV_5898" ">ASV_5899" ">ASV_5900"
## [5901] ">ASV_5901" ">ASV_5902" ">ASV_5903" ">ASV_5904" ">ASV_5905" ">ASV_5906" ">ASV_5907" ">ASV_5908" ">ASV_5909" ">ASV_5910" ">ASV_5911" ">ASV_5912" ">ASV_5913" ">ASV_5914" ">ASV_5915" ">ASV_5916" ">ASV_5917" ">ASV_5918" ">ASV_5919" ">ASV_5920"
## [5921] ">ASV_5921" ">ASV_5922" ">ASV_5923" ">ASV_5924" ">ASV_5925" ">ASV_5926" ">ASV_5927" ">ASV_5928" ">ASV_5929" ">ASV_5930" ">ASV_5931" ">ASV_5932" ">ASV_5933" ">ASV_5934" ">ASV_5935" ">ASV_5936" ">ASV_5937" ">ASV_5938" ">ASV_5939" ">ASV_5940"
## [5941] ">ASV_5941" ">ASV_5942" ">ASV_5943" ">ASV_5944" ">ASV_5945" ">ASV_5946" ">ASV_5947" ">ASV_5948" ">ASV_5949" ">ASV_5950" ">ASV_5951" ">ASV_5952" ">ASV_5953" ">ASV_5954" ">ASV_5955" ">ASV_5956" ">ASV_5957" ">ASV_5958" ">ASV_5959" ">ASV_5960"
## [5961] ">ASV_5961" ">ASV_5962" ">ASV_5963" ">ASV_5964" ">ASV_5965" ">ASV_5966" ">ASV_5967" ">ASV_5968" ">ASV_5969" ">ASV_5970" ">ASV_5971" ">ASV_5972" ">ASV_5973" ">ASV_5974" ">ASV_5975" ">ASV_5976" ">ASV_5977" ">ASV_5978" ">ASV_5979" ">ASV_5980"
## [5981] ">ASV_5981" ">ASV_5982" ">ASV_5983" ">ASV_5984" ">ASV_5985" ">ASV_5986" ">ASV_5987" ">ASV_5988" ">ASV_5989" ">ASV_5990" ">ASV_5991" ">ASV_5992" ">ASV_5993" ">ASV_5994" ">ASV_5995" ">ASV_5996" ">ASV_5997" ">ASV_5998" ">ASV_5999" ">ASV_6000"
## [6001] ">ASV_6001" ">ASV_6002" ">ASV_6003" ">ASV_6004" ">ASV_6005" ">ASV_6006" ">ASV_6007" ">ASV_6008" ">ASV_6009" ">ASV_6010" ">ASV_6011" ">ASV_6012" ">ASV_6013" ">ASV_6014" ">ASV_6015" ">ASV_6016" ">ASV_6017" ">ASV_6018" ">ASV_6019" ">ASV_6020"
## [6021] ">ASV_6021" ">ASV_6022" ">ASV_6023" ">ASV_6024" ">ASV_6025" ">ASV_6026" ">ASV_6027" ">ASV_6028" ">ASV_6029" ">ASV_6030" ">ASV_6031" ">ASV_6032" ">ASV_6033" ">ASV_6034" ">ASV_6035" ">ASV_6036" ">ASV_6037" ">ASV_6038" ">ASV_6039" ">ASV_6040"
## [6041] ">ASV_6041" ">ASV_6042" ">ASV_6043" ">ASV_6044" ">ASV_6045" ">ASV_6046" ">ASV_6047" ">ASV_6048" ">ASV_6049" ">ASV_6050" ">ASV_6051" ">ASV_6052" ">ASV_6053" ">ASV_6054" ">ASV_6055" ">ASV_6056" ">ASV_6057" ">ASV_6058" ">ASV_6059" ">ASV_6060"
## [6061] ">ASV_6061" ">ASV_6062" ">ASV_6063" ">ASV_6064" ">ASV_6065" ">ASV_6066" ">ASV_6067" ">ASV_6068" ">ASV_6069" ">ASV_6070" ">ASV_6071" ">ASV_6072" ">ASV_6073" ">ASV_6074" ">ASV_6075" ">ASV_6076" ">ASV_6077" ">ASV_6078" ">ASV_6079" ">ASV_6080"
## [6081] ">ASV_6081" ">ASV_6082" ">ASV_6083" ">ASV_6084" ">ASV_6085" ">ASV_6086" ">ASV_6087" ">ASV_6088" ">ASV_6089" ">ASV_6090" ">ASV_6091" ">ASV_6092" ">ASV_6093" ">ASV_6094" ">ASV_6095" ">ASV_6096" ">ASV_6097" ">ASV_6098" ">ASV_6099" ">ASV_6100"
## [6101] ">ASV_6101" ">ASV_6102" ">ASV_6103" ">ASV_6104" ">ASV_6105" ">ASV_6106" ">ASV_6107" ">ASV_6108" ">ASV_6109" ">ASV_6110" ">ASV_6111" ">ASV_6112" ">ASV_6113" ">ASV_6114" ">ASV_6115" ">ASV_6116" ">ASV_6117" ">ASV_6118" ">ASV_6119" ">ASV_6120"
## [6121] ">ASV_6121" ">ASV_6122" ">ASV_6123" ">ASV_6124" ">ASV_6125" ">ASV_6126" ">ASV_6127" ">ASV_6128" ">ASV_6129" ">ASV_6130" ">ASV_6131" ">ASV_6132" ">ASV_6133" ">ASV_6134" ">ASV_6135" ">ASV_6136" ">ASV_6137" ">ASV_6138" ">ASV_6139" ">ASV_6140"
## [6141] ">ASV_6141" ">ASV_6142" ">ASV_6143" ">ASV_6144" ">ASV_6145" ">ASV_6146" ">ASV_6147" ">ASV_6148" ">ASV_6149" ">ASV_6150" ">ASV_6151" ">ASV_6152" ">ASV_6153" ">ASV_6154" ">ASV_6155" ">ASV_6156" ">ASV_6157" ">ASV_6158" ">ASV_6159" ">ASV_6160"
## [6161] ">ASV_6161" ">ASV_6162" ">ASV_6163" ">ASV_6164" ">ASV_6165" ">ASV_6166" ">ASV_6167" ">ASV_6168" ">ASV_6169" ">ASV_6170" ">ASV_6171" ">ASV_6172" ">ASV_6173" ">ASV_6174" ">ASV_6175" ">ASV_6176" ">ASV_6177" ">ASV_6178" ">ASV_6179" ">ASV_6180"
## [6181] ">ASV_6181" ">ASV_6182" ">ASV_6183" ">ASV_6184" ">ASV_6185" ">ASV_6186" ">ASV_6187" ">ASV_6188" ">ASV_6189" ">ASV_6190" ">ASV_6191" ">ASV_6192" ">ASV_6193" ">ASV_6194" ">ASV_6195" ">ASV_6196" ">ASV_6197" ">ASV_6198" ">ASV_6199" ">ASV_6200"
## [6201] ">ASV_6201" ">ASV_6202" ">ASV_6203" ">ASV_6204" ">ASV_6205" ">ASV_6206" ">ASV_6207" ">ASV_6208" ">ASV_6209" ">ASV_6210" ">ASV_6211" ">ASV_6212" ">ASV_6213" ">ASV_6214" ">ASV_6215" ">ASV_6216" ">ASV_6217" ">ASV_6218" ">ASV_6219" ">ASV_6220"
## [6221] ">ASV_6221" ">ASV_6222" ">ASV_6223" ">ASV_6224" ">ASV_6225" ">ASV_6226" ">ASV_6227" ">ASV_6228" ">ASV_6229" ">ASV_6230" ">ASV_6231" ">ASV_6232" ">ASV_6233" ">ASV_6234" ">ASV_6235" ">ASV_6236" ">ASV_6237" ">ASV_6238" ">ASV_6239" ">ASV_6240"
## [6241] ">ASV_6241" ">ASV_6242" ">ASV_6243" ">ASV_6244" ">ASV_6245" ">ASV_6246" ">ASV_6247" ">ASV_6248" ">ASV_6249" ">ASV_6250" ">ASV_6251" ">ASV_6252" ">ASV_6253" ">ASV_6254" ">ASV_6255" ">ASV_6256" ">ASV_6257" ">ASV_6258" ">ASV_6259" ">ASV_6260"
## [6261] ">ASV_6261" ">ASV_6262" ">ASV_6263" ">ASV_6264" ">ASV_6265" ">ASV_6266" ">ASV_6267" ">ASV_6268" ">ASV_6269" ">ASV_6270" ">ASV_6271" ">ASV_6272" ">ASV_6273" ">ASV_6274" ">ASV_6275" ">ASV_6276" ">ASV_6277" ">ASV_6278" ">ASV_6279" ">ASV_6280"
## [6281] ">ASV_6281" ">ASV_6282" ">ASV_6283" ">ASV_6284" ">ASV_6285" ">ASV_6286" ">ASV_6287" ">ASV_6288" ">ASV_6289" ">ASV_6290" ">ASV_6291" ">ASV_6292" ">ASV_6293" ">ASV_6294" ">ASV_6295" ">ASV_6296" ">ASV_6297" ">ASV_6298" ">ASV_6299" ">ASV_6300"
## [6301] ">ASV_6301" ">ASV_6302" ">ASV_6303" ">ASV_6304" ">ASV_6305" ">ASV_6306" ">ASV_6307" ">ASV_6308" ">ASV_6309" ">ASV_6310" ">ASV_6311" ">ASV_6312" ">ASV_6313" ">ASV_6314" ">ASV_6315" ">ASV_6316" ">ASV_6317" ">ASV_6318" ">ASV_6319" ">ASV_6320"
## [6321] ">ASV_6321" ">ASV_6322" ">ASV_6323" ">ASV_6324" ">ASV_6325" ">ASV_6326" ">ASV_6327" ">ASV_6328" ">ASV_6329" ">ASV_6330" ">ASV_6331" ">ASV_6332" ">ASV_6333" ">ASV_6334" ">ASV_6335" ">ASV_6336" ">ASV_6337" ">ASV_6338" ">ASV_6339" ">ASV_6340"
## [6341] ">ASV_6341" ">ASV_6342" ">ASV_6343" ">ASV_6344" ">ASV_6345" ">ASV_6346" ">ASV_6347" ">ASV_6348" ">ASV_6349" ">ASV_6350" ">ASV_6351" ">ASV_6352" ">ASV_6353" ">ASV_6354" ">ASV_6355" ">ASV_6356" ">ASV_6357" ">ASV_6358" ">ASV_6359" ">ASV_6360"
## [6361] ">ASV_6361" ">ASV_6362" ">ASV_6363" ">ASV_6364" ">ASV_6365" ">ASV_6366" ">ASV_6367" ">ASV_6368" ">ASV_6369" ">ASV_6370" ">ASV_6371" ">ASV_6372" ">ASV_6373" ">ASV_6374" ">ASV_6375" ">ASV_6376" ">ASV_6377" ">ASV_6378" ">ASV_6379" ">ASV_6380"
## [6381] ">ASV_6381" ">ASV_6382" ">ASV_6383" ">ASV_6384" ">ASV_6385" ">ASV_6386" ">ASV_6387" ">ASV_6388" ">ASV_6389" ">ASV_6390" ">ASV_6391" ">ASV_6392" ">ASV_6393" ">ASV_6394" ">ASV_6395" ">ASV_6396" ">ASV_6397" ">ASV_6398" ">ASV_6399" ">ASV_6400"
## [6401] ">ASV_6401" ">ASV_6402" ">ASV_6403" ">ASV_6404" ">ASV_6405" ">ASV_6406" ">ASV_6407" ">ASV_6408" ">ASV_6409" ">ASV_6410" ">ASV_6411" ">ASV_6412" ">ASV_6413" ">ASV_6414" ">ASV_6415" ">ASV_6416" ">ASV_6417" ">ASV_6418" ">ASV_6419" ">ASV_6420"
## [6421] ">ASV_6421" ">ASV_6422" ">ASV_6423" ">ASV_6424" ">ASV_6425" ">ASV_6426" ">ASV_6427" ">ASV_6428" ">ASV_6429" ">ASV_6430" ">ASV_6431" ">ASV_6432" ">ASV_6433" ">ASV_6434" ">ASV_6435" ">ASV_6436" ">ASV_6437" ">ASV_6438" ">ASV_6439" ">ASV_6440"
## [6441] ">ASV_6441" ">ASV_6442" ">ASV_6443" ">ASV_6444" ">ASV_6445" ">ASV_6446" ">ASV_6447" ">ASV_6448" ">ASV_6449" ">ASV_6450" ">ASV_6451" ">ASV_6452" ">ASV_6453" ">ASV_6454" ">ASV_6455" ">ASV_6456" ">ASV_6457" ">ASV_6458" ">ASV_6459" ">ASV_6460"
## [6461] ">ASV_6461" ">ASV_6462" ">ASV_6463" ">ASV_6464" ">ASV_6465" ">ASV_6466" ">ASV_6467" ">ASV_6468" ">ASV_6469" ">ASV_6470" ">ASV_6471" ">ASV_6472" ">ASV_6473" ">ASV_6474" ">ASV_6475" ">ASV_6476" ">ASV_6477" ">ASV_6478" ">ASV_6479" ">ASV_6480"
## [6481] ">ASV_6481" ">ASV_6482" ">ASV_6483" ">ASV_6484" ">ASV_6485" ">ASV_6486" ">ASV_6487" ">ASV_6488" ">ASV_6489" ">ASV_6490" ">ASV_6491" ">ASV_6492" ">ASV_6493" ">ASV_6494" ">ASV_6495" ">ASV_6496" ">ASV_6497" ">ASV_6498" ">ASV_6499" ">ASV_6500"
## [6501] ">ASV_6501" ">ASV_6502" ">ASV_6503" ">ASV_6504" ">ASV_6505" ">ASV_6506" ">ASV_6507" ">ASV_6508" ">ASV_6509" ">ASV_6510" ">ASV_6511" ">ASV_6512" ">ASV_6513" ">ASV_6514" ">ASV_6515" ">ASV_6516" ">ASV_6517" ">ASV_6518" ">ASV_6519" ">ASV_6520"
## [6521] ">ASV_6521" ">ASV_6522" ">ASV_6523" ">ASV_6524" ">ASV_6525" ">ASV_6526" ">ASV_6527" ">ASV_6528" ">ASV_6529" ">ASV_6530" ">ASV_6531" ">ASV_6532" ">ASV_6533" ">ASV_6534" ">ASV_6535" ">ASV_6536" ">ASV_6537" ">ASV_6538" ">ASV_6539" ">ASV_6540"
## [6541] ">ASV_6541" ">ASV_6542" ">ASV_6543" ">ASV_6544" ">ASV_6545" ">ASV_6546" ">ASV_6547" ">ASV_6548" ">ASV_6549" ">ASV_6550" ">ASV_6551" ">ASV_6552" ">ASV_6553" ">ASV_6554" ">ASV_6555" ">ASV_6556" ">ASV_6557" ">ASV_6558" ">ASV_6559" ">ASV_6560"
## [6561] ">ASV_6561" ">ASV_6562" ">ASV_6563" ">ASV_6564" ">ASV_6565" ">ASV_6566" ">ASV_6567" ">ASV_6568" ">ASV_6569" ">ASV_6570" ">ASV_6571" ">ASV_6572" ">ASV_6573" ">ASV_6574" ">ASV_6575" ">ASV_6576" ">ASV_6577" ">ASV_6578" ">ASV_6579" ">ASV_6580"
## [6581] ">ASV_6581" ">ASV_6582" ">ASV_6583" ">ASV_6584" ">ASV_6585" ">ASV_6586" ">ASV_6587" ">ASV_6588" ">ASV_6589" ">ASV_6590" ">ASV_6591" ">ASV_6592" ">ASV_6593" ">ASV_6594" ">ASV_6595" ">ASV_6596" ">ASV_6597" ">ASV_6598" ">ASV_6599" ">ASV_6600"
## [6601] ">ASV_6601" ">ASV_6602" ">ASV_6603" ">ASV_6604" ">ASV_6605" ">ASV_6606" ">ASV_6607" ">ASV_6608" ">ASV_6609" ">ASV_6610" ">ASV_6611" ">ASV_6612" ">ASV_6613" ">ASV_6614" ">ASV_6615" ">ASV_6616" ">ASV_6617" ">ASV_6618" ">ASV_6619" ">ASV_6620"
## [6621] ">ASV_6621" ">ASV_6622" ">ASV_6623" ">ASV_6624" ">ASV_6625" ">ASV_6626" ">ASV_6627" ">ASV_6628" ">ASV_6629" ">ASV_6630" ">ASV_6631" ">ASV_6632" ">ASV_6633" ">ASV_6634" ">ASV_6635" ">ASV_6636" ">ASV_6637" ">ASV_6638" ">ASV_6639" ">ASV_6640"
## [6641] ">ASV_6641" ">ASV_6642" ">ASV_6643" ">ASV_6644" ">ASV_6645" ">ASV_6646" ">ASV_6647" ">ASV_6648" ">ASV_6649" ">ASV_6650" ">ASV_6651" ">ASV_6652" ">ASV_6653" ">ASV_6654" ">ASV_6655" ">ASV_6656" ">ASV_6657" ">ASV_6658" ">ASV_6659" ">ASV_6660"
## [6661] ">ASV_6661" ">ASV_6662" ">ASV_6663" ">ASV_6664" ">ASV_6665" ">ASV_6666" ">ASV_6667" ">ASV_6668" ">ASV_6669" ">ASV_6670" ">ASV_6671" ">ASV_6672" ">ASV_6673" ">ASV_6674" ">ASV_6675" ">ASV_6676" ">ASV_6677" ">ASV_6678" ">ASV_6679" ">ASV_6680"
## [6681] ">ASV_6681" ">ASV_6682" ">ASV_6683" ">ASV_6684" ">ASV_6685" ">ASV_6686" ">ASV_6687" ">ASV_6688" ">ASV_6689" ">ASV_6690" ">ASV_6691" ">ASV_6692" ">ASV_6693" ">ASV_6694" ">ASV_6695" ">ASV_6696" ">ASV_6697" ">ASV_6698" ">ASV_6699" ">ASV_6700"
## [6701] ">ASV_6701" ">ASV_6702" ">ASV_6703" ">ASV_6704" ">ASV_6705" ">ASV_6706" ">ASV_6707" ">ASV_6708" ">ASV_6709" ">ASV_6710" ">ASV_6711" ">ASV_6712" ">ASV_6713" ">ASV_6714" ">ASV_6715" ">ASV_6716" ">ASV_6717" ">ASV_6718" ">ASV_6719" ">ASV_6720"
## [6721] ">ASV_6721" ">ASV_6722" ">ASV_6723" ">ASV_6724" ">ASV_6725" ">ASV_6726" ">ASV_6727" ">ASV_6728" ">ASV_6729" ">ASV_6730" ">ASV_6731" ">ASV_6732" ">ASV_6733" ">ASV_6734" ">ASV_6735" ">ASV_6736" ">ASV_6737" ">ASV_6738" ">ASV_6739" ">ASV_6740"
## [6741] ">ASV_6741" ">ASV_6742" ">ASV_6743" ">ASV_6744" ">ASV_6745" ">ASV_6746" ">ASV_6747" ">ASV_6748" ">ASV_6749" ">ASV_6750" ">ASV_6751" ">ASV_6752" ">ASV_6753" ">ASV_6754" ">ASV_6755" ">ASV_6756" ">ASV_6757" ">ASV_6758" ">ASV_6759" ">ASV_6760"
## [6761] ">ASV_6761" ">ASV_6762" ">ASV_6763" ">ASV_6764" ">ASV_6765" ">ASV_6766" ">ASV_6767" ">ASV_6768" ">ASV_6769" ">ASV_6770" ">ASV_6771" ">ASV_6772" ">ASV_6773" ">ASV_6774" ">ASV_6775" ">ASV_6776" ">ASV_6777" ">ASV_6778" ">ASV_6779" ">ASV_6780"
## [6781] ">ASV_6781" ">ASV_6782" ">ASV_6783" ">ASV_6784" ">ASV_6785" ">ASV_6786" ">ASV_6787" ">ASV_6788" ">ASV_6789" ">ASV_6790" ">ASV_6791" ">ASV_6792" ">ASV_6793" ">ASV_6794" ">ASV_6795" ">ASV_6796" ">ASV_6797" ">ASV_6798" ">ASV_6799" ">ASV_6800"
## [6801] ">ASV_6801" ">ASV_6802" ">ASV_6803" ">ASV_6804" ">ASV_6805" ">ASV_6806" ">ASV_6807" ">ASV_6808" ">ASV_6809" ">ASV_6810" ">ASV_6811" ">ASV_6812" ">ASV_6813" ">ASV_6814" ">ASV_6815" ">ASV_6816" ">ASV_6817" ">ASV_6818" ">ASV_6819" ">ASV_6820"
## [6821] ">ASV_6821" ">ASV_6822" ">ASV_6823" ">ASV_6824" ">ASV_6825" ">ASV_6826" ">ASV_6827" ">ASV_6828" ">ASV_6829" ">ASV_6830" ">ASV_6831" ">ASV_6832" ">ASV_6833" ">ASV_6834" ">ASV_6835" ">ASV_6836" ">ASV_6837" ">ASV_6838" ">ASV_6839" ">ASV_6840"
## [6841] ">ASV_6841" ">ASV_6842" ">ASV_6843" ">ASV_6844" ">ASV_6845" ">ASV_6846" ">ASV_6847" ">ASV_6848" ">ASV_6849" ">ASV_6850" ">ASV_6851" ">ASV_6852" ">ASV_6853" ">ASV_6854" ">ASV_6855" ">ASV_6856" ">ASV_6857" ">ASV_6858" ">ASV_6859" ">ASV_6860"
## [6861] ">ASV_6861" ">ASV_6862" ">ASV_6863" ">ASV_6864" ">ASV_6865" ">ASV_6866" ">ASV_6867" ">ASV_6868" ">ASV_6869" ">ASV_6870" ">ASV_6871" ">ASV_6872" ">ASV_6873" ">ASV_6874" ">ASV_6875" ">ASV_6876" ">ASV_6877" ">ASV_6878" ">ASV_6879" ">ASV_6880"
## [6881] ">ASV_6881" ">ASV_6882" ">ASV_6883" ">ASV_6884" ">ASV_6885" ">ASV_6886" ">ASV_6887" ">ASV_6888" ">ASV_6889" ">ASV_6890" ">ASV_6891" ">ASV_6892" ">ASV_6893" ">ASV_6894" ">ASV_6895" ">ASV_6896" ">ASV_6897" ">ASV_6898" ">ASV_6899" ">ASV_6900"
## [6901] ">ASV_6901" ">ASV_6902" ">ASV_6903" ">ASV_6904" ">ASV_6905" ">ASV_6906" ">ASV_6907" ">ASV_6908" ">ASV_6909" ">ASV_6910" ">ASV_6911" ">ASV_6912" ">ASV_6913" ">ASV_6914" ">ASV_6915" ">ASV_6916" ">ASV_6917" ">ASV_6918" ">ASV_6919" ">ASV_6920"
## [6921] ">ASV_6921" ">ASV_6922" ">ASV_6923" ">ASV_6924" ">ASV_6925" ">ASV_6926" ">ASV_6927" ">ASV_6928" ">ASV_6929" ">ASV_6930" ">ASV_6931" ">ASV_6932" ">ASV_6933" ">ASV_6934" ">ASV_6935" ">ASV_6936" ">ASV_6937" ">ASV_6938" ">ASV_6939" ">ASV_6940"
## [6941] ">ASV_6941" ">ASV_6942" ">ASV_6943" ">ASV_6944" ">ASV_6945" ">ASV_6946" ">ASV_6947" ">ASV_6948" ">ASV_6949" ">ASV_6950" ">ASV_6951" ">ASV_6952" ">ASV_6953" ">ASV_6954" ">ASV_6955" ">ASV_6956" ">ASV_6957" ">ASV_6958" ">ASV_6959" ">ASV_6960"
## [6961] ">ASV_6961" ">ASV_6962" ">ASV_6963" ">ASV_6964" ">ASV_6965" ">ASV_6966" ">ASV_6967" ">ASV_6968" ">ASV_6969" ">ASV_6970" ">ASV_6971" ">ASV_6972" ">ASV_6973" ">ASV_6974" ">ASV_6975" ">ASV_6976" ">ASV_6977" ">ASV_6978" ">ASV_6979" ">ASV_6980"
## [6981] ">ASV_6981" ">ASV_6982" ">ASV_6983" ">ASV_6984" ">ASV_6985" ">ASV_6986" ">ASV_6987" ">ASV_6988" ">ASV_6989" ">ASV_6990" ">ASV_6991" ">ASV_6992" ">ASV_6993" ">ASV_6994" ">ASV_6995" ">ASV_6996" ">ASV_6997" ">ASV_6998" ">ASV_6999" ">ASV_7000"
## [7001] ">ASV_7001" ">ASV_7002" ">ASV_7003" ">ASV_7004" ">ASV_7005" ">ASV_7006" ">ASV_7007" ">ASV_7008" ">ASV_7009" ">ASV_7010" ">ASV_7011" ">ASV_7012" ">ASV_7013" ">ASV_7014" ">ASV_7015" ">ASV_7016" ">ASV_7017" ">ASV_7018" ">ASV_7019" ">ASV_7020"
## [7021] ">ASV_7021" ">ASV_7022" ">ASV_7023" ">ASV_7024" ">ASV_7025" ">ASV_7026" ">ASV_7027" ">ASV_7028" ">ASV_7029" ">ASV_7030" ">ASV_7031" ">ASV_7032" ">ASV_7033" ">ASV_7034" ">ASV_7035" ">ASV_7036" ">ASV_7037" ">ASV_7038" ">ASV_7039" ">ASV_7040"
## [7041] ">ASV_7041" ">ASV_7042" ">ASV_7043" ">ASV_7044" ">ASV_7045" ">ASV_7046" ">ASV_7047" ">ASV_7048" ">ASV_7049" ">ASV_7050" ">ASV_7051" ">ASV_7052" ">ASV_7053" ">ASV_7054" ">ASV_7055" ">ASV_7056" ">ASV_7057" ">ASV_7058" ">ASV_7059" ">ASV_7060"
## [7061] ">ASV_7061" ">ASV_7062" ">ASV_7063" ">ASV_7064" ">ASV_7065" ">ASV_7066" ">ASV_7067" ">ASV_7068" ">ASV_7069" ">ASV_7070" ">ASV_7071" ">ASV_7072" ">ASV_7073" ">ASV_7074" ">ASV_7075" ">ASV_7076" ">ASV_7077" ">ASV_7078" ">ASV_7079" ">ASV_7080"
## [7081] ">ASV_7081" ">ASV_7082" ">ASV_7083" ">ASV_7084" ">ASV_7085" ">ASV_7086" ">ASV_7087" ">ASV_7088" ">ASV_7089" ">ASV_7090" ">ASV_7091" ">ASV_7092" ">ASV_7093" ">ASV_7094" ">ASV_7095" ">ASV_7096" ">ASV_7097" ">ASV_7098" ">ASV_7099" ">ASV_7100"
## [7101] ">ASV_7101" ">ASV_7102" ">ASV_7103" ">ASV_7104" ">ASV_7105" ">ASV_7106" ">ASV_7107" ">ASV_7108" ">ASV_7109" ">ASV_7110" ">ASV_7111" ">ASV_7112" ">ASV_7113" ">ASV_7114" ">ASV_7115" ">ASV_7116" ">ASV_7117" ">ASV_7118" ">ASV_7119" ">ASV_7120"
## [7121] ">ASV_7121" ">ASV_7122" ">ASV_7123" ">ASV_7124" ">ASV_7125" ">ASV_7126" ">ASV_7127" ">ASV_7128" ">ASV_7129" ">ASV_7130" ">ASV_7131" ">ASV_7132" ">ASV_7133" ">ASV_7134" ">ASV_7135" ">ASV_7136" ">ASV_7137" ">ASV_7138" ">ASV_7139" ">ASV_7140"
## [7141] ">ASV_7141" ">ASV_7142" ">ASV_7143" ">ASV_7144" ">ASV_7145" ">ASV_7146" ">ASV_7147" ">ASV_7148" ">ASV_7149" ">ASV_7150" ">ASV_7151" ">ASV_7152" ">ASV_7153" ">ASV_7154" ">ASV_7155" ">ASV_7156" ">ASV_7157" ">ASV_7158" ">ASV_7159" ">ASV_7160"
## [7161] ">ASV_7161" ">ASV_7162" ">ASV_7163" ">ASV_7164" ">ASV_7165" ">ASV_7166" ">ASV_7167" ">ASV_7168" ">ASV_7169" ">ASV_7170" ">ASV_7171" ">ASV_7172" ">ASV_7173" ">ASV_7174" ">ASV_7175" ">ASV_7176" ">ASV_7177" ">ASV_7178" ">ASV_7179" ">ASV_7180"
## [7181] ">ASV_7181" ">ASV_7182" ">ASV_7183" ">ASV_7184" ">ASV_7185" ">ASV_7186" ">ASV_7187" ">ASV_7188" ">ASV_7189" ">ASV_7190" ">ASV_7191" ">ASV_7192" ">ASV_7193" ">ASV_7194" ">ASV_7195" ">ASV_7196" ">ASV_7197" ">ASV_7198" ">ASV_7199" ">ASV_7200"
## [7201] ">ASV_7201" ">ASV_7202" ">ASV_7203" ">ASV_7204" ">ASV_7205" ">ASV_7206" ">ASV_7207" ">ASV_7208" ">ASV_7209" ">ASV_7210" ">ASV_7211" ">ASV_7212" ">ASV_7213" ">ASV_7214" ">ASV_7215" ">ASV_7216" ">ASV_7217" ">ASV_7218" ">ASV_7219" ">ASV_7220"
## [7221] ">ASV_7221" ">ASV_7222" ">ASV_7223" ">ASV_7224" ">ASV_7225" ">ASV_7226" ">ASV_7227" ">ASV_7228" ">ASV_7229" ">ASV_7230" ">ASV_7231" ">ASV_7232" ">ASV_7233" ">ASV_7234" ">ASV_7235" ">ASV_7236" ">ASV_7237" ">ASV_7238" ">ASV_7239" ">ASV_7240"
## [7241] ">ASV_7241" ">ASV_7242" ">ASV_7243" ">ASV_7244" ">ASV_7245" ">ASV_7246" ">ASV_7247" ">ASV_7248" ">ASV_7249" ">ASV_7250" ">ASV_7251" ">ASV_7252" ">ASV_7253" ">ASV_7254" ">ASV_7255" ">ASV_7256" ">ASV_7257" ">ASV_7258" ">ASV_7259" ">ASV_7260"
## [7261] ">ASV_7261" ">ASV_7262" ">ASV_7263" ">ASV_7264" ">ASV_7265" ">ASV_7266" ">ASV_7267" ">ASV_7268" ">ASV_7269" ">ASV_7270" ">ASV_7271" ">ASV_7272" ">ASV_7273" ">ASV_7274" ">ASV_7275" ">ASV_7276" ">ASV_7277" ">ASV_7278" ">ASV_7279" ">ASV_7280"
## [7281] ">ASV_7281" ">ASV_7282" ">ASV_7283" ">ASV_7284" ">ASV_7285" ">ASV_7286" ">ASV_7287" ">ASV_7288" ">ASV_7289" ">ASV_7290" ">ASV_7291" ">ASV_7292" ">ASV_7293" ">ASV_7294" ">ASV_7295" ">ASV_7296" ">ASV_7297" ">ASV_7298" ">ASV_7299" ">ASV_7300"
## [7301] ">ASV_7301" ">ASV_7302" ">ASV_7303" ">ASV_7304" ">ASV_7305" ">ASV_7306" ">ASV_7307" ">ASV_7308" ">ASV_7309" ">ASV_7310" ">ASV_7311" ">ASV_7312" ">ASV_7313" ">ASV_7314" ">ASV_7315" ">ASV_7316" ">ASV_7317" ">ASV_7318" ">ASV_7319" ">ASV_7320"
## [7321] ">ASV_7321" ">ASV_7322" ">ASV_7323" ">ASV_7324" ">ASV_7325" ">ASV_7326" ">ASV_7327" ">ASV_7328" ">ASV_7329" ">ASV_7330" ">ASV_7331" ">ASV_7332" ">ASV_7333" ">ASV_7334" ">ASV_7335" ">ASV_7336" ">ASV_7337" ">ASV_7338" ">ASV_7339" ">ASV_7340"
## [7341] ">ASV_7341" ">ASV_7342" ">ASV_7343" ">ASV_7344" ">ASV_7345" ">ASV_7346" ">ASV_7347" ">ASV_7348" ">ASV_7349" ">ASV_7350" ">ASV_7351" ">ASV_7352" ">ASV_7353" ">ASV_7354" ">ASV_7355" ">ASV_7356" ">ASV_7357" ">ASV_7358" ">ASV_7359" ">ASV_7360"
## [7361] ">ASV_7361" ">ASV_7362" ">ASV_7363" ">ASV_7364" ">ASV_7365" ">ASV_7366" ">ASV_7367" ">ASV_7368" ">ASV_7369" ">ASV_7370" ">ASV_7371" ">ASV_7372" ">ASV_7373" ">ASV_7374" ">ASV_7375" ">ASV_7376" ">ASV_7377" ">ASV_7378" ">ASV_7379" ">ASV_7380"
## [7381] ">ASV_7381" ">ASV_7382" ">ASV_7383" ">ASV_7384" ">ASV_7385" ">ASV_7386" ">ASV_7387" ">ASV_7388" ">ASV_7389" ">ASV_7390" ">ASV_7391" ">ASV_7392" ">ASV_7393" ">ASV_7394" ">ASV_7395" ">ASV_7396" ">ASV_7397" ">ASV_7398" ">ASV_7399" ">ASV_7400"
## [7401] ">ASV_7401" ">ASV_7402" ">ASV_7403" ">ASV_7404" ">ASV_7405" ">ASV_7406" ">ASV_7407" ">ASV_7408" ">ASV_7409" ">ASV_7410" ">ASV_7411" ">ASV_7412" ">ASV_7413" ">ASV_7414" ">ASV_7415" ">ASV_7416" ">ASV_7417" ">ASV_7418" ">ASV_7419" ">ASV_7420"
## [7421] ">ASV_7421" ">ASV_7422" ">ASV_7423" ">ASV_7424" ">ASV_7425" ">ASV_7426" ">ASV_7427" ">ASV_7428" ">ASV_7429" ">ASV_7430" ">ASV_7431" ">ASV_7432" ">ASV_7433" ">ASV_7434" ">ASV_7435" ">ASV_7436" ">ASV_7437" ">ASV_7438" ">ASV_7439" ">ASV_7440"
## [7441] ">ASV_7441" ">ASV_7442" ">ASV_7443" ">ASV_7444" ">ASV_7445" ">ASV_7446" ">ASV_7447" ">ASV_7448" ">ASV_7449" ">ASV_7450" ">ASV_7451" ">ASV_7452" ">ASV_7453" ">ASV_7454" ">ASV_7455" ">ASV_7456" ">ASV_7457" ">ASV_7458" ">ASV_7459" ">ASV_7460"
## [7461] ">ASV_7461" ">ASV_7462" ">ASV_7463" ">ASV_7464" ">ASV_7465" ">ASV_7466" ">ASV_7467" ">ASV_7468" ">ASV_7469" ">ASV_7470" ">ASV_7471" ">ASV_7472" ">ASV_7473" ">ASV_7474" ">ASV_7475" ">ASV_7476" ">ASV_7477" ">ASV_7478" ">ASV_7479" ">ASV_7480"
## [7481] ">ASV_7481" ">ASV_7482" ">ASV_7483" ">ASV_7484" ">ASV_7485" ">ASV_7486" ">ASV_7487" ">ASV_7488" ">ASV_7489" ">ASV_7490" ">ASV_7491" ">ASV_7492" ">ASV_7493" ">ASV_7494" ">ASV_7495" ">ASV_7496" ">ASV_7497" ">ASV_7498" ">ASV_7499" ">ASV_7500"
## [7501] ">ASV_7501" ">ASV_7502" ">ASV_7503" ">ASV_7504" ">ASV_7505" ">ASV_7506" ">ASV_7507" ">ASV_7508" ">ASV_7509" ">ASV_7510" ">ASV_7511" ">ASV_7512" ">ASV_7513" ">ASV_7514" ">ASV_7515" ">ASV_7516" ">ASV_7517" ">ASV_7518" ">ASV_7519" ">ASV_7520"
## [7521] ">ASV_7521" ">ASV_7522" ">ASV_7523" ">ASV_7524" ">ASV_7525" ">ASV_7526" ">ASV_7527" ">ASV_7528" ">ASV_7529" ">ASV_7530" ">ASV_7531" ">ASV_7532" ">ASV_7533" ">ASV_7534" ">ASV_7535" ">ASV_7536" ">ASV_7537" ">ASV_7538" ">ASV_7539" ">ASV_7540"
## [7541] ">ASV_7541" ">ASV_7542" ">ASV_7543" ">ASV_7544" ">ASV_7545" ">ASV_7546" ">ASV_7547" ">ASV_7548" ">ASV_7549" ">ASV_7550" ">ASV_7551" ">ASV_7552" ">ASV_7553" ">ASV_7554" ">ASV_7555" ">ASV_7556" ">ASV_7557" ">ASV_7558" ">ASV_7559" ">ASV_7560"
## [7561] ">ASV_7561" ">ASV_7562" ">ASV_7563" ">ASV_7564" ">ASV_7565" ">ASV_7566" ">ASV_7567" ">ASV_7568" ">ASV_7569" ">ASV_7570" ">ASV_7571" ">ASV_7572" ">ASV_7573" ">ASV_7574" ">ASV_7575" ">ASV_7576" ">ASV_7577" ">ASV_7578" ">ASV_7579" ">ASV_7580"
## [7581] ">ASV_7581" ">ASV_7582" ">ASV_7583" ">ASV_7584" ">ASV_7585" ">ASV_7586" ">ASV_7587" ">ASV_7588" ">ASV_7589" ">ASV_7590" ">ASV_7591" ">ASV_7592" ">ASV_7593" ">ASV_7594" ">ASV_7595" ">ASV_7596" ">ASV_7597" ">ASV_7598" ">ASV_7599" ">ASV_7600"
## [7601] ">ASV_7601" ">ASV_7602" ">ASV_7603" ">ASV_7604" ">ASV_7605" ">ASV_7606" ">ASV_7607" ">ASV_7608" ">ASV_7609" ">ASV_7610" ">ASV_7611" ">ASV_7612" ">ASV_7613" ">ASV_7614" ">ASV_7615" ">ASV_7616" ">ASV_7617" ">ASV_7618" ">ASV_7619" ">ASV_7620"
## [7621] ">ASV_7621" ">ASV_7622" ">ASV_7623" ">ASV_7624" ">ASV_7625" ">ASV_7626" ">ASV_7627" ">ASV_7628" ">ASV_7629" ">ASV_7630" ">ASV_7631" ">ASV_7632" ">ASV_7633" ">ASV_7634" ">ASV_7635" ">ASV_7636" ">ASV_7637" ">ASV_7638" ">ASV_7639" ">ASV_7640"
## [7641] ">ASV_7641" ">ASV_7642" ">ASV_7643" ">ASV_7644" ">ASV_7645" ">ASV_7646" ">ASV_7647" ">ASV_7648" ">ASV_7649" ">ASV_7650" ">ASV_7651" ">ASV_7652" ">ASV_7653" ">ASV_7654" ">ASV_7655" ">ASV_7656" ">ASV_7657" ">ASV_7658" ">ASV_7659" ">ASV_7660"
## [7661] ">ASV_7661" ">ASV_7662" ">ASV_7663" ">ASV_7664" ">ASV_7665" ">ASV_7666" ">ASV_7667" ">ASV_7668" ">ASV_7669" ">ASV_7670" ">ASV_7671" ">ASV_7672" ">ASV_7673" ">ASV_7674" ">ASV_7675" ">ASV_7676" ">ASV_7677" ">ASV_7678" ">ASV_7679" ">ASV_7680"
## [7681] ">ASV_7681" ">ASV_7682" ">ASV_7683" ">ASV_7684" ">ASV_7685" ">ASV_7686" ">ASV_7687" ">ASV_7688" ">ASV_7689" ">ASV_7690" ">ASV_7691" ">ASV_7692" ">ASV_7693" ">ASV_7694" ">ASV_7695" ">ASV_7696" ">ASV_7697" ">ASV_7698" ">ASV_7699" ">ASV_7700"
## [7701] ">ASV_7701" ">ASV_7702" ">ASV_7703" ">ASV_7704" ">ASV_7705" ">ASV_7706" ">ASV_7707" ">ASV_7708" ">ASV_7709" ">ASV_7710" ">ASV_7711" ">ASV_7712" ">ASV_7713" ">ASV_7714" ">ASV_7715" ">ASV_7716" ">ASV_7717" ">ASV_7718" ">ASV_7719" ">ASV_7720"
## [7721] ">ASV_7721" ">ASV_7722" ">ASV_7723" ">ASV_7724" ">ASV_7725" ">ASV_7726" ">ASV_7727" ">ASV_7728" ">ASV_7729" ">ASV_7730" ">ASV_7731" ">ASV_7732" ">ASV_7733" ">ASV_7734" ">ASV_7735" ">ASV_7736" ">ASV_7737" ">ASV_7738" ">ASV_7739" ">ASV_7740"
## [7741] ">ASV_7741" ">ASV_7742" ">ASV_7743" ">ASV_7744" ">ASV_7745" ">ASV_7746" ">ASV_7747" ">ASV_7748" ">ASV_7749" ">ASV_7750" ">ASV_7751" ">ASV_7752" ">ASV_7753" ">ASV_7754" ">ASV_7755" ">ASV_7756" ">ASV_7757" ">ASV_7758" ">ASV_7759" ">ASV_7760"
## [7761] ">ASV_7761" ">ASV_7762" ">ASV_7763" ">ASV_7764" ">ASV_7765" ">ASV_7766" ">ASV_7767" ">ASV_7768" ">ASV_7769" ">ASV_7770" ">ASV_7771" ">ASV_7772" ">ASV_7773" ">ASV_7774" ">ASV_7775" ">ASV_7776" ">ASV_7777" ">ASV_7778" ">ASV_7779" ">ASV_7780"
## [7781] ">ASV_7781" ">ASV_7782" ">ASV_7783" ">ASV_7784" ">ASV_7785" ">ASV_7786" ">ASV_7787" ">ASV_7788" ">ASV_7789" ">ASV_7790" ">ASV_7791" ">ASV_7792" ">ASV_7793" ">ASV_7794" ">ASV_7795" ">ASV_7796" ">ASV_7797" ">ASV_7798" ">ASV_7799" ">ASV_7800"
## [7801] ">ASV_7801" ">ASV_7802" ">ASV_7803" ">ASV_7804" ">ASV_7805" ">ASV_7806" ">ASV_7807" ">ASV_7808" ">ASV_7809" ">ASV_7810" ">ASV_7811" ">ASV_7812" ">ASV_7813" ">ASV_7814" ">ASV_7815" ">ASV_7816" ">ASV_7817" ">ASV_7818" ">ASV_7819" ">ASV_7820"
## [7821] ">ASV_7821" ">ASV_7822" ">ASV_7823" ">ASV_7824" ">ASV_7825" ">ASV_7826" ">ASV_7827" ">ASV_7828" ">ASV_7829" ">ASV_7830" ">ASV_7831" ">ASV_7832" ">ASV_7833" ">ASV_7834" ">ASV_7835" ">ASV_7836" ">ASV_7837" ">ASV_7838" ">ASV_7839" ">ASV_7840"
## [7841] ">ASV_7841" ">ASV_7842" ">ASV_7843" ">ASV_7844" ">ASV_7845" ">ASV_7846" ">ASV_7847" ">ASV_7848" ">ASV_7849" ">ASV_7850" ">ASV_7851" ">ASV_7852" ">ASV_7853" ">ASV_7854" ">ASV_7855" ">ASV_7856" ">ASV_7857" ">ASV_7858" ">ASV_7859" ">ASV_7860"
## [7861] ">ASV_7861" ">ASV_7862" ">ASV_7863" ">ASV_7864" ">ASV_7865" ">ASV_7866" ">ASV_7867" ">ASV_7868" ">ASV_7869" ">ASV_7870" ">ASV_7871" ">ASV_7872" ">ASV_7873" ">ASV_7874" ">ASV_7875" ">ASV_7876" ">ASV_7877" ">ASV_7878" ">ASV_7879" ">ASV_7880"
## [7881] ">ASV_7881" ">ASV_7882" ">ASV_7883" ">ASV_7884" ">ASV_7885" ">ASV_7886" ">ASV_7887" ">ASV_7888" ">ASV_7889" ">ASV_7890" ">ASV_7891" ">ASV_7892" ">ASV_7893" ">ASV_7894" ">ASV_7895" ">ASV_7896" ">ASV_7897" ">ASV_7898" ">ASV_7899" ">ASV_7900"
## [7901] ">ASV_7901" ">ASV_7902" ">ASV_7903" ">ASV_7904" ">ASV_7905" ">ASV_7906" ">ASV_7907" ">ASV_7908" ">ASV_7909" ">ASV_7910" ">ASV_7911" ">ASV_7912" ">ASV_7913" ">ASV_7914" ">ASV_7915" ">ASV_7916" ">ASV_7917" ">ASV_7918" ">ASV_7919" ">ASV_7920"
## [7921] ">ASV_7921" ">ASV_7922" ">ASV_7923" ">ASV_7924" ">ASV_7925" ">ASV_7926" ">ASV_7927" ">ASV_7928" ">ASV_7929" ">ASV_7930" ">ASV_7931" ">ASV_7932" ">ASV_7933" ">ASV_7934" ">ASV_7935" ">ASV_7936" ">ASV_7937" ">ASV_7938" ">ASV_7939" ">ASV_7940"
## [7941] ">ASV_7941" ">ASV_7942" ">ASV_7943" ">ASV_7944" ">ASV_7945" ">ASV_7946" ">ASV_7947" ">ASV_7948" ">ASV_7949" ">ASV_7950" ">ASV_7951" ">ASV_7952" ">ASV_7953" ">ASV_7954" ">ASV_7955" ">ASV_7956" ">ASV_7957" ">ASV_7958" ">ASV_7959" ">ASV_7960"
## [7961] ">ASV_7961" ">ASV_7962" ">ASV_7963" ">ASV_7964" ">ASV_7965" ">ASV_7966" ">ASV_7967" ">ASV_7968" ">ASV_7969" ">ASV_7970" ">ASV_7971" ">ASV_7972" ">ASV_7973" ">ASV_7974" ">ASV_7975" ">ASV_7976" ">ASV_7977" ">ASV_7978" ">ASV_7979" ">ASV_7980"
## [7981] ">ASV_7981" ">ASV_7982" ">ASV_7983" ">ASV_7984" ">ASV_7985" ">ASV_7986" ">ASV_7987" ">ASV_7988" ">ASV_7989" ">ASV_7990" ">ASV_7991" ">ASV_7992" ">ASV_7993" ">ASV_7994" ">ASV_7995" ">ASV_7996" ">ASV_7997" ">ASV_7998" ">ASV_7999" ">ASV_8000"
## [8001] ">ASV_8001" ">ASV_8002" ">ASV_8003" ">ASV_8004" ">ASV_8005" ">ASV_8006" ">ASV_8007" ">ASV_8008" ">ASV_8009" ">ASV_8010" ">ASV_8011" ">ASV_8012" ">ASV_8013" ">ASV_8014" ">ASV_8015" ">ASV_8016" ">ASV_8017" ">ASV_8018" ">ASV_8019" ">ASV_8020"
## [8021] ">ASV_8021" ">ASV_8022" ">ASV_8023" ">ASV_8024" ">ASV_8025" ">ASV_8026" ">ASV_8027" ">ASV_8028" ">ASV_8029" ">ASV_8030" ">ASV_8031" ">ASV_8032" ">ASV_8033" ">ASV_8034" ">ASV_8035" ">ASV_8036" ">ASV_8037" ">ASV_8038" ">ASV_8039" ">ASV_8040"
## [8041] ">ASV_8041" ">ASV_8042" ">ASV_8043" ">ASV_8044" ">ASV_8045" ">ASV_8046" ">ASV_8047" ">ASV_8048" ">ASV_8049" ">ASV_8050" ">ASV_8051" ">ASV_8052" ">ASV_8053" ">ASV_8054" ">ASV_8055" ">ASV_8056" ">ASV_8057" ">ASV_8058" ">ASV_8059" ">ASV_8060"
## [8061] ">ASV_8061" ">ASV_8062" ">ASV_8063" ">ASV_8064" ">ASV_8065" ">ASV_8066" ">ASV_8067" ">ASV_8068" ">ASV_8069" ">ASV_8070" ">ASV_8071" ">ASV_8072" ">ASV_8073" ">ASV_8074" ">ASV_8075" ">ASV_8076" ">ASV_8077" ">ASV_8078" ">ASV_8079" ">ASV_8080"
## [8081] ">ASV_8081" ">ASV_8082" ">ASV_8083" ">ASV_8084" ">ASV_8085" ">ASV_8086" ">ASV_8087" ">ASV_8088" ">ASV_8089" ">ASV_8090" ">ASV_8091" ">ASV_8092" ">ASV_8093" ">ASV_8094" ">ASV_8095" ">ASV_8096" ">ASV_8097" ">ASV_8098" ">ASV_8099" ">ASV_8100"
## [8101] ">ASV_8101" ">ASV_8102" ">ASV_8103" ">ASV_8104" ">ASV_8105" ">ASV_8106" ">ASV_8107" ">ASV_8108" ">ASV_8109" ">ASV_8110" ">ASV_8111" ">ASV_8112" ">ASV_8113" ">ASV_8114" ">ASV_8115" ">ASV_8116" ">ASV_8117" ">ASV_8118" ">ASV_8119" ">ASV_8120"
## [8121] ">ASV_8121" ">ASV_8122" ">ASV_8123" ">ASV_8124" ">ASV_8125" ">ASV_8126" ">ASV_8127" ">ASV_8128" ">ASV_8129" ">ASV_8130" ">ASV_8131" ">ASV_8132" ">ASV_8133" ">ASV_8134" ">ASV_8135" ">ASV_8136" ">ASV_8137" ">ASV_8138" ">ASV_8139" ">ASV_8140"
## [8141] ">ASV_8141" ">ASV_8142" ">ASV_8143" ">ASV_8144" ">ASV_8145" ">ASV_8146" ">ASV_8147" ">ASV_8148" ">ASV_8149" ">ASV_8150" ">ASV_8151" ">ASV_8152" ">ASV_8153" ">ASV_8154" ">ASV_8155" ">ASV_8156" ">ASV_8157" ">ASV_8158" ">ASV_8159" ">ASV_8160"
## [8161] ">ASV_8161" ">ASV_8162" ">ASV_8163" ">ASV_8164" ">ASV_8165" ">ASV_8166" ">ASV_8167" ">ASV_8168" ">ASV_8169" ">ASV_8170" ">ASV_8171" ">ASV_8172" ">ASV_8173" ">ASV_8174" ">ASV_8175" ">ASV_8176" ">ASV_8177" ">ASV_8178" ">ASV_8179" ">ASV_8180"
## [8181] ">ASV_8181" ">ASV_8182" ">ASV_8183" ">ASV_8184" ">ASV_8185" ">ASV_8186" ">ASV_8187" ">ASV_8188" ">ASV_8189" ">ASV_8190" ">ASV_8191" ">ASV_8192" ">ASV_8193" ">ASV_8194" ">ASV_8195" ">ASV_8196" ">ASV_8197" ">ASV_8198" ">ASV_8199" ">ASV_8200"
## [8201] ">ASV_8201" ">ASV_8202" ">ASV_8203" ">ASV_8204" ">ASV_8205" ">ASV_8206" ">ASV_8207" ">ASV_8208" ">ASV_8209" ">ASV_8210" ">ASV_8211" ">ASV_8212" ">ASV_8213" ">ASV_8214" ">ASV_8215" ">ASV_8216" ">ASV_8217" ">ASV_8218" ">ASV_8219" ">ASV_8220"
## [8221] ">ASV_8221" ">ASV_8222" ">ASV_8223" ">ASV_8224" ">ASV_8225" ">ASV_8226" ">ASV_8227" ">ASV_8228" ">ASV_8229" ">ASV_8230" ">ASV_8231" ">ASV_8232" ">ASV_8233" ">ASV_8234" ">ASV_8235" ">ASV_8236" ">ASV_8237" ">ASV_8238" ">ASV_8239" ">ASV_8240"
## [8241] ">ASV_8241" ">ASV_8242" ">ASV_8243" ">ASV_8244" ">ASV_8245" ">ASV_8246" ">ASV_8247" ">ASV_8248" ">ASV_8249" ">ASV_8250" ">ASV_8251" ">ASV_8252" ">ASV_8253" ">ASV_8254" ">ASV_8255" ">ASV_8256" ">ASV_8257" ">ASV_8258" ">ASV_8259" ">ASV_8260"
## [8261] ">ASV_8261" ">ASV_8262" ">ASV_8263" ">ASV_8264" ">ASV_8265" ">ASV_8266" ">ASV_8267" ">ASV_8268" ">ASV_8269" ">ASV_8270" ">ASV_8271" ">ASV_8272" ">ASV_8273" ">ASV_8274" ">ASV_8275" ">ASV_8276" ">ASV_8277" ">ASV_8278" ">ASV_8279" ">ASV_8280"
## [8281] ">ASV_8281" ">ASV_8282" ">ASV_8283" ">ASV_8284" ">ASV_8285" ">ASV_8286" ">ASV_8287" ">ASV_8288" ">ASV_8289" ">ASV_8290" ">ASV_8291" ">ASV_8292" ">ASV_8293" ">ASV_8294" ">ASV_8295" ">ASV_8296" ">ASV_8297" ">ASV_8298" ">ASV_8299" ">ASV_8300"
## [8301] ">ASV_8301" ">ASV_8302" ">ASV_8303" ">ASV_8304" ">ASV_8305" ">ASV_8306" ">ASV_8307" ">ASV_8308" ">ASV_8309" ">ASV_8310" ">ASV_8311" ">ASV_8312" ">ASV_8313" ">ASV_8314" ">ASV_8315" ">ASV_8316" ">ASV_8317" ">ASV_8318" ">ASV_8319" ">ASV_8320"
## [8321] ">ASV_8321" ">ASV_8322" ">ASV_8323" ">ASV_8324" ">ASV_8325" ">ASV_8326" ">ASV_8327" ">ASV_8328" ">ASV_8329" ">ASV_8330" ">ASV_8331" ">ASV_8332" ">ASV_8333" ">ASV_8334" ">ASV_8335" ">ASV_8336" ">ASV_8337" ">ASV_8338" ">ASV_8339" ">ASV_8340"
## [8341] ">ASV_8341" ">ASV_8342" ">ASV_8343" ">ASV_8344" ">ASV_8345" ">ASV_8346" ">ASV_8347" ">ASV_8348" ">ASV_8349" ">ASV_8350" ">ASV_8351" ">ASV_8352" ">ASV_8353" ">ASV_8354" ">ASV_8355" ">ASV_8356" ">ASV_8357" ">ASV_8358" ">ASV_8359" ">ASV_8360"
## [8361] ">ASV_8361" ">ASV_8362" ">ASV_8363" ">ASV_8364" ">ASV_8365" ">ASV_8366" ">ASV_8367" ">ASV_8368" ">ASV_8369" ">ASV_8370" ">ASV_8371" ">ASV_8372" ">ASV_8373" ">ASV_8374" ">ASV_8375" ">ASV_8376" ">ASV_8377" ">ASV_8378" ">ASV_8379" ">ASV_8380"
## [8381] ">ASV_8381" ">ASV_8382" ">ASV_8383" ">ASV_8384" ">ASV_8385" ">ASV_8386" ">ASV_8387" ">ASV_8388" ">ASV_8389" ">ASV_8390" ">ASV_8391" ">ASV_8392" ">ASV_8393" ">ASV_8394" ">ASV_8395" ">ASV_8396" ">ASV_8397" ">ASV_8398" ">ASV_8399" ">ASV_8400"
## [8401] ">ASV_8401" ">ASV_8402" ">ASV_8403" ">ASV_8404" ">ASV_8405" ">ASV_8406" ">ASV_8407" ">ASV_8408" ">ASV_8409" ">ASV_8410" ">ASV_8411" ">ASV_8412" ">ASV_8413" ">ASV_8414" ">ASV_8415" ">ASV_8416" ">ASV_8417" ">ASV_8418" ">ASV_8419" ">ASV_8420"
## [8421] ">ASV_8421" ">ASV_8422" ">ASV_8423" ">ASV_8424" ">ASV_8425" ">ASV_8426" ">ASV_8427" ">ASV_8428" ">ASV_8429" ">ASV_8430" ">ASV_8431" ">ASV_8432" ">ASV_8433" ">ASV_8434" ">ASV_8435" ">ASV_8436" ">ASV_8437" ">ASV_8438" ">ASV_8439" ">ASV_8440"
## [8441] ">ASV_8441" ">ASV_8442" ">ASV_8443" ">ASV_8444" ">ASV_8445" ">ASV_8446" ">ASV_8447" ">ASV_8448" ">ASV_8449" ">ASV_8450" ">ASV_8451" ">ASV_8452" ">ASV_8453" ">ASV_8454" ">ASV_8455" ">ASV_8456" ">ASV_8457" ">ASV_8458" ">ASV_8459" ">ASV_8460"
## [8461] ">ASV_8461" ">ASV_8462" ">ASV_8463" ">ASV_8464" ">ASV_8465" ">ASV_8466" ">ASV_8467" ">ASV_8468" ">ASV_8469" ">ASV_8470" ">ASV_8471" ">ASV_8472" ">ASV_8473" ">ASV_8474" ">ASV_8475" ">ASV_8476" ">ASV_8477" ">ASV_8478" ">ASV_8479" ">ASV_8480"
## [8481] ">ASV_8481" ">ASV_8482" ">ASV_8483" ">ASV_8484" ">ASV_8485" ">ASV_8486" ">ASV_8487" ">ASV_8488" ">ASV_8489" ">ASV_8490" ">ASV_8491" ">ASV_8492" ">ASV_8493" ">ASV_8494" ">ASV_8495" ">ASV_8496" ">ASV_8497" ">ASV_8498" ">ASV_8499" ">ASV_8500"
## [8501] ">ASV_8501" ">ASV_8502" ">ASV_8503" ">ASV_8504" ">ASV_8505" ">ASV_8506" ">ASV_8507" ">ASV_8508" ">ASV_8509" ">ASV_8510" ">ASV_8511" ">ASV_8512" ">ASV_8513" ">ASV_8514" ">ASV_8515" ">ASV_8516" ">ASV_8517" ">ASV_8518" ">ASV_8519" ">ASV_8520"
## [8521] ">ASV_8521" ">ASV_8522" ">ASV_8523" ">ASV_8524" ">ASV_8525" ">ASV_8526" ">ASV_8527" ">ASV_8528" ">ASV_8529" ">ASV_8530" ">ASV_8531" ">ASV_8532" ">ASV_8533" ">ASV_8534" ">ASV_8535" ">ASV_8536" ">ASV_8537" ">ASV_8538" ">ASV_8539" ">ASV_8540"
## [8541] ">ASV_8541" ">ASV_8542" ">ASV_8543" ">ASV_8544" ">ASV_8545" ">ASV_8546" ">ASV_8547" ">ASV_8548" ">ASV_8549" ">ASV_8550" ">ASV_8551" ">ASV_8552" ">ASV_8553" ">ASV_8554" ">ASV_8555" ">ASV_8556" ">ASV_8557" ">ASV_8558" ">ASV_8559" ">ASV_8560"
## [8561] ">ASV_8561" ">ASV_8562" ">ASV_8563" ">ASV_8564" ">ASV_8565" ">ASV_8566" ">ASV_8567" ">ASV_8568" ">ASV_8569" ">ASV_8570" ">ASV_8571" ">ASV_8572" ">ASV_8573" ">ASV_8574" ">ASV_8575" ">ASV_8576" ">ASV_8577" ">ASV_8578" ">ASV_8579" ">ASV_8580"
## [8581] ">ASV_8581" ">ASV_8582" ">ASV_8583" ">ASV_8584" ">ASV_8585" ">ASV_8586" ">ASV_8587" ">ASV_8588" ">ASV_8589" ">ASV_8590" ">ASV_8591" ">ASV_8592" ">ASV_8593" ">ASV_8594" ">ASV_8595" ">ASV_8596" ">ASV_8597" ">ASV_8598" ">ASV_8599" ">ASV_8600"
## [8601] ">ASV_8601" ">ASV_8602" ">ASV_8603" ">ASV_8604" ">ASV_8605" ">ASV_8606" ">ASV_8607" ">ASV_8608" ">ASV_8609" ">ASV_8610" ">ASV_8611" ">ASV_8612" ">ASV_8613" ">ASV_8614" ">ASV_8615" ">ASV_8616" ">ASV_8617" ">ASV_8618" ">ASV_8619" ">ASV_8620"
## [8621] ">ASV_8621" ">ASV_8622" ">ASV_8623" ">ASV_8624" ">ASV_8625" ">ASV_8626" ">ASV_8627" ">ASV_8628" ">ASV_8629" ">ASV_8630" ">ASV_8631" ">ASV_8632" ">ASV_8633" ">ASV_8634" ">ASV_8635" ">ASV_8636" ">ASV_8637" ">ASV_8638" ">ASV_8639" ">ASV_8640"
## [8641] ">ASV_8641" ">ASV_8642" ">ASV_8643" ">ASV_8644" ">ASV_8645" ">ASV_8646" ">ASV_8647" ">ASV_8648" ">ASV_8649" ">ASV_8650" ">ASV_8651" ">ASV_8652" ">ASV_8653" ">ASV_8654" ">ASV_8655" ">ASV_8656" ">ASV_8657" ">ASV_8658" ">ASV_8659" ">ASV_8660"
## [8661] ">ASV_8661" ">ASV_8662" ">ASV_8663" ">ASV_8664" ">ASV_8665" ">ASV_8666" ">ASV_8667" ">ASV_8668" ">ASV_8669" ">ASV_8670" ">ASV_8671" ">ASV_8672" ">ASV_8673" ">ASV_8674" ">ASV_8675" ">ASV_8676" ">ASV_8677" ">ASV_8678" ">ASV_8679" ">ASV_8680"
## [8681] ">ASV_8681" ">ASV_8682" ">ASV_8683" ">ASV_8684" ">ASV_8685" ">ASV_8686" ">ASV_8687" ">ASV_8688" ">ASV_8689" ">ASV_8690" ">ASV_8691" ">ASV_8692" ">ASV_8693" ">ASV_8694" ">ASV_8695" ">ASV_8696" ">ASV_8697" ">ASV_8698" ">ASV_8699" ">ASV_8700"
## [8701] ">ASV_8701" ">ASV_8702" ">ASV_8703" ">ASV_8704" ">ASV_8705" ">ASV_8706" ">ASV_8707" ">ASV_8708" ">ASV_8709" ">ASV_8710" ">ASV_8711" ">ASV_8712" ">ASV_8713" ">ASV_8714" ">ASV_8715" ">ASV_8716" ">ASV_8717" ">ASV_8718" ">ASV_8719" ">ASV_8720"
## [8721] ">ASV_8721" ">ASV_8722" ">ASV_8723" ">ASV_8724" ">ASV_8725" ">ASV_8726" ">ASV_8727" ">ASV_8728" ">ASV_8729" ">ASV_8730" ">ASV_8731" ">ASV_8732" ">ASV_8733" ">ASV_8734" ">ASV_8735" ">ASV_8736" ">ASV_8737" ">ASV_8738" ">ASV_8739" ">ASV_8740"
## [8741] ">ASV_8741" ">ASV_8742" ">ASV_8743" ">ASV_8744" ">ASV_8745" ">ASV_8746" ">ASV_8747" ">ASV_8748" ">ASV_8749" ">ASV_8750" ">ASV_8751" ">ASV_8752" ">ASV_8753" ">ASV_8754" ">ASV_8755" ">ASV_8756" ">ASV_8757" ">ASV_8758" ">ASV_8759" ">ASV_8760"
## [8761] ">ASV_8761" ">ASV_8762" ">ASV_8763" ">ASV_8764" ">ASV_8765" ">ASV_8766" ">ASV_8767" ">ASV_8768" ">ASV_8769" ">ASV_8770" ">ASV_8771" ">ASV_8772"
```

```r
##### Rename ASVs in table then write out our ASV fasta file! 
#View(seqtab_nochim)
asv_tab <- t(seqtab_nc_len)
#View(asv_tab)

## Rename our asvs! 
row.names(asv_tab) <- sub(">", "", asv_headers)
#View(asv_tab)

# Write the count table to a file! 
# You may need to create the 01_dada2_exports folder
write.table(asv_tab, "data/01_dada2_exports/ASV_counts.tsv", sep = "\t", quote = FALSE, col.names = NA)

# Write out the fasta file for reference later on for what seq matches what ASV
asv_fasta <- c(rbind(asv_headers, asv_seqs))

# Save to a file!
write(asv_fasta, "data/01_dada2_exports/ASVs.fasta")

# Also write to TaxAss folder
# You may need to make the TaxAss_analysis folder
write(asv_fasta, "data/TaxAss_analysis/ASVs.fasta")

# And save our seqtab_nc_len R Object for taxonomic assignment
save(seqtab_nc_len, file = "data/01_dada2_exports/seqtab_nc_len.RData")

# And save our meta_track R object
save(meta_track, file = "data/01_dada2_exports/meta_track.RData")
```


# Errata (not included in final outputs)

## Testing more stringent trimming - don't do!


```r
filtered_forward_reads_ag <- file.path("data", "filtered_aggressive", paste0(sample_names, "_R1_filtered.fastq.gz"))
filtered_reverse_reads_ag <- file.path("data", "filtered_aggressive", paste0(sample_names, "_R2_filtered.fastq.gz"))


filtered_out_ag <- filterAndTrim(forward_reads, filtered_forward_reads_ag,
                              reverse_reads, filtered_reverse_reads_ag,
                              truncLen = c(220,200), trimLeft = c(25,25),
                              maxN = 0, maxEE = c(1,1), truncQ = 2, 
                              rm.phix = TRUE, compress = TRUE, 
                              multithread = TRUE)

err_forward_reads_ag <- learnErrors(filtered_forward_reads_ag, multithread = TRUE) 

err_reverse_reads_ag <- learnErrors(filtered_reverse_reads_ag, multithread = TRUE)

dada_forward_ag <- dada(filtered_forward_reads_ag, err = err_forward_reads_ag, multithread = TRUE)

# run dada2 on the reverse sequences 
dada_reverse_ag <- dada(filtered_reverse_reads_ag, err = err_reverse_reads_ag, multithread = TRUE)


merged_amplicons_ag <- mergePairs(dada_forward_ag, filtered_forward_reads_ag, 
                               dada_reverse_ag, filtered_reverse_reads_ag,
                               verbose = TRUE)

seqtab_ag <- makeSequenceTable(merged_amplicons_ag)

dim(seqtab_ag) # 151 21008

data.frame(Seq_Length_ag = nchar(getSequences(seqtab_ag))) %>%
  ggplot(aes(x = Seq_Length_ag )) + 
  geom_histogram()

table(nchar(getSequences(seqtab_ag))) %>% sort()

seqtab_nochim_ag <- removeBimeraDenovo(seqtab_ag, verbose = TRUE) 
# Identified 1167 bimeras out of 21008 input sequences.

asv_keeps_ag <- nchar(getSequences(seqtab_nochim_ag)) %in% c(242, 243)


seqtab_nc_len_ag <- seqtab_nochim_ag[,asv_keeps_ag]

track_ag <- cbind(filtered_out_ag, 
               sapply(dada_forward_ag, getN),
               sapply(dada_reverse_ag, getN),
               sapply(merged_amplicons_ag, getN),
               rowSums(seqtab_nc_len_ag))


colnames(track_ag) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochim")
rownames(track_ag) <- sample_names


track_ag %>%
  # make it a dataframe
  as.data.frame() %>%
  rownames_to_column(var = "sample_name") %>%
  pivot_longer(input:nochim, names_to = "read_type", values_to = "num_reads") %>%
  #left_join(metadata, by = "sample_name") %>% 
  mutate(read_type = fct_relevel(read_type, 
                                 "input", "filtered", "denoisedF", "denoisedR", "merged", "nochim"),
         is_blank = sample_name%in%c("PCR_Blanks","AP_D61","AP_D62","AP_D153","AP_D154")) %>%
  ggplot(aes(x = read_type, y = num_reads, fill = read_type, color = is_blank)) + 
  #facet_grid(~strata) + 
  geom_line(aes(group = sample_name)) + 
  geom_point(shape = 21, size = 3, alpha = 0.8, color = "grey") + 
  scale_fill_brewer(palette = "Spectral") + 
  theme_bw() + 
  labs(x = "Filtering Step", y = "Number of Sequences") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


# Check the mock commmunity 
mock_sample_ag <- seqtab_nc_len_ag["Zymo_Mock_R1_filtered.fastq.gz",]

length(mock_sample_ag)

# Drop ASVs absent from mock community 
length(mock_sample_ag[mock_sample_ag > 0])

mock_sample_sub_ag <- sort(mock_sample_ag[mock_sample_ag > 0], decreasing = TRUE)

length(mock_sample_sub_ag)

cat("DADA2 inferred", length(mock_sample_sub_ag), "ASVs present in the Mock Community.")

#Who are they in the mock community? 

#### Compare our ASVs from the mock community to the reference fasta!
mock_reference_ag <- getSequences(file.path("data/mock_fastas/", "mock_amplicons.fasta"))

matches_ag <- sapply(names(mock_sample_sub),
                             function(x) any(grepl(x, mock_reference)))

match_mock_ref_ag <- sum(matches_ag)

total_mock_reads_ag <- sum(mock_sample_sub_ag)
matched_reads_ag <- sum(mock_sample_sub_ag[matches_ag])
err_reads_ag <- sum(mock_sample_sub_ag[!matches_ag])


cat(sum(match_mock_ref_ag), "ASVs were exact matches to the expected reference sequences.")

cat("DADA2 inferred ", length(mock_sample_sub_ag) - match_mock_ref_ag, "erroneous ASVs out of the mock community with aggressive filtering")

cat("Erroneous matches represented ", 100*err_reads_ag/total_mock_reads_ag, "% of total reads in the mock with aggressive filtering")
```


# Session Info


```r
sessioninfo::session_info()
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
##  ! package              * version    date (UTC) lib source
##  P abind                  1.4-5      2016-07-21 [?] CRAN (R 4.3.2)
##  P ade4                   1.7-22     2023-02-06 [?] CRAN (R 4.3.2)
##  P ape                    5.7-1      2023-03-13 [?] CRAN (R 4.3.2)
##  P Biobase                2.62.0     2023-10-24 [?] Bioconductor
##  P BiocGenerics         * 0.48.1     2023-11-01 [?] Bioconductor
##  P BiocManager            1.30.22    2023-08-08 [?] CRAN (R 4.3.2)
##  P BiocParallel           1.36.0     2023-10-24 [?] Bioconductor
##  P biomformat             1.30.0     2023-10-24 [?] Bioconductor
##  P Biostrings           * 2.70.1     2023-10-25 [?] Bioconductor
##  P bit                    4.0.5      2022-11-15 [?] CRAN (R 4.3.2)
##  P bit64                  4.0.5      2020-08-30 [?] CRAN (R 4.3.2)
##  P bitops                 1.0-7      2021-04-24 [?] CRAN (R 4.3.2)
##  P bslib                  0.5.1      2023-08-11 [?] CRAN (R 4.3.2)
##  P cachem                 1.0.8      2023-05-01 [?] CRAN (R 4.3.2)
##  P cli                    3.6.1      2023-03-23 [?] CRAN (R 4.3.2)
##  P cluster                2.1.4      2022-08-22 [?] CRAN (R 4.3.2)
##  P codetools              0.2-19     2023-02-01 [?] CRAN (R 4.3.3)
##  P colorspace             2.1-0      2023-01-23 [?] CRAN (R 4.3.2)
##  P crayon                 1.5.2      2022-09-29 [?] CRAN (R 4.3.2)
##  P cytolib                2.14.1     2024-01-18 [?] Bioconduc~
##  P dada2                * 1.30.0     2023-10-24 [?] Bioconductor
##  P data.table             1.15.2     2024-02-29 [?] CRAN (R 4.3.2)
##  P DelayedArray           0.28.0     2023-10-24 [?] Bioconductor
##  P deldir                 2.0-4      2024-02-28 [?] CRAN (R 4.3.2)
##  P digest                 0.6.33     2023-07-07 [?] CRAN (R 4.3.2)
##  P dplyr                * 1.1.3      2023-09-03 [?] CRAN (R 4.3.2)
##  P ellipsis               0.3.2      2021-04-29 [?] CRAN (R 4.3.2)
##  P evaluate               0.23       2023-11-01 [?] CRAN (R 4.3.2)
##  P fansi                  1.0.5      2023-10-08 [?] CRAN (R 4.3.2)
##  P farver                 2.1.1      2022-07-06 [?] CRAN (R 4.3.2)
##  P fastmap                1.1.1      2023-02-24 [?] CRAN (R 4.3.2)
##  P flowCore               2.14.2     2024-03-18 [?] Bioconduc~
##  P forcats              * 1.0.0      2023-01-29 [?] CRAN (R 4.3.2)
##  P foreach                1.5.2      2022-02-02 [?] CRAN (R 4.3.2)
##  P fun.gus              * 0.3.1      2025-06-26 [?] Github (MarschmiLab/fun.gus@7daa3fa)
##  P furrr                  0.3.1      2022-08-15 [?] CRAN (R 4.3.2)
##  P future                 1.33.1     2023-12-22 [?] CRAN (R 4.3.3)
##  P generics               0.1.3      2022-07-05 [?] CRAN (R 4.3.2)
##  P GenomeInfoDb         * 1.38.0     2023-10-24 [?] Bioconductor
##  P GenomeInfoDbData       1.2.11     2023-11-07 [?] Bioconductor
##  P GenomicAlignments      1.38.2     2024-01-16 [?] Bioconduc~
##  P GenomicRanges          1.54.1     2023-10-29 [?] Bioconductor
##  P ggplot2              * 3.5.0      2024-02-23 [?] CRAN (R 4.3.2)
##  P globals                0.16.2     2022-11-21 [?] CRAN (R 4.3.3)
##  P glue                   1.6.2      2022-02-24 [?] CRAN (R 4.3.2)
##  P gtable                 0.3.4      2023-08-21 [?] CRAN (R 4.3.2)
##  P highr                  0.10       2022-12-22 [?] CRAN (R 4.3.2)
##  P hms                    1.1.3      2023-03-21 [?] CRAN (R 4.3.2)
##  P htmltools              0.5.7      2023-11-03 [?] CRAN (R 4.3.2)
##  P httpuv                 1.6.12     2023-10-23 [?] CRAN (R 4.3.2)
##  P hwriter                1.3.2.1    2022-04-08 [?] CRAN (R 4.3.2)
##  P igraph                 1.5.1      2023-08-10 [?] CRAN (R 4.3.2)
##  P interp                 1.1-6      2024-01-26 [?] CRAN (R 4.3.2)
##  P IRanges              * 2.36.0     2023-10-24 [?] Bioconductor
##  P iterators              1.0.14     2022-02-05 [?] CRAN (R 4.3.2)
##  P jpeg                   0.1-10     2022-11-29 [?] CRAN (R 4.3.2)
##  P jquerylib              0.1.4      2021-04-26 [?] CRAN (R 4.3.2)
##  P jsonlite               1.8.7      2023-06-29 [?] CRAN (R 4.3.2)
##  P knitr                  1.45       2023-10-30 [?] CRAN (R 4.3.2)
##  P labeling               0.4.3      2023-08-29 [?] CRAN (R 4.3.2)
##  P later                  1.3.1      2023-05-02 [?] CRAN (R 4.3.2)
##  P lattice                0.21-9     2023-10-01 [?] CRAN (R 4.3.2)
##  P latticeExtra           0.6-30     2022-07-04 [?] CRAN (R 4.3.2)
##  P lifecycle              1.0.3      2022-10-07 [?] CRAN (R 4.3.2)
##  P listenv                0.9.1      2024-01-29 [?] CRAN (R 4.3.2)
##  P lubridate            * 1.9.3      2023-09-27 [?] CRAN (R 4.3.2)
##  P magrittr               2.0.3      2022-03-30 [?] CRAN (R 4.3.2)
##  P MASS                   7.3-60     2023-05-04 [?] CRAN (R 4.3.2)
##  P Matrix                 1.6-1.1    2023-09-18 [?] CRAN (R 4.3.2)
##  P MatrixGenerics         1.14.0     2023-10-24 [?] Bioconductor
##  P matrixStats            1.2.0      2023-12-11 [?] CRAN (R 4.3.2)
##  P mgcv                   1.9-0      2023-07-11 [?] CRAN (R 4.3.2)
##  P mime                   0.12       2021-09-28 [?] CRAN (R 4.3.2)
##  P multtest               2.58.0     2023-10-24 [?] Bioconductor
##  P munsell                0.5.0      2018-06-12 [?] CRAN (R 4.3.2)
##  P NatParksPalettes     * 0.2.0      2022-10-09 [?] CRAN (R 4.3.2)
##  P nlme                   3.1-163    2023-08-09 [?] CRAN (R 4.3.2)
##  P pacman                 0.5.1      2019-03-11 [?] CRAN (R 4.3.2)
##  P parallelly             1.37.1     2024-02-29 [?] CRAN (R 4.3.3)
##  P patchwork            * 1.2.0.9000 2025-06-26 [?] Github (thomasp85/patchwork@d943757)
##  P permute                0.9-7      2022-01-27 [?] CRAN (R 4.3.2)
##  P phyloseq             * 1.46.0     2023-10-24 [?] Bioconductor
##  P pillar                 1.9.0      2023-03-22 [?] CRAN (R 4.3.2)
##  P pkgconfig              2.0.3      2019-09-22 [?] CRAN (R 4.3.2)
##  P plyr                   1.8.9      2023-10-02 [?] CRAN (R 4.3.2)
##  P png                    0.1-8      2022-11-29 [?] CRAN (R 4.3.2)
##  P promises               1.2.1      2023-08-10 [?] CRAN (R 4.3.2)
##  P purrr                * 1.0.2      2023-08-10 [?] CRAN (R 4.3.2)
##  P R6                     2.5.1      2021-08-19 [?] CRAN (R 4.3.2)
##  P RColorBrewer           1.1-3      2022-04-03 [?] CRAN (R 4.3.2)
##  P Rcpp                 * 1.0.11     2023-07-06 [?] CRAN (R 4.3.2)
##  P RcppParallel           5.1.7      2023-02-27 [?] CRAN (R 4.3.2)
##  P RCurl                  1.98-1.13  2023-11-02 [?] CRAN (R 4.3.2)
##  P readr                * 2.1.5      2024-01-10 [?] CRAN (R 4.3.2)
##    renv                   1.0.5      2024-02-29 [1] CRAN (R 4.3.2)
##  P reshape2               1.4.4      2020-04-09 [?] CRAN (R 4.3.2)
##  P rhdf5                  2.46.1     2023-11-29 [?] Bioconduc~
##  P rhdf5filters           1.14.1     2023-11-06 [?] Bioconductor
##  P Rhdf5lib               1.24.2     2024-02-07 [?] Bioconduc~
##  P rlang                  1.1.2      2023-11-04 [?] CRAN (R 4.3.2)
##  P rmarkdown              2.25       2023-09-18 [?] CRAN (R 4.3.2)
##  P RProtoBufLib           2.14.1     2024-03-18 [?] Bioconduc~
##  P Rsamtools              2.18.0     2023-10-24 [?] Bioconductor
##  P rstudioapi             0.15.0     2023-07-07 [?] CRAN (R 4.3.2)
##  P S4Arrays               1.2.0      2023-10-24 [?] Bioconductor
##  P S4Vectors            * 0.40.1     2023-10-26 [?] Bioconductor
##  P sass                   0.4.7      2023-07-15 [?] CRAN (R 4.3.2)
##  P scales                 1.3.0      2023-11-28 [?] CRAN (R 4.3.2)
##  P sessioninfo            1.2.2      2021-12-06 [?] CRAN (R 4.3.2)
##  P shiny                  1.7.5.1    2023-10-14 [?] CRAN (R 4.3.2)
##  P ShortRead              1.60.0     2023-10-24 [?] Bioconductor
##  P SparseArray            1.2.4      2024-02-11 [?] Bioconduc~
##  P stringi                1.7.12     2023-01-11 [?] CRAN (R 4.3.2)
##  P stringr              * 1.5.0      2022-12-02 [?] CRAN (R 4.3.2)
##  P SummarizedExperiment   1.32.0     2023-10-24 [?] Bioconductor
##  P survival               3.5-8      2024-02-14 [?] CRAN (R 4.3.3)
##  P tibble               * 3.2.1      2023-03-20 [?] CRAN (R 4.3.2)
##  P tidyr                * 1.3.1      2024-01-24 [?] CRAN (R 4.3.2)
##  P tidyselect             1.2.0      2022-10-10 [?] CRAN (R 4.3.2)
##  P tidyverse            * 2.0.0      2023-02-22 [?] CRAN (R 4.3.2)
##  P timechange             0.3.0      2024-01-18 [?] CRAN (R 4.3.2)
##  P tzdb                   0.4.0      2023-05-12 [?] CRAN (R 4.3.2)
##  P utf8                   1.2.4      2023-10-22 [?] CRAN (R 4.3.2)
##  P vctrs                  0.6.4      2023-10-12 [?] CRAN (R 4.3.2)
##  P vegan                  2.6-4      2022-10-11 [?] CRAN (R 4.3.2)
##  P vroom                  1.6.5      2023-12-05 [?] CRAN (R 4.3.2)
##  P withr                  2.5.2      2023-10-30 [?] CRAN (R 4.3.2)
##  P xfun                   0.52       2025-04-02 [?] CRAN (R 4.3.3)
##  P xtable                 1.8-4      2019-04-21 [?] CRAN (R 4.3.2)
##  P XVector              * 0.42.0     2023-10-24 [?] Bioconductor
##  P yaml                   2.3.7      2023-01-23 [?] CRAN (R 4.3.2)
##  P zlibbioc               1.48.0     2023-10-24 [?] Bioconductor
## 
##  [1] /local/workdir/arp277/Pendleton_2025_Ontario_Publication_Repo/renv/library/R-4.3/x86_64-pc-linux-gnu
##  [2] /home/arp277/.cache/R/renv/sandbox/R-4.3/x86_64-pc-linux-gnu/fd835031
## 
##  P ── Loaded and on-disk path mismatch.
## 
## ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
```
