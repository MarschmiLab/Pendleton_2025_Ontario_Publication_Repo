---
title: "Functional_Prediction" 
author: "Augustus Pendleton"
date: "27 June, 2025"
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




# Load packages 

```r
# Efficiently load packages 
pacman::p_load(phyloseq, tidyverse, install = FALSE)

knitr::write_bib(file = "data/11_functional_exports/packages.bib")

# load in functions and color preferences
source("code/R/plotting_aesthetics.R")
```


# Load Data


```r
load("data/08_compositional_exports/full_abs_physeq.RData")
```


## Functional annotation with FAPROTAX

Remember for bash chunks, always begin in the project root directory - I'll be explicit with cd if I'm navigating elsewhere


```bash

cd data/11_functional_exports/faprotax

wget https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/SECTION_Download/MODULE_Downloads/CLASS_Latest%20release/UNIT_FAPROTAX_1.2.10/FAPROTAX_1.2.10.zip

unzip FAPROTAX_1.2.10.zip
```


```r
# We will do this with cell-abundance normalized counts, going in. 
# We are using the silva assignments; as it is the normalized taxonomy used by FAPROTAX (though it makes almost no difference)

both <- read_tsv("data/02_taxonomy_exports/ASV_both_taxonomy.tsv")

tax_names <- both %>%
  mutate(taxonomy = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = "; ")) %>%
  select(ASV, taxonomy)

input_table <- full_abs_physeq %>%
  otu_table %>%
  psmelt() %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>%
  dplyr::rename(ASV = OTU) %>%
  left_join(tax_names)


write_tsv(input_table, file = "data/11_functional_exports/faprotax/input_table.tsv")
```


```bash

# CONFIRM YOU ARE IN THE faprotax DIRECTORY!

mkdir -p faprotax_outputs

conda create -n "python_for_faprotax" python=3.7.6 ipython

conda activate python_for_faprotax
conda install numpy

FAPROTAX_1.2.10/collapse_table.py -i input_table.tsv -g FAPROTAX_1.2.10/FAPROTAX.txt -d taxonomy --omit_columns 0 -f -o faprotax_outputs/func_table.tsv -r faprotax_outputs/fapro_report.txt

cd ../../../

```

How many taxa were assigned at least one function?


```r
readLines("data/11_functional_exports/faprotax/faprotax_outputs/fapro_report.txt")[107]
```

```
## [1] "# 1764 out of 7278 records (24.2374 %) were assigned to at least one group"
```

```r
readLines("data/11_functional_exports/faprotax/faprotax_outputs/fapro_report.txt")[108]
```

```
## [1] "# 5514 out of 7278 records (75.7626 %) could not be assigned to any group (leftovers)"
```

This is an important caveat here - so many taxa are unannotated!

# Writing out select FAPROTAX functions for plotting in QGIS

We don't do any visualization in R here - hence the writing out. 


```r
fapro_res <- read_tsv(file = "data/11_functional_exports/faprotax/faprotax_outputs/func_table.tsv")

by_rep_fapro <- fapro_res %>%
  pivot_longer(!group, names_to = "Rep_ID", values_to = "Count") %>%
  group_by(group) %>% # Group by functional annotation
  dplyr::filter(any(Count!=0)) %>% # Get rid of any groups that are all zero
  ungroup()

for_join <- sample_data(full_abs_physeq) %>%
  data.frame() %>%
  select(Rep_ID, Comp_Group_Hier)

semi_clean_fapro <- by_rep_fapro %>%
  left_join(sample_data(full_abs_physeq) %>%
  data.frame() %>%
  select(Rep_ID, month, Depth_Class, Latitude, Longitude)) %>%
  filter(Depth_Class!= "M") %>%
  filter(group %in% c("methanotrophy","aerobic_ammonia_oxidation", "sulfate_respiration","photoautotrophy")) %>%
  pivot_wider(names_from = "group", values_from = "Count")

semi_clean_fapro %>%
  sf::st_as_sf(coords = c("Longitude","Latitude"), crs = "EPSG:4326") %>%
  filter(month == "September", Depth_Class == "B") %>%
  sf::st_transform(crs = "EPSG:3174")%>%
  sf::st_write("analysis/QGIS_Work/sep_b_mt.gpkg",
               append = FALSE)
```

```
## Deleting layer `sep_b_mt' using driver `GPKG'
## Writing layer `sep_b_mt' to data source `analysis/QGIS_Work/sep_b_mt.gpkg' using driver `GPKG'
## Writing 15 features with 7 fields and geometry type Point.
```

```r
semi_clean_fapro %>%
  sf::st_as_sf(coords = c("Longitude","Latitude"), crs = "EPSG:4326") %>%
  filter(month == "September", Depth_Class == "E") %>%
  sf::st_transform(crs = "EPSG:3174")%>%
  sf::st_write("analysis/QGIS_Work/sep_e_mt.gpkg",
               append = FALSE)
```

```
## Deleting layer `sep_e_mt' using driver `GPKG'
## Writing layer `sep_e_mt' to data source `analysis/QGIS_Work/sep_e_mt.gpkg' using driver `GPKG'
## Writing 15 features with 7 fields and geometry type Point.
```

```r
semi_clean_fapro %>%
  sf::st_as_sf(coords = c("Longitude","Latitude"), crs = "EPSG:4326") %>%
  filter(month == "May", Depth_Class == "B") %>%
  sf::st_transform(crs = "EPSG:3174")%>%
  sf::st_write("analysis/QGIS_Work/may_b_mt.gpkg",
               append = FALSE)
```

```
## Deleting layer `may_b_mt' using driver `GPKG'
## Writing layer `may_b_mt' to data source `analysis/QGIS_Work/may_b_mt.gpkg' using driver `GPKG'
## Writing 9 features with 7 fields and geometry type Point.
```

```r
semi_clean_fapro %>%
  sf::st_as_sf(coords = c("Longitude","Latitude"), crs = "EPSG:4326") %>%
  filter(month == "May", Depth_Class == "E") %>%
  sf::st_transform(crs = "EPSG:3174")%>%
  sf::st_write("analysis/QGIS_Work/may_e_mt.gpkg",
               append = FALSE)
```

```
## Deleting layer `may_e_mt' using driver `GPKG'
## Writing layer `may_e_mt' to data source `analysis/QGIS_Work/may_e_mt.gpkg' using driver `GPKG'
## Writing 15 features with 7 fields and geometry type Point.
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
##  date     2025-06-27
##  pandoc   3.1.1 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/ (via rmarkdown)
## 
## ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
##  ! package          * version   date (UTC) lib source
##  P ade4               1.7-22    2023-02-06 [?] CRAN (R 4.3.2)
##  P ape                5.7-1     2023-03-13 [?] CRAN (R 4.3.2)
##  P Biobase            2.62.0    2023-10-24 [?] Bioconductor
##  P BiocGenerics       0.48.1    2023-11-01 [?] Bioconductor
##  P BiocManager        1.30.22   2023-08-08 [?] CRAN (R 4.3.2)
##  P biomformat         1.30.0    2023-10-24 [?] Bioconductor
##  P Biostrings         2.70.1    2023-10-25 [?] Bioconductor
##  P bit                4.0.5     2022-11-15 [?] CRAN (R 4.3.2)
##  P bit64              4.0.5     2020-08-30 [?] CRAN (R 4.3.2)
##  P bitops             1.0-7     2021-04-24 [?] CRAN (R 4.3.2)
##  P bslib              0.5.1     2023-08-11 [?] CRAN (R 4.3.2)
##  P cachem             1.0.8     2023-05-01 [?] CRAN (R 4.3.2)
##  P class              7.3-22    2023-05-03 [?] CRAN (R 4.3.3)
##  P classInt           0.4-10    2023-09-05 [?] CRAN (R 4.3.2)
##  P cli                3.6.1     2023-03-23 [?] CRAN (R 4.3.2)
##  P cluster            2.1.4     2022-08-22 [?] CRAN (R 4.3.2)
##  P codetools          0.2-19    2023-02-01 [?] CRAN (R 4.3.3)
##  P colorspace         2.1-0     2023-01-23 [?] CRAN (R 4.3.2)
##  P crayon             1.5.2     2022-09-29 [?] CRAN (R 4.3.2)
##  P data.table         1.15.2    2024-02-29 [?] CRAN (R 4.3.2)
##  P DBI                1.2.2     2024-02-16 [?] CRAN (R 4.3.2)
##  P digest             0.6.33    2023-07-07 [?] CRAN (R 4.3.2)
##  P dplyr            * 1.1.3     2023-09-03 [?] CRAN (R 4.3.2)
##  P e1071              1.7-14    2023-12-06 [?] CRAN (R 4.3.2)
##  P evaluate           0.23      2023-11-01 [?] CRAN (R 4.3.2)
##  P fansi              1.0.5     2023-10-08 [?] CRAN (R 4.3.2)
##  P fastmap            1.1.1     2023-02-24 [?] CRAN (R 4.3.2)
##  P forcats          * 1.0.0     2023-01-29 [?] CRAN (R 4.3.2)
##  P foreach            1.5.2     2022-02-02 [?] CRAN (R 4.3.2)
##  P generics           0.1.3     2022-07-05 [?] CRAN (R 4.3.2)
##  P GenomeInfoDb       1.38.0    2023-10-24 [?] Bioconductor
##  P GenomeInfoDbData   1.2.11    2023-11-07 [?] Bioconductor
##  P ggplot2          * 3.5.0     2024-02-23 [?] CRAN (R 4.3.2)
##  P glue               1.6.2     2022-02-24 [?] CRAN (R 4.3.2)
##  P gtable             0.3.4     2023-08-21 [?] CRAN (R 4.3.2)
##  P hms                1.1.3     2023-03-21 [?] CRAN (R 4.3.2)
##  P htmltools          0.5.7     2023-11-03 [?] CRAN (R 4.3.2)
##  P igraph             1.5.1     2023-08-10 [?] CRAN (R 4.3.2)
##  P IRanges            2.36.0    2023-10-24 [?] Bioconductor
##  P iterators          1.0.14    2022-02-05 [?] CRAN (R 4.3.2)
##  P jquerylib          0.1.4     2021-04-26 [?] CRAN (R 4.3.2)
##  P jsonlite           1.8.7     2023-06-29 [?] CRAN (R 4.3.2)
##  P KernSmooth         2.23-22   2023-07-10 [?] CRAN (R 4.3.3)
##  P knitr              1.45      2023-10-30 [?] CRAN (R 4.3.2)
##  P lattice            0.21-9    2023-10-01 [?] CRAN (R 4.3.2)
##  P lifecycle          1.0.3     2022-10-07 [?] CRAN (R 4.3.2)
##  P lubridate        * 1.9.3     2023-09-27 [?] CRAN (R 4.3.2)
##  P magrittr           2.0.3     2022-03-30 [?] CRAN (R 4.3.2)
##  P MASS               7.3-60    2023-05-04 [?] CRAN (R 4.3.2)
##  P Matrix             1.6-1.1   2023-09-18 [?] CRAN (R 4.3.2)
##  P mgcv               1.9-0     2023-07-11 [?] CRAN (R 4.3.2)
##  P multtest           2.58.0    2023-10-24 [?] Bioconductor
##  P munsell            0.5.0     2018-06-12 [?] CRAN (R 4.3.2)
##  P NatParksPalettes * 0.2.0     2022-10-09 [?] CRAN (R 4.3.2)
##  P nlme               3.1-163   2023-08-09 [?] CRAN (R 4.3.2)
##  P pacman             0.5.1     2019-03-11 [?] CRAN (R 4.3.2)
##  P permute            0.9-7     2022-01-27 [?] CRAN (R 4.3.2)
##  P phyloseq         * 1.46.0    2023-10-24 [?] Bioconductor
##  P pillar             1.9.0     2023-03-22 [?] CRAN (R 4.3.2)
##  P pkgconfig          2.0.3     2019-09-22 [?] CRAN (R 4.3.2)
##  P plyr               1.8.9     2023-10-02 [?] CRAN (R 4.3.2)
##  P proxy              0.4-27    2022-06-09 [?] CRAN (R 4.3.2)
##  P purrr            * 1.0.2     2023-08-10 [?] CRAN (R 4.3.2)
##  P R6                 2.5.1     2021-08-19 [?] CRAN (R 4.3.2)
##  P Rcpp               1.0.11    2023-07-06 [?] CRAN (R 4.3.2)
##  P RCurl              1.98-1.13 2023-11-02 [?] CRAN (R 4.3.2)
##  P readr            * 2.1.5     2024-01-10 [?] CRAN (R 4.3.2)
##    renv               1.0.5     2024-02-29 [1] CRAN (R 4.3.2)
##  P reshape2           1.4.4     2020-04-09 [?] CRAN (R 4.3.2)
##  P rhdf5              2.46.1    2023-11-29 [?] Bioconduc~
##  P rhdf5filters       1.14.1    2023-11-06 [?] Bioconductor
##  P Rhdf5lib           1.24.2    2024-02-07 [?] Bioconduc~
##  P rlang              1.1.2     2023-11-04 [?] CRAN (R 4.3.2)
##  P rmarkdown          2.25      2023-09-18 [?] CRAN (R 4.3.2)
##  P rstudioapi         0.15.0    2023-07-07 [?] CRAN (R 4.3.2)
##  P S4Vectors          0.40.1    2023-10-26 [?] Bioconductor
##  P sass               0.4.7     2023-07-15 [?] CRAN (R 4.3.2)
##  P scales             1.3.0     2023-11-28 [?] CRAN (R 4.3.2)
##  P sessioninfo        1.2.2     2021-12-06 [?] CRAN (R 4.3.2)
##  P sf                 1.0-15    2023-12-18 [?] CRAN (R 4.3.2)
##  P stringi            1.7.12    2023-01-11 [?] CRAN (R 4.3.2)
##  P stringr          * 1.5.0     2022-12-02 [?] CRAN (R 4.3.2)
##  P survival           3.5-8     2024-02-14 [?] CRAN (R 4.3.3)
##  P tibble           * 3.2.1     2023-03-20 [?] CRAN (R 4.3.2)
##  P tidyr            * 1.3.1     2024-01-24 [?] CRAN (R 4.3.2)
##  P tidyselect         1.2.0     2022-10-10 [?] CRAN (R 4.3.2)
##  P tidyverse        * 2.0.0     2023-02-22 [?] CRAN (R 4.3.2)
##  P timechange         0.3.0     2024-01-18 [?] CRAN (R 4.3.2)
##  P tzdb               0.4.0     2023-05-12 [?] CRAN (R 4.3.2)
##  P units              0.8-5     2023-11-28 [?] CRAN (R 4.3.2)
##  P utf8               1.2.4     2023-10-22 [?] CRAN (R 4.3.2)
##  P vctrs              0.6.4     2023-10-12 [?] CRAN (R 4.3.2)
##  P vegan              2.6-4     2022-10-11 [?] CRAN (R 4.3.2)
##  P vroom              1.6.5     2023-12-05 [?] CRAN (R 4.3.2)
##  P withr              2.5.2     2023-10-30 [?] CRAN (R 4.3.2)
##  P xfun               0.52      2025-04-02 [?] CRAN (R 4.3.3)
##  P XVector            0.42.0    2023-10-24 [?] Bioconductor
##  P yaml               2.3.7     2023-01-23 [?] CRAN (R 4.3.2)
##  P zlibbioc           1.48.0    2023-10-24 [?] Bioconductor
## 
##  [1] /local/workdir/arp277/Pendleton_2025_Ontario_Publication_Repo/renv/library/R-4.3/x86_64-pc-linux-gnu
##  [2] /home/arp277/.cache/R/renv/sandbox/R-4.3/x86_64-pc-linux-gnu/fd835031
## 
##  P ── Loaded and on-disk path mismatch.
## 
## ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
```
