## Code and Data for "Upwelling periodically disturbs the ecological assembly of microbial communities in the Laurentian Great Lakes"

Augustus Pendleton, Mathew G. Wells, and Marian Schmidt

#### DOI Link for this Repo: 

Please find the citeable DOI for this repo on [Zenodo](https://doi.org/10.5281/zenodo.15793805).

#### Link to paper:

Read the manuscript on [BioRxiv](https://www.biorxiv.org/content/10.1101/2025.01.17.633667v1) or access the Manuscript.pdf and SI_Appendix.pdf directly from this repo.

#### Abstract

The Laurentian Great Lakes hold 21% of the world’s surface freshwater and supply drinking water to nearly 40 million people. Here we provide the first evidence that wind-driven upwelling fundamentally restructures microbial communities in Lake Ontario, with its effects extended and redistributed by an internal Kelvin wave propagating along the shoreline for >2 weeks. While thermal stratification is known to organize microbial communities by depth and season, we show that this vertical structure arises from contrasting mechanisms: homogenizing selection in surface waters and dispersal limitation with drift in the hypolimnion. Kelvin wave-driven upwelling disrupts this scaffold, displacing rare taxa into the surface and creating novel coastal communities enriched in methane oxidation and sulfur metabolism genes—functional traits absent elsewhere in the lake. We observed a Kelvin wave lasting over two weeks and propagating eastward at \~60 km day⁻¹. Given the \~10–12 day recurrence of wind events, at any time, at least one segment of Lake Ontario’s coastline is typically experiencing upwelling, producing pulses frequent and sustained enough to remodel microbial communities on ecologically relevant timescales. Recurrent upwellings, sustained and redistributed by Kelvin waves, act as a biological disturbance that overrides stratification, mobilizes rare functional potential, and assembles novel coastal microbial communities. As climate change lengthens stratified periods and reshapes large-lake circulation, understanding how physical forcing governs microbial assembly is essential for forecasting the biogeochemical future of Earth’s great lakes, especially in shoreline zones where ecological shifts directly affect human communities.

#### Significance Statement

The Laurentian Great Lakes hold 21% of the world’s surface freshwater and nearly 85% of North America’s—supplying drinking water to over 40 million people. Yet the microbial life that underpins their ecological function remains poorly understood. We show that wind-driven upwelling, followed by internal Kelvin waves that redistribute and sustain upwelling, acts as a recurring physical disturbance that overrides thermal stratification, redistributing rare microbes and assembling novel shoreline communities. These shifts bring unexpected functional changes, including enrichment of methane oxidation and sulfur metabolism genes, where people most interact with lake water. As climate change reshapes large lake circulation and intensifies thermal layering, understanding how microbial ecosystems respond is essential for forecasting transformations in water quality, ecosystem function, and biogeochemical resilience.

#### Description of this repo:

This directory contains the data and code for the first manuscript coming out of the 2023 Ontario CSMI sampling trip. Herein, I amplified the 16S V4 region for May and September samples across Lake Ontario (and over three depths). Note, raw sequencing files and flow cytometry .fcs files are not present within this repo, though can be downloaded from SRA (PRJNA1212049) and FlowRepository (FR-FCM-Z8SJ), respectively. The analysis folder contains .Rmd files wherein I analyze the data and produce figures, and the knitted .md and .html outputs. Note that html outputs won't display in Github - clone this repo onto your device and open from them. The data folder holds exports (and inputs) of each step in the analysis pipeline. Large outputs may not be present (check the .gitignore), but the code to produce them is present. The code folder holds stand-alone R scripts that either load functions or plotting aesthetics. The figures folder holds exported figures.
