# This file indicates that colors, shapes, and themes across my documents
# let's keep the visualizations consistent! :) 
library(ggplot2)
library(NatParksPalettes)


ont_theme <- theme_classic() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size=12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        strip.background = element_rect(colour = NA, fill = 'transparent'),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )

theme_set(ont_theme)


# Making touching the x-axis step shorter
y_ex <- scale_y_continuous(expand = c(0.0,0))

# Defining some colors
acadia <- natparks.pals(name = "Acadia", n = 9)

arches <- natparks.pals(name = "Arches", n = 9)

bryce <- natparks.pals(name = "BryceCanyon", n = 9)

triglav <- natparks.pals(name = "Triglav", n = 6)

capitol_reef <- natparks.pals(name = "CapitolReef", n = 6)

volcanoes <- natparks.pals(name = "Volcanoes", n = 5)

month_shapes <- c(
  "May" = 21,
  "September" = 24)

upwelling_shapes <- c(
  "Other" = 21,
  "Upwelling" = 24,
  "Downwelling" = 25,
  "Welland Canal" = 22
)
# Month colors from 
month_colors <- c(
  "May" = arches[2],
  "September" = arches[8]
)

sampletype_colors <- c(
  "Sample" = bryce[2],
  "Blank" = bryce[5])


depth_colors <- c(
  "E" = triglav[2],
  "M" = triglav[1],
  "B" = triglav[6])

comp_group_colors <- c(
  "Shallow_May" = triglav[5],
  "Shallow_September" = triglav[4],
  "Deep_May" = triglav[6],
  "Deep_September" = triglav[3])

comp_group_colors_hier <- c(
  "Shallow_May" = triglav[5],
  "Shallow_September" = triglav[4],
  "Deep (May)" = triglav[6],
  "Deep (September)" = triglav[3])

comp_three_colors <- c(
  "Shallow_May" = triglav[5],
  "Shallow_September" = triglav[4],
  "Deep" = triglav[6]
)
station_depth_colors <- c(
  "Nearshore" = triglav[3],
  "Deep" = triglav[1]
)

transect_colors <- volcanoes
names(transect_colors) <- c("T1","T2","T3","T4","T5")

north_vs_south_colors <- c(
  "North" = volcanoes[2],
  "Center" = volcanoes[3],
  "South" = volcanoes[4],
  "Weird Stn 74" = volcanoes[5]
)
# Set station colors
station_colors <- data.frame(Station = c(8,  12,  17,  29,  33,  35,  38,  41,  43,  48 , 55,  62,  64,  66 , 74, 717),
                             Station_Color = c(natparks.pals("DeathValley", n = 8), natparks.pals("CraterLake", n = 8)))

# Set the phylum colors
phylum_colors <- c(
  Acidobacteriota = "navy", 
  Actinobacteriota = "darkslategray2", 
  Armatimonadota = "deeppink1",
  Alphaproteobacteria = "plum2", 
  Bacteroidota = "gold", 
  Betaproteobacteria = "plum1", 
  Bdellovibrionota = "red1",
  Chloroflexi="black", 
  Crenarchaeota = "plum1",
  Cyanobacteria = "limegreen",
  Campylobacterota = "grey30", 
  Desulfobacterota="magenta",
  Firmicutes = "#3E9B96",
  Gammaproteobacteria = "greenyellow",
  Patescibacteria = "olivedrab",
  Myxococcota = "#B5D6AA",
  Nitrospirota = "palevioletred1",
  Proteobacteria = "royalblue",
  Planctomycetota = "darkorange", 
  "SAR324 clade(Marine group B)" = "olivedrab",
  Gemmatimonadota = "greenyellow",
  Thermoplasmatota = "green",
  Verrucomicrobiota = "darkorchid1",
  Other = "grey50",
  Unidentified = "grey90")

rare_phylum_colors <- c(
  Acidobacteriota = "navy", 
  Actinobacteriota = "darkslategray2", 
  Dependentiae = "royalblue",
  Alphaproteobacteria = "plum2", 
  Bacteroidota = "gold", 
  Bdellovibrionota = "red1",
  Chloroflexi="black", 
  Cyanobacteria = "limegreen",
  Spirochaetota = "grey30", 
  Desulfobacterota="magenta",
  Firmicutes = "#3E9B96",
  Gammaproteobacteria = "greenyellow",
  Patescibacteria = "olivedrab",
  Myxococcota = "#B5D6AA",
  Nitrospirota = "palevioletred1",
  Proteobacteria = "royalblue",
  Planctomycetota = "darkorange", 
  "SAR324 clade(Marine group B)" = "olivedrab",
  Elusimicrobiota = "greenyellow",
  Nanoarchaeota = "green",
  Verrucomicrobiota = "darkorchid1",
  Campylobacterota = "#6A3D9A",
  Gemmatimonadota = "steelblue4",
  Rare = "grey50",
  Unidentified = "white")

# Taken from:
# https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)


genus_colors <- c(
  "Rare" = "grey90",
  "OM27 clade" = "gold1",
  "betI-A" = "green4",
  "betI-B" = "green1",
  "betIV-A" = "palegreen2",
  "Candidatus Nitrotoga" = "#FB9A99",
  "acIV-C" = "blue1",
  "acIV-A" = "steelblue4",
  "acI-A" = "dodgerblue2",
  "acI-B" = "skyblue2",
  "hgcI clade" = "darkturquoise",
  "acI-C" = "turquoise1",
  "Cyanobium PCC-6307" = "olivedrab1",
  "bacIII-A" = "#FF7F00",
  "bacIII-B" = "#FDBF6F",
  "Fluviicola" = "khaki2",
  "bacI-A" = "yellow4",
  "bacII-A" = "yellow3",
  "Lacunisphaera" = "#CAB2D6",
  "LD19" = "#6A3D9A",
  "DEV114" = "orchid1",
  "SH3-11" = "deeppink1",
  "Chthoniobacter" = "maroon",
  "CL500-3" = "tan4",
  "Nitrospira" = "black",
  "alfV-A" = "#E31A1C"
)

genus_colors_2 <- c(
  "Rare" = "grey90",
  "acSTL-A" = "gold1",
  "betI-A" = "green4",
  "betI-B" = "green1",
  "betIV-A" = "palegreen2",
  "Hydrogenophaga" = "#FB9A99",
  "acIV-C" = "blue1",
  "acIV-A" = "steelblue4",
  "acI-A" = "dodgerblue2",
  "acI-B" = "skyblue2",
  "hgcI clade" = "darkturquoise",
  "acI-C" = "turquoise1",
  "Cyanobium PCC-6307" = "olivedrab1",
  "bacIII-A" = "#FF7F00",
  "bacIII-B" = "#FDBF6F",
  "Fluviicola" = "khaki2",
  "bacI-A" = "yellow4",
  "bacII-A" = "yellow3",
  "Flavobacterium" = "#CAB2D6",
  "LD19" = "#6A3D9A",
  "DEV114" = "orchid1",
  "SH3-11" = "deeppink1",
  "Luna1-A" = "maroon",
  "CL500-3" = "tan4",
  "Pnec" = "black",
  "alfV-A" = "#E31A1C"
)

class_colors <- c("Rare" = "grey80",
                  "Alphaproteobacteria" = "green4",
                  "Gammaproteobacteria"= "gold1",
                  "Bacteroidia" = "#FF7F00",
                  "Kapabacteria"= "#FDBF6F",
                  "Actinobacteria"      = "turquoise1",
                  "Verrucomicrobiae" = "orchid1",
                  "Chlamydiae"= "deeppink1",
                  "Cyanobacteriia"     = "olivedrab1",
                  "Acidimicrobiia" = "blue1",
                  "Vicinamibacteria" = "blue4",
                  "Holophagae"=  "steelblue2",
                  "Anaerolineae"        = "black",
                  "SL56 marine group"   = "grey20",
                  "JG30-KF-CM66"= "grey40",
                  "Chloroflexia"= "grey60",
                  "Phycisphaerae" = "tan4",
                  "Planctomycetes"= "tan3",
                  "OM190"= "tan",
                  "Nitrospiria" = "#6A3D9A",
                  "Bdellovibrionia"= "yellow4",
                  "Oligoflexia"=  "steelblue4",
                  "Armatimonadia"= "red",
                  "Gemmatimonadetes"= "cornflowerblue",
                  "Nitrososphaeria"= "maroon",
                  "Polyangia"  = "brown"
)

ggpreview <- function(...) {
  fname <- tempfile(fileext = ".png")
  ggsave(filename = fname, ...)
  system2("open", fname)
  invisible(NULL)
}

archaea_class_colors <- 
  c(
    "Methanobacteria" = "#1A3D82",
    "Methanosarcinia" = "#4499F5",
    "Aenigmarchaeia" = "#F0AC7D",
    "Nanoarchaeia" = "#832B0F",
    "Iainarchaeia" = "#426737",
    "Nitrososphaeria" = "#864F8F",
    "Rare" = "grey80",
    "Thermoplasmata" = "green4",
    "Methanomicrobia"= "steelblue4",
    "Micrarchaeia" = "#FF7F00",
    "Bathyarchaeia" = "gold1"
  )


process_colors <- 
  c(
    "DL" = "#1A3D82",
    "DR/HoS" = "#4499F5",
    "DL/HeS" = "#F0AC7D",
    "DR" = "#832B0F",
    "HoS" = "#426737",
    "DL/DR" = "#864F8F"
  )