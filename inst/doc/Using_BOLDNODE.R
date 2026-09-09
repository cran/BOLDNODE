## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE,
  warning = FALSE,
  message = FALSE
)

## ----install packages,message=F,warning=F,echo=FALSE,include=F----------------
# CRAN packages
library(BOLDNODE)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)

## ----package_install_for_vignette,message=F,warning=F,echo=FALSE--------------
# library(BOLDNODE)

## ----setup-parquet_toshow-----------------------------------------------------
# Path to the downloaded BOLD Parquet release file
# parquet_file <- "G:/usr/path/file.parquet"

## ----setup-parquet,include=F--------------------------------------------------
# Path to the downloaded BOLD Parquet release file

parquet_file <- system.file(
  "extdata",
  "test_data.parquet",
  package = "BOLDNODE"
)

## ----search,message=F,warning=F-----------------------------------------------
# Search the BOLD dataset
cerambycidae_search <- bold_parquet_search(
  input.parquet = parquet_file,
  taxonomy = "Cerambycidae",
  geography = "Canada",
  marker = "COI-5P"
)

## ----search-records,include=F-------------------------------------------------
tot_records <- cerambycidae_search %>%
  dplyr::summarise(Total_records = dplyr::n()) %>%
  dplyr::collect() %>%
  dplyr::pull(Total_records)

## ----print_records,echo=FALSE-------------------------------------------------
cat("The search has", tot_records, "records in the dataset")

## ----collect,warning=F,message=F----------------------------------------------
# Collect the search results into memory
cerambycidae_data <- bold_search_collect(
  cerambycidae_search,
  chunk.size = 50000,
  export = FALSE
)

## ----display-tabular_output,message=F,warning=F-------------------------------
# Inspect the collected data
DT::datatable(
  head(cerambycidae_data, 100) %>% select(processid, sampleid, bin_uri, family, genus, species),
  options = list(pageLength = 10, scrollX = TRUE)
)

## ----summary,message=F,warning=F----------------------------------------------
# Generate concise summary
cerambycidae_summary <- get_concise_summary(cerambycidae_search)
DT::datatable(cerambycidae_summary)

## ----step4-bin-consensus,message=F,warning=F----------------------------------
# Get BIN consensus with strict consensus threshold (1.0)
bin_consensus <- get_bin_consensus(
  cerambycidae_search,
  threshold = 1.0,
  min.ids = 1
)

# View consensus at different taxonomic ranks
table(bin_consensus$concordant_rank)

# Examine the results
DT::datatable(head(bin_consensus, 50),
  options = list(
    pageLength = 10,
    scrollX = TRUE
  )
)

## ----step5-occurrence-matrix,message=F,warning=F------------------------------
# Generate occurrence matrix at species level, grouped by province/state
occ_matrix <- bcdm_to_occmatrix(
  cerambycidae_search,
  kingdom = "Animalia",
  taxon.rank = "species",
  site.cat = "province.state"
)

# Visualize the output
DT::datatable(
  head(occ_matrix, 20),
  options = list(
    pageLength = 10,
    scrollX = TRUE
  )
)

## ----vegan-example-richness,message=F,warning=F-------------------------------
library(vegan)
species_richness_est <- poolaccum(occ_matrix)
species_richness_est

## ----column-to-rowname, include=F,message=F,warning=F-------------------------
occ_matrix <- occ_matrix %>%
  column_to_rownames("province.state")

## ----vegan-example-beta-diversity,message=F,warning=F-------------------------
# Calculate pairwise beta diversity
beta_diversity <- vegdist(occ_matrix, method = "bray") |> round(2)

# Visualize
DT::datatable(
  data.frame(as.matrix(beta_diversity)),
  options = list(
    pageLength = 10,
    scrollX = TRUE
  )
)

## ----nMDS,,message=F,warning=F,fig.width=12, fig.height=12, dpi=300-----------
library(ggrepel)
# Using the occurrence matrix generated above for nMDS
nmds <- metaMDS(
  occ_matrix,
  distance = "bray",
  k = 2,
  trymax = 100,
  trace = FALSE
)

# Site scores
sites <- as.data.frame(scores(nmds, display = "sites"))
sites$Site <- rownames(sites)

# Species scores
species <- as.data.frame(scores(nmds, display = "species"))
species$Species <- rownames(species)

# Keep only the 40 species furthest from the origin
species$dist <- sqrt(species$NMDS1^2 + species$NMDS2^2)
species20 <- species[order(species$dist, decreasing = TRUE), ][1:40, ]

# Plot
ggplot() +
  geom_point(
    data = sites,
    aes(NMDS1, NMDS2),
    size = 4
  ) +
  geom_text(
    data = sites,
    aes(NMDS1, NMDS2, label = Site),
    size = 4,
    alpha = 0.6,
    vjust = -0.4
  ) +
  geom_text_repel(
    data = species20,
    aes(NMDS1, NMDS2, label = Species),
    colour = "red",
    fontface = "bold",
    size = 3,
    max.overlaps = Inf
  ) +
  coord_equal() +
  theme_bw() +
  labs(
    title = "nMDS (Bray-Curtis)",
    x = "NMDS1",
    y = "NMDS2"
  )

## ----step5-bin-reps,message=F,warning=F---------------------------------------
# Sample 3 records per BIN-taxon combination
bin_reps <- get_bin_reps(
  bold.search.res = cerambycidae_search,
  Nreps = 3,
  by.tax = TRUE,
  enforce.scientific = TRUE,
  criteria = list(
    seq_length = 658,
    id_method = "Morphology",
    vouchered = TRUE
  )
)

# Compare sample to full dataset
n_bins <- length(unique(bin_reps$bin_uri))
cat(
  "Sampled", nrow(bin_reps), "representatives from", n_bins, "BINs",
  "(out of", tot_records, "total records)."
)

## ----step5-bin-rep-DT,message=F,warning=F-------------------------------------
DT::datatable(
  head(bin_reps, 50),
  options = list(
    pageLength = 10,
    scrollX = TRUE
  )
)

## ----DNAStringSet-object,message=F,warning=F----------------------------------
library(Biostrings)
# Filter for genus Clytus
clytus_search <- bold_parquet_search(
  input.parquet = parquet_file,
  taxonomy = "Clytus",
  geography = "Canada",
  marker = "COI-5P",
  basecount = 658,
  ambi.base.cutoff = "<1%"
)

# Collect the filtered data
clytus_data <- bcdm_to_dnastringset(clytus_search,
  cols_for_seq_names = c("processid", "bin_uri", "species")
)
head(clytus_data, 3)

## ----basepair-freq,message=F,warning=F----------------------------------------
library(Biostrings)
# Count bases across all sequences
base_counts <- colSums(alphabetFrequency(clytus_data, baseOnly = TRUE))
# Proportions
base_props <- base_counts / sum(base_counts)
base_props

## ----sequence_alignment,message=F,warning=F-----------------------------------
library(Biostrings)
library(muscle)
alignment_muscle <- muscle(clytus_data)

## ----NJ tree,message=F,warning=F,fig.width=12, fig.height=12, dpi=300---------
library(ape)
library(phangorn)
# Convert the alignment to DNABin
dna_bin <- as.DNAbin(alignment_muscle)

# Distance matrix using K80
dist_matrix <- dist.dna(
  dna_bin,
  model = "K80"
)

# Neighbor Joining tree
nj_tree <- nj(dist_matrix)

# midpoint rooting
nj_tree <- midpoint(nj_tree)

# Plot tree
plot(
  nj_tree,
  cex = 0.6,
  main = "Neighbor-Joining Tree"
)

## ----sf_object----------------------------------------------------------------
# Filter for genus Monochamus
monochamus_search <- bold_parquet_search(
  input.parquet = parquet_file,
  taxonomy = "Monochamus",
  geography = "Canada",
  marker = "COI-5P"
)

# Collect the filtered data
monochamus_data <- bcdm_to_sf(monochamus_search)

## ----occurrence-map,message=F,warning=F---------------------------------------
library(sf)
# Creating a background map using the maps package; some map_data country names (ID column) are changed to suit the BCDM country.ocean names
map_data <- st_as_sf(maps::map("world",
  plot = FALSE,
  fill = TRUE
)) %>%
  filter(ID == "Canada")
# Convert the data to WGS84
map_data <- st_transform(
  map_data,
  4326
)
# Plot
map_plot <- ggplot() +
  geom_sf(
    data = map_data,
    alpha = 0.3,
    linewidth = 0.4
  ) +
  geom_point(
    data = monochamus_data,
    mapping = aes(
      x = lon,
      y = lat
    ),
    colour = "#011B26",
    fill = "#F78E1E",
    size = 3,
    pch = 21
  ) +
  theme_bw(base_size = 15) +
  theme(panel.grid.major = element_line(
    colour = "grey50",
    size = 0.3,
    linetype = 3
  )) +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_sf(expand = FALSE) +
  ggtitle("Distribution map")

map_plot

