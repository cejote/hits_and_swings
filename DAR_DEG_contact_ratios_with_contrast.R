#
# Script to decide on how to merge DAR/DEG windows
# Produces randomized background data for enrichment tests
#

options(readr.num_threads = 1)
options(scipen = 99999)
options(stringsAsFactors = FALSE)
options(digits = 7)

library(tidyverse)
library(ggrepel)
library(ggbeeswarm)
library(patchwork)

outdir = "DARDEG"

###################
# Get DAR/DEG regions before overlapping with HITS
dars = read_tsv("./data/DARs_uq_windows.bed") %>% arrange(chrom, start)
degs = read_tsv("/data/DEGS_uq_windows.bed") %>% arrange(chrom, start)

chrOrder = paste0("chr", c(1:19))

dars = dars %>% filter(chrom %in% chrOrder)
degs = degs %>% filter(chrom %in% chrOrder)

###################
# Determine merge distance effect by incrementing window join size
regioncounts = c()
for(dmax in 1:10) {
  BINS = 30000 * dmax
  
  dardeg_windows = bind_rows(dars, degs) %>%
    distinct(chrom, start, end) %>%
    arrange(chrom, start) %>%
    mutate(
      grp = cumsum(start - lag(start, default = 0) > BINS),
      cgrp = cumsum(start - lag(start, default = 0)),
      dist = start - lag(start)
    ) %>%
    group_by(chrom, grp) %>%
    mutate(
      grpstart = min(start),
      grpend = max(start),
      grplength = max(start) - min(start),
      grpstr = sprintf("%s:%s-%s", chrom, grpstart, grpend)
    ) %>%
    arrange(chrom, grp) %>%
    distinct(grpstr)
  
  regioncounts = regioncounts %>% bind_rows(data.frame(dmax = dmax, regions = nrow(dardeg_windows)))
}

ggplot(regioncounts %>% mutate(loss = lag(regions, default = 0) - regions), aes(x = as.factor(dmax), y = abs(loss))) +
  geom_point() +
  scale_y_log10() +
  ggtitle("Number of regions lost joining [dmax] windows") +
  theme_bw()

bind_rows(dars %>% mutate(g = "dar"), degs %>% mutate(g = "deg")) %>%
  filter(chrom %in% paste0("chr", 1:6)) %>%
  ggplot(aes(xmin = start, xmax = pmax(end, start + 50000), ymin = 0, ymax = 1, fill = g)) +
  geom_rect() +
  facet_grid(chrom ~ .) +
  ggtitle("Positions of DAR and DEG loci (chroms 1-6)") +
  theme_bw()

###################
# Create randomized regions and define merged DAR/DEG blocks
BINS = 30000

dardeg_windows = bind_rows(dars, degs) %>%
  distinct(chrom, start, end) %>%
  group_by(chrom) %>%
  arrange(start) %>%
  mutate(
    grp = cumsum(start - lag(start, default = 0) > BINS),
    dist = start - lag(start)
  ) %>%
  group_by(chrom, grp) %>%
  mutate(
    grpstart = min(start),
    grpend = max(start),
    grplength = max(start) - min(start) + 30000,
    grpstr = sprintf("%s:%s-%s", chrom, grpstart, grpend)
  ) %>%
  ungroup() %>%
  arrange(chrom, grp)

# Inspect grouping
dardeg_windows %>% distinct(chrom, grp) %>% count()
dardeg_windows %>% distinct(grpstr) %>% count()
dars %>% count()
dars %>% count(chrom)
degs %>% count()
degs %>% count(chrom)
dardeg_windows %>% group_by(chrom) %>% summarize(gpm = max(grp))
dardeg_windows %>% distinct(chrom, grp)

# Plot distribution of distances between unmerged regions
dardeg_windows.dists = dardeg_windows %>% filter(!is.na(dist), dist > 30000) %>% count(dist = log10(dist))
ggplot(dardeg_windows.dists, aes(x = 10^dist, y = n)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  ggtitle("Distance between DAR/DEG regions before merging") +
  theme_bw()

###################
# Randomize regions for background generation

g.background = read_tsv("genome_windows_30000.bed", col_names = c("chrom", "start", "end")) %>%
  group_by(chrom) %>%
  summarise(glen = plyr::round_any(max(end), 30000, floor))

g.background %>% count(chrom) %>% print(n = 22)

# Join genomic background size info
dardeg_windows = dardeg_windows %>% left_join(g.background, by = c("chrom"))

system(sprintf("mkdir -p ./results/%s/randomregions/", outdir))

numrand = 100
rands = unique(plyr::round_any(runif(400, 30000 * 10, 30000 * 10000), 30000))[1:numrand]
shift = c(0, sort(rands))

i = 0
for(k in shift) {
  i = i + 1
  print(sprintf("%s: %s", i, k))
  
  shifted_windows = dardeg_windows %>%
    mutate(
      start.adj = (start + k) %% glen,
      end.adj = (end + k) %% glen,
      grpstart.shifted = grpstart + k,
      grpend.shifted = grpend + k,
      g2 = grpend.shifted > glen,
      g1 = grpstart.shifted < glen & g2
    ) %>%
    filter(!g1)  # Exclude regions that wrap around chromosome ends
  
  shifted_windows %>%
    write_tsv(sprintf("./results/%s/randomregions/dardegs_%d.tsv", outdir, k))
  
  shifted_windows %>%
    select(chrom, start = start.adj, end = end.adj, grp, grplength) %>%
    write_tsv(sprintf("./results/%s/randomregions/dardegs_%d.bed", outdir, k))
}

length(shift)
sum(duplicated(shift))
