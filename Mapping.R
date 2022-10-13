library(onemap)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(gdata)
library(here)
​
​
# user-defined objects. Output file names are automatically attached
fam <- "CM7" # correctly specify the family
l <- 10 # chromosome on which genetic map is build
dir.in <- "RESULT/ONEMAP" # dir of the input file
file.in <- "F2_CM7_QCFI_final4.onemap.txt" # file name of the input file
map_results_dir <- "RESULT/FinalMaps" # folder to save genetic map
​
# make dir to save
#dir.create(map_results_dir, recursive = TRUE)
​
# import genotype data as onemap object
markers_onemap <- read_onemap(dir = dir.in, inputfile = file.in)
​
# all markers
n.mrk.all <- markers_onemap$n.mar
​
# twoopt analysis (by Dan -- I am not sure if this is critical...)
LOD_sug <- suggest_lod(markers_onemap)
twopts_f2 <- rf_2pts(markers_onemap, LOD = LOD_sug)
​
# extract pairwise recombination
pairwise_recombination <- twopts_f2$analysis %>%
  as_tibble() %>%
  mutate(marker1 = names(twopts_f2$analysis %>% as_tibble())) %>%
  select(marker1, everything()) %>%
  gather(marker2, recomb_rate, -marker1, convert = TRUE) %>%
  filter(marker1 != marker2)
​
# segregation test
seg_test <- test_segregation(markers_onemap)
png(paste0(map_results_dir, "/Segreg.Test_", fam, "_chr", l, ".png"))
plot(seg_test) # some markers are significantly distorted
dev.off()
​
# markers to be removed
mrk.name.all <- colnames(twopts_f2$analysis)
mrk.name.rm.dist <- select_segreg(seg_test, distorted = T, numbers = F)
mrk.rm.tf <- mrk.name.all %in% mrk.name.rm.dist
mrk.retain.num <- which(!mrk.rm.tf)
​
# make sequence from those markers
mark_no_dist_filtered_f2 <- make_seq(twopts_f2, mrk.retain.num)
​
# make a linkage group (Modified by Ryokei, 09/29/2021, given an error caused by weak linkage within a chromosome)
LGs_f2 <- group(mark_no_dist_filtered_f2) # NO LONGER USE THIS FUNCTION
snpname.retained <- colnames(markers_onemap$geno)[mark_no_dist_filtered_f2$seq.num]
myfun <- function(x){as.integer(gsub("chr", "", strsplit(x, "_")[[1]][1]))}
LGs_f2$groups <- sapply(snpname.retained, myfun) # OVERWRITE
LG_summary <- LGs_f2$groups %>%
  as_tibble() %>%
  count(value) %>%
  rename(LG = value) %>%
  rename(nr_markers = n)
​
# set map function
set_map_fun(type = "kosambi")
​
# make seq
LG_f2 <- make_seq(LGs_f2, l)
​
# de novo ordered map
LG_f2_ord <- order_seq(input.seq = LG_f2, n.init = 5, subset.search = "twopt", twopt.alg = "rcd", THRES = 6, touchdown = TRUE)
lg_safe_map <- make_seq(LG_f2_ord, "safe")
lg_all_map <- make_seq(LG_f2_ord, "force")
​
# rearrange safe markers in B73 order & make map
markers_name_num_pairs <- markers_onemap$geno %>%
  as_tibble() %>%
  names() %>%
  as_tibble() %>%
  rowid_to_column() %>%
  rename(numeric_id = rowid, marker = value)
safe_known_order <- markers_name_num_pairs %>%
  filter(numeric_id %in% lg_safe_map$seq.num) %>%
  separate(marker, c("chr", "pos"), sep = "_", remove = FALSE, convert = TRUE) %>%
  arrange(chr, pos) %>%
  select(numeric_id) %>%
  unlist() %>%
  as_vector()
all_known_order <- markers_name_num_pairs %>%
  filter(numeric_id %in% lg_all_map$seq.num) %>%
  separate(marker, c("chr", "pos"), sep = "_", remove = FALSE, convert = TRUE) %>%
  arrange(chr, pos) %>%
  select(numeric_id) %>%
  unlist() %>%
  as_vector()
lg_b73ord_safe_map <- make_seq(twopts_f2, safe_known_order) %>% onemap::map()
lg_b73ord_all_map <- make_seq(twopts_f2, all_known_order) %>% onemap::map()
​
# write map
write_map(lg_safe_map, here(map_results_dir, paste0("x", fam, "_chr", l, "_safe_map3.txt")))
write_map(lg_all_map, here(map_results_dir, paste0("x", fam, "_chr", l, "_all_map3.txt")))
write_map(lg_b73ord_safe_map, here(map_results_dir, paste0("x", fam, "_chr", l, "_safe_b73ord_map3.txt")))
write_map(lg_b73ord_all_map, here(map_results_dir, paste0("x", fam, "_chr", l, "_all_b73ord_map3.txt")))




