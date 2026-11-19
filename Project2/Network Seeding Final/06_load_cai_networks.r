
# ============================================================
# 06_load_cai_networks.r
# Load and construct Cai et al. village networks
# ============================================================

library(tidyverse)
library(igraph)
library(haven)


# ------------------------------------------------------------
# 1. Load the original Cai network data
# ------------------------------------------------------------

cai_all <- read_dta("cai_data/0422allinforawnet.dta")


# ------------------------------------------------------------
# 2. Identify villages
# ------------------------------------------------------------

village_names <- unique(cai_all$address)

village_names <- village_names[
  !is.na(village_names) & village_names != ""
]


# ------------------------------------------------------------
# 3. Keep only valid within-village connections
# ------------------------------------------------------------

cai_edges <- cai_all %>%
  filter(
    address == network_address,
    !is.na(address),
    address != "",
    !is.na(network_id),
    id != 99,
    network_id != 99
  )


# ------------------------------------------------------------
# 4. Count retained edge records per village
# ------------------------------------------------------------

edge_counts <- cai_edges %>%
  count(address, name = "n_edges")


# ------------------------------------------------------------
# 5. Retain villages with more than 25 edge records
# ------------------------------------------------------------

village_names <- village_names[
  village_names %in%
    edge_counts$address[edge_counts$n_edges > 25]
]


# ------------------------------------------------------------
# 6. Construct a directed graph for each village
# ------------------------------------------------------------

networks_cai <- lapply(village_names, function(v) {

  edges_v <- cai_edges %>%
    filter(address == v) %>%
    select(id, network_id)

  edge_matrix <- as.matrix(
    data.frame(
      from = as.character(edges_v$id),
      to = as.character(edges_v$network_id)
    )
  )

  graph_from_edgelist(
    edge_matrix,
    directed = TRUE
  )
})


# ------------------------------------------------------------
# 7. Assign village names to the network list
# ------------------------------------------------------------

names(networks_cai) <- village_names

