# Progetto Tesi: Importance of Protected Areas in Italy
# This analysis aims to ascertain how many IUCN Red List (Italian) species are contained within Italy's protected areas (Euap Sites and Network-Nature 00 Sites).

# Install packages ----
# install.packages("rgbif")
# install.packages("mapview")
# install.packages("CoordinateCleaner")
# install.packages("countrycode")
# install.packages("sp")
# install.packages("maps")
# install.packages("TNRS")
# install.packages("sf")
# install.packages("styler")

# Load packages ----
library(rgbif)
library(mapview)
library(countrycode)
library(CoordinateCleaner)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sp)
library(maps)
#mapdata
library(TNRS)
library(sf)
library(styler)

# Get Data from Gbif ----
gbif_occurrences <- occ_download(
  pred_in("basisOfRecord", c("OBSERVATION", "HUMAN_OBSERVATION")),
  pred("country", "IT"),
  pred("taxonKey", 7707728),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred_gte("year", 1985),
  pred_gte("coordinateUncertaintyInMeters", 0),
  pred_lte("coordinateUncertaintyInMeters", 5000),
  format = "SIMPLE_CSV",
  user = "user", pwd = "pwdw", email = "email"
)

gbif_status <- occ_download_wait(gbif_occurrences[1])

gbif_it <- occ_download_get(gbif_occurrences[1]) |>
  occ_download_import()

write.csv(gbif_it, "gbif_tracheofite.csv", row.names = FALSE)
gbif_it <- read.csv("gbif_tracheofite.csv")

# Select columns of interest
gbif_it <- gbif_it |>
  dplyr::select(
    species, decimalLongitude, decimalLatitude, countryCode, individualCount,
    gbifID, family, taxonRank, coordinateUncertaintyInMeters, year,
    basisOfRecord, institutionCode, scientificName, stateProvince
  )

# Standardization taxonomy
gbif_it <- dplyr::filter(gbif_it, species != "")

names_species <- unique(gbif_it$species)
is.vector(names_species)

names_correct <- TNRS(names_species, sources = "wcvp")
?TNRS
TNRS_sources()

# accepted_name

names_correct <- names_correct |>
  dplyr::select(Name_submitted, Accepted_name)

# write.csv(names_correct, "dateset_nomi_corretti_Gbif.csv")
# names_correct <- read.csv("dateset_nomi_corretti_Gbif.csv", sep = ",")

colnames(names_correct)[1] <- "species"

gbif_it <- merge(
  x = gbif_it, y = names_correct,
  by = "species", all.x = TRUE
)
gbif_it <- gbif_it[, -1]

gbif_it <- gbif_it |>
  dplyr::select(Accepted_name, everything())


# Red List ----
red_list<- read.csv("assessments.csv", sep=",")
names(red_list)
summary(red_list)

# Taxonomy standardisation 
names_species <- unique(red_list$scientificName)
names_correct <- TNRS(names_species)
names_confront <- names_correct %>%
  select(Name_submitted, Name_matched)
write.csv(names_confront, "nomi_stand_RL.csv")

colnames(names_confront)[1] <- "scientificName"

red_list <- merge(x=red_list,y=names_confront, 
                  by="scientificName", all.x=TRUE)
red_list <- red_list[,-1]

red_list <- red_list %>%  
  select(Name_matched, everything())  
colnames(red_list)[1] <- "Name_species"
summary(red_list)

# Rimozione di alcune colonne non rilevanti
# red_list_pulito <- red_list %>%
# select(-language, -yearLastSeen)
write.csv(red_list, "red_list_finale.csv")

red_list <- read.csv("red_list_finale.csv")

# Check which species from the GBIF dataset are present in the Red List
colnames(gbif_it)[1] <- "Name_species"
species_gbif_rl <- dplyr::left_join(gbif_it, red_list, by = "Name_species")
species_gbif_rl <- species_gbif_rl[!is.na(species_gbif_rl$redlistCategory), ]
class(species_gbif_rl)
summary(species_gbif_rl)
sp_gbif_rl <- unique(species_gbif_rl$Name_species)

## Gbif Species Catalog in the Red List ----
lista_sp_cat_rischio <- split(species_gbif_rl, species_gbif_rl$redlistCategory)

names(lista_sp_cat_rischio)

sp_CE <- lista_sp_cat_rischio[[1]]
sp_DD <- lista_sp_cat_rischio[[2]]
sp_E <- lista_sp_cat_rischio[[3]]
sp_LC <- lista_sp_cat_rischio[[4]]
sp_NT <- lista_sp_cat_rischio[[5]]
sp_V <- lista_sp_cat_rischio[[6]]

unique(sp_CE$Name_species) # Count the number of species for each Risk Category

#rm(sp_CE,sp_DD,sp_E,sp_LC,sp_NT,sp_V)

# Transform 'gbif_it' into an 'sf' (spatial) object
species_gbif_rl <- st_as_sf(species_gbif_rl, coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +datum=WGS84")

# Display the occurrences in 'gbif_rl'
ita_map <- map_data("italy")
ggplot() +
  # coord_fixed() +
  geom_polygon(data = ita_map, aes(x = long, y = lat, group = group), colour = "gray50", fill = "gray50", alpha = 0.5) +
  geom_sf(data = species_gbif_rl, colour = "darkred", size = 0.01, alpha = 0.5) +
  theme_minimal()
dev.off()


# DIFFERENZE SPECIE_GBIF_RL FIRST AND LAST DATASET ----
species_gbif_rl_first <- st_read("species_gbif_rl.shp")
species_gbif_rl_first <- species_gbif_rl_first |> 
  filter(year >= 1985 | is.na(year)) # seleziono solo le specie dal 1985 ad OGGI ed elimino gli NA 

species_gbif_rl_no_geom <- st_drop_geometry(species_gbif_rl)

## specie gbif_rl che sono presenti solo nel dataset FIRST ----
colnames(species_gbif_rl_first)[1] <- "Name_species"
species_gbif_rl_miss <-  anti_join(species_gbif_rl_first, species_gbif_rl_no_geom, by = "Name_species")
unique(species_gbif_rl_miss$Name_species) #168 specie perse dal 1° dataset...

species_gbif_rl_miss_cat_risk <- species_gbif_rl_miss |> 
  group_by(rdlstCt) |>
  summarise(
    num_species = n_distinct(Name_species)
  ) |> 
  st_drop_geometry()

species_gbif_rl_miss_basis_of_Rc <- species_gbif_rl_miss |> 
  group_by(bssOfRc) |>
  summarise(
    num_species = n_distinct(Name_species)
  ) |> 
  st_drop_geometry()

species_gbif_rl_miss_ho <- species_gbif_rl_miss |> 
  filter(bssOfRc == "HUMAN_OBSERVATION" | bssOfRc == "OCCURRENCE")

unique(species_gbif_rl_miss_ho$Name_species)

# VI Elenco Ufficiale Aree Naturali Protette (EUAP) ----
# Load the .shp file for protected areas
siti_protet <- st_read("C:/Project_tirocinio/data/Siti_protetti_EUAP.shp")

# Fix the geometry of the protected areas' polygons
siti_protet <- st_make_valid(siti_protet)

# st_crs(species_gbif_rl) == st_crs(siti_protet) # Verify the coordinate reference systems (CRS)

unique(siti_protet$tipo)
class(siti_protet$tipo)

# Remove Marine Protected Areas and Underwater Natural Parks
siti_protet <- siti_protet |>
  filter(tipo != "MAR" & tipo != "GAPN")

# Convert the 'tipo' column to an ordered factor
livelli_tipo <- c("RNS", "PNZ", "PNR", "RNR", "AANP")
siti_protet$tipo <- factor(siti_protet$tipo, levels = livelli_tipo)

# mapview(siti_protet) # open street map

# Plot of protected areas categorized by type
ggplot() +
  coord_fixed() +
  geom_polygon(data = ita_map, aes(x = long, y = lat, group = group), colour = "gray50", fill = "grey70", alpha = 0.5) +
  geom_sf(data = siti_protet, aes(fill = tipo), colour = NA) +
  scale_fill_manual(
    values = c("RNS" = "darkred", "PNZ" = "red", "PNR" = "orange", "RNR" = "yellow", "AANP" = "green"),
    name = "Tipologia Area Protetta",
    labels = c(
      "Riserve Naturali Statali (RNS)", "Parchi Naturali Nazionali (PNZ)",
      "Parchi Naturali Regionali (PNR)", "Riserve Naturali Regionali (RNR)",
      "Altre Aree Naturali Protette (AANP)"
    )
  ) +
  theme_minimal()
dev.off()

## Total number of areas and total area size for each type of protected area ----
# Convert the 'superficie' column to numeric
siti_protet$superficie <- as.numeric(siti_protet$superficie)

# Group by 'tipo', count the number of areas, and calculate total area size, then drop geometry
siti_protet_type <- siti_protet |>
  group_by(tipo) |>
  summarise(
    num = n(),
    superficie_tot = sum(superficie)
  ) |>
  st_drop_geometry()


# SPECIES COUNT ----
## Count occurrences WITHIN the protected areas ----

species_gbif_rl$index_sp <- 1:nrow(species_gbif_rl)
sp_in_aree_pr <- st_within(species_gbif_rl, siti_protet)
sp_in_aree_pr_df <- as.data.frame(sp_in_aree_pr)
sp_in_aree_pr_df <- subset(sp_in_aree_pr_df, !duplicated(row.id))

colnames(sp_in_aree_pr_df)[1] <- "index_sp"
sp_in_aree_pr_df$col.id <- siti_protet$nome_gazze[sp_in_aree_pr_df$col.id]
colnames(sp_in_aree_pr_df)[2] <- "protected_area"
# species_cat_risk$index_sp <- as.numeric(species_cat_risk$index_sp)
sp_in_aree_pr_df <- left_join(sp_in_aree_pr_df, species_gbif_rl, by = "index_sp")
sp_in_aree_pr_df <- sp_in_aree_pr_df[, -1]
#colnames(sp_in_aree_pr_df)[5] <- "occurence_gbifID"
sp_in_aree_pr_df <- sp_in_aree_pr_df[, c(5, 1:4, 6:ncol(sp_in_aree_pr_df))]
sp_in_aree_pr_df <- sp_in_aree_pr_df[, c(1, 3, 2, 4:ncol(sp_in_aree_pr_df))]
st_write(sp_in_aree_pr_df, "C:/Project_tirocinio/OUTPUT_TIROCINIO_2/sp_in_aree_pr_df.shp")


unique(sp_in_aree_pr_df$Name_species) # number of species

# Species divided by risk category, summing the number of species and occurrences
sp_in_aree_pr_cat_risk <- sp_in_aree_pr_df |>
  group_by(redlistCategory) |>
  summarise(
    num_species = n_distinct(Name_species),
    num_occorrenze = n()
  )

# Division by risk categories
list_sp_in_cat_risk <- split(sp_in_aree_pr_df, sp_in_aree_pr_df$rdlstCt)
sp_in_DD <- list_sp_in_cat_risk[[2]]
sp_in_LC <- list_sp_in_cat_risk[[4]]
sp_in_NT <- list_sp_in_cat_risk[[5]]
sp_in_V <- list_sp_in_cat_risk[[6]]
sp_in_E <- list_sp_in_cat_risk[[3]]
sp_in_CE <- list_sp_in_cat_risk[[1]]

#rm(sp_in_DD, sp_in_LC, sp_in_NT, sp_in_V, sp_in_E, sp_in_CE)

## Count occurrences OUTSIDE the protected areas ----
sp_out_aree_pr_df <- as.data.frame(sp_in_aree_pr)
colnames(sp_out_aree_pr_df)[1] <- "index_sp"
sp_out_aree_pr_df$col.id <- siti_protet$nome_gazze[sp_out_aree_pr_df$col.id]
colnames(sp_out_aree_pr_df)[2] <- "protected_area"
sp_out <- pull(sp_out_aree_pr_df, index_sp)
sp_out_aree_pr_df <- species_gbif_rl |>
  filter(!(index_sp %in% sp_out))
sp_out_aree_pr_df <- sp_out_aree_pr_df[, -ncol(sp_out_aree_pr_df)]
st_write(sp_out_aree_pr_df, "sp_out_aree_pr_df.shp")
# nrow(sp_in_aree_pr_df) + nrow(sp_out_aree_pr_df)  # check

# number of species
unique(sp_out_aree_pr_df$Name_species)

# Species categorized by risk category, summing the number of species and occurrences
sp_out_aree_pr_cat_risk <- sp_out_aree_pr_df |>
  group_by(redlistCategory) |>
  summarise(
    num_species = n_distinct(Name_species),
    num_occorrenze = n()
  ) |>
  st_drop_geometry()

# Categorization by risk categories
list_sp_out_cat_risk <- split(sp_out_aree_pr_df, sp_out_aree_pr_df$rdlstCt)
sp_out_DD <- list_sp_out_cat_risk[[2]]
sp_out_LC <- list_sp_out_cat_risk[[4]]
sp_out_NT <- list_sp_out_cat_risk[[5]]
sp_out_V <- list_sp_out_cat_risk[[6]]
sp_out_E <- list_sp_out_cat_risk[[3]]
sp_out_CE <- list_sp_out_cat_risk[[1]]


# plot delle occorrenze fuori dalle aree protette #nel caso lo metta devo costruirmi
# i vari sotto dataset per le cat di rischio
# ita_map <- map_data("italy")
# ggplot() +
#   coord_fixed() +
#   geom_polygon(data = ita_map, aes(x = long, y = lat, group = group), colour = "gray50", fill = "gray50", alpha = 0.5) +
#   geom_sf(data = siti_protet, color = "darkgreen", fill = "green1", alpha = 0.3) +
#   geom_sf(data = sp_CE_out_areepr, aes(color = "Critically Endangered (CE)"), size = 0.04) +
#   geom_sf(data = sp_E_out_areepr, aes(color = "Endangered (EN)"), size = 0.04) +
#   geom_sf(data = sp_V_out_areepr, aes(color = "Vulnerable (VU)"), size = 0.04) +
#   geom_sf(data = sp_NT_out_areepr, aes(color = "Near Threatened (NT)"), size = 0.04) +
#   scale_color_manual(
#     name = "RedList Category",
#     values = c("darkred", "red", "orange2", "yellow"),
#     labels = c(
#       "Critically Endangered (CE)", "Endangered (EN)", "Vulnerable(VU)",
#       "Near Threatened (NT)"
#     )
#   ) +
#   theme_minimal()
# dev.off()


# per un po' di pulizia:
# rm(sp_in_aree_pr, sp_in_aree_pr_cat_risk, sp_in_aree_pr_df)
# rm(sp_out_aree_pr_cat_risk, sp_out_aree_pr_df)


## SHARED SPECIES ----
sp_condivise_out <- semi_join(sp_out_aree_pr_df, sp_in_aree_pr_df, by = "Name_species")
sp_condivise_in <- semi_join(sp_in_aree_pr_df, sp_out_aree_pr_df, by = "Name_species")
# or
# sp_condivise <- inner_join(sp_out_aree_pr_df, sp_in_aree_pr_df, by = "Nm_spcs")
nrow(sp_condivise_in) + nrow(sp_condivise_out)

length(unique(sp_condivise_in$Name_species))
length(unique(sp_condivise_out$Name_species))
# length(unique(species_gbif_rl$Nm_spcs))

# Group shared species ("inside and outside") generically
sp_condivise_cat_risk <- sp_condivise_in |> # not change if I use sp_condivise_out
  group_by(redlistCategory) |>
  summarise(num_species = n_distinct(Name_species)) |>
  st_drop_geometry()

## Count the species that are ONLY WITHIN the protected areas----
df_sp_only_in <- anti_join(sp_in_aree_pr_df, sp_out_aree_pr_df, by = "Name_species")
unique(df_sp_only_in$Name_species)

# Remove occurrences with risk level DD (Data Deficient) or LC (Least Concern)
# df_sp_only_in <- df_sp_only_in[!(df_sp_only_out$rdlstCt %in% c("Least Concern", "Data Deficient")), ]

# write.csv(df_sp_only_in,"df_sp_only_in.csv")
df_sp_only_in <- read.csv("df_sp_only_in.csv")

# Group species and occurrences by risk category
sp_only_in_cat_risk <- df_sp_only_in |>
  group_by(redlistCategory) |>
  summarise(
    num_species = n_distinct(Name_species),
    num_occorrenze = n()
  )

# Division by risk categories
list_sp_only_in_cat_risk <- split(df_sp_only_in, df_sp_only_in$rdlstCt)
sp_only_in_DD <- list_sp_only_in_cat_risk[[2]]
sp_only_in_LC <- list_sp_only_in_cat_risk[[4]]
sp_only_in_NT <- list_sp_only_in_cat_risk[[5]]
sp_only_in_V <- list_sp_only_in_cat_risk[[6]]
sp_only_in_E <- list_sp_only_in_cat_risk[[3]]
sp_only_in_CE <- list_sp_only_in_cat_risk[[1]]

unique(sp_only_in_CE$Name_species) # to count the number of species for each category


# some check
nrow(sp_in_aree_pr_df) == nrow(sp_condivise_in) + nrow(df_sp_only_in)
length(unique(sp_in_aree_pr_df$Name_species)) == length(unique(sp_condivise_in$Name_species)) + length(unique(df_sp_only_in$Name_species))


## Count the species that are ONLY OUTSIDE the protected areas ----
df_sp_only_out <- anti_join(sp_out_aree_pr_df, sp_in_aree_pr_df, by = "Name_species")
unique(df_sp_only_out$Name_species)

# Remove occurrences with risk level DD (Data Deficient) or LC (Least Concern)
# df_sp_only_out <- df_sp_only_out[!(df_sp_only_out$rdlstCt %in% c("Least Concern", "Data Deficient")), ]

#write.csv(df_sp_only_out, "df_sp_only_out.csv")
df_sp_only_out <- read.csv("df_sp_only_out.csv")

# Group species and occurrences by risk category
sp_only_out_cat_risk <- df_sp_only_out |>
  group_by(redlistCategory) |>
  summarise(
    num_species = n_distinct(Name_species),
    num_occorrenze = n()
  ) |>
  st_drop_geometry()

# Categorization by risk categories
list_sp_only_out_cat_risk <- split(df_sp_only_out, df_sp_only_out$rdlstCt)
# sp_only_out_DD
# sp_only_out_LC
sp_only_out_NT <- list_sp_only_out_cat_risk[[5]]
sp_only_out_V <- list_sp_only_out_cat_risk[[6]]
sp_only_out_E <- list_sp_only_out_cat_risk[[3]]
sp_only_out_CE <- list_sp_only_out_cat_risk[[1]]

unique(sp_only_out_CE$Nm_spcs)

# some checks
nrow(sp_out_aree_pr_df) == nrow(sp_condivise_out) + nrow(df_sp_only_out)
length(unique(sp_out_aree_pr_df$Name_species)) == length(unique(sp_condivise_out$Name_species)) + length(unique(df_sp_only_out$Name_species))


#### PLOT of occurrences belonging to species that are ONLY OUTSIDE the protected areas ----
# Exclude species with risk levels DD (Data Deficient) and LC (Least Concern) as they are not very relevant for this analysis
df_sp_only_out_most_risk <- df_sp_only_out |>
  filter(!(redlistCategory %in% c("Data Deficient", "Least Concern")))

unique(df_sp_only_out_most_risk$Name_species)

ita_map <- map_data("italy")
df_sp_only_out_most_risk$redlistCategory <- factor(df_sp_only_out_most_risk$redlistCategory,
  levels = c(
    "Critically Endangered", "Endangered",
    "Vulnerable", "Near Threatened"
  )
)
ggplot() +
  # coord_fixed() +
  geom_polygon(data = ita_map, aes(x = long, y = lat, group = group), colour = "gray50", fill = "gray70", alpha = 0.5) +
  geom_sf(data = siti_protet, color = "darkgreen", fill = "green3", alpha = 0.3, size = 0.1) +
  geom_sf(data = df_sp_only_out_most_risk, aes(color = redlistCategory), size = 0.6) +
  scale_color_manual(
    name = "RedList Category",
    values = c(
      "Critically Endangered" = "darkred", "Endangered" = "red",
      "Vulnerable" = "orange2", "Near Threatened" = "yellow"
    ),
    labels = c(
      "Critically Endangered (CE)", "Endangered (EN)", "Vulnerable(VU)",
      "Near Threatened (NT)"
    ),
    guide = guide_legend(override.aes = list(size = 1)) # cambiare dimensione legenda
  ) +
  theme_minimal()
dev.off()

# per pulire tutto l'environment
rm(list = ls())

# capire se serve o meno ----
## QUALI SP SONO CONTENUTE IN MENO AREE PROTETTE? ----
### Sp DENTRO ----
sp_in_low_protect <- sp_in_aree_pr_df %>%
  group_by(Nm_spcs) %>%
  summarise(
    num_areepr = n_distinct(protected_area),
    num_occorrenze = n()
  ) %>%
  arrange(num_areepr)

### Sp SOLO DENTRO ----
sp_only_in_low_protect <- df_sp_only_in %>%
  group_by(Nm_spcs) %>%
  summarise(
    num_areepr = n_distinct(protected_area),
    num_occorrenze = n()
  ) %>%
  arrange(num_areepr)

rm(names_correct, sp_condivise_cat_risk, sp_in_aree_pr_cat_risk,
   sp_only_in_cat_risk, sp_only_out_cat_risk, sp_out_aree_pr_cat_risk)

# ANALYSIS FOR EACH TYPE OF PROTECTED AREA ----
lista_aree_pr <- split(siti_protet, siti_protet$tipo)

lista_aree_pr[[2]] <- lista_aree_pr[[2]] |>
  filter(nome_gazze != "Parco Nazionale dell'Arcipelago di La Maddalena")
pnz <- lista_aree_pr[[2]]

species_gbif_rl$index_sp <- 1:nrow(species_gbif_rl)

# Create a function to analyze the types of protected areas

function_type_areepr <- function(tipo_areepr) {
  sp_in <- st_within(species_gbif_rl, tipo_areepr)
  sp_in <- as.data.frame(sp_in)
  colnames(sp_in)[1] <- "index_sp"
  colnames(sp_in)[2] <- "protected_area"
  sp_in$protected_area <- tipo_areepr$nome_gazze[sp_in$protected_area]
  sp_in <- left_join(sp_in, species_gbif_rl, by = "index_sp")
  sp_in <- sp_in[, -1]
  colnames(sp_in)[5] <- "occurence_gbifID"
  sp_in <- sp_in[, c(5, 1:4, 6:ncol(sp_in))]
  sp_in <- sp_in[, c(1, 3, 2, 4:ncol(sp_in))]
  return(sp_in)
}

results_typearee_sp_in <- lapply(lista_aree_pr, function_type_areepr)

# Dataset containing occurrences and species for each type of protected area
sp_in_rns <- results_typearee_sp_in[[1]]
sp_in_pnz <- results_typearee_sp_in[[2]]
sp_in_pnr <- results_typearee_sp_in[[3]]
sp_in_rnr <- results_typearee_sp_in[[4]]
sp_in_aanp <- results_typearee_sp_in[[5]]

unique(sp_in_rns$Nm_spcs) # Get the total number of species for each type of protected area


# create function that groups by Red List risk category. For each category, it calculates how many distinct species there are
# and how many total observations or occurrences there are for those species
function_sp_in_cat_risk <- function(specie_in_areepr) {
  df_sp_in_aggreg <- specie_in_areepr |>
    group_by(redlistCategory) |>
    summarise(
      num_species = n_distinct(Name_species),
      num_occorrenze = n()
    ) |>
    st_drop_geometry()
  return(df_sp_in_aggreg)
}


list_df_sp_in_aggr_catrisk <- lapply(results_typearee_sp_in, function_sp_in_cat_risk)

sp_in_rns_aggr_catrisk <- list_df_sp_in_aggr_catrisk[[1]]
sp_in_pnz_aggr_catrisk <- list_df_sp_in_aggr_catrisk[[2]]
sp_in_pnr_aggr_catrisk <- list_df_sp_in_aggr_catrisk[[3]]
sp_in_rnr_aggr_catrisk <- list_df_sp_in_aggr_catrisk[[4]]
sp_in_aanp_aggr_catrisk <- list_df_sp_in_aggr_catrisk[[5]]


# function that aggregates the data by protected area and counts how many distinct species are contained in each area.
count_specie <- function(specie_in_areepr) {
  count_specie_areepr <- specie_in_areepr |>
    group_by(protected_area) |>
    summarise(num_species = n_distinct(Name_species)) |>
    arrange(desc(num_species))
  return(count_specie_areepr)
}


list_df_sp_in_num_type <- lapply(results_typearee_sp_in, count_specie)

num_sp_in_rns <- list_df_sp_in_num_type[[1]]
num_sp_in_pnz <- list_df_sp_in_num_type[[2]]
num_sp_in_pnr <- list_df_sp_in_num_type[[3]]
num_sp_in_rnr <- list_df_sp_in_num_type[[4]]
num_sp_in_aanp <- list_df_sp_in_num_type[[5]]


# creates function that show for each protected area, the total number of species
# and the number of species for each Red List risk category.
full_function <- function(specie_in_areepr) {
  species_count <- specie_in_areepr |>
    group_by(protected_area) |>
    summarise(total_species = n_distinct(Name_species))

  # Count of species by risk category
  risk_count <- specie_in_areepr |>
    group_by(protected_area, redlistCategory) |>
    summarise(risk_species = n_distinct(Name_species)) # %>%
  # spread(key = rdlstCt, value = risk_species, fill = 0)

  # Merge the two dataframes
  result <- left_join(risk_count, species_count, by = "protected_area")
  result <- result |>
    arrange(desc(total_species))
  return(result)
}


lista_full <- lapply(results_typearee_sp_in, full_function)
rns_full <- lista_full[[1]]
pnr_full <- lista_full[[3]]
rnr_full <- lista_full[[4]]
aanp_full <- lista_full[[5]]
pnz_full <- lista_full[[2]]

# Parco Nazionale dell'Arcipelago di La Maddalena does not contain species listed in the Red List

# Some protected areas do not contain species from the Red List...


# For the barplots of RNS, PNR, RNR, and AANP, I could select the top areas with the most species and create two barplots

## plots Parchi Nazionali ----

write.csv(num_sp_in_pnz, "num_sp_in_pnz.csv")

# Create a vector with the acronyms
sigle <- c(
  "PN-AS", "PN-GOG", "PN-GP", "PN-VG", "PN-ST", "PN-DB", "PN-CT", "PN-ATE",
  "PN-FMC", "PN-MS", "PN-GSML", "PN-MAI", "PN-ALM", "PN-CIR", "PN-VES", "PN-CVD",
  "PN-GAR", "PN-AM ", "PN-ALL", "PN-POL", "PN-SIL", "PN-ASP", "PN-ALM", "PN-AT"
)

# pnz$sigle <- sigle
lista_aree_pr[[2]]$sigle <- sigle

colnames(pnz)[4] <- "protected_area"
num_sp_in_pnz <- merge(num_sp_in_pnz, pnz[, c("protected_area", "sigle")], by = "protected_area", all.x = TRUE)

num_sp_in_pnz <- num_sp_in_pnz |>
  arrange(desc(num_species)) |>
  st_drop_geometry()


# Barplot showing the number of species categorized by risk category
pnz_full$redlistCategory <- factor(pnz_full$redlistCategory,
  levels = c(
    "Critically Endangered", "Endangered",
    "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"
  )
)

pnz_full <- merge(pnz_full, pnz[, c("protected_area", "sigle")], by = "protected_area", all.x = TRUE)

ggplot(pnz_full, aes(x = reorder(sigle, -total_species), y = risk_species, fill = redlistCategory)) +
  geom_bar(stat = "identity", width = 0.9) +
  labs(x = "National Park", y = "Number of species", fill = "Categoria Lista Rossa") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "whitesmoke", colour = NA)
  ) +
  scale_fill_manual(values = c(
    "Least Concern" = "green4", "Near Threatened" = "yellow",
    "Vulnerable" = "orange", "Endangered" = "red",
    "Critically Endangered" = "darkred", "Data Deficient" = "darkgrey"
  ))
dev.off()


## plots Riserve naturali statali Nazionali ----

# ggplot(num_sp_in_rns, aes(x = reorder(protected_area, -num_species), y = num_species)) +
#  geom_bar(stat = "identity")+
#  labs(x = "Riserve Naturali Statali", y = "Num species") +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels for better readability
#   panel.background = element_rect(fill = "lightgrey", colour = NA))
# dev.off()

num_sp_in_rns <- num_sp_in_rns |>
  filter(num_species >= 10)


# Barplot displaying the number of species divided by risk category from the dataset
rns_full <- factor(rns_full$redlistCategory,
  levels = c(
    "Critically Endangered", "Endangered",
    "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"
  )
)

rns_full <- rns_full |>
  filter(total_species >= 10)

ggplot(rns_full, aes(x = reorder(protected_area, -total_species), y = risk_species, fill = redlistCategory)) +
  geom_bar(stat = "identity") +
  labs(x = "Protected Area", y = "Number of species ", fill = "Categoria Lista Rossa") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "whitesmoke", colour = NA)
  ) +
  scale_fill_manual(values = c(
    "Least Concern" = "green4", "Near Threatened" = "yellow",
    "Vulnerable" = "orange", "Endangered" = "red",
    "Critically Endangered" = "darkred", "Data Deficient" = "darkgrey"
  ))
dev.off()

## plots Parchi naturali Regionali ----

# ggplot(num_sp_in_pnr, aes(x = reorder(protected_area, -num_species), y = num_species)) +
#  geom_bar(stat = "identity") +
#  labs(x = "Parchi Naturali Regionali", y = "Num species") +
# theme(
#    axis.text.x = element_text(angle = 45, hjust = 1), # Rotazione etichette sull'asse x per una migliore leggibilità
#    panel.background = element_rect(fill = "lightgrey", colour = NA)
#  )
# dev.off()

num_sp_in_pnr <- num_sp_in_pnr |>
  filter(num_species >= 50)


# Barplot showing the number of species categorized by risk category
pnr_full$redlistCategory <- factor(pnr_full$redlistCategory,
  levels = c(
    "Critically Endangered", "Endangered",
    "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"
  )
)

pnr_full <- pnr_full |>
  filter(total_species >= 50)

ggplot(pnr_full, aes(x = reorder(protected_area, -total_species), y = risk_species, fill = redlistCategory)) +
  geom_bar(stat = "identity") +
  labs(x = "Parco Regionale", y = "Numero di specie", fill = "Categoria Lista Rossa") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "whitesmoke", colour = NA)
  ) +
  scale_fill_manual(values = c(
    "Least Concern" = "green4", "Near Threatened" = "yellow",
    "Vulnerable" = "orange", "Endangered" = "red",
    "Critically Endangered" = "darkred", "Data Deficient" = "darkgrey"
  ))
dev.off()


## plots Riserve naturali Regionali ----

# ggplot(num_sp_in_rnr, aes(x = reorder(protected_area, -num_species), y = num_species)) +
#  geom_bar(stat = "identity") +
#  labs(x = "Riserve Naturali Regionali", y = "Num species") +
#  theme(
#    axis.text.x = element_text(angle = 45, hjust = 1), # Rotazione etichette sull'asse x per una migliore leggibilità
#    panel.background = element_rect(fill = "lightgrey", colour = NA)
#  )
# dev.off()

num_sp_in_rnr <- num_sp_in_rnr |>
  filter(num_species >= 20)


# Barplot displaying the number of species divided by risk category from the dataset
rnr_full$redlistCategory <- factor(rnr_full$redlistCategory,
  levels = c(
    "Critically Endangered", "Endangered",
    "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"
  )
)

rnr_full <- rnr_full |>
  filter(total_species >= 20)

ggplot(rnr_full, aes(x = reorder(protected_area, -total_species), y = risk_species, fill = redlistCategory)) +
  geom_bar(stat = "identity") +
  labs(x = "Riserva Naturale Regionale", y = "Numero di specie", fill = "Categoria Lista Rossa") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "whitesmoke", colour = NA)
  ) +
  scale_fill_manual(values = c(
    "Least Concern" = "green4", "Near Threatened" = "yellow",
    "Vulnerable" = "orange", "Endangered" = "red",
    "Critically Endangered" = "darkred", "Data Deficient" = "darkgrey"
  ))
dev.off()

## plots Altre Aree Naturali Protette ----

# ggplot(num_sp_in_aanp, aes(x = reorder(protected_area, -num_species), y = num_species)) +
#  geom_bar(stat = "identity") +
#  labs(x = "Altre Riserve Naturali", y = "Num species") +
#  theme(
#    axis.text.x = element_text(angle = 45, hjust = 1), # Rotazione etichette sull'asse x per una migliore leggibilità
#    panel.background = element_rect(fill = "lightgrey", colour = NA)
#  )
# dev.off()

num_sp_in_aanp <- num_sp_in_aanp |>
  filter(num_species >= 10)


# Barplot showing the number of species categorized by risk category
aanp_full$redlistCategory <- factor(aanp_full$redlistCategory,
  levels = c(
    "Critically Endangered", "Endangered",
    "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"
  )
)

aanp_full <- aanp_full |>
  filter(total_species >= 10)

ggplot(aanp_full, aes(x = reorder(protected_area, -total_species), y = risk_species, fill = redlistCategory)) +
  geom_bar(stat = "identity") +
  labs(x = "Altre Riserve Naturali", y = "Numero di specie", fill = "Categoria Lista Rossa") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "whitesmoke", colour = NA)
  ) +
  scale_fill_manual(values = c(
    "Least Concern" = "green4", "Near Threatened" = "yellow",
    "Vulnerable" = "orange", "Endangered" = "red",
    "Critically Endangered" = "darkred", "Data Deficient" = "darkgrey"
  ))
dev.off()


##### functions ----
# Identifying functions
functions <- ls()[sapply(ls(), function(x) is.function(get(x)))]


# Remove functions
rm(list = functions)

#rm(list_df_sp_in_aggr_catrisk, list_df_sp_in_num_type,
#              lista_full, pnz, pnz_full, results_typearee_sp_in)


# Siti Rete Natura 00 ----
sf_use_s2(FALSE) # 'sf' package S2 for handling spherical geometries
rete_nat_00 <- st_read("C:/Project_tirocinio/data/Rete_Nat_2000.shp")
rete_nat_00 <- st_make_valid(rete_nat_00)
# write.csv(rete_nat_00, "rete_nat_00.csv")
# st_write(rete_nat_00,"rete_nat_00.shp")

st_crs(rete_nat_00) == st_crs(siti_protet)

summary(rete_nat_00)
unique(rete_nat_00$tipo_sito)

# ZPS = site type A; SIC-ZSC = site type B; SIC-ZSC coinciding with ZPS = site type C
# ZPS: Special Protection Areas established by the Birds Directive
# ZSC: Special Areas of Conservation from the Habitats Directive

unique(rete_nat_00$reg_biog)

# Total area
rete_nat_00$hectares <- as.numeric(rete_nat_00$hectares)
sup_totale <- round(sum(rete_nat_00$hectares), 3)

rete_nat_aggr <- rete_nat_00 |>
  st_drop_geometry() |>
  group_by(tipo_sito) |>
  summarise(num_sit = n(), tot_area = round(sum(hectares), 3))

# grafico siti rete nat 2000 suddivisi per tipologia
rete_nat_00$tipo_sito <- factor(rete_nat_00$tipo_sito,
  levels = c("A", "B", "C")
)
ggplot() +
  coord_fixed() +
  geom_polygon(data = ita_map, aes(x = long, y = lat, group = group), colour = "gray50", fill = "gray70", alpha = 0.5) +
  geom_sf(data = rete_nat_00, aes(fill = tipo_sito), color = NA) +
  scale_fill_manual(
    values = c("A" = "red", "B" = "orange", "C" = "yellow"),
    name = "Tipologia Sito Protetto",
    labels = c("ZPS", "SIC-ZSC", "SIC-ZSC ≡ ZPS")
  ) +
  theme_minimal()
dev.off()


lista_siti_rete_nat <- split(rete_nat_00, rete_nat_00$tipo_sito)

# per dataframe della "lista_siti_rete_nat"
# mi riesco a calcolare quanti sono i siti per ogni reg biologica
# .... <- .... %>%
#  st_drop_geometry() %>%
#  group_by(reg_biog) %>%
#  summarise(num_sit = n())


## QUANTE SPECIE ----
species_gbif_rl$index_sp <- 1:nrow(species_gbif_rl)

function_sp_sititype <- function(sito_nat){
  sp_in <- st_within(species_gbif_rl, sito_nat)
  sp_in <- as.data.frame(sp_in)
  colnames(sp_in)[1] <- "index_sp"
  colnames(sp_in)[2] <- "protected_area"
  sp_in$protected_area <- sito_nat$denominazi[sp_in$protected_area]
  sp_in <- left_join(sp_in, species_gbif_rl, by = "index_sp")
  sp_in <- sp_in[,-1]
  colnames(sp_in)[5] <- "occurence_gbifID"
  sp_in <- sp_in[, c(5, 1:4, 6:ncol(sp_in))]
  sp_in <- sp_in[, c(1, 3, 2, 4:ncol(sp_in))]
  return(sp_in)
}

results_typesiti_sp_in <- lapply(lista_siti_rete_nat, function_sp_sititype)


### siti a ----
siti_a <- rete_nat_00 |>
  filter(tipo_sito == "A")

siti_a_reg_bio <- siti_a |>
  st_drop_geometry() |>
  group_by(reg_biog) |>
  summarise(num_sit = n())

sp_in_a <- results_typesiti_sp_in[[1]]

sp_in_a <- sp_in_a |>
  rename(denominazi = protected_area)

sp_in_a <- merge(sp_in_a, siti_a[, c("denominazi", "reg_biog")], by = "denominazi", all.x = T)

unique(sp_in_a$Name_species) # 578

a_num_sp_reg_bio <- sp_in_a |>
  st_drop_geometry() |>
  group_by(reg_biog) |>
  summarise(num_species = n_distinct(Name_species))

a_num_sp_cat_risk <- sp_in_a |>
  st_drop_geometry() |>
  group_by(redlistCategory) |>
  summarise(num_species = n_distinct(Name_species))


# Count of species per protected area
a_species_count <- sp_in_a |>
  group_by(protected_area) |>
  summarise(total_species = n_distinct(Name_species)) |>
  arrange(desc(total_species))

# Count of species by risk category
a_risk_count <- sp_in_a |>
  group_by(protected_area, redlistCategory) |>
  summarise(risk_species = n_distinct(Name_species)) # %>%
# spread(key = rdlstCt, value = risk_species, fill = 0)

# Merge the two dataframes
a_result <- left_join(a_risk_count, a_species_count, by = "protected_area")
a_result <- a_result |> 
  arrange(desc(total_species))

# Plot of the top protected sites with 50+ species

a_first_result <- a_result[1:58, ]

a_first_result$redlistCategory <- factor(a_first_result$redlistCategory,
  levels = c(
    "Critically Endangered", "Endangered",
    "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"
  )
)

ggplot(a_first_result, aes(x = reorder(protected_area, -total_species), y = risk_species, fill = redlistCategory)) +
  geom_bar(stat = "identity") +
  labs(x = "Siti Protetti tipo A", y = "Numero di specie", fill = "Categoria Lista Rossa") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "whitesmoke", colour = NA)
  ) +
  scale_fill_manual(values = c(
    "Least Concern" = "green4", "Near Threatened" = "yellow",
    "Vulnerable" = "orange", "Endangered" = "red",
    "Critically Endangered" = "darkred", "Data Deficient" = "darkgrey"
  ))
dev.off()



### siti b ----
siti_b <- rete_nat_00 |>
  filter(tipo_sito == "B")

siti_b_reg_bio <- siti_b |>
  st_drop_geometry() |>
  group_by(reg_biog) |>
  summarise(num_sit = n())

sp_in_b <- results_typesiti_sp_in[[2]]

sp_in_b <- sp_in_b |>
  rename(denominazi = protected_area)

sp_in_b <- merge(sp_in_b, siti_b[, c("denominazi", "reg_biog")], by = "denominazi", all.x = T)

unique(sp_in_b$Name_species) # 650

b_num_sp_reg_bio <- sp_in_b |>
  st_drop_geometry() |>
  group_by(reg_biog) |>
  summarise(num_species = n_distinct(Name_species))

b_num_sp_cat_risk <- sp_in_b |>
  st_drop_geometry() |>
  group_by(redlistCategory) |>
  summarise(num_species = n_distinct(Name_species))

# Count of species per protected area
b_species_count <- sp_in_b |>
  group_by(protected_area) |>
  summarise(total_species = n_distinct(Name_species)) |>
  arrange(desc(total_species))

# Count of species per risk category
b_risk_count <- sp_in_b |>
  group_by(protected_area, redlistCategory) |>
  summarise(risk_species = n_distinct(Name_species)) # %>%
# spread(key = rdlstCt, value = risk_species, fill = 0)

# Merge the two dataframes
b_result <- left_join(b_risk_count, b_species_count, by = "protected_area")
b_result <- b_result |>
  arrange(desc(total_species))

# Plot of the top protected sites with 50+ species

b_first_result <- b_result[1:59, ]

b_first_result$redlistCategory <- factor(b_first_result$redlistCategory,
  levels = c(
    "Critically Endangered", "Endangered",
    "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"
  )
)

ggplot(b_first_result, aes(x = reorder(protected_area, -total_species), y = risk_species, fill = redlistCategory)) +
  geom_bar(stat = "identity") +
  labs(x = "Siti Protetti tipo A", y = "Numero di specie", fill = "Categoria Lista Rossa") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "whitesmoke", colour = NA)
  ) +
  scale_fill_manual(values = c(
    "Least Concern" = "green4", "Near Threatened" = "yellow",
    "Vulnerable" = "orange", "Endangered" = "red",
    "Critically Endangered" = "darkred", "Data Deficient" = "darkgrey"
  ))
dev.off()


### siti c ----
siti_c <- rete_nat_00 |>
  filter(tipo_sito == "C")

siti_c_reg_bio <- siti_c |>
  st_drop_geometry() |>
  group_by(reg_biog) |>
  summarise(num_sit = n())

sp_in_c <- results_typesiti_sp_in[[3]]

sp_in_c <- sp_in_c |>
  rename(denominazi = protected_area)

sp_in_c <- merge(sp_in_c, siti_c[, c("denominazi", "reg_biog")], by = "denominazi", all.x = T)

unique(sp_in_c$Name_species) # 492

c_num_sp_reg_bio <- sp_in_c |>
  st_drop_geometry() |>
  group_by(reg_biog) |>
  summarise(num_species = n_distinct(Name_species))

c_num_sp_cat_risk <- sp_in_c |>
  st_drop_geometry() |>
  group_by(redlistCategory) |>
  summarise(num_species = n_distinct(Name_species))

# Count the number of species for each protected area
c_species_count <- sp_in_c |>
  group_by(protected_area) |>
  summarise(total_species = n_distinct(Name_species)) |>
  arrange(desc(total_species))

# Count the number of species for each risk category
c_risk_count <- sp_in_c |>
  group_by(protected_area, redlistCategory) |>
  summarise(risk_species = n_distinct(Name_species)) # %>%
# spread(key = rdlstCt, value = risk_species, fill = 0)

# Merge the two dataframes
c_result <- left_join(c_risk_count, c_species_count, by = "protected_area")
c_result <- c_result |>
  arrange(desc(total_species))

# Plot of the top protected sites with 50+ species

c_first_result <- c_result[1:26, ]

c_first_result$redlistCategory <- factor(c_first_result$redlistCategory,
  levels = c(
    "Critically Endangered", "Endangered",
    "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"
  )
)

ggplot(c_first_result, aes(x = reorder(protected_area, -total_species), y = risk_species, fill = redlistCategory)) +
  geom_bar(stat = "identity") +
  labs(x = "Siti Protetti tipo A", y = "Numero di specie", fill = "Categoria Lista Rossa") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "whitesmoke", colour = NA)
  ) +
  scale_fill_manual(values = c(
    "Least Concern" = "green4", "Near Threatened" = "yellow",
    "Vulnerable" = "orange", "Endangered" = "red",
    "Critically Endangered" = "darkred", "Data Deficient" = "darkgrey"
  ))
dev.off()


sp_only_out_euap <- unique(df_sp_only_out$Name_species)

# Species in sp_only_out are located within site A?
sp_only_out_in_a <- sp_only_out_euap %in% sp_in_a$Name_species

sp_only_out_in_a <- sp_only_out_euap[sp_only_out_in_a]

# Species in sp_only_out are located within site B?
sp_only_out_in_b <- sp_only_out_euap %in% sp_in_b$Name_species

sp_only_out_in_b <- sp_only_out_euap[sp_only_out_in_b]

# Species in sp_only_out are located within site C?
sp_only_out_in_c <- sp_only_out_euap %in% sp_in_c$Name_species

sp_only_out_in_c <- sp_only_out_euap[sp_only_out_in_c]

sp_only_out_in_nat_00 <- unique(c(sp_only_out_in_a, sp_only_out_in_b, sp_only_out_in_c))

# Out of the 188 species in sp_only_out euap, 85 are located within Natura 2000 sites!
# Which species are they?

sp_only_out_nat_00 <- setdiff(sp_only_out_euap, sp_only_out_in_nat_00)

species_gbif_rl_total_out <- species_gbif_rl |>
  filter(Name_species %in% sp_only_out_nat_00)

class(species_gbif_rl_total_out)

species_gbif_rl_total_out <- species_gbif_rl_total_out |>
  distinct()

unique(species_gbif_rl_total_out$Name_species)

unique(species_gbif_rl_total_out$redlistCategory)

# Plot occurrences of species both outside protected areas and from Natura 2000 sites ----

ita_map <- map_data("italy")

species_gbif_rl_total_out$redlistCategory <- factor(species_gbif_rl_total_out$redlistCategory,
  levels = c(
    "Critically Endangered", "Endangered",
    "Vulnerable", "Near Threatened",
    "Least Concern", "Data Deficient"
  )
)


ggplot() +
  # coord_fixed() +
  geom_polygon(data = ita_map, aes(x = long, y = lat, group = group), colour = "gray50", fill = "gray70", alpha = 0.5) +
  geom_sf(data = siti_protet, color = NA, fill = "darkgreen", alpha = 0.4, size = 0.1) +
  geom_sf(data = rete_nat_00, color = NA, fill = "chartreuse2", alpha = 0.4, size = 0.1) +
  geom_sf(data = species_gbif_rl_total_out, aes(color = redlistCategory), size = 0.8) +
  scale_color_manual(
    name = "RedList Category",
    values = c(
      "Critically Endangered" = "darkred", "Endangered" = "red",
      "Vulnerable" = "orange2", "Near Threatened" = "yellow",
      "Least Concern" = "burlywood2", "Data Deficient" = "cyan2"
    ),
    labels = c(
      "Critically Endangered (CE)", "Endangered (EN)", "Vulnerable(VU)",
      "Near Threatened (NT)", "Least Concern (LC)", "Data Deficient (DD)"
    ),
    guide = guide_legend(override.aes = list(size = 1)) # cambiare dimensione legenda
  ) +
  theme_minimal()
dev.off()
                         
