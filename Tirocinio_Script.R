# Progetto Tesi: Importance of Protected Areas in Italy
# This analysis aims to ascertain how many IUCN Red List (Italian) species are contained within Italy's protected areas (Euap Sites).

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
  user = "user", pwd = "pwd", email = "email"
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

gbif_it <- gbif_it |>
  filter(!is.na(decimalLatitude) & !is.na(decimalLongitude))

# Flag problems
flags <- clean_coordinates(x = gbif_it, 
                           lon = "decimalLongitude", 
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                           tests = c("capitals", "centroids", "equal","gbif", "institutions",
                                     "zeros", "countries")) 

summary(flags)
plot_flags <- plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")
plot_flags

# Exlcude problematic records
gbif_it <- gbif_it[flags$.summary,]

# The flagged records
#dat_fl <- gbif_it[!flags$.summary,]


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

# Gbif Species Catalog in the Red List ----
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

gbif_it_sf <- st_as_sf(gbif_it, coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +datum=WGS84")

# DISPLAY the occurrences in 'species_gbif_rl' ----
ita_map <- map_data("italy")

# plot_1 <- ggplot() +
#   # coord_fixed() +
#   geom_polygon(data = ita_map, aes(x = long, y = lat, group = group), colour = "gray50", fill = "gray50", alpha = 0.5) +
#   geom_sf(data = species_gbif_rl, colour = "darkred", size = 0.01, alpha = 0.5) +
#   theme_minimal()
# dev.off()


# Combine and preprocess both datasets

gbif_coords <- st_coordinates(gbif_it_sf$geometry)
gbif_rl_coords <- st_coordinates(species_gbif_rl$geometry)

# Create hexbin plots to get counts

gbif_hex <- ggplot() +
  geom_polygon(
    data = ita_map, aes(x = long, y = lat, group = group),
    fill = "gray30", size = 0.2, alpha = 0.3
  ) +
  stat_bin_hex(data = as.data.frame(gbif_coords), aes(x = X, y = Y), bins = 50)

gbif_rl_hex <- ggplot() +
  geom_polygon(
    data = ita_map, aes(x = long, y = lat, group = group),
    fill = "gray30", size = 0.2, alpha = 0.3
  ) +
  stat_bin_hex(data = as.data.frame(gbif_rl_coords), aes(x = X, y = Y), bins = 50)

# Extract hexbin counts

gbif_counts <- ggplot_build(gbif_hex)$data[[2]]$count
gbif_rl_counts <- ggplot_build(gbif_rl_hex)$data[[2]]$count

# Get the maximum count

global_max_count <- max(gbif_counts)
global_max_count_rl <- max(gbif_rl_counts)

# Build the final plots

gbif_hex_plot <- ggplot() +
  geom_polygon(
    data = ita_map, aes(x = long, y = lat, group = group),
    fill = "gray30", size = 0.2, alpha = 0.3
  ) +
  stat_bin_hex(
    data = as.data.frame(gbif_coords), aes(x = X, y = Y),
    bins = 50, color = "gray20", size = 0.2
  ) +
  scale_fill_viridis(
    option = "G", limits = c(0, global_max_count),
    direction = -1, transform = "sqrt"
  ) +
  theme_void() +
  labs(
    title = "GBIF",
    x = "Longitude",
    y = "Latitude",
    fill = "No of\noccurrences"
  ) +
  theme(
    legend.text = element_text(size = 4), legend.title = element_text(
      size = 8, #7
      vjust = 0.8
    ),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

gbif_rl_hex_plot <- ggplot() +
  geom_polygon(
    data = ita_map, aes(x = long, y = lat, group = group),
    fill = "gray30", size = 0.2, alpha = 0.3
  ) +
  stat_bin_hex(
    data = as.data.frame(gbif_rl_coords), aes(x = X, y = Y),
    bins = 50, color = "black", size = 0.2
  ) +
  scale_fill_viridis(
    option = "G", limits = c(0, global_max_count_rl),
    direction = -1, transform = "sqrt"
  ) +
  theme_void() +
  labs(
    title = "GBIF_RED_LIST",
    x = "Longitude",
    y = "Latitude",
    fill = "Number of\noccurrences"
  ) +
  theme(
    legend.text = element_text(size = 4), legend.title = element_text(
      size = 8, #7
      vjust = 0.8
    ),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Arrange the plots side by side for comparison

combined_plot <- ggarrange(gbif_hex_plot, gbif_rl_hex_plot,
                           ncol = 2, nrow = 1,
                           common.legend = FALSE, legend = "bottom"
)

# Display the combined plot

print(combined_plot)


# DIFFERENZE SPECIE_GBIF_RL FIRST AND LAST DATASET ----
species_gbif_rl_first <- st_read("species_gbif_rl.shp")
species_gbif_rl_first <- species_gbif_rl_first |> 
  filter(year >= 1985 & !is.na(year)) # seleziono solo le specie dal 1985 ad OGGI 

species_gbif_rl_no_geom <- st_drop_geometry(species_gbif_rl)

## specie gbif_rl che sono presenti solo nel dataset FIRST ----
colnames(species_gbif_rl_first)[1] <- "Name_species"

#species_gbif_rl_miss_prova <-  anti_join(species_gbif_rl_first_no_geom, species_gbif_rl_no_geom, by = "Name_species")
# 4 specie miss; 

species_gbif_rl_miss <-  anti_join(species_gbif_rl_first, species_gbif_rl_no_geom, by = "Name_species")
unique(species_gbif_rl_miss$Name_species) #140 specie perse dal 1° dataset...

# suddivsione delle 140 specie per categorie di rischio
species_gbif_rl_miss_cat_risk <- species_gbif_rl_miss |> 
  group_by(rdlstCt) |>
  summarise(
    num_species = n_distinct(Name_species)
  ) |> 
  st_drop_geometry()

# suddivisione delle 140 specie per Basis of Record
species_gbif_rl_miss_basis_of_Rc <- species_gbif_rl_miss |> 
  group_by(bssOfRc) |>
  summarise(
    num_species = n_distinct(Name_species)
  ) |> 
  st_drop_geometry() # non tornano le 140 perchè alcune specie posso avere più tipologie di basis of Record

# cercare di capire perchè tot specie sono perse dal 1° dataset
# species_gbif_rl_first_no_geom <- st_drop_geometry(species_gbif_rl_first)
# 
# species_gbif_rl_first_no_geom <- species_gbif_rl_first_no_geom |> 
#   filter(!is.na(year) & !is.na(crdnUIM))

# delle specie che ho perso seleziono quelle con la giusta tipologia di Basis of Record
species_gbif_rl_miss_ho <- species_gbif_rl_miss |> 
  filter(bssOfRc == "HUMAN_OBSERVATION" | bssOfRc == "OCCURRENCE")
unique(species_gbif_rl_miss_ho$Name_species) # 61 specie che hanno crdnUIM = "NA"

# le raggruppo in base all'anno per capire quali eliminare 
group_year <- species_gbif_rl_miss_ho |> 
  group_by(year) |> 
  summarise(
    num_species = n_distinct(Name_species)) |> 
  st_drop_geometry()


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
plot_aree_pr <- ggplot() +
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
## Count occurrences IN the protected areas ----

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
sp_condivise <- inner_join(sp_out_aree_pr_df, sp_in_aree_pr_df, by = "Name_species")
#nrow(sp_condivise_in) + nrow(sp_condivise_out)

length(unique(sp_condivise_in$Name_species))
length(unique(sp_condivise_out$Name_species))
length(unique(species_gbif_rl$Name_species))

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

# write.csv(df_sp_only_in,"C:/Project_tirocinio/OUTPUT_TIROCINIO_2/df_sp_only_in.csv")
# read.csv

# Group species and occurrences by risk category
sp_only_in_cat_risk <- df_sp_only_in |>
  group_by(redlistCategory) |>
  summarise(
    num_species = n_distinct(Name_species),
    num_occorrenze = n()
  )

# Division by risk categories
list_sp_only_in_cat_risk <- split(df_sp_only_in, df_sp_only_in$redlistCategory)
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

#write.csv(df_sp_only_out, "C:/Project_tirocinio/OUTPUT_TIROCINIO_2/df_sp_only_out.csv")
# read.csv

# Group species and occurrences by risk category
sp_only_out_cat_risk <- df_sp_only_out |>
  group_by(redlistCategory) |>
  summarise(
    num_species = n_distinct(Name_species),
    num_occorrenze = n()
  ) |>
  st_drop_geometry()

# Categorization by risk categories
list_sp_only_out_cat_risk <- split(df_sp_only_out, df_sp_only_out$redlistCategory)
# sp_only_out_DD
# sp_only_out_LC
sp_only_out_NT <- list_sp_only_out_cat_risk[[5]]
sp_only_out_V <- list_sp_only_out_cat_risk[[6]]
sp_only_out_E <- list_sp_only_out_cat_risk[[3]]
sp_only_out_CE <- list_sp_only_out_cat_risk[[1]]

unique(sp_only_out_CE$Name_species)

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
plot_sp_only_out <- ggplot() +
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
plot_3
dev.off()

# per pulire tutto l'environment
rm(list = ls())


### Sp SOLO DENTRO ----
sp_only_in_low_protect <- df_sp_only_in %>%
  group_by(Nm_spcs) %>%
  summarise(
    num_areepr = n_distinct(protected_area),
    num_occorrenze = n()
  ) %>%
  arrange(num_areepr)

# rm(plot_3, list_sp_only_in_cat_risk, list_sp_only_out_cat_risk, lista_sp_cat_rischio,
#    names_correct, sp_condivise_cat_risk, sp_condivise_in, sp_condivise_out,
#    sp_only_in_cat_risk, sp_only_out_cat_risk, sp_only_out_CE,
#    sp_only_out_E, sp_only_out_NT, sp_only_out_V, sp_out)


# ANALYSIS FOR EACH TYPE OF PROTECTED AREA ----
lista_aree_pr <- split(siti_protet, siti_protet$tipo)

lista_aree_pr[[2]] <- lista_aree_pr[[2]] |>
  filter(nome_gazze != "Parco Nazionale dell'Arcipelago di La Maddalena")

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

# Datasets  containing occurrences and species for each type of protected area
sp_in_rns <- results_typearee_sp_in[[1]]
sp_in_pnz <- results_typearee_sp_in[[2]]
sp_in_pnr <- results_typearee_sp_in[[3]]
sp_in_rnr <- results_typearee_sp_in[[4]]
sp_in_aanp <- results_typearee_sp_in[[5]]

unique(sp_in_aanp$Name_species) # Get the total number of species for each type of protected area

unique(sp_in_aanp$protected_area) # Get the total number of protected area containing at least 1 red list species

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
# count_specie <- function(specie_in_areepr) {
#   count_specie_areepr <- specie_in_areepr |>
#     group_by(protected_area) |>
#     summarise(num_species = n_distinct(Name_species)) |>
#     arrange(desc(num_species))
#   return(count_specie_areepr)
# }
# 
# 
# list_df_sp_in_num_type <- lapply(results_typearee_sp_in, count_specie)
# 
# num_sp_in_rns <- list_df_sp_in_num_type[[1]]
# num_sp_in_pnz <- list_df_sp_in_num_type[[2]]
# num_sp_in_pnr <- list_df_sp_in_num_type[[3]]
# num_sp_in_rnr <- list_df_sp_in_num_type[[4]]
# num_sp_in_aanp <- list_df_sp_in_num_type[[5]]


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
pnz_full <- lista_full[[2]]
pnr_full <- lista_full[[3]]
rnr_full <- lista_full[[4]]
aanp_full <- lista_full[[5]]


# Some protected areas do not contain species from the Red List...


## plots Parchi Nazionali ----

#write.csv(num_sp_in_pnz, "num_sp_in_pnz.csv")

unique(sp_in_pnz$Name_species) # Get the total number of species for each type of protected area

unique(sp_in_pnz$protected_area) # Get the total number of protected area containing at least 1 red list species

# Create a vector with the acronyms
sigle <- c(
  "PN-AS", "PN-GOG", "PN-GP", "PN-VG", "PN-ST", "PN-DB", "PN-CT", "PN-ATE",
  "PN-FMC", "PN-MS", "PN-GSML", "PN-MAI", "PN-ALM", "PN-CIR", "PN-VES", "PN-CVD",
  "PN-GAR", "PN-AM ", "PN-ALL", "PN-POL", "PN-SIL", "PN-ASP", "PN-ALM", "PN-AT"
)

# pnz$sigle <- sigle
lista_aree_pr[[2]]$sigle <- sigle

pnz <- lista_aree_pr[[2]]
colnames(pnz)[4] <- "protected_area"

#num_sp_in_pnz <- merge(num_sp_in_pnz, pnz[, c("protected_area", "sigle")], by = "protected_area", all.x = TRUE)

# num_sp_in_pnz <- num_sp_in_pnz |>
#   arrange(desc(num_species)) |>
#   st_drop_geometry()


# Barplot showing the number of species categorized by risk category
pnz_full$redlistCategory <- factor(pnz_full$redlistCategory,
  levels = c(
    "Critically Endangered", "Endangered",
    "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"
  )
)

pnz_full <- merge(pnz_full, pnz[, c("protected_area", "sigle")], by = "protected_area", all.x = TRUE)

geom_bar_pnz <- ggplot(pnz_full, aes(x = reorder(sigle, -total_species), y = risk_species, fill = redlistCategory)) +
  geom_bar(stat = "identity", width = 0.9) +
  labs(x = "National Park", y = "Number of species", fill = "Red List Category") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "whitesmoke", colour = NA)
  ) +
  scale_fill_manual(values = c(
    "Least Concern" = "green4", "Near Threatened" = "yellow",
    "Vulnerable" = "orange", "Endangered" = "red",
    "Critically Endangered" = "darkred", "Data Deficient" = "darkgrey"
  ))

geom_bar_pnz
dev.off()


# For the barplots of RNS, PNR, RNR, and AANP, I could select the top areas with the most species and create two barplots

## plots Riserve naturali statali Nazionali ----

# ggplot(num_sp_in_rns, aes(x = reorder(protected_area, -num_species), y = num_species)) +
#  geom_bar(stat = "identity")+
#  labs(x = "Riserve Naturali Statali", y = "Num species") +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels for better readability
#   panel.background = element_rect(fill = "lightgrey", colour = NA))
# dev.off()

# num_sp_in_rns <- num_sp_in_rns |>
#   filter(num_species >= 10)


# Barplot displaying the number of species divided by risk category from the dataset
rns_full$redlistCategory <- factor(rns_full$redlistCategory,
  levels = c(
    "Critically Endangered", "Endangered",
    "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"
  )
)

rns_full_fltr <- rns_full |>
  filter(total_species >= 10)

geom_bar_rns <- ggplot(rns_full_fltr, aes(x = reorder(protected_area, -total_species), y = risk_species, fill = redlistCategory)) +
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

geom_bar_rns
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

# num_sp_in_pnr <- num_sp_in_pnr |>
#   filter(num_species >= 50)


# Barplot showing the number of species categorized by risk category
pnr_full$redlistCategory <- factor(pnr_full$redlistCategory,
  levels = c(
    "Critically Endangered", "Endangered",
    "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"
  )
)

pnr_full_fltr <- pnr_full |>
  filter(total_species >= 50)

geom_bar_pnr <- ggplot(pnr_full_fltr, aes(x = reorder(protected_area, -total_species), y = risk_species, fill = redlistCategory)) +
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


rnr_full_fltr <- rnr_full |>
  filter(total_species >= 20 | redlistCategory == "Critically Endangered")

geom_bar_rnr <- ggplot(rnr_full_fltr, aes(x = reorder(protected_area, -total_species), y = risk_species, fill = redlistCategory)) +
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

geom_bar_rnr
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

# num_sp_in_aanp <- num_sp_in_aanp |>
#   filter(num_species >= 10)


# Barplot showing the number of species categorized by risk category
aanp_full$redlistCategory <- factor(aanp_full$redlistCategory,
  levels = c(
    "Critically Endangered", "Endangered",
    "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"
  )
)

aanp_full_fltr <- aanp_full |>
  filter(total_species >= 10)

geom_bar_aanp <- ggplot(aanp_full_fltr, aes(x = reorder(protected_area, -total_species), y = risk_species, fill = redlistCategory)) +
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

geom_bar_aanp
dev.off()


##### functions ----
# Identifying functions
functions <- ls()[sapply(ls(), function(x) is.function(get(x)))]

# Remove functions
rm(list = functions)

                         
# Statistical ANALYSIS ----
mapview_species_rl = mapview(species_gbif_rl, cex = 3, alpha = .5, popup = NULL)

mapview_species_rl

# confini dell'italia

italy_boundary <- ne_countries(country = "Italy", scale = "medium", returnclass = "sf")
ggplot(data = italy_boundary) +
  geom_sf()
dev.off()

# Convertire il CRS in un sistema di coordinate in metri (es. EPSG:3035)

italy_boundary <- st_transform(italy_boundary, crs = 3035)

st_crs(italy_boundary)

# Creazione griglia di quadrati

net_grid <- st_make_grid(italy_boundary, cellsize = c(10000, 10000), square = TRUE)
plot(net_grid)
dev.off()

net_grid_sf <- st_sf(net_grid)

# Uso sparse = FALSE per ottenere una matrice booleana

intersections <- st_intersects(net_grid_sf, italy_boundary, sparse = FALSE)

# Filtrare net_grid utilizzando la matrice delle intersezioni

net_grid_sf_italy <- net_grid_sf[rowSums(intersections) > 0, ]

plot(net_grid_sf_italy[1])
dev.off()

# Indicizzazione dei quadrati della griglia per facilitare analisi 

net_grid_sf_italy <- net_grid_sf_italy |> 
  mutate(grid_id = 1:n())

st_crs(net_grid_sf_italy)

net_grid_sf_italy <- st_make_valid(net_grid_sf_italy)

siti_protet_union <- st_union(siti_protet)

siti_protet_union <- st_make_valid(siti_protet_union)

siti_protet_union <- st_transform(siti_protet_union, crs = st_crs(net_grid_sf_italy))

st_crs(net_grid_sf_italy) == st_crs(siti_protet_union)

# Intersezione tra griglia e aree protette per trovare la superficie protetta
intersezioni <- st_intersection(net_grid_sf_italy, siti_protet_union)

str(intersezioni)

# Rendi valide le geometrie
intersezioni <- st_make_valid(intersezioni)

# # Rimuovi le intersezione duplicate basandoti sulla geometria
# dddddd <- intersezioni |> 
#   distinct(net_grid, .keep_all = TRUE)
# rm(dddddd)


# intersezioni_sf <- st_sf(intersezioni)

# Unisci la geometria della griglia e delle intersezioni
# intersezioni_plot <- ggplot() +
#   geom_sf(data = net_grid_sf_italy, fill = "red", color = "darkred", alpha = 0.3) +  # Griglia di quadrati
#   geom_sf(data = intersezioni_sf, fill = "blue", color = "darkblue", alpha = 0.5) +  # Aree protette
#   labs(title = "Intersezioni tra Griglia e Aree Protette", x = "Longitudine", y = "Latitudine") +
#   theme_minimal() +
#   theme(legend.position = "none")
# intersezioni_plot
# dev.off()

# Sistema di riferimento adeguato per st_area in m^2

st_crs(intersezioni) # 3035 va bene, invece crs 4326 no, perchè è in gradi

# intersezioni <- st_transform(intersezioni, crs = 3035) # crs in metri


# Calcola l'area delle intersezioni

intersezioni$area_intersecata <- st_area(intersezioni)


# Aggrega e conta l'area protetta (delle intersezioni) per ciascuna cella della griglia
area_protetta_per_cella <- intersezioni |> 
  group_by(grid_id) |> 
  summarize(area_protetta = sum(area_intersecata, na.rm = TRUE))

# A questo punto, ho un dataframe con l'area protetta per ogni quadrato

# Calcola la superficie totale di ciascuna cella della griglia

net_grid_sf_italy$area_totale <- st_area(net_grid_sf_italy)


net_grid_italy <- st_drop_geometry(net_grid_sf_italy)

area_protetta_per_cella_nogeom <- st_drop_geometry(area_protetta_per_cella)                                  

net_grid_italy <- left_join(net_grid_italy, area_protetta_per_cella_nogeom, by = "grid_id")

net_grid_italy <- net_grid_italy |> 
  mutate(percentuale_protetta = round((area_protetta / area_totale) * 100,1))

net_grid_italy$area_protetta[is.na(net_grid_italy$area_protetta)] <- 0
net_grid_italy$percentuale_protetta[is.na(net_grid_italy$percentuale_protetta)] <- 0



# lavoro con le occ. totali "species_gbif_rl" e con le occ. "sp_in_aree_pr_df"

species_gbif_rl <- st_transform(species_gbif_rl, crs = st_crs(net_grid_sf_italy))

species_gbif_rl_w_grid <- st_join(species_gbif_rl, net_grid_sf_italy, join = st_within)

# Calcola il numero di occorrenze totali per cella

occ_tot_per_cella <- species_gbif_rl_w_grid |> 
  group_by(grid_id) |> 
  summarise(occ_totali = n(), .groups = 'drop') |> 
  st_drop_geometry()

# Aggiungi una colonna 'grid_id' al dataframe delle occorrenze protette

sp_in_aree_pr_sf <- st_as_sf(sp_in_aree_pr_df)

sp_in_aree_pr_sf <- st_transform(sp_in_aree_pr_sf, crs = 3035)

# st_crs(sp_in_aree_pr_sf) == st_crs(net_grid_sf_italy)

sp_in_aree_pr_w_grid <- st_join(sp_in_aree_pr_sf, net_grid_sf_italy, join = st_within)

# Calcola il numero di occorrenze protette per cella

occ_prot_per_cella <- sp_in_aree_pr_w_grid |> 
  group_by(grid_id) |> 
  summarise(occ_protette = n(), .groups = 'drop') |> 
  st_drop_geometry()

# Unisci i conteggi delle occorrenze totali e protette nel dataframe della griglia

net_grid_italy <- net_grid_italy |> 
  left_join(occ_tot_per_cella, by = "grid_id") |> 
  left_join(occ_prot_per_cella, by = "grid_id")

# 

net_grid_italy <- net_grid_italy |> 
  mutate(percentuale_occ_protette = round((occ_protette / occ_totali) * 100,1))


net_grid_italy <- net_grid_italy |> 
  mutate(
    occ_totali = replace_na(occ_totali, 0),
    occ_protette = replace_na(occ_protette, 0),
    percentuale_occ_protette = replace_na(percentuale_occ_protette, 0)
  )

write.csv(net_grid_italy, "linear_mod.csv")

# Calcola la regressione lineare
mod <- lm(percentuale_occ_protette ~ as.numeric(percentuale_protetta), data = net_grid_italy)

summary(mod)

# Aggiungi i valori predetti e calcola lo scarto

net_grid_italy <- net_grid_italy |> 
  mutate(predicted = predict(mod), 
         scarto = percentuale_occ_protette - predicted)

# plot
scatter_plot_final <- ggplot(net_grid_italy, aes(x = as.numeric(percentuale_protetta), y = percentuale_occ_protette))+
  geom_point(size = 0.9) +
  geom_line(aes(y = predicted), color = "blue", linetype = "dashed") +  # retta di regressione
  #geom_segment(aes(xend = as.numeric(percentuale_protetta), yend = predicted), color = "red", linetype = "dotted") +  # scarti
  labs(title = "Scatter Plot EUAP sites",
       x = "Protected area %",
       y = "Protected occurrences %") +
  annotate("text", x = 90, y = 20, 
           label = paste("Slope: ", round(slope, 2), "\n",
           "t-value: 76.77", "\n",
           "p-value: ", p_value, "\n",
           "R^2: ", round(r_squared, 4)),
            size = 4, color = "black", hjust = 0) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
    # panel.grid.major = element_blank(),  # Removes major grid lines
    # panel.grid.minor = element_blank()   # Removes minor grid lines
  )

scatter_plot_final
dev.off()


scarto <- net_grid_italy[,c("grid_id","scarto")]
net_grid_sf_italy <- left_join(net_grid_sf_italy, scarto, by = "grid_id")

st_crs(net_grid_sf_italy)

# italy_shp_region <- st_read("C:/Project_tirocinio/data/gadm41_ITA_1.shp")
italy_shp_region <- st_transform(italy_shp_region, crs = st_crs(net_grid_sf_italy))

color <- colorRampPalette(c("darkred", "white", "darkgreen"))(50)
plot_finale <- ggplot()+
  geom_sf(data = net_grid_sf_italy, aes(fill = scarto), colour = "black") +
  # geom_sf(data = italy_shp_region, fill = NA, color = "grey20", alpha = 0.1) +
  scale_fill_gradientn(colors = color, name = "Residuals") +
  labs(title = "Residuals Sites EUAP") +
  theme_minimal() +
  theme( 
    legend.title = element_text(size = 12, face = "plain"),  # Titolo della legenda
    legend.text = element_text(size = 10),  # Testo delle etichette
    legend.key.size = unit(0.6, "cm"),# Dimensione dei simboli nella legenda
    # legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_blank(),     # Remove x-axis labels
    axis.text.y = element_blank(),      # Remove y-axis labels  
    # panel.background = element_blank(),    # Remove panel background
    plot.background = element_blank(),     # Remove plot background
    panel.grid = element_blank()           # Remove grid lines
  )
plot_finale
dev.off()

# Estrai i coefficienti e il R^2
slope <- coef(mod)[2]  # pendenza
intercept <- coef(mod)[1]  # intercetta
r_squared <- summary(mod)$r.squared  # R^2
p_value <- "<2e-16***"# summary(mod)$coefficients[2, 4]  # p-value per la pendenza                
