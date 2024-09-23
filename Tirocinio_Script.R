#Progetto Tesi: Importanza delle Aree Protette in Italia
#Quest' analisi ha l'obiettivo di verificare quante specie della Red List IUCN (italiana) sono contenute all'interno delle aree protette d'Italia (Siti Euap e Siti Rete-Natura 00).

#install.packages("CoordinateCleaner")
#install.packages("countrycode")
#install.packages("TNRS")
#install.packages("dtrackr")
#install.packages("DescTools")
#install.packages("dbscan")
#install.packages("mapview")
library(mapview)
library(countrycode)
library(CoordinateCleaner)
library(dplyr)
library(ggplot2)
library(rgbif)
library(sp)
library(maps)
library(TNRS)
library(sf)
#llibrary(dbscan)
#library(DescTools)
#library(dtrackr)

dat <- read.csv("data/Tracheofite_Gbif_nuovo.csv", sep = "\t")
names(dat) #a lot of columns
glimpse(dat)
str(dat)

# Pulizia del dataset Gbif ----
# Select columns of interest
dat <- dat %>%
  dplyr::select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
                gbifID, family, taxonRank, coordinateUncertaintyInMeters, year,
                basisOfRecord, institutionCode)

# Remove records without coordinates
dat <- dat%>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude))

# Plot data to get an overview
wm <- map_data("world")
ggplot() + 
  coord_fixed() + 
  geom_polygon(data = wm, aes(x = long, y = lat, group = group), colour = "gray50", fill = "gray50") +
  geom_point(data = dat, aes(x = decimalLongitude, y = decimalLatitude), colour = "darkred", size = 1) +
  theme_bw()

# Select records only in Italian coordinates
dat <- dat %>%
  filter(!(decimalLongitude < 6 | decimalLongitude > 18)) %>%
  filter(!(decimalLatitude < 36 | decimalLatitude > 47))

# Convert country code from ISO2c to ISO3c
dat$countryCode <-  countrycode(dat$countryCode, origin =  'iso2c', destination = 'iso3c')

# Check that all records are in the ITA country
prova <-filter(dat, countryCode == "ITA")

# dat$decimalLatitude <- as.numeric(dat$decimalLatitude)
summary(dat)

# Flag problems
flags <- clean_coordinates(x = dat, 
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
dat_cl <- dat[flags$.summary,]

# The flagged records
dat_fl <- dat[!flags$.summary,]

write.csv(dat_cl, "dataset_gbif_pulito_pt1.csv")
dat_cl <- read.csv("dataset_gbif_pulito_pt1.csv")

# write.table(dat_cl, "dataset_pulito_pt1.csv", sep="\t")
# dat <- read.csv("dataset_pulito_pt1.csv", sep= "\t")

# Remove records with low coordinate precision
summary(dat_cl$coordinateUncertaintyInMeters)
hist(dat_cl$coordinateUncertaintyInMeters / 1000, breaks = 20)

# some attempts
dat_cl <- dat_cl %>%
  filter(coordinateUncertaintyInMeters / 1000 <= 100 | is.na(coordinateUncertaintyInMeters))
dat_cl <- dat_cl %>%
  filter(coordinateUncertaintyInMeters <= 20 | is.na(coordinateUncertaintyInMeters))
dat_cl <- dat_cl %>%
  filter(coordinateUncertaintyInMeters <= 5 | is.na(coordinateUncertaintyInMeters))

hist(dat_cl$coordinateUncertaintyInMeters, breaks = 5)

# Remove unsuitable data sources, especially fossils Which are responsible for the majority of problems in this case
table(dat$basisOfRecord)

dat_cl <- filter(dat_cl, basisOfRecord == "HUMAN_OBSERVATION" | 
                   basisOfRecord == "OBSERVATION" |
                   basisOfRecord == "PRESERVED_SPECIMEN" |
                   basisOfRecord == "LIVING_SPECIMEN" |
                   basisOfRecord == "MATERIAL_SAMPLE" |
                   basisOfRecord == "OCCURRENCE")

write.csv(dat_cl, "dataset_gbif_pulito_pt2.csv")

dat <- read.csv("dataset_gbif_pulito_pt2.csv")

# Taxonomy standardisation 
dat <- filter(dat, species != "")

names_species <- unique(dat$species)
is.vector(names_species)

names_correct <- TNRS(names_species)

names_correct <- names_correct %>%
  dplyr::select(Name_submitted, Name_matched)

write.csv(names_correct, "dateset_nomi_corretti_Gbif.csv")
names_correct <- read.csv("dateset_nomi_corretti_Gbif.csv", sep=",")

colnames(names_correct)[1] <- "species"

dat <- merge(x=dat,y=names_correct, 
             by="species", all.x=TRUE)
dat <- dat[,-1]

dat <- dat %>%   select(Name_matched, everything())  
colnames(dat)[1] <- "Name species"
# dat <- dat[ , -c(2,3)]

# dataset pulito finale
# write.table(dat, "dataset_pulito_finale.csv", sep=",")
dat_pulito <- read.csv("dataset_pulito_finale.csv", sep=",")


## Filtrazione delle occorrenze Gbif in base al dato temporale----
dat_pulito <- read.csv("dataset_pulito_finale.csv", sep=",")
summary(dat_pulito)
dat_pulito <- dat_pulito %>%
  filter(year >= 1985 | is.na(year))

write.csv(dat_pulito, "dat_pulito_tempo.csv")

# Quante specie ci sono
species <- unique(dat_pulito$Name.species)
rm(species)

frequenze_sp <- Freq(dat_pulito$Name.species, ord = "desc")

# prendo le prime 26 specie che sono distribuite/coprono il 10% degli individui, il restante 90% è 
# rappresentato dalle restanti 6813 specie
# freq_sp_maggior <- frequenze_sp[1:26,]
# barplot(height = freq_sp_maggior$perc *100, names.arg = freq_sp_maggior$level, horiz = T, las = 1)

# Vedere come si distribuiscono le specie in Italia
ita_map <- map_data("italy")
ggplot() + 
  coord_fixed() + 
  geom_polygon(data = ita_map, aes(x = long, y = lat, group = group), colour = "gray50", fill = "gray50") +
  geom_point(data = dat_pulito, aes(x = decimalLongitude, y = decimalLatitude), colour = "darkred", size = 0.1) +
  theme_bw()
dev.off()



# RED LIST e Dataset Gbif ----
# Loading redlist
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

red_list <-read.csv("red_list_finale.csv")

# Verifico quali specie del dataset gbif sono presenti nella lista rossa
colnames(dat_pulito_sf)[1] <- "Name_species"
species_gbif_rl <- left_join(dat_pulito_sf, red_list, by = "Name_species")
species_gbif_rl <- species_gbif_rl[!is.na(species_gbif_rl$redlistCategory),]
class(species_gbif_rl)
summary(species_gbif_rl)
unique(species_gbif_rl$Name_species)

st_write(species_gbif_rl, "species_gbif_rl.shp")
species_gbif_rl <- st_read("species_gbif_rl.shp")

## Catalogo Specie Gbif nella Red List ----
species_gbif_rl <- species_gbif_rl %>%
  filter(year >= 1985 | is.na(year)) # seleziono solo le specie dal 1985 ad OGGI ed elimino gli NA 
lista_sp_cat_rischio <- split(species_gbif_rl, species_gbif_rl$rdlstCt)

names(lista_sp_cat_rischio)

sp_CE <- lista_sp_cat_rischio[[1]]
sp_DD <- lista_sp_cat_rischio[[2]]
sp_E <- lista_sp_cat_rischio[[3]]
sp_LC <- lista_sp_cat_rischio[[4]]
sp_NT <- lista_sp_cat_rischio[[5]]
sp_V <- lista_sp_cat_rischio[[6]]

unique(sp_DD$Nm_spcs) # per contare il numero di specie per ciascuna Categoria di Rischio

# rm(sp_CE,sp_DD,sp_E,sp_LC,sp_NT,sp_V)

# Trasformazione dat_pulito in elemento sf
dat_pulito_sf <- st_as_sf(dat_pulito, coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +datum=WGS84")
# verifico che i sistemi di riferimento dei due oggetti_sf siano gli stessi
st_crs(dat_pulito_sf) == st_crs(species_gbif_rl)

# Visualizzo le occorrenze gbif_rl 
ggplot() +
  coord_fixed() + 
  geom_polygon(data = ita_map, aes(x = long, y = lat, group = group), colour = "gray50", fill = "gray50", alpha = 0.5) +
  geom_sf(data = species_gbif_rl, colour = "darkred", size = 0.01, alpha = 0.5)+
  theme_minimal()
dev.off()



# VI Elenco Ufficiale Aree Naturali Protette (EUAP) ----
# Carico file .shp aree protette
siti_protet <- st_read("C:/Project_tirocinio/data/Siti_protetti_EUAP.shp")

# Correggo la geometria dei poligoni delle aree protette
siti_protet <- st_make_valid(siti_protet)

# st_crs(dat_pulito_sf) == st_crs(siti_protet) #verifca dei sistemi di riferimento

unique(siti_protet$tipo)
class(siti_protet$tipo)

# Elimino le Aree marine protette e Parchi Naturali Sommersi
siti_protet <- siti_protet %>%
  filter(tipo != "MAR" & tipo != "GAPN")

# Conversione della colonna tipo in un fattore ordinato
livelli_tipo <- c("RNS", "PNZ", "PNR", "RNR", "AANP")
siti_protet$tipo <- factor(siti_protet$tipo, levels = livelli_tipo)

# mapview(siti_protet) # open street map

# Plot delle aree protette suddivise per tipologia
ggplot() +
  coord_fixed() + 
  geom_polygon(data = ita_map, aes(x = long, y = lat, group = group), colour = "gray50", fill = "grey70", alpha = 0.5) +
  geom_sf(data = siti_protet, aes(fill = tipo), colour = NA)+
  scale_fill_manual(values = c("RNS" = "darkred", "PNZ" = "red", "PNR" = "orange", "RNR" = "yellow", "AANP" = "green"),
                    name = "Tipologia Area Protetta",
                    labels = c("Riserve Naturali Statali (RNS)", "Parchi Naturali Nazionali (PNZ)", 
                               "Parchi Naturali Regionali (PNR)", "Riserve Naturali Regionali (RNR)", 
                               "Altre Aree Naturali Protette (AANP)")) +
  theme_minimal()
dev.off()

## Numero totale aree e superficie totale per ogni tipologia di area protetta ----
siti_protet$superficie <- as.numeric(siti_protet$superficie)

siti_protet_type <- siti_protet %>% 
  group_by(tipo) %>%
  summarise(num = n(),
            superficie_tot = sum(superficie))%>% 
  st_drop_geometry()


## Contare le occorrenze DENTRO le aree protette----

species_gbif_rl$index_sp <- 1:nrow(species_gbif_rl)
sp_in_aree_pr <- st_within(species_gbif_rl, siti_protet)
sp_in_aree_pr_df <- as.data.frame(sp_in_aree_pr)
sp_in_aree_pr_df <- subset(sp_in_aree_pr_df, !duplicated(row.id))

colnames(sp_in_aree_pr_df)[1] <- "index_sp"
sp_in_aree_pr_df$col.id <- siti_protet$nome_gazze[sp_in_aree_pr_df$col.id]
colnames(sp_in_aree_pr_df)[2] <- "protected_area"
# species_cat_risk$index_sp <- as.numeric(species_cat_risk$index_sp)
sp_in_aree_pr_df <- left_join(sp_in_aree_pr_df,species_gbif_rl, by = "index_sp")
sp_in_aree_pr_df <- sp_in_aree_pr_df[,-1]
colnames(sp_in_aree_pr_df)[5] <- "occurence_gbifID"
sp_in_aree_pr_df <- sp_in_aree_pr_df[, c(5, 1:4, 6:ncol(sp_in_aree_pr_df))]
sp_in_aree_pr_df <- sp_in_aree_pr_df[, c(1, 3, 2, 4:ncol(sp_in_aree_pr_df))]
st_write(sp_in_aree_pr_df, "sp_in_aree_pr_df.shp")

# numero di specie 
unique(sp_in_aree_pr_df$Nm_spcs)

# specie suddivise per categoria di rischio, e sommo numero di specie e occorrenze
sp_in_aree_pr_cat_risk <- sp_in_aree_pr_df %>%
  group_by(rdlstCt) %>%
  summarise(num_species = n_distinct(Nm_spcs),
            num_occorrenze = n())

# suddivisione per categorie di rischio
list_sp_in_cat_risk <- split(sp_in_aree_pr_df, sp_in_aree_pr_df$rdlstCt)
sp_in_DD <- list_sp_in_cat_risk[[2]]
sp_in_LC <- list_sp_in_cat_risk[[4]]
sp_in_NT <- list_sp_in_cat_risk[[5]]
sp_in_V <- list_sp_in_cat_risk[[6]]
sp_in_E <- list_sp_in_cat_risk[[3]]
sp_in_CE <- list_sp_in_cat_risk[[1]]


## Contare le occorrenze FUORI dalle aree protette ----
sp_out_aree_pr_df <- as.data.frame(sp_in_aree_pr)
colnames(sp_out_aree_pr_df)[1] <- "index_sp"
sp_out_aree_pr_df$col.id <- siti_protet$nome_gazze[sp_out_aree_pr_df$col.id]
colnames(sp_out_aree_pr_df)[2] <- "protected_area"
sp_out <- pull(sp_out_aree_pr_df, index_sp)
sp_out_aree_pr_df <- species_gbif_rl %>% 
  filter(!(index_sp %in% sp_out))
sp_out_aree_pr_df <- sp_out_aree_pr_df[ , -ncol(sp_out_aree_pr_df)]
st_write(sp_out_aree_pr_df, "sp_out_aree_pr_df.shp")
# nrow(sp_in_aree_pr_df) + nrow(sp_out_aree_pr_df)  check

# numero di specie 
num_species_out <- unique(sp_out_aree_pr_df$Nm_spcs)

# specie suddivise per categoria di rischio, e sommo numero di specie e occorrenze
sp_out_aree_pr_cat_risk <- sp_out_aree_pr_df %>%
  group_by(rdlstCt) %>%
  summarise(num_species = n_distinct(Nm_spcs),
            num_occorrenze = n()) %>% 
  st_drop_geometry()

# suddivisione per categorie di rischio
list_sp_out_cat_risk <- split(sp_out_aree_pr_df, sp_out_aree_pr_df$rdlstCt)
sp_in_DD <- list_sp_out_cat_risk[[2]]
sp_in_LC <- list_sp_out_cat_risk[[4]]
sp_in_NT <- list_sp_out_cat_risk[[5]]
sp_in_V <- list_sp_out_cat_risk[[6]]
sp_in_E <- list_sp_out_cat_risk[[3]]
sp_in_CE <- list_sp_out_cat_risk[[1]]


# plot delle occorrenze fuori dalle aree protette #nel caso lo metta devo costruirmi 
# i vari sotto dataset per le cat di rischio 
ita_map <- map_data("italy")
ggplot() +
  coord_fixed() + 
  geom_polygon(data = ita_map, aes(x = long, y = lat, group = group), colour = "gray50", fill = "gray50", alpha = 0.5) +
  geom_sf(data = siti_protet, color = "darkgreen", fill = "green1", alpha = 0.3) +
  geom_sf(data = sp_CE_out_areepr, aes(color = "Critically Endangered (CE)"), size = 0.04) +
  geom_sf(data = sp_E_out_areepr, aes(color = "Endangered (EN)"), size = 0.04) +
  geom_sf(data = sp_V_out_areepr, aes(color = "Vulnerable (VU)"), size = 0.04) +
  geom_sf(data = sp_NT_out_areepr, aes(color = "Near Threatened (NT)"), size = 0.04) +
  scale_color_manual(name = "RedList Category",
                     values = c("darkred", "red", "orange2","yellow"),
                     labels = c("Critically Endangered (CE)", "Endangered (EN)", "Vulnerable(VU)",
                                "Near Threatened (NT)")) +
  theme_minimal()
dev.off()
# per la presentazione in latex meglio esportare il plot come pdf

# per un po' di pulizia: 
# rm(sp_in_aree_pr, sp_in_aree_pr_cat_risk, sp_in_aree_pr_df)
# rm(sp_out_aree_pr_cat_risk, sp_out_aree_pr_df)


## Specie CONDIVISE ----
sp_condivise_out <- semi_join(sp_out_aree_pr_df, sp_in_aree_pr_df, by = "Nm_spcs")
sp_condivise_in <- semi_join(sp_in_aree_pr_df, sp_out_aree_pr_df, by = "Nm_spcs") 
# or
# sp_condivise <- inner_join(sp_out_aree_pr_df, sp_in_aree_pr_df, by = "Nm_spcs")
nrow(sp_condivise_in) + nrow(sp_condivise_out)

length(unique(sp_condivise_in$Nm_spcs)) 
length(unique(sp_condivise_out$Nm_spcs))
# length(unique(species_gbif_rl$Nm_spcs)) 

# Raggruppo per categorie specie condivise ("dentro e fuori") generico
sp_condivise_cat_risk <- sp_condivise_in %>%  # non cambia se uso sp_condivise_out
  group_by(rdlstCt) %>%
  summarise(num_species = n_distinct(Nm_spcs)) %>% 
  st_drop_geometry()


### Quante sono le specie che sono SOLO DENTRO alle aree protette ----
df_sp_only_in <- anti_join(sp_in_aree_pr_df, sp_out_aree_pr_df, by = "Nm_spcs")
unique(df_sp_only_in$Nm_spcs)

# elimino le occorrenze che hanno come livello di rischio DD o LC
# df_sp_only_in <- df_sp_only_in[!(df_sp_only_out$rdlstCt %in% c("Least Concern", "Data Deficient")), ]

# write.csv(df_sp_only_in,"df_sp_only_in.csv")
df_sp_only_in <- read.csv("df_sp_only_in.csv")

# raggruppo specie e occorrenze per categoria di rischio
sp_only_in_cat_risk <- df_sp_only_in %>%
  group_by(rdlstCt) %>%
  summarise(num_species = n_distinct(Nm_spcs),
            num_occorrenze = n())

# suddivisione per categorie di rischio
list_sp_only_in_cat_risk <-split(df_sp_only_in, df_sp_only_in$rdlstCt)
sp_only_in_DD <- list_sp_only_in_cat_risk[[2]]
sp_only_in_LC <- list_sp_only_in_cat_risk[[4]]
sp_only_in_NT <- list_sp_only_in_cat_risk[[5]]
sp_only_in_V <- list_sp_only_in_cat_risk[[6]]
sp_only_in_E <- list_sp_only_in_cat_risk[[3]]
sp_only_in_CE <- list_sp_only_in_cat_risk[[1]]

unique(sp_only_in_CE$Name_species)  # per contare il numero di specie per ciascuna categoria

# qualche controllo 
nrow(sp_in_aree_pr_df) == nrow(sp_condivise_in) + nrow(df_sp_only_in)
length(unique(sp_in_aree_pr_df$Nm_spcs)) == length(unique(sp_condivise_in$Nm_spcs)) + length(unique(df_sp_only_in$Nm_spcs))



### Quante sono le specie che sono SOLO FUORI dalle aree protette ----
df_sp_only_out <- anti_join(sp_out_aree_pr_df, sp_in_aree_pr_df, by = "Nm_spcs")
sp_only_out_euap <- unique(df_sp_only_out$Nm_spcs)

# elimino le occorrenze che hanno come livello di rischio DD o LC
# df_sp_only_out <- df_sp_only_out[!(df_sp_only_out$rdlstCt %in% c("Least Concern", "Data Deficient")), ]

write.csv(df_sp_only_out,"df_sp_only_out.csv")
df_sp_only_out <- read.csv("df_sp_only_out.csv")

# raggruppo specie e occorrenze per categoria di rischio
sp_only_out_cat_risk <- df_sp_only_out %>%
  group_by(rdlstCt) %>%
  summarise(num_species = n_distinct(Nm_spcs),
            num_occorrenze = n()) %>% 
  st_drop_geometry()

# suddivisione per categorie di rischio
list_sp_only_out_cat_risk <-split(df_sp_only_out, df_sp_only_out$rdlstCt)
# sp_only_out_DD
# sp_only_out_LC
sp_only_out_NT <- list_sp_only_out_cat_risk[[5]]
sp_only_out_V <- list_sp_only_out_cat_risk[[6]]
sp_only_out_E <- list_sp_only_out_cat_risk[[3]]
sp_only_out_CE <- list_sp_only_out_cat_risk[[1]]

unique(sp_only_out_CE$Nm_spcs)

# qualche controllo
nrow(sp_out_aree_pr_df) == nrow(sp_condivise_out) + nrow(df_sp_only_out)
length(unique(sp_out_aree_pr_df$Nm_spcs)) == length(unique(sp_condivise_out$Nm_spcs)) + length(unique(df_sp_only_out$Nm_spcs))


#### PLOT delle OCCORRENZE appartanenti alle specie che sono SOLAMENTE FUORI dalle aree protette ----
# Non considero le sp DD e LC siccome non sono troppo rilevanti per quest'analisi
df_sp_only_out_most_risk <- df_sp_only_out %>% 
  filter(!(rdlstCt %in% c("Data Deficient", "Least Concern")))


ita_map <- map_data("italy")
df_sp_only_out_most_risk$rdlstCt <- factor(df_sp_only_out_most_risk$rdlstCt, 
                                           levels = c("Critically Endangered", "Endangered", 
                                                      "Vulnerable", "Near Threatened"))
ggplot() +
  #coord_fixed() + 
  geom_polygon(data = ita_map, aes(x = long, y = lat, group = group), colour = "gray50", fill = "gray70", alpha = 0.5) +
  geom_sf(data = siti_protet, color = "darkgreen", fill = "green3", alpha = 0.3, size = 0.1) +
  geom_sf(data = df_sp_only_out_most_risk, aes(color = rdlstCt), size = 0.6)+
  scale_color_manual(name = "RedList Category",
                     values = c("Critically Endangered" = "darkred", "Endangered" = "red", 
                                "Vulnerable" = "orange2", "Near Threatened" = "yellow"),
                     labels = c("Critically Endangered (CE)", "Endangered (EN)", "Vulnerable(VU)",
                                "Near Threatened (NT)"),
                     guide = guide_legend(override.aes = list(size = 1)) # cambiare dimensione legenda
  ) + 
  theme_minimal()
dev.off()

# per pulire tutto l'environment 
rm(list = ls())


## QUALI SP SONO CONTENUTE IN MENO AREE PROTETTE? ----
### Sp DENTRO ---- 
sp_in_low_protect <- sp_in_aree_pr_df %>%
  group_by(Nm_spcs) %>%
  summarise(num_areepr = n_distinct(protected_area),
            num_occorrenze = n()) %>% 
  arrange(num_areepr)

### Sp SOLO DENTRO ----
sp_only_in_low_protect <- df_sp_only_in %>%
  group_by(Nm_spcs) %>%
  summarise(num_areepr = n_distinct(protected_area),
            num_occorrenze = n()) %>% 
  arrange(num_areepr)



# ANALISI PER OGNI TIPOLOGIA DI AREA PROTETTA ----
lista_aree_pr <- split(siti_protet, siti_protet$tipo)

lista_aree_pr[[2]] <- lista_aree_pr[[2]] %>% 
  filter(nome_gazze != "Parco Nazionale dell'Arcipelago di La Maddalena")
pnz <- lista_aree_pr[[2]]

species_gbif_rl$index_sp <- 1:nrow(species_gbif_rl)

# creare una funzione per analizzare le tipologie di aree protette

function_type_areepr <- function(tipo_areepr){
  sp_in <- st_within(species_gbif_rl, tipo_areepr)
  sp_in <- as.data.frame(sp_in)
  colnames(sp_in)[1] <- "index_sp"
  colnames(sp_in)[2] <- "protected_area"
  sp_in$protected_area <- tipo_areepr$nome_gazze[sp_in$protected_area]
  sp_in <- left_join(sp_in, species_gbif_rl, by = "index_sp")
  sp_in <- sp_in[,-1]
  colnames(sp_in)[5] <- "occurence_gbifID"
  sp_in <- sp_in[, c(5, 1:4, 6:ncol(sp_in))]
  sp_in <- sp_in[, c(1, 3, 2, 4:ncol(sp_in))]
  return(sp_in)
}

results_typearee_sp_in <- lapply(lista_aree_pr, function_type_areepr)

# dataset contenente occorrenze e quindi specie per ogni tipologia di area protetta
sp_in_rns <- results_typearee_sp_in[[1]]
sp_in_pnz <- results_typearee_sp_in[[2]]
sp_in_pnr <- results_typearee_sp_in[[3]]
sp_in_rnr <- results_typearee_sp_in[[4]]
sp_in_aanp <- results_typearee_sp_in[[5]]

# function unique per trovare il totale di specie per ogni tipologia di area protetta!
unique(sp_in_aanp$Nm_spcs)


function_sp_in_cat_risk <- function(specie_in_areepr){
  df_sp_in_aggreg <- specie_in_areepr %>%
    group_by(rdlstCt) %>%
    summarise(num_species = n_distinct(Nm_spcs),
              num_occorrenze = n()) %>% 
    st_drop_geometry()
  return(df_sp_in_aggreg)
}


list_df_sp_in_aggr_catrisk <- lapply(results_typearee_sp_in, function_sp_in_cat_risk)

# per ogni tipologia di area protetta suddivisione delle specie per categoria di rischio
# conta numero specie e num occorrenze per ogni categoria
sp_in_rns_aggr_catrisk <- list_df_sp_in_aggr_catrisk[[1]]
sp_in_pnz_aggr_catrisk <- list_df_sp_in_aggr_catrisk[[2]]
sp_in_pnr_aggr_catrisk <- list_df_sp_in_aggr_catrisk[[3]]
sp_in_rnr_aggr_catrisk <- list_df_sp_in_aggr_catrisk[[4]]
sp_in_aanp_aggr_catrisk <- list_df_sp_in_aggr_catrisk[[5]]


count_specie <- function(specie_in_areepr){
  count_specie_areepr <- specie_in_areepr %>% 
    group_by(protected_area) %>%
    summarise(num_species = n_distinct(Nm_spcs)) %>% 
    arrange(desc(num_species))
  return(count_specie_areepr)
}


# lista delle tipologie di area protetta e per ciascuna area del dataframe sono indicate il numero di specie contenente
list_df_sp_in_num_type <- lapply(results_typearee_sp_in,count_specie)

num_sp_in_rns <- list_df_sp_in_num_type[[1]]
num_sp_in_pnz <- list_df_sp_in_num_type[[2]]
num_sp_in_pnr <- list_df_sp_in_num_type[[3]]
num_sp_in_rnr <- list_df_sp_in_num_type[[4]]
num_sp_in_aanp <- list_df_sp_in_num_type[[5]]


library(tidyr)
prova_funct <- function(specie_in_areepr){
  species_count <- specie_in_areepr %>%
    group_by(protected_area) %>%
    summarise(total_species = n_distinct(Nm_spcs))
  
  # Conteggio delle specie per categoria di rischio
  risk_count <- specie_in_areepr %>%
    group_by(protected_area, rdlstCt) %>%
    summarise(risk_species = n_distinct(Nm_spcs)) #%>%
  # spread(key = rdlstCt, value = risk_species, fill = 0)
  
  # Uniamo i due DataFrame
  result <- left_join(risk_count, species_count, by = "protected_area")
  result <- result %>% 
    arrange(desc(total_species))
  return(result)
  
}

lista_prova <- lapply(results_typearee_sp_in, prova_funct)
rns_prova <- lista_prova[[1]]
pnr_prova <- lista_prova[[3]]
rnr_prova <- lista_prova[[4]]
aanp_prova <- lista_prova[[5]]
pnz_prova <- lista_prova[[2]]

# Parco Nazionale dell'Arcipelago di La Maddalena non ha all'interno
# specie presenti nella Red List


# ci sono molte aree protette che non contengono specie della lista Rossa...


########### per i barplot delle RNS, PNR, RNR e AANP potrei prendere le prime con più specie e fare i 2 barplot ###########

## qualche grafico barplot per Parchi Nazionali ----

write.csv(num_sp_in_pnz, "num_sp_in_pnz.csv")

# Crea un vettore con le sigle manualmente
sigle <- c("PN-AS", "PN-GOG", "PN-GP", "PN-VG", "PN-ST","PN-DB","PN-CT","PN-ATE",
           "PN-FMC","PN-MS","PN-GSML","PN-MAI","PN-ALM","PN-CIR","PN-VES","PN-CVD",
           "PN-GAR","PN-AM ","PN-ALL","PN-POL","PN-SIL","PN-ASP","PN-ALM","PN-AT")

#pnz$sigle <- sigle

# barplot num specie Parchi Nazionali

lista_aree_pr[[2]]$sigle <- sigle

colnames(pnz)[4] <- "protected_area"
num_sp_in_pnz <- merge(num_sp_in_pnz, pnz[, c("protected_area","sigle")], by = "protected_area", all.x = TRUE)

num_sp_in_pnz <- num_sp_in_pnz %>% 
  arrange(desc(num_species)) %>% 
  st_drop_geometry()

ggplot(num_sp_in_pnz, aes(x = reorder(sigle, -num_species), y = num_species)) +
  geom_bar(stat = "identity")+
  labs(x = "National Park", y = "Num species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotazione etichette sull'asse x per una migliore leggibilità
        panel.background = element_rect(fill = "lightgrey", colour = NA))
dev.off()

# barplot più specifico numero di specie suddiviso per cat di rischio, (capire come aggregare i dati di sp_pnz_df)
pnz_prova$rdlstCt <- factor(pnz_prova$rdlstCt,
                            levels = c("Critically Endangered", "Endangered", 
                                       "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"))

pnz_prova <- merge(pnz_prova, pnz[, c("protected_area","sigle")], by = "protected_area", all.x = TRUE)

ggplot(pnz_prova, aes(x = reorder(sigle, -total_species), y = risk_species, fill = rdlstCt)) +
  geom_bar(stat = "identity", width = 0.9) +
  labs(x = "Parco Nazionale", y = "Numero di specie", fill = "Categoria Lista Rossa") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "whitesmoke", colour = NA)) +
  scale_fill_manual(values = c("Least Concern" = "green4", "Near Threatened" = "yellow", 
                               "Vulnerable" = "orange", "Endangered" = "red", 
                               "Critically Endangered" = "darkred", "Data Deficient" = "darkgrey"))
dev.off()

 
### barplot num specie Riserve naturali statali Nazionali ----

# ggplot(num_sp_in_rns, aes(x = reorder(protected_area, -num_species), y = num_species)) +
#  geom_bar(stat = "identity")+
#  labs(x = "Riserve Naturali Statali", y = "Num species") +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotazione etichette sull'asse x per una migliore leggibilità
#        panel.background = element_rect(fill = "lightgrey", colour = NA))
# dev.off()

num_sp_in_rns <- num_sp_in_rns %>% 
  filter(num_species >= 10)

ggplot(num_sp_in_rns, aes(x = reorder(protected_area, -num_species), y = num_species)) +
  geom_bar(stat = "identity")+
  labs(x = "Riserve Naturali Statali", y = "Num species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotazione etichette sull'asse x per una migliore leggibilità
        panel.background = element_rect(fill = "lightgrey", colour = NA))
dev.off()

# barplot più specifico numero di specie suddiviso per cat di rischio

rns_prova$rdlstCt <- factor(rns_prova$rdlstCt,
                            levels = c("Critically Endangered", "Endangered", 
                                       "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"))

rns_prova <- rns_prova %>% 
  filter(total_species >= 10)

ggplot(rns_prova, aes(x = reorder(protected_area, -total_species), y = risk_species, fill = rdlstCt)) +
  geom_bar(stat = "identity") +
  labs(x = "Parco Nazionale", y = "Numero di specie", fill = "Categoria Lista Rossa") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "whitesmoke", colour = NA)) +
  scale_fill_manual(values = c("Least Concern" = "green4", "Near Threatened" = "yellow", 
                               "Vulnerable" = "orange", "Endangered" = "red", 
                               "Critically Endangered" = "darkred", "Data Deficient" = "darkgrey"))
dev.off()

### barplot num specie Parchi naturali Regionali ----

ggplot(num_sp_in_pnr, aes(x = reorder(protected_area, -num_species), y = num_species)) +
  geom_bar(stat = "identity")+
  labs(x = "Parchi Naturali Regionali", y = "Num species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotazione etichette sull'asse x per una migliore leggibilità
        panel.background = element_rect(fill = "lightgrey", colour = NA))
dev.off()

num_sp_in_pnr <- num_sp_in_pnr %>% 
  filter(num_species >= 50)

ggplot(num_sp_in_pnr, aes(x = reorder(protected_area, -num_species), y = num_species)) +
  geom_bar(stat = "identity")+
  labs(x = "Parchi Naturali Regionali", y = "Num species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotazione etichette sull'asse x per una migliore leggibilità
        panel.background = element_rect(fill = "lightgrey", colour = NA))
dev.off()

# barplot più specifico numero di specie suddiviso per cat di rischio 

pnr_prova$rdlstCt <- factor(pnr_prova$rdlstCt,
                            levels = c("Critically Endangered", "Endangered", 
                                       "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"))

pnr_prova <- pnr_prova %>% 
  filter(total_species >= 50)

ggplot(pnr_prova, aes(x = reorder(protected_area, -total_species), y = risk_species, fill = rdlstCt)) +
  geom_bar(stat = "identity") +
  labs(x = "Parco Regionale", y = "Numero di specie", fill = "Categoria Lista Rossa") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "whitesmoke", colour = NA)) +
  scale_fill_manual(values = c("Least Concern" = "green4", "Near Threatened" = "yellow", 
                               "Vulnerable" = "orange", "Endangered" = "red", 
                               "Critically Endangered" = "darkred", "Data Deficient" = "darkgrey"))
dev.off()


### barplot num specie Riserve naturali Regionali ----

ggplot(num_sp_in_rnr, aes(x = reorder(protected_area, -num_species), y = num_species)) +
  geom_bar(stat = "identity")+
  labs(x = "Riserve Naturali Regionali", y = "Num species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotazione etichette sull'asse x per una migliore leggibilità
        panel.background = element_rect(fill = "lightgrey", colour = NA))
dev.off()

num_sp_in_rnr <- num_sp_in_rnr %>% 
  filter(num_species >= 20)

ggplot(num_sp_in_rnr, aes(x = reorder(protected_area, -num_species), y = num_species)) +
  geom_bar(stat = "identity")+
  labs(x = "Riserve Naturali Regionali", y = "Num species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotazione etichette sull'asse x per una migliore leggibilità
        panel.background = element_rect(fill = "lightgrey", colour = NA))
dev.off()

# barplot più specifico numero di specie suddiviso per cat di rischio 

rnr_prova$rdlstCt <- factor(rnr_prova$rdlstCt,
                            levels = c("Critically Endangered", "Endangered", 
                                       "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"))

rnr_prova <- rnr_prova %>% 
  filter(total_species >= 20)

ggplot(rnr_prova, aes(x = reorder(protected_area, -total_species), y = risk_species, fill = rdlstCt)) +
  geom_bar(stat = "identity") +
  labs(x = "Riserva Naturale Regionale", y = "Numero di specie", fill = "Categoria Lista Rossa") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "whitesmoke", colour = NA)) +
  scale_fill_manual(values = c("Least Concern" = "green4", "Near Threatened" = "yellow", 
                               "Vulnerable" = "orange", "Endangered" = "red", 
                               "Critically Endangered" = "darkred", "Data Deficient" = "darkgrey"))
dev.off()

### barplot num specie Altree aree naturali Protette ----

ggplot(num_sp_in_aanp, aes(x = reorder(protected_area, -num_species), y = num_species)) +
  geom_bar(stat = "identity")+
  labs(x = "Altre Riserve Naturali", y = "Num species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotazione etichette sull'asse x per una migliore leggibilità
        panel.background = element_rect(fill = "lightgrey", colour = NA))
dev.off()

num_sp_in_aanp <- num_sp_in_aanp %>% 
  filter(num_species >= 10)

ggplot(num_sp_in_aanp, aes(x = reorder(protected_area, -num_species), y = num_species)) +
  geom_bar(stat = "identity")+
  labs(x = "Altre Riserve Naturali", y = "Num species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotazione etichette sull'asse x per una migliore leggibilità
        panel.background = element_rect(fill = "lightgrey", colour = NA))
dev.off()

# barplot più specifico numero di specie suddiviso per cat di rischio

aanp_prova$rdlstCt <- factor(aanp_prova$rdlstCt,
                             levels = c("Critically Endangered", "Endangered", 
                                        "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"))

aanp_prova <- aanp_prova %>% 
  filter(total_species >= 10)

ggplot(aanp_prova, aes(x = reorder(protected_area, -total_species), y = risk_species, fill = rdlstCt)) +
  geom_bar(stat = "identity") +
  labs(x = "Altre Riserve Naturali", y = "Numero di specie", fill = "Categoria Lista Rossa") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "whitesmoke", colour = NA)) +
  scale_fill_manual(values = c("Least Concern" = "green4", "Near Threatened" = "yellow", 
                               "Vulnerable" = "orange", "Endangered" = "red", 
                               "Critically Endangered" = "darkred", "Data Deficient" = "darkgrey"))
dev.off()


##### funzioni ----
# Identificare le funzioni 
functions <- ls()[sapply(ls(), function(x) is.function(get(x)))]

# Rimuovere le funzioni
rm(list = functions)



# Siti Rete Natura 00 ----
sf_use_s2(FALSE) #pacchetto sf S2 gestione delle geometrie sferiche 
rete_nat_00 <- st_read("C:/Project_tirocinio/data/Rete_Nat_2000.shp")
rete_nat_00 <- st_make_valid(rete_nat_00)
# write.csv(rete_nat_00, "rete_nat_00.csv")
# st_write(rete_nat_00,"rete_nat_00.shp")

st_crs(rete_nat_00) == st_crs(siti_protet)

summary(rete_nat_00)
unique(rete_nat_00$tipo_sito)

# ZPS = sito di tipo A, SIC-ZSC = siti di tipo B e SIC-ZSC coincidenti con ZPS = siti di tipo C
# ZPS zone di protezione speciale istituite dalla Direttiva Uccelli
# ZSC zone speciali di conservazione della Direttiva Habitat

unique(rete_nat_00$reg_biog)

# area totale 
rete_nat_00$hectares <- as.numeric(rete_nat_00$hectares)
sup_totale <- round(sum(rete_nat_00$hectares), 3)

rete_nat_aggr <- rete_nat_00 %>% 
  st_drop_geometry() %>% 
  group_by(tipo_sito) %>% 
  summarise(num_sit = n(), tot_area = round(sum(hectares), 3))

# grafico siti rete nat 2000 suddivisi per tipologia 
rete_nat_00$tipo_sito <- factor(rete_nat_00$tipo_sito, 
                                levels = c("A","B","C"))
ggplot() +
  coord_fixed() +
  geom_polygon(data = ita_map, aes(x = long, y = lat, group = group), colour = "gray50", fill = "gray70", alpha = 0.5) +
  geom_sf(data = rete_nat_00, aes(fill = tipo_sito), color = NA)+
  scale_fill_manual(values = c("A" = "red", "B" = "orange", "C" = "yellow"),
                    name = "Tipologia Sito Protetto",
                    labels = c("ZPS", "SIC-ZSC", "SIC-ZSC ≡ ZPS")) +
  theme_minimal()
dev.off()


lista_siti_rete_nat <- split(rete_nat_00, rete_nat_00$tipo_sito)

# per dataframe della "lista_siti_rete_nat" 
# mi riesco a calcolare quanti sono i siti per ogni reg biologica
#.... <- .... %>% 
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
siti_a <- rete_nat_00 %>% 
  filter(tipo_sito == "A")

siti_a_reg_bio <- siti_a %>% 
  st_drop_geometry() %>% 
  group_by(reg_biog) %>% 
  summarise(num_sit = n())

sp_in_a <- results_typesiti_sp_in[[1]]

sp_in_a <- sp_in_a %>%
  rename(denominazi = protected_area)

sp_in_a <- merge(sp_in_a, siti_a[,c("denominazi","reg_biog")], by = "denominazi", all.x = T)

unique(sp_in_a$Nm_spcs) # 578

a_num_sp_reg_bio <- sp_in_a %>% 
  st_drop_geometry() %>% 
  group_by(reg_biog) %>% 
  summarise(num_species = n_distinct(Nm_spcs))

a_num_sp_cat_risk <- sp_in_a %>% 
  st_drop_geometry() %>% 
  group_by(rdlstCt) %>% 
  summarise(num_species = n_distinct(Nm_spcs))


# Conteggio delle specie per area protetta
a_species_count <- sp_in_a %>%
  group_by(protected_area) %>%
  summarise(total_species = n_distinct(Nm_spcs)) %>% 
  arrange(desc(num_species))

# Conteggio delle specie per categoria di rischio
a_risk_count <- sp_in_a %>%
  group_by(protected_area, rdlstCt) %>%
  summarise(risk_species = n_distinct(Nm_spcs)) #%>%
# spread(key = rdlstCt, value = risk_species, fill = 0)

# Uniamo i due DataFrame
a_result <- left_join(a_risk_count, a_species_count, by = "protected_area")
a_result <- a_result %>% 
  arrange(desc(total_species))

# plot primi siti protetti con 50+ specie 

a_first_result <- a_result[1:58,]

a_first_result$rdlstCt <- factor(a_first_result$rdlstCt,
                                 levels = c("Critically Endangered", "Endangered", 
                                            "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"))

ggplot(a_first_result, aes(x = reorder(protected_area, -total_species), y = risk_species, fill = rdlstCt)) +
  geom_bar(stat = "identity") +
  labs(x = "Siti Protetti tipo A", y = "Numero di specie", fill = "Categoria Lista Rossa") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "whitesmoke", colour = NA)) +
  scale_fill_manual(values = c("Least Concern" = "green4", "Near Threatened" = "yellow", 
                               "Vulnerable" = "orange", "Endangered" = "red", 
                               "Critically Endangered" = "darkred", "Data Deficient" = "darkgrey"))
dev.off()



### siti b : ----
siti_b <- rete_nat_00 %>% 
  filter(tipo_sito == "B")

siti_b_reg_bio <- siti_b %>% 
  st_drop_geometry() %>% 
  group_by(reg_biog) %>% 
  summarise(num_sit = n())

sp_in_b <- results_typesiti_sp_in[[2]]

sp_in_b <- sp_in_b %>%
  rename(denominazi = protected_area)

sp_in_b <- merge(sp_in_b, siti_b[,c("denominazi","reg_biog")], by = "denominazi", all.x = T)

unique(sp_in_b$Nm_spcs) # 650

b_num_sp_reg_bio <- sp_in_b %>% 
  st_drop_geometry() %>% 
  group_by(reg_biog) %>% 
  summarise(num_species = n_distinct(Nm_spcs))

b_num_sp_cat_risk <- sp_in_b %>% 
  st_drop_geometry() %>% 
  group_by(rdlstCt) %>% 
  summarise(num_species = n_distinct(Nm_spcs))

# Conteggio delle specie per area protetta
b_species_count <- sp_in_b %>%
  group_by(protected_area) %>%
  summarise(total_species = n_distinct(Nm_spcs)) %>% 
  arrange(desc(total_species))

# Conteggio delle specie per categoria di rischio
b_risk_count <- sp_in_b %>%
  group_by(protected_area, rdlstCt) %>%
  summarise(risk_species = n_distinct(Nm_spcs)) #%>%
# spread(key = rdlstCt, value = risk_species, fill = 0)

# Uniamo i due DataFrame
b_result <- left_join(b_risk_count, b_species_count, by = "protected_area")
b_result <- b_result %>% 
  arrange(desc(total_species))

# plot primi siti protetti con 50+ specie 

b_first_result <- b_result[1:59,]

b_first_result$rdlstCt <- factor(b_first_result$rdlstCt,
                                 levels = c("Critically Endangered", "Endangered", 
                                            "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"))

ggplot(b_first_result, aes(x = reorder(protected_area, -total_species), y = risk_species, fill = rdlstCt)) +
  geom_bar(stat = "identity") +
  labs(x = "Siti Protetti tipo A", y = "Numero di specie", fill = "Categoria Lista Rossa") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "whitesmoke", colour = NA)) +
  scale_fill_manual(values = c("Least Concern" = "green4", "Near Threatened" = "yellow", 
                               "Vulnerable" = "orange", "Endangered" = "red", 
                               "Critically Endangered" = "darkred", "Data Deficient" = "darkgrey"))
dev.off()


### siti c : ----
siti_c <- rete_nat_00 %>% 
  filter(tipo_sito == "C")

siti_c_reg_bio <- siti_c %>% 
  st_drop_geometry() %>% 
  group_by(reg_biog) %>% 
  summarise(num_sit = n())

sp_in_c <- results_typesiti_sp_in[[3]]

sp_in_c <- sp_in_c %>%
  rename(denominazi = protected_area)

sp_in_c <- merge(sp_in_c, siti_c[,c("denominazi","reg_biog")], by = "denominazi", all.x = T)

unique(sp_in_c$Nm_spcs) # 492

c_num_sp_reg_bio <- sp_in_c %>% 
  st_drop_geometry() %>% 
  group_by(reg_biog) %>% 
  summarise(num_species = n_distinct(Nm_spcs))

c_num_sp_cat_risk <- sp_in_c %>% 
  st_drop_geometry() %>% 
  group_by(rdlstCt) %>% 
  summarise(num_species = n_distinct(Nm_spcs))

# Conteggio delle specie per area protetta
c_species_count <- sp_in_c %>%
  group_by(protected_area) %>%
  summarise(total_species = n_distinct(Nm_spcs)) %>% 
  arrange(desc(total_species))

# Conteggio delle specie per categoria di rischio
c_risk_count <- sp_in_c %>%
  group_by(protected_area, rdlstCt) %>%
  summarise(risk_species = n_distinct(Nm_spcs)) #%>%
# spread(key = rdlstCt, value = risk_species, fill = 0)

# Uniamo i due DataFrame
c_result <- left_join(c_risk_count, c_species_count, by = "protected_area")
c_result <- c_result %>% 
  arrange(desc(total_species))

# plot primi siti protetti con 50+ specie 

c_first_result <- c_result[1:26,]

c_first_result$rdlstCt <- factor(c_first_result$rdlstCt,
                                 levels = c("Critically Endangered", "Endangered", 
                                            "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient"))

ggplot(c_first_result, aes(x = reorder(protected_area, -total_species), y = risk_species, fill = rdlstCt)) +
  geom_bar(stat = "identity") +
  labs(x = "Siti Protetti tipo A", y = "Numero di specie", fill = "Categoria Lista Rossa") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "whitesmoke", colour = NA)) +
  scale_fill_manual(values = c("Least Concern" = "green4", "Near Threatened" = "yellow", 
                               "Vulnerable" = "orange", "Endangered" = "red", 
                               "Critically Endangered" = "darkred", "Data Deficient" = "darkgrey"))
dev.off()



# sp_only out sono dentro i siti a?
sp_only_out_in_a <- sp_only_out_euap %in% sp_in_a$Nm_spcs

sp_only_out_in_a <- sp_only_out_euap[sp_only_out_in_a]

# sp_only out sono dentro i siti b?
sp_only_out_in_b <- sp_only_out_euap %in% sp_in_b$Nm_spcs

sp_only_out_in_b <- sp_only_out_euap[sp_only_out_in_b]

# sp_only out sono dentro i siti c?
sp_only_out_in_c <- sp_only_out_euap %in% sp_in_c$Nm_spcs

sp_only_out_in_c <- sp_only_out_euap[sp_only_out_in_c]

sp_only_out_in_nat_00 <- unique(c(sp_only_out_in_a, sp_only_out_in_b, sp_only_out_in_c))

# delle 188 sp only out euap, 85 sono all'interno dei siti nat 00!
# Quali sono?

sp_only_out_nat_00 <- setdiff(sp_only_out_euap, sp_only_out_in_nat_00)

species_gbif_rl_total_out <- species_gbif_rl %>%
  filter(Nm_spcs %in% sp_only_out_nat_00)

class(species_gbif_rl_total_out)

species_gbif_rl_total_out <- species_gbif_rl_total_out %>%
  distinct()

unique(species_gbif_rl_total_out$Nm_spcs)




## plot occorrenze delle specie sia fuori dalle aree protette e sia dai siti rete nat 00 ----

ita_map <- map_data("italy")

species_gbif_rl_total_out$rdlstCt <- factor(species_gbif_rl_total_out$rdlstCt, 
                                            levels = c("Critically Endangered", "Endangered", 
                                                       "Vulnerable", "Near Threatened", 
                                                       "Least Concern", "Data Deficient"))



ggplot() +
  #coord_fixed() + 
  geom_polygon(data = ita_map, aes(x = long, y = lat, group = group), colour = "gray50", fill = "gray70", alpha = 0.5) +
  geom_sf(data = siti_protet, color = NA, fill = "darkgreen", alpha = 0.4, size = 0.1) +
  geom_sf(data = rete_nat_00, color = NA, fill = "chartreuse2", alpha = 0.4, size = 0.1)+
  geom_sf(data = species_gbif_rl_total_out, aes(color = rdlstCt), size = 0.8)+
  scale_color_manual(name = "RedList Category",
                     values = c("Critically Endangered" = "darkred", "Endangered" = "red", 
                                "Vulnerable" = "orange2", "Near Threatened" = "yellow",
                                "Least Concern" = "burlywood2", "Data Deficient" = "cyan2"),
                     labels = c("Critically Endangered (CE)", "Endangered (EN)", "Vulnerable(VU)",
                                "Near Threatened (NT)", "Least Concern (LC)", "Data Deficient (DD)"),
                     guide = guide_legend(override.aes = list(size = 1)) # cambiare dimensione legenda
  ) + 
  theme_minimal()
dev.off()
