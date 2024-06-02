
# Spatial Analysis and Statistical Modeling with R and spmodel ------------


# Spatial Linear Models in spmodel ----------------------------------------

# Moss Data ---------------------------------------------------------------
library(spdata.sfs24)
library(spmodel)
library(ggplot2)

data('moss')

ggplot(moss, aes(color = log_Zn)) +
  geom_sf(size = 2) +
  scale_color_viridis_c() +
  scale_x_continuous(breaks = seq(-163, -164, length.out = 2)) +
  theme_gray(base_size = 18)

spmod <- splm(formula = log_Zn ~ log_dist2road, data = moss,
              spcov_type = "exponential")
summary(spmod)

varcomp(spmod)

tidy(spmod)

tidy(spmod, effects = "spcov")

glance(spmod)

augment(spmod)

plot(spmod, which = 4)

none <- splm(formula = log_Zn ~ log_dist2road, data = moss,
             spcov_type = "none")

glances(spmod, none)

loocv(spmod)

loocv(none)

# Moose Data --------------------------------------------------------------

data('moose')

ggplot(data = moose, aes(colour = count)) +
  geom_sf(size = 2) +
  scale_colour_viridis_c(limits = c(0, 40)) +
  theme_minimal(base_size = 18)

ggplot(data = moose_preds) +
  geom_sf(size = 2) +
  theme_minimal(base_size = 18)

moosemod <- splm(count ~ elev * strat, data = moose,
                 spcov_type = "spherical")
tidy(moosemod)

predict(moosemod, newdata = moose_preds)

moose_aug <- augment(moosemod, newdata = moose_preds)
moose_aug

ggplot(data = moose, aes(colour = count)) +
  geom_sf(alpha = 0.4) +
  geom_sf(data = moose_aug, aes(colour = .fitted)) +
  scale_colour_viridis_c(limits = c(0, 40)) +
  theme_minimal()

# Spatial Generalized Linear Models in spmodel ----------------------------

poismod <- spglm(count ~ elev * strat, data = moose,
                 family = poisson, spcov_type = "spherical")

summary(poismod)

predict(poismod, newdata = moose_preds)

augment(poismod, newdata = moose_preds)


# GIS in R ----------------------------------------------------------------


# Points, lines, and polygons ---------------------------------------------

library(tidyverse)
library(ggplot2)

id <- c(1:5)
cities <- c('Ashland','Corvallis','Bend','Portland','Newport')
longitude <- c(-122.699, -123.275, -121.313, -122.670, -124.054)
latitude <- c(42.189, 44.57, 44.061, 45.523, 44.652)
population <- c(20062, 50297, 61362, 537557, 9603)

oregon_cities <- data.frame(id, cities, longitude, latitude, population)

ggplot(
  data = oregon_cities,
  aes(x = longitude,
      y = latitude,
      size = population,
      label = cities)) +
  geom_point() +
  geom_text(hjust = 1, vjust = 1) +
  theme_bw()

library(sf)
ls("package:sf")

oregon_cities <- oregon_cities %>%
  st_as_sf(coords = c('longitude', 'latitude'), crs = 4269)

print(oregon_cities)

st_crs(oregon_cities) == st_crs(5070)

oregon_cities <-
  oregon_cities %>%
  st_transform(crs = 5070)

cities_buffer <-
  oregon_cities %>%
  st_buffer(100000)

ggplot(data = cities_buffer) +
  geom_sf(aes(fill = cities), alpha = 0.5) +
  geom_sf(data = st_centroid(oregon_cities)) +
  theme_bw()

cities_buffer <- cities_buffer %>%
  st_intersection() %>%
  mutate(area = st_area(.) %>%
           units::drop_units(),
         id = as.factor(1:nrow(.)))

ggplot(data = cities_buffer) +
  geom_sf(aes(fill = id), alpha = 0.5) +
  theme_bw()


# Raster Data -------------------------------------------------------------

library(terra)
ls("package:terra")

r <- rast(ncol=10, nrow = 10)

r[] <- runif(n=ncell(r))

r

plot(r)

r[12]

r[2, 2]

r2 <- r * 50
r3 <- sqrt(r * 5)

s <- c(r, r2, r3)
names(s) <- c('r1', 'r2', 'r3')

plot(s)

# Working with Real Data --------------------------------------------------

library(FedData)

# Select just Corvallis and calculate a 10,000-m buffer
corvallis <-
  oregon_cities %>%
  filter(cities == 'Corvallis') %>%
  st_buffer(10000)

# Download national elevation data (ned)
ned <- FedData::get_ned(
  template = corvallis,
  label = "corvallis")

ned <- terra::project(ned,
                      'epsg:5070',
                      method = 'bilinear')

# zonal function in terra to calculate zonal statistics
terra::zonal(ned,

             # Need to convert corvallis `sf` object to terra vector
             terra::vect(corvallis),

             # Metric to be calculated
             mean, na.rm = T)


# Solution to excercise 1 ----------------------------------------------------------------

library(maps)
data('us.cities')

my_city <- us.cities %>%
  filter(name == 'Idaho Falls ID') %>%
  st_as_sf(coords = c('long', 'lat'), crs = 4269) %>%
  st_transform(crs = 5070) %>%
  st_buffer(10000)

ned <- FedData::get_ned(
  template = my_city,
  label = "Idaho Falls")

ned <- terra::project(ned,
                      'epsg:5070',
                      method = 'bilinear')

terra::zonal(ned,
             terra::vect(my_city),
             mean, na.rm = T)

# Watershed Delineation ---------------------------------------------------

# USGS StreamStats

streamstats_ws = function(state, longitude, latitude){
  p1 = 'https://streamstats.usgs.gov/streamstatsservices/watershed.geojson?rcode='
  p2 = '&xlocation='
  p3 = '&ylocation='
  p4 = '&crs=4269&includeparameters=false&includeflowtypes=false&includefeatures=true&simplify=true'
  query <-  paste0(p1, state, p2, toString(longitude), p3, toString(latitude), p4)
  mydata <- jsonlite::fromJSON(query, simplifyVector = FALSE, simplifyDataFrame = FALSE)
  poly_geojsonsting <- jsonlite::toJSON(mydata$featurecollection[[2]]$feature, auto_unbox = TRUE)
  poly <- geojsonio::geojson_sf(poly_geojsonsting)
  poly
}

# Define location for delineation (Calapooia Watershed)
state <- 'OR'
latitude <- 44.62892
longitude <- -123.13113

# Delineate watershed
cal_ws <- streamstats_ws('OR', longitude, latitude) %>%
  st_transform(crs = 5070)

library(mapview)
mapview::mapviewOptions(fgb=FALSE)
mapview(cal_ws, alpha.regions = .08)

# nhdplusTools & NLDI service

library(nhdplusTools)

# Simple feature option to generate point without any other attributes
cal_pt <- st_sfc(st_point(c(longitude, latitude)), crs = 4269)

# Identify the network location (NHDPlus common ID or COMID)
start_comid <- nhdplusTools::discover_nhdplus_id(cal_pt)

# Combine info into list (required by NLDI basin function)
ws_source <- list(featureSource = "comid", featureID = start_comid)

cal_ws2 <- nhdplusTools::get_nldi_basin(nldi_feature = ws_source)

mapview(cal_ws, alpha.regions = .08) +
  mapview(cal_ws2, alpha.regions = .08, col.regions = 'red')


# Solution to excercise 2 ----------------------------------------------------------------

# Logan River Watershed
latitude <- 41.707
longitude <- -111.855

# Define the lat/lon
start_point <- st_sfc(st_point(c(longitude, latitude)), crs = 4269)
# Find COMID of this point
start_comid <- nhdplusTools::discover_nhdplus_id(start_point)
# Create a list object that defines the feature source and starting COMID
ws_source <- list(featureSource = "comid", featureID = start_comid)
# Delineate basin
logan_ws <- nhdplusTools::get_nldi_basin(nldi_feature = ws_source)

mapview::mapview(logan_ws)

# StreamCatTools ----------------------------------------------------------

# StreamCat

nlcd <- get_nlcd(
  template = cal_ws2,
  year = 2019,
  label = 'Calapooia') %>%
  terra::project('epsg:5070', method = 'near')

terra::extract(nlcd,
               terra::vect(cal_ws2)) %>%
  group_by(Class) %>%
  summarise(count = n()) %>%
  mutate(percent = (count / sum(count)) * 100) %>%
  filter(Class == 'Cultivated Crops') %>%
  print()

library(StreamCatTools)

comid <- '23763521'
sc_get_data(comid = comid,
            metric = 'PctCrop2019',
            aoi = 'watershed')

sc_get_data(comid = comid,
            metric = 'PctCrop2019',
            aoi = 'catchment,watershed')

sc_get_data(comid = comid,
            metric = 'PctCrop2019',
            aoi = 'catchment,riparian_catchment,watershed,riparian_watershed') %>%
  as.data.frame()

iowa_crop <- sc_get_data(state = 'IA',
                         metric = 'PctCrop2019',
                         aoi = 'watershed')

ggplot() +
  geom_histogram(data = iowa_crop,
                 aes(x = PCTCROP2019WS)) +
  theme_bw()


# LakeCAt

library(nhdplusTools)

# Pelican Lake, WI
latitude <- 45.502840
longitude <- -89.198694

pelican_pt <- data.frame(site_name = 'Pelican Lake',
                         latitude = latitude,
                         longitude = longitude) %>%
  st_as_sf(coords = c('longitude', 'latitude'), crs = 4326)

pelican_lake <- nhdplusTools::get_waterbodies(pelican_pt)

comid <- pelican_lake %>%
  pull(comid)

lc_get_data(metric = 'elev, cao, sand, om, pctdecid2019',
            aoi = 'watershed',
            comid = comid)


# Solution to excercise 3 ----------------------------------------------------------------

library(StreamCatTools)

comids <- "9028333,9025609,9025611,9025635"

lc_get_data(metric = 'pctdecid2019, pctconif2019',
            aoi = 'watershed',
            comid = comids)


# Lake Conductivity -------------------------------------------------------

library(tigris)
library(StreamCatTools)
library(spmodel)
library(data.table)
library(spdata.sfs24)
#library(data.table)

data('cond_nla_data')

# Read in states to give some context
states <- tigris::states(cb = TRUE, progress_bar = FALSE)  %>%
  filter(!STUSPS %in% c('HI', 'PR', 'AK', 'MP', 'GU', 'AS', 'VI'))  %>%
  st_transform(crs = 5070)

ggplot() +
  geom_sf(data = states,
          fill = NA) +
  geom_sf(data = lake_cond,
          aes(color = DSGN_CYCLE)) +
  scale_color_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a")) +
  theme_bw() +
  theme(legend.position="bottom")

MN <- states %>%
  filter(STUSPS == 'MN')

cond_mn <- lake_cond %>%
  st_filter(MN) %>%
  rename(year = DSGN_CYCLE)

ggplot() +
  geom_sf(data = MN,
          fill = NA) +
  geom_sf(data = cond_mn,
          aes(color = log(COND_RESULT))) +
  scale_color_distiller(palette = 'YlOrRd', direction = 1) +
  theme_bw() +
  theme(legend.position="bottom")


comids <- cond_mn$COMID

mn_lakecat <- lc_get_data(comid = comids,
                          metric = 'Tmean8110, Precip8110, S') %>%
  select(COMID, TMEAN8110CAT, PRECIP8110WS, SWS)

mn_lakecat


crop <-

  # Grab LakeCat crop data
  lc_get_data(comid = comids,
              aoi = 'watershed',
              metric = 'pctcrop2006, pctcrop2011, pctcrop2016') %>%

  # Remove watershed area from data
  select(-WSAREASQKM) %>%

  # Pivot table to long to create "year" column
  pivot_longer(!COMID, names_to = 'tmpcol', values_to = 'PCTCROPWS') %>%

  # Remove PCTCROP and WS to make "year" column
  mutate(year = as.integer(
    str_replace_all(tmpcol, 'PCTCROP|WS', ''))) %>%

  # Add 1 to each year to match NLA years
  mutate(year = factor(year + 1)) %>%

  # Remove the tmp column
  select(-tmpcol)


model_data <- cond_mn %>%
  left_join(mn_lakecat, join_by(COMID)) %>%
  left_join(crop, join_by(COMID, year))

formula <-
  log(COND_RESULT) ~
  AREA_HA +
  year +
  TMEAN8110CAT +
  PRECIP8110WS +
  PCTCROPWS +
  SWS

cond_mod <- splm(formula = formula,
                 data = model_data,
                 spcov_type = 'none')

cond_spmod <- splm(formula = formula,
                   data = model_data,
                   spcov_type = 'exponential')

summary(cond_spmod)

glances(cond_mod, cond_spmod)

prd_mod <- spmodel::loocv(cond_mod, se.fit = TRUE, cv_predict = TRUE)

prd_spmod <- spmodel::loocv(cond_spmod, se.fit = TRUE, cv_predict = TRUE)

rbind(prd_mod %>% pluck('stats'),
      prd_spmod %>% pluck('stats'))

# Combine predictions with model data (spatial points)
model_data <-
  model_data %>%
  mutate(prd_cond = prd_spmod %>%
           pluck('cv_predict'),
         se_fit = prd_spmod %>%
           pluck('se.fit'))
# Plot predictions
ggplot() +
  geom_sf(data = MN,
          fill = NA) +
  geom_sf(data = model_data,
          aes(color = prd_cond)) +
  scale_color_distiller(palette = 'YlOrRd', direction = 1) +
  theme_bw() +
  theme(legend.position="bottom")

# Plot standard errors
ggplot() +
  geom_sf(data = MN,
          fill = NA) +
  geom_sf(data = model_data,
          aes(color = se_fit)) +
  scale_color_distiller(palette = 'Reds', direction = 1) +
  theme_bw() +
  theme(legend.position="bottom")

library(finsyncR)
library(pROC)

macros <- getInvertData(dataType = "occur",
                        taxonLevel = "Genus",
                        agency = "EPA",
                        lifestage = FALSE,
                        rarefy = TRUE,
                        rarefyCount = 300,
                        sharedTaxa = FALSE,
                        seed = 1,
                        boatableStreams = T)

# Flexible code so we could model another taxon
genus <- 'Argia'

taxon = macros %>%
  dplyr::select(SampleID,
                ProjectLabel,
                CollectionDate,
                Latitude_dd,
                Longitude_dd,
                all_of(genus))  %>%
  #filter(ProjectLabel != 'WSA') %>%
  mutate(CollectionDate = date(CollectionDate),
         presence =
           as.factor(pull(., genus)))  %>%
  st_as_sf(coords = c('Longitude_dd', 'Latitude_dd'), crs = 4269)  %>%
  st_transform(crs = 5070)


ggplot() +
  geom_sf(data = states, fill = NA) +
  geom_sf(data = taxon,
          aes(color = presence),
          size = 1.5,
          alpha = 0.65) +
  scale_color_manual(values=c("#d9d9d9", "#08519c")) +
  theme_bw() +
  theme(legend.position="bottom")

ggplot() +
  geom_sf(data = states, fill = NA) +
  geom_sf(data = taxon,
          aes(color = ProjectLabel),
          size = 1.5,
          alpha = 0.75) +
  scale_color_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")) +
  theme_bw() +
  theme(legend.position="bottom")

region <- states %>%
  filter(STUSPS %in% c('VT', 'NH', 'ME', 'NY', 'RI',
                       'MA', 'CT', 'NJ', 'PA', 'DE'))

# Use region as spatial filter (sf::st_filter()) for taxon of interest
taxon_rg <- taxon %>%
  st_filter(region) %>%
  filter(ProjectLabel %in% c('NRSA1314', 'NRSA1819')) %>%
  mutate(year = year(ymd(CollectionDate))) %>%
  select(SampleID:CollectionDate, presence:year)

taxon_rg %>%
  pull(presence) %>%
  table()

ggplot() +
  geom_sf(data = region, fill = NA) +
  geom_sf(data = taxon_rg,
          aes(color = presence)) +
  scale_color_manual(values=c("#d9d9d9", "#08519c")) +
  theme_bw() +
  theme(legend.position="bottom")

comids <- sc_get_comid(taxon_rg)

# data('comids')

comid_vect <-
  comids %>%
  str_split(',') %>%
  unlist() %>%
  as.integer()

taxon_rg <-
  taxon_rg %>%
  mutate(COMID = comid_vect)

sc <-
  sc_get_data(comid = comids,
              aoi = 'watershed',
              metric = 'bfi, precip8110, wetindex, elev',
              showAreaSqKm = TRUE)


library(prism)

# Get these years of PRISM
years <- c(2013, 2014, 2018, 2019)

# Set the PRISM directory (creates directory in not present)
prism_set_dl_dir("./data/prism_data", create = TRUE)

# Download monthly PRISM rasters (tmean)
get_prism_monthlys('tmean',
                   years = years,
                   mon = 7:8,
                   keepZip = FALSE)


# Create stack of downloaded PRISM rasters
tmn <- pd_stack((prism_archive_subset("tmean","monthly",
                                      years = years,
                                      mon = 7:8)))

# Extract tmean at sample points and massage data
tmn <- terra::extract(tmn,
                      # Transform taxon_rg to CRS of PRISM on the fly
                      taxon_rg %>%
                        st_transform(crs = st_crs(tmn))) %>%

  # Add COMIDs to extracted values
  data.frame(COMID = comid_vect, .) %>%

  # Remove front and back text from PRISM year/month in names
  rename_with( ~ stringr::str_replace_all(., 'PRISM_tmean_stable_4kmM3_|_bil', '')) %>%

  # Pivot to long table and calle column TMEANPRISMXXXXPT, XXXX indicates year
  pivot_longer(!COMID, names_to = 'year_month',
               values_to = 'TMEANPRISMXXXXPT') %>%

  # Create new column of year
  mutate(year = year(ym(year_month))) %>%

  # Average July and August temperatures
  summarise(TMEANPRISMXXXXPT = mean(TMEANPRISMXXXXPT, na.rm = TRUE),
            .by = c(COMID, year))


model_data <-
  taxon_rg %>%
  left_join(sc, join_by(COMID)) %>%
  left_join(tmn, join_by(COMID, year)) %>%
  drop_na()

formula <-
  presence ~
  I(log10(WSAREASQKM)) +
  ELEVWS +
  WETINDEXWS +
  BFIWS +
  PRECIP8110WS +
  TMEANPRISMXXXXPT

bin_mod <- spglm(formula = formula,
                 data = model_data,
                 family = 'binomial',
                 spcov_type = 'none')

bin_spmod <- spglm(formula = formula,
                   data = model_data,
                   family = 'binomial',
                   spcov_type = 'exponential')

summary(bin_spmod)

glances(bin_mod, bin_spmod)

# Function to convert from log odds to probability
to_prob <- function(x) exp(x)/(1+exp(x))

# loocv of non-spatial model
loocv_mod <- loocv(bin_mod, cv_predict = TRUE, se.fit = TRUE)

# loocv of spatial model
loocv_spmod <- loocv(bin_spmod, cv_predict = TRUE, se.fit = TRUE)

pROC::auc(model_data$presence, loocv_mod$cv_predict)

prob_spmod <- to_prob(loocv_spmod$cv_predict)
sefit_spmod <- loocv_spmod$se.fit

model_data <- model_data %>%
  mutate(prob_spmod = prob_spmod,
         sefit_spmod = sefit_spmod)

# Plot predicted values
ggplot() +
  geom_sf(data = region, fill = NA) +
  geom_sf(data = model_data,
          aes(color = prob_spmod)) +
  scale_color_viridis_b() +
  theme_bw() +
  theme(legend.position="bottom")

# Plot standard errors
ggplot() +
  geom_sf(data = region, fill = NA) +
  geom_sf(data = model_data,
          aes(color = sefit_spmod)) +
  scale_color_distiller(palette = 'Reds', direction = 1) +
  theme_bw() +
  theme(legend.position="bottom")

