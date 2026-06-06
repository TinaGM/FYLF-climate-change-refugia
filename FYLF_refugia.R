#Code for foothill yellow-legged frog modeling
#Written by Tina Mozelewski during University of Massachusetts-Amherst + Northeast Climate Adaptation Science Center postdoc (2021-2022)
#Use species distribution modeling to map habitat suitability at 2020, 2040, and 2080
#Places on the landscape that are suitable habitat now AND suitable habitat in the future are CC refugia

library(raster)
library(rgdal)
library(dismo)
library(tidyr)
library(tidyverse)
library(dplyr)
library(stats)
library(sf)
library(sp)
library(rgeos)
library(car)
#library(ggplot2)
library(usdm)
library(enmSdmX) #remotes::install_github('adamlilith/enmSdm', dependencies=TRUE)
library(sdm) #devtools::install_github("babaknaimi/sdm")
#library(regclass)
library(corpcor)
library(rJava)
library(terra)
library(exactextractr)
library(gam)
#library(safeBinaryRegression) #used as another check for collinearity
#remove.packages("safeBinaryRegression")

#set working directory
setwd("D:/Frog_datasets/")

#read in shapefile of FYLF occurrence points; points were converted from occupied stream segments in ArcGIS
ca_shp<-shapefile('D:/Frog_datasets/Occur_Sierra_270.shp')

#read in static variables (data downloaded from LANDFIRE https://landfire.gov/topographic.php)
aspect <- raster("Raster_Crop/Asp_clip.tif")
elevation <- raster("Raster_Crop/Elev_clip.tif")
slope <- raster("Raster_Crop/SLP_clip.tif")
static_v <- stack(aspect, elevation, slope)

#extract data to FYLF points
occsAES <- raster::extract(static_v, ca_shp)
occsAES <- as.data.frame(occsAES)
occsAES$ID <- (1:nrow(occsAES))

#extract veg data to FYLF points; this is coded out because each iteration of the code showed veg was not a significant predictor
#and FWS listing proposal echoed that FYLF found in streams with diverse surrounding vegetation
#veg <- raster("veg_reclass.tif")
#veg <- veg + 1
#veg <- as.factor(veg)
#veg_d <- raster::extract(veg, ca_shp)
#veg_d <- as.data.frame(veg_d)
#veg_d$ID <- (1:nrow(veg_d))

#extract stream order to FYLF points (downloaded from NHD+ https://www.usgs.gov/national-hydrography/nhdplus-high-resolution)
streamorder <- raster("Hydrography/streamo.tif")
streampts <- raster::extract(streamorder, ca_shp)
streampts <- as.data.frame(streampts)
streampts$ID <- (1:nrow(streampts))
sum(is.na(streampts$streampts))

#extract soil K factor (a proxy for soil erosion https://www.arcgis.com/home/item.html?id=ac1bc7c30bd4455e85f01fc51055e586) to FYLF points
#opted not to use this layer bc missing a ton of data, especially in national parks
#soil <- raster("SoilKfactor.tif")
#soilpts <- raster::extract(soil, ca_shp)
#soilpts <- as.data.frame(soilpts)
#soilpts$ID <- (1:nrow(soilpts))
#sum(is.na(soilpts$soilpts))

#extract stream temp (downloaded data from https://www.fs.usda.gov/rm/boise/AWAE/projects/NorWeST.html) to FYLF points
stemp <- raster("Raster_Crop_mess/Stemp_hist.tif")
stemp_hist <- raster::extract(stemp, ca_shp)
stemp_hist <- as.data.frame(stemp_hist)
stemp_hist$ID <- (1:nrow(stemp_hist))
sum(is.na(stemp_hist$stemp_hist))

#extract mean summer flows to FYLF points (downloaded from https://www.fs.usda.gov/rm/boise/AWAE/projects/VIC_streamflowmetrics/archived_modeled_stream_flow_metrics.shtml)
sflow <- raster("MS_flowhist1.tif")
sflow_hist <- raster::extract(sflow, ca_shp)
sflow_hist <- as.data.frame(sflow_hist)
sflow_hist$ID <- (1:nrow(sflow_hist))
sum(is.na(sflow_hist$sflow_hist))

#extract winter flows above 95th percentile to FYLF points (downloaded from https://www.fs.usda.gov/rm/boise/AWAE/projects/VIC_streamflowmetrics/archived_modeled_stream_flow_metrics.shtml)
wflow<- raster("W95_flowhist1.tif")
wflow_hist <- raster::extract(wflow, ca_shp)
wflow_hist <- as.data.frame(wflow_hist)
wflow_hist$ID <- (1:nrow(wflow_hist))
sum(is.na(wflow_hist$wflow_hist))

#extract stream alteration (downloaded from McManamay, R.A., A dataset of modeled hydrologic alteration and ecological consequences in stream
#reaches of the conterminous United States. Zenodo 10.5281/zenodo.5839011. (2022))
streamalt<- raster("stream_alt/Stream_alt.tif")
streamalt_hist <- raster::extract(streamalt, ca_shp)
streamalt_hist <- as.data.frame(streamalt_hist)
streamalt_hist$ID <- (1:nrow(streamalt_hist))
sum(is.na(streamalt_hist$streamalt_hist))

#merge all extracted environmental predictor variables
static_noveg<-merge(occsAES, stemp_hist, by="ID")
static_flow<-merge(static_noveg, sflow_hist, by="ID")
static_order<-merge(static_flow, streampts, by="ID")
#static_soil<-merge(static_order, soilpts, by="ID")
static_most<-merge(static_order, wflow_hist, by="ID")
static_all <- merge(static_most, streamalt_hist, by="ID")
#static_veg <- merge(static_noveg, veg_d, by="ID")
#static_veg_f<-static_veg[complete.cases(static_veg), ]

#human modification data received from FWS, D. Theobald data product
hm <- raster("HM_US_v3_2019_Pacific_90m.tif")

#read in FYLF occupied stream segments
#convert polylines to points
ptsFYLF = as(ca_shp, "SpatialPointsDataFrame")
ptsFYLF$ID <- (1:nrow(ptsFYLF))
#remove streamsegments outside of California
ptsFYLF <- ptsFYLF[ptsFYLF$state == 'CA',]
#get coordinates for geoThinning
pc<-coordinates(spTransform(ptsFYLF, CRS("+proj=longlat +datum=WGS84")))
ptsFYLF$long <- pc[,1]
ptsFYLF$lat <- pc[,2]
#convert to dataframe for geoThinning
ptsdf<-as.data.frame(ptsFYLF)

#if want to remove 2020 occupancy data to be used later for validation (I did not do this; validation from eDNA or egg mass):
#ptsdf <- ptsdf[ptsdf$last_year_ != '2020',]

#remove coastal clades... will eventually run 2 analyses, one for all clades and one for just mountain clades
ptsdf <- ptsdf[ptsdf$Clade != 'Northwest/North Coast',]
ptsdf <- ptsdf[ptsdf$Clade != 'West/Central Coast',]
ptsdf <- ptsdf[ptsdf$Clade != 'Southwest/South Coast',]


#get data from BCM model historical (2000-2010)
table <- NULL

timestep <- as.character(seq(from=2000, to=2010, by=1))

for (t in timestep) {
  #t=2000
print(t)

#set observation year
ptsyr <- ptsdf[ptsdf$last_year_ == t,]

#create for loop for each year from 2000 to 2010 to geoThin occupancy points and assign environmental variables to each point
thinned <- geoThin(ptsyr, longLat=c(16,17), 270)

#read in annual predictors (downloaded from CA Basin Characterization Model https://ca.water.usgs.gov/projects/reg_hydro/basin-characterization-model.html)
soil_water <- raster(paste0("BCM_hist/str",t,".tif"))
excess_water <- raster(paste0("BCM_hist/exc",t,".tif"))
cwd <- raster(paste0("BCM_hist/cwd",t,".tif"))
tmax <- raster(paste0("BCM_hist/tmx",t,".tif"))
tmin <- raster(paste0("BCM_hist/tmn",t,".tif"))
ppt <- raster(paste0("BCM_hist/ppt",t,".tif"))
snow <- raster(paste0("BCM_hist/snw",t,".tif"))

#prj <-projection(tmax)

#create stack of biophysical variables with which to overlay occurrence data
env = stack(soil_water, excess_water, cwd, tmax, tmin, ppt, snow)

FYLF_t <- ptsFYLF[ptsFYLF@data$ID %in% thinned$ID,]
FYLF_t <- spTransform(FYLF_t, CRS("+proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs "))
FYLF_t@bbox <- as.matrix(extent(soil_water))

occsEnv <- raster::extract(env, FYLF_t)
occsEnv <- as.data.frame(occsEnv) # need to do this for prediction later
occsEnv$ID <- FYLF_t$ID

#wanted to test whether buffered or unbuffered human modification had greater predictive power
#extract unbuffered human mod to each FYLF point
occsHM <- raster::extract(hm, FYLF_t)
occsHM <- as.data.frame(occsHM)
occsHM$ID <- FYLF_t$ID

#extract buffered human mod to each FYLF point
occsHM_buf <- raster::extract(hm,             # raster layer
                            FYLF_t,   # SPDF with centroids for buffer
                            buffer = 500,     # buffer size, units depend on CRS
                            fun=mean)         # what to value to extract
occsHM_buf <- as.data.frame(occsHM_buf)
occsHM_buf$ID <- FYLF_t$ID

colnames(occsEnv) <- c("soil_water", "excess_water", "cwd", "tmax", "tmin", "ppt", "snow", "ID")

occs_Env_Hm <- merge(occsEnv, occsHM, by="ID")
occs_all <- merge(occs_Env_Hm, occsHM_buf, by="ID")

table <- rbind(table, occs_all)
}


#get data from BCM model future (2011-2099)
table2 <- NULL

timestep2 <- as.character(seq(from=2011, to=2020, by=1))

for (i in timestep2) {
  #t=2000
  print(i)
  
  #set observation year
  ptsyr <- ptsdf[ptsdf$last_year_ == i,]
  
  #create for loop for each year from 2000 to 2010 to geoThin occupancy points and assign environmental variables to each point
  thinned <- geoThin(ptsyr, longLat=c(16,17), 270)
  
  #read in annual predictors (downloaded from CA Basin Characterization Model https://ca.water.usgs.gov/projects/reg_hydro/basin-characterization-model.html)
  soil_water <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85str",i,".tif"))
  excess_water <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85exc",i,".tif"))
  cwd <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85cwd",i,".tif"))
  tmax <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85tmx",i,".tif"))
  tmin <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85tmn",i,".tif"))
  ppt <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85ppt",i,".tif"))
  snow <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85snw",i,".tif"))
  
  #prj <-projection(tmax)
  
  #create stack of biophysical variables with which to overlay occurrence data
  env = stack(soil_water, excess_water, cwd, tmax, tmin, ppt, snow)

  FYLF_t <- ptsFYLF[ptsFYLF@data$ID %in% thinned$ID,]
  FYLF_t <- spTransform(FYLF_t, CRS("+proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs "))
  FYLF_t@bbox <- as.matrix(extent(soil_water))
  
  occsEnv <- raster::extract(env, FYLF_t)
  occsEnv <- as.data.frame(occsEnv) # need to do this for prediction later
  occsEnv$ID <- FYLF_t$ID
  
  occsHM <- raster::extract(hm, FYLF_t)
  occsHM <- as.data.frame(occsHM)
  occsHM$ID <- FYLF_t$ID
  
  occsHM_buf <- raster::extract(hm,             # raster layer
                                FYLF_t,   # SPDF with centroids for buffer
                                buffer = 500,     # buffer size, units depend on CRS
                                fun=mean)         # what to value to extract
  occsHM_buf <- as.data.frame(occsHM_buf)
  occsHM_buf$ID <- FYLF_t$ID
  
  colnames(occsEnv) <- c("soil_water", "excess_water", "cwd", "tmax", "tmin", "ppt", "snow", "ID")
  
  occs_Env_Hm <- merge(occsEnv, occsHM, by="ID")
  occs_all <- merge(occs_Env_Hm, occsHM_buf, by="ID")
  
  table2 <- rbind(table2, occs_all)
  
}

var <- rbind(table, table2)

#create table of environmental predictor variables for each FYLF point
env_var <- merge(var, static_all, by="ID")
env_var_f <- env_var[,c(1:8, 10:18)] #remove hm no buf
env_var_f<-env_var_f[complete.cases(env_var_f), ]
write.csv(env_var_f,"occurrences_10092023.csv", row.names = FALSE, col.names = TRUE)


#check for colinearity among variables using usdm package
v1 <- vifstep(env_var_f)
v2 <- vifcor(env_var_f)

v1
v2

#add column to indicate whether presence/absence point
env_var_f$species <- 1

pres_coords <- ptsdf[,c(15:17)]
env_var_f <- merge(env_var_f, pres_coords, by="ID")


#made random points in r for background points in species distribution modeling since no true absences
absences <- readOGR("Back270_rand.shp")
#absences <- spsample(sierra,n=5000,"random")


#absences<-coordinates(spTransform(absences, CRS("+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")))
absences_df <- as.data.frame(absences)
absences_df$ID <- (1:nrow(absences_df))

occsAES_ab <- raster::extract(static_v, absences)
occsAES_ab <- as.data.frame(occsAES_ab)
occsAES_ab$ID <- (1:nrow(occsAES_ab))

#again, didn't wind up using veg data
#veg_d_ab <- raster::extract(veg, absences)
#veg_d_ab <- as.data.frame(veg_d_ab)
#veg_d_ab$ID <- (1:nrow(veg_d_ab))

stemp_hist_ab <- raster::extract(stemp, absences)
stemp_hist_ab <- as.data.frame(stemp_hist_ab)
stemp_hist_ab$ID <- (1:nrow(stemp_hist_ab))
colnames(stemp_hist_ab) <- c("stemp_hist", "ID")

sflow_hist_ab <- raster::extract(sflow, absences)
sflow_hist_ab <- as.data.frame(sflow_hist_ab)
sflow_hist_ab$ID <- (1:nrow(sflow_hist_ab))

streampts_ab <- raster::extract(streamorder, absences)
streampts_ab <- as.data.frame(streampts_ab)
streampts_ab$ID <- (1:nrow(streampts_ab))
#sum(is.na(streampts$streampts))

wflow_hist_ab <- raster::extract(wflow, absences)
wflow_hist_ab <- as.data.frame(wflow_hist_ab)
wflow_hist_ab$ID <- (1:nrow(wflow_hist_ab))

streamalt_ab <- raster::extract(streamalt, absences)
streamalt_ab <- as.data.frame(streamalt_ab)
streamalt_ab$ID <- (1:nrow(streamalt_ab))
sum(is.na(streamalt_ab$streamalt_ab))

occsHM_ab <- raster::extract(hm, absences)
occsHM_ab <- as.data.frame(occsHM_ab)
occsHM_ab$ID <- (1:nrow(occsHM_ab))
colnames(occsHM_ab) <- c("occsHM", "ID")

occsHM_buf_ab <- raster::extract(hm,             # raster layer
                              absences,   # SPDF with centroids for buffer
                              buffer = 500,     # buffer size, units depend on CRS
                              fun=mean)         # what to value to extract
occsHM_buf_ab <- as.data.frame(occsHM_buf_ab)
occsHM_buf_ab$ID <- (1:nrow(occsHM_buf_ab))
colnames(occsHM_buf_ab) <- c("occsHM_buf", "ID")


static_hm_ab<- merge(occsHM_ab, occsHM_buf_ab, by="ID")
static_layers_ab <- merge(static_hm_ab, occsAES_ab, by="ID")
static_stemp_ab<- merge(static_layers_ab, stemp_hist_ab, by="ID")
static_flow_ab <- merge(static_stemp_ab, sflow_hist_ab, by="ID")
static_flow_pts <- merge(static_flow_ab, streampts_ab)
static_flow_w <- merge(static_flow_pts, wflow_hist_ab)
static_all_ab <- merge(static_flow_w, streamalt_ab, by="ID")
#static_all_ab_f<-static_all_ab[complete.cases(static_all_ab), ]

#need to assign random year to background points for data with temporal component (e.g., Basin Characterization Model data)
year <- rep(c(2000:2020), each = 184)
year <- as.data.frame(year)
colnames(year) <- "year"
year2 <- rep(c(2000:2006))
year2 <- as.data.frame(year2)
colnames(year2) <- "year"
years <- rbind(year, year2)

absences_df$year <- years
absences$ID <- (1:nrow(absences_df))

table_ab <- NULL

timestep <- as.character(seq(from=2000, to=2010, by=1))

for (t in timestep) {
  #t=2000
  print(t)
  
  #set observation year
  ptsyr_ab <- absences_df[absences_df$year == t,]
  
  #create_ab for loop for each year from 2000 to 2010 to geoThin occupancy points and assign environmental variables to each point
  thinned <- geoThin(ptsyr_ab, longLat=c(2,3), 270)
  
  #read in annual predictors (downloaded from CA Basin Characterization Model https://ca.water.usgs.gov/projects/reg_hydro/basin-characterization-model.html)
  soil_water <- raster(paste0("BCM_hist/str",t,".tif"))
  excess_water <- raster(paste0("BCM_hist/exc",t,".tif"))
  cwd <- raster(paste0("BCM_hist/cwd",t,".tif"))
  tmax <- raster(paste0("BCM_hist/tmx",t,".tif"))
  tmin <- raster(paste0("BCM_hist/tmn",t,".tif"))
  ppt <- raster(paste0("BCM_hist/ppt",t,".tif"))
  snow <- raster(paste0("BCM_hist/snw",t,".tif"))
  
  #prj <-projection(tmax)
  
  #create stack of biophysical variables with which to overlay occurrence data
  env = stack(soil_water, excess_water, cwd, tmax, tmin, ppt, snow)
  
  FYLF_ab <- absences[absences@data$ID %in% thinned$ID,]
  FYLF_ab <- spTransform(FYLF_ab, CRS("+proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs "))
  FYLF_ab@bbox <- as.matrix(extent(soil_water))
  
  occsEnv_ab <- raster::extract(env, FYLF_ab)
  occsEnv_ab <- as.data.frame(occsEnv_ab) # need to do this for prediction later
  occsEnv_ab$ID <- FYLF_ab$ID
  
  colnames(occsEnv_ab) <- c("soil_water", "excess_water", "cwd", "tmax", "tmin", "ppt", "snow", "ID")
  
  table_ab <- rbind(table_ab, occsEnv_ab)
}


#get data from BCM model future (2011-2099)

table2_ab <- NULL

timestep2 <- as.character(seq(from=2011, to=2020, by=1))

for (i in timestep2) {
  #t=2000
  print(i)
  
  #set observation year
  ptsyr <- ptsdf[ptsdf$last_year_ == i,]
  
    #create_ab for loop for each year from 2000 to 2010 to geoThin occupancy points and assign environmental variables to each point
    thinned <- geoThin(ptsyr_ab, longLat=c(2,3), 270)
    
    #read in annual predictors (downloaded from CA Basin Characterization Model https://ca.water.usgs.gov/projects/reg_hydro/basin-characterization-model.html)
    soil_water <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85str",i,".tif"))
    excess_water <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85exc",i,".tif"))
    cwd <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85cwd",i,".tif"))
    tmax <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85tmx",i,".tif"))
    tmin <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85tmn",i,".tif"))
    ppt <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85ppt",i,".tif"))
    snow <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85snw",i,".tif"))
    
    #prj <-projection(tmax)
    
    #create stack of biophysical variables with which to overlay occurrence data
    env = stack(soil_water, excess_water, cwd, tmax, tmin, ppt, snow)
    
    FYLF_ab <- absences[absences@data$ID %in% thinned$ID,]
    FYLF_ab <- spTransform(FYLF_ab, CRS("+proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs "))
    FYLF_ab@bbox <- as.matrix(extent(soil_water))
    
    occsEnv_ab <- raster::extract(env, FYLF_ab)
    occsEnv_ab <- as.data.frame(occsEnv_ab) # need to do this for prediction later
    occsEnv_ab$ID <- FYLF_ab$ID
    
    colnames(occsEnv_ab) <- c("soil_water", "excess_water", "cwd", "tmax", "tmin", "ppt", "snow", "ID")
    
    table2_ab <- rbind(table2_ab, occsEnv_ab)
  }

var_ab <- rbind(table_ab, table2_ab)


env_var_ab <- merge(var_ab, static_all_ab, by="ID")
env_var_ab_f <- env_var_ab[,c(1:8, 10:18)] #remove hm no buf
env_var_ab_f<-env_var_ab_f[complete.cases(env_var_ab_f), ]
write.csv(env_var_ab_f,"background_10092023.csv", row.names = FALSE, col.names = TRUE)


#add column to indicate whether presence/absence point
env_var_ab_f$species <- 0

#get species points for sdm modeling
abs_coords <- absences_df[,c(4, 2:3)]
colnames(abs_coords) <- c("ID", "long", "lat")
env_var_ab_f <- merge (env_var_ab_f, abs_coords, by="ID")


colnames(env_var_f) <- c("ID","soil_water","excess_water","cwd","tmax","tmin","ppt","snow","HM_buf","asp","elev",'slp',"stemp_hist","sflow_hist","streampts","wflow_hist","alt","species","long","lat")
colnames(env_var_ab_f) <- c("ID","soil_water","excess_water","cwd","tmax","tmin","ppt","snow","HM_buf","asp","elev",'slp',"stemp_hist","sflow_hist","streampts","wflow_hist","alt","species","long","lat")


all_data <- rbind (env_var_f, env_var_ab_f)
#colnames(all_data) <- c("ID","soil_water","excess_water","cwd","tmax","tmin","ppt","snow","HM","HM_buf","asp","elev",'slp',"stemp_hist","sflow_hist","wflow_hist","species","long","lat")
write.csv(all_data, "all_data_10092023.csv", row.names = FALSE, col.names = TRUE)

all_data <- read.csv("all_data_10092023.csv")
all_data <- all_data[complete.cases(all_data), ]
#all_data$veg <- as.factor(all_data$veg)

pres <- all_data[all_data$species ==1,]

abs <- all_data[all_data$species ==0,]

#linear regression to determine predictive power of each environmental variable
#also test whether greater predictive power buffering or not buffering around points for human modification with linear regression

#run global model to find AIC of all variables
#unbuffered human mod
glmALL_nobuf <-stats::glm(species ~ soil_water+excess_water+cwd+tmax+tmin+ppt+snow+HM+asp+slp+stemp_hist+sflow_hist+wflow_hist, data= all_data, family ="binomial")
summary(glmALL_nobuf)
#do stepwise regression to find combination of predictor variables with lowest AIC + help reduce number of variables used to reduce model overfitting
slmALL_nobuf <- step(glmALL_nobuf)
summary(slmALL_nobuf)

#buffered human mod
glmALL_buf <-stats::glm(species ~ soil_water+excess_water+cwd+tmax+tmin+ppt+snow+HM_buf+asp+slp+stemp_hist+sflow_hist+streampts+wflow_hist+alt, data= all_data, family ="binomial")
summary(glmALL_buf)
slmALL_buf <- step(glmALL_buf)
summary(slmALL_buf)

glmALL_nocolin <-stats::glm(species ~ soil_water+cwd+tmax+ppt+snow+HM_buf+asp+slp+stemp_hist+sflow_hist+streampts+wflow_hist+alt, data= all_data, family ="binomial")
summary(glmALL_nocolin)
slmALL_nocolin <- step(glmALL_nocolin)
summary(slmALL_nocolin)


#check again for multicollinearity
car::vif(glmALL_buf)

pres_data_var <-all_new[all_new$species==1,]
pres_data_var <- pres_data_var[,c(3:13, 18:21)]
cor2pcor(cov(pres_data_var))
slmALL_buf <- step(glmALL_buf)
summary(slmALL_buf)

#got warning when running glm that the model did not coverge
#ran with package safeBinaryRegression to figure out what was causing the warning:
#1: glm.fit: algorithm did not converge 
#2: glm.fit: fitted probabilities numerically 0 or 1 occurred
#glm_test<-glm(species ~ soil_water+excess_water+cwd+tmax+tmin+ppt+snow+asp+elev+slp+stemp_hist+veg, data= all_data, family ="binomial",
    #separation = c("find", "test"))
#After running glm through the package safeBinaryResgression, the following terms are causing separation among the sample points: excess_water, snow

#also scaled data then did linear regression to see if variable range was inflating variable importance
scaled_data <- all_data

#normalize/standardize variables
#have not updated this part of code!!
scaled_data$soil_water<-scale(scaled_data$soil_water)
scaled_data$excess_water<-scale(scaled_data$excess_water)
scaled_data$cwd<-scale(scaled_data$cwd)
scaled_data$tmax<-scale(scaled_data$tmax)
scaled_data$tmin<-scale(scaled_data$tmin)
scaled_data$ppt<-scale(scaled_data$ppt)
scaled_data$snow<-scale(scaled_data$snow)
scaled_data$HM<-scale(scaled_data$HM)
scaled_data$HM_buf<-scale(scaled_data$HM_buf)
scaled_data$asp<-scale(scaled_data$asp)
scaled_data$elev<-scale(scaled_data$elev)
scaled_data$slp<-scale(scaled_data$slp)
scaled_data$stemp_hist<-scale(scaled_data$stemp_hist)
scaled_data$sflow_hist<-scale(scaled_data$sflow_hist)
scaled_data$wflow_hist<-scale(scaled_data$wflow_hist)

str(all_data)

write.csv(scaled_data, "scaled_data_flow2.csv", row.names = FALSE, col.names = TRUE)

scaled_data <- read.csv("scaled_data_flow2.csv")


#test how a few variable performs independently to get reference AIC
V1 <- glm(species ~ cwd_scaled, data = scaled_data, family = "binomial")
V2 <- glm(species ~ ppt_scaled, data = scaled_data, family = "binomial")
#V3 <- glm(species ~ veg, data = scaled_data, family = "binomial")
V4 <- glm(species ~ tmax_scaled, data = scaled_data, family = "binomial")
V5 <- glm(species ~ tmin_scaled, data = scaled_data, family = "binomial")
V6 <- glm(species ~ elev_scaled, data = scaled_data, family = "binomial")
AIC(V1,V2,V3,V4,V5,V6)

#run global model to find AIC of scaled variables
glmscaled_nobuf <-stats::glm(species ~ soil_water+excess_water+cwd+tmax+tmin+ppt+snow+HM+asp+elev+slp+stemp_hist+sflow_hist+wflow_hist, data= scaled_data, family ="binomial")
summary(glmscaled_nobuf)
slmscaled_nobuf <- step(glmscaled_nobuf)
summary(slmscaled_nobuf)

glmscaled_buf <-stats::glm(species ~ soil_water+excess_water+cwd+tmax+tmin+ppt+snow+HM_buf+asp+elev+slp+stemp_hist+sflow_hist+wflow_hist, data= scaled_data, family ="binomial")
summary(glmscaled_buf)
slmscaled_buf <- step(glmscaled_buf)
summary(slmscaled_buf)


#SDM modeling
#use significant predictors identified in stepwise glm
installAll(sdm)

d <- sdmData(species~soil_water+cwd+tmax+ppt+snow+HM_buf+asp+slp+sflow_hist+streampts+wflow_hist+coords(long+lat),train=all_data)

getmethodNames()

#if need to include veg and get error about factor (vector name length mismatch error) try names(m) <- envvar
m <- sdm::sdm(species~., d, methods = c("glm","rf","brt","maxent"), replication = c("boot"), n=500)
write.sdm(m,'sdm')

#m <- read.sdm("sdm.sdm")

#remove some data to clear memory

rm(streamalt)
#rm(glmALL_ns)
rm(ca_shp)
rm(elevation)
rm(excess_water)
rm(tmin)
rm(absences)
rm(env)
rm(FYLF_ab)
rm(FYLF_t)
rm(ptsFYLF)
rm(static_v)
rm(stemp)
rm(test_r)
#rm(slmALL_ns)


######################2020
#predict potential current habitat suitability

#read in annual predictors
soil_water <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85str2020.tif"))
#excess_water <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85exc2020.tif"))
cwd <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85cwd2020.tif"))
tmax <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85tmx2020.tif"))
#tmin <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85tmn2020.tif"))
ppt <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85ppt2020.tif"))
snow <- raster(paste0("BCM_future/CNRM_rcp85/CNRM_rcp85snw2020.tif"))

#reproject static layers
slope_r <- projectRaster(slope, cwd)
slope_r <- setExtent(slope_r, cwd)
slp <- slope_r

aspect_r <- projectRaster(aspect, cwd)
aspect_r <- setExtent(aspect_r, cwd)
asp <- aspect_r

hm_r <- projectRaster(hm, cwd)
hm_r <- setExtent(hm_r, cwd)

sflow_hist_r <- projectRaster(sflow, cwd)
sflow_hist_r <- setExtent(sflow_hist_r, cwd)

wflow_hist_r <- projectRaster(wflow, cwd)
wflow_hist_r <- setExtent(wflow_hist_r, cwd)

streamorder_r <- projectRaster(streamorder, cwd)
streamorder_r <- setExtent(wflow_hist_r, cwd)

#soil_r <- projectRaster(soil, cwd)
#soil_r <- setExtent(soil_r, cwd)

bio_curr<-raster::stack(soil_water,cwd,tmax,ppt,snow,hm_r,asp,slp,sflow_hist_r,streamorder_r,wflow_hist_r)
names(bio_curr) <- c("soil_water","cwd","tmax","ppt","snow","HM_buf","asp","slp","sflow_hist","streampts","wflow_hist")
#names(bio_curr) <- paste0(c("soil_water","cwd","tmax","ppt","snow","HM_buf","asp","slp","sflow_hist","streampts","wflow_hist"))


bio_curr_r<-projectRaster(bio_curr, crs = ("+proj=longlat +datum=WGS84"))

ens_2020 <- ensemble(m, newdata=bio_curr_r, filename='hs_2020_10312023.tif',setting=list(method='weighted',stat='AUC'))

######################2080
#predict future habitat suitability in 2080

#read in annual predictors
predictor <- c("tmx","snw","ppt","str","cwd")
#mods <- c("MIROC_rcp85", "CCSM4_rcp85", "CNRM_rcp85", "FGOALS_rcp85")
mods <- c("CCSM4_rcp85", "CNRM_rcp85", "FGOALS_rcp85")

sflow <- raster("MS_flow20801.tif")
sflow_2080_r <- projectRaster(sflow, cwd)
sflow_2080_r <- setExtent(sflow_2080_r, cwd)

wflow <- raster("W95_flow20801.tif")
wflow_2080_r <- projectRaster(sflow, cwd)
wflow_2080_r <- setExtent(wflow_2080_r, cwd)

for (o in mods) {
  #o="MIROC_rcp45"
  print(o)

for (p in predictor) {
  #p="tmn"
  #read in annual predictors
  p2080 <- raster(paste0("BCM_future/",o,"/",o,p,"2080.tif"))
  p2081 <- raster(paste0("BCM_future/",o,"/",o,p,"2081.tif"))
  p2082 <- raster(paste0("BCM_future/",o,"/",o,p,"2082.tif"))
  p2083 <- raster(paste0("BCM_future/",o,"/",o,p,"2083.tif"))
  p2084 <- raster(paste0("BCM_future/",o,"/",o,p,"2084.tif"))
  p2085 <- raster(paste0("BCM_future/",o,"/",o,p,"2085.tif"))
  p2086 <- raster(paste0("BCM_future/",o,"/",o,p,"2086.tif"))
  p2087 <- raster(paste0("BCM_future/",o,"/",o,p,"2087.tif"))
  p2088 <- raster(paste0("BCM_future/",o,"/",o,p,"2088.tif"))
  p2089 <- raster(paste0("BCM_future/",o,"/",o,p,"2089.tif"))
  
  stacky <-stack(p2080,p2081,p2082,p2083,p2084,p2085,p2086,p2087,p2088,p2089)
  
  mean_s <- mean(stacky)
  
  assign(paste0(p,"2080"), mean_s)
  
  #writeRaster(mean_s, filename = paste0("BCM_future/MIROC_rcp85/",p,"2080.tif"))
  
}

#stemp_2080_r <- raster(paste0("Stream_temp_future/stemp_2080v2.tif"))
#stemp_2080_r <- projectRaster(stemp_hist_r, cwd)


bio_curr<-raster::stack(str2080,cwd2080,tmx2080,ppt2080,snw2080,hm_r,asp,slp,sflow_2080_r,streamorder_r,wflow_2080_r)
names(bio_curr) <- c("soil_water","cwd","tmax","ppt","snow","HM_buf","asp","slp","sflow_hist","streampts","wflow_hist")

bio_curr_2080<-projectRaster(bio_curr, crs = ("+proj=longlat +datum=WGS84"))

ens_2080 <- ensemble(m, newdata=bio_curr_2080, filename=paste0('hs_2080_',o,'_10312023.tif'),setting=list(method='weighted',stat='AUC'))
}

######################2040
#predict future habitat suitability in 2040

#read in annual predictors
predictor <- c("tmx","snw","ppt","str","cwd")
#mods <- c("MIROC_rcp85", "CCSM4_rcp85", "CNRM_rcp85", "FGOALS_rcp85")
mods <- c("MIROC_rcp85","CCSM4_rcp85", "CNRM_rcp85", "FGOALS_rcp85")

sflow <- raster("MS_flow20401.tif")
sflow_2040_r <- projectRaster(sflow, cwd)
sflow_2040_r <- setExtent(sflow_2040_r, cwd)

wflow <- raster("W95_flow20401.tif")
wflow_2040_r <- projectRaster(sflow, cwd)
wflow_2040_r <- setExtent(wflow_2040_r, cwd)

for (o in mods) {
  #o="MIROC_rcp45"
  print(o)
  
  for (p in predictor) {
    #p="tmn"
    #read in annual predictors
    p2040 <- raster(paste0("BCM_future/",o,"/",o,p,"2040.tif"))
    p2081 <- raster(paste0("BCM_future/",o,"/",o,p,"2081.tif"))
    p2082 <- raster(paste0("BCM_future/",o,"/",o,p,"2082.tif"))
    p2083 <- raster(paste0("BCM_future/",o,"/",o,p,"2083.tif"))
    p2084 <- raster(paste0("BCM_future/",o,"/",o,p,"2084.tif"))
    p2085 <- raster(paste0("BCM_future/",o,"/",o,p,"2085.tif"))
    p2086 <- raster(paste0("BCM_future/",o,"/",o,p,"2086.tif"))
    p2087 <- raster(paste0("BCM_future/",o,"/",o,p,"2087.tif"))
    p2088 <- raster(paste0("BCM_future/",o,"/",o,p,"2088.tif"))
    p2089 <- raster(paste0("BCM_future/",o,"/",o,p,"2089.tif"))
    
    stacky <-stack(p2040,p2081,p2082,p2083,p2084,p2085,p2086,p2087,p2088,p2089)
    
    mean_s <- mean(stacky)
    
    assign(paste0(p,"2040"), mean_s)
    
    #writeRaster(mean_s, filename = paste0("BCM_future/MIROC_rcp85/",p,"2040.tif"))
    
  }
  
  #stemp_2040_r <- raster(paste0("Stream_temp_future/stemp_2040v2.tif"))
  #stemp_2040_r <- projectRaster(stemp_hist_r, cwd)
  
  
  bio_curr<-raster::stack(str2040,cwd2040,tmx2040,ppt2040,snw2040,hm_r,asp,slp,sflow_2040_r,streamorder_r,wflow_2040_r)
  names(bio_curr) <- c("soil_water","cwd","tmax","ppt","snow","HM_buf","asp","slp","sflow_hist","streampts","wflow_hist")
  
  bio_curr_2040<-projectRaster(bio_curr, crs = ("+proj=longlat +datum=WGS84"))
  
  ens_2040 <- ensemble(m, newdata=bio_curr_2040, filename=paste0('hs_2040_',o,'_10312023.tif'),setting=list(method='weighted',stat='AUC'))
}

#find CC refugia + transition to newly suitable habitat

###FYLF refugia

#read in habitat suitability rasters for 2010 and 2080 (can also do same steps for 2040 time step)
suit_2020<- raster('hs_2020_10312023.tif')

mods <- c("CCSM4_rcp85", "CNRM_rcp85", "FGOALS_rcp85","MIROC_rcp85")

for (o in mods) {
  #o="CCSM4_rcp85"
  print(o)

suit_2020<- raster('hs_2020_10312023.tif')
suit_2040 <- raster(paste0('hs_2040_',o,'.tif'))
suit_2080 <- raster(paste0('hs_2080_',o,'.tif'))

#suitability threshold set to 0.5
#set conditions to give any cell less than 0.5 an NA value
suit_2020[suit_2020 < 0.5] <- NA
suit_2040[suit_2040 < 0.5] <- NA
suit_2080[suit_2080 < 0.5] <- NA

#subtract historic habitat suitability from projected future FYLF for conditional change in habitat suitability map
FYLF_r40 <- (suit_2040 - suit_2020)
FYLF_r80 <- (suit_2080 - suit_2020)

#set any cell that was NA under our conditions that changed in the raster math back to NA (R is weird about this,
#ArcGIS did it automatically.)
FYLF_r40[is.na(suit_2020)] <- NA
FYLF_r40[is.na(suit_2040)] <- NA
FYLF_r80[is.na(suit_2020)] <- NA
FYLF_r80[is.na(suit_2080)] <- NA

#plot(FYLF_r)
writeRaster(FYLF_r40, paste0('ref_2040_',o,'.tif'), overwrite = TRUE)
writeRaster(FYLF_r80, paste0('ref_2080_',o,'.tif'), overwrite = TRUE)
}

###FYLF transition

#read in habitat suitability rasters for 2010 and 2080
suit_2020<- raster('hs_2020_10092023.tif')

for (o in mods) {
  #o="MIROC_rcp85"
  print(o)
  
suit_2040 <- raster(paste0('hs_2040_',o,'.tif'))
suit_2080 <- raster(paste0('hs_2080_',o,'.tif'))

#set conditions to give any cell less than 0.5 an NA value
suit_2020[suit_2020 >= 0.5] <- NA
suit_2040[suit_2040 < 0.5] <- NA
suit_2080[suit_2080 < 0.5] <- NA

#subtract historic habitat suitability from projected future FYLF for conditional change in habitat suitability map
FYLF_t40 <- (suit_2040 - suit_2020)
FYLF_t80 <- (suit_2080 - suit_2020)

#set any cell that was NA under our conditions that changed in the raster math back to NA (R is weird about this,
#ArcGIS did it automatically.)
FYLF_t40[is.na(suit_2020)] <- NA
FYLF_t40[is.na(suit_2040)] <- NA
FYLF_t80[is.na(suit_2020)] <- NA
FYLF_t80[is.na(suit_2080)] <- NA

#plot(FYLF_t)
writeRaster(FYLF_t40, paste0('tran_2040_',o,'.tif'), overwrite = TRUE)
writeRaster(FYLF_t80, paste0('tran_2080_',o,'.tif'), overwrite = TRUE)
}

#find means
#2080
MIROC_rcp80 <- raster("hs_2080_MIROC_rcp85.tif")
CCSM4_rcp80 <- raster("hs_2080_CCSM4_rcp85.tif")
CNRM_rcp80 <- raster("hs_2080_CNRM_rcp85.tif")
FGOALS_rcp80 <- raster("hs_2080_FGOALS_rcp85.tif")
stacky80 <- stack(MIROC_rcp80,CCSM4_rcp80,CNRM_rcp80,FGOALS_rcp80)
mean80 <- mean(stacky80)

writeRaster(mean80, paste0('hs_mean_2080.tif'), overwrite = TRUE)

#2040
MIROC_rcp40 <- raster("hs_2040_MIROC_rcp85.tif")
CCSM4_rcp40 <- raster("hs_2040_CCSM4_rcp85.tif")
CNRM_rcp40 <- raster("hs_2040_CNRM_rcp85.tif")
FGOALS_rcp40 <- raster("hs_2040_FGOALS_rcp85.tif")
stacky40 <- stack(MIROC_rcp40,CCSM4_rcp40,CNRM_rcp40,FGOALS_rcp40)
mean40 <- mean(stacky40)

writeRaster(mean40, paste0('hs_mean_2040.tif'), overwrite = TRUE)


###find refugia and transition using means rasters
suit_2020<- raster('hs_2020_10092023.tif')
suit_2040 <- mean40
suit_2080 <- mean80

#suitability threshold set to 0.5
#set conditions to give any cell less than 0.5 an NA value
suit_2020[suit_2020 < 0.5] <- NA
suit_2040[suit_2040 < 0.5] <- NA
suit_2080[suit_2080 < 0.5] <- NA

#subtract historic habitat suitability from projected future FYLF for conditional change in habitat suitability map
FYLF_r40 <- (suit_2040 - suit_2020)
FYLF_r80 <- (suit_2080 - suit_2020)

#set any cell that was NA under our conditions that changed in the raster math back to NA (R is weird about this,
#ArcGIS did it automatically.)
FYLF_r40[is.na(suit_2020)] <- NA
FYLF_r40[is.na(suit_2040)] <- NA
FYLF_r80[is.na(suit_2020)] <- NA
FYLF_r80[is.na(suit_2080)] <- NA

#plot(FYLF_r)
writeRaster(FYLF_r40, paste0('ref_2040_mean.tif'), overwrite = TRUE)
writeRaster(FYLF_r80, paste0('ref_2080_mean.tif'), overwrite = TRUE)


###FYLF transition

#read in habitat suitability rasters for 2010 and 2080
suit_2020<- raster('hs_2020_10092023.tif')
suit_2040 <- mean40
suit_2080 <- mean80

#set conditions to give any cell less than 0.5 an NA value
suit_2020[suit_2020 >= 0.5] <- NA
suit_2040[suit_2040 < 0.5] <- NA
suit_2080[suit_2080 < 0.5] <- NA

#subtract historic habitat suitability from projected future FYLF for conditional change in habitat suitability map
FYLF_t40 <- (suit_2040 - suit_2020)
FYLF_t80 <- (suit_2080 - suit_2020)

#set any cell that was NA under our conditions that changed in the raster math back to NA (R is weird about this,
#ArcGIS did it automatically.)
FYLF_t40[is.na(suit_2020)] <- NA
FYLF_t40[is.na(suit_2040)] <- NA
FYLF_t80[is.na(suit_2020)] <- NA
FYLF_t80[is.na(suit_2080)] <- NA

#plot(FYLF_t)
writeRaster(FYLF_t40, paste0('tran_2040_mean.tif'), overwrite = TRUE)
writeRaster(FYLF_t80, paste0('tran_2080_mean.tif'), overwrite = TRUE)

###############################################
#find km of suitable habitat
sierras <- readOGR('D:/Frog_datasets/sierra_clades.shp')
sierras_t <- spTransform(sierras, CRS("+proj=longlat +datum=WGS84"))
#sierras_t@bbox <- as.matrix(extent(suit_2020))

suit_2020<- raster('hs_2020_10092023.tif')
suit_2020[suit_2020 < 0.5] <- NA
suit_2020[suit_2020 >= 0.5] <- 1

masked <- mask(x = suit_2020, mask = sierras_t)
hs_2020s <- raster::crop(masked, extent(sierras_t))
hs_2020s[hs_2020s < 1] <- NA

hs_2020r<- rast(hs_2020s)
test <- cellSize(hs_2020r, unit ="km")
#test_sum <- sum(test, na.rm=TRUE)
hs2020_sum<- exact_extract(test, sierras_t, fun="sum")


plot(hs_2020s)
lines(sierras_t)

mods <- c("CCSM4_rcp85", "CNRM_rcp85", "FGOALS_rcp85","MIROC_rcp85")

sum_2040 <- as.data.frame(c(1:6))
sum_2080 <- as.data.frame(c(1:6))

for (o in mods) {
  #o="MIROC_rcp85"
  print(o)
  
  suit_2040 <- raster(paste0('hs_2040_',o,'.tif'))
  suit_2080 <- raster(paste0('hs_2080_',o,'.tif'))
  
  suit_2040[suit_2040 < 0.5] <- NA
  suit_2040[suit_2040 >= 0.5] <- 1
  
  masked <- mask(x = suit_2040, mask = sierras_t)
  hs_2040s <- raster::crop(masked, extent(sierras_t))
  hs_2040s[hs_2040s < 1] <- NA
  
  hs_2040r<- rast(hs_2040s)
  test <- cellSize(hs_2040r, unit ="km")
  hs2040_sum<- exact_extract(test, sierras_t, fun="sum")
  hs2040_sum<- as.data.frame(hs2040_sum)
  sum_2040 <- cbind(sum_2040, hs2040_sum)
  
  suit_2080[suit_2080 < 0.5] <- NA
  suit_2080[suit_2080 >= 0.5] <- 1
  
  masked <- mask(x = suit_2080, mask = sierras_t)
  hs_2080s <- raster::crop(masked, extent(sierras_t))
  hs_2080s[hs_2080s < 1] <- NA
  
  hs_2080r<- rast(hs_2080s)
  test <- cellSize(hs_2080r, unit ="km")
  hs2080_sum<- exact_extract(test, sierras_t, fun="sum")
  hs2080_sum<- as.data.frame(hs2080_sum)
  sum_2080 <- cbind(sum_2080, hs2080_sum)
  }

colSums(sum_2040, na.rm = FALSE)
colSums(sum_2080, na.rm = FALSE)

hs_mean_2040<- raster('hs_mean_2040.tif')
hs_mean_2080<- raster('hs_mean_2080.tif')

hs_mean_2040[hs_mean_2040 < 0.5] <- NA
hs_mean_2040[hs_mean_2040 >= 0.5] <- 1

masked <- mask(x = hs_mean_2040, mask = sierras_t)
hs_mean_2040s <- raster::crop(masked, extent(sierras_t))
hs_mean_2040s[hs_mean_2040s < 1] <- NA

hs_mean_2040r<- rast(hs_mean_2040s)
test <- cellSize(hs_mean_2040r, unit ="km")
hsmean_2040_sum<- exact_extract(test, sierras_t, fun="sum")
hsmean_2040_sum<- as.data.frame(hsmean_2040_sum)

hs_mean_2080[hs_mean_2080 < 0.5] <- NA
hs_mean_2080[hs_mean_2080 >= 0.5] <- 1

masked <- mask(x = hs_mean_2080, mask = sierras_t)
hs_mean_2080s <- raster::crop(masked, extent(sierras_t))
hs_mean_2080s[hs_mean_2080s < 1] <- NA

hs_mean_2080r<- rast(hs_mean_2080s)
test <- cellSize(hs_mean_2080r, unit ="km")
hsmean_2080_sum<- exact_extract(test, sierras_t, fun="sum")
hsmean_2080_sum<- as.data.frame(hsmean_2080_sum)

colSums(hsmean_2040_sum, na.rm = FALSE)
colSums(hsmean_2080_sum, na.rm = FALSE)


#find refugia river km 
sum_2040 <- as.data.frame(c(1:6))
sum_2080 <- as.data.frame(c(1:6))

mods <- c("CCSM4_rcp85", "CNRM_rcp85", "FGOALS_rcp85","MIROC_rcp85")

for (o in mods) {
  #o="CCSM4_rcp85"
  print(o)

  ref_2040 <- raster(paste0('ref_2040_',o,'.tif'))
  ref_2080 <- raster(paste0('ref_2080_',o,'.tif'))

  #ref_2040[ref_2040 >= 0] <- 1
  
  masked <- mask(x = ref_2040, mask = sierras_t)
  ref_2040s <- raster::crop(masked, extent(sierras_t))
  #ref_2040s[ref_2040s < 1] <- NA
  
  ref_2040r<- rast(ref_2040s)
  test <- cellSize(ref_2040r, unit ="km")
  ref2040_sum<- exact_extract(test, sierras_t, fun="sum")
  ref2040_sum<- as.data.frame(ref2040_sum)
  sum_2040 <- cbind(sum_2040, ref2040_sum)
  
  masked <- mask(x = ref_2080, mask = sierras_t)
  ref_2080s <- raster::crop(masked, extent(sierras_t))
  #ref_2080s[ref_2080s < 1] <- NA
  
  ref_2080r<- rast(ref_2080s)
  test <- cellSize(ref_2080r, unit ="km")
  ref2080_sum<- exact_extract(test, sierras_t, fun="sum")
  ref2080_sum<- as.data.frame(ref2080_sum)
  sum_2080 <- cbind(sum_2080, ref2080_sum)
}
colSums(sum_2040, na.rm = FALSE)
colSums(sum_2080, na.rm = FALSE)
#mean
ref_2040 <- raster(paste0('ref_2040_mean.tif'))
ref_2080 <- raster(paste0('ref_2080_mean.tif'))

masked <- mask(x = ref_2040, mask = sierras_t)
ref_2040s <- raster::crop(masked, extent(sierras_t))

ref_2040r<- rast(ref_2040s)
test <- cellSize(ref_2040r, unit ="km")
ref2040_sum<- exact_extract(test, sierras_t, fun="sum")
ref2040_sum<- as.data.frame(ref2040_sum)
sum_2040 <- cbind(sum_2040, ref2040_sum)

masked <- mask(x = ref_2080, mask = sierras_t)
ref_2080s <- raster::crop(masked, extent(sierras_t))

ref_2080r<- rast(ref_2080s)
test <- cellSize(ref_2080r, unit ="km")
ref2080_sum<- exact_extract(test, sierras_t, fun="sum")
ref2080_sum<- as.data.frame(ref2080_sum)

colSums(ref2040_sum, na.rm = FALSE)
colSums(ref2080_sum, na.rm = FALSE)


#find transition river km 
sum_2040 <- as.data.frame(c(1:6))
sum_2080 <- as.data.frame(c(1:6))

mods <- c("CCSM4_rcp85", "CNRM_rcp85", "FGOALS_rcp85","MIROC_rcp85")

for (o in mods) {
  #o="CCSM4_rcp85"
print(o)

tran_2040 <- raster(paste0('tran_2040_',o,'.tif'))
tran_2080 <- raster(paste0('tran_2080_',o,'.tif'))

#tran_2040[tran_2040 >= 0] <- 1

masked <- mask(x = tran_2040, mask = sierras_t)
tran_2040s <- raster::crop(masked, extent(sierras_t))
#tran_2040s[tran_2040s < 1] <- NA

tran_2040r<- rast(tran_2040s)
test <- cellSize(tran_2040r, unit ="km")
tran2040_sum<- exact_extract(test, sierras_t, fun="sum")
tran2040_sum<- as.data.frame(tran2040_sum)
sum_2040 <- cbind(sum_2040, tran2040_sum)

masked <- mask(x = tran_2080, mask = sierras_t)
tran_2080s <- raster::crop(masked, extent(sierras_t))
#tran_2080s[tran_2080s < 1] <- NA

tran_2080r<- rast(tran_2080s)
test <- cellSize(tran_2080r, unit ="km")
tran2080_sum<- exact_extract(test, sierras_t, fun="sum")
tran2080_sum<- as.data.frame(tran2080_sum)
sum_2080 <- cbind(sum_2080, tran2080_sum)
}
colSums(sum_2040, na.rm = FALSE)
colSums(sum_2080, na.rm = FALSE)
#mean
tran_2040 <- raster(paste0('tran_2040_mean.tif'))
tran_2080 <- raster(paste0('tran_2080_mean.tif'))

masked <- mask(x = tran_2040, mask = sierras_t)
tran_2040s <- raster::crop(masked, extent(sierras_t))

tran_2040r<- rast(tran_2040s)
test <- cellSize(tran_2040r, unit ="km")
tran2040_sum<- exact_extract(test, sierras_t, fun="sum")
tran2040_sum<- as.data.frame(tran2040_sum)

masked <- mask(x = tran_2080, mask = sierras_t)
tran_2080s <- raster::crop(masked, extent(sierras_t))

tran_2080r<- rast(tran_2080s)
test <- cellSize(tran_2080r, unit ="km")
tran2080_sum<- exact_extract(test, sierras_t, fun="sum")
tran2080_sum<- as.data.frame(tran2080_sum)

colSums(tran2040_sum, na.rm = FALSE)
colSums(tran2080_sum, na.rm = FALSE)

###############################################
#find km of suitable habitat in merced/tuolumne
MT <- readOGR('D:/Frog_datasets/Merced_Tuolumne_dis.shp')
MT_t <- spTransform(MT, CRS("+proj=longlat +datum=WGS84"))
#MT_t@bbox <- as.matrix(extent(suit_2020))

suit_2020<- raster('hs_2020_10092023.tif')
suit_2020[suit_2020 < 0.5] <- NA
suit_2020[suit_2020 >= 0.5] <- 1

masked <- mask(x = suit_2020, mask = MT_t)
hs_2020s <- raster::crop(masked, extent(MT_t))
hs_2020s[hs_2020s < 1] <- NA

hs_2020r<- rast(hs_2020s)
test <- cellSize(hs_2020r, unit ="km")
#test_sum <- sum(test, na.rm=TRUE)
hs2020_sum<- exact_extract(test, MT_t, fun="sum")


plot(hs_2020s)
lines(MT_t)

hs_mean_2040<- raster('hs_mean_2040.tif')
hs_mean_2080<- raster('hs_mean_2080.tif')

hs_mean_2040[hs_mean_2040 < 0.5] <- NA
hs_mean_2040[hs_mean_2040 >= 0.5] <- 1

masked <- mask(x = hs_mean_2040, mask = MT_t)
hs_mean_2040s <- raster::crop(masked, extent(MT_t))
hs_mean_2040s[hs_mean_2040s < 1] <- NA

hs_mean_2040r<- rast(hs_mean_2040s)
test <- cellSize(hs_mean_2040r, unit ="km")
hsmean_2040_sum<- exact_extract(test, MT_t, fun="sum")
hsmean_2040_sum<- as.data.frame(hsmean_2040_sum)

hs_mean_2080[hs_mean_2080 < 0.5] <- NA
hs_mean_2080[hs_mean_2080 >= 0.5] <- 1

masked <- mask(x = hs_mean_2080, mask = MT_t)
hs_mean_2080s <- raster::crop(masked, extent(MT_t))
hs_mean_2080s[hs_mean_2080s < 1] <- NA

hs_mean_2080r<- rast(hs_mean_2080s)
test <- cellSize(hs_mean_2080r, unit ="km")
hsmean_2080_sum<- exact_extract(test, MT_t, fun="sum")
hsmean_2080_sum<- as.data.frame(hsmean_2080_sum)

colSums(hsmean_2040_sum, na.rm = FALSE)
colSums(hsmean_2080_sum, na.rm = FALSE)


#find refugia river km 
#mean
ref_2040 <- raster(paste0('ref_2040_mean.tif'))
ref_2080 <- raster(paste0('ref_2080_mean.tif'))

masked <- mask(x = ref_2040, mask = MT_t)
ref_2040s <- raster::crop(masked, extent(MT_t))

ref_2040r<- rast(ref_2040s)
test <- cellSize(ref_2040r, unit ="km")
ref2040_sum<- exact_extract(test, MT_t, fun="sum")
ref2040_sum<- as.data.frame(ref2040_sum)


masked <- mask(x = ref_2080, mask = MT_t)
ref_2080s <- raster::crop(masked, extent(MT_t))

ref_2080r<- rast(ref_2080s)
test <- cellSize(ref_2080r, unit ="km")
ref2080_sum<- exact_extract(test, MT_t, fun="sum")
ref2080_sum<- as.data.frame(ref2080_sum)



#find transition river km 
#mean
tran_2040 <- raster(paste0('tran_2040_mean.tif'))
tran_2080 <- raster(paste0('tran_2080_mean.tif'))

masked <- mask(x = tran_2040, mask = MT_t)
tran_2040s <- raster::crop(masked, extent(MT_t))

tran_2040r<- rast(tran_2040s)
test <- cellSize(tran_2040r, unit ="km")
tran2040_sum<- exact_extract(test, MT_t, fun="sum")
tran2040_sum<- as.data.frame(tran2040_sum)

masked <- mask(x = tran_2080, mask = MT_t)
tran_2080s <- raster::crop(masked, extent(MT_t))

tran_2080r<- rast(tran_2080s)
test <- cellSize(tran_2080r, unit ="km")
tran2080_sum<- exact_extract(test, MT, fun="sum")
tran2080_sum<- as.data.frame(tran2080_sum)

colSums(tran2040_sum, na.rm = FALSE)
colSums(tran2080_sum, na.rm = FALSE)

