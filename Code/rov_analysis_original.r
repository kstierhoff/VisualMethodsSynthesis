# Code for analyzing data from ROV surveys during the comparative rockfish survey (COAST 11)
# 8 July 2012	Kevin L. Stierhoff (NOAA Fisheries-SWFSC)

rm(list=ls())
# set system time zone to GMT
Sys.setenv(TZ='GMT')
# load libraries
library(RODBC);library(reshape2);library(plyr);library(ggplot2);library(vegan);
library(maps);library(scales);library(forecast)
# determine the OS for setting the working directory
if(.Platform$OS.type == "unix") {
	root.dir <- "/Users/kevinstierhoff/Dropbox/R/coast11/"
} else {
	root.dir <- "C:/Users/kls/Documents/Projects/2011 Comparative Rockfish Survey/"
}
# define working directories and file paths
dbdir  <- "C:/Users/kls/Documents/Data/rov_data/"
dbname <- "ROV_Master.accdb"
wkdir  <- root.dir
FOVdir <- paste(root.dir,"3Beam_Logs",sep="")

# #read data from ROV database
# setwd(dbdir)
# channel		<- odbcConnectAccess2007(dbname)
# nav 		<- sqlFetch(channel,"dbo_tbl_NAV_DATA_COAST11")
# slope 		<- sqlFetch(channel,"dbo_tbl_NAV_SLOPE_COAST11")
# depth 		<- sqlFetch(channel,"dbo_tbl_NAV_DEPTH_COAST11")
# obs 		<- sqlFetch(channel,"dbo_tbl_VIDEO_OBS_COAST11")
# spp 		<- sqlFetch(channel,"dbo_tlu_SPECIES_CODES")
# site.names 	<- sqlFetch(channel,"dbo_tlu_SITE_GEN_COAST11") 
# lw.data 	<- sqlFetch(channel,"dbo_tlu_LENGTH_WEIGHT_DATA") 
# close(channel)

# setwd(wkdir)
# # read survey area data (use survey area defined by Yok and Watters)
# survey.area <- read.csv("survey_area_Yok_coast11.txt") 
# # save data
# save(nav,slope,depth,obs,spp,site.names,survey.area,lw.data,file = "coast11.Rdata")

#load data from file
setwd(wkdir)
load("coast11.Rdata")
# read good nav data from GIS analysis
good.nav <- read.csv("good_nav_COAST11.csv")
nav <- subset(nav,nav$nav_id %in% good.nav$nav_id)
nav <- nav[order(nav$nav_id),]
#add slope data to nav
nav <- merge(nav,slope[,c("nav_id","slope_dem")],by = "nav_id")
nav <- merge(nav,depth[,c("nav_id","depth_dem")],by = "nav_id")
nav <- merge(nav,site.names,by="dive_name")
# add species names to obs
obs <- droplevels(merge(obs,spp[c("species_code","sci_name_full","common_name","genus_t")],by = "species_code",all.x = TRUE))
obs <- merge(obs,site.names,by="dive_name")
# subset species list to include only those in the video tapes
spp <- droplevels(subset(spp[ ,c("species_code","sci_name","sci_name_full","common_name","genus_t","subgenus_t","species_t")],
	spp$species_code %in% unique(obs$species_code)))
spp <- spp[order(spp$species_code), ]
spp[order(spp$sci_name_full),c("species_code","sci_name_full","common_name")]
# subset L/W data to include only species of interest
lw.data <- merge(lw.data,spp[,c("sci_name_full","genus_t")],by = "sci_name_full")
lw.data <- droplevels(subset(lw.data, lw.data$genus_t %in% c("Sebastes","Sebastomus","Sebastolobus","Chromis","Caulolatilus",
	"Oxyjulis","Trachurus","Merluccius","Ophiodon") & lw.data$to_use == 1))
lw.table <- merge(lw.data,spp[ ,c("sci_name_full","common_name")],by = "sci_name_full")
write.csv(lw.table[ ,c("sci_name_full","common_name","a","b","sex","reference","comment")],file = "lw_table.csv",row.names = FALSE,na="")

# #############################################################################################
# Process nav data #####
# #############################################################################################
# get list of dives from nav data
dives.nav <- unique(nav$dive_name)	
#create empty data frame to store smoothed data	
nav.smooth <- data.frame()
# clean CTD and nav data from the Phantom ROV
pb1 <- txtProgressBar(min = 0, max = length(dives.nav), style = 3)
for (i in 1:length(dives.nav)){
# subset data
	temp.nav <- subset(nav[order(nav$nav_id),],nav$dive_name == dives.nav[i])
# fit a Loess smoother to depth data
	depth.loess <- loess(depth~nav_id,temp.nav, span = 0.02)
	# plot(temp.nav$depth)
	# lines(depth.loess$fitted,col = 'red',lwd = 2)
	temp.nav$depth <- depth.loess$fitted
# fit a Loess smoother to temperature data
	temp.loess <- loess(temperature~nav_id,temp.nav, span = 0.1)
	# plot(temp.nav$date_time,temp.nav$temperature, type = "l")
	# lines(temp.nav$date_time,temp.loess$fitted,col = 'red',lwd = 2)
	temp.nav$temperature <- temp.loess$fitted
# fit a Loess smoother to conductivity data
	cond.loess <- loess(conductivity~nav_id,temp.nav, span = 0.025)
	# plot(temp.nav$date_time,temp.nav$conductivity, type = "l")
	# lines(temp.nav$date_time,cond.loess$fitted,col = 'red',lwd = 2)
	temp.nav$conductivity <- cond.loess$fitted
# fit a Loess smoother to pressure data
	press.loess <- loess(pressure~nav_id,temp.nav, span = 0.025)
	# plot(temp.nav$date_time,temp.nav$pressure, type = "l")
	# lines(temp.nav$date_time,press.loess$fitted,col = 'red',lwd = 2)
	temp.nav$pressure <- press.loess$fitted
# fit a Loess smoother to salinity data
	sal.loess <- loess(salinity~nav_id,temp.nav, span = 0.025)
	# plot(temp.nav$date_time,temp.nav$salinity, type = "l")
	# lines(temp.nav$date_time,sal.loess$fitted,col = 'red',lwd = 2)
	temp.nav$salinity <- sal.loess$fitted
# fit a Loess smoother to sound velocity data
	sound.loess <- loess(sound_vel~nav_id,temp.nav, span = 0.02)
	# plot(temp.nav$date_time,temp.nav$sound_vel, type = "l")
	# lines(temp.nav$date_time,sound.loess$fitted,col = 'red',lwd = 2)
	temp.nav$sound_vel <- sound.loess$fitted
# fit a Loess smoother to oxygen conc data
	oc.loess <- loess(oxygen_conc~nav_id,temp.nav, span = 0.02)
	# plot(temp.nav$date_time,temp.nav$oxygen_conc, type = "l")
	# lines(temp.nav$date_time,oc.loess$fitted,col = 'red',lwd = 2)
	temp.nav$oxygen_conc <- oc.loess$fitted
# fit a Loess smoother to oxygen sat data
	os.loess <- loess(oxygen_sat~nav_id,temp.nav, span = 0.02)
	# plot(temp.nav$date_time,temp.nav$oxygen_sat, type = "l")
	# lines(temp.nav$date_time,os.loess$fitted,col = 'red',lwd = 2)
	temp.nav$oxygen_sat <- os.loess$fitted
# fit a Loess smoother to altitude data
	alt.loess <- loess(altitude~nav_id,temp.nav, span = 0.01)
	# plot(temp.nav$date_time,temp.nav$altitude, type = "l")
	# lines(temp.nav$date_time,alt.loess$fitted,col = 'red',lwd = 2)
	temp.nav$altitude <- alt.loess$fitted
# # fit a Loess smoother to lat/lon data
	# lat.loess <- loess(lat~nav_id,temp.nav, span = 0.005)
	# lon.loess <- loess(lon~nav_id,temp.nav, span = 0.005)
	# plot(temp.nav$lon,temp.nav$lat, type = "l")
	# lines(lon.loess$fitted,lat.loess$fitted,col = 'red',lwd = 2)
	# temp.nav$lat <- lat.loess$fitted
	# temp.nav$lon <- lon.loess$fitted
# # fit a Loess smoother to northing/easting data
	# n.loess <- loess(northing_r~nav_id,temp.nav, span = 0.005)
	# e.loess <- loess(easting_r~nav_id,temp.nav, span = 0.005)
	# # plot(temp.nav$easting_r,temp.nav$northing_r, type = "l")
	# # lines(e.loess$fitted,n.loess$fitted,col = 'red',lwd = 2)
	# temp.nav$northing_r <- n.loess$fitted
	# temp.nav$easting_r <- e.loess$fitted
# fit a Loess smoother to ROV speed data
	speed.loess <- loess(speed~nav_id,temp.nav, span = 0.01)
	# plot(temp.nav$date_time,temp.nav$speed_r, type = "l")
	# lines(temp.nav$date_time,speed.loess$fitted,col = 'red',lwd = 2)
	temp.nav$speed_r <- speed.loess$fitted
# add present dive's data to the temporary nav data frame
	nav.smooth <- rbind(nav.smooth,temp.nav)
# update progress bar
	setTxtProgressBar(pb1, i)
}
close(pb1)
# replace nav with smoothed data
nav <- nav.smooth

# create vector with unique dive names; used for subsetting larger datasets
dives.obs <- droplevels(sort(unique(obs$dive_name)))
# remove nav_id, lag_s, and slope_dem from data frame
obs <- obs[,-which(names(obs) %in% c("nav_id","lag_s","slope_dem"))]
#create geol_ps varible for nav and obs
obs$geol_ps <- as.factor(paste(obs$geol_prim,obs$geol_sec,sep = ""))
# subset obs to include time 'on effort', select species, and dive segments
obs.sub <- droplevels(subset(obs, obs$genus_t %in% c("Sebastes","Sebastomus","Sebastolobus","Chromis","Caulolatilus",
	"Oxyjulis","Trachurus","Merluccius","Ophiodon")))

# summarize species abundance for all dives
spp.summ <- as.data.frame(t(sapply(sort(unique(obs.sub$species_code)),function(x){
	data <- subset(obs.sub,species_code == x)
	species <- as.character(data$species_code[1])
	total.obs <- as.numeric(sum(data$counts, na.rm = TRUE))
	c(species_code = species, sum = total.obs)
})))
spp.summ$sum <- as.numeric(levels(spp.summ$sum)[spp.summ$sum])
spp.summ$pct.total <- spp.summ$sum/sum(spp.summ$sum)*100 
spp.summ <- merge(spp.summ,spp[,c("species_code","sci_name_full","common_name")],by = "species_code", all.x = TRUE)
spp.summ <- spp.summ[order(spp.summ$sum, decreasing = TRUE),]
# print results to screen
spp.summ[,c("sci_name_full","common_name","sum","pct.total")]
write.csv(spp.summ[,c("sci_name_full","common_name","sum","pct.total")],file = "species_obs_summary.csv",row.names = FALSE)

# restrict nav data to dives where video data are present
nav <- droplevels(subset(nav,nav$dive_name %in% dives.obs))

# #############################################################################################
# summarize nav data by dive #####
# #############################################################################################
nav.summ <- as.data.frame(t(sapply(sort(unique(nav$dive_name)),function(x){
	data 		<- subset(nav,dive_name == x)
	dive 		<- as.character(data$dive_name[1])
	lat.start 	<- data$lat[1]
	lon.start 	<- data$lon[1]
	lat.end		<- data$lat[length(data$lat)]
	lon.end		<- data$lon[length(data$lat)]
	date.start 	<- as.character(data$date_time[1])
	date.end 	<- as.character(data$date_time[length(data$date_time)])
	z.mean 		<- mean(abs(data$depth_dem))
	z.min 		<- min(abs(data$depth_dem))
	z.max 		<- max(abs(data$depth_dem))
	z.range 	<- max(abs(data$depth_dem)) - min(abs(data$depth_dem))
	slope.mean 	<- mean(abs(data$slope_dem))
	distance 	<- sum(data$disp)
	c(dive.name = dive, lat.start=lat.start,lon.start=lon.start,lat.end=lat.end,lon.end=lon.end,
		date.start=date.start,date.end=date.end,z.mean = z.mean, z.min = z.min, z.max = z.max, 
		z.range = z.range,slope.mean = slope.mean, distance=distance)
	# c(dive.name = dive, z.mean = z.mean, z.min = z.min, z.max = z.max, z.range = z.range, slope.mean = slope.mean, disp.sum = disp.sum)
})))
nav.summ
write.csv(nav.summ,file = "nav_summary.csv",row.names = FALSE)

# #############################################################################################
# summarize species data by dive #####
# #############################################################################################
obs.sub$key.summ <- paste(obs.sub$dive_name,obs.sub$sci_name_full)
obs.summ <- as.data.frame(t(sapply(sort(unique(obs.sub$key.summ)),function(x){
	data 		<- subset(obs.sub,key.summ == x)
	dive 		<- as.character(data$dive_name[1])
	species		<- as.character(data$sci_name_full[1])
	counts	 	<- sum(data$counts)
	c(dive.name = dive, species=species,counts=counts)
	# c(dive.name = dive, z.mean = z.mean, z.min = z.min, z.max = z.max, z.range = z.range, slope.mean = slope.mean, disp.sum = disp.sum)
})))
obs.summ

# #############################################################################################
# merge summarized species data by dive #####
# #############################################################################################
rov.summary <- merge(obs.summ,nav.summ,by="dive.name")
# write data to file for Rick Methot et al.
write.csv(rov.summary,file = "obs_summary_by_dive.csv",row.names = FALSE)

# #############################################################################################
# ADD NAV DATA TO VIDEO OBSERVATIONS #####
# #############################################################################################
#create variables
video_id   	<- numeric()
nav_id 		<- video_id
lag_s    	<- video_id
#create status bar
pb1 <- txtProgressBar(min = 0, max = length(dives.obs), style = 3)
for (i in 1:length(dives.obs)){
# get nav ID and time lag for each video observation
# subset nav and obs data to only include values for one dive at a time
temp.nav <- subset(nav[ ,c("nav_id","date_time")],nav$dive_name == dives.obs[i])
temp.obs <- subset(obs[ ,c("video_id","date_time")],obs$dive_name == dives.obs[i])
# create temporary variables to hold sync results
temp.vid.id 		<- numeric()
temp.vid.nav.id 	<- temp.vid.id
temp.vid.lag.s		<- temp.vid.id
	# get habitat info for each nav record
	for (j in 1:length(temp.obs$date_time)){
		# calculate the time difference between jth video obs and each nav record
		time.diff 		<- abs(difftime(temp.obs$date_time[j],temp.nav$date_time,units = "secs"))
		temp.vid.id 	<- c(temp.vid.id,temp.obs$video_id[j])
		temp.vid.nav.id <- c(temp.vid.nav.id,temp.nav$nav_id[which.min(time.diff)])
		temp.vid.lag.s 	<- c(temp.vid.lag.s,min(time.diff))
	}
	hist(temp.vid.lag.s, xlab = "Time lag (s)", main = dives.obs[i])
	# hist(temp.vid.lag.s)
	# add sync results from jth dive to previous sync results
	video_id 	<- c(video_id,temp.vid.id)
	nav_id 		<- c(nav_id,  temp.vid.nav.id)
	lag_s 		<- c(lag_s,   temp.vid.lag.s)
	
	# update progress bar
	setTxtProgressBar(pb1, i)
}
close(pb1)
# store video sync results in a data frame
vid.nav.sync <- data.frame(video_id,nav_id,lag_s)
str(vid.nav.sync)
hist(vid.nav.sync$lag_s, breaks = seq(0,max(vid.nav.sync$lag_s),max(vid.nav.sync$lag_s)/30))

# #############################################################################################
# add nav data to video observations #####
# #############################################################################################
str(obs)
obs <- droplevels(merge(obs,vid.nav.sync,by = "video_id",all.x = TRUE))
# add nav data to obs
obs <- merge(obs,nav[,c("nav_id","depth","temperature","salinity","oxygen_conc","oxygen_sat","depth_dem","slope_dem","lat","lon")],by="nav_id")
write.csv(obs[,c("video_id","nav_id","lag_s","slope_dem","depth_dem")],file = "obs_nav_sync_data.csv",row.names = FALSE)
# write observation data to csv file for Juan Zwolinski (for COAST analysis)
write.csv(obs[obs$genus %in% c("Sebastes","Sebastomus","Sebastolobus","Chromis","Caulolatilus","Oxyjulis","Trachurus","Merluccius","Ophiodon"),
	c("video_id","dive_name","site_name","date_time","sci_name_full","counts","lat","lon","depth")],file="coast11_video_obs.csv",row.names=FALSE)
# write file for plotting in GIS
write.csv(obs[,c("video_id","dive_name","sci_name_full","counts","size_code","geol_prim","lat","lon")],file = "all_observations.csv",row.names=FALSE)
# quick plot of obs positions
# qplot(lon,lat,data = nav,colour = depth_dem)
# qplot(lon,lat,data = obs,colour = depth_dem)

# #############################################################################################
# ADD VIDEO OBSERVATIONS TO NAV DATA #####
# #############################################################################################
# find unique date/times to cut nav data
obs.unique <- subset(obs,!duplicated(obs$date_time) == TRUE)
#create variable for video_id
video_id   		<- numeric()
nav_count 		<- numeric()
obs_id_count	<- numeric()
#create status bar
pb1 <- txtProgressBar(min = 0, max = length(dives.obs), style = 3)
for (i in 1:length(dives.obs)){
	#subset nav and obs data to only include values for one dive at a time
	temp.nav <- subset(nav[ ,c("nav_id","date_time")],nav$dive_name == dives.obs[i])
	temp.obs <- subset(obs.unique[ ,c("video_id","date_time")],obs.unique$dive_name == dives.obs[i])
	cuts <- cut(as.numeric(temp.nav$date_time),unique(as.numeric(temp.obs$date_time)), labels = FALSE,include.lowest = TRUE, right = FALSE)
	obs_id_count <- c(obs_id_count,length(temp.obs$video_id[cuts]))
	nav_count 	 <- c(nav_count,nrow(temp.nav))
	video_id  	 <- c(video_id,temp.obs$video_id[cuts])
	# update progress bar
	setTxtProgressBar(pb1, i)
}
close(pb1)
nav$video_id <- video_id
nav <- merge(nav,obs[,c("video_id","search_code","geol_prim","geol_sec")],by = "video_id", all.x = TRUE, sort = FALSE)
nav$geol_ps <- as.factor(paste(nav$geol_prim,nav$geol_sec,sep = ""))
# remove nav data collected before and after video analysis began
nav <- droplevels(subset(nav,is.na(nav$geol_prim)==FALSE))
nav <- nav[order(nav$nav_id), ]
# write nav data to csv file for Cutter
write.csv(nav[,c("video_id","nav_id","object_id","dive_name","date_time","lat","lon","depth",
"site_name","geol_prim","geol_sec","geol_ps")],row.names=FALSE,file="coast11_seabed_obs.csv")
# quick plot lat/lon for nav and obs
# qplot(lon,lat,data = nav,colour = geol_prim)
# qplot(lon,lat,data = obs,colour = geol_prim)

# #############################################################################################
#create variable depth.bin to cut survey effort and observations by depth for each dive #####
depth.bin <- seq(-600,0,100)	#create depth bins
depth.lab <- seq(-550,50,100)	#create depth labels (midpoints of depth bins, in this case)

# #############################################################################################
#create depth_stratum factor for nav and obs records #####
nav$depth_stratum <- cut(nav$depth_dem,depth.bin)
nav$key.z <- as.factor(paste(nav$dive_name,nav$depth_stratum))
# create depth_stratum factor for obs records
obs$depth_stratum <- cut(obs$depth_dem,depth.bin)
obs$key.z <- as.factor(paste(obs$dive_name,obs$depth_stratum))

# #############################################################################################
#calculate distance within each dive and depth stratum #####
nav.depth.cut <- tapply(nav$disp,list(nav$dive_name,nav$depth_stratum),sum)
nav.hab.cut <- tapply(nav$disp,list(nav$dive_name,nav$geol_prim),sum)
#combine primary and secondary geol classes
nav$geol_ps <- paste(nav$geol_prim,nav$geol_sec,sep = "")
nav.geolps.cut <- tapply(nav$disp,list(nav$dive_name,nav$geol_ps),sum)
#melt nav.cut
nav.melt <- melt(nav.depth.cut)
names(nav.melt) <- c("dive_name","depth_stratum","disp.sum")
nav.melt$key.z <- as.factor(paste(nav.melt$dive_name,nav.melt$depth_stratum))
#create dive list for transects greater than 200m
nav.list <- droplevels(subset(nav.melt,nav.melt$disp.sum > 100))
#select nav data from dive list
nav.sub <- droplevels(subset(nav,nav$key.z %in% nav.list$key.z))
# summarize nav data by site and depth
nav.summ.site <- dcast(melt(nav.sub[,c("site_name","depth_stratum","disp")]),site_name+depth_stratum~.,sum)
names(nav.summ.site) <- c("site_name","depth_stratum","disp")
write.csv(nav.summ.site,file = "nav_summary_site.csv",row.names = FALSE)
# summarize nav data by transect and site
nav.summ.trans <- dcast(melt(nav.sub[,c("site_name","key.z","depth_stratum","disp")]),key.z+site_name+depth_stratum~.,sum)
names(nav.summ.trans) <- c("transect","site_name","depth_stratum","disp")
# summarize depth by transect
nav.z.trans <- dcast(melt(nav.sub[,c("site_name","key.z","depth_stratum","depth")]),key.z+site_name+depth_stratum~.,mean)
names(nav.z.trans) <- c("transect","site_name","depth_stratum","depth")
write.csv(nav.summ.trans,file = "nav_summary_site.csv",row.names = FALSE)
# summarize nav data by transect and site
nav.summ.hab <- dcast(melt(nav.sub[,c("site_name","geol_prim","disp")]),site_name~geol_prim,sum)
write.csv(nav.summ.hab,file = "nav_summary_hab.csv",row.names = FALSE)

# #############################################################################################
# ANALYSIS OF FIELD OF VIEW FOR ALL TRANSECTS #####
# #############################################################################################
setwd(FOVdir)
dir <- dir(pattern = "*_ProcessingLog.txt", ignore.case = TRUE)
cat("Available 3Beam analyses:",length(dir),"files \n")
dir
#create output data.frame
fov.output <- data.frame()
temp.fov.id <- 1

for (i in 1:length(dir)){
cat("Processing",i,"of",length(dir),"; SERIES name:",dir[i],"\n")
flush.console()
data <- read.csv(dir[i])
name <- substr(dir[i],1,7)						# extract the dive name from the file name
series <- substr(dir[i],1,nchar(dir[i])-18)		# extract the series name from the file name

#convert factors to character
data$ImageFilename 	<- as.character(data$ImageFilename)
data$ImageDate 		<- as.character(data$ImageDate)
data$ImageTime 		<- as.character(data$ImageTime)
data$NavDate 		<- as.character(data$NavDate)
data$NavTime 		<- as.character(data$NavTime)
#create file paths for images
data$ImageFilename 	<- substr(data$ImageFilename,nchar(data$ImageFilename)-12,nchar(data$ImageFilename))
data$ImagePath 		<- paste("\\\\abulon\\ROV\\Data\\3Beam_data\\3Beam_ImageData\\",series,"\\",data$ImageFilename, sep = "")
#convert date/time fields to POSIX (on which you can do math)
data$img_datetime 	<- as.POSIXct(paste(data$ImageDate,data$ImageTime),format = "%m/%d/%Y %H:%M:%S") 	#"10/08/2010 18:29:03"
data$nav_datetime 	<- as.POSIXct(paste(data$NavDate,data$NavTime),format = "%m/%d/%Y %H:%M:%S") 		#"10/08/2010 18:29:03"
data$ImageDate 		<- NULL
data$ImageTime 		<- NULL
data$NavDate 		<- NULL
data$NavTime 		<- NULL
#clean unused data
data$ImageFrame 	<- NULL
data$NavRoll 		<- NULL
data$NavPitch 		<- NULL
data$NavTrueHead 	<- NULL
data$NavAltitude 	<- NULL
data$NavNorthVel 	<- NULL
data$NavEastVel 	<- NULL
data$NavDepth 		<- NULL
data$NavLong 		<- NULL
data$NavLat 		<- NULL
data$ProcessArea 	<- NULL
#rename RGB data
data$LeftR <- data$R; data$R 		<- NULL
data$LeftG <- data$G; data$G 		<- NULL
data$LeftB <- data$B; data$B 		<- NULL
data$RightR <- data$R.1; data$R.1 	<- NULL
data$RightG <- data$G.1; data$G.1 	<- NULL
data$RightB <- data$B.1; data$B.1 	<- NULL
data$CrossR <- data$R.2; data$R.2 	<- NULL
data$CrossG <- data$G.2; data$G.2 	<- NULL
data$CrossB <- data$B.2; data$B.2 	<- NULL

# # plot results
# plot(data$CenterWidth, type = "l", xlab = "Image number",ylab ="Center Width")

# compute distance at each time step
data$dist <- rep(0,length(data$NavMeanDist))
data$dist[2:end(data$dist)[1]] <- data$NavMeanDist[2:end(data$NavMeanDist)[1]]-data$NavMeanDist[1:end(data$NavMeanDist)[1]-1]
plot(data$dist,type = "l",xlab = "Image number",ylab ="Distance")

# replace unfound lasers with the mean center width
# make a copy of original data
data.corr <- data
# find bad data or outliers in CenterWidth(CenterWidth >8m)
data.corr$CenterWidth[data.corr$CenterWidth == 0 | data.corr$CenterWidth > 8] <- NA
# interpolate missing data
data.corr$CenterWidth <- na.interp(data.corr$CenterWidth)
data.corr$CenterWidth[is.na(data.corr$CenterWidth)] <- mean(data.corr$CenterWidth, na.rm = T)
# comparison plot of results
layout(matrix(c(1,1,2,3),2,2,byrow = TRUE))
plot(data$CenterWidth, main = c("3Beam Analysis Results:",series),ylab = "Center width (m)", pch = 20, ylim = c(0,10))
par(new = T)
plot(data.corr$CenterWidth,ylab = "",ylim = c(0,10), type = "l", col = "red")
plot(data.corr$LeftX,data.corr$LeftY,col = "red",pch = 20, xlim = c(150,500),ylim = c(225,325), 
	xlab = "X Coordinate",ylab = "Y Coordinate", main = "Laser Coordinates")
par(new = T)
plot(data.corr$RightX,data.corr$RightY,col = "blue",pch = 20, xlim = c(150,500),ylim = c(225,325), 
	xlab = "",ylab = "", main = "")
par(new = T)
plot(data.corr$CrossX,data.corr$CrossY,col = "green",pch = 20, xlim = c(150,500),ylim = c(225,325), 
	xlab = "",ylab = "", main = "")
hist(data.corr$CenterWidth,breaks = seq(0,10,0.2), main = c("Mean Center Width:",format(mean(data.corr$CenterWidth),digits = 3)),
	xlab = "Center width(m)")
#save plot of processing results as JPG
dev.copy(jpeg,paste("FOV_plot_",series,".jpg", sep = ""))
dev.off()
#calculate mean fields of view for raw and corrected results
raw.mean <- mean(data$CenterWidth, na.rm = T)
corr.mean <- mean(data.corr$CenterWidth)
#print results
cat("Results for series:",series,"\n")
cat("Raw mean:",raw.mean,",Corrected mean:",corr.mean,"\n")
#str(data.corr)
#create index and series name
id <- seq(temp.fov.id,temp.fov.id+(length(data.corr$CenterWidth)-1),1)
#create new temp.fov.id for next series
temp.fov.id <- max(id) + 1
# create dive and series name vectors
dive_name <- rep(name,length(data.corr$CenterWidth))
series_name <- rep(series,length(data.corr$CenterWidth))
#assemble results
temp.data <- data.frame(
	id, 
	dive_name,
	series_name,
	# nav_id = rep(NA,length(data.corr$CenterWidth)),
	# lag_s = rep(NA,length(data.corr$CenterWidth)),
	filepath = data.corr$ImagePath,
	date_time = data.corr$img_datetime,
	found_lasers = data.corr$FoundLasers,
	center_width = data.corr$CenterWidth,
	slant_range = rep(-999,length(data.corr$CenterWidth)),
	area_below_center = rep(-999,length(data.corr$CenterWidth)),
	left_x = data.corr$LeftX,
	left_y = data.corr$LeftY,
	left_r = data.corr$LeftR,
	left_g = data.corr$LeftG,
	left_b = data.corr$LeftB,
	right_x = data.corr$RightX,
	right_y = data.corr$RightY,
	right_r = data.corr$RightR,
	right_g = data.corr$RightG,
	right_b = data.corr$RightB,
	cross_x = data.corr$CrossX,
	cross_y = data.corr$CrossY,
	cross_r = data.corr$CrossR,
	cross_g = data.corr$CrossG,
	cross_b = data.corr$CrossB,
	disp = data.corr$dist,
	total_disp = data.corr$NavMeanDist
	
	)
#add results of series i with previously analyzed series
fov.output <- rbind(fov.output,temp.data)
}
str(fov.output)
unique(fov.output$series_name)
mean(fov.output$center_width)

# #############################################################################################
# ADD NAV DATA TO FOV ESTIMATES #####
# #############################################################################################
#create variables
fov_id   	<- numeric()
nav_id 		<- fov_id
lag_s    	<- fov_id
# create vector with unique dive names; used for subsetting larger datasets
dives.fov <- as.character(droplevels(sort(unique(fov.output$dive_name))))
#create status bar
pb1 <- txtProgressBar(min = 0, max = length(dives.fov), style = 3)
for (i in 1:length(dives.fov)){
# get nav ID and time lag for each video observation
# subset nav and obs data to only include values for one dive at a time
temp.nav <- subset(nav[ ,c("nav_id","date_time")],as.character(nav$dive_name) %in% dives.fov[i])
temp.fov <- subset(fov.output[ ,c("id","date_time")],as.character(fov.output$dive_name) %in% dives.fov[i])
# create temporary variables to hold sync results
temp.fov.id 		<- numeric()
temp.fov.nav.id 	<- temp.fov.id
temp.fov.lag.s		<- temp.fov.id
	# get habitat info for each nav record
	for (j in 1:length(temp.fov$date_time)){
	# calculate the time difference between jth video obs and each nav record
	time.diff 		<- abs(difftime(temp.fov$date_time[j],temp.nav$date_time,units = "secs"))
	temp.fov.id 	<- c(temp.fov.id,temp.fov$id[j])
	temp.fov.nav.id <- c(temp.fov.nav.id,temp.nav$nav_id[which.min(time.diff)])
	temp.fov.lag.s 	<- c(temp.fov.lag.s,min(time.diff))
	}
	layout(1)
	hist(temp.fov.lag.s, xlab = "Time lag (s)", main = dives.fov[i])
	# hist(temp.fov.lag.s)
	# add sync results from jth dive to previous sync results
	fov_id 	<- c(fov_id,temp.fov.id)
	nav_id 	<- c(nav_id,temp.fov.nav.id)
	lag_s 	<- c(lag_s,temp.fov.lag.s)
	
	# update progress bar
	setTxtProgressBar(pb1, i)
}
close(pb1)
# store video sync results in a data frame
fov.nav.sync <- data.frame(fov_id,nav_id,lag_s)
# rename variables in data frame
names(fov.nav.sync) <- c("id","nav_id","lag_s")
hist(fov.nav.sync$lag_s, breaks = seq(0,max(fov.nav.sync$lag_s),max(fov.nav.sync$lag_s)/30))
# add nav data to FOV estimates
fov.output <- droplevels(merge(fov.output,fov.nav.sync,by = "id"))

# #############################################################################################
# summarize FOV results and save to file #####
# #############################################################################################
setwd(wkdir)
#format date for database import
fov.output$date_time <- format(fov.output$date_time, format = "%m/%d/%Y %H:%M:%S")
fov.output$dive_name <- as.character(fov.output$dive_name)
fov.output$series_name <- as.character(fov.output$series_name)
#write results to text file
write.csv(fov.output,"coast11_FOV.txt", row.names=FALSE,quote = FALSE, na = "")

#summarize results by series name
fov.summ.series <- t(sapply(sort(unique(fov.output$series_name)),function(x){
	fov.sub  <- subset(fov.output,series_name == x)
	# series 	 <- as.character(fov.sub$series_name[1])
	samples  <- length(fov.sub$series_name)
	fov.mean <- mean(fov.sub$center_width)
	fov.sd 	 <- sd(fov.sub$center_width)
	
	c(n_images = samples,mean = fov.mean, sd = fov.sd)
}))
fov.summ.series <- data.frame(fov.summ.series)
series <- unique(fov.output$series_name)
fov.summ.series <- cbind(series,fov.summ.series)
fov.summ.series
write.csv(fov.summ.series,"coast11_FOV_summ_series.csv", row.names=FALSE,quote = FALSE, na = "")

#summarize results by dive name
fov.summ.dive 	<- t(sapply(sort(unique(fov.output$dive_name)),function(x){
	fov.sub 	<- subset(fov.output,dive_name == x)
	samples 	<- length(fov.sub$dive_name)
	fov.mean 	<- mean(fov.sub$center_width)
	fov.sd 		<- sd(fov.sub$center_width)
	
	c(n_images = samples,mean = fov.mean, sd = fov.sd)
}))
fov.summ.dive <- data.frame(fov.summ.dive)
dive <- unique(fov.output$dive_name)
fov.summ.dive <- cbind(dive,fov.summ.dive)
fov.summ.dive
write.csv(fov.summ.dive,"coast11_FOV_summ_dive.csv", row.names=FALSE,quote = FALSE, na = "")

#save fov.output to Rdata file
save(fov.output,file = "coast11_FOV.Rdata")

# #############################################################################################
# ADD FOV DATA TO NAV #####
# #############################################################################################
#create variables
nav_id   	 <- numeric()
fov_id 		 <- nav_id
lag_s    	 <- nav_id
center_width <- nav_id

#create new data frame for syncing FOV data with nav
fov.sync <- fov.output
fov.sync$date_time 	<- as.POSIXct(fov.sync$date_time,format = "%m/%d/%Y %H:%M:%S") #09/21/2011 14:59:24

#create status bar
pb1 <- txtProgressBar(min = 0, max = length(dives.obs), style = 3)
for (i in 1:length(dives.obs)){
# get FOV ID and time lag for each nav record
# subset nav and FOV data to only include values for one dive at a time
temp.nav <- subset(nav[ ,c("nav_id","date_time")],nav$dive_name == dives.obs[i])
temp.fov <- subset(fov.sync[ ,c("id","date_time","center_width")],fov.output$dive_name == dives.obs[i])
# create temporary variables to hold sync results
temp.nav.id 	<- numeric()
temp.fov.id 	<- temp.nav.id
temp.fov.lag.s	<- temp.nav.id
	# get FOV info for each nav record
	for (j in 1:length(temp.nav$date_time)){
	# calculate the time difference between j-th nav record and each FOV estimate
	time.diff 		<- abs(difftime(temp.nav$date_time[j],temp.fov$date_time,units = "secs"))
	temp.nav.id 	<- c(temp.nav.id,temp.nav$nav_id[j])
	temp.fov.id 	<- c(temp.fov.id,temp.fov$id[which.min(time.diff)])
	temp.fov.lag.s 	<- c(temp.fov.lag.s,min(time.diff))
	}
	hist(temp.fov.lag.s, xlab = "Time lag (s)", main = dives.obs[i])
	# add sync results from j-th dive to previous sync results
	nav_id 		<- c(nav_id,temp.nav.id)
	fov_id 		<- c(fov_id,temp.fov.id)
	lag_s 		<- c(lag_s,temp.fov.lag.s)
	
	# update progress bar
	setTxtProgressBar(pb1, i)
}
close(pb1)
# store video sync results in a data frame
fov.nav.sync <- data.frame(nav_id,fov_id,lag_s)
names(fov.nav.sync) <- c("nav_id","id","fov_lag_s")
str(fov.nav.sync)
hist(fov.nav.sync$fov_lag_s, breaks = seq(0,max(fov.nav.sync$fov_lag_s),max(fov.nav.sync$fov_lag_s)/30))
# add center width data
fov.nav.sync <- merge(fov.nav.sync,fov.sync[ ,c("id","center_width")],by="id")
# add center width and fov lag data to nav
nav <- merge(nav,fov.nav.sync[ ,c("nav_id","fov_lag_s","center_width")],by="nav_id")
# save nav data to Rdata file
save(nav,file="nav_data.Rdata")
# save subset of nav data for Di Watters
nav.diwatters <- nav[ ,c("key.z","center_width")]
save(nav.diwatters,file="nav_data_boot.Rdata")

# #############################################################################################
# EXPLORATORY PLOTS OF DATA #####
# #############################################################################################

# #############################################################################################
# plot abundance vs. height for each species #####
# #############################################################################################
species.list <- c("Sebastes","Sebastomus","Sebastolobus")
obs.gg <- droplevels(subset(obs,obs$genus_t %in% species.list & obs$key.z %in% nav.list$key))
height.temp <- (dcast(obs.gg[ ,c("sci_name_full","height_code","counts")],sci_name_full ~ height_code,sum,margins="height_code"))
height.pct <- height.temp[,8]
height.temp[ ,2:7] <- (height.temp[ ,2:7]/height.pct)*100
height.gg <- melt(height.temp[ ,1:7])
names(height.gg) <- c("species","height_code","counts")

H <- ggplot(height.gg, aes(x=height_code,y=counts))
H + geom_bar(stat="identity",binwidth = 1) + coord_flip() + facet_wrap(~species) +
#   H + geom_bar(stat="identity",binwidth = 1) + coord_flip() + facet_wrap(~species, scales = "free") +
	scale_x_discrete("Height above the seabed (m)\n",labels = c("on","0.5","1","2","3",">3")) +
  scale_y_continuous("Percentage",expand=c(0,0)) + theme_bw() +
	theme(title = element_text("Rockfish height above the seabed\n"),
	axis.text.x=element_text(angle=45))
write.csv(height.gg,file="height_gg.csv",quote=FALSE,row.names=FALSE)
ggsave(filename = "rockfish_heights.pdf",width=13.3,height=9.8)
ggsave(filename = "rockfish_heights.png",width=13.3,height=9.8)

# #############################################################################################
# Plot encounter rate by seabed class and depth  #####
# #############################################################################################
#summarize effort by depth and seabed type (sz)
nav.sz <- dcast(melt(nav.sub[,c("geol_prim","depth_stratum","disp")]),depth_stratum+geol_prim~variable,sum)
nav.sz$key.sz <- paste(nav.sz$depth_stratum,nav.sz$geol_prim)
# summarize observations by depth and seabed type (sz)
obs.sz <- dcast(melt(obs.gg[,c("sci_name_full","geol_prim","depth_stratum","counts")]),depth_stratum+geol_prim~sci_name_full,sum)
obs.sz$key.sz <- paste(obs.sz$depth_stratum,obs.sz$geol_prim)
# add survey distance to observation data
obs.sz <- merge(obs.sz,nav.sz[,c("key.sz","disp")],by = "key.sz")
# normalize observations by effort (encounter rate)
obs.sz[,4:39] <- obs.sz[,4:39]/(obs.sz[,"disp"]/1000)
# sort seabed factors for plotting
hab.fac <- as.data.frame(cbind(levels(obs$geol_prim),c(6,5,7,4,1,3,2)))
names(hab.fac) <- c("geol_prim","geol_code")
hab.fac$geol_code <- as.numeric(hab.fac$geol_code)
# summarize species abundance by seabed type
obs.sz.gg <- melt(obs.sz[ ,-40])
# add ordered seabed categories to plot data
obs.sz.gg$geol_code <- as.factor(hab.fac$geol_code[obs.sz.gg$geol_prim])
# RAW ABUNDANCES (PLOT DOMINATED BY A FEW HIGHLY ABUNDANCE SPECIES)
SZ <- ggplot(obs.sz.gg,aes(geol_code,depth_stratum,z=value))
SZ + geom_point(aes(size = value)) + facet_wrap(~variable) +
	scale_x_discrete("Seabed type (low to high complexity)",labels = c("M","S","P","L","C","B","H")) + 
	scale_y_discrete("Depth (m)\n",labels = c("300-400","200-300","100-200","0-100")) + theme_bw() +
	opts(title = "Encounter rate (fish/km) of rockfishes by depth and seabed type\n",
		panel.grid.major = theme_blank(), # switch off major gridlines
		panel.grid.minor = theme_blank()) # switch off minor gridlines
# RESCALED ABUNDANCE VALUES (EACH SPECIES ON A 0-1 SCALE)
obs.sz.rescaled <- ddply(obs.sz.gg,.(variable),function(d){d$value <- rescale(d$value); d})
SZr <- ggplot(obs.sz.rescaled,aes(geol_code,depth_stratum,z=value))
SZr + geom_point(aes(size = value, colour = value),legend=FALSE) + facet_wrap(~variable) +
	scale_colour_gradient(low = "gray70",high = "black")  +
	scale_x_discrete("\nSeabed type (low to high complexity)",labels = c("M","S","P","L","C","B","H")) + 
	scale_y_discrete("Depth (m)\n",labels = c("300-400","200-300","100-200","0-100")) + theme_bw() +
	# annotate("text", size=5,x=1, y=6, label="n") +
	opts(title = "Encounter rate (fish/km) of rockfishes by depth and seabed type\n",
		panel.grid.major = theme_blank(), # switch off major gridlines
		panel.grid.minor = theme_blank()) # switch off minor gridlines
	# scale_y_discrete("Depth (m)") + theme_bw()
ggsave(filename = "rockfish_depth_hab.pdf",width=13.3,height=9.8)
ggsave(filename = "rockfish_depth_hab.png",width=13.3,height=9.8)

# #############################################################################################
# Plot fish lengths by species  #####
# #############################################################################################
length.gg <- dcast(obs.gg[ ,c("sci_name_full","size_code","counts")],sci_name_full ~ size_code,sum)
length.gg <- melt(length.gg)
names(length.gg) <- c("species","size_code","counts")
#plot abundance vs. height for each species
L <- ggplot(length.gg, aes(x=size_code,y=counts))
L + geom_bar(binwidth = 1) + facet_wrap(~species, scales = "free") +
	scale_x_discrete("Total length (cm)",labels = c("<10","10","20","30","40","50","60+")) +
	scale_y_continuous("Frequency") + theme_bw() +
	opts(title = "Rockfish length distributions\n")
ggsave(filename = "rockfish_lengths.pdf",width=13.3,height=9.8)
ggsave(filename = "rockfish_lengths.png",width=13.3,height=9.8)

# #############################################################################################
# Plot species encounter rate by depth  #####
# #############################################################################################
# summarize effort by depth (z)
nav.z <- dcast(melt(nav.sub[,c("depth_stratum","disp")]),depth_stratum~variable,sum)
# summarize observations by depth (z)
obs.z <- dcast(melt(obs.gg[,c("sci_name_full","depth_stratum","counts")]),depth_stratum~sci_name_full,sum)
obs.z <- merge(obs.z,nav.z,by = "depth_stratum")
# normalize observations by effort (encounter rate)
obs.z[,2:37] <- obs.z[,2:37]/(obs.z[,"disp"]/1000)
obs.z.gg <- melt(obs.z[ ,-38])
names(obs.z.gg) <- c("depth","species","counts")
#plot abundance vs. height for each species
Z <- ggplot(obs.z.gg, aes(x=depth,y=counts))
Z + geom_bar(binwidth = 1) + coord_flip() + facet_wrap(~species, scales = "free") +
	scale_x_discrete("Depth (m)\n",labels = c("300-400","200-300","100-200","0-100")) +
	scale_y_continuous("Encounter rate (fish/km)",expand=c(0,0)) + theme_bw() +
	opts(title = "Rockfish encounter rate by depth\n", axis.text.x=theme_text(angle=45))
ggsave(filename = "rockfish_depths.pdf",width=13.3,height=9.8)
ggsave(filename = "rockfish_depths.png",width=13.3,height=9.8)

# #############################################################################################
# Plot species encounter rate by seabed type  #####
# #############################################################################################
# summarize effort by depth (z)
nav.s <- dcast(melt(nav.sub[,c("geol_prim","disp")]),geol_prim~variable,sum)
# summarize observations by depth (z)
obs.s <- dcast(melt(obs.gg[,c("sci_name_full","geol_prim","counts")]),geol_prim~sci_name_full,sum)
obs.s <- merge(obs.s,nav.s,by = "geol_prim")
# normalize observations by effort (encounter rate)
obs.s[,2:37] <- obs.s[,2:37]/(obs.s[,"disp"]/1000)
obs.s.gg <- melt(obs.s[ ,-38])
names(obs.s.gg) <- c("geol_prim","species","counts")
# add size ranks to seabed types
hab.fac.obs <- as.data.frame(cbind(levels(obs.gg$geol_prim),c(6,5,7,4,1,3,2)))
names(hab.fac.obs) <- c("geol_prim","geol_code")
hab.fac.obs$geol_code <- as.numeric(hab.fac.obs$geol_code)
#add ordered seabed categories to plot data
obs.s.gg$seabed_code <- as.factor(hab.fac.obs$geol_code[obs.s.gg$geol_prim])
S <- ggplot(obs.s.gg, aes(x=seabed_code,y=counts))
S + geom_bar(binwidth = 1) + facet_wrap(~species, scales = "free") +
	scale_x_discrete("Seabed type (low to high complexity)",labels = c("M","S","P","L","C","B","H")) + 
	scale_y_continuous("Encounter rate (fish/km)\n") + theme_bw() +
	opts(title = "Rockfish encounter rate by seabed type\n")
ggsave(filename = "rockfish_seabed_type.pdf",width=13.3,height=9.8)
ggsave(filename = "rockfish_seabed_type.png",width=13.3,height=9.8)

# #############################################################################################
# Plot species distribution and abundance by location  #####
# #############################################################################################
obs.gg.rescaled <- ddply(obs.gg,.(sci_name_full),function(d){d$counts <- rescale(d$counts); d})
D <- ggplot(obs.gg.rescaled,aes(x=lon,y=lat))
D + geom_point(aes(size = counts,colour = counts),legend=FALSE) + facet_wrap(~sci_name_full) +
	scale_x_continuous("Longitude") + scale_y_continuous("Latitude\n") + theme_bw() +
	scale_colour_gradient(low = "gray40",high = "black")  +
	opts(title = "Species Abundance at Footprint and Piggy Bank\n",
		axis.text.x=theme_text(angle=45)) 
ggsave(filename = "abundance_dist_by_species.pdf",width=13.3,height=9.8)
ggsave(filename = "abundance_dist_by_species.png",width=13.3,height=9.8)

# #############################################################################################
# Plot seabed observations by location  #####
# #############################################################################################
nav.gg <- nav.sub
hab.fac.nav <- as.data.frame(cbind(levels(nav.gg$geol_prim),c(6,5,7,4,1,3,2)))
names(hab.fac.nav) <- c("geol_prim","geol_code")
hab.fac.nav$geol_code <- as.numeric(hab.fac.nav$geol_code)
#add ordered seabed categories to plot data
nav.gg$seabed_code <- as.factor(hab.fac.nav$geol_code[nav.gg$geol_prim])
H <- ggplot(nav.gg,aes(x=lon,y=lat))
H + geom_point(aes(colour = seabed_code)) +
	scale_x_continuous("Longitude") + scale_y_continuous("Latitude\n") + theme_bw() +
	scale_colour_hue(h=c(0, 180),'Seabed type', labels=c('Mud', 'Sand', 'Pebble','Low-complexity','Cobble','Boulder','High-complexity')) +
	opts(title = "Seabed type at Footprint and Piggy Bank\n")
ggsave(filename = "seabed_observations.pdf",width=13.3,height=9.8)
ggsave(filename = "seabed_observations.png",width=13.3,height=9.8)
write.csv(nav.gg[ ,c("nav_id","object_id","dive_name","lat","lon","slope_dem","geol_prim","geol_sec","geol_ps")],file = "nav_seabed.csv",row.names = FALSE)

# #############################################################################################
# Plot standardized environmental variables by location  #####
# #############################################################################################
vars <- c("lat","lon","depth","salinity","temperature","oxygen_conc","oxygen_sat","slope_dem")
nav.gg <- melt(nav[order(nav$nav_id),vars], id = c("lat","lon"))
nav.gg.rescaled <- ddply(nav.gg,.(variable),function(d){d$value <- rescale(d$value); d})
# # plot actual data
# m <- ggplot(nav.gg,aes(lon,lat))
# m + geom_point(aes(colour = value)) + facet_wrap(~variable, scales = "free")
# plot reclassified data
m1 <- ggplot(nav.gg.rescaled,aes(lon,lat))
m1 + geom_point(aes(colour = value)) + facet_wrap(~variable)
ggsave(filename = "env_data.pdf",width=13.3,height=9.8)
ggsave(filename = "env_data.png",width=13.3,height=9.8)
# #############################################################################################
# Species accumulation (rarefaction) plots for each bank #####
# #############################################################################################
# remove unidentified species from species accumulation curve data
genus.list <- c("Sebastes","Sebastomus","Sebastolobus","Chromis","Caulolatilus","Oxyjulis","Trachurus","Merluccius","Ophiodon")
# species list for accum curves
spp.list.a <- unique(obs$sci_name_full[obs$genus_t %in% genus.list])
spp.list.a <- droplevels(spp.list.a[-grep("sp.",spp.list.a)])
obs.accum.fp <- dcast(obs[obs$sci_name_full %in% spp.list.a & obs$site_name == "The Footprint",c("dive_name","species_code","counts")],dive_name ~ species_code,sum)
obs.accum.fp <- obs.accum.fp[ ,-1]
obs.accum.pb <- dcast(obs[obs$sci_name_full %in% spp.list.a & obs$site_name == "Piggy Bank",c("dive_name","species_code","counts")],dive_name ~ species_code,sum)
obs.accum.pb <- obs.accum.pb[ ,-1]
# calculate curves
accum.exact.fp  <- specaccum(obs.accum.fp, method = "exact", permutations = 100, conditioned =TRUE, gamma = "jack1")
accum.random.fp <- specaccum(obs.accum.fp, method = "random", permutations = 100, conditioned =TRUE, gamma = "jack1")
accum.coll.fp   <- specaccum(obs.accum.fp, method = "collector", permutations = 100, conditioned =TRUE, gamma = "jack1")
accum.exact.pb  <- specaccum(obs.accum.pb, method = "exact", permutations = 100, conditioned =TRUE, gamma = "jack1")
accum.random.pb <- specaccum(obs.accum.pb, method = "random", permutations = 100, conditioned =TRUE, gamma = "jack1")
accum.coll.pb   <- specaccum(obs.accum.pb, method = "collector", permutations = 100, conditioned =TRUE, gamma = "jack1")
# plot results of exact method
layout(1)
plot(accum.exact.fp,xlab="ROV samples",ylab="Number of species")
plot(accum.exact.pb,xlab="ROV samples",ylab="Number of species")	

# plot results of random method
boxplot(accum.random.fp)	
boxplot(accum.random.pb)	
# convert accum list to a data frame for plotting in ggplot2
accum.exact.fp.gg <- data.frame(accum.exact.fp[[3]],accum.exact.fp[[4]],accum.exact.fp[[5]])
names(accum.exact.fp.gg) <- c("samples","richness","stdev")
accum.exact.fp.gg$site_name <- "The Footprint"
accum.exact.pb.gg <- data.frame(accum.exact.pb[[3]],accum.exact.pb[[4]],accum.exact.pb[[5]])
names(accum.exact.pb.gg) <- c("samples","richness","stdev")
accum.exact.pb.gg$site_name <- "Piggy Bank"
accum.exact.all.gg <- rbind(accum.exact.fp.gg,accum.exact.pb.gg)

# plot exact results as points
ae <- ggplot(accum.exact.all.gg,aes(group = site_name,y = richness,x = samples))
limits <- aes(ymax = richness+stdev,ymin = richness-stdev)
ae + geom_errorbar(limits, colour = "black",width=0.3) + geom_point(aes(colour = site_name),size=3) + 
	scale_x_continuous('\nNumber of transects') + scale_y_continuous('Species richness\n') + theme_bw() + labs("Site name") +
	opts(title="Species Accumulation Curves for \nthe Footprint and Piggy Banks\n",
	legend.position = c(0.85,0.5),
	legend.title = theme_blank(),
	panel.grid.major = theme_blank(), # switch off major gridlines
	panel.grid.minor = theme_blank())		
ggsave(filename = "rockfish_accum_exact.pdf", width=13.3, height=9.8)
ggsave(filename = "rockfish_accum_exact.png", width=13.3, height=9.8)

# plot random results as box plots
# format data for The Footprint
accum.rand.fp.gg <- data.frame(accum.random.fp[[6]])
accum.rand.fp.gg$samples <- seq(1,nrow(accum.rand.fp.gg),1)
accum.rand.fp.gg$site_name <- as.factor("The Footprint")
accum.rand.fp.gg <- melt(accum.rand.fp.gg,id=c("site_name","samples"))
# format data for Piggy Bank
accum.rand.pb.gg <- data.frame(accum.random.pb[[6]])
accum.rand.pb.gg$samples <- seq(1,nrow(accum.rand.pb.gg),1)
accum.rand.pb.gg$site_name <- as.factor("Piggy Bank")
accum.rand.pb.gg <- melt(accum.rand.pb.gg,id=c("site_name","samples"))
# combine data from both banks
accum.rand.all.gg <- rbind(accum.rand.fp.gg,accum.rand.pb.gg)
# plot Footprint data
arand.fp <- ggplot(accum.rand.fp.gg, aes(factor(samples), value))
arand.fp + geom_jitter(alpha = 0.25) + geom_boxplot(fill = "grey80",colour = "grey40",outlier.colour = "red", outlier.size = 2) + 
	xlab("\nSamples") + ylab("Species richness\n") + opts(title="Species Accumulation at \nFootprint and Piggy Banks\n") + theme_bw() +
	scale_x_discrete("samples") + scale_y_continuous(limits = c(0, 35))
# plot Piggy Bank data
arand.pb <- ggplot(accum.rand.pb.gg, aes(factor(samples), value))
arand.pb + geom_jitter(alpha = 0.25) + geom_boxplot(fill = "grey80",colour = "grey40",outlier.colour = "red", outlier.size = 2) + 
	xlab("\nSamples") + ylab("Species richness\n") + opts(title="Species Accumulation at \nFootprint and Piggy Banks\n") + theme_bw() +
	scale_x_discrete("samples") + scale_y_continuous(limits = c(0, 35))
# plot data from both banks
arand.all <- ggplot(accum.rand.all.gg, aes(x=as.factor(samples), y=value, fill=site_name))
arand.all + geom_jitter(alpha = 0.25) + geom_boxplot(outlier.colour = "red", outlier.size = 2) +  
	xlab("\nSamples") + ylab("Species richness\n") + theme_bw() +
	opts(title="Species Accumulation at \nFootprint and Piggy Banks\n",
		legend.position = c(0.85,0.5),
		legend.title = theme_blank(),
		panel.grid.major = theme_blank(), # switch off major gridlines
		panel.grid.minor = theme_blank()) 
# save box plots
ggsave(filename = "rockfish_accum_random.pdf", width=13.3, height=9.8)
ggsave(filename = "rockfish_accum_random.png", width=13.3, height=9.8)

# #############################################################################################
# Bootstrap analysis of species abundance #####
# #############################################################################################
# add biomass estimates to obs data
length.bin 	<- sort(unique(obs$size_code[is.na(obs$size_code)==FALSE]))
length.cm 	<- seq(5,65,10) 
l.lookup <- data.frame(length.bin,length.cm)
obs2 <- droplevels(merge(obs,l.lookup,by.x="size_code",by.y="length.bin",all.x = TRUE))
obs2 <- merge(obs2,lw.data[,c("sci_name_full","a","b")],by="sci_name_full",all.x=TRUE)
obs2$biomass <- obs2$a*(obs2$length.cm^obs2$b)*obs2$counts

# summarize length data by species for binned length estimates
length.summ.spp <- as.data.frame(t(sapply(sort(unique(obs2$sci_name_full)),function(x){
	data <- subset(obs2,sci_name_full == x)
	species <- as.character(data$sci_name_full[1])
	n.meas	<- as.numeric(dim(data)[1]) #good
	l.mean 	<- as.numeric(mean(data$length.cm))
	l.sd 	<- as.numeric(sd(data$length.cm))
	c(spp = species, n.meas=n.meas,l.mean = l.mean,l.sd=l.sd)
})))
length.summ.spp
write.csv(length.summ.spp, file = "rockfish_length_summary.csv",row.names = FALSE, na = "")

# list of species to be analyzed
genus.list <- c("Sebastes","Sebastomus","Sebastolobus","Chromis","Caulolatilus","Oxyjulis","Trachurus","Merluccius","Ophiodon")
# summarize fish abundance by transect
obs.counts <- droplevels(melt(dcast(obs2[obs2$key.z %in% nav.list$key.z & obs2$genus_t %in% genus.list,c("key.z","sci_name_full","counts")],key.z~sci_name_full,sum)))
names(obs.counts) <- c("transect","species","counts")
obs.counts$key.ts <- paste(obs.counts$transect,obs.counts$species)
# summarize fish biomass by transect
obs.biomass <- droplevels(melt(dcast(obs2[obs2$key.z %in% nav.list$key.z & obs2$genus_t %in% genus.list,c("key.z","sci_name_full","biomass")],key.z~sci_name_full,sum)))
names(obs.biomass) <- c("transect","species","biomass")
obs.biomass$key.ts <- paste(obs.biomass$transect,obs.biomass$species)
# merge obs data
obs.boot <- merge(obs.counts,obs.biomass[,c("key.ts","biomass")],by = "key.ts")
# extract dive name from transect name
obs.boot$dive_name <- as.factor(substr(obs.boot$transect,1,7))
# add site names to obs.boot
obs.boot <- merge(obs.boot,site.names,by = "dive_name")

# summarize transect length by transect
nav.boot <-  droplevels(dcast(nav.sub[nav.sub$key.z %in% nav.list$key.z ,c("key.z","disp")],key.z~.,sum))
names(nav.boot) <- c("transect","length")
# merge obs and nav data
boot.data <- merge(obs.boot,nav.boot,by = "transect")
boot.data <- droplevels(merge(boot.data,nav.list[ ,c("key.z","depth_stratum")],by.x="transect",by.y="key.z"))
boot.data$key.ts <- NULL
boot.data$key <- paste(boot.data$dive_name,boot.data$depth_stratum)
spp.list <- unique(boot.data$species)
depth.list <- droplevels(unique(boot.data$depth_stratum))
# save data for bootstrap analysis
write.csv(boot.data, file = "boot_data.csv",row.names = FALSE, na = "")
save(boot.data,file="boot_data.Rdata")

# initialize variables for bootstrapped results
# stratum data
spp.boot 	<- character()
site.boot 	<- character()
depth.boot 	<- character()
n.transects	<- numeric()
length.t	<- numeric()
e.rate		<- numeric()
n.obs		<- numeric()
# abundance results
D.boot		<- numeric()
N.boot.mean <- numeric()
N.boot.cv 	<- numeric()
N.boot.ciL 	<- numeric()
N.boot.ciH 	<- numeric()
N.anal.mean <- numeric()
N.anal.cv	<- numeric()
N.anal.ciL	<- numeric()
N.anal.ciH	<- numeric()
# biomass results
B.boot.mean <- numeric()
B.boot.cv 	<- numeric()
B.boot.ciL 	<- numeric()
B.boot.ciH 	<- numeric()

# TEMP variables for testing the loops
# setwd(wkdir)
# i = "Sebastes levis"
# j = "The Footprint"
# k = "(-100,0]"
# k = "(-200,-100]"
# k = "(-300,-200]"
# k = "(-400,-300]"
# compute abundance and variance estimates for each species
pb1 <- txtProgressBar(min = 0, max = length(spp.list), style = 3)
for (i in spp.list){
	# subset data by species
	spp.temp <- droplevels(subset(boot.data,boot.data$species == i))
	for (j in unique(spp.temp$site_name)){
		# subset data by site (Footprint or Piggy Bank)
		site.temp <- droplevels(subset(spp.temp,spp.temp$site_name == j))
		for (k in unique(site.temp$depth_stratum)){
			# subset data by depth stratum (100m bins from 0-400m)
			depth.temp <- droplevels(subset(site.temp,site.temp$depth_stratum == k))
			# subset nav to include data from dives in this stratum
			nav.boot <- subset(nav,nav$key.z %in% depth.temp$key)
			# Total survey area (km sq, as defined by Yok and Watters, excluding the NE sliver)
			A <- sum(survey.area$area[survey.area$site_name %in% depth.temp$site_name[1] & survey.area$depth_stratum %in% depth.temp$depth_stratum[1]])
			n <- sum(depth.temp$counts)		 # total n of sightings
			L <- sum(depth.temp$length)/1000 # total transect distance (km)
			T <- length(depth.temp$transect) # number of transects
			g.mean 	<- 1					 # detection probability mean
			g.cv 	<- 0					 # detection probablility CV
			# w.mean <- 1.5/1000			 # strip width mean (km); assumed to be 3m prior to 3Beam analysis
			# w.cv 	 <- 0.5					 # strip width CV; assumed to be 50% prior to 3Beam analysis
			w.mean 	<- mean(nav.boot$center_width)/(1000) 	 # strip width mean (km) from 3Beam analysis
			w.cv 	  <- sd(nav.boot$center_width)/(1000)/w.mean # strip width CV from 3Beam analysis
			meanN 	<- (A*n)/(L*w.mean)	 	  # mean total abundance assuming strip transect methodology
			meanD 	<- meanN/A				      # mean density
								
			# calculate mean (theta-hat) and SEM (SE-hat) for all transects in this stratum
			# calculate transect mean abundance (N.tr)
			n.tr <- depth.temp$counts		    # species abundance for each transect in this stratum
			b.tr <- depth.temp$biomass		  # species biomass for each transect in this stratum
			L.tr <- depth.temp$length/1000	# length of each transect in this stratum
			w.tr <- numeric()
			for (kk in unique(nav.boot$key.z)){
				nav.boot.sub <- subset(nav.boot,nav.boot$key.z %in% kk)
				w.tr <- c(w.tr,mean(nav.boot.sub$center_width)/(1000)) # mean strip width of each transect
			}
			
			# initialize transect abundance (N.tr) and biomass (B.tr)
			N.tr <- rep(0,length(n.tr))
			B.tr <- rep(0,length(n.tr))			
			# calculate N for each transect (N.tr)
			for (m in 1:length(N.tr)){ 
				N.tr[m] <- (A*n.tr[m])/(L.tr[m]*w.tr[m])		  # abundance per transect
				B.tr[m] <- (A*b.tr[m])/(L.tr[m]*w.tr[m])/1000	# biomass per transect (kg)
				}
			# calculate the CV of the encounter rate of each species using analytical approx.
			nL	 <- n.tr/L.tr									# Encounter rate (fish per km) for each transect
			D	 <- n.tr/w.tr*L.tr
			cvnL <- sd(nL)/(mean(nL)*sqrt(dim(depth.temp)[1]))	# CV of encounter rate in all transects
			cvD  <- sd(D)/(mean(D)*sqrt(dim(depth.temp)[1]))	# CV of density in all transects
			cvN  <- sqrt((cvnL)^2 + (w.cv)^2 + (g.cv)^2) 		# CV of N (total abundance) in all transects
			cvN.d<- sqrt((cvD)^2 + (g.cv)^2) 		# CV of N (total abundance) in all transects
			
			# calculate CIs for normally distributed data
			Z <- 1.645
			sdN.tr <- cvN*mean(N.tr)
			seN.tr <- sdN.tr/sqrt(length(N.tr))
			upperCIN.norm <- mean(N.tr) + (Z*seN.tr)
			lowerCIN.norm <- mean(N.tr) - (Z*seN.tr)
			
			# calculate the upper and lower abundance CIs for Log-normal data
			c <- exp(Z*sqrt(log(1+(sdN.tr^2/mean(N.tr)^2))))
			lowerCIN.logn <- mean(N.tr)/c
			upperCIN.logn <- mean(N.tr)*c
						
			# Calculate N and CV(N) using a non-parametric bootstrap
			B <- 1000
			boot.samples <- length(N.tr)
			N.boot <- numeric()
			while (length(N.boot) < B){
				bN.boot <- mean(sample(N.tr, boot.samples, replace = TRUE))
				N.boot <- c(N.boot,bN.boot)
			}
			# Calculate B and CV(B) using a non-parametric bootstrap
			B.boot <- numeric()
			# if biomass info not available, return NA
			if(is.na(mean(B.tr))) {			
				B.boot  <- NA
			# else, compute bootstrapped mean biomass
			} else {
				while (length(B.boot) < B){
					bB.boot <- mean(sample(B.tr, boot.samples, replace = TRUE))
					B.boot <- c(B.boot,bB.boot)
					}
			}
			# abundance estimates
			bN <- mean(N.boot)								# mean of bootstrapped abundance estimates
			bSD.n <- sqrt(sum((N.boot-bN)^2)/(B-1))			# SD of bootstrapped abundance estimates (even though it says SE)
			bCV.n <- bSD.n/bN 								# CV of bootstrapped abundance estimates
			bCI.n <- quantile(N.boot, probs=c(0.05, 0.95),na.rm = TRUE) 	# quantile 90% conficence interval
			
			# biomass estimates
			bB <- mean(B.boot, rm.na = TRUE)				# mean of bootstrapped biomass estimates
			bSD.b <- sqrt(sum((B.boot-bB)^2)/(B-1))			# SE of bootstrapped biomass estimates
			bCV.b <- bSD.b/bB 								# CV of bootstrapped biomass estimates
			bCI.b <- quantile(B.boot, probs=c(0.05, 0.95),na.rm = TRUE) 	# quantile 90% conficence interval
			
			# histogram of bootstrapped mean samples	
			hist(N.boot, main = "Distribution of bootstrapped means",xlab="Mean abundance")
			abline(v = mean(N.boot),col = "gray70", lty=4)
			
			spp.boot 	<- c(spp.boot,i) 	 			# add species name to list
			site.boot 	<- c(site.boot,j) 	 			# add species name to list
			depth.boot 	<- c(depth.boot,k)	 			# add depth strata to list
			n.transects <- c(n.transects,T)  			# add number of transects
			length.t	<- c(length.t,L)	 			# add total survey distance
			e.rate		<- c(e.rate,mean(nL))			# add encounter rate (no. per km)
			n.obs		<- c(n.obs,n)		 			# add number of fish
			# density extimate
			D.boot		<- c(D.boot,bN/A)		 		# add density
			# abundance estimates
			N.boot.mean <- c(N.boot.mean,bN) 	 		# add bootstrapped mean
			N.boot.cv 	<- c(N.boot.cv,bCV.n)		 		# add bootstrapped CV
			N.boot.ciL 	<- c(N.boot.ciL,bCI.n[1])  		# add bootstrapped CI-lower
			N.boot.ciH 	<- c(N.boot.ciH,bCI.n[2]) 		# add bootstrapped CI-upper
			N.anal.mean <- c(N.anal.mean,meanN)  	  	# add analytical mean
			N.anal.cv	<- c(N.anal.cv,cvN)		 	  	# add analytical CV
			N.anal.ciL	<- c(N.anal.ciL,lowerCIN.norm)	# add analysital CI-lower
			N.anal.ciH	<- c(N.anal.ciH,upperCIN.norm)	# add analysital CI-upper
			# biomass estimates
			B.boot.mean <- c(B.boot.mean,bB) 	 		# add bootstrapped mean
			B.boot.cv 	<- c(B.boot.cv,bCV.b)	 		# add bootstrapped CV
			B.boot.ciL 	<- c(B.boot.ciL,bCI.b[1])		# add bootstrapped CI-lower
			B.boot.ciH 	<- c(B.boot.ciH,bCI.b[2])		# add bootstrapped CI-upper
		}
	}
	setTxtProgressBar(pb1, i)
}
close(pb1)
# combine results from all species
boot.summ <- data.frame(spp.boot,site.boot,depth.boot,n.transects,length.t,n.obs,e.rate,
		D.boot,N.boot.mean,N.boot.cv,N.boot.ciL,N.boot.ciH,N.anal.mean,N.anal.cv,N.anal.ciL,N.anal.ciH,
		B.boot.mean,B.boot.cv,B.boot.ciL,B.boot.ciH)
# rename columns
names(boot.summ) <- c("spp","site","depth","n.transects","length.t","n.obs","e.rate",
	"D.boot","N.boot.mean","N.boot.cv","N.boot.ciL","N.boot.ciH","N.anal.mean","N.anal.cv","N.anal.ciL","N.anal.ciH",
	"B.boot.mean","B.boot.cv","B.boot.ciL","B.boot.ciH")
# replace negative CIs with zeros		
boot.summ <- do.call(data.frame, lapply(boot.summ, function(x) {
  if(is.numeric(x)) {x[x < 0] <- 0} 
  x
}))
boot.summ <- do.call(data.frame, lapply(boot.summ, function(x) {
  if(is.numeric(x)) {x[is.nan(x)] <- NA} 
  x
}))
# sort results by species, site, and depth
write.csv(boot.summ, file = "boot_summary_all.csv",row.names = FALSE, na = "")

# #############################################################################################
# Summarize bootstrap results by bank and species #####
# #############################################################################################
# summarize species abundance and CV for all banks
boot.summ$key <- as.character(paste(boot.summ$spp,boot.summ$site))
boot.summ.site <- as.data.frame(t(sapply(sort(unique(boot.summ$key)),function(x){
	data 	<- subset(boot.summ,key == x)
	species <- as.character(data$spp[1])
	site	<- as.character(data$site[1])
	N.mean 	<- as.numeric(sum(data$N.boot.mean, na.rm = TRUE)) #good
	N.cv 	<- as.numeric(sqrt(sum(data$N.boot.cv^2 * data$N.boot.mean^2, na.rm = TRUE))/sum(data$N.boot.mean), na.rm = TRUE)
	B.mean 	<- as.numeric(sum(data$B.boot.mean, na.rm = TRUE)) #good
	B.cv 	<- as.numeric(sqrt(sum(data$B.boot.cv^2 * data$B.boot.mean^2, na.rm = TRUE))/sum(data$B.boot.mean), na.rm = TRUE)
	c(spp = species, site = site,N.mean = N.mean, N.cv = N.cv,B.mean = B.mean,B.cv = B.cv)
})))
write.csv(boot.summ.site, file = "boot_summary_bank.csv",row.names = FALSE, na = "")

# summarize species abundance and CV for all banks
boot.summ.spp <- as.data.frame(t(sapply(sort(unique(boot.summ$spp)),function(x){
	data <- subset(boot.summ,spp == x)
	species <- as.character(data$spp[1])
	N.mean 	<- as.numeric(sum(data$N.boot.mean, na.rm = TRUE)) #good
	N.cv 	<- as.numeric(sqrt(sum(data$N.boot.cv^2 * data$N.boot.mean^2, na.rm = TRUE))/sum(data$N.boot.mean), na.rm = TRUE)
	B.mean 	<- as.numeric(sum(data$B.boot.mean, na.rm = TRUE)) #good
	B.cv 	<- as.numeric(sqrt(sum(data$B.boot.cv^2 * data$B.boot.mean^2, na.rm = TRUE))/sum(data$B.boot.mean), na.rm = TRUE)
	c(spp = species, N.mean = N.mean, N.cv = N.cv,B.mean = B.mean,B.cv = B.cv)
})))
boot.summ.spp
write.csv(boot.summ.spp, file = "boot_summary_spp.csv",row.names = FALSE, na = "")

# #############################################################################################
# Plots of fish lengths of target species measured using the reference lasers #####
# #############################################################################################
lengths.l <- read.csv("fish_lengths_lasers.csv")
LL.gg <- ggplot(lengths.l,aes(x=total_length))
LL.gg + geom_histogram() + facet_wrap(~sci_name_full, scales="free_y") + theme_bw() +
	scale_x_continuous("\nTotal length (cm)\n") + scale_y_continuous("Frequency\n") + 
	opts(title = "Rockfish lenghts using lasers\n",legend.position="none")
# LL.gg + geom_histogram() + facet_wrap(~sci_name_full, scales="free")  + theme_bw()
# LL.gg + geom_histogram(aes(y = ..density.., alpha=0.5)) + geom_density() + facet_wrap(~sci_name_full, scales="free_y") +
	# scale_x_continuous("\nTotal length (cm)\n") + scale_y_continuous("Density\n") + theme_bw() + 
	# opts(title = "Rockfish lenghts using lasers\n",legend.position="none")
ggsave(filename = "rockfish_lengths_lasers.pdf",width=13.3,height=9.8)
ggsave(filename = "rockfish_lengths_lasers.png",width=13.3,height=9.8)

# summarize length estimates by species for laser measurements
length.summ.l.spp <- as.data.frame(t(sapply(sort(unique(lengths.l$sci_name_full)),function(x){
	data <- subset(lengths.l,sci_name_full == x)
	species <- as.character(data$sci_name_full[1])
	n.meas	<- as.numeric(dim(data)[1]) #good
	l.mean 	<- as.numeric(mean(data$total_length))
	l.sd 	<- as.numeric(sd(data$total_length))
	c(spp = species, n.meas=n.meas,l.mean = l.mean,l.sd=l.sd)
})))
write.csv(length.summ.l.spp, file = "rockfish_length_summary_lasers.csv",row.names = FALSE, na = "")

# #############################################################################################
# Analysis of optimal sample allocation #####
# #############################################################################################
# Neyman sampling formula
# h = stratum
# n = sample size
# N = population size
# S = standard deviation
# n.h <- n*(N.h * S.h)/(sum(N.h[i]*S.h[i]))

# Determine optimal sample allocation across both banks
neyman.all <- data.frame()
for (ii in unique(boot.summ$spp)){
neyman.temp <- subset(boot.summ,boot.summ$spp == ii)
spp.na <- neyman.temp$spp
site.na <- neyman.temp$site
depth.na <- neyman.temp$depth
t.act <- neyman.temp$n.transects
t.opt <- round((sum(neyman.temp$n.transects)*(neyman.temp$N.boot.mean^2 * neyman.temp$N.boot.cv))/sum(neyman.temp$N.boot.mean^2 * neyman.temp$N.boot.cv, na.rm=TRUE))
n.h <- data.frame(spp.na,site.na,depth.na,t.act,t.opt)
neyman.all <- rbind(neyman.all,n.h)
}
# create key to merge with neyman.site
neyman.all$key <- paste(neyman.all$site.na,neyman.all$depth.na,neyman.all$spp.na)
# Determine optimal sample allocation within each bank
neyman.site <- data.frame()
for (ii in unique(boot.summ$spp)){
	neyman.temp <- subset(boot.summ,boot.summ$spp == ii)
	for (jj in unique(neyman.temp$site)){
		n.h <- data.frame()
		neyman.site.temp <- subset(neyman.temp,site == jj)
		spp.ns <- neyman.site.temp$spp
		site.ns <- neyman.site.temp$site
		depth.ns <- neyman.site.temp$depth
		t.act <- neyman.site.temp$n.transects
		t.opt <- round((sum(neyman.site.temp$n.transects)*(neyman.site.temp$N.boot.mean^2 * 
			neyman.site.temp$N.boot.cv))/sum(neyman.site.temp$N.boot.mean^2 * neyman.site.temp$N.boot.cv, na.rm=TRUE))
		n.h <- data.frame(spp.ns,site.ns,depth.ns,t.act,t.opt)
		neyman.site <- rbind(neyman.site,n.h)
		}		
}
# create a key to merge with neyman.all
neyman.site$key <- paste(neyman.site$site.ns,neyman.site$depth.ns,neyman.site$spp.ns)
# merge with neyman.all
neyman.summ <- merge(neyman.all,neyman.site[,c("key","t.opt")],by="key")
neyman.summ <- neyman.summ[ ,-1]
write.csv(neyman.summ,file="neyman_summary.csv", row.names = FALSE,na="")

# #############################################################################################
# Analysis of biodiversity #####
# #############################################################################################
# compute biodiversity statistics for each transect
obs.div <- dcast(obs.boot[ ,c("transect","species","counts")],transect~species,sum)
div.mat <- data.matrix(obs.div)
obs.div$H.shannon <- diversity(div.mat,index="shannon")
obs.div$H.simpson <- diversity(div.mat,index="simpson")
obs.div$S.rich <- specnumber(div.mat)
obs.div$S.raref <- rarefy(div.mat, sample = 22)
obs.div <- merge(obs.div,nav.summ.trans[,c("transect","depth_stratum","site_name")],by="transect")
obs.div.transect <- obs.div[ ,c("transect","site_name","depth_stratum","H.shannon","H.simpson","S.rich","S.raref")]
write.csv(obs.div.transect,file="diversity_transect.csv",row.names=FALSE)
# summarize biodiversity statistics for each site and depth
obs.div.summ <- dcast(melt(obs.div.transect),site_name+depth_stratum~variable,mean)
write.csv(obs.div.summ[ ,c("site_name","depth_stratum","S.rich","S.raref","H.shannon","H.simpson")],file="diversity_stratum.csv",row.names=FALSE)

# #############################################################################################
# Analysis of fish reactions #####
# #############################################################################################
# summarize rockfish reactions by species
rxn.codes <- read.csv("rxn_codes.csv")
# subset obs to get only reaction data
obs.rxn <- obs[obs$genus_t %in% genus.list,c("sci_name_full","reaction","counts")]
# add description of reactions
obs.rxn <- merge(obs.rxn,rxn.codes,by = "reaction")
# summarize counts for each species and reaction type
rxn <- dcast(melt(obs.rxn[ ,c("sci_name_full","descr","counts")], value.name="counts"),sci_name_full~descr,margins = TRUE,sum)
# convert counts to proportions
rxn[,2:8] <- rxn[,2:8]/rxn[,9]
# add common names to the table and write to CSV file
rxn <- merge(rxn,spp[,c("sci_name_full","common_name")], by = "sci_name_full")
write.csv(rxn,file="reaction_summary.csv",row.names=FALSE)

#######################################################################################################
# OLD OR UNUSED CODE #####
#######################################################################################################
#for exploring long time lags
# fov.sub <- subset(fov.nav.sync,fov.nav.sync$lag_s > 200)
# temp.nav <- subset(nav[ ,c("nav_id","date_time")],nav$dive_name == "11-338A")
# temp.fov <- subset(fov.sync[ ,c("id","date_time","center_width")],fov.output$dive_name == "11-338A")
# range(temp.fov$date_time)
# range(temp.nav$date_time)
