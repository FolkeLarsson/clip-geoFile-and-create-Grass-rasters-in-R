## runs from Grass command prompt, called by Rscript
## calculate mean values from different circular areas and testing them against values from random points 
## two different means are calculated, one from original map-file and one setting a minimum value removing low values, 
##  in this case open areas
## for testing, non-parametric methods are used, ks.test and wilcox.test 
## tested with rt90 and GTiff 
## tested on Ubuntu 14 and windows 10 
## arguments from command line 
## uses existing location and mapset
## paths to executables hardcoded in script for windows and else
## Folke Larsson Boden, Sweden 2017 


library(XML)
library(sp)
library(rgdal)
library(maptools)
library(spgrass6)
library(raster)
library(MASS)
library(stringr)

# init
args<-commandArgs(TRUE)
options(warn=0)


# constants, table columns and shortened version of rt90 proj4 
rt90 =  "+proj=tmerc +lat_0=0 +lon_0=15.80827777777778 +k=1 +x_0=1500000 +y_0=0 +ellps=bessel +units=m +no_defs"
tab_cols <- c("mean r-250", "mean r-500", "mean r-750", "mean r-1000", "mean r-250 nocc", "mean r-500 nocc", "mean r-750 nocc", "mean r-1000 nocc")
matr_cols = c("var1", "nVar1", "var2", "nVar2", "radius", "D/W", "p", "method", "alternative")


# arguments from command line
in_MinLimit     <- args[1]  # for adjusting min-value in maps from 0 to a defined one, in this case excluding open areas, clearcuts etc
in_Prefix       <- args[2]  # directory name for both logging dir and csv-files etc
in_GrassRaster  <- args[3]  # name of raster map from Grass 
in_ShapeFile    <- args[4]  # path and name of shape file 
in_Location     <- args[5]  # existing Grass location
in_Mapset       <- args[6]  # existing Grass mapset
in_ShpAttribute <- args[7]  # name of attribute in shape layer, not mandatory if name is in shape file 


#variables
shp_mean_vect =  c()  # contain mean values from shape file for different rediuses
p100_mean_vect = c()  # contain mean values from random points for different rediuses
filename_divider <- "_-_"


# directory for log- and error-files and some initial system parameters
work_dir <- paste("dir_", in_Prefix, sep="")
file_sep <- .Platform$file.sep
file_psep   <- .Platform$path.sep
system_name <- Sys.info()['sysname']
system_info <- str(Sys.info)
OS_type     <- .Platform$OS.type
os_version  <- version$os

# for command prompt 
print("-------- command line arguments ------------- ")
print(args)
print("-------- end command line arguments --------- ")

# setting file separator if OS is Windows
# then init Grass paths depending on OS, only tested on windows 10 and Ubuntu 14 
if (OS_type == "windows") {
	file_sep <- "\\\\" 
	curr_gisBase  <- "C:\\OSGeo4W64\\apps\\grass\\grass-6.4.3"
	curr_gisDbase <- "C:\\GIS\\grass\\grassdata"	
} else {
	curr_gisBase  <- "/usr/lib/grass64/"
	curr_gisDbase <- "/home/folke/grass/grassdata/"	
}


# create dir for logging 
could_create_dir <- dir.create(work_dir)
print(paste("  logging directory path: ", could_create_dir))
dir_path = paste(work_dir, file_sep, sep="")
tmp_path <- paste(dir_path, paste(file_sep, "tmp", sep=""), sep="")
print(paste(" dir: ", dir_path))

# directing output to files 
msg <- file(paste(dir_path, "errors_and_warnings.log", sep=""), open="wt")
out <- file(paste(dir_path, "logging_output.log" , sep=""), open="wt")
sink(msg, type="message")
sink(out, type="output")


## functions 
# replacing values lower then given nim value with NA
set_min_value_to_zero <- function(x) { 
 x[x<as.numeric(in_MinLimit)] <- 0 ; return(x)
}

set_zero_to_NA <- function(x) { 
 x[x<1] <- NA; return(x)
}


# returning rownames from attribute in shapefile, attribute "name" choosen or second attribute selected if no value given
# ToDo: checck for valid shapefile object
get_shape_rownames = function(v_shp, v_AttributeName) {
	row_names <- "uninitialized"
	shp_attr_vect <- unlist(names(v_shp))
	l_shp_attr_vect <- length(shp_attr_vect)
	if(is.na(v_AttributeName)) { # no given attribute in layer
		row_names_attr <- match("name", shp_attr_vect) # check for standard attribute "name"
		if(!is.na(row_names_attr))
			row_names <- paste("v_shp@data$", "name", sep="")
		else { # returning first attribute
			shp_attr_vect <- unlist(names(v_shp))
			l_shp_attr_vect <- length(shp_attr_vect)
			if((l_shp_attr_vect > 0) && (!is.na(shp_attr_vect)) ) {
				row_names <- paste("v_shp@data$", shp_attr_vect[1], sep="") # choose first attribute
			} else {
				print(paste("v_shp@data$", shp_attr_vect[1], sep=""))
				stop(" no attribute in shapefile ")
			}
		}
	} else { # name of attribute given 
		row_names_attr <- match(v_AttributeName, shp_attr_vect)
		if(!is.na(row_names_attr))
			row_names <- paste("v_shp@data$", v_AttributeName, sep="")
		else {
			print(paste(" in func: ", "v_shp@data$", v_AttributeName, sep=""))
			stop(" given name doesnt match any attribute in shapefile, or no valid shapefile")
		}
	}
	return(eval(parse(text=row_names)))
} # get_shape_rownames


# wrapper for statistical calculations 
stat_test = function(var1, var2, in_alternative, radius, method) {
 stat_vect = c()  # return vector
 stat_vect = c()
 #test_vect = c()

 var1_name <- var1$variable[1]
 var2_name <- var2$variable[1]
 var1_length <- length(var1[,1])
 var2_length <- length(var2[,1])
 var1_noNA <- na.omit(var1)
 var2_noNA <- na.omit(var2)
 #test_vect = var1_noNA[,1]

 var1_noNA_length <- length(var1_noNA[,1])
 var2_noNA_length <- length(var2_noNA[,1])
 v_alternative = paste("alternative = \"", in_alternative, "\"", sep="")
 v_alternative2 = paste("alternative=", '\"', in_alternative, '\"', sep="")
 v_alternative3 = as.character(paste("alternative=", "\"", sep=""))
 v_alternative4 = paste(' alternative = "', in_alternative, '"', sep="")

 if(method == "ks") {
  v_ks <- ks.test(var1_noNA[,1], var2_noNA[,1], alternative=in_alternative)
 } else if(method == "wc") {
  v_ks <- wilcox.test(var1_noNA[,1], var2_noNA[,1], alternative=in_alternative)
 } else {
   print(in_alternative)
   print(method)
   print(" NA NA NA no such method! ")
   stop()
   #v_ks <- NA  # should be caught
 }

 v_dataname <- v_ks$data.name
 v_dataname_l <- str_length(v_dataname)
 v_cut_pos <- str_locate(v_dataname, " and ")
 v_dataname_1 <- substr(v_dataname, 1, v_cut_pos-1)
 v_dataname_2 <- substr(v_dataname, v_cut_pos+5, v_dataname_l)

 stat_vect[1] = var1_name
 stat_vect[2] = var1_noNA_length
 stat_vect[3] = var2_name
 stat_vect[4] = var2_noNA_length
 stat_vect[5] = radius
 stat_vect[6] = round(v_ks$statistic, 2)
 stat_vect[7] = round(v_ks$p.value, 2)
 stat_vect[8] = v_ks$method
 stat_vect[9] = v_ks$alternative

 return(stat_vect)
} #stat_test



# wrapper for boxplot function
f_boxplot = function(v_shp_var, v_p100_var, v_radius, v_cc_type, v_mean_crop, v_mean_nocrop, v_prefix, v_dir_path) {
	mean_p100 <- mean(v_p100_var, na.rm = TRUE)
	p100_var_noNA <- na.omit(v_p100_var)
	if(v_cc_type == "nocc") {
		ylab_txt <- "mean values from raster without open areas"
		bp_name <- ", no clearcuts "
	} else {
		ylab_txt <- "mean values from raster including open areas"
		bp_name <- ""
	}
	main_txt_row1 <- paste( "shape vs random points, radius ", v_radius, sep="") 
	v_prefix_array <- unlist(str_split(v_prefix, "_"))
	main_txt_row2 <- ""
	for(ind in 2:length(v_prefix_array)  ) {
		main_txt_row2 <- paste(main_txt_row2, v_prefix_array[ind], sep=" ")
	}
	main_txt <- paste(main_txt_row1, main_txt_row2, sep="\n")
	bp_name <- paste(v_radius, bp_name, ".png", sep="")
	png(paste(v_dir_path, "boxplot_", v_radius, "_", v_cc_type, "_" , v_prefix, ".png", sep=""))
	plot.new()
	
	boxplot(v_shp_var, p100_var_noNA, ylab=ylab_txt, xlab=v_prefix_array[1], names=c("shapefile points", "random points"), main=main_txt, col="green", plot=TRUE)
	abline(h=v_mean_nocrop, col="red")
	abline(h=v_mean_crop, col="green")
	dev.off()
	mean_p100 <- mean(v_p100_var, na.rm = TRUE)
} # f_boxplot


# mean value for raster 
f_cs_mean = function(x) {
 cellStats(x, 'mean', na.rm=TRUE) 
}


# setting extent from geo-file in R-raster because readGrass6 returns raster with location extent 
set_extent = function(v_grassRaster, v_raster) {
 regexp_extent_N <- "north=[0-9]+"
 regexp_extent_S <- "south=[0-9]+"
 regexp_extent_W <- "west=[0-9]+"
 regexp_extent_E <- "east=[0-9]+"
	stat_str  <- execGRASS("r.info",  flags=c("g"), map=v_grassRaster, intern=TRUE)
	extent_W_vect <- unlist(str_split(regmatches(stat_str, regexpr(regexp_extent_E, stat_str)), "="))
	extent_E_vect <- unlist(str_split(regmatches(stat_str, regexpr(regexp_extent_W, stat_str)),  "="))
	extent_S_vect <- unlist(str_split(regmatches(stat_str, regexpr(regexp_extent_N, stat_str)), "=")) 
	extent_N_vect <- unlist(str_split(regmatches(stat_str, regexpr(regexp_extent_S, stat_str)), "="))
	temp_ext <- extent(as.numeric(extent_E_vect[2]), as.numeric(extent_W_vect[2]), as.numeric(extent_N_vect[2]), as.numeric(extent_S_vect[2]))

	is_equal <- all.equal(extent(v_raster), temp_ext)
	if( is_equal == FALSE) {	
	  extent(v_raster) <- temp_ext
	}

	return(v_raster)
} # set_extent


# setting bounding box from geo-file in Grass-raster because readGrass6 returns raster with location bounding box 
set_bbox = function(v_grassRaster, v_gridDataFrame) { 
	# ToDo: find an equal function
	regexp_extent_N <- "north=[0-9]+"
	regexp_extent_S <- "south=[0-9]+"
	regexp_extent_W <- "west=[0-9]+"
	regexp_extent_E <- "east=[0-9]+"
	stat_str  <- execGRASS("r.info",  flags=c("g"), map=v_grassRaster, intern=TRUE)
	extent_W_vect <- unlist(str_split(regmatches(stat_str, regexpr(regexp_extent_E, stat_str)), "="))
	extent_E_vect <- unlist(str_split(regmatches(stat_str, regexpr(regexp_extent_W, stat_str)),  "="))
	extent_S_vect <- unlist(str_split(regmatches(stat_str, regexpr(regexp_extent_N, stat_str)), "=")) 
	extent_N_vect <- unlist(str_split(regmatches(stat_str, regexpr(regexp_extent_S, stat_str)), "="))	
	
	v_gridDataFrame@bbox[1,2] <- as.numeric(extent_W_vect[2])
	v_gridDataFrame@bbox[1,1] <- as.numeric(extent_E_vect[2])
	v_gridDataFrame@bbox[2,2] <- as.numeric(extent_S_vect[2])
	v_gridDataFrame@bbox[2,1] <- as.numeric(extent_N_vect[2])

	return(v_gridDataFrame)
} # set_bbox


# logging system parameters and arguments 
print(paste(" file sep: ", file_sep))
print(paste(" file psep: ",file_psep))
print(paste(" OS type: ", OS_type))
print(paste(" OS version: ", os_version))
print(paste(" sysname: ", system_name))
print(paste(" sysinfo: ", system_info))     


# from command prompt
print(paste("input min value  ", in_MinLimit ))
print(paste("input map prefix: ", in_Prefix))
print(paste("input raster map: ", in_GrassRaster))
print(paste("input vector map, shapefile: ", in_ShapeFile))
print(paste("input vector map layer name, shapefile: ", in_ShpAttribute))
print(paste("input location: ", in_Location))
print(paste("input mapset: ", in_Mapset))


# setting Grass environment and logging result 
cat("\n")
init_grass <- initGRASS(gisBase=curr_gisBase, home = tempdir(), gisDbase=curr_gisDbase, location=in_Location, mapset=in_Mapset, override = TRUE)
cat(" ---  init ----------\n")
print(init_grass) 
cat(" ---  init end ---------- \n\n")


#making objects from Grass raster-map
rast_cla <- readRAST6(in_GrassRaster, useGDAL=TRUE)
cat(" \n -- set bbox --- \n")
print(rast_cla@bbox)
rast_cla <- set_bbox(in_GrassRaster, rast_cla)
print(rast_cla@bbox)
cat("  ------- end set bbox ----- \n")

rast_cla_r <- raster(rast_cla)

cat("\n")
print(paste("  ----- before set extent :  ", str(extent(rast_cla_r))))
rast_cla_r <- set_extent(in_GrassRaster, rast_cla_r)
print(paste(" ----- after set extent :  ", str(extent(rast_cla_r))))
cat("\n")


# making objects from the shape-file and an extent to crop the raster from raster-file 
shp <- readShapePoints(in_ShapeFile, proj4string = CRS(rt90), verbose = TRUE, repair = TRUE)

#shp extent
shp_ext <- extent(shp)

#creating rasters cropped after shape-file 
rast_cla_r_cr <- crop(rast_cla_r, shp_ext)


# mean values for raters, original and cropped to shapefile, used in boxplots
mean_nocrop <- f_cs_mean(rast_cla_r)
mean_crop <- f_cs_mean(rast_cla_r_cr)
cat(" \n ------ mean values for raters, original and cropped to shapefile --- \n")
print(paste("mean value from geo-file: ", mean_nocrop))
print(paste("mean value from geo-file cropped to shapefile", mean_crop))
cat(" ------- end mean values for raters, original and cropped to shapefile --- \n")


# removing values lover then given min-level
rast_cla_r_nocc_temp <- calc(rast_cla_r, set_min_value_to_zero)  # curr did crop
rast_cla_r_nocc <- calc(rast_cla_r_nocc_temp, set_zero_to_NA)

rast_cla_r_cr_nocc_temp <- calc(rast_cla_r_cr, set_min_value_to_zero)  # curr did crop
rast_cla_r_cr_nocc <- calc(rast_cla_r_cr_nocc_temp, set_zero_to_NA)


cat(" \n ------ histograms geo-file area and area cropped to shapefile and finally with open areas removed ------- \n")
png(paste(dir_path, "histogram_geofile_area _", in_Prefix, ".png", sep=""))
plot.new()
h0 <- hist(rast_cla_r, plot=TRUE, freq=TRUE)
print(h0)
dev.off()

png(paste(dir_path, "histogram_geofile_area_cropped_tp_shapefile_", in_Prefix, ".png", sep=""))
plot.new()
h1 <- hist(rast_cla_r_cr, plot=TRUE, freq=TRUE)
print(h1)
dev.off()

png(paste(dir_path, "histogram_cropped_after_remove_openareas_", in_Prefix, ".png", sep=""))
plot.new()
h2 <- hist(rast_cla_r_cr_nocc, plot=TRUE, freq=TRUE)
print(h2)
dev.off()

cat(" \n ------ end histograms geo-file area and area cropped to shapefile and finally with open areas removed ------- \n\n")

# objects from a number of random points, less then 100 
random_points <- spsample(rast_cla, 100, "random") # makes random  number of random points, less then 100 
random_spatial_points = SpatialPoints(coordinates(random_points), CRS(rt90))


# to save coordinates for making other diagramms from several runs
write.matrix(random_spatial_points@coords, paste(dir_path, "p100_coordinates.csv", sep=""), sep=",")


# plot with shape points versus random points saved as a .png
png(paste(dir_path, "shp_and_random_points_mean_", in_Prefix, ".png", sep=""))
plot.new()
plot(random_spatial_points, col="red")
plot(shp, col="green", add=TRUE)
dev.off()


# mean from shapefile points with different radiuses 
shp_mean_250m <-  round(extract(rast_cla_r_cr, shp, buffer=250,  fun=mean, na.rm=TRUE), 1)  # earlier 0 
shp_mean_500m <-  round(extract(rast_cla_r_cr, shp, buffer=500,  fun=mean, na.rm=TRUE), 1)
shp_mean_750m <-  round(extract(rast_cla_r_cr, shp, buffer=750,  fun=mean, na.rm=TRUE), 1)
shp_mean_1000m <- round(extract(rast_cla_r_cr, shp, buffer=1000, fun=mean, na.rm=TRUE), 1)


# mean valuea from  shapefile points without clearcuts, cc
shp_mean_250m_nocc <-  round(extract(rast_cla_r_cr_nocc, shp, buffer=250,  fun=mean, na.rm=TRUE), 1)
shp_mean_500m_nocc <-  round(extract(rast_cla_r_cr_nocc, shp, buffer=500,  fun=mean, na.rm=TRUE), 1)
shp_mean_750m_nocc <-  round(extract(rast_cla_r_cr_nocc, shp, buffer=750,  fun=mean, na.rm=TRUE), 1)
shp_mean_1000m_nocc <- round(extract(rast_cla_r_cr_nocc, shp, buffer=1000, fun=mean, na.rm=TRUE), 1)


# mean values from random number of sample points
p100_mean_250m <-  round(extract(rast_cla_r, random_spatial_points, buffer=250,  fun=mean, na.rm=TRUE), 1)
p100_mean_500m <-  round(extract(rast_cla_r, random_spatial_points, buffer=500,  fun=mean, na.rm=TRUE), 1)
p100_mean_750m <-  round(extract(rast_cla_r, random_spatial_points, buffer=750,  fun=mean, na.rm=TRUE), 1)
p100_mean_1000m <- round(extract(rast_cla_r, random_spatial_points, buffer=1000, fun=mean, na.rm=TRUE), 1)


# mean values from random number of sample points without clearcuts, nocc
p100_mean_250m_nocc <-  round(extract(rast_cla_r_nocc, random_spatial_points, buffer=250,  fun=mean, na.rm=TRUE), 1)
p100_mean_500m_nocc <-  round(extract(rast_cla_r_nocc, random_spatial_points, buffer=500,  fun=mean, na.rm=TRUE), 1)
p100_mean_750m_nocc <-  round(extract(rast_cla_r_nocc, random_spatial_points, buffer=750,  fun=mean, na.rm=TRUE), 1)
p100_mean_1000m_nocc <- round(extract(rast_cla_r_nocc, random_spatial_points, buffer=1000, fun=mean, na.rm=TRUE), 1)


# calculating means for different radiuses 
colmeans_shp_250m <- round(mean(shp_mean_250m, na.rm=TRUE),0)
shp_mean_vect[1] <- colmeans_shp_250m
colmeans_shp_500m <- round(mean(shp_mean_500m, na.rm=TRUE),0)
shp_mean_vect[2] <- colmeans_shp_500m
colmeans_shp_750m <- round(mean(shp_mean_750m, na.rm=TRUE),0)
shp_mean_vect[3] <- colmeans_shp_750m
colmeans_shp_1000m <- round(mean(shp_mean_1000m, na.rm=TRUE),0)
shp_mean_vect[4] <- colmeans_shp_1000m

colmeans_shp_250m_nocc <- round(mean(shp_mean_250m_nocc, na.rm=TRUE),0)
shp_mean_vect[5] <- colmeans_shp_250m_nocc
colmeans_shp_500m_nocc <- round(mean(shp_mean_500m_nocc, na.rm=TRUE),0)
shp_mean_vect[6] <- colmeans_shp_500m_nocc
colmeans_shp_750m_nocc <- round(mean(shp_mean_750m_nocc, na.rm=TRUE),0)
shp_mean_vect[7] <- colmeans_shp_750m_nocc
colmeans_shp_1000m_nocc <- round(mean(shp_mean_1000m_nocc, na.rm=TRUE),0)
shp_mean_vect[8] <- colmeans_shp_1000m_nocc
shp_mean_vect[is.na(shp_mean_vect)] <- 0  # 0 for missing values

colmeans_p100_250m <- round(mean(p100_mean_250m, na.rm=TRUE),0)
p100_mean_vect[1] <- colmeans_p100_250m
colmeans_p100_500m <- round(mean(p100_mean_500m, na.rm=TRUE),0)
p100_mean_vect[2] <- colmeans_p100_500m
colmeans_p100_750m <-round( mean(p100_mean_750m, na.rm=TRUE),0)
p100_mean_vect[3] <- colmeans_p100_750m
colmeans_p100_1000m <- round(mean(p100_mean_1000m, na.rm=TRUE),0)
p100_mean_vect[4] <- colmeans_p100_1000m

colmeans_p100_250m_nocc <- round(mean(p100_mean_250m_nocc, na.rm=TRUE),0)
p100_mean_vect[5] <- colmeans_p100_250m_nocc
colmeans_p100_500m_nocc <- round(mean(p100_mean_500m_nocc, na.rm=TRUE),0)
p100_mean_vect[6] <- colmeans_p100_500m_nocc
colmeans_p100_750m_nocc <- round(mean(p100_mean_750m_nocc, na.rm=TRUE),0)
p100_mean_vect[7] <- colmeans_p100_750m_nocc
colmeans_p100_1000m_nocc <- round(mean(p100_mean_1000m_nocc, na.rm=TRUE),0)
p100_mean_vect[8] <- colmeans_p100_1000m_nocc

p100_mean_vect[is.na(p100_mean_vect)] <- 0   # check


# create matrix with mean values for the random points and writes it to a csv-file
p100_mean_matrix <- matrix(p100_mean_vect, 1, length(p100_mean_vect))
colnames(p100_mean_matrix) <- tab_cols
shp_mean_matrix <- matrix(shp_mean_vect, 1, length(shp_mean_vect))
colnames(shp_mean_matrix) <- tab_cols
shp_vs_p100_mean_matrix <- rbind(shp_mean_matrix, p100_mean_matrix)
rownames(shp_vs_p100_mean_matrix) <- c("shp means", "p100 means")
write.csv(shp_vs_p100_mean_matrix, paste(dir_path, "means_shp_vs_p100_allradiuses",filename_divider,  in_Prefix, ".csv", sep=""), quote=FALSE)


# creating csv from a matrix of mean values from shape points
# ToDo: create variable 
shp_matrix_rows <- length(shp_mean_250m)
shp_mean_matrix <- matrix("", shp_matrix_rows, 0)
shp_mean_matrix <- cbind(shp_mean_matrix, shp_mean_250m)
shp_mean_matrix <- cbind(shp_mean_matrix, shp_mean_500m)
shp_mean_matrix <- cbind(shp_mean_matrix, shp_mean_750m)
shp_mean_matrix <- cbind(shp_mean_matrix, shp_mean_1000m)
shp_mean_matrix <- cbind(shp_mean_matrix, shp_mean_250m_nocc)
shp_mean_matrix <- cbind(shp_mean_matrix, shp_mean_500m_nocc)
shp_mean_matrix <- cbind(shp_mean_matrix, shp_mean_750m_nocc)
shp_mean_matrix <- cbind(shp_mean_matrix, shp_mean_1000m_nocc)

colnames(shp_mean_matrix) <- tab_cols
rownames(shp_mean_matrix) <- get_shape_rownames(shp, in_ShpAttribute)
shp_mean_matrix <- rbind(shp_mean_matrix, shp_mean_vect)
rownames(shp_mean_matrix)[nrow(shp_mean_matrix)] <- "mean value"
write.csv(shp_mean_matrix, paste(dir_path,  "shp_mean_matrix", filename_divider, in_Prefix, ".csv", sep=""), quote=FALSE)


#creating csv from a matrix of random points after creating a dataframe to create rownumbers
# the parameter stringsAsFactors=FALSE is needed to avoid random NA values in converting from matrix to dataframe
#  converting a matrix to a dataframe is a way to get rownumbers
p100_matrix_rows <- length(p100_mean_500m)
p100_mean_matrix <- matrix("", p100_matrix_rows, 0)
p100_mean_matrix <- cbind(p100_mean_matrix, p100_mean_250m)
p100_mean_matrix <- cbind(p100_mean_matrix, p100_mean_500m)
p100_mean_matrix <- cbind(p100_mean_matrix, p100_mean_750m)
p100_mean_matrix <- cbind(p100_mean_matrix, p100_mean_1000m)
p100_mean_matrix <- cbind(p100_mean_matrix, p100_mean_250m_nocc)
p100_mean_matrix <- cbind(p100_mean_matrix, p100_mean_500m_nocc)
p100_mean_matrix <- cbind(p100_mean_matrix, p100_mean_750m_nocc)
p100_mean_matrix <- cbind(p100_mean_matrix, p100_mean_1000m_nocc)

colnames(p100_mean_matrix) <- tab_cols
p100_mean_matrix <- rbind(p100_mean_matrix, p100_mean_vect)
df_p100_mean_matrix <- data.frame(p100_mean_matrix, stringsAsFactors=FALSE)  # dup rows
rownames(df_p100_mean_matrix)[nrow(df_p100_mean_matrix)] <- "mean value"
write.csv(df_p100_mean_matrix, paste(dir_path, "p100_mean_matrix_df", filename_divider, in_Prefix, ".csv", sep=""), quote=FALSE)


# variables with filenames 
shp_mean_tab_250m <-   paste(in_Prefix, filename_divider, "rast_mean_r_cr_shp_250_t",  ".tab", sep="")
shp_mean_tab_500m <-   paste(in_Prefix, filename_divider, "rast_mean_r_cr_shp_500_t",  ".tab", sep="")
shp_mean_tab_750m <-   paste(in_Prefix, filename_divider, "rast_mean_r_cr_shp_750_t",  ".tab", sep="")
shp_mean_tab_1000m <-  paste(in_Prefix, filename_divider, "rast_mean_r_cr_shp_1000_t", ".tab", sep="")

p100_mean_tab_250m <-  paste(in_Prefix, filename_divider, "p100_mean_250_m",  ".tab", sep="")
p100_mean_tab_500m <-  paste(in_Prefix, filename_divider, "p100_mean_500_m",  ".tab", sep="")
p100_mean_tab_750m <-  paste(in_Prefix, filename_divider, "p100_mean_750_m",  ".tab", sep="")
p100_mean_tab_1000m <- paste(in_Prefix, filename_divider, "p100_mean_1000_m", ".tab", sep="")


# creating dataframes from shape points as are easier to use later in statistical tests etc
df_shp_mean_250m  <- data.frame(shp_mean_250m)
df_shp_mean_500m  <- data.frame(shp_mean_500m)
df_shp_mean_750m  <- data.frame(shp_mean_750m)
df_shp_mean_1000m <- data.frame(shp_mean_1000m)
df_shp_mean_250m_nocc  <- data.frame(shp_mean_250m_nocc)
df_shp_mean_500m_nocc  <- data.frame(shp_mean_500m_nocc)
df_shp_mean_750m_nocc  <- data.frame(shp_mean_750m_nocc)
df_shp_mean_1000m_nocc <- data.frame(shp_mean_1000m_nocc)


# creating dataframes from random points as are easier to use later in statistical tests etc
df_p100_mean_250m  <- data.frame(p100_mean_250m)
df_p100_mean_500m  <- data.frame(p100_mean_500m)
df_p100_mean_750m  <- data.frame(p100_mean_750m)
df_p100_mean_1000m <- data.frame(p100_mean_1000m)
df_p100_mean_250m_nocc  <- data.frame(p100_mean_250m_nocc)
df_p100_mean_500m_nocc  <- data.frame(p100_mean_500m_nocc)
df_p100_mean_750m_nocc  <- data.frame(p100_mean_750m_nocc)
df_p100_mean_1000m_nocc <- data.frame(p100_mean_1000m_nocc)


# adding columns for the statistic result matrix to separaate random and shape points from eachother
df_shp_mean_250m$variable  <- "shp"
df_shp_mean_500m$variable  <- "shp"
df_shp_mean_750m$variable  <- "shp"
df_shp_mean_1000m$variable <- "shp"
df_shp_mean_250m_nocc$variable  <- "shp nocc"
df_shp_mean_500m_nocc$variable  <- "shp nocc"
df_shp_mean_750m_nocc$variable  <- "shp nocc"
df_shp_mean_1000m_nocc$variable <- "shp nocc"

df_p100_mean_250m$variable  <- "p100"
df_p100_mean_500m$variable  <- "p100"
df_p100_mean_750m$variable  <- "p100"
df_p100_mean_1000m$variable <- "p100"
df_p100_mean_250m_nocc$variable  <- "p100 nocc"
df_p100_mean_500m_nocc$variable  <- "p100 nocc"
df_p100_mean_750m_nocc$variable  <- "p100 nocc"
df_p100_mean_1000m_nocc$variable <- "p100 nocc"


#  statistical tests comparing series of shape points with random points 
# ks: two-sample Kolmogorov-Smirnov test , wc: two-sample Wilcoxon test
ks_2s_mean_250m <-   stat_test(df_shp_mean_250m, df_p100_mean_250m, "two.sid", "250m", "ks")
ks_gr_mean_250m <-   stat_test(df_shp_mean_250m, df_p100_mean_250m, "greater", "250m", "ks")
ks_less_mean_250m <- stat_test(df_shp_mean_250m, df_p100_mean_250m, "less",    "250m", "ks")

ks_2s_mean_500m <-   stat_test(df_shp_mean_500m, df_p100_mean_500m, "two.sid", "500m", "ks")
ks_gr_mean_500m <-   stat_test(df_shp_mean_500m, df_p100_mean_500m, "greater", "500m", "ks")
ks_less_mean_500m <- stat_test(df_shp_mean_500m, df_p100_mean_500m, "less",    "500m", "ks")

ks_2s_mean_750m <-   stat_test(df_shp_mean_750m, df_p100_mean_750m, "two.sid", "750m", "ks")
ks_gr_mean_750m <-   stat_test(df_shp_mean_750m, df_p100_mean_750m, "greater", "750m", "ks")
ks_less_mean_750m <- stat_test(df_shp_mean_750m, df_p100_mean_750m, "less",    "750m", "ks")

ks_2s_mean_1000m <-   stat_test(df_shp_mean_1000m, df_p100_mean_1000m, "two.sid", "1000m", "ks")
ks_gr_mean_1000m <-   stat_test(df_shp_mean_1000m, df_p100_mean_1000m, "greater", "1000m", "ks")
ks_less_mean_1000m <- stat_test(df_shp_mean_1000m, df_p100_mean_1000m, "less",    "1000m", "ks")

ks_2s_mean_250m_nocc <-   stat_test(df_shp_mean_250m_nocc, df_p100_mean_250m_nocc, "two.sid", "250m", "ks")
ks_gr_mean_250m_nocc <-   stat_test(df_shp_mean_250m_nocc, df_p100_mean_250m_nocc, "greater", "250m", "ks")
ks_less_mean_250m_nocc <- stat_test(df_shp_mean_250m_nocc, df_p100_mean_250m_nocc, "less",    "250m", "ks")

ks_2s_mean_500m_nocc <-   stat_test(df_shp_mean_500m_nocc, df_p100_mean_500m_nocc, "two.sid", "500m", "ks")
ks_gr_mean_500m_nocc <-   stat_test(df_shp_mean_500m_nocc, df_p100_mean_500m_nocc, "greater", "500m", "ks")
ks_less_mean_500m_nocc <- stat_test(df_shp_mean_500m_nocc, df_p100_mean_500m_nocc, "less",    "500m", "ks")

ks_2s_mean_750m_nocc <-   stat_test(df_shp_mean_750m_nocc, df_p100_mean_750m_nocc, "two.sid", "750m", "ks")
ks_gr_mean_750m_nocc <-   stat_test(df_shp_mean_750m_nocc, df_p100_mean_750m_nocc, "greater", "750m", "ks")
ks_less_mean_750m_nocc <- stat_test(df_shp_mean_750m_nocc, df_p100_mean_750m_nocc, "less",    "750m", "ks")

ks_2s_mean_1000m_nocc <-   stat_test(df_shp_mean_1000m_nocc, df_p100_mean_1000m_nocc, "two.sid", "1000m", "ks")
ks_gr_mean_1000m_nocc <-   stat_test(df_shp_mean_1000m_nocc, df_p100_mean_1000m_nocc, "greater", "1000m", "ks")
ks_less_mean_1000m_nocc <- stat_test(df_shp_mean_1000m_nocc, df_p100_mean_1000m_nocc, "less",    "1000m", "ks")


wc_2s_mean_250m <-   stat_test(df_shp_mean_250m, df_p100_mean_250m, "two.sid", "250m", "wc")
wc_gr_mean_250m <-   stat_test(df_shp_mean_250m, df_p100_mean_250m, "greater", "250m", "wc")
wc_less_mean_250m <- stat_test(df_shp_mean_250m, df_p100_mean_250m, "less",    "250m", "wc")

wc_2s_mean_500m <-   stat_test(df_shp_mean_500m, df_p100_mean_500m, "two.sid", "500m", "wc")
wc_gr_mean_500m <-   stat_test(df_shp_mean_500m, df_p100_mean_500m, "greater", "500m", "wc")
wc_less_mean_500m <- stat_test(df_shp_mean_500m, df_p100_mean_500m, "less",    "500m", "wc")

wc_2s_mean_750m <-   stat_test(df_shp_mean_750m, df_p100_mean_750m, "two.sid", "750m", "wc")
wc_gr_mean_750m <-   stat_test(df_shp_mean_750m, df_p100_mean_750m, "greater", "750m", "wc")
wc_less_mean_750m <- stat_test(df_shp_mean_750m, df_p100_mean_750m, "less",    "750m", "wc")

wc_2s_mean_1000m <-   stat_test(df_shp_mean_1000m, df_p100_mean_1000m, "two.sid", "1000m", "wc")
wc_gr_mean_1000m <-   stat_test(df_shp_mean_1000m, df_p100_mean_1000m, "greater", "1000m", "wc")
wc_less_mean_1000m <- stat_test(df_shp_mean_1000m, df_p100_mean_1000m,  "less",    "1000m", "wc")

wc_2s_mean_250m_nocc <-   stat_test(df_shp_mean_250m_nocc, df_p100_mean_250m_nocc, "two.sid", "250m", "wc")
wc_gr_mean_250m_nocc <-   stat_test(df_shp_mean_250m_nocc, df_p100_mean_250m_nocc, "greater", "250m", "wc")
wc_less_mean_250m_nocc <- stat_test(df_shp_mean_250m_nocc, df_p100_mean_250m_nocc, "less",    "250m", "wc")

wc_2s_mean_500m_nocc <-   stat_test(df_shp_mean_500m_nocc, df_p100_mean_500m_nocc, "two.sid", "500m", "wc")
wc_gr_mean_500m_nocc <-   stat_test(df_shp_mean_500m_nocc, df_p100_mean_500m_nocc, "greater", "500m", "wc")
wc_less_mean_500m_nocc <- stat_test(df_shp_mean_500m_nocc, df_p100_mean_500m_nocc, "less",    "500m", "wc")

wc_2s_mean_750m_nocc <-   stat_test(df_shp_mean_750m_nocc, df_p100_mean_750m_nocc, "two.sid", "750m", "wc")
wc_gr_mean_750m_nocc <-   stat_test(df_shp_mean_750m_nocc, df_p100_mean_750m_nocc, "greater", "750m", "wc")
wc_less_mean_750m_nocc <- stat_test(df_shp_mean_750m_nocc, df_p100_mean_750m_nocc, "less",    "750m", "wc")

wc_2s_mean_1000m_nocc <-   stat_test(df_shp_mean_1000m_nocc, df_p100_mean_1000m_nocc, "two.sid", "1000m", "wc")
wc_gr_mean_1000m_nocc <-   stat_test(df_shp_mean_1000m_nocc, df_p100_mean_1000m_nocc, "greater", "1000m", "wc")
wc_less_mean_1000m_nocc <- stat_test(df_shp_mean_1000m_nocc, df_p100_mean_1000m_nocc, "less",    "1000m", "wc")


#creating matrix with statistical testresults and saving them to a csv
stat_matr = matrix(ks_2s_mean_250m, 1, 9)
colnames(stat_matr) <- matr_cols

stat_matr <- rbind(stat_matr, ks_2s_mean_250m)
stat_matr <- rbind(stat_matr, ks_gr_mean_250m)
stat_matr <- rbind(stat_matr, ks_less_mean_250m)

stat_matr <- rbind(stat_matr, ks_2s_mean_500m)
stat_matr <- rbind(stat_matr, ks_gr_mean_500m)
stat_matr <- rbind(stat_matr, ks_less_mean_500m)

stat_matr <- rbind(stat_matr, ks_2s_mean_750m)
stat_matr <- rbind(stat_matr, ks_gr_mean_750m)
stat_matr <- rbind(stat_matr, ks_less_mean_750m)

stat_matr <- rbind(stat_matr, ks_2s_mean_1000m)
stat_matr <- rbind(stat_matr, ks_gr_mean_1000m)
stat_matr <- rbind(stat_matr, ks_less_mean_1000m)

stat_matr <- rbind(stat_matr, ks_2s_mean_250m_nocc)
stat_matr <- rbind(stat_matr, ks_gr_mean_250m_nocc)
stat_matr <- rbind(stat_matr, ks_less_mean_250m_nocc)

stat_matr <- rbind(stat_matr, ks_2s_mean_500m_nocc)
stat_matr <- rbind(stat_matr, ks_gr_mean_500m_nocc)
stat_matr <- rbind(stat_matr, ks_less_mean_500m_nocc)

stat_matr <- rbind(stat_matr, ks_2s_mean_750m_nocc)
stat_matr <- rbind(stat_matr, ks_gr_mean_750m_nocc)
stat_matr <- rbind(stat_matr, ks_less_mean_750m_nocc)

stat_matr <- rbind(stat_matr, ks_2s_mean_1000m_nocc)
stat_matr <- rbind(stat_matr, ks_gr_mean_1000m_nocc)
stat_matr <- rbind(stat_matr, ks_less_mean_1000m_nocc)

stat_matr <- rbind(stat_matr, wc_2s_mean_250m)
stat_matr <- rbind(stat_matr, wc_gr_mean_250m)
stat_matr <- rbind(stat_matr, wc_less_mean_250m)

stat_matr <- rbind(stat_matr, wc_2s_mean_500m)
stat_matr <- rbind(stat_matr, wc_gr_mean_500m)
stat_matr <- rbind(stat_matr, wc_less_mean_500m)

stat_matr <- rbind(stat_matr, wc_2s_mean_750m)
stat_matr <- rbind(stat_matr, wc_gr_mean_750m)
stat_matr <- rbind(stat_matr, wc_less_mean_750m)

stat_matr <- rbind(stat_matr, wc_2s_mean_1000m)
stat_matr <- rbind(stat_matr, wc_gr_mean_1000m)
stat_matr <- rbind(stat_matr, wc_less_mean_1000m)

stat_matr <- rbind(stat_matr, wc_2s_mean_250m_nocc)
stat_matr <- rbind(stat_matr, wc_gr_mean_250m_nocc)
stat_matr <- rbind(stat_matr, wc_less_mean_250m_nocc)

stat_matr <- rbind(stat_matr, wc_2s_mean_500m_nocc)
stat_matr <- rbind(stat_matr, wc_gr_mean_500m_nocc)
stat_matr <- rbind(stat_matr, wc_less_mean_500m_nocc)

stat_matr <- rbind(stat_matr, wc_2s_mean_750m_nocc)
stat_matr <- rbind(stat_matr, wc_gr_mean_750m_nocc)
stat_matr <- rbind(stat_matr, wc_less_mean_750m_nocc)

stat_matr <- rbind(stat_matr, wc_2s_mean_1000m_nocc)
stat_matr <- rbind(stat_matr, wc_gr_mean_1000m_nocc)
stat_matr <- rbind(stat_matr, wc_less_mean_1000m_nocc)


# matrix with statistical data
matrix_name <- paste("stattest_means_", in_Prefix, ".csv", sep="")
write.matrix(stat_matr, paste(dir_path, matrix_name, sep=""), sep=",")

# matrix with statistical data, ordered
stat_matr2 <- stat_matr[order(stat_matr[,8], stat_matr[,7], decreasing=FALSE, stat_matr[,6]),]
matrix_name2 <- paste("stattest_means_sorted_", in_Prefix, ".csv", sep="")
write.matrix(stat_matr2, paste(dir_path, matrix_name2, sep=""), sep=",")

# matrix with statistical data, only significant 
stat_matr3 <- matrix(stat_matr[1,], 1, 9)
colnames(stat_matr3) <- matr_cols
for(ind in 2:length(stat_matr[,1])) {
	if( (stat_matr[ind,7] < 0.06) && ((stat_matr[ind,9] == "two-sided") || (stat_matr[ind,9] == "two.sided")) ) {  
		stat_matr3 <- rbind(stat_matr3, stat_matr[ind,])
	}
} # for
matrix_name3 <- paste("stattest_means_relevant_", in_Prefix, ".csv", sep="")
write.matrix(stat_matr3, paste(dir_path, matrix_name3, sep=""), sep=",")


# boxplots comparing shp and random points 
f_boxplot(shp_mean_250m,  p100_mean_250m,  "250m",  "cc", mean_crop, mean_nocrop, in_Prefix, dir_path)
f_boxplot(shp_mean_500m,  p100_mean_500m,  "500m",  "cc", mean_crop, mean_nocrop, in_Prefix, dir_path)
f_boxplot(shp_mean_750m,  p100_mean_750m,  "750m",  "cc", mean_crop, mean_nocrop, in_Prefix, dir_path)
f_boxplot(shp_mean_1000m, p100_mean_1000m, "1000m", "cc", mean_crop, mean_nocrop, in_Prefix, dir_path)
f_boxplot(shp_mean_250m_nocc,  p100_mean_250m_nocc,   "250m", "nocc", mean_crop, mean_nocrop, in_Prefix, dir_path)
f_boxplot(shp_mean_500m_nocc,  p100_mean_500m_nocc,   "500m", "nocc", mean_crop, mean_nocrop, in_Prefix, dir_path)
f_boxplot(shp_mean_750m_nocc,  p100_mean_750m_nocc,   "750m", "nocc", mean_crop, mean_nocrop, in_Prefix, dir_path)
f_boxplot(shp_mean_1000m_nocc, p100_mean_1000m_nocc, "1000m", "nocc", mean_crop, mean_nocrop, in_Prefix, dir_path)

v_warn <- warnings()
cat(" last warnings: \n")
print(v_warn)
rm(list = ls())



