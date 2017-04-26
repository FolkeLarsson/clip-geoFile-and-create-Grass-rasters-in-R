## runs from Grass command prompt, called by Rscript
## calculate parts of classes in circles with different radiuses around shape-points and testing them against random points
## for testing, the non-parametric method wilcox.test is used 
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


# command line aruments
in_Prefix              <- args[1] #directory name for both logging dir and csv-files etc
in_GrassRasterClasses  <- args[2] # name of raster map with classes from Grass 
in_ShapeFile           <- args[3] # # path and name of shape file 
in_Location            <- args[4]  # new or existing Grass location
in_Mapset              <- args[5]  # existing Grass mapset, if new localion always PERMANENT
in_ShpAttribute        <- args[6] # name of attribute in shape layer, not mandatory if name is in shape file 

#variables
stat_matr = c() # statistic data


#constants
regexp_classes <- "\\|[0-9]\\|[0-9a-zA-Z]{1,}"
rt90 =   "+proj=tmerc +lat_0=0 +lon_0=15.80827777777778 +k=1 +x_0=1500000 +y_0=0 +ellps=bessel +units=m +no_defs"
matr_cols = c("var1", "nVar1", "var2", "nVar2", "D/W", "p", "alternative", "method", "radius")
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


# extracting a vector of column fields of classes from the classified raster, classes originated from classfile 
# extracting headers from original classfile is easier and gives the same result if not modified or no errors occurs, 
# therefore this method is safer 
get_header_vect = function( v_header_text, v_regexpr) {
 v_header_vect <- c()
 v_reg_exp <- regexpr(v_regexpr, v_header_text)
 v_str_mat <- regmatches(v_header_text, v_reg_exp)
 v_header_l <- length(v_str_mat)
 for (ind in 1:v_header_l) {
  v_curr_head <- substr(v_str_mat[ind], 2, str_length(v_str_mat[ind]))
  v_curr_head_l <- str_length(v_curr_head)
  v_cut_pos <- str_locate(v_curr_head, "\\|")
  v_curr_head <- substr(v_curr_head, v_cut_pos+1, v_curr_head_l)
  v_header_vect[ind] <- v_curr_head
 }
 return(v_header_vect)
} # get_header_vect

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



# wrapper function for statistical tests, now Kolmogorov-Smirnov(ks) and Wilcoxon-Smith(ws) two-sample tests
stat_test = function(v_var1, v_var2, v_alternative, v_radius, v_method, v_limit) {
 stat_vect = c()  # return vector
 v_var1_name <- v_var1$variable[1]
 v_var2_name <- v_var2$variable[1]
 v_var1_noNA <- na.omit(v_var1)
 v_var2_noNA <- na.omit(v_var2)
 v_var1_length <- length(v_var1[,1])
 v_var2_length <- length(v_var2[,1])
 v_var1_noNA_length <- length(v_var1_noNA[,1])
 v_var2_noNA_length <- length(v_var2_noNA[,1])
 
 if ((v_var1_noNA_length > 1) && (v_var2_noNA_length > 1)) {
  if(v_method == "ks") {
   stat_res <- ks.test(v_var1_noNA[,1], v_var2_noNA[,1], alternative=v_alternative)
  } else if(v_method == "wc") {
   stat_res <- wilcox.test(v_var1_noNA[,1], v_var2_noNA[,1], alternative=v_alternative)
  } else {
    print(v_alternative)
    print(v_method)
    print(" NA NA NA no such method! ")
    stat_vect <- NA  # should be caught
  }

  v_dataname <- stat_vect$data.name
  v_dataname_l <- str_length(v_dataname)
  v_cut_pos  <- str_locate(v_dataname, " and ")
  v_dataname_1 <- substr(v_dataname, 1, v_cut_pos-1)
  v_dataname_2 <- substr(v_dataname, v_cut_pos+5, v_dataname_l)

  stat_vect[1] = v_var1_name
  stat_vect[2] = v_var1_noNA_length
  stat_vect[3] = v_var2_name
  stat_vect[4] = v_var2_noNA_length
  stat_vect[5] = v_radius
  stat_vect[6] = round(stat_res$statistic, 2)
  stat_vect[7] = round(stat_res$p.value, 2)
  stat_vect[8] = stat_res$method
  stat_vect[9] = stat_res$alternative
  return(stat_vect)
 } else {
  stop(" STOP and bel�gg, to few variables not NA for statistical test ")
 }
} #stat_test


# setting bounding box from geo-file in Grass-raster because readGrass6 returns raster with location bounding box 
set_bbox = function(v_grassRaster, v_gridDataFrame) { 
	# find an equal function
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


# logging system parameters and arguments 
print(paste(" file sep: ", file_sep))
print(paste(" file psep: ",file_psep))
print(paste(" OS type: ", OS_type))
print(paste(" OS version: ", os_version))
print(paste(" sysname: ", system_name))
print(paste(" sysinfo: ", system_info))     


# from command prompt
print(paste("input map prefix: ", in_Prefix))
print(paste("input reclassed raster map: ", in_GrassRasterClasses))
print(paste("input vector map, shapefile: ", in_ShapeFile))
print(paste("input vector map layer name, shapefile: ", in_ShpAttribute))
print(paste("input location: ", in_Location))
print(paste("input mapset: ", in_Mapset))

cat("\n")
init_grass <- initGRASS(gisBase=curr_gisBase, home = tempdir(), gisDbase=curr_gisDbase, location=in_Location, mapset=in_Mapset, override = TRUE)
cat(" ---  init ----------\n")
print(init_grass) 
cat(" ---  init end ----------\n\n")

# making objects from reclassed Grass raster-map
rast_cla <- readRAST6(in_GrassRasterClasses, useGDAL=TRUE)
rast_cla <- set_bbox(in_GrassRasterClasses, rast_cla)
rast_cla_r <- raster(rast_cla)
rast_cla_r <- set_extent(in_GrassRasterClasses, rast_cla_r)

# statistics about categories and their distribution
rep_str <- execGRASS("r.report", flags=c("h", "C"), map=in_GrassRasterClasses, intern=TRUE)
stat_str <- execGRASS("r.stats", flags=c("c"), input=in_GrassRasterClasses, intern=TRUE)
cat("\n--- statistics -------\n")
print(rep_str)
print(stat_str)
cat("\n--- end statistics ---\n")

# getting headers( class-names ) from reclassed raster using regular expressions
headers <- get_header_vect(as.character(rep_str), regexp_classes)

# making objects from the shape-file and an extent to crop the raster from raster-file 
shp <- readShapePoints(in_ShapeFile, proj4string = CRS(rt90))
shp_ext <- extent(shp)

#creating rasters cropped after shape-file 
rast_cla_r_cr <- crop(rast_cla_r, shp_ext)

# create random  number of random points, less then 100 
sp_100p <- spsample(rast_cla, 100, "random") 

# objects from a number of random points, less then 100 
spatial_points <- SpatialPointsDataFrame(coordinates(sp_100p),data.frame(id=1:nrow(sp_100p@coords)))
random_spatial_points = SpatialPoints(coordinates(sp_100p), CRS(rt90))

# plot with shape points versus random points saved as a .png
png(paste(dir_path, "shp_and_random_points.png", sep=""))
plot.new()
plot(spatial_points, type="p", col="red")
plot(shp, type="p", col="green", add=TRUE)
dev.off()

# class tables from the random points 
#p100_cla_50 <-  extract(rast_cla_r, random_spatial_points, buffer=125) # not with KNN-data
p100_cla_250 <-  extract(rast_cla_r, random_spatial_points, buffer=250)
p100_cla_500 <-  extract(rast_cla_r, random_spatial_points, buffer=500)
p100_cla_750 <-  extract(rast_cla_r, random_spatial_points, buffer=750)
p100_cla_1000 <- extract(rast_cla_r, random_spatial_points, buffer=1000)

#p100_cla_50_t <-  lapply(p100_cla_50, table)  # not with KNN-data
p100_cla_250_t <-  lapply(p100_cla_250, table)
p100_cla_500_t <-  lapply(p100_cla_500, table)
p100_cla_750_t <-  lapply(p100_cla_750, table)
p100_cla_1000_t <- lapply(p100_cla_1000, table)


cat(" ------ histograms geo-file area and cropped area ------- \n")
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
cat(" \n ------ end histograms geo-file area and cropped area  ------ \n")


# classes table names, shp and random points
#shp_classes_tab_50m <-   paste(dir_path,  "rast_cla_r_cr_shp_50m",  filename_divider,   in_Prefix, ".csv", sep="")  # not with KNN-data
shp_classes_tab_250m <-   paste(dir_path,  "rast_cla_r_cr_shp_250m",   filename_divider,   in_Prefix, ".csv", sep="")
shp_classes_tab_500m <-   paste(dir_path,  "rast_cla_r_cr_shp_500m",   filename_divider,   in_Prefix, ".csv", sep="")
shp_classes_tab_750m <-   paste(dir_path, "rast_cla_r_cr_shp_750m",    filename_divider,   in_Prefix, ".csv", sep="")
shp_classes_tab_1000m <-  paste(dir_path,  "rast_cla_r_cr_shp_1000m",  filename_divider,   in_Prefix, ".csv", sep="")

#p100_classes_tab_50m <-  paste(dir_path,  "rast_cla_r_cr_p100_50m",  filename_divider,   in_Prefix, ".csv", sep="")   # not with KNN-data
p100_classes_tab_250m <-  paste(dir_path,  "rast_cla_r_cr_p100_250m",   filename_divider,   in_Prefix, ".csv", sep="")
p100_classes_tab_500m <-  paste(dir_path,  "rast_cla_r_cr_p100_500m",   filename_divider,   in_Prefix, ".csv", sep="")
p100_classes_tab_750m <-  paste(dir_path,  "rast_cla_r_cr_p100_750m",   filename_divider,   in_Prefix, ".csv", sep="")
p100_classes_tab_1000m <- paste(dir_path,  "rast_cla_r_cr_p100_1000m",  filename_divider,   in_Prefix, ".csv", sep="")


# class tables from the shape file 
#rast_cla_r_cr_shp_50 <-  extract(rast_cla_r_cr, shp, buffer=125)  # not with KNN-data
rast_cla_r_cr_shp_250 <-  extract(rast_cla_r_cr, shp, buffer=250)
rast_cla_r_cr_shp_500 <-  extract(rast_cla_r_cr, shp, buffer=500)
rast_cla_r_cr_shp_750 <-  extract(rast_cla_r_cr, shp, buffer=750)
rast_cla_r_cr_shp_1000 <- extract(rast_cla_r_cr, shp, buffer=1000)


#rast_cla_r_cr_shp_50_t <-  lapply(rast_cla_r_cr_shp_50, table)  # not with KNN-data
rast_cla_r_cr_shp_250_t <-  lapply(rast_cla_r_cr_shp_250, table)
rast_cla_r_cr_shp_500_t <-  lapply(rast_cla_r_cr_shp_500, table)
rast_cla_r_cr_shp_750_t <-  lapply(rast_cla_r_cr_shp_750, table)
rast_cla_r_cr_shp_1000_t <- lapply(rast_cla_r_cr_shp_1000, table)

write.matrix(rast_cla_r_cr_shp_250_t, file=as.character(shp_classes_tab_250m), sep=",")
table_classes_shp_250m <- read.table(as.character(shp_classes_tab_250m), fill=TRUE, col.names = headers, sep=",")
write.matrix(rast_cla_r_cr_shp_500_t, file=as.character(shp_classes_tab_500m), sep=",")
table_classes_shp_500m <- read.table(as.character(shp_classes_tab_500m), fill=TRUE, col.names = headers, sep=",")
write.matrix(rast_cla_r_cr_shp_750_t, file=as.character(shp_classes_tab_750m), sep=",")
table_classes_shp_750m <- read.table(as.character(shp_classes_tab_750m), fill=TRUE, col.names = headers, sep=",")
write.matrix(rast_cla_r_cr_shp_1000_t, file=as.character(shp_classes_tab_1000m), sep=",")
table_classes_shp_1000m <- read.table(as.character(shp_classes_tab_1000m), fill=TRUE, col.names = headers, sep=",")


# dataframes/tables for classes from random points for radius 500, 750, 1000 m 
#write.matrix(p100_cla_50_t, file=as.character(p100_classes_tab_50m), sep=",")  # not with KNN-data
#table_classes_p100_50m <- read.table(as.character(p100_classes_tab_50m), fill=TRUE, col.names = headers, sep=",")  # not with KNN-data
write.matrix(p100_cla_250_t, file=as.character(p100_classes_tab_250m), sep=",")
table_classes_p100_250m <- read.table(as.character(p100_classes_tab_250m), fill=TRUE, col.names = headers, sep=",")
write.matrix(p100_cla_500_t, file=as.character(p100_classes_tab_500m), sep=",")
table_classes_p100_500m <- read.table(as.character(p100_classes_tab_500m), fill=TRUE, col.names = headers, sep=",")
write.matrix(p100_cla_750_t, file=as.character(p100_classes_tab_750m), sep=",")
table_classes_p100_750m <- read.table(as.character(p100_classes_tab_750m), fill=TRUE, col.names = headers, sep=",")
write.matrix(p100_cla_1000_t, file=as.character(p100_classes_tab_1000m), sep=",")
table_classes_p100_1000m <- read.table(as.character(p100_classes_tab_1000m), fill=TRUE, col.names = headers, sep=",")


# adding column names from class file to dataframes for shapefile 
#colnames(table_classes_shp_50m)  <-  headers  # not with KNN-data
colnames(table_classes_shp_250m)  <-  headers
colnames(table_classes_shp_500m)  <-  headers
colnames(table_classes_shp_750m)  <-  headers
colnames(table_classes_shp_1000m) <-  headers

# adding column names from class file to dataframes for random points
#colnames(table_classes_p100_50m) <-  headers   # not with KNN-data
colnames(table_classes_p100_250m) <-  headers
colnames(table_classes_p100_500m) <-  headers
colnames(table_classes_p100_750m) <-  headers
colnames(table_classes_p100_1000m) <- headers

#rownames(table_classes_shp_50m) <- get_shape_rownames(shp, in_ShpAttribute)   # eval(parse(text=shp_name_full[2]))  # not with KNN-data
rownames(table_classes_shp_250m) <-  get_shape_rownames(shp, in_ShpAttribute)  #eval(parse(text=shp_name_full[2]))
rownames(table_classes_shp_500m) <-  get_shape_rownames(shp, in_ShpAttribute)  # eval(parse(text=shp_name_full[2]))   #eval(parse(text=shp_name_full))
rownames(table_classes_shp_750m) <-  get_shape_rownames(shp, in_ShpAttribute)  # eval(parse(text=shp_name_full[2])) #eval(parse(text=shp_name_full))
rownames(table_classes_shp_1000m) <- get_shape_rownames(shp, in_ShpAttribute) # eval(parse(text=shp_name_full[2])) #eval(parse(text=shp_name_full))

# adding a new column telling source of data for the statistic test matrix 
#table_classes_shp_50m$variable   <-  "shape "   # not with KNN-data
table_classes_shp_250m$variable   <-  "shape "
table_classes_shp_500m$variable   <-  "shape "
table_classes_shp_750m$variable   <-  "shape "
table_classes_shp_1000m$variable  <-  "shape "

#table_classes_p100_50m$variable  <-  "100p " # not with KNN-data
table_classes_p100_250m$variable  <-  "100p "
table_classes_p100_500m$variable  <-  "100p "
table_classes_p100_750m$variable  <-  "100p "
table_classes_p100_1000m$variable <-  "100p "


# statistical tests comparing series of shape points with random points 
# ks: two-sample Kolmogorov-Smirnov test , wc: two-sample Wilcoxon test
# iterating with ind representing the individual class from the classification
# saved as CSV-files 
for(ind in 1:length(headers)) {
 #shp_var_50m <- table_classes_shp_50m[ind]   # not with KNN-data	
 shp_var_250m <-  table_classes_shp_250m[ind]
 shp_var_500m <-  table_classes_shp_500m[ind]
 shp_var_750m <-  table_classes_shp_750m[ind]
 shp_var_1000m <- table_classes_shp_1000m[ind]

 #p100_var_50m <- table_classes_p100_50m[ind]  # not with KNN-data
 p100_var_250m <-  table_classes_p100_250m[ind]
 p100_var_500m <-  table_classes_p100_500m[ind]
 p100_var_750m <-  table_classes_p100_750m[ind]
 p100_var_1000m <- table_classes_p100_1000m[ind]
 
 #shp_var_50m$variable <- paste(table_classes_shp_500m$variable, names(shp_var_50m))  # not with KNN-data
 shp_var_250m$variable  <- paste(table_classes_shp_500m$variable, names(shp_var_250m))
 shp_var_500m$variable  <- paste(table_classes_shp_500m$variable, names(shp_var_500m))
 shp_var_750m$variable  <- paste(table_classes_shp_500m$variable, names(shp_var_750m))
 shp_var_1000m$variable <- paste(table_classes_shp_500m$variable, names(shp_var_1000m))
 #p100_var_50m$variable <- paste(table_classes_p100_500m[1]$variable, names(p100_var_50m))  # not with KNN-data
 p100_var_250m$variable  <- paste(table_classes_p100_500m[1]$variable, names(p100_var_250m))
 p100_var_500m$variable  <- paste(table_classes_p100_500m[1]$variable, names(p100_var_500m))
 p100_var_750m$variable  <- paste(table_classes_p100_500m[1]$variable, names(p100_var_750m))
 p100_var_1000m$variable <- paste(table_classes_p100_500m[1]$variable, names(p100_var_1000m))

 ks_2s_cla2 <- stat_test(shp_var_500m, p100_var_500m, "two.sid", "500m", "ks", 0)
 ks_less_cla2 <- stat_test(shp_var_500m, p100_var_500m, "less", "500m", "ks", 0)
 ks_greater_cla2 <- stat_test(shp_var_500m, p100_var_500m, "greater", "500m", "ks", 0)
 wc_2s_cla2_250m <- stat_test(shp_var_250m, p100_var_250m, "two.sid", "250m", "wc", 0)
 wc_less_cla2_250m <- stat_test(shp_var_250m, p100_var_250m, "less", "250m", "wc", 0)
 wc_greater_cla2_250m <- stat_test(shp_var_250m, p100_var_250m, "greater", "250m", "wc", 0)
 wc_2s_cla2 <- stat_test(shp_var_500m, p100_var_500m, "two.sid", "500m", "wc", 0)
 wc_less_cla2 <- stat_test(shp_var_500m, p100_var_500m, "less", "500m", "wc", 0)
 wc_greater_cla2 <- stat_test(shp_var_500m, p100_var_500m, "greater", "500m", "wc", 0)
 wc_2s_cla2_750m <- stat_test(shp_var_750m, p100_var_750m, "two.sid", "750m", "wc", 0)
 wc_less_cla2_750m <- stat_test(shp_var_750m, p100_var_750m, "less", "750m", "wc", 0)
 wc_greater_cla2_750m <- stat_test(shp_var_750m, p100_var_750m, "greater", "750m", "wc", 0)
 wc_2s_cla2_1000m <- stat_test(shp_var_1000m, p100_var_1000m, "two.sid", "1000m", "wc", 0)
 wc_less_cla2_1000m <- stat_test(shp_var_1000m, p100_var_1000m, "less", "1000m", "wc", 0)
 wc_greater_cla2_1000m <- stat_test(shp_var_1000m, p100_var_1000m, "greater", "1000m", "wc", 0)

 if (ind == 1) {
  stat_matr = matrix(ks_2s_cla2, 1, 9)
  dimnames(stat_matr) <- list(NULL, matr_cols)
  stat_matr <- rbind(stat_matr, ks_greater_cla2)
  stat_matr <- rbind(stat_matr, ks_less_cla2)
  stat_matr <- rbind(stat_matr, wc_2s_cla2_250m)
  stat_matr <- rbind(stat_matr, wc_less_cla2_250m)
  stat_matr <- rbind(stat_matr, wc_greater_cla2_250m)
  stat_matr <- rbind(stat_matr, wc_2s_cla2)
  stat_matr <- rbind(stat_matr, wc_greater_cla2)
  stat_matr <- rbind(stat_matr, wc_less_cla2)
  stat_matr <- rbind(stat_matr, wc_2s_cla2_750m)
  stat_matr <- rbind(stat_matr, wc_less_cla2_750m)
  stat_matr <- rbind(stat_matr, wc_greater_cla2_750m)
  stat_matr <- rbind(stat_matr, wc_2s_cla2_1000m)
  stat_matr <- rbind(stat_matr, wc_less_cla2_1000m)
  stat_matr <- rbind(stat_matr, wc_greater_cla2_1000m)

 } else {
  stat_matr <- rbind(stat_matr, wc_2s_cla2_250m)
  stat_matr <- rbind(stat_matr, wc_less_cla2_250m)
  stat_matr <- rbind(stat_matr, wc_greater_cla2_250m)
  stat_matr <- rbind(stat_matr, wc_2s_cla2)
  stat_matr <- rbind(stat_matr, wc_greater_cla2)
  stat_matr <- rbind(stat_matr, wc_less_cla2)
  stat_matr <- rbind(stat_matr, wc_2s_cla2_750m)
  stat_matr <- rbind(stat_matr, wc_less_cla2_750m)
  stat_matr <- rbind(stat_matr, wc_greater_cla2_750m)
  stat_matr <- rbind(stat_matr, wc_2s_cla2_1000m)
  stat_matr <- rbind(stat_matr, wc_less_cla2_1000m)
  stat_matr <- rbind(stat_matr, wc_greater_cla2_1000m)
 }
 #print(kruskal.test(shp_var_50m))  # not with KNN-data
 print(kruskal.test(shp_var_250m))
 print(kruskal.test(shp_var_500m))
 print(kruskal.test(shp_var_750m))
 print(kruskal.test(shp_var_1000m))
 

 # finally writing the matrix to a file
 if(ind == length(headers)) {
  matrix_name <- paste("stattest_classes_", in_Prefix, ".csv", sep="")
  matrix_name2 <- paste("stattest_classes_sorted_", in_Prefix, ".csv", sep="")
  write.matrix(stat_matr, paste(dir_path, matrix_name, sep=""), sep=",")

stat_matr2 <- stat_matr[order(stat_matr[,8], stat_matr[,7], decreasing=FALSE, stat_matr[,5]),]
  write.matrix(stat_matr2, paste(dir_path, matrix_name2, sep=""), sep=",")
 }
}


# removong the "shape", "100p" column
#table_classes_shp_50m <-   subset(table_classes_shp_50m, select = -variable)  # not with KNN-data
table_classes_shp_250m <-   subset(table_classes_shp_250m, select = -variable)
table_classes_shp_500m <-   subset(table_classes_shp_500m, select = -variable)
table_classes_shp_750m <-   subset(table_classes_shp_750m, select = -variable)
table_classes_shp_1000m <-  subset(table_classes_shp_1000m, select = -variable)

#table_classes_p100_50m <-  subset(table_classes_p100_50m, select = -variable)  # not with KNN-data
table_classes_p100_250m <-  subset(table_classes_p100_250m, select = -variable)
table_classes_p100_500m <-  subset(table_classes_p100_500m, select = -variable)
table_classes_p100_750m <-  subset(table_classes_p100_750m, select = -variable)
table_classes_p100_1000m <- subset(table_classes_p100_1000m, select = -variable)


# create new column: colsums with calculated values 
#table_classes_shp_50m$colsums <- apply(table_classes_shp_50m, 1, function(x) round(sum(x, na.rm=TRUE), digits=0))  # not with KNN-data
table_classes_shp_250m$colsums  <- apply(table_classes_shp_250m, 1, function(x) round(sum(x, na.rm=TRUE), digits=0))
table_classes_shp_500m$colsums  <- apply(table_classes_shp_500m, 1, function(x) round(sum(x, na.rm=TRUE), digits=0))
table_classes_shp_750m$colsums  <- apply(table_classes_shp_750m, 1, function(x) round(sum(x, na.rm=TRUE), digits=0))
table_classes_shp_1000m$colsums <- apply(table_classes_shp_1000m, 1, function(x) round(sum(x, na.rm=TRUE), digits=0))

#table_classes_p100_50m$colsums <- apply(table_classes_p100_50m, 1, function(x) round(sum(x, na.rm=TRUE), digits=0))  # not with KNN-data
table_classes_p100_250m$colsums  <- apply(table_classes_p100_250m, 1, function(x) round(sum(x, na.rm=TRUE), digits=0))
table_classes_p100_500m$colsums  <- apply(table_classes_p100_500m, 1, function(x) round(sum(x, na.rm=TRUE), digits=0))
table_classes_p100_750m$colsums  <- apply(table_classes_p100_750m, 1, function(x) round(sum(x, na.rm=TRUE), digits=0))
table_classes_p100_1000m$colsums <- apply(table_classes_p100_1000m, 1, function(x) round(sum(x, na.rm=TRUE), digits=0))

# as headers but with sum column 
headers2 <- colnames(table_classes_shp_250m)


#v_colmeans_shp_50m <- round(colMeans(table_classes_shp_50m, na.rm =TRUE),0)  # not with KNN-data
v_colmeans_shp_250m  <- round(colMeans(table_classes_shp_250m,  na.rm =TRUE),0)
v_colmeans_shp_500m  <- round(colMeans(table_classes_shp_500m,  na.rm =TRUE),0)
v_colmeans_shp_750m  <- round(colMeans(table_classes_shp_750m,  na.rm =TRUE), 0)
v_colmeans_shp_1000m <- round(colMeans(table_classes_shp_1000m, na.rm =TRUE), 0)

#v_colmeans_p100_50m <-  round(colMeans(table_classes_p100_50m, na.rm =TRUE), digits=0)  # not with KNN-data
v_colmeans_p100_250m <-  round(colMeans(table_classes_p100_250m,  na.rm =TRUE), digits=0) 
v_colmeans_p100_500m <-  round(colMeans(table_classes_p100_500m,  na.rm =TRUE), digits=0) 
v_colmeans_p100_750m <-  round(colMeans(table_classes_p100_750m,  na.rm =TRUE), digits=0)
v_colmeans_p100_1000m <- round(colMeans(table_classes_p100_1000m, na.rm =TRUE), digits=0)


#table_classes_shp_500m <- rbind(table_classes_shp_500m, v_colmeans_shp_500m)
table_classes_shp_250m  <- rbind(table_classes_shp_250m,  v_colmeans_shp_250m)
table_classes_shp_500m  <- rbind(table_classes_shp_500m,  v_colmeans_shp_500m)
table_classes_shp_750m  <- rbind(table_classes_shp_750m,  v_colmeans_shp_750m)
table_classes_shp_1000m <- rbind(table_classes_shp_1000m, v_colmeans_shp_1000m)

# setting name in last row
#row.names(table_classes_shp_50m)[nrow(table_classes_shp_50m)]    <- "AVG"   # not with KNN-data
row.names(table_classes_shp_250m)[nrow(table_classes_shp_250m)]   <- "AVG"
row.names(table_classes_shp_500m)[nrow(table_classes_shp_500m)]   <- "AVG"
row.names(table_classes_shp_750m)[nrow(table_classes_shp_750m)]   <- "AVG"
row.names(table_classes_shp_1000m)[nrow(table_classes_shp_1000m)] <- "AVG"


# saving CSV's of shp and p100 tables with all points
#write.csv(table_classes_shp_50m,   paste(dir_path, "table_classes_shp_50m.csv"), quote=FALSE)   # not with KNN-data
write.csv(table_classes_shp_250m,   paste(dir_path, "table_classes_shp_250m",  filename_divider,   in_Prefix, ".csv", sep=""), quote=FALSE)
write.csv(table_classes_shp_500m,   paste(dir_path, "table_classes_shp_500m",  filename_divider,   in_Prefix, ".csv", sep=""), quote=FALSE)
write.csv(table_classes_shp_750m,   paste(dir_path, "table_classes_shp_750m",  filename_divider,   in_Prefix, ".csv", sep=""), quote=FALSE)
write.csv(table_classes_shp_1000m,  paste(dir_path, "table_classes_shp_1000m", filename_divider,   in_Prefix, ".csv", sep=""), quote=FALSE)

#write.csv(table_classes_p100_50m,  paste(dir_path, "table_classes_p100_50m",  filename_divider,   in_Prefix, ".csv", sep=""), quote=FALSE)   # not with KNN-data
write.csv(table_classes_p100_250m,  paste(dir_path, "table_classes_p100_250m",  filename_divider,   in_Prefix, ".csv", sep=""), quote=FALSE)
write.csv(table_classes_p100_500m,  paste(dir_path, "table_classes_p100_500m",  filename_divider,   in_Prefix, ".csv", sep=""), quote=FALSE)
write.csv(table_classes_p100_750m,  paste(dir_path, "table_classes_p100_750m",  filename_divider,   in_Prefix, ".csv", sep=""), quote=FALSE)
write.csv(table_classes_p100_1000m, paste(dir_path, "table_classes_p100_1000m", filename_divider,   in_Prefix, ".csv", sep=""), quote=FALSE)


# creating summary tables for shp and p100 tables for later merging
#rowsum_matrix_shp_50m   <- matrix(v_colmeans_shp_50m,  1, length(headers2) #  not with KNN-data  6
rowsum_matrix_shp_250m  <- matrix(v_colmeans_shp_250m,  1, length(headers2)) 
rowsum_matrix_shp_500m  <- matrix(v_colmeans_shp_500m,  1, length(headers2)) 
rowsum_matrix_shp_750m  <- matrix(v_colmeans_shp_750m,  1, length(headers2))
rowsum_matrix_shp_1000m <- matrix(v_colmeans_shp_1000m, 1, length(headers2))

colnames(rowsum_matrix_shp_250m)  <-  headers2 
colnames(rowsum_matrix_shp_500m)  <-  headers2 
colnames(rowsum_matrix_shp_750m)  <-  headers2 
colnames(rowsum_matrix_shp_1000m) <-  headers2 


#rowsum_matrix_p100_50m   <-  matrix(v_colmeans_p100_50m,  1, length(headers2))   # not with KNN-data
rowsum_matrix_p100_250m  <-  matrix(v_colmeans_p100_250m,  1, length(headers2))
rowsum_matrix_p100_500m  <-  matrix(v_colmeans_p100_500m,  1, length(headers2))
rowsum_matrix_p100_750m  <-  matrix(v_colmeans_p100_750m,  1, length(headers2))
rowsum_matrix_p100_1000m <-  matrix(v_colmeans_p100_1000m, 1,  length(headers2))

colnames(rowsum_matrix_p100_250m)  <- headers2  
colnames(rowsum_matrix_p100_500m)  <- headers2  
colnames(rowsum_matrix_p100_750m)  <- headers2  
colnames(rowsum_matrix_p100_1000m) <- headers2  

 
# merging summaries for shp and p100 tables
#rowsum_matrix_50m <-   matrix(v_colmeans_shp_50m, 1, length(colnames(table_classes_shp_50m)))  # not with KNN-data
#rowsum_matrix_50m <-   rbind(rowsum_matrix_50m, v_colmeans_p100_50m)  # not with KNN-data
rowsum_matrix_250m <-  matrix(v_colmeans_shp_250m, 1, length(headers2))
rowsum_matrix_250m <-  rbind(rowsum_matrix_250m, v_colmeans_p100_250m)
rowsum_matrix_500m <-  matrix(v_colmeans_shp_500m, 1, length(headers2))
rowsum_matrix_500m <-  rbind(rowsum_matrix_500m, v_colmeans_p100_500m)
rowsum_matrix_750m <-  matrix(v_colmeans_shp_750m, 1, length(headers2))
rowsum_matrix_750m <-  rbind(rowsum_matrix_750m, v_colmeans_p100_750m)
rowsum_matrix_1000m <- matrix(v_colmeans_shp_1000m, 1, length(headers2))
rowsum_matrix_1000m <- rbind(rowsum_matrix_1000m, v_colmeans_p100_1000m)

#rownames(rowsum_matrix_50m)[1]   <- "shp means"   # not with KNN-data
#rownames(rowsum_matrix_50m)[2]   <- "p100 means"  # not with KNN-data
rownames(rowsum_matrix_250m)[1]  <- "shp means"
rownames(rowsum_matrix_250m)[2]  <- "p100 means"
rownames(rowsum_matrix_500m)[1]  <- "shp means"
rownames(rowsum_matrix_500m)[2]  <- "p100 means"
rownames(rowsum_matrix_750m)[1]  <- "shp means"
rownames(rowsum_matrix_750m)[2]  <- "p100 means"
rownames(rowsum_matrix_1000m)[1] <- "shp means"
rownames(rowsum_matrix_1000m)[2] <- "p100 means"


# calculating values for percent comparisin tables
#rowsum_matrix_50m_pc   <- matrix("", 2,  length(headers)+1) # not with KNN-data
rowsum_matrix_250m_pc  <- matrix("", 2, length(headers)+1)
rowsum_matrix_500m_pc  <- matrix("", 2, length(headers)+1)
rowsum_matrix_750m_pc  <- matrix("", 2, length(headers)+1)
rowsum_matrix_1000m_pc <- matrix("", 2, length(headers)+1)

for(ind in 1:length(headers)) {
 #rowsum_matrix_50m_pc[1,ind]  <- round((rowsum_matrix_50m[1, ind] /  rowsum_matrix_50m[1, length(headers)+1])*100, 0) #2  # not with KNN-data
 #rowsum_matrix_50m_pc[2,ind]  <- round((rowsum_matrix_50m[2, ind] /  rowsum_matrix_50m[2, length(headers)+1])*100, 0)	# not with KNN-data
 rowsum_matrix_250m_pc[1,ind]  <- round((rowsum_matrix_250m[1, ind] /  rowsum_matrix_250m[1,  length(headers)+1])*100, 0) #2
 rowsum_matrix_250m_pc[2,ind]  <- round((rowsum_matrix_250m[2, ind] /  rowsum_matrix_250m[2,  length(headers)+1])*100, 0)
 rowsum_matrix_500m_pc[1,ind]  <- round((rowsum_matrix_500m[1, ind] /  rowsum_matrix_500m[1,  length(headers)+1])*100, 0)
 rowsum_matrix_500m_pc[2,ind]  <- round((rowsum_matrix_500m[2, ind] /  rowsum_matrix_500m[2,  length(headers)+1])*100, 0)
 rowsum_matrix_750m_pc[1,ind]  <- round((rowsum_matrix_750m[1, ind] /  rowsum_matrix_750m[1,  length(headers)+1])*100, 0)
 rowsum_matrix_750m_pc[2,ind]  <- round((rowsum_matrix_750m[2, ind] /  rowsum_matrix_750m[2,  length(headers)+1])*100, 0)
 rowsum_matrix_1000m_pc[1,ind] <- round((rowsum_matrix_1000m[1, ind] / rowsum_matrix_1000m[1, length(headers)+1])*100, 0)
 rowsum_matrix_1000m_pc[2,ind] <- round((rowsum_matrix_1000m[2, ind] / rowsum_matrix_1000m[2, length(headers)+1])*100, 0)
}


# decorating shp- and p100 percent comparison tables with column- and rovnames and removing colsums colomn 
#colnames(rowsum_matrix_50m_pc)   <- headers2    # not with KNN-data
#rownames(rowsum_matrix_50m_pc)   <- rownames(rowsum_matrix_50m)  # not with KNN-data
colnames(rowsum_matrix_250m_pc)  <- headers2  
rownames(rowsum_matrix_250m_pc)  <- rownames(rowsum_matrix_250m)
colnames(rowsum_matrix_500m_pc)  <- headers2  
rownames(rowsum_matrix_500m_pc)  <- rownames(rowsum_matrix_500m)
colnames(rowsum_matrix_750m_pc)  <- headers2  
rownames(rowsum_matrix_750m_pc)  <- rownames(rowsum_matrix_750m)
colnames(rowsum_matrix_1000m_pc) <- headers2  
rownames(rowsum_matrix_1000m_pc) <- rownames(rowsum_matrix_1000m)

rowsum_matrix_250m_pc   <- subset(rowsum_matrix_250m_pc,   select = -colsums)
rowsum_matrix_500m_pc   <- subset(rowsum_matrix_500m_pc,   select = -colsums)
rowsum_matrix_750m_pc   <- subset(rowsum_matrix_750m_pc,   select = -colsums)
rowsum_matrix_1000m_pc  <- subset(rowsum_matrix_1000m_pc,  select = -colsums)

 
# saving CSV's of percent coomparison tables between shape and p100
#write.csv(rowsum_matrix_50m_pc, paste(dir_path,   "classes_shp_vs_p100_50m_pc.csv",   sep=""), quote=FALSE)  # not with KNN-data
write.csv(rowsum_matrix_250m_pc, paste(dir_path,  "classes_shp_vs_p100_250m_pc",  filename_divider,  in_Prefix, ".csv",  sep=""), quote=FALSE)
write.csv(rowsum_matrix_500m_pc, paste(dir_path,  "classes_shp_vs_p100_500m_pc",  filename_divider,  in_Prefix, ".csv",  sep=""), quote=FALSE)
write.csv(rowsum_matrix_750m_pc, paste(dir_path,  "classes_shp_vs_p100_750m_pc",  filename_divider,  in_Prefix, ".csv",  sep=""), quote=FALSE)
write.csv(rowsum_matrix_1000m_pc, paste(dir_path, "classes_shp_vs_p100_1000m_pc", filename_divider,  in_Prefix, ".csv", sep=""), quote=FALSE)

#write.csv(rowsum_matrix_50m, paste(dir_path,   "classes_shp_vs_p100_50m.csv",   sep=""), quote=FALSE)  # not with KNN-data
write.csv(rowsum_matrix_250m, paste(dir_path,  "classes_shp_vs_p100_250m",  filename_divider,  in_Prefix, ".csv",  sep=""), quote=FALSE)
write.csv(rowsum_matrix_500m, paste(dir_path,  "classes_shp_vs_p100_500m",  filename_divider,  in_Prefix, ".csv",  sep=""), quote=FALSE)
write.csv(rowsum_matrix_750m, paste(dir_path,  "classes_shp_vs_p100_750m",  filename_divider,  in_Prefix, ".csv",  sep=""), quote=FALSE)
write.csv(rowsum_matrix_1000m, paste(dir_path, "classes_shp_vs_p100_1000m", filename_divider,  in_Prefix, ".csv", sep=""), quote=FALSE)

v_nrow_shp  <- nrow(table_classes_shp_250m)
v_nrow_p100 <- nrow(table_classes_p100_250m)
v_ncol <- ncol(table_classes_shp_500m)


# creating percent-tables
#table_classes_shp_50m_pc <- matrix("",  v_nrow_shp-1, v_ncol)  # not with KNN-data
table_classes_shp_250m_pc  <- matrix("",  v_nrow_shp-1, v_ncol) 
table_classes_shp_500m_pc  <- matrix("",  v_nrow_shp-1, v_ncol) # -1
table_classes_shp_750m_pc  <- matrix("",  v_nrow_shp-1, v_ncol)
table_classes_shp_1000m_pc <- matrix("", v_nrow_shp-1, v_ncol)

#table_classes_p100_50m_pc <- matrix("",  v_nrow_p100-1, v_ncol)  # not with KNN-data
table_classes_p100_250m_pc  <- matrix("",  v_nrow_p100-1, v_ncol)
table_classes_p100_500m_pc  <- matrix("",  v_nrow_p100-1, v_ncol)
table_classes_p100_750m_pc  <- matrix("",  v_nrow_p100-1, v_ncol)
table_classes_p100_1000m_pc <- matrix("", v_nrow_p100-1, v_ncol)


# calculating values in percnt-tables
for (irow in 1:(v_nrow_shp-1)) { #-1
 for (icol in 1:(v_ncol-1)) {
  #table_classes_shp_50m_pc[irow, icol] <- round((table_classes_shp_50m[irow, icol] /    table_classes_shp_50m[irow, v_ncol]), 2)*100 #2 # not with KNN-data	 
  table_classes_shp_250m_pc[irow, icol]  <- round((table_classes_shp_250m[irow, icol] /   table_classes_shp_250m[irow, v_ncol]), 2)*100 #2
  table_classes_shp_500m_pc[irow, icol]  <- round((table_classes_shp_500m[irow, icol] /   table_classes_shp_500m[irow, v_ncol]), 2)*100
  table_classes_shp_750m_pc[irow, icol]  <- round((table_classes_shp_750m[irow, icol] /   table_classes_shp_750m[irow, v_ncol]), 2)*100
  table_classes_shp_1000m_pc[irow, icol] <- round((table_classes_shp_1000m[irow, icol] / table_classes_shp_1000m[irow, v_ncol]), 2)*100
 }
}

for (irow in 1:(v_nrow_p100-1)) { #-1
 for (icol in 1:(v_ncol-1)) {
  #table_classes_p100_50m_pc[irow, icol]  <- round((table_classes_p100_50m[irow, icol] /    table_classes_p100_50m[irow, v_ncol]), 2)*100	 # not with KNN-data
  table_classes_p100_250m_pc[irow, icol]  <- round((table_classes_p100_250m[irow, icol] /   table_classes_p100_250m[irow, v_ncol]), 2)*100
  table_classes_p100_500m_pc[irow, icol]  <- round((table_classes_p100_500m[irow, icol] /   table_classes_p100_500m[irow, v_ncol]), 2)*100
  table_classes_p100_750m_pc[irow, icol]  <- round((table_classes_p100_750m[irow, icol] /   table_classes_p100_750m[irow, v_ncol]), 2)*100
  table_classes_p100_1000m_pc[irow, icol] <- round((table_classes_p100_1000m[irow, icol] / table_classes_p100_1000m[irow, v_ncol]), 2)*100
 }
}


# decorating shp- and p100 percent-tables with column- and rovnames. Addind summay at bottom
#colnames(table_classes_shp_50m_pc) <- headers2 #  list(NULL, colnames(table_classes_shp_50m))  # not with KNN-data
colnames(table_classes_shp_250m_pc)  <- headers2  
colnames(table_classes_shp_500m_pc)  <- headers2  
colnames(table_classes_shp_750m_pc)  <- headers2  
colnames(table_classes_shp_1000m_pc) <- headers2  

#rownames(table_classes_shp_50m_pc)  <- get_shape_rownames(shp, in_ShpAttribute)  #  eval(parse(text=shp_name_full[2]))  # not with KNN-data
rownames(table_classes_shp_250m_pc)  <- get_shape_rownames(shp, in_ShpAttribute)  # eval(parse(text=shp_name_full[2]))
rownames(table_classes_shp_500m_pc)  <- get_shape_rownames(shp, in_ShpAttribute)  # eval(parse(text=shp_name_full[2]))
rownames(table_classes_shp_750m_pc)  <- get_shape_rownames(shp, in_ShpAttribute)  # eval(parse(text=shp_name_full[2]))
rownames(table_classes_shp_1000m_pc) <- get_shape_rownames(shp, in_ShpAttribute)  # eval(parse(text=shp_name_full[2]))

#colnames(table_classes_p100_50m_pc) <-  headers2    # not with KNN-data
colnames(table_classes_p100_250m_pc)  <- headers2  
colnames(table_classes_p100_500m_pc)  <- headers2  
colnames(table_classes_p100_750m_pc)  <- headers2  
colnames(table_classes_p100_1000m_pc) <- headers2  


# saving CSV's of shape-points and random points for classes of different raduis
#write.csv(table_classes_shp_50m_pc,    paste(dir_path, "classes_shp_50m_pc",   filename_divider,   in_Prefix, ".csv",  sep=""), quote=FALSE)  # not with KNN-data
write.csv(table_classes_shp_250m_pc,   paste(dir_path, "classes_shp_250m_pc",  filename_divider,   in_Prefix, ".csv",  sep=""), quote=FALSE)
write.csv(table_classes_shp_500m_pc,   paste(dir_path, "classes_shp_500m_pc",  filename_divider,   in_Prefix, ".csv",  sep=""), quote=FALSE)
write.csv(table_classes_shp_750m_pc,   paste(dir_path, "classes_shp_750m_pc",  filename_divider,   in_Prefix, ".csv",  sep=""), quote=FALSE)
write.csv(table_classes_shp_1000m_pc,  paste(dir_path, "classes_shp_1000m_pc", filename_divider,  in_Prefix, ".csv",  sep=""), quote=FALSE)

#write.csv(table_classes_p100_50m_pc,   paste(dir_path, "classes_p100_50m_pc",   filename_divider,   in_Prefix, ".csv", sep=""), quote=FALSE)  # not with KNN-data
write.csv(table_classes_p100_250m_pc,  paste(dir_path, "classes_p100_250m_pc",  filename_divider,   in_Prefix, ".csv", sep=""), quote=FALSE)
write.csv(table_classes_p100_500m_pc,  paste(dir_path, "classes_p100_500m_pc",  filename_divider,   in_Prefix, ".csv", sep=""), quote=FALSE)
write.csv(table_classes_p100_750m_pc,  paste(dir_path, "classes_p100_750m_pc",  filename_divider,   in_Prefix, ".csv", sep=""), quote=FALSE)
write.csv(table_classes_p100_1000m_pc, paste(dir_path, "classes_p100_1000m_pc", filename_divider,  in_Prefix, ".csv", sep=""), quote=FALSE)

stat_matr3 = matrix(stat_matr[1,], 1, 9)
colnames(stat_matr3) <- matr_cols

for(ind in 2:length(stat_matr[,1])) {
	#if(stat_matr[ind ,7] < 0.06) {
	if( (stat_matr[ind,7] < 0.06) && ((stat_matr[ind,9] == "two-sided") || (stat_matr[ind,9] == "two.sided")) ) { 
		stat_matr3 <- rbind(stat_matr3, stat_matr[ind,])
	}
} # for
stat_matr4 <- stat_matr3[order(stat_matr[,7], decreasing=FALSE)]
matrix_name3 <- paste("stattest_classes_relevant_", in_Prefix, ".csv", sep="")
write.matrix(stat_matr3, paste(dir_path, matrix_name3, sep=""), sep=",")
write.matrix(stat_matr4, paste(dir_path, matrix_name3, "_test_", sep=""), sep=",")

v_warn <- warnings()
cat("\n warnings \n")
print(v_warn)
rm(list = ls())


