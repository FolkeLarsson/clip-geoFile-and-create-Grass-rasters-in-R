## runs from Grass command prompt, called by Rscript
## cutting a GDAL imagefile and creating a new, creating GRASS raster and a reclassing of it based on a rules file  
## tested with rt90 and GTiff 
## tested on Ubuntu 14 and windows 10 
## arguments from command line 
## can create new location based on settings of input imagefile or use existing


# init
library(sp)
library(gdalUtils)
library(rgdal)
library(spgrass6)
library(raster)
library(grid)
library(stringr)
require(methods)

args<-commandArgs(TRUE)
options(warn=0)


# command line parameters 
in_MapFile        <- args[1]  # source GTiff map file
in_ulx            <- args[2]  # Upper-Left X-coordinate 
in_uly            <- args[3]  # Upper-Left Y-coordinate 
in_lrx            <- args[4]  # Lower-Left X-coordinate 
in_lry            <- args[5]  # Lower-Left Y-coordinate 
in_Prefix         <- args[6]  # directory name for log/error files
in_Suffix         <- args[7]  # specific name of the new map, like area-name, added to given source filename 
in_ClassRules     <- args[8]  # ASCII-file with rules for creating classes in Grass raster 
in_location       <- args[9]  # new or existing Grass location
in_mapset         <- args[10] # existing Grass mapset, if new localion always PERMANENT
in_isNewLocation  <- args[11] # TRUE or FALSE if a new location should be created 


# area of destination file, to cut out from given Gtiff map
projwin_str1 = paste(in_ulx, ", ", in_uly, ", ", in_lrx, ", ", in_lry, sep="")
print(" projwin ")
print(projwin_str1)


# variables 
imagefile_type <- "GTiff"
imagefile_class <- "RasterBrick"
#curr_location <- "rt90"
#curr_mapset   <- "AKH"
#new_location <- "rt90_AKH"


## functions

# a fix if using older versions of R
paste0 <- function(..., collapse = NULL) {
    paste(..., sep = "", collapse = collapse)
}

# check if valid installation of gdal exists
gdal_OK = function() {
	gdal_setInstallation()
	valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
	print( paste("----- valid install of gdal ? ------", valid_install) )
	return(valid_install) # TRUE or FALSE 
} # gdal_OK

# from given file, GTiff so far, and given clip-window, creating new smaller file.
# returns name of new file if its a raster, else returns the other type
cut_mapfile = function(v_in_mapfile, v_file_sep, v_suffix, v_lrx, v_lry, v_ulx, v_uly) {
	projwin_str = paste(v_ulx, ", ", v_uly, ", ", v_lrx, ", ", v_lry, sep="")
	in_file_length <- str_length(as.character(v_in_mapfile))
	cut_pos <- as.numeric(str_locate(v_in_mapfile, '\\.'))
	in_file_name <- substr(v_in_mapfile, 1, cut_pos-1)
	in_file_type <- substr(v_in_mapfile, cut_pos+1, in_file_length)
	dest_file_name <- paste(substr(v_in_mapfile, 1, nchar(v_in_mapfile)-4), v_suffix, sep="")
	dest_file <- paste(dest_file_name, paste(".", in_file_type, sep=""), sep="")
	can_cut_mapfile <- tryCatch( {
				gdal_translate(v_in_mapfile, dest_file, of=imagefile_type, projwin=projwin_str, output_Raster=TRUE,  verbose=TRUE) 
				#return_value <- dest_file
			},
			error=function(err) {
				stop(err)
			},
			warning=function( war) {
				message(" could not cut  off GTiff file with gdal_translate ")
			},
			finally={ 
				print( paste(" new  GTiff created: ", dest_file) )
			}
	) # can_cut_mapfile
	if( class(can_cut_mapfile) == imagefile_class ) { 
		return(dest_file)
	} else {
		return(class(can_cut_mapfile))
	} # class
} # cut_mapfile 

# create Grass raster from given image file, returns name of rastermap
create_grass_rastermap = function(v_mapfile, v_file_sep, v_location, v_mapset, v_is_new_location) {
	vect_filepath <- unlist(str_split(v_mapfile, v_file_sep))  # unlist to create vector from list
	l_filepath <- length(vect_filepath)
	raster_name <- paste(vect_filepath[l_filepath], sep="")
	raster_name_t <- substr(raster_name, 1, nchar(raster_name)-4)
	if(v_is_new_location == TRUE) {
		err_exec <- execGRASS("r.in.gdal", input=v_mapfile, output=raster_name_t, location=new_location, flags=c("overwrite", "e"))
	} else if(v_is_new_location == FALSE) {
		init_grass <- initGRASS(gisBase=curr_gisBase, home = tempdir(), gisDbase=curr_gisDbase, location=v_location, mapset=v_mapset, override = TRUE)
		print(" ---  init ----------")
		print(init_grass) 	
		err_exec <- execGRASS("r.in.gdal", input=v_mapfile, output=raster_name_t, flags=c("overwrite", "o"))
	} else {
		stop(" in create_grass_rastermap, v_is_new_location must be TRUE or FALSE !! ")
	} # v_is_new_location 
	 
	if (err_exec == 0) {
		return(raster_name_t)
	} else {
		stop( paste(" error creating grass raster  map", as.character(err_exec)) )
	} # err_exec
} # create_grass_rastermap

# reclassify a Grass rastermap from given rules file
create_grass_classes_rastermap = function(v_raster_name, v_rules_file) {
	raster_name_classes <- paste("classes_", v_raster_name, sep="")
	classes_title <- paste("Reclassed_", raster_name_classes, sep="")
	err_exec <- execGRASS("r.reclass", input=raster_name, output=raster_name_classes, rules=v_rules_file, title=classes_title, flags=("overwrite"))
	if (err_exec == 0) {
		return(raster_name_classes)
	} else {
		print( paste(" error creating grass classes raster map", as.character(err_exec)) )
		stop( paste(" error creating grass classes raster map", as.character(err_exec)) )
	}
}# create_grass_classes_rastermap


# directory for log- and error-files and some initial logging
work_dir <- paste("dir_", in_Prefix, sep="")
file_sep <-  .Platform$file.sep
file_psep <- .Platform$path.sep
OS_type <- .Platform$OS.type
print("----- file sep -----")
print(file_sep)
print(file_psep)
print(OS_type)
print(Sys.info()["sysname"])

could_create_dir <- dir.create(work_dir)
dir_path = paste(work_dir, file_sep, sep="")
print(" ********* dir create for logging and error/warnings ******** ")
print(could_create_dir)
print("----- logging directory path -----")
print(dir_path)

# directing output to files 
msg <- file(paste(dir_path, "errors_and_warnings.log", sep=""), open="wt")
out <- file(paste(dir_path, "logging_output.log" , sep=""), open="wt")
sink(msg, type="message")
sink(out, type="output")

print("command line arguments: ")
print(args)

# setting file separator if OS is Windows
# then init Grass paths depending on OS, only tested on windows 10 and Ubuntu 14 
print(file_sep)
if (OS_type == "windows") {
	file_sep <- "\\\\" 
	curr_gisBase  <- "C:\\OSGeo4W64\\apps\\grass\\grass-6.4.3"
	curr_gisDbase <- "C:\\GIS\\grass\\grassdata"	
} else {
	curr_gisBase  <- "/usr/lib/grass64/"
	curr_gisDbase <- "/home/folke/grass/grassdata/"	
}
print(file_sep)



if(in_isNewLocation == TRUE) {
	print(paste("new loc:  ", in_isNewLocation))
} else if(in_isNewLocation == FALSE) {
	print(paste("new loc:  ", in_isNewLocation))
} else {
	print(paste("new loc:  ", in_isNewLocation))
	stop("  boolean trouble in_isNewLocation")
}


# the "script itself"
print(" ---  GDAL installation OK ?----------")
is_gdal_OK <- gdal_OK()
print(is_gdal_OK)
if (is_gdal_OK == TRUE) {
	new_mapfile <- cut_mapfile(in_MapFile, file_sep, in_Suffix, in_lrx, in_lry, in_ulx, in_uly)
	print(" ---  new map ----------")
	print(new_mapfile)
	raster_name <- create_grass_rastermap(new_mapfile, file_sep, in_location, in_mapset, in_isNewLocation)
	print(raster_name)
 if (in_isNewLocation == FALSE) {
	#init_grass <- initGRASS(gisBase=curr_gisBase, home = tempdir(), gisDbase=curr_gisDbase, location=new_location, mapset="PERMANENT", override = TRUE)
	print(" ---  init ----------")
	#print(init_grass) 	
 }# is_new_location	
raster_name_classes <- create_grass_classes_rastermap(raster_name, in_ClassRules) 
	print(raster_name_classes)
} else {
	stop(" found no valid GDAL installation ")
} # is_gdal_OK




