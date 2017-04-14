## runs from Grass command prompt, called by Rscript
## cutting a GDAL imagefile and creating a new, creating GRASS raster and a reclassing of it based on a rules file  
## tested with rt90 and GTiff 
## tested on Ubuntu 14 and windows 10 
## arguments from command line 
## can create new location based on settings of input imagefile or use existing
## paths to executables, imagefile type and class hardcoded in script for windows and else
## Folke Larsson, Boden Sweden, 2017


library(sp)
library(gdalUtils)
library(rgdal)
library(spgrass6)
library(raster)
library(grid)
library(stringr)


# init 
args<-commandArgs(TRUE)
options(warn=0) # 


# command line arguments 
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


## variables
# projwin_str: area of destination file, to cut out from given image
projwin_str     <- paste(in_ulx, ", ", in_uly, ", ", in_lrx, ", ", in_lry, sep="")
imagefile_type  <- "GTiff"
imagefile_class <- "RasterBrick"


## functions
# a fix if using older versions of R
paste0 <- function(..., collapse = NULL) {
    paste(..., sep = "", collapse = collapse)
}

# check if valid installation of gdal exists
gdal_OK = function() {
	gdal_setInstallation()
	valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
	return(valid_install) # TRUE or FALSE 
} # gdal_OK

# from given file, like a GTiff and given clip-window, creating new smaller file.
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
		print(" ---  end init ----------")
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
		stop( paste(" error creating grass classes raster map", as.character(err_exec)) )
	}
}# create_grass_classes_rastermap


# directory for log- and error-files and some initial logging of system parameters
work_dir    <- paste("dir_", in_Prefix, sep="")
file_sep    <- .Platform$file.sep
file_psep   <- .Platform$path.sep
OS_type     <- .Platform$OS.type
os_version  <- version$os
system_name <- Sys.info()['sysname']
system_info <- str(Sys.info)

# for command prompt 
print("-------- command line arguments ------------- ")
print(args)
if (is.na(in_isNewLocation)) {
	stop(" no in_isNewLocation argument !! ")
}
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

# logging system parameters and arguments
print(paste(" projwin ",    projwin_str))
print(paste(" file sep: ",  file_sep))
print(paste(" file psep: ", file_psep))
print(paste(" OS type: ",   OS_type))
print(paste(" sysname: ",   system_name))
print(paste(" sysinfo: ",   system_info))


# the "script itself"
is_gdal_OK <- gdal_OK()
print( paste("----- valid install of gdal  :", is_gdal_OK) )
if (is_gdal_OK == TRUE) {
	new_mapfile <- cut_mapfile(in_MapFile, file_sep, in_Suffix, in_lrx, in_lry, in_ulx, in_uly)
	print(paste("new map file: ", new_mapfile))
	raster_name <- create_grass_rastermap(new_mapfile, file_sep, in_location, in_mapset, in_isNewLocation)
	print(paste("raster name: ", raster_name))
    raster_name_classes <- create_grass_classes_rastermap(raster_name, in_ClassRules) 
	print(paste("raster name classes: ", raster_name_classes))
} else {
	stop(" found no valid GDAL installation ")
} # is_gdal_OK




