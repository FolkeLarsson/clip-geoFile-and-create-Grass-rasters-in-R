


## runs from Grass in ordinary command prompt, called by Rscript
## cutting a GDAL imagefile and creating a new, creating GRASS raster and a reclassing of it based on a rules file  
## can also create another imagefile, like the first one but transformed to another coordinate system
## tested with rt90, sweref99 and GTiff 
## tested on Ubuntu 14-17, FreeBSD 11 and windows 10 
## arguments from command line or an ini-file with more parameters, so using command line one has to hardcode some arguments in script 
## can create new location based on settings of input imagefile or use existing
## paths to executables, imagefile type and class hardcoded in script for different OS
## Folke Larsson, Boden Sweden, 2018. Version 1.0


library(sp)
library(gdalUtils)
library(rgdal)
library(rgrass7)
library(raster)
library(grid)
library(stringr)
library(methods)
library(getArgsflaGrassStat) # for the arguments,  install_github("FolkeLarsson/getArgsFlaGrassStat")

# init 
options(warn=0) # 
read_from_ini_file <- TRUE
task <- "clip_and_create_rasters" # only for command line arguments
#options(error=recover)


##constants
rt90     <- "+proj=tmerc +lat_0=0 +lon_0=15.80827777777778 +k=1 +x_0=1500000 +y_0=0 +ellps=bessel +towgs84=414.1,41.3,603.1,-0.855,2.141,-7.023,0 +units=m +no_defs"
sweref99 <- "+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
utm34N   <- "+proj=utm +no_defs +zone=34 +a=6378137 +rf=298.257223563 +towgs84=0.000,0.000,0.000 +to_meter=1"
wgs84    <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"


### variables ###
# projwin_str: area of destination file, to cut out from given image
imagefile_type   <- "GTiff"
imagefile_class  <- "RasterBrick"
raster_nullvalue1 <- 65535   # like in forest maps, skogskartor, fom SLU
raster_nullvalue2 <- -1      # like in laser data from SVS
max_mem_size <- "500M" # for gdal_translate and gdalwarp
inifile_name <- "RGrass_parameters.ini"
inifile_path <- "/home/follar/grass7/maps/parameter_files"
df_arguments <- data.frame("arguments") # to hold parameters either from an ini-file or command line



### functions ###
# a fix if using older versions of R
paste0 <- function(..., collapse = NULL) {
	paste(..., sep = "", collapse = collapse)
}# paste0


# check if valid installation of gdal exists
gdal_OK = function() {
	gdal_setInstallation()
	valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
	return(valid_install) # TRUE or FALSE 
} # gdal_OK


# from given file, like a GTiff and given clip-window, creating new smaller file.
# returns name of new file if its a raster, else returns the other type
cut_mapfile = function(v_in_mapfile, v_file_sep, v_suffix, v_lrx, v_lry, v_ulx, v_uly, v_imagefile_type, v_imagefile_class) {
	projwin_str = paste(v_ulx, ", ", v_uly, ", ", v_lrx, ", ", v_lry, sep="")
	cat(projwin_str)
	cat("\n")
	in_file_length <- str_length(as.character(v_in_mapfile))
	cut_pos <- as.numeric(str_locate(v_in_mapfile, '\\.'))
	in_file_name <- substr(v_in_mapfile, 1, cut_pos-1)
	in_file_type <- substr(v_in_mapfile, cut_pos+1, in_file_length)
	dest_file_name <- paste(substr(v_in_mapfile, 1, nchar(v_in_mapfile)-4), v_suffix, sep="")
	dest_file <- paste(dest_file_name, paste(".", in_file_type, sep=""), sep="")
	can_cut_mapfile <- tryCatch( {
				gdal_translate(v_in_mapfile, dest_file, of=v_imagefile_type, projwin=projwin_str, output_Raster=TRUE,  verbose=TRUE) 
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
	if( class(can_cut_mapfile) == v_imagefile_class ) { 
		return(dest_file)
	} else {
		return(class(can_cut_mapfile))
	} # class
} # cut_mapfile 


# converting coordinate system and creating new raster-file 
# todo srcnodata, overwrite option  
warp_mapfile = function (v_in_srs, v_out_srs, v_srcnodata_value, v_dstnodata_value, v_imagefile_type, v_file_sep, v_in_mapfile, v_out_mapfilename) {
	cs_name <- as.character(v_out_srs)
	in_crs  <- eval(parse(text=v_in_srs))
	out_crs  <- eval(parse(text=v_out_srs))
	cut_pos <- as.numeric(str_locate(v_in_mapfile, '\\.'))
	print(cut_pos)
	in_file_name <- substr(v_in_mapfile, 1, cut_pos-1)
	in_file_length <- str_length(as.character(v_in_mapfile))
	in_file_type <- substr(v_in_mapfile, cut_pos+1, in_file_length)
	in_file_type_length <- nchar(in_file_type)
	dest_file_name <- paste(substr(v_in_mapfile, 1, nchar(v_in_mapfile)-4), "_", cs_name, sep="")
	dest_file <- paste(dest_file_name, paste(".", in_file_type, sep=""), sep="")		
	err_exc <- gdalwarp(v_in_mapfile, dest_file, t_srs=out_crs, s_srs=in_crs, of=v_imagefile_type, srcnodata = v_srcnodata_value, dstnodata=v_dstnodata_value, output_Raster=TRUE, wm=500, verbose=TRUE)
	print(err_exc)
	return(dest_file)
} # warp_mapfile


# create Grass raster from given image file, returns name of rastermap
create_grass_rastermap = function(v_mapfile, v_file_sep, v_location, v_mapset, v_is_new_location) {
	vect_filepath <- unlist(str_split(v_mapfile, v_file_sep))  # unlist to create vector from list
	l_filepath <- length(vect_filepath)
	raster_name <- paste(vect_filepath[l_filepath], sep="")
	raster_name_t <- substr(raster_name, 1, nchar(raster_name)-4)
	if(v_is_new_location == TRUE) {
		print(" ----location TRUE ----")
		curr_gisBase <- "/usr/lib/grass722/"
		init_grass <- initGRASS(gisBase=curr_gisBase, home = tempdir(), gisDbase=curr_gisDbase, location="temp", mapset="PERMANENT", override = TRUE)
		print(init_grass)
		err_exec <- execGRASS("r.in.gdal", input=v_mapfile, output=raster_name_t, location=v_location, flags=c("overwrite", "e"))
		print(err_exec)
		cat("\n")
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
create_grass_classes_rastermap = function(v_raster_name, v_rules_file, v_location, v_mapset) {
	raster_name_classes <- paste("classes_", v_raster_name, sep="")
	classes_title <- paste("Reclassed_", raster_name_classes, sep="")
	curr_gisBase <- "/usr/lib/grass722/"
	init_grass <- initGRASS(gisBase=curr_gisBase, home = tempdir(), gisDbase=curr_gisDbase, location=v_location, mapset=v_mapset, override = TRUE)
	print(init_grass)
	# todo catch errors
	err_exec <- execGRASS("r.reclass", input=raster_name, output=raster_name_classes, rules=v_rules_file, title=classes_title, flags=("overwrite"))
	if (err_exec == 0) {
		return(raster_name_classes)
	} else {
		stop( paste(" error creating grass classes raster map", as.character(err_exec)) )
	}
}# create_grass_classes_rastermap



### ---- the script itself ---- ### 


## setting parameters coming either from a ini-file or command prompt
#curr_task <- getArgsGstat$task_clip_create_rasters
df_arguments <- getArgsflaGrassStat::read_arguments(inifile_path, inifile_name, read_from_ini_file, task)
in_MapFile        <- df_arguments$mapfile       # source GTiff map file
in_ulx            <- df_arguments$ulx           # Upper-Left X-coordinate 
in_uly            <- df_arguments$uly           # Upper-Left Y-coordinate 
in_lrx            <- df_arguments$lrx           # Lower-Left X-coordinate 
in_lry            <- df_arguments$lry           # Lower-Left Y-coordinate 
in_Prefix         <- df_arguments$prefix        # directory name for log/error files
in_Suffix         <- df_arguments$suffix        # specific name of the new map, like area-name, added to given source filename 
in_ClassRules     <- df_arguments$classfile     # ASCII-file with rules for creating classes in Grass raster 
in_location       <- df_arguments$location      # new or existing Grass location
in_mapset         <- df_arguments$mapset        # existing Grass mapset, if new localion always PERMANENT
in_isNewLocation  <- df_arguments$isnewlocation # TRUE or FALSE if a new location should be created 
in_convert_crs    <- df_arguments$convert_crs   #  TRUE or FALSE if the cutted tif-file should be converted to IUT_CRS too
in_inCRS          <- df_arguments$incrs
in_outCRS         <- df_arguments$outcrs

## area of destination file, to cut out from given image
projwin_str     <- paste(in_ulx, ", ", in_uly, ", ", in_lrx, ", ", in_lry, sep="")

# directory for log- and error-files and some initial logging of system parameters
work_dir    <- paste("dir_", in_Prefix, sep="")
file_sep    <- .Platform$file.sep
file_psep   <- .Platform$path.sep
OS_type     <- .Platform$OS.type
os_version  <- version$os
system_name <- Sys.info()['sysname']
system_info <- str(Sys.info)


# listing received arguments  
print("-------- command line or ini-file arguments ------------- ")
print(df_arguments)
if (is.null(in_isNewLocation)) {
	stop(" no in_isNewLocation argument !! ")
}
print("-------- end command line arguments --------- ")


## create dir for logging 
could_create_dir <- dir.create(work_dir)
print(paste("  logging directory path: ", could_create_dir))
dir_path = paste(work_dir, file_sep, sep="")
tmp_path <- paste(dir_path, paste(file_sep, "tmp", sep=""), sep="")
print(paste(" dir: ", dir_path))


## directing output to files 
msg <- file(paste(dir_path, "errors_and_warnings.log", sep=""), open="wt")
out <- file(paste(dir_path, "logging_output.log" , sep=""), open="wt")
sink(msg, type="message")
sink(out, type="output")


## logging system parameters and arguments
print(paste(" projwin ",    projwin_str))
print(paste(" file sep: ",  file_sep))
print(paste(" file psep: ", file_psep))
print(paste(" OS type: ",   OS_type))
print(paste(" sysname: ",   system_name))
print(paste(" sysinfo: ",   system_info))


## Grass environmemt
if (read_from_ini_file == FALSE) {
	if (OS_type == "windows") {
		file_sep <- "\\\\" 
		curr_gisBase  <- "C:\\OSGeo4W64\\apps\\grass\\grass-6.4.3"
		curr_gisDbase <- "C:\\GIS\\grass\\grassdata"	
	} else {
		curr_gisBase  <- "/usr/local/grass722"
		curr_gisDbase <- "/home/folke/grass/grassdata/"	
	}# if-else
} else {
	if (OS_type == "windows") {
		file_sep <- "\\\\" 
		curr_gisBase  <- df_arguments$gisbase
		curr_gisDbase <- df_arguments$gisdbase	
	} else {
		curr_gisBase  <- df_arguments$gisbase
		curr_gisDbase <- df_arguments$gisdbase	
	}# if-else
}# if-else


### the script logic section 
is_gdal_OK <- gdal_OK()
print( paste("----- valid install of gdal  :", is_gdal_OK) )
if (is_gdal_OK == TRUE) {
	new_mapfile <- cut_mapfile(in_MapFile, file_sep, in_Suffix, in_lrx, in_lry, in_ulx, in_uly, imagefile_type, imagefile_class)
	print(paste("new map file: ", new_mapfile))
	raster_name <- create_grass_rastermap(new_mapfile, file_sep, in_location, in_mapset, in_isNewLocation)
	print(paste("raster name: ", raster_name))
	in_crs <- eval(parse(text=in_inCRS))
	out_crs <- eval(parse(text=in_outCRS))
	if ((in_convert_crs == TRUE) && (in_inCRS != in_outCRS)) {
		converted_file_name <- warp_mapfile(in_inCRS, in_outCRS, raster_nullvalue2, raster_nullvalue1, imagefile_type, file_sep, new_mapfile, in_outCRS)
		print(paste(" converted file name", converted_file_name))
	} # if
	raster_name_classes <- create_grass_classes_rastermap(raster_name, in_ClassRules, in_location, in_mapset) 
	print(paste("raster name classes: ", raster_name_classes))
} else {
	stop(" found no valid GDAL installation ")
} # is_gdal_OK








