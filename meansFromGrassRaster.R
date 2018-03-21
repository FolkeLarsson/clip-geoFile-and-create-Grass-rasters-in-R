# TODO: Add comment
# 
# Author: Folke Larsson Boden, Sweden 2018 
###############################################################################


## runs from Grass command prompt, or the ordinary one, called by Rscript
## calculate mean values from different circular areas and testing them against values from random points 
## two different means are calculated, one from original map-file and one setting a minimum value removing low values, 
## in this case open areas
## for statistical testing, non-parametric methods are used, ks.test and wilcox.test 
## tested with rt90, sweref99 and GTiff 
## tested on Ubuntu 14-17, FreeBSD10-11 and windows 10 
## arguments from command line
## uses selected location and mapset
## paths to executables hardcoded in script for windows and else
## Folke Larsson Boden, Sweden 2018, version 1.0 

library(methods)
library(XML)
library(sp)
library(rgdal)
library(rgeos)
library(maptools)
library(rgrass7)
library(raster)
library(MASS)
library(stringr)
library(getArgsflaGrassStat) # arguments from command line or ini-file ## install_github("FolkeLarsson/getArgsFlaGrassStat")


## init
options(warn=0)
#options(error=recover)

## constants, table columns and shortened version of rt90 proj4 
rt90 =  "+proj=tmerc +lat_0=0 +lon_0=15.80827777777778 +k=1 +x_0=1500000 +y_0=0 +ellps=bessel +units=m +no_defs"
sweref99 <- "+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
wgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
utm34N   <- "+proj=utm +no_defs +zone=34 +a=6378137 +rf=298.257223563 +towgs84=0.000,0.000,0.000 +to_meter=1"
test_method1 <- "ks"
test_method2 <- "wc"
curr_unit <- "m"



read_from_inifile <- TRUE
task <- "calculate_means"
inifile_name <- "RGrass_parameters_means.ini"
inifile_path <- "/home/follar/grass7/maps/parameter_files"
df_arguments <- data.frame("arguments") # to hold parameters either from an ini-file or command line
df_arguments <- getArgsflaGrassStat::read_arguments(inifile_path, inifile_name, read_from_inifile, task)


## arguments from command line
in_MinLimit     <- df_arguments$minlimit        # for adjusting min-value in maps from 0 to a defined one, in this case excluding open areas, clearcuts etc
in_Prefix       <- df_arguments$prefix          # directory name for both logging dir and csv-files etc
in_GrassRaster  <- df_arguments$grassraster     # name of raster map from Grass 
in_ShapeFile    <- df_arguments$shapefile       # path and name of shape file 
in_Location     <- df_arguments$location        # existing Grass location
in_Mapset       <- df_arguments$mapset          # existing Grass mapset
in_ShpAttribute <- df_arguments$shapeattribute  # name of attribute in shape layer, not mandatory if name is in shape file 
curr_gisBase    <- df_arguments$gisbase         # Grass libraries
curr_gisDbase   <- df_arguments$gisdbase        # Grass home for locations
in_Radiusseq    <- df_arguments$radiusseq       # start, step and end for radiuses
in_Radiuses     <- df_arguments$radiuses        # list with radiuses


# variables
filename_divider <- "_-_"
current_CRS <- sweref99
point_type1 <- "shape"
point_type2 <- "p100"
area_type1 <- "OA"    # calculations including Open Areas like agricultural land, clear cuttings etc
area_type2 <- "noOA"  # calculations without those Open Areas
matrix_type_sum <- "shape_vs_p100_means"
matrix_type <- "mean_matrix"
matrix_type_percentiles <- "percentiles"
radiuses <- c("125m", "250m", "500m", "750m", "1000m")

## directory for log- and error-files and some initial system parameters
work_dir <- paste("dir_", in_Prefix, sep="")
file_sep <- .Platform$file.sep
file_psep   <- .Platform$path.sep
system_name <- Sys.info()['sysname']
system_info <- str(Sys.info)
OS_type     <- .Platform$OS.type
os_version  <- version$os


## listing received arguments 
print("-------- ini-file or command line arguments ------------- ")
print(df_arguments)
print("-------- end ini-file or command line arguments --------- ")


## setting file separator if OS is Windows
## then init Grass paths depending on OS, tested on windows 10 and Ubuntu 14-17 and Freebsd10-11 

if( read_from_inifile == FALSE ) {
	if (OS_type == "windows") {
		file_sep <- "\\\\" 
		curr_gisBase  <- "C:\\OSGeo4W64\\apps\\grass\\grass-6.4.3"
		curr_gisDbase <- "C:\\GIS\\grass\\grassdata"	
	} else {
		curr_gisBase  <- "/usr/lib/grass722/"
		curr_gisDbase <- "/home/follar/grass7/grassdata/"	
	}# if else
	
}# if


## create dir for logging and naming csv-files
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


###   functions    ###

## replacing values lower then given min value with NA
set_min_value_to_zero <- function(x) { 
	x[x<as.numeric(in_MinLimit)] <- 0 ; return(x)
}

set_zero_to_NA <- function(x) { 
	x[x<1] <- NA; return(x)
}


# returna a vector with radiuses and units, like "250m" from either a vector with radiuses or sequence parameters 
get_radius_vect <- function(v_radius_row, v_unit, v_is_sequence) {
	r_vect <- unlist(strsplit(v_radius_row, ","))
	if(v_is_sequence == TRUE) {
		r_vect <- unlist(r_vect)
		curr_seq <- seq(as.numeric(r_vect[1]), as.numeric(r_vect[2]), as.numeric(r_vect[3]))
		n_rows <- length(curr_seq)
		r_vect <- c()
		for(ind in 1:n_rows) {
			r_vect[ind] <- paste(curr_seq[ind], v_unit, sep="")	
		}# for	
	} else if(v_is_sequence == FALSE) {
		n_rows <- length(r_vect)
		for(ind in 1:n_rows) {
			r_vect[ind] <- paste(r_vect[ind], v_unit, sep="")
		}
		r_vect <- unlist(r_vect)
	} else {
		stop("-- v_is_sequence has to be TRUE or FALSE--")
	}
	return(r_vect)
}# get_radius_vect


## returning rownames from attribute in shapefile, attribute "name" choosen or second attribute selected if no value given
## ToDo: checck for valid shapefile object
get_shape_rownames = function(v_shp, v_AttributeName) {
if(read_from_inifile == FALSE) { # temporary fix
	if(is.na(v_AttributeName)) {
		v_AttributeName <- NULL
	}
}# if
	row_names <- "uninitialized"
	shp_attr_vect <- unlist(names(v_shp))
	l_shp_attr_vect <- length(shp_attr_vect)
	if(is.null(v_AttributeName)) { # no given attribute in layer
		row_names_attr <- match("name", shp_attr_vect) # check for standard attribute "name"
		if(!is.null(row_names_attr))
			row_names <- paste("v_shp@data$", "name", sep="")
		else { # returning first attribute
			shp_attr_vect <- unlist(names(v_shp))
			l_shp_attr_vect <- length(shp_attr_vect)
			if((l_shp_attr_vect > 0) && (!is.null(shp_attr_vect)) ) {
				row_names <- paste("v_shp@data$", shp_attr_vect[1], sep="") # choose first attribute
			} else {
				print(paste("v_shp@data$", shp_attr_vect[1], sep=""))
				stop(" no attribute in shapefile ")
			}
		}
	} else { # name of attribute given 
		row_names_attr <- match(v_AttributeName, shp_attr_vect)
		if(!is.null(row_names_attr))
			row_names <- paste("v_shp@data$", v_AttributeName, sep="")
		else {
			print(paste(" in func: ", "v_shp@data$", v_AttributeName, sep=""))
			stop(" given name doesnt match any attribute in shapefile, or no valid shapefile")
		}
	}
	return(eval(parse(text=row_names)))
} # get_shape_rownames



## wrapper for statistical calculations to make single test of two variables
stat_test = function(v_var1, v_var2, v_alternative, v_radius, v_method, v_area_type) {
	stat_vect = c()  # return vector
	var1_name <- v_var1$variable[1]
	var2_name <- v_var2$variable[1]
	var1_length <- length(v_var1[,1])
	var2_length <- length(v_var2[,1])
	var1_noNA <- na.omit(v_var1)
	var2_noNA <- na.omit(v_var2)
	var1_noNA_length <- length(var1_noNA[,1])
	var2_noNA_length <- length(var2_noNA[,1])
	
	if(v_method == "ks") {
		v_ks <- ks.test(var1_noNA[,1], var2_noNA[,1], alternative=v_alternative)
	} else if(v_method == "wc") {
		v_ks <- wilcox.test(var1_noNA[,1], var2_noNA[,1], alternative=v_alternative)
	} else {
		print(v_alternative)
		print(method)
		print(" NA NA NA no such method! ")
		stop()
		#v_ks <- NA  # should be caught
	}
	
	stat_vect[1] = var1_name
	stat_vect[2] = var1_noNA_length
	stat_vect[3] = var2_name
	stat_vect[4] = var2_noNA_length
	stat_vect[5] = paste(v_radius, "_", v_area_type)
	stat_vect[6] = round(v_ks$statistic, 2)
	stat_vect[7] = round(v_ks$p.value, 2)
	stat_vect[8] = v_ks$method
	stat_vect[9] = v_ks$alternative
	
	return(stat_vect)
} #stat_test



## wrapper for boxplot function to plot one object
f_boxplot = function(v_shp_var, v_p100_var, v_radius, v_cc_type, v_mean_crop, v_mean_nocrop, v_prefix, v_dir_path) {
	OA_type <- "noOA"  # todo parameter, NO
	line_color1 <- "red"
	line_color2 <- "green"
	mean_p100 <- mean(v_p100_var, na.rm = TRUE)
	p100_var_noNA <- na.omit(v_p100_var)
	#print(v_p100_var)
	#print(p100_var_noNA)
	if(v_cc_type == OA_type) {
		ylab_txt <- "mean values from raster without open areas"
		bp_name <- ", no clearcuts etc"
	} else {  # todo if-else
		ylab_txt <- "mean values from raster including open areas"
		bp_name <- ""
	}
	main_txt_row1 <- paste( "shape vs random points: ", v_radius, sep="") 
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


## calling f_boxplot to plot from two matrixes 
write_boxplots <- function(v_in_Prefix, v_dir_path, v_mean_nocrop, v_mean_crop, v_matrix1, v_matrix2) {
	nr_cols1 <- ncol(v_matrix1)
	nr_cols2 <- ncol(v_matrix2)
	nr_rows1 <- nrow(v_matrix1)
	nr_rows2 <- nrow(v_matrix2)
	# todo check if equal lengths in columns
	headers <- colnames(v_matrix1)
	
	for(ind in 1:nr_cols1) {
		col_factors <- unlist(str_split(headers[ind], " "))
		method <- col_factors[length(col_factors)]
		colname <- col_factors[length(col_factors)-1]
		v_matrix_row1 <- v_matrix1[, ind]
		v_matrix_row2 <- v_matrix2[, ind]
		f_boxplot(v_matrix_row1,  v_matrix_row2,  colname, method, mean_crop, mean_nocrop, v_in_Prefix, v_dir_path)
	}# for
	
}# write_boxplots 

## mean value for raster 
f_cs_mean = function(x) {
	cellStats(x, 'mean', na.rm=TRUE) 
}


## setting extent from geo-file in R-raster because readGrass6 returns raster with location extent 
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


## setting bounding box from geo-file in Grass-raster because readGrass6 returns raster with location bounding box 
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

## setting rownums in random points mean matrix
get_random_points_rownums <- function(v_random_points) {
	rownum_vect = c()
	n_coord_rows <- nrow(v_random_points@coords)
	row_seq <- seq(1, n_coord_rows, 1)
	return(row_seq)
}# get_random_points_rownums


## setting row index in stattest matrix that now uses 3 test methods, two-sided, greater and less 
f_row_index <- function(v_ind) {
	rowind <- 0
	suff <- v_ind-2
	if(suff >= 0){
		rowind <- 2*v_ind+suff
	}
	else {
		rowind <- v_ind*v_ind 
	}
	return(rowind)
}# index_factor


## create and return martix with different percentiles using quantile function 
calculate_percentiles <- function(v_matrix) {
	q_vect = c()
	p1 <- 0.10
	p2 <- 0.20
	p3 <- 0.80
	p4 <- 0.90
	p_vect <- c(p1, p2, p3, p4, 0, 0, 0, 0, 0)
	p_vect0 <- c(p1, p2, p3, p4)
	row_vect <- c(paste("p1:", p1, sep=""), paste("p2:", p2, sep=""), paste("p3:", p3, sep=""), paste("p4:", p4, sep=""))
	row_vect <- c(paste("p1:", p1, sep=""), paste("p2:", p2, sep=""), paste("p3:", p3, sep=""), paste("p4:", p4, sep=""))
#row_vect <- c(paste("p1:", p1, sep=""), paste("p2:", p2, sep=""), paste("p3:", p3, sep=""), paste("p4:", p4, sep=""))
	row_vect2 <- c("p1", "p2", "p3", "p4", "min", "max", "mean", "min_max_span", "q_span")
	colnames <- colnames(v_matrix)
	nr_cols <- ncol(v_matrix)
	nr_rows <- nrow(v_matrix)
	q_df <- data.frame(row.names=row_vect2, p_vect)
	for(ind in 1:nr_cols) {
		max_value   <- round(max(v_matrix[,ind], na.rm=TRUE), 0)
		min_value   <- round(min(v_matrix[1:(nr_rows-1),ind], na.rm=TRUE), 0)
		mean_value  <- round( mean(v_matrix[1:(nr_rows-1),ind], na.rm=TRUE)  , 0)
		minmax_span <- max_value - min_value
		q_list <- round(quantile(v_matrix[,ind], p_vect0, na.rm=TRUE), 0)
		q_span <- round(abs(as.numeric(q_list[1]) - as.numeric(q_list[4])), 0)
		q_list1 <- list(min_value, max_value, mean_value, minmax_span, q_span)
		q_df <- cbind(q_df,unlist(append(q_list, q_list1)))
	}
	q_df <- q_df[,2:(nr_cols+1)]
	colnames(q_df) <- colnames(v_matrix)
	return(q_df)
}# calculate_percentiles

## create and return matrix after extracting values from a raster in circles with different radiuses
calculate_means_matrix_with_different_radiuses <- function(v_radius_vector, v_raster, v_points, v_point_type, v_matrix_name, v_row_attribute) {
	mean_vect = c()
	column_vect = c()
	if(v_point_type == "shape") {
		curr_df <- data.frame(row.names = get_shape_rownames(v_points, v_row_attribute))
	} else {  # todo if-else
		n_rows <- nrow(v_points@coords)
		row_seq <- seq(1, n_rows, 1)
		curr_df <- data.frame(row.names = row_seq)
	} # if else
	
	for (ind in 1:length(v_radius_vector)) {
		column_name <- paste("r-", v_radius_vector[ind], " ", v_matrix_name, sep="")
		column_vect[ind] <- column_name
		curr_regexpr_pattern <- "[0-9]+"
		curr_regexpr <- regexpr(curr_regexpr_pattern, v_radius_vector[ind])
		curr_string_match <- as.numeric(regmatches(v_radius_vector[ind], curr_regexpr))
		curr_extracted_means <- round(extract(v_raster, v_points, buffer=curr_string_match, fun = mean, na.rm = TRUE), 1)
		curr_df <- cbind(curr_df, unlist(curr_extracted_means))
		mean_vect[ind] <- round(mean(unlist(curr_extracted_means), na.rm=TRUE),1)
		tull <- 1
	} # for
	colnames(curr_df) <- column_vect
	curr_df <- rbind(curr_df, mean_vect)
	rownames(curr_df)[nrow(curr_df)] <- "mean values"
	return(curr_df)
}# calculate_means_matrix_with_different_radiuses


## 
write_matrix_as_csv <- function(v_matrix_type, v_area_type, v_point_type, v_matrix, v_dir_path, v_filename_divider, v_in_Prefix) {
	write.csv(v_matrix, paste(v_dir_path, v_matrix_type, "_", v_area_type, "_", v_point_type, v_filename_divider,  v_in_Prefix, ".csv", sep=""), quote=FALSE)
}# write_matrix_as_csv


## creates matrix with rows from vector returned from stat_test
calculate_stattest_matrix <- function(v_radius_vector, v_point_vector1, v_point_vector2, v_point_type1, v_point_type2, v_test_method, v_area_type) {
	matr_cols = c("var1", "nVar1", "var2", "nVar2", "radius", "D/W", "p", "method", "alternative")
	alternative1 <- "two.sid"
	alternative2 <- "greater"
	alternative3 <- "less"
	test_v1 <- v_point_vector1[,1]
	test_v2 <- v_point_vector2[,1]
	test_v1 <- data.frame(test_v1)
	test_v2 <- data.frame(test_v2)
	
	test_v1$variable <- v_point_type1
	test_v2$variable <- v_point_type2
	dummy_row <- stat_test(test_v1, test_v2, alternative1, v_radius_vector[1], v_test_method, v_area_type)
	
	stat_matr = matrix(dummy_row, 1, length(dummy_row))
	colnames(stat_matr) <- matr_cols
	
	for (ind in 1:length(v_radius_vector)) {
		v_point_vector1 <- v_point_vector1[1:(nrow(v_point_vector1)-1),]
		v_point_vector2 <- v_point_vector2[1:(nrow(v_point_vector2)-1),]
		v_point_vector1_column <- v_point_vector1[,ind]
		v_point_vector2_column <- v_point_vector2[,ind]
		v_point_vector1_column <- data.frame(v_point_vector1_column)
		v_point_vector2_column <- data.frame(v_point_vector2_column)
		v_point_vector1_column$variable <- v_point_type1
		v_point_vector2_column$variable <- v_point_type2
		
		test_result <- stat_test(v_point_vector1_column, v_point_vector2_column, alternative1, v_radius_vector[ind], v_test_method, v_area_type)
		stat_matr <- rbind(stat_matr, test_result)
		test_result <- stat_test(v_point_vector1_column, v_point_vector2_column, alternative2, v_radius_vector[ind], v_test_method, v_area_type)
		stat_matr <- rbind(stat_matr, test_result)
		test_result <- stat_test(v_point_vector1_column, v_point_vector2_column, alternative3, v_radius_vector[ind], v_test_method, v_area_type)
		stat_matr <- rbind(stat_matr, test_result)
		
	}# for
	stat_matr <- stat_matr[2:nrow(stat_matr),]
	tull <- 1
	return(stat_matr)
	
} # calculate_stattest_matrix


## writing given stattest matrix in 3 versions
## one original, one ordered after p-value and one limited after given alternative and p-value
write_stattest_matrixes <- function(v_in_Prefix, v_dir_path, v_stattest_matrix, v_limit) {
	matrix_name <- paste("stattest_means_", v_in_Prefix, ".csv", sep="")
	write.matrix(v_stattest_matrix, paste(v_dir_path, matrix_name, sep=""), sep=",")
	
	# matrix ordered after increasing p-value
	matrix_name <- paste("stattest_means_", v_in_Prefix, "_ordered", ".csv", sep="")
	v_stattest_matrix <- v_stattest_matrix[order(v_stattest_matrix[,8], v_stattest_matrix[,7], decreasing=FALSE, v_stattest_matrix[,6]),]	
	write.matrix(v_stattest_matrix, paste(v_dir_path, paste(matrix_name,  sep=""), sep=""), sep=",")
	
	# matrix with p-value lower then given limit and with given alternative 
	v_stattest_matrix_rel <- matrix(v_stattest_matrix[1,], 1, ncol(v_stattest_matrix))
	for(ind in 2:length(v_stattest_matrix[,1])) {
		if( (v_stattest_matrix[ind,7] < as.numeric(v_limit)) && ((v_stattest_matrix[ind,9] == "two-sided") || (v_stattest_matrix[ind,9] == "two.sided")) ) {  
			v_stattest_matrix_rel <- rbind(v_stattest_matrix_rel, v_stattest_matrix[ind,])
		}# if
	} # for	
	matrix_name <- paste("stattest_means_", v_in_Prefix, "_relevant", ".csv", sep="")
	write.matrix(v_stattest_matrix_rel, paste(v_dir_path, paste(matrix_name, sep=""), sep=""), sep=",")
	
}# write_stattest_matrix


## adding one stattest matirx after another
concat_stattest_matrixes <- function(v_matrix1, v_matrix2) {
	col_names <- colnames(v_matrix1)
	nr_rows_matrix1 <- nrow(v_matrix1)
	nr_cols_matrix1 <- ncol(v_matrix1)
	new_matrix <- matrix(v_matrix1, nr_rows_matrix1, nr_cols_matrix1)
	new_matrix <- rbind(new_matrix, v_matrix2)
	return(new_matrix)
}# concat_stattest_matrixes

## adding one mean matirx after another
concat_mean_matrixes <- function(v_matrix1, v_matrix2) {
	nr_rows_matrix1 <- nrow(v_matrix1)
	nr_cols <- ncol(v_matrix2)
	columns1 <- as.character(colnames(v_matrix1))
	columns2 <- as.character(colnames(v_matrix2))
	columns <- rbind(unlist(columns1), unlist(columns2))
	columns_c <- cbind(unlist(columns1), unlist(columns2))
	return_matrix <- data.frame(cbind(v_matrix1, v_matrix2))
	colnames(return_matrix) <- columns_c
	return(return_matrix)
}# concat_mean_matrixes


## final matrix with rows for shape and random means for different radiuses 
create_matrix_shape_p100_means <- function(v_shape_matrix, v_p100_matrix, v_point_type1, v_point_type2) {
	nr_rows1 <- nrow(v_shape_matrix)
	nr_rows2 <- nrow(v_p100_matrix)
	return_matrix <- v_shape_matrix[nr_rows1,]
	return_matrix <- rbind(return_matrix, v_p100_matrix[nr_rows2,])
	rownames(return_matrix) <- c(paste(v_point_type1, " means", sep=""), paste(v_point_type2, " means", sep=""))
	return(return_matrix)
}# create_matrix_shape_p100_means


## logging system parameters and arguments 
print(paste(" file sep: ", file_sep))
print(paste(" file psep: ",file_psep))
print(paste(" OS type: ", OS_type))
print(paste(" OS version: ", os_version))
print(paste(" sysname: ", system_name))
print(paste(" sysinfo: ", system_info))     


## from command prompt
print(paste("input min value  ", in_MinLimit ))
print(paste("input map prefix: ", in_Prefix))
print(paste("input raster map: ", in_GrassRaster))
print(paste("input vector map, shapefile: ", in_ShapeFile))
print(paste("input vector map layer name, shapefile: ", in_ShpAttribute))
print(paste("input location: ", in_Location))
print(paste("input mapset: ", in_Mapset))


## setting Grass environment and logging the result 
init_grass <- initGRASS(gisBase=curr_gisBase, home = tempdir(), gisDbase=curr_gisDbase, location=in_Location, mapset=in_Mapset, override = TRUE)
cat(" ---  init ----------\n")
print(init_grass) 
cat(" ---  init end ---------- \n\n")


## making objects from Grass raster-map
rast_cla <- readRAST(in_GrassRaster, useGDAL=TRUE)  # readRAST6 for grass 6
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


## making objects from the shape-file and an extent to crop the raster from raster-file 
shp <- readShapePoints(in_ShapeFile, proj4string = CRS(current_CRS), verbose = TRUE, repair = TRUE)

## shp extent
shp_ext <- extent(shp)

## creating rasters cropped after shape-file 
rast_cla_r_cr <- crop(rast_cla_r, shp_ext)


## mean values for raters, original and cropped to shapefile, used in boxplots
mean_nocrop <- round(f_cs_mean(rast_cla_r), 0)
mean_crop <- round(f_cs_mean(rast_cla_r_cr), 0)
cat(" \n ------ mean values for raters, original and cropped to shapefile --- \n")
print(paste("mean value from geo-file: ", mean_nocrop))
print(paste("mean value from geo-file cropped to shapefile", mean_crop))
cat(" ------- end mean values for raters, original and cropped to shapefile --- \n")


## removing values lover then given min-level
rast_cla_r_nocc_temp <- calc(rast_cla_r, set_min_value_to_zero)  # curr did crop
rast_cla_r_nocc <- calc(rast_cla_r_nocc_temp, set_zero_to_NA)

rast_cla_r_cr_nocc_temp <- calc(rast_cla_r_cr, set_min_value_to_zero)  # curr did crop
rast_cla_r_cr_nocc <- calc(rast_cla_r_cr_nocc_temp, set_zero_to_NA)

## plotting raster histograms 
cat(" \n ------ histograms geo-file area and area cropped to shapefile and finally with open areas removed ------- \n")
png(paste(dir_path, "histogram_geofile_area _", in_Prefix, ".png", sep=""))
plot.new()
h0 <- hist(rast_cla_r, plot=TRUE, freq=TRUE)
print(h0)
dev.off()

png(paste(dir_path, "histogram_geofile_area_cropped_to_shapefile_", in_Prefix, ".png", sep=""))
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

## objects from a number of random points, less then 100 or 100 depending on OS 
random_points <- spsample(rast_cla, 100, "random") # makes random  number of random points, less then 100 
random_spatial_points = SpatialPoints(coordinates(random_points), CRS(current_CRS))


## to save coordinates for making other diagramms from several runs
write.matrix(random_spatial_points@coords, paste(dir_path, "p100_coordinates.csv", sep=""), sep=",")


## plot with shape points versus random points saved as a .png
png(paste(dir_path, "shp_and_random_points_mean_", in_Prefix, ".png", sep=""))
plot.new()
plot(random_spatial_points, col="red")
plot(shp, col="green", add=TRUE)
dev.off()

#radius_vect <- testrad2

if (read_from_inifile == TRUE) {
	#browser()
	if(!is.null(df_arguments$radiusseq)) {
		radius_vect <- get_radius_vect(in_Radiusseq, curr_unit, TRUE)
	} else if(!is.null(df_arguments$radiuses)) {
		radius_vect <- get_radius_vect(in_Radiuses, curr_unit, FALSE)
	} else {
		stop("-- both the vector with radiuses and the sequence cannot be NULL")
	} # if else
} else {
	radius_vect <- radiuses
} # if else	


## create matrixes with means from shapefile points with different radiuses , save as csv's
df_shape_matrix_OA   <- calculate_means_matrix_with_different_radiuses(radius_vect, rast_cla_r_cr, shp, point_type1, area_type1, in_ShpAttribute)
df_shape_matrix_noOA <- calculate_means_matrix_with_different_radiuses(radius_vect, rast_cla_r_cr_nocc, shp, point_type1, area_type2, in_ShpAttribute)
df_shape_matrix <- concat_mean_matrixes(df_shape_matrix_OA, df_shape_matrix_noOA)
write_matrix_as_csv(matrix_type, area_type1, point_type1, df_shape_matrix_OA, dir_path, filename_divider, in_Prefix)
write_matrix_as_csv(matrix_type, area_type2, point_type1, df_shape_matrix_noOA, dir_path, filename_divider, in_Prefix)
write_matrix_as_csv(matrix_type, paste(area_type1, "_and_", area_type2, sep=""), point_type1, df_shape_matrix, dir_path, filename_divider, in_Prefix)

## create matrixes with means from shapefile points with different radiuses , save as csv's
df_p100_matrix_OA   <- calculate_means_matrix_with_different_radiuses(radius_vect, rast_cla_r, random_spatial_points, point_type2, area_type1, in_ShpAttribute)
df_p100_matrix_noOA <- calculate_means_matrix_with_different_radiuses(radius_vect, rast_cla_r_nocc, random_spatial_points, point_type2, area_type2, in_ShpAttribute)
df_p100_matrix <- concat_mean_matrixes(df_p100_matrix_OA, df_p100_matrix_noOA)
write_matrix_as_csv(matrix_type, area_type1, point_type2, df_p100_matrix_OA, dir_path, filename_divider, in_Prefix)
write_matrix_as_csv(matrix_type, area_type2, point_type2, df_p100_matrix_noOA, dir_path, filename_divider, in_Prefix)
write_matrix_as_csv(matrix_type, paste(area_type1, "_and_", area_type2, sep=""), point_type2, df_p100_matrix, dir_path, filename_divider, in_Prefix)

## writing boxplots with shape vs random points
write_boxplots(in_Prefix, dir_path, mean_nocrop, mean_crop, df_shape_matrix_OA, df_p100_matrix_OA)
write_boxplots(in_Prefix, dir_path, mean_nocrop, mean_crop, df_shape_matrix_noOA, df_p100_matrix_noOA)

## statistical tests for shape vs random points
stattest_matrix_ks_OA <- calculate_stattest_matrix(radius_vect, df_shape_matrix_OA, df_p100_matrix_OA, point_type1, point_type2, test_method1, area_type1)
stattest_matrix_ks_noOA <- calculate_stattest_matrix(radius_vect, df_shape_matrix_noOA, df_p100_matrix_noOA, point_type1, point_type2, test_method1, area_type2)
stattest_matrix_wc_OA <- calculate_stattest_matrix(radius_vect, df_shape_matrix_OA, df_p100_matrix_OA, point_type1, point_type2, test_method2, area_type1)
stattest_matrix_wc_noOA <- calculate_stattest_matrix(radius_vect, df_shape_matrix_noOA, df_p100_matrix_noOA, point_type1, point_type2, test_method2, area_type2)

stattest_matrix <- concat_stattest_matrixes(stattest_matrix_ks_OA, stattest_matrix_ks_noOA)
stattest_matrix <- concat_stattest_matrixes(stattest_matrix, stattest_matrix_wc_OA)
stattest_matrix <- concat_stattest_matrixes(stattest_matrix, stattest_matrix_wc_noOA)

write_stattest_matrixes(in_Prefix, dir_path, stattest_matrix, 0.06)



## calculate matrixes from shape points with percentiles, min, max etc to see boundaries
percentiles_shape_OA <- calculate_percentiles(df_shape_matrix_OA)
percentiles_shape_noOA <- calculate_percentiles(df_shape_matrix_noOA)
percentiles_shape_matrix <- concat_mean_matrixes(percentiles_shape_OA, percentiles_shape_noOA)
write_matrix_as_csv(matrix_type_percentiles, area_type1, point_type1, percentiles_shape_OA, dir_path, filename_divider, in_Prefix)
write_matrix_as_csv(matrix_type_percentiles, area_type2, point_type1, percentiles_shape_noOA, dir_path, filename_divider, in_Prefix)
write_matrix_as_csv(matrix_type_percentiles, paste(area_type1, "_and_", area_type2, sep=""), point_type1, percentiles_shape_matrix, dir_path, filename_divider, in_Prefix)


## calculate matrixes from random points with percentiles, min, max etc to see boundaries 
percentiles_p100_OA <- calculate_percentiles(df_p100_matrix_OA)
percentiles_p100_noOA <- calculate_percentiles(df_p100_matrix_noOA)
percentiles_p100_matrix <- concat_mean_matrixes(percentiles_p100_OA, percentiles_p100_noOA)
write_matrix_as_csv(matrix_type_percentiles, area_type1, point_type2, percentiles_p100_OA, dir_path, filename_divider, in_Prefix)
write_matrix_as_csv(matrix_type_percentiles, area_type2,point_type2, percentiles_p100_noOA, dir_path, filename_divider, in_Prefix)
write_matrix_as_csv(matrix_type_percentiles, paste(area_type1, "_and_", area_type2, sep=""), point_type2, percentiles_p100_matrix, dir_path, filename_divider, in_Prefix)


## mean matrix with one row for shape-points and another for random points 
sum_matrix <- create_matrix_shape_p100_means(df_shape_matrix, df_p100_matrix, point_type1, point_type2)
write_matrix_as_csv(matrix_type_sum, paste(area_type1, "_and_", area_type2, sep=""), paste(point_type1, "_", point_type2, sep=""), sum_matrix, dir_path, filename_divider, in_Prefix)


v_warn <- warnings()
cat(" last warnings: \n")
print(v_warn)
rm(list = ls())




