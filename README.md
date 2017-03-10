# clip-geoFile-and-create-Grass-rasters-in-R
cuts a georeferenced imagefile into a smaller one, creates a Grass raster and a reclassification using categories from a rules file



## runs from Grass command prompt, called by Rscript
## cutting a GDAL imagefile and creating a new, creating GRASS raster and a reclassing of it based on a rules file  
## tested with rt90 and GTiff 
## tested on Ubuntu 14 and windows 10 
## arguments from command line 
## can create new location based on settings of input imagefile or use existing

# I have used maps from   "ftp://salix.slu.se/download/skogskarta/"
# One rules-file downloaded to repository, age_classes 

ToDo's: 
 ini-file instead or as an alternative to command line arguments
  check if the files exists
  check if imagefile is valid and of right type, GTiff for example
  check if clipped area less then original
 geo-statistics and histogram info from original and new imagefile in logging output 
