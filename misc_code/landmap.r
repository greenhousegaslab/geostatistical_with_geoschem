#-----------------------------------------------------------------------------------------------#
# FUNCTION: landmap.r										#
# PURRPOSE: Create a land map with land and ocean for inversions using OCO-2.  			#
# S. Miller, Nov 27, 2021									#
#												#
#-----------------------------------------------------------------------------------------------#

#----------#
# NOTES:   #
#----------#

	# Use a value of 1 for land and 0 for ocean/water

#----------------------#
# Required libraries   #
#----------------------#

	library(ncdf4)
	library(R.matlab)
	library(fields)

#--------------------------------------------------#
# Read in MERRA-2 file with land/water variables   #
#--------------------------------------------------#

	# Read in the MERRA2 file from netcdf
        merra2path <- "/home-4/smill191@jhu.edu/work/GEOSChem1/ExtData/GEOS_4x5/MERRA2/2015/01/"
        ncfile <- nc_open(paste(merra2path,"MERRA2.20150101.CN.4x5.nc4",sep=""))

	# Read in individual land and water cover variables
	frland  <- ncvar_get(ncfile,"FRLAND") # Fraction land cover
	# frocean <- ncvar_get(ncfile,"FROCEAN") # Fraction ocean cover
	# frlandic<- ncvar_get(ncfile,"FRLANDIC") # Fraction of land covered in ice


#------------------------------#
# Create land-water mask       #
#------------------------------#

	# Set a threshold of 50% land to include in land cover mask
	landmap <- matrix(0,nrow=nrow(frland),ncol=ncol(frland))

	landmap[frland >= 0.3] <- 1

	print(dim(landmap))

	# Plot
	image.plot(landmap)
	
#-----------------#
# Write to file   #
#-----------------#

        writeMat("~/work/smiller/data/OCO2_MIP/misc_inputs/landmap.mat",landmap=landmap);

#-----------------------------------------------------------------------------------------------#
# END OF SCRIPT 										#
#-----------------------------------------------------------------------------------------------#

