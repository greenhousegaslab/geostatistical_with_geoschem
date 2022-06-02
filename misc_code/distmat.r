#-----------------------------------------------------------------------------------------------------------------------#
# FUNCTION: distmat.r 													#
# PURPOSE: Create a distance matrix that lists the distannces among grid points on the GEOS-Chem 4x5 grid		#
# S. Miller, Nov 27, 2021
#															#
#-----------------------------------------------------------------------------------------------------------------------#


#-------------------#
# FUNCTION NOTES:   #
#-------------------#

	# This script is constructed assuming that the rows of the flux matrix are longitudes and the columns are latitudes

#-------------------------------#
# Required libraries and paths  #
#-------------------------------#

	library(ncdf4)
	library(R.matlab)
	merra2path <- "/home-4/smill191@jhu.edu/work/GEOSChem1/ExtData/GEOS_4x5/MERRA2/2015/01/"
	source("~/scripts/utilities/greatcircledist.r")


#---------------------------------------------------------------------------------#
# Read in the land mask and select out the longitudes/latitudes of the land grid  #
#---------------------------------------------------------------------------------#

	ncfile <- nc_open(paste(merra2path,"MERRA2.20150101.CN.4x5.nc4",sep=""))
	lon1   <- ncvar_get(ncfile,"lon")
	lat1   <- ncvar_get(ncfile,"lat")

	# The first and last latitude value in the netcdf file is -90 and 90. That produces NA values in the distance matrix.
	# By contrast, the GEOS-Chem website lists -89 and 89 as the latitude limits.
	lat1[1]            <- -89
	lat1[length(lat1)] <- 89

        # Create two vectors that give the latitude and longitude, respectively, of the 
        # middle of each geographic grid box in the inversion domain.
        lat <- rep(lat1,each=length(lon1))
        lon <- rep(lon1,times=length(lat1))

        # Create distance matrix
                lat1 <- lat %*% t(rep(1,length(lat)))
                lat2 <- t(lat1)

                # Put zeros in the diagonals (otherwise, greatcircledist() will produce NaNs in some locations)
                diag(lat1)<-0
                diag(lat2)<-0

                # Convert to vectors
                lat1<-as.vector(lat1)
                lat2<-as.vector(lat2)

                lon1 <- lon %*% t(rep(1,length(lon)))
                lon2<-t(lon1)

                # Put zeros in diagonals
                diag(lon1)<-0
                diag(lon2)<-0

                # Convert to vectors
                lon1<-as.vector(lon1)
                lon2<-as.vector(lon2)

                deltamat<-greatcircledist(lat1=lat1,lat2=lat2,lon1=lon1,lon2=lon2)
                deltamat<-matrix(data=deltamat,nrow=length(lat),ncol=length(lat))
                deltamat<-deltamat/1000 # Convert distances from m to km
                # Deltamat contains the distances (in km) between every pair of grid cells
                # in the inverse modeling  domain.
                # Note that the function is not fooled by the dateline (and the switch from pos. to neg. longitudes)

		# Sometimes the greatcircledist function will produce NA values if lat1=42 and lat2=-42

                rm(lat1);rm(lat2);rm(lon1);rm(lon2);

	print("Distmat dimensions")
	print(dim(deltamat))

	# Write the distance matrix as a Matlab file
	# writeMat("~/data/smiller/outputs/xiaoling/distmat.mat",deltamat=deltamat);
	writeMat("~/work/smiller/data/OCO2_MIP/misc_inputs/deltamat.mat",deltamat=deltamat);


#-----------------------------------------------------------------------------------------------------------------------#
# END OF SCRIPT
#-----------------------------------------------------------------------------------------------------------------------#
	
